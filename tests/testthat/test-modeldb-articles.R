# Every article the site links to, and every article a model names, must be
# an .Rmd that pkgdown builds. pkgdown renders vignettes/<x>.Rmd and
# vignettes/articles/<x>.Rmd alike as articles/<x>.html, so both trees count.
#
# A model with no article must be listed, with a reason, in
# models-without-article.csv: `template` for generic building blocks that have
# no source to validate against, `article-pending` for published models that
# should have one. The list is explicit rather than a rule (such as "reference
# is 'nlmixr2lib template'") so that a NEW article-less model fails here until
# someone decides which it is.

.articleExemptStatus <- c("template", "article-pending")

# Stems of every .Rmd pkgdown builds into articles/.
.builtArticles <- function(root) {
  vign <- file.path(root, "vignettes")
  c(
    tools::file_path_sans_ext(list.files(vign, pattern = "\\.Rmd$")),
    tools::file_path_sans_ext(list.files(file.path(vign, "articles"), pattern = "\\.Rmd$"))
  )
}

# Every `href:` value anywhere in a parsed _pkgdown.yml. `[[` rather than `$`,
# which matches names partially on a list.
.collectHrefs <- function(x) {
  if (!is.list(x)) {
    return(character())
  }
  own <- if (is.character(x[["href"]])) x[["href"]] else character()
  c(own, unlist(lapply(x, .collectHrefs), use.names = FALSE))
}

# Article links in _pkgdown.yml that pkgdown cannot build: navbar hrefs into
# articles/, and the plain entries of each `articles:` group's `contents`
# (`articles/<x>` lives under vignettes/articles, a bare `<x>` under
# vignettes/). pkgdown selectors such as `starts_with("x")` are not file names
# and are skipped.
.danglingArticleLinks <- function(root) {
  ymlPath <- file.path(root, "_pkgdown.yml")
  yml <- yaml::read_yaml(ymlPath)
  built <- .builtArticles(root)
  hrefs <- grep("^articles/", .collectHrefs(yml), value = TRUE)
  out <- data.frame(where = character(), target = character(), stringsAsFactors = FALSE)
  for (h in hrefs) {
    if (!(sub("^articles/(.*)\\.html$", "\\1", h) %in% built)) {
      out[nrow(out) + 1L, ] <- list("navbar href", h)
    }
  }
  vign <- file.path(root, "vignettes")
  indexed <- character()
  for (grp in yml[["articles"]]) {
    for (entry in unlist(grp[["contents"]])) {
      if (grepl("(", entry, fixed = TRUE)) {
        next
      }
      indexed <- c(indexed, entry)
      rmd <- if (startsWith(entry, "articles/")) {
        file.path(vign, "articles", paste0(sub("^articles/", "", entry), ".Rmd"))
      } else {
        file.path(vign, paste0(entry, ".Rmd"))
      }
      if (!file.exists(rmd)) {
        out[nrow(out) + 1L, ] <- list(paste("articles group", grp[["title"]]), entry)
      }
    }
  }
  attr(out, "nHref") <- length(hrefs)
  attr(out, "indexed") <- indexed
  out
}

# Plain list items (`- <name>`, no `:` and no selector call) in the top-level
# `articles:` section of _pkgdown.yml, counted from the raw lines. This is an
# independent count of what .danglingArticleLinks() walks.
.articleEntryLineCount <- function(lines) {
  start <- grep("^articles:", lines)
  if (length(start) != 1L) {
    return(0L)
  }
  keys <- grep("^[A-Za-z_]", lines)
  end <- c(keys[keys > start], length(lines) + 1L)[[1]] - 1L
  sum(grepl("^[[:space:]]*-[[:space:]]+[^:(]+$", lines[seq_len(end - start) + start]))
}

# One message per model whose `vignette` names no built article, per
# article-less model missing from the exemption list, and per exemption that is
# stale or malformed. `built = NULL` skips the first check, which needs the
# source tree.
.modelArticleIssues <- function(db, exempt, built = NULL) {
  msg <- character()
  hasArticle <- !is.na(db$vignette)
  if (!is.null(built)) {
    bad <- hasArticle & !(db$vignette %in% built)
    msg <- c(msg, sprintf("%s: vignette '%s' is not an .Rmd under vignettes/", db$name[bad], db$vignette[bad]))
  }
  for (m in setdiff(db$name[!hasArticle], exempt$model)) {
    msg <- c(msg, sprintf("%s: no article and not listed in models-without-article.csv", m))
  }
  for (m in setdiff(exempt$model, db$name)) {
    msg <- c(msg, sprintf("%s: listed in models-without-article.csv but not in the modeldb", m))
  }
  for (m in intersect(exempt$model, db$name[hasArticle])) {
    msg <- c(msg, sprintf("%s: listed in models-without-article.csv but has an article; remove the row", m))
  }
  for (m in exempt$model[duplicated(exempt$model)]) {
    msg <- c(msg, sprintf("%s: listed twice in models-without-article.csv", m))
  }
  for (i in which(!(exempt$status %in% .articleExemptStatus))) {
    msg <- c(
      msg,
      sprintf(
        "%s: status '%s' is not one of: %s",
        exempt$model[[i]],
        exempt$status[[i]],
        paste(.articleExemptStatus, collapse = ", ")
      )
    )
  }
  for (m in exempt$model[!nzchar(trimws(exempt$reason))]) {
    msg <- c(msg, sprintf("%s: listed in models-without-article.csv without a reason", m))
  }
  msg
}

.readArticleExemptions <- function(path = testthat::test_path("models-without-article.csv")) {
  utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = character())
}

test_that("every _pkgdown.yml article link resolves to an article pkgdown builds", {
  skip_if_not_installed("yaml")
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  skip_if(!file.exists(file.path(root, "_pkgdown.yml")), "_pkgdown.yml not present (installed package)")
  skip_if(!dir.exists(file.path(root, "vignettes", "articles")), "vignettes/articles not present (installed package)")

  dangling <- .danglingArticleLinks(root)
  # The YAML walk must see every href and every `articles:` entry the file
  # holds; a parse that silently found none would pass the check below on an
  # empty set.
  lines <- readLines(file.path(root, "_pkgdown.yml"), warn = FALSE)
  expect_equal(attr(dangling, "nHref"), sum(grepl("^[[:space:]]*-?[[:space:]]*href:[[:space:]]*articles/", lines)))
  indexed <- attr(dangling, "indexed")
  expect_equal(length(indexed), .articleEntryLineCount(lines))
  # pkgdown aborts the site build when an article is missing from a non-empty
  # `articles:` index, so every built article must be listed, and listed once.
  expect_setequal(sub("^articles/", "", indexed), .builtArticles(root))
  expect_equal(anyDuplicated(indexed), 0L)
  expect_equal(
    nrow(dangling),
    0L,
    info = paste0("\n", paste(utils::capture.output(print(dangling)), collapse = "\n"))
  )
})

test_that("every model has an article, or is listed with a reason", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  built <- if (dir.exists(file.path(root, "vignettes", "articles"))) .builtArticles(root)
  # Without the source tree the article files cannot be checked, but the
  # exemption list can: an article-less model still has to be listed.
  issues <- .modelArticleIssues(nlmixr2lib::modeldb, .readArticleExemptions(), built)
  expect_equal(issues, character(), info = paste0("\n", paste(issues, collapse = "\n")))
})

test_that("the article checks detect each defect class they claim to", {
  skip_if_not_installed("yaml")
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "vignettes", "articles"), recursive = TRUE)
  writeLines("", file.path(tmp, "vignettes", "guide.Rmd"))
  writeLines("", file.path(tmp, "vignettes", "articles", "Real_2020_drug.Rmd"))
  writeLines(
    c(
      "navbar:",
      "  components:",
      "    drugs:",
      "      menu:",
      "      - text: Real",
      "        href: articles/Real_2020_drug.html",
      "      - text: Ghost",
      "        href: articles/Ghost_2018_drug.html",
      "    guide:",
      "      href: articles/guide.html",
      "    news:",
      "      href: news/index.html",
      "articles:",
      "- title: General",
      "  contents:",
      "  - guide",
      "  - missing_guide",
      "- title: internal",
      "  contents:",
      "  - articles/Real_2020_drug",
      "  - articles/Gone_2019_drug",
      "  - starts_with(\"x\")"
    ),
    file.path(tmp, "_pkgdown.yml")
  )
  dangling <- .danglingArticleLinks(tmp)
  expect_equal(attr(dangling, "nHref"), 3L)
  expect_equal(
    attr(dangling, "indexed"),
    c("guide", "missing_guide", "articles/Real_2020_drug", "articles/Gone_2019_drug")
  )
  expect_equal(.articleEntryLineCount(readLines(file.path(tmp, "_pkgdown.yml"))), 4L)
  # A section the walk cannot see is counted as empty, not as an error.
  expect_equal(.articleEntryLineCount(c("url: x", "navbar:", "  left: []")), 0L)
  expect_equal(dangling$target, c("articles/Ghost_2018_drug.html", "missing_guide", "articles/Gone_2019_drug"))
  expect_equal(dangling$where, c("navbar href", "articles group General", "articles group internal"))

  db <- data.frame(
    name = c("Real_2020_drug", "Ghost_2018_drug", "Template_pk", "Unlisted_pk", "Shared_2020_sibling"),
    vignette = c("Real_2020_drug", "Ghost_2018_drug", NA, NA, "Real_2020_drug"),
    stringsAsFactors = FALSE
  )
  exempt <- data.frame(
    model = c("Template_pk", "Real_2020_drug", "Removed_pk", "Template_pk", "Bad_status"),
    status = c("template", "template", "template", "template", "wibble"),
    reason = c("a building block", "stale", "stale", "dup", ""),
    stringsAsFactors = FALSE
  )
  issues <- .modelArticleIssues(db, exempt, built = .builtArticles(tmp))
  expect_equal(
    issues,
    c(
      "Ghost_2018_drug: vignette 'Ghost_2018_drug' is not an .Rmd under vignettes/",
      "Unlisted_pk: no article and not listed in models-without-article.csv",
      "Removed_pk: listed in models-without-article.csv but not in the modeldb",
      "Bad_status: listed in models-without-article.csv but not in the modeldb",
      "Real_2020_drug: listed in models-without-article.csv but has an article; remove the row",
      "Template_pk: listed twice in models-without-article.csv",
      "Bad_status: status 'wibble' is not one of: template, article-pending",
      "Bad_status: listed in models-without-article.csv without a reason"
    )
  )
  # A sibling sharing another model's article is legitimate, and without the
  # source tree the article files are not checked at all.
  expect_false(any(grepl("Shared_2020_sibling", issues)))
  expect_false(any(grepl("is not an .Rmd", .modelArticleIssues(db, exempt, built = NULL), fixed = TRUE)))
})
