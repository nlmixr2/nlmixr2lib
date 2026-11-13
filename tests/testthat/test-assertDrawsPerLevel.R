# assertDrawsPerLevel() and the `nesting_assertion` convention rule.
#
# An `eta ~ variance | G` eta is drawn per subject and level of G only when some
# subject has two or more levels of G; otherwise rxode2 shares one draw per level
# without a message. The helper asserts the draws on a simulation, and the rule
# requires the article of every model with a below-subject eta to call it.

# One row per (subject, level, time), with one draw per (subject, level) pair
# unless `shared`, which gives every subject on a level the level's draw.
.drawsSim <- function(levels = list(`1` = 1:3, `2` = 4:6), shared = FALSE, times = c(1, 2)) {
  rows <- lapply(names(levels), function(lev) {
    ids <- levels[[lev]]
    draw <- if (shared) rep(as.numeric(lev) / 10, length(ids)) else as.numeric(lev) + ids / 100
    data.frame(
      id = rep(ids, each = length(times)),
      time = rep(times, length(ids)),
      OCC = as.integer(lev),
      etaiov = rep(draw, each = length(times))
    )
  })
  do.call(rbind, rows)
}

test_that("a draw per subject and level passes and the simulation comes back", {
  sim <- .drawsSim()
  expect_identical(assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC"), sim)
  expect_invisible(assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC"))
})

test_that("a draw shared by the subjects of a level is an error naming the counts", {
  sim <- .drawsSim(shared = TRUE)
  expect_error(
    assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC"),
    "`etaiov` is shared: 6 (subject, `OCC`) pairs carry 2 distinct draws; level 1: 1 distinct of 3 subject(s); level 2: 1 distinct of 3 subject(s)",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
})

test_that("a single level is checked the same way", {
  one <- list(`1` = 1:5)
  expect_silent(assertDrawsPerLevel(.drawsSim(one), eta = "etaiov", by = "OCC"))
  expect_error(
    assertDrawsPerLevel(.drawsSim(one, shared = TRUE), eta = "etaiov", by = "OCC"),
    "level 1: 1 distinct of 5 subject(s)",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
})

test_that("a subject that carries the same draw on two levels shares it", {
  sim <- .drawsSim(list(`1` = 1:3, `2` = 1:3))
  sim$etaiov <- sim$id / 10
  expect_error(
    assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC"),
    "6 (subject, `OCC`) pairs carry 3 distinct draws",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
})

test_that("every named eta is checked and the shared one is named", {
  sim <- .drawsSim()
  sim$etaother <- 1
  expect_error(
    assertDrawsPerLevel(sim, eta = c("etaiov", "etaother"), by = "OCC"),
    "`etaother` is shared",
    class = "nlmixr2lib_drawsPerLevel"
  )
  sim$etaother <- sim$etaiov * 2
  expect_silent(assertDrawsPerLevel(sim, eta = c("etaiov", "etaother"), by = "OCC"))
})

test_that("a draw that changes within a subject and level is an error", {
  sim <- .drawsSim()
  sim$etaiov[1L] <- 99
  expect_error(
    assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC"),
    "`etaiov` takes more than one value within 1 (subject, `OCC`) pair(s)",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
})

test_that("n requires exactly that many distinct draws on every level", {
  sim <- .drawsSim()
  expect_silent(assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC", n = 3L))
  expect_error(
    assertDrawsPerLevel(sim, eta = "etaiov", by = "OCC", n = 4L),
    "`etaiov` should have 4 distinct draws on each `OCC` level: level 1: 3 distinct of 3 subject(s)",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
  expect_error(
    assertDrawsPerLevel(.drawsSim(shared = TRUE), eta = "etaiov", by = "OCC", n = 3L),
    "should have 3 distinct draws",
    class = "nlmixr2lib_drawsPerLevel"
  )
})

test_that("NA, a missing column and a non-numeric eta are typed errors", {
  sim <- .drawsSim()
  na_eta <- sim
  na_eta$etaiov[2L] <- NA
  expect_error(assertDrawsPerLevel(na_eta, "etaiov", "OCC"), "1 NA value", class = "nlmixr2lib_drawsPerLevel")
  na_level <- sim
  na_level$OCC[2L] <- NA
  expect_error(assertDrawsPerLevel(na_level, "etaiov", "OCC"), "has NA values", class = "nlmixr2lib_drawsPerLevel")
  na_id <- sim
  na_id$id[2L] <- NA
  expect_error(assertDrawsPerLevel(na_id, "etaiov", "OCC"), "has NA values", class = "nlmixr2lib_drawsPerLevel")
  expect_error(
    assertDrawsPerLevel(sim, c("etaiov", "etanone"), "OCC"),
    "no column `etanone`",
    class = "nlmixr2lib_drawsPerLevel"
  )
  expect_error(assertDrawsPerLevel(sim, "etaiov", "NOLEVEL"), "no column `NOLEVEL`", class = "nlmixr2lib_drawsPerLevel")
  chr <- sim
  chr$etaiov <- as.character(chr$etaiov)
  expect_error(assertDrawsPerLevel(chr, "etaiov", "OCC"), "is not numeric", class = "nlmixr2lib_drawsPerLevel")
})

test_that("a degenerate call is an error, not a pass", {
  sim <- .drawsSim()
  expect_error(assertDrawsPerLevel(sim[0L, ], "etaiov", "OCC"), "rows")
  expect_error(assertDrawsPerLevel(sim, character(), "OCC"), "length")
  expect_error(assertDrawsPerLevel(sim, "etaiov", c("OCC", "id")), "length 1")
  expect_error(assertDrawsPerLevel(sim, c("etaiov", "etaiov"), "OCC"), "duplicated")
  expect_error(assertDrawsPerLevel(sim, "etaiov", "OCC", n = 0L), "n")
  expect_error(assertDrawsPerLevel(list(a = 1), "etaiov", "OCC"), "data.frame")
})

test_that("the subject is `id`, `sim.id`, both, or the column given", {
  sim <- .drawsSim()
  names(sim)[names(sim) == "id"] <- "sim.id"
  expect_silent(assertDrawsPerLevel(sim, "etaiov", "OCC"))
  expect_error(
    assertDrawsPerLevel(sim, "etaiov", "OCC", id = "id"),
    "no column `id`",
    class = "nlmixr2lib_drawsPerLevel"
  )
  none <- .drawsSim()
  names(none)[names(none) == "id"] <- "who"
  expect_error(
    assertDrawsPerLevel(none, "etaiov", "OCC"),
    "no `id` or `sim.id` column",
    class = "nlmixr2lib_drawsPerLevel"
  )
  expect_silent(assertDrawsPerLevel(none, "etaiov", "OCC", id = "who"))
  # Two replicates of one subject table: the same subject number is a different
  # subject in each, so its draws may differ between them.
  two <- rbind(transform(.drawsSim(), sim.id = 1L), transform(.drawsSim(), sim.id = 2L, etaiov = etaiov + 0.5))
  expect_silent(assertDrawsPerLevel(two, "etaiov", "OCC"))
})

test_that("rxode2 shares the draw when every subject has one level, and the helper says so", {
  skip_on_cran()
  fun <- function() {
    ini({
      lcl <- log(5)
      etacl ~ 0.09
      etaiov ~ 0.25 | OCC
    })
    model({
      cl <- exp(lcl + etacl + etaiov)
      cp <- 100 / cl
    })
  }
  mod <- suppressWarnings(suppressMessages(rxode2::rxode2(fun)))
  events <- data.frame(
    id = rep(1:6, each = 2L),
    time = rep(c(1, 2), 6L),
    evid = 0L,
    OCC = rep(1:2, each = 6L)
  )
  rxode2::rxSetSeed(20260101)
  shared <- as.data.frame(rxode2::rxSolve(mod, events))
  expect_error(
    assertDrawsPerLevel(shared, "etaiov", "OCC"),
    "level 1: 1 distinct of 3 subject(s); level 2: 1 distinct of 3 subject(s)",
    fixed = TRUE,
    class = "nlmixr2lib_drawsPerLevel"
  )
  # One record that gives subject 1 a third level makes rxode2 draw per subject
  # and level.
  extra <- rbind(events, data.frame(id = 1L, time = 3, evid = 2L, OCC = 3L))
  extra <- extra[order(extra$id, extra$time), ]
  per_level <- as.data.frame(rxode2::rxSolve(mod, extra))
  expect_silent(assertDrawsPerLevel(per_level, "etaiov", "OCC"))
})

# --- the rule -----------------------------------------------------------------

# A model function whose ini() has `eta ~ 0.25 | group` and the `extra` ini()
# lines, as an rxUi; every eta they declare is used in model(). `vignette` NULL
# leaves the field out.
.nestFixtureUi <- function(group = "OCC", eta = "etaiov", vignette = "Fixture_article", extra = character()) {
  extra_etas <- trimws(unlist(strsplit(sub("~.*$", "", extra[grepl("~", extra, fixed = TRUE)]), "+", fixed = TRUE)))
  text <- c(
    "function() {",
    if (!is.null(vignette)) sprintf("  vignette <- \"%s\"", vignette),
    "  ini({",
    "    lcl <- log(5)",
    "    etacl ~ 0.09",
    sprintf("    %s ~ 0.25 | %s", eta, group),
    paste0("    ", extra),
    "    addSd <- 0.1",
    "  })",
    "  model({",
    sprintf("    cl <- exp(%s)", paste(c("lcl", "etacl", eta, extra_etas), collapse = " + ")),
    "    cp <- 100 / cl",
    "    cp ~ add(addSd)",
    "  })",
    "}"
  )
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr(eval(parse(text = paste(text, collapse = "\n"))))))
}

# An article directory with `Fixture_article.Rmd` holding `code` in one chunk.
.nestFixtureArticles <- function(code, header = "```{r sim}", name = "Fixture_article", env = parent.frame()) {
  dir <- file.path(withr::local_tempdir(.local_envir = env), "articles")
  dir.create(dir)
  writeLines(c("# Title", "", header, code, "```"), file.path(dir, paste0(name, ".Rmd")))
  dir
}

.nestRuleIssues <- function(ui, dir, model = "Fixture_occ") {
  .checkNestingAssertion(ui, model, articleDir = dir)
}

.callOk <- "nlmixr2lib:::assertDrawsPerLevel(sim, eta = \"etaiov\", by = \"OCC\")"

test_that("a nested eta with the helper called on its eta and group in the article is clean", {
  skip_if_not_installed("nlmixr2est")
  ui <- .nestFixtureUi()
  cases <- list(
    full = .callOk,
    bare = "assertDrawsPerLevel(sim, eta = \"etaiov\", by = \"OCC\")",
    positional = "assertDrawsPerLevel(sim, \"etaiov\", \"OCC\")",
    reordered = "assertDrawsPerLevel(by = \"OCC\", sim = sim, eta = \"etaiov\", n = 3L)",
    exported = "nlmixr2lib::assertDrawsPerLevel(sim, eta = \"etaiov\", by = \"OCC\")",
    piped = "sim |> nlmixr2lib:::assertDrawsPerLevel(eta = \"etaiov\", by = \"OCC\")",
    several = "assertDrawsPerLevel(sim, eta = c(\"other\", \"etaiov\"), by = \"OCC\")"
  )
  for (nm in names(cases)) {
    dir <- .nestFixtureArticles(c("sim <- run()", cases[[nm]]))
    expect_equal(nrow(.nestRuleIssues(ui, dir)), 0L, label = nm)
  }
})

test_that("a nested eta without a matching call is an error naming the eta and the group", {
  skip_if_not_installed("nlmixr2est")
  ui <- .nestFixtureUi()
  none <- .nestRuleIssues(ui, .nestFixtureArticles("sim <- run()"))
  expect_equal(nrow(none), 1L)
  expect_identical(none$category, "nesting_assertion")
  expect_identical(none$severity, "error")
  expect_identical(none$name, "etaiov")
  expect_match(
    none$message,
    "`etaiov ~ ... | OCC` is drawn per level of OCC: the article has no assertDrawsPerLevel() call",
    fixed = TRUE
  )
  expect_match(none$suggestion, "assertDrawsPerLevel(sim, eta = ", fixed = TRUE)

  wrong <- list(
    group = "assertDrawsPerLevel(sim, eta = \"etaiov\", by = \"STUDY\")",
    eta = "assertDrawsPerLevel(sim, eta = \"etacl\", by = \"OCC\")",
    twogroups = "assertDrawsPerLevel(sim, eta = \"etaiov\", by = c(\"OCC\", \"id\"))"
  )
  for (nm in names(wrong)) {
    v <- .nestRuleIssues(ui, .nestFixtureArticles(wrong[[nm]]))
    expect_equal(nrow(v), 1L, label = nm)
    expect_match(
      v$message,
      "no assertDrawsPerLevel() call in the article has `eta` containing \"etaiov\" and `by = \"OCC\"`",
      fixed = TRUE
    )
  }
  not_literal <- .nestRuleIssues(
    ui,
    .nestFixtureArticles(c("nm <- \"etaiov\"", "assertDrawsPerLevel(sim, eta = nm, by = \"OCC\")"))
  )
  expect_match(not_literal$message, "1 call(s) have an `eta` or `by` that is not a string literal", fixed = TRUE)
})

test_that("a call in a chunk that does not run does not count, however eval is written", {
  skip_if_not_installed("nlmixr2est")
  ui <- .nestFixtureUi()
  headers <- c("```{r sim, eval = FALSE}", "```{r sim, eval=F}")
  for (h in headers) {
    expect_equal(nrow(.nestRuleIssues(ui, .nestFixtureArticles(.callOk, header = h))), 1L, label = h)
  }
  for (opt in c("#| eval: false", "#| eval = FALSE")) {
    expect_equal(nrow(.nestRuleIssues(ui, .nestFixtureArticles(c(opt, .callOk)))), 1L, label = opt)
  }
  expect_equal(nrow(.nestRuleIssues(ui, .nestFixtureArticles(.callOk, header = "```{r sim, eval = TRUE}"))), 0L)
})

test_that("an article that cannot be read, or has a chunk that does not parse, is an error", {
  skip_if_not_installed("nlmixr2est")
  ui <- .nestFixtureUi()
  empty <- .nestFixtureArticles("1")
  expect_match(.nestRuleIssues(.nestFixtureUi(vignette = NULL), empty)$message, "has no `vignette` field", fixed = TRUE)
  expect_match(
    .nestRuleIssues(.nestFixtureUi(vignette = "Missing_article"), empty)$message,
    "the article `Missing_article.Rmd` named by the model's `vignette` field is not in vignettes/articles or vignettes",
    fixed = TRUE
  )
  broken <- .nestFixtureArticles(c(.callOk, "x <- (1 +"))
  v <- .nestRuleIssues(ui, broken)
  expect_equal(nrow(v), 1L)
  expect_match(
    v$message,
    "the article has no assertDrawsPerLevel() call; a running chunk at line 3 does not parse",
    fixed = TRUE
  )
  # The call is in a chunk that parses and a later chunk does not.
  later <- .nestFixtureArticles(c(.callOk, "```", "", "```{r later}", "x <- (1 +"))
  v <- .nestRuleIssues(ui, later)
  expect_equal(nrow(v), 1L)
  expect_match(
    v$message,
    "a running chunk at line 7 does not parse, so the article cannot be read in full",
    fixed = TRUE
  )
})

test_that("an article in vignettes/ beside articles/ is found", {
  skip_if_not_installed("nlmixr2est")
  articles <- .nestFixtureArticles("1")
  writeLines(c("```{r}", .callOk, "```"), file.path(dirname(articles), "Fixture_top.Rmd"))
  expect_equal(nrow(.nestRuleIssues(.nestFixtureUi(vignette = "Fixture_top"), articles)), 0L)
})

test_that("the rule does nothing without an articles directory, a nested eta, or for an evidenced above-subject group", {
  skip_if_not_installed("nlmixr2est")
  expect_equal(nrow(.nestRuleIssues(.nestFixtureUi(), NULL)), 0L)
  plain <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr(function() {
    ini({
      lcl <- log(5)
      etacl ~ 0.09
      addSd <- 0.1
    })
    model({
      cl <- exp(lcl + etacl)
      cp <- 100 / cl
      cp ~ add(addSd)
    })
  })))
  expect_equal(nrow(.nestRuleIssues(plain, .nestFixtureArticles("1"))), 0L)
  study <- .nestFixtureUi(group = "SIDN", eta = "etalcl_study", vignette = NULL)
  expect_equal(nrow(.nestRuleIssues(study, .nestFixtureArticles("1"), model = "Qi_2024_vosoritide")), 0L)
  expect_equal(nrow(.nestRuleIssues(study, .nestFixtureArticles("1"), model = "Other_2030_model")), 1L)
  # The allowlist is per group: a below-subject eta of an allowlisted model is checked.
  mixed <- .nestFixtureUi(extra = "etalcl_study ~ 0.04 | SIDN")
  v <- .nestRuleIssues(mixed, .nestFixtureArticles("1"), model = "Qi_2024_vosoritide")
  expect_identical(v$name, "etaiov")
})

test_that("each uncovered eta of a model is reported, and a correlated block yields both etas", {
  skip_if_not_installed("nlmixr2est")
  ui <- .nestFixtureUi(extra = c("etaiov2 + etaiov3 ~ c(0.1, 0.01, 0.1) | OCC", "etaother ~ 0.2 | MD_CATHETER"))
  expect_identical(
    .nestedEtaRows(ui),
    data.frame(
      eta = c("etaiov", "etaiov2", "etaiov3", "etaother"),
      group = c("OCC", "OCC", "OCC", "MD_CATHETER"),
      stringsAsFactors = FALSE
    )
  )
  code <- c(
    "assertDrawsPerLevel(sim, eta = c(\"etaiov\", \"etaiov3\"), by = \"OCC\")",
    "assertDrawsPerLevel(sim, eta = \"etaother\", by = \"MD_CATHETER\")"
  )
  v <- .nestRuleIssues(ui, .nestFixtureArticles(code))
  expect_identical(v$name, "etaiov2")
})

test_that("checkModelConventions() runs the rule on a model with a nested eta", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2est")
  skip_if(is.null(.nestingArticleDir()), "no vignettes/articles directory (installed package)")
  ui <- .nestFixtureUi(vignette = "No_such_article_in_the_source_tree")
  res <- suppressWarnings(suppressMessages(.checkOneModel(ui, verbose = FALSE)))
  hit <- res[res$category == "nesting_assertion", ]
  expect_equal(nrow(hit), 1L)
  expect_identical(hit$severity, "error")
})

test_that("the above-subject allowlist names only nested groups that still exist, with evidence", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2est")
  above <- .nestingAboveSubjectGroups
  expect_true(all(nzchar(above$evidence)))
  for (i in seq_len(nrow(above))) {
    ui <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr(readModelDb(above$model[[i]]))))
    expect_true(above$group[[i]] %in% .nestedEtaRows(ui)$group, label = above$model[[i]])
  }
})

# The models of the library that have a below-subject `| G` eta, with their
# articles' calls, and the proof that the rule goes red when the call is removed.
.nestedLibraryModels <- c(
  "Bindellini_2025_cefazolin",
  "ElKhateeb_2023_tacrolimus",
  "ElKhateeb_2023_tacrolimus_ckdstage",
  "ElKhateeb_2023_tacrolimus_egfr",
  "Lindauer_2025_risperidone_oral"
)

test_that("the library's below-subject models have their draws asserted, and lose that when the call goes", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2est")
  articles <- .nestingArticleDir()
  skip_if(is.null(articles), "no vignettes/articles directory (installed package)")
  scratch <- file.path(withr::local_tempdir(), "articles")
  dir.create(scratch)
  for (m in .nestedLibraryModels) {
    ui <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr(readModelDb(m))))
    expect_gt(nrow(.nestedEtaRows(ui)), 0L, label = m)
    expect_equal(nrow(.checkNestingAssertion(ui, m)), 0L, label = m)
    # The same article with the call disarmed: assertDrawsPerLevel( -> invisible(.
    path <- file.path(articles, paste0(ui$meta$vignette, ".Rmd"))
    txt <- readLines(path, warn = FALSE)
    expect_true(any(grepl("assertDrawsPerLevel(", txt, fixed = TRUE)), label = m)
    writeLines(gsub("assertDrawsPerLevel(", "invisible(", txt, fixed = TRUE), file.path(scratch, basename(path)))
    red <- .checkNestingAssertion(ui, m, articleDir = scratch)
    expect_equal(sort(red$name), sort(.nestedEtaRows(ui)$eta), label = m)
    expect_true(all(red$severity == "error"), label = m)
  }
})

test_that("a model with a nested eta is not kept in the build cache, one without is", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2est")
  models <- system.file("modeldb", package = "nlmixr2lib")
  nested <- list.files(models, "^Lindauer_2025_risperidone_oral\\.R$", recursive = TRUE, full.names = TRUE)
  plain <- list.files(models, "^PK_1cmt\\.R$", recursive = TRUE, full.names = TRUE)
  skip_if(length(nested) != 1L || length(plain) != 1L, "fixture models not found")
  fixtureDir <- withr::local_tempdir()
  file.copy(c(nested, plain), fixtureDir)
  res <- suppressMessages(suppressWarnings(nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = list())))
  expect_setequal(res$modeldb$name, c("Lindauer_2025_risperidone_oral", "PK_1cmt"))
  expect_named(res$entries, "PK_1cmt.R")
})
