# Gate: an article that simulates a model with a below-subject `eta ~ var | G`
# eta must assert that the eta is drawn per subject and level. See
# helper-articleNesting.R for why a silent share is possible.

.nestRoot <- function() {
  testthat::test_path("..", "..")
}

.nestSkipWithoutSource <- function() {
  testthat::skip_on_cran()
  root <- .nestRoot()
  testthat::skip_if(
    !dir.exists(file.path(root, "inst", "modeldb")) || !dir.exists(file.path(root, "vignettes")),
    "package source tree not present (installed package)"
  )
  root
}

# A model file with one `| G` eta, and the articles that exercise the gate.
.nestFixtureDir <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  dir.create(file.path(dir, "models"))
  dir.create(file.path(dir, "articles"))
  writeLines(
    c(
      "Fix_occ <- function() {",
      "  ini({",
      "    lcl <- log(5)",
      "    etalcl ~ 0.09",
      "    etaiov ~ 0.25 | OCC # per-occasion draw",
      "  })",
      "  model({",
      "    cl <- exp(lcl + etalcl + etaiov)",
      "    cp ~ add(addSd) | cmt",
      "  })",
      "}"
    ),
    file.path(dir, "models", "Fix_occ.R")
  )
  writeLines(
    c(
      "Fix_study <- function() {",
      "  ini({",
      "    lcl <- log(5)",
      "    etastudy ~ 0.04 | SIDN",
      "  })",
      "  model({",
      "    cl <- exp(lcl + etastudy)",
      "  })",
      "}"
    ),
    file.path(dir, "models", "Fix_study.R")
  )
  dir
}

.nestWriteArticle <- function(dir, name, code) {
  path <- file.path(dir, "articles", paste0(name, ".Rmd"))
  writeLines(c("# Title", "", "```{r sim}", code, "```"), path)
  path
}

test_that("the scanner finds the `| G` etas of a model and nothing else", {
  dir <- .nestFixtureDir()
  expect_equal(
    .nestedEtas(file.path(dir, "models", "Fix_occ.R")),
    data.frame(eta = "etaiov", group = "OCC", stringsAsFactors = FALSE)
  )
  # the `cp ~ add(addSd) | cmt` endpoint form in model() is not an eta
  expect_false("cmt" %in% .nestedEtas(file.path(dir, "models", "Fix_occ.R"))$group)
  nested <- .nestedModels(
    file.path(dir, "models"),
    above = data.frame(model = "Fix_study", group = "SIDN", stringsAsFactors = FALSE)
  )
  expect_equal(nested$model[order(nested$model)], c("Fix_occ", "Fix_study"))
  expect_equal(nested$above_subject[order(nested$model)], c(FALSE, TRUE))
})

test_that("an article that simulates a below-subject model without an assertion is flagged", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "NoAssertion",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), art)
  expect_equal(v$article, "NoAssertion.Rmd")
  expect_equal(v$models, "Fix_occ")
  expect_equal(v$etas, "etaiov")
})

test_that("an assertion on something other than per-level distinct draws does not count", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "WrongAssertion",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(nrow(sim) == 3)",
      "# stopifnot(dplyr::n_distinct(sim$etaiov) == 3)"
    )
  )
  expect_equal(.nestDrawViolations(file.path(dir, "models"), art)$article, "WrongAssertion.Rmd")
})

test_that("a distinct-draw assertion, directly or through a summary, clears the article", {
  dir <- .nestFixtureDir()
  direct <- .nestWriteArticle(
    dir,
    "Direct",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(length(unique(round(sim$etaiov, 8))) == 3L)"
    )
  )
  summarised <- .nestWriteArticle(
    dir,
    "Summarised",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "per_occ <- sim |>",
      "  dplyr::group_by(OCC) |>",
      "  dplyr::summarise(n_distinct_iov = dplyr::n_distinct(round(etaiov, 8)))",
      "stopifnot(all(per_occ$n_distinct_iov == 3L))"
    )
  )
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), c(direct, summarised))), 0L)
})

test_that("an article that only solves with zeroRe(), or never solves, is not flagged", {
  dir <- .nestFixtureDir()
  typical <- .nestWriteArticle(
    dir,
    "Typical",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "sim <- rxode2::rxSolve(rxode2::zeroRe(mod), events)"
    )
  )
  no_solve <- .nestWriteArticle(
    dir,
    "NoSolve",
    c(
      "mod <- readModelDb(\"Fix_occ\")",
      "mod$theta"
    )
  )
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), c(typical, no_solve))), 0L)
})

test_that("an above-subject grouping on the allowlist is not flagged, and unlisted ones are", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "Study",
    c(
      "mod <- readModelDb(\"Fix_study\")",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  listed <- data.frame(model = "Fix_study", group = "SIDN", stringsAsFactors = FALSE)
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), art, above = listed)), 0L)
  unlisted <- data.frame(model = character(), group = character(), stringsAsFactors = FALSE)
  expect_equal(.nestDrawViolations(file.path(dir, "models"), art, above = unlisted)$models, "Fix_study")
})

test_that("the above-subject allowlist names only `| G` etas that still exist, with evidence", {
  root <- .nestSkipWithoutSource()
  nested <- .nestedModels(file.path(root, "inst", "modeldb"))
  listed <- .nestAboveSubjectGroups
  expect_false(anyDuplicated(paste(listed$model, listed$group)) > 0L)
  expect_true(all(nzchar(listed$evidence)))
  for (i in seq_len(nrow(listed))) {
    hit <- nested[nested$model == listed$model[i] & nested$group == listed$group[i], ]
    expect_gt(nrow(hit), 0L)
  }
})

test_that("every article that simulates a below-subject `| G` model asserts per-level distinct draws", {
  root <- .nestSkipWithoutSource()
  articles <- list.files(
    file.path(root, "vignettes"),
    pattern = "[.]Rmd$",
    recursive = TRUE,
    full.names = TRUE
  )
  expect_gt(length(articles), 100L)
  v <- .nestDrawViolations(file.path(root, "inst", "modeldb"), articles)
  expect_equal(
    nrow(v),
    0L,
    info = paste(
      "Articles simulating a below-subject `| G` model without asserting that the eta is drawn per subject and level:",
      paste0(v$article, " (", v$models, "; ", v$etas, ")", collapse = "; "),
      "rxode2 shares one draw per level when no subject has two or more levels of G.",
      "Give one subject a second level with an evid = 2 record and assert that",
      "n_distinct(round(<eta>, 8)) per level equals the subjects on that level;",
      "see ElKhateeb_2023_tacrolimus.Rmd."
    )
  )
})
