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
      "Fixture_occ <- function() {",
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
    file.path(dir, "models", "Fixture_occ.R")
  )
  writeLines(
    c(
      "Fixture_study <- function() {",
      "  ini({",
      "    lcl <- log(5)",
      "    etastudy ~ 0.04 | SIDN",
      "  })",
      "  model({",
      "    cl <- exp(lcl + etastudy)",
      "  })",
      "}"
    ),
    file.path(dir, "models", "Fixture_study.R")
  )
  dir
}

.nestWriteArticle <- function(dir, name, code, header = "```{r sim}") {
  path <- file.path(dir, "articles", paste0(name, ".Rmd"))
  writeLines(c("# Title", "", header, code, "```"), path)
  path
}

# A model file whose ini() block is `ini_lines`.
.nestWriteModel <- function(dir, name, ini_lines) {
  path <- file.path(dir, paste0(name, ".R"))
  writeLines(
    c(
      paste0(name, " <- function() {"),
      "  ini({",
      "    lcl <- log(5)",
      paste0("    ", ini_lines),
      "  })",
      "  model({",
      "    cl <- exp(lcl)",
      "  })",
      "}"
    ),
    path
  )
  path
}

test_that("the scanner finds the `| G` etas of a model and nothing else", {
  dir <- .nestFixtureDir()
  expect_equal(
    .nestedEtas(file.path(dir, "models", "Fixture_occ.R")),
    data.frame(eta = "etaiov", group = "OCC", stringsAsFactors = FALSE)
  )
  # the `cp ~ add(addSd) | cmt` endpoint form in model() is not an eta
  expect_false("cmt" %in% .nestedEtas(file.path(dir, "models", "Fixture_occ.R"))$group)
  nested <- .nestedModels(
    file.path(dir, "models"),
    above = data.frame(model = "Fixture_study", group = "SIDN", stringsAsFactors = FALSE)
  )
  expect_equal(nested$model[order(nested$model)], c("Fixture_occ", "Fixture_study"))
  expect_equal(nested$above_subject[order(nested$model)], c(FALSE, TRUE))
})

test_that("an article that simulates a below-subject model without an assertion is flagged", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "NoAssertion",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), art)
  expect_equal(v$article, "NoAssertion.Rmd")
  expect_equal(v$models, "Fixture_occ")
  expect_equal(v$etas, "etaiov")
})

test_that("an assertion on something other than per-level distinct draws does not count", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "WrongAssertion",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
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
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(length(unique(round(sim$etaiov[sim$OCC == 1], 8))) == 3L)"
    )
  )
  summarised <- .nestWriteArticle(
    dir,
    "Summarised",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
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
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(rxode2::zeroRe(mod), events)"
    )
  )
  no_solve <- .nestWriteArticle(
    dir,
    "NoSolve",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
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
      "mod <- readModelDb(\"Fixture_study\")",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  listed <- data.frame(model = "Fixture_study", group = "SIDN", stringsAsFactors = FALSE)
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), art, above = listed)), 0L)
  unlisted <- data.frame(model = character(), group = character(), stringsAsFactors = FALSE)
  expect_equal(.nestDrawViolations(file.path(dir, "models"), art, above = unlisted)$models, "Fixture_study")
})

test_that("an eta whose `~` and `|` sit on different lines, or behind a comment, is still found", {
  dir <- withr::local_tempdir()
  split <- .nestWriteModel(dir, "Fixture_split", c("etaiov ~", "  0.25 | OCC"))
  commented <- .nestWriteModel(dir, "Fixture_commented", c("etaiov ~ # between-occasion variance", "  0.25 | OCC"))
  both <- data.frame(eta = "etaiov", group = "OCC", stringsAsFactors = FALSE)
  expect_equal(.nestedEtas(split), both)
  expect_equal(.nestedEtas(commented), both)
})

test_that("a correlated block over several etas yields one row per eta", {
  dir <- withr::local_tempdir()
  path <- .nestWriteModel(dir, "Fixture_block", "etaiov_a + etaiov_b ~ c(0.04, 0.01, 0.09) | OCC")
  expect_equal(
    .nestedEtas(path),
    data.frame(eta = c("etaiov_a", "etaiov_b"), group = "OCC", stringsAsFactors = FALSE)
  )
})

test_that("a grouping that is not a plain column name is a typed error, not a silent skip", {
  dir <- withr::local_tempdir()
  path <- .nestWriteModel(dir, "Fixture_interaction", "etaiov ~ 0.25 | interaction(a, b)")
  expect_error(.nestedEtas(path), class = "nestUnsupportedGroup")
  expect_error(.nestedEtas(path), "Fixture_interaction.R")
  expect_error(.nestedModels(dir), class = "nestUnsupportedGroup")
})

test_that("a running chunk that does not parse is reported, with its line", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "Broken",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod,"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), art)
  expect_equal(v$article, "Broken.Rmd")
  expect_equal(v$reason, "unverifiable: R chunk at line 3 does not parse")
})

test_that("chunks that do not run are skipped, however the option is written", {
  dir <- .nestFixtureDir()
  code <- c("mod <- readModelDb(\"Fixture_occ\")", "sim <- rxode2::rxSolve(mod, events")
  header <- .nestWriteArticle(dir, "HeaderOff", code, header = "```{r sim, eval = FALSE}")
  yaml <- .nestWriteArticle(dir, "YamlOff", c("#| eval: false", code))
  rstyle <- .nestWriteArticle(dir, "RStyleOff", c("#| eval = FALSE", code))
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), c(header, yaml, rstyle))), 0L)
  # the same code in a running chunk is reported
  running <- .nestWriteArticle(dir, "Running", code)
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), running)), 1L)
})

test_that("an assertion inside a chunk that does not run does not count", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "AssertOff",
    c(
      "#| eval: false",
      "stopifnot(dplyr::n_distinct(sim$etaiov[sim$OCC == 1]) == 100)"
    )
  )
  writeLines(
    c(
      readLines(art),
      "",
      "```{r run}",
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
      "```"
    ),
    art
  )
  expect_equal(.nestDrawViolations(file.path(dir, "models"), art)$article, "AssertOff.Rmd")
})

test_that("a count of subjects and a finiteness check on the eta do not make an assertion on the eta per level", {
  dir <- .nestFixtureDir()
  loose <- .nestWriteArticle(
    dir,
    "Loose",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(length(unique(sim$id)) == 100, all(is.finite(sim$etaiov)))"
    )
  )
  split_args <- .nestWriteArticle(
    dir,
    "SplitArgs",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(dplyr::n_distinct(sim$id) == 100, all(sim$OCC %in% 1:3), all(is.finite(sim$etaiov)))"
    )
  )
  per_level <- .nestWriteArticle(
    dir,
    "PerLevel",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(dplyr::n_distinct(sim$etaiov[sim$OCC == 1]) == 100)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(loose, split_args, per_level))
  expect_equal(v$article, c("Loose.Rmd", "SplitArgs.Rmd"))
  expect_equal(v$etas, c("etaiov", "etaiov"))
})

test_that("a model name built at run time that could name a below-subject model is unverifiable", {
  dir <- .nestFixtureDir()
  built <- .nestWriteArticle(
    dir,
    "Built",
    c(
      "kind <- \"occ\"",
      "mod <- readModelDb(paste0(\"Fixture_\", kind))",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  through_symbol <- .nestWriteArticle(
    dir,
    "ThroughSymbol",
    c(
      "nm <- sprintf(\"Fixture_%s\", \"occ\")",
      "mod <- nlmixr2lib::readModelDb(nm)",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  from_data <- .nestWriteArticle(
    dir,
    "FromData",
    c(
      "mod <- readModelDb(paste0(prefix, \"_\", suffix))",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(built, through_symbol, from_data))
  expect_equal(v$article, c("Built.Rmd", "ThroughSymbol.Rmd", "FromData.Rmd"))
  expect_true(all(grepl("cannot be resolved", v$reason)))
})

test_that("a built model name that cannot name a below-subject model, or is only used with zeroRe(), is not flagged", {
  dir <- .nestFixtureDir()
  unrelated <- .nestWriteArticle(
    dir,
    "Unrelated",
    c(
      "mod <- readModelDb(paste0(\"Unrelated_2020_\", drug))",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  typical <- .nestWriteArticle(
    dir,
    "BuiltTypical",
    c(
      "mod <- readModelDb(paste0(\"Fixture_\", kind))",
      "sim <- rxode2::rxSolve(rxode2::zeroRe(mod), events)"
    )
  )
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), c(unrelated, typical))), 0L)
})

test_that("a count that only builds the event table, or sits in an unrelated argument, does not cover the eta", {
  dir <- .nestFixtureDir()
  one_argument <- .nestWriteArticle(
    dir,
    "OneArgument",
    c(
      "mod <- readModelDb(\"Fixture_occ\")",
      "sim <- rxode2::rxSolve(mod, events)",
      "stopifnot(dplyr::n_distinct(sim$id) == 100 && all(is.finite(sim$etaiov)) && all(sim$OCC > 0))"
    )
  )
  through_events <- .nestWriteArticle(
    dir,
    "ThroughEvents",
    c(
      "ids <- 1:100",
      "ev <- expand.grid(id = unique(ids), OCC = 1:3)",
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), ev)",
      "stopifnot(all(is.finite(sim$etaiov[sim$OCC > 0])))"
    )
  )
  summary_other <- .nestWriteArticle(
    dir,
    "SummaryOther",
    c(
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
      "chk <- sim |>",
      "  dplyr::group_by(OCC) |>",
      "  dplyr::summarise(n_ids = dplyr::n_distinct(id), spread = sd(etaiov))",
      "stopifnot(all(chk$n_ids == 100))"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(one_argument, through_events, summary_other))
  expect_equal(v$article, c("OneArgument.Rmd", "ThroughEvents.Rmd", "SummaryOther.Rmd"))
  expect_equal(v$reason, rep("no per-level distinct-draw assertion", 3L))
})

test_that("a distinct-draw count over the eta and its group clears the article however it is written", {
  dir <- .nestFixtureDir()
  forms <- list(
    across = c(
      "per_occ <- sim |>",
      "  dplyr::group_by(OCC) |>",
      "  dplyr::summarise(dplyr::across(c(etaiov), \\(x) dplyr::n_distinct(round(x, 8)), .names = \"n_{.col}\"))",
      "stopifnot(per_occ$n_etaiov == 100)"
    ),
    by_argument = c(
      "per_occ <- dplyr::summarise(sim, n = dplyr::n_distinct(etaiov), .by = OCC)",
      "stopifnot(per_occ$n == 100)"
    ),
    magrittr = c(
      "per_occ <- sim %>% dplyr::group_by(OCC) %>% dplyr::summarise(n = dplyr::n_distinct(etaiov))",
      "stopifnot(per_occ$n == 100)"
    ),
    tapply = c(
      "per_occ <- tapply(sim$etaiov, sim$OCC, function(x) length(unique(x)))",
      "stopifnot(all(per_occ == 100))"
    ),
    expect = c(
      "testthat::expect_equal(dplyr::n_distinct(sim$etaiov[sim$OCC == 2]), 100)"
    )
  )
  arts <- vapply(
    names(forms),
    function(nm) {
      .nestWriteArticle(dir, nm, c("sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)", forms[[nm]]))
    },
    character(1L)
  )
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), arts)), 0L)
})

test_that("a model name that comes from data, a loop over the library, params or lapply() is unverifiable", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    FromLibrary = c(
      "nm <- modeldb$name[grepl(\"Fixture\", modeldb$name)]",
      "mod <- readModelDb(nm[1])",
      sim
    ),
    LoopOverLibrary = c(
      "for (m in modeldb$name) {",
      "  mod <- readModelDb(m)",
      "  sim <- rxode2::rxSolve(mod, events)",
      "}"
    ),
    FromParams = c("mod <- readModelDb(params$model)", sim),
    LapplyLoader = c(
      "nm <- modeldb$name[grepl(\"Fixture\", modeldb$name)]",
      "mods <- lapply(nm, readModelDb)",
      "sim <- rxode2::rxSolve(mods[[1]], events)"
    ),
    PurrrLoader = c(
      "nm <- readLines(\"models.txt\")",
      "mods <- purrr::map(nm, nlmixr2lib::readModelDb)",
      "sim <- rxode2::rxSolve(mods[[1]], events)"
    ),
    FunctionParameter = c(
      "load_one <- function(nm) readModelDb(nm)",
      "mods <- lapply(read.csv(\"models.csv\")$model, load_one)",
      "sim <- rxode2::rxSolve(mods[[1]], events)"
    )
  )
  arts <- vapply(
    names(cases),
    function(nm) .nestWriteArticle(dir, nm, cases[[nm]]),
    character(1L)
  )
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, paste0(names(cases), ".Rmd"))
  expect_true(all(grepl("cannot be resolved", v$reason)))
})

test_that("a name that passes through a string transformer is unverifiable whatever its fragments", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    Sub = c("mod <- readModelDb(sub(\"x$\", \"\", \"Fixture_occx\"))", sim),
    Lower = c("mod <- readModelDb(tolower(\"Fixture_OCC\"))", sim),
    ThroughSymbol = c("nm <- toupper(\"unrelated_a\")", "mod <- readModelDb(nm)", sim),
    AroundConcatenation = c("nm <- gsub(\"-\", \"_\", paste0(\"Unrelated_\", id))", "mod <- readModelDb(nm)", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, paste0(names(cases), ".Rmd"))
})

test_that("names that resolve to literals or to a pattern no below-subject model matches are not flagged", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    LapplyLiterals = c("mods <- lapply(c(\"Unrelated_a\", \"Unrelated_b\"), readModelDb)", "mod <- mods[[1]]", sim),
    LambdaLiterals = c("mods <- sapply(c(\"Unrelated_a\"), function(nm) readModelDb(nm))", "mod <- mods[[1]]", sim),
    UserFunction = c("load_one <- function(nm) readModelDb(nm)", "mod <- load_one(\"Unrelated_a\")", sim),
    UserFunctionIterated = c(
      "load_one <- function(nm) readModelDb(nm)",
      "mods <- lapply(c(\"Unrelated_a\", \"Unrelated_b\"), load_one)",
      "mod <- mods[[1]]",
      sim
    ),
    TribbleColumn = c(
      "tab <- tibble::tribble(~label, ~model, \"first\", \"Unrelated_a\", \"second\", \"Unrelated_b\")",
      "mod <- readModelDb(tab$model[1])",
      sim
    ),
    TibbleColumnWithLabelTransform = c(
      "tab <- tibble::tibble(model = c(\"Unrelated_a\", \"Unrelated_b\"), label = sub(\"^Unrelated_\", \"\", model))",
      "mod <- readModelDb(tab$model[2])",
      sim
    ),
    Branches = c("mod <- readModelDb(ifelse(dose == 1, \"Unrelated_a\", \"Unrelated_b\"))", sim),
    BuiltFromUnrelatedPrefix = c("mod <- readModelDb(sprintf(\"Unrelated_2020_%s_mouse\", drug))", sim),
    TransformedPieceOfAName = c("mod <- readModelDb(paste0(\"Unrelated_\", gsub(\"-\", \"_\", id)))", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), arts)), 0L)
})

test_that("a group that only filters the events, or only the data, does not group the count", {
  dir <- .nestFixtureDir()
  filtered_events <- .nestWriteArticle(
    dir,
    "FilteredEvents",
    c(
      "ev <- dplyr::filter(events, OCC <= 2)",
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), ev)",
      "stopifnot(dplyr::n_distinct(sim$etaiov) > 1)"
    )
  )
  filtered_data <- .nestWriteArticle(
    dir,
    "FilteredData",
    c(
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
      "chk <- sim |>",
      "  dplyr::filter(OCC <= 2) |>",
      "  dplyr::summarise(n = dplyr::n_distinct(etaiov))",
      "stopifnot(chk$n > 1)"
    )
  )
  subset_symbol <- .nestWriteArticle(
    dir,
    "SubsetSymbol",
    c(
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
      "x <- sim$etaiov[sim$OCC == 1]",
      "stopifnot(length(unique(x)) == 100)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(filtered_events, filtered_data, subset_symbol))
  expect_equal(v$article, c("FilteredEvents.Rmd", "FilteredData.Rmd", "SubsetSymbol.Rmd"))
  expect_equal(v$reason, rep("no per-level distinct-draw assertion", 3L))
})

test_that("a count of the group alone, or of something else, does not cover the eta", {
  dir <- .nestFixtureDir()
  group_only <- .nestWriteArticle(
    dir,
    "GroupOnly",
    c(
      "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
      "chk <- sim |>",
      "  dplyr::filter(is.finite(etaiov)) |>",
      "  dplyr::group_by(OCC) |>",
      "  dplyr::summarise(n = dplyr::n_distinct(OCC))",
      "stopifnot(all(chk$n == 1))"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), group_only)
  expect_equal(v$article, "GroupOnly.Rmd")
})

test_that("a loader reached through a pipe, a formula, an alias or do.call() is followed", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  unseen <- "nm <- params$model"
  cases <- list(
    NativePipe = c(unseen, "mod <- nm |> readModelDb()", sim),
    MagrittrPipe = c(unseen, "mod <- nm %>% readModelDb()", sim),
    MagrittrPipeNamespaced = c(unseen, "mod <- nm %>% nlmixr2lib::readModelDb()", sim),
    PipedIterator = c(unseen, "mods <- nm %>% purrr::map(readModelDb)", "mod <- mods[[1]]", sim),
    FormulaIterator = c(unseen, "mods <- purrr::map(nm, ~ readModelDb(.x))", "mod <- mods[[1]]", sim),
    LambdaCalledInPlace = c(unseen, "mod <- (\\(x) readModelDb(x))(nm)", sim),
    Alias = c(unseen, "rd <- readModelDb", "mod <- rd(nm)", sim),
    AliasByName = c(unseen, "rd <- get(\"readModelDb\")", "mod <- rd(nm)", sim),
    AliasIterated = c(unseen, "rd <- nlmixr2lib::readModelDb", "mods <- lapply(nm, rd)", "mod <- mods[[1]]", sim),
    DoCall = c(unseen, "mod <- do.call(readModelDb, list(nm))", sim),
    DoCallByName = c(unseen, "mod <- do.call(\"readModelDb\", list(nm))", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, paste0(names(cases), ".Rmd"))
  expect_true(all(grepl("cannot be resolved", v$reason)))
})

test_that("a loader reached the same ways on literal names is not flagged", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    NativePipe = c("mod <- \"Unrelated_a\" |> readModelDb()", sim),
    MagrittrPipe = c("nm <- \"Unrelated_a\"", "mod <- nm %>% readModelDb()", sim),
    FormulaIterator = c("mods <- purrr::map(c(\"Unrelated_a\"), ~ readModelDb(.x))", "mod <- mods[[1]]", sim),
    Alias = c("rd <- readModelDb", "mod <- rd(\"Unrelated_a\")", sim),
    DoCall = c("mod <- do.call(readModelDb, list(\"Unrelated_a\"))", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), arts)), 0L)
})

test_that("one value that cannot be seen makes a name unverifiable, whatever else is bound to it", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    Reassigned = c("nm <- \"Unrelated_a\"", "nm <- params$model", "mod <- readModelDb(nm)", sim),
    HelperLiteralAndData = c(
      "load_one <- function(n) readModelDb(n)",
      "mod <- load_one(\"Unrelated_a\")",
      "mod <- load_one(params$model)",
      sim
    ),
    IteratedLiteralAndData = c(
      "load_one <- function(n) readModelDb(n)",
      "mods <- lapply(c(\"Unrelated_a\"), load_one)",
      "mods <- lapply(params$models, load_one)",
      "mod <- mods[[1]]",
      sim
    ),
    HelperNeverCalled = c("load_one <- function(n) readModelDb(n)", "mod <- NULL", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, paste0(names(cases), ".Rmd"))
})

test_that("a parameter is looked up in its own function, not among the article's other variables", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    LoopVariableElsewhere = c(
      "fill <- function(covs) for (nm in names(covs)) ev[[nm]] <- covs[[nm]]",
      "mods <- lapply(c(\"Unrelated_a\", \"Unrelated_b\"), function(nm) readModelDb(nm))",
      "mod <- mods[[1]]",
      sim
    ),
    LocalAssignmentElsewhere = c(
      "f <- function() { nm <- params$model; nm }",
      "nm <- \"Unrelated_a\"",
      "mod <- readModelDb(nm)",
      sim
    ),
    FormulaParameter = c(
      "mods <- purrr::map(c(\"Unrelated_a\"), ~ readModelDb(.x))",
      "x <- params$value",
      "mod <- mods[[1]]",
      sim
    )
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), arts)), 0L)
})

test_that("predict() and nlmixr2(est = \"rxSolve\") count as simulating, zeroRe() does not", {
  dir <- .nestFixtureDir()
  predicted <- .nestWriteArticle(
    dir,
    "Predicted",
    c("mod <- readModelDb(\"Fixture_occ\")", "pr <- predict(mod, events)")
  )
  fitted <- .nestWriteArticle(
    dir,
    "Fitted",
    c("mod <- readModelDb(\"Fixture_occ\")", "fit <- nlmixr2::nlmixr2(mod, data, est = \"rxSolve\")")
  )
  typical <- .nestWriteArticle(
    dir,
    "TypicalPredict",
    c("mod <- readModelDb(\"Fixture_occ\")", "pr <- predict(rxode2::zeroRe(mod), events)")
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(predicted, fitted, typical))
  expect_equal(v$article, c("Predicted.Rmd", "Fitted.Rmd"))
})

test_that("an article with unresolved loader names is excused only by an evidenced entry for its model family", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "FromData",
    c("mod <- readModelDb(read.csv(\"models.csv\")$model[1])", "sim <- rxode2::rxSolve(mod, events)")
  )
  none <- .nestDynamicLoaderExceptions[0L, ]
  family <- data.frame(
    article = "FromData.Rmd",
    prefix = "Unrelated_",
    evidence = "all Unrelated_* models",
    stringsAsFactors = FALSE
  )
  overlapping <- data.frame(
    article = "FromData.Rmd",
    prefix = "Fixture_",
    evidence = "wrong family",
    stringsAsFactors = FALSE
  )
  models <- file.path(dir, "models")
  expect_equal(nrow(.nestDrawViolations(models, art, exceptions = none)), 1L)
  expect_equal(nrow(.nestDrawViolations(models, art, exceptions = family)), 0L)
  expect_equal(nrow(.nestDrawViolations(models, art, exceptions = overlapping)), 1L)
})

test_that("a filter inside the data of group_by() does not group the count, whatever the pipe", {
  dir <- .nestFixtureDir()
  head <- "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)"
  native <- .nestWriteArticle(
    dir,
    "FilterThenGroupNative",
    c(
      head,
      "chk <- sim |>",
      "  dplyr::filter(OCC == 1) |>",
      "  dplyr::group_by(id) |>",
      "  dplyr::summarise(n = dplyr::n_distinct(etaiov))",
      "stopifnot(all(chk$n == 1))"
    )
  )
  magrittr <- .nestWriteArticle(
    dir,
    "FilterThenGroupMagrittr",
    c(
      head,
      "chk <- sim %>%",
      "  dplyr::filter(OCC == 1) %>%",
      "  dplyr::group_by(id) %>%",
      "  dplyr::summarise(n = dplyr::n_distinct(etaiov))",
      "stopifnot(all(chk$n == 1))"
    )
  )
  grouped <- .nestWriteArticle(
    dir,
    "GroupThenFilterNotNeeded",
    c(
      head,
      "chk <- sim %>%",
      "  dplyr::group_by(OCC) %>%",
      "  dplyr::summarise(n = dplyr::n_distinct(etaiov))",
      "stopifnot(all(chk$n == 100))"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(native, magrittr, grouped))
  expect_equal(v$article, c("FilterThenGroupNative.Rmd", "FilterThenGroupMagrittr.Rmd"))
})

test_that("grouping does not survive extracting a column or ungroup()", {
  dir <- .nestFixtureDir()
  head <- c(
    "sim <- rxode2::rxSolve(readModelDb(\"Fixture_occ\"), events)",
    "g <- sim |> dplyr::group_by(OCC)"
  )
  column <- .nestWriteArticle(dir, "GroupedColumn", c(head, "stopifnot(dplyr::n_distinct(g$etaiov) > 1)"))
  ungrouped <- .nestWriteArticle(
    dir,
    "Ungrouped",
    c(head, "chk <- g |> dplyr::ungroup() |> dplyr::summarise(n = dplyr::n_distinct(etaiov))", "stopifnot(chk$n > 1)")
  )
  ungrouped_symbol <- .nestWriteArticle(
    dir,
    "UngroupedSymbol",
    c(
      head,
      "u <- dplyr::ungroup(g)",
      "chk <- dplyr::summarise(u, n = dplyr::n_distinct(etaiov))",
      "stopifnot(chk$n > 1)"
    )
  )
  kept <- .nestWriteArticle(
    dir,
    "GroupKept",
    c(head, "chk <- dplyr::summarise(g, n = dplyr::n_distinct(etaiov))", "stopifnot(all(chk$n == 100))")
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(column, ungrouped, ungrouped_symbol, kept))
  expect_equal(v$article, c("GroupedColumn.Rmd", "Ungrouped.Rmd", "UngroupedSymbol.Rmd"))
})

test_that("zeroRe() that keeps the omegas still simulates", {
  dir <- .nestFixtureDir()
  make <- function(name, call) {
    .nestWriteArticle(
      dir,
      name,
      c("mod <- readModelDb(\"Fixture_occ\")", paste0("sim <- rxode2::rxSolve(", call, ", events)"))
    )
  }
  arts <- c(
    SigmaOnly = make("SigmaOnly", "rxode2::zeroRe(mod, which = \"sigma\")"),
    SigmaPositional = make("SigmaPositional", "rxode2::zeroRe(mod, \"sigma\")"),
    UnseenWhich = make("UnseenWhich", "rxode2::zeroRe(mod, which = w)"),
    Default = make("Default", "rxode2::zeroRe(mod)"),
    Omega = make("Omega", "rxode2::zeroRe(mod, which = \"omega\")"),
    Both = make("Both", "rxode2::zeroRe(mod, which = c(\"omega\", \"sigma\"))")
  )
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, c("SigmaOnly.Rmd", "SigmaPositional.Rmd", "UnseenWhich.Rmd"))
})

test_that("a solver passed as a function, not called, still simulates", {
  dir <- .nestFixtureDir()
  make <- function(name, code) .nestWriteArticle(dir, name, c("mod <- readModelDb(\"Fixture_occ\")", code))
  arts <- c(
    DoCall = make("DoCall", "sim <- do.call(\"rxSolve\", list(mod, events))"),
    DoCallSymbol = make("DoCallSymbol", "sim <- do.call(rxode2::rxSolve, list(mod, events))"),
    Lapply = make("Lapply", "sims <- lapply(list(events), rxSolve, object = mod)"),
    Pmap = make("Pmap", "sims <- purrr::pmap(list(events), rxSolve, object = mod)"),
    MapMany = make("MapMany", "sims <- Map(rxSolve, list(mod), list(events))"),
    Typical = make("Typical", "sims <- lapply(list(events), rxode2::zeroRe)")
  )
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, c("DoCall.Rmd", "DoCallSymbol.Rmd", "Lapply.Rmd", "Pmap.Rmd", "MapMany.Rmd"))
})

test_that("a default of a parameter is one of its values", {
  dir <- .nestFixtureDir()
  sim <- "sim <- rxode2::rxSolve(mod, events)"
  cases <- list(
    BuiltDefault = c("load_one <- function(nm = paste0(\"Fixture_\", z)) readModelDb(nm)", "mod <- load_one()", sim),
    UnseenDefault = c(
      "load_one <- function(nm = params$model) readModelDb(nm)",
      "mod <- load_one(\"Unrelated_a\")",
      sim
    ),
    LiteralDefault = c("load_one <- function(nm = \"Unrelated_a\") readModelDb(nm)", "mod <- load_one()", sim)
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, c("BuiltDefault.Rmd", "UnseenDefault.Rmd"))
})

test_that("a below-subject model named only in a parameter default is found", {
  dir <- .nestFixtureDir()
  cases <- list(
    NameDefault = c(
      "run <- function(nm = \"Fixture_occ\") rxode2::rxSolve(readModelDb(nm), events)",
      "sim <- run()"
    ),
    LoaderDefault = c(
      "run <- function(mod = readModelDb(\"Fixture_occ\")) rxode2::rxSolve(mod, events)",
      "sim <- run()"
    ),
    LambdaDefault = c(
      "run <- \\(nm = \"Fixture_occ\") rxode2::rxSolve(readModelDb(nm), events)",
      "sim <- run()"
    )
  )
  arts <- vapply(names(cases), function(nm) .nestWriteArticle(dir, nm, cases[[nm]]), character(1L))
  v <- .nestDrawViolations(file.path(dir, "models"), arts)
  expect_equal(v$article, c("NameDefault.Rmd", "LoaderDefault.Rmd", "LambdaDefault.Rmd"))
  expect_equal(unique(v$models), "Fixture_occ")

  asserted <- .nestWriteArticle(
    dir,
    "AssertedDefault",
    c(
      cases$NameDefault,
      "stopifnot(length(unique(round(sim$etaiov[sim$OCC == 1], 8))) == 3L)"
    )
  )
  expect_equal(nrow(.nestDrawViolations(file.path(dir, "models"), asserted)), 0L)
})

test_that("a matrix or a diagonal can hold model names", {
  dir <- .nestFixtureDir()
  art <- .nestWriteArticle(
    dir,
    "MatrixNames",
    c(
      "nm <- matrix(c(\"Unrelated_a\", params$model), 1L)",
      "mod <- readModelDb(nm[1, 2])",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  expect_equal(.nestDrawViolations(file.path(dir, "models"), art)$article, "MatrixNames.Rmd")
})

test_that("an exception never excuses a name that can be a below-subject model", {
  dir <- .nestFixtureDir()
  data_only <- .nestWriteArticle(
    dir,
    "DataOnly",
    c("mod <- readModelDb(read.csv(\"models.csv\")$model[1])", "sim <- rxode2::rxSolve(mod, events)")
  )
  built <- .nestWriteArticle(
    dir,
    "AlsoBuilt",
    c(
      "mod <- readModelDb(read.csv(\"models.csv\")$model[1])",
      "other <- readModelDb(paste0(\"Fixture_\", kind))",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  hidden <- .nestWriteArticle(
    dir,
    "AlsoHidden",
    c(
      "mod <- readModelDb(read.csv(\"models.csv\")$model[1])",
      "other <- readModelDb(sub(\"x$\", \"\", \"Fixture_occ\"))",
      "sim <- rxode2::rxSolve(mod, events)"
    )
  )
  exceptions <- data.frame(
    article = c("DataOnly.Rmd", "AlsoBuilt.Rmd", "AlsoHidden.Rmd"),
    prefix = "Unrelated_",
    evidence = "all Unrelated_* models",
    stringsAsFactors = FALSE
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(data_only, built, hidden), exceptions = exceptions)
  expect_equal(unique(v$article), c("AlsoBuilt.Rmd", "AlsoHidden.Rmd"))
  expect_true(any(grepl("a name built from pieces can be Fixture_occ", v$reason)))
  expect_true(any(grepl("a string inside an unseen value is part of Fixture_occ", v$reason)))
  expect_false("DataOnly.Rmd" %in% v$article)
})

test_that("a model name held in a symbol is matched, and an assertion clears it", {
  dir <- .nestFixtureDir()
  flagged <- .nestWriteArticle(
    dir,
    "SymbolFlagged",
    c(
      "nm <- \"Fixture_occ\"",
      "sim <- rxode2::rxSolve(readModelDb(nm), events)"
    )
  )
  cleared <- .nestWriteArticle(
    dir,
    "SymbolCleared",
    c(
      "nm <- \"Fixture_occ\"",
      "sim <- rxode2::rxSolve(readModelDb(nm), events)",
      "stopifnot(dplyr::n_distinct(sim$etaiov[sim$OCC == 1]) == 100)"
    )
  )
  v <- .nestDrawViolations(file.path(dir, "models"), c(flagged, cleared))
  expect_equal(v$article, "SymbolFlagged.Rmd")
  expect_equal(v$reason, "no per-level distinct-draw assertion")
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
  nested <- .nestedModels(file.path(root, "inst", "modeldb"))
  # without the exceptions, so the table can be checked against what the gate says
  v <- .nestDrawViolations(
    file.path(root, "inst", "modeldb"),
    articles,
    exceptions = .nestDynamicLoaderExceptions[0L, ]
  )
  excepted <- v[v$article %in% .nestDynamicLoaderExceptions$article, ]
  rest <- v[!(v$article %in% .nestDynamicLoaderExceptions$article), ]
  expect_equal(
    nrow(rest),
    0L,
    info = paste(
      "Articles the gate cannot clear for a below-subject `| G` model:",
      paste0(rest$article, " (", rest$reason, "; ", rest$models, "; ", rest$etas, ")", collapse = "; "),
      "rxode2 shares one draw per level when no subject has two or more levels of G.",
      "Give one subject a second level with an evid = 2 record and assert that",
      "n_distinct(round(<eta>, 8)) per level equals the subjects on that level;",
      "see ElKhateeb_2023_tacrolimus.Rmd."
    )
  )
  # each excepted article still needs its row, for the unresolved loader names only
  expect_setequal(excepted$article, .nestDynamicLoaderExceptions$article)
  expect_true(all(grepl("cannot be resolved", excepted$reason)))
  # and no below-subject model is in an excepted family
  below <- nested$model[!nested$above_subject]
  for (prefix in .nestDynamicLoaderExceptions$prefix) {
    expect_false(any(startsWith(below, prefix)), info = prefix)
  }
  expect_true(all(nzchar(.nestDynamicLoaderExceptions$evidence)))
})
