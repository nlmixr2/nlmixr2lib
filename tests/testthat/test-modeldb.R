test_that("modeldb is correctly generated", {
  expect_true("PK_3cmt" %in% modeldb$name)
})

# Fixture model files for the build tests. The clean ones are the shipped PK
# templates, or PK_1cmt under another name; each broken one differs from a
# clean one in a single way, so it fails for exactly one reason.
.templateModelLines <- function(name) {
  readLines(system.file("modeldb", paste0(name, ".R"), package = "nlmixr2lib"))
}

# PK_1cmt with its function renamed to `newName`; `dropDescription` removes the
# required `description`, an error-severity convention issue.
.renamedTemplate <- function(newName, dropDescription = FALSE) {
  x <- sub("^PK_1cmt <- function", paste0(newName, " <- function"), .templateModelLines("PK_1cmt"))
  if (dropDescription) {
    x <- x[!grepl("^\\s*description <- ", x)]
  }
  x
}

# Five files in C-collation order: three clean, one with a convention error
# and one that cannot load (its function name does not match its file name).
.buildFixtureFiles <- function() {
  list(
    "PK_1cmt.R" = .templateModelLines("PK_1cmt"),
    "PK_2cmt.R" = .templateModelLines("PK_2cmt"),
    "sub/Fixture_2026_clean.R" = .renamedTemplate("Fixture_2026_clean"),
    "sub/Fixture_2026_mismatch.R" = .renamedTemplate("Other_2026_name"),
    "sub/Fixture_2026_nodesc.R" = .renamedTemplate("Fixture_2026_nodesc", dropDescription = TRUE)
  )
}

.writeModelFixtures <- function(dir, files) {
  for (nm in names(files)) {
    path <- file.path(dir, nm)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[nm]], path)
  }
  invisible(dir)
}

# A minimal package tree that buildModelDb() can build: it finds the root by
# its DESCRIPTION and hashes R/modeldb.R into the cache signature.
.fixturePackage <- function(root, files = .buildFixtureFiles()) {
  writeLines("Package: fixture", file.path(root, "DESCRIPTION"))
  dir.create(file.path(root, "R"))
  writeLines("# fixture", file.path(root, "R", "modeldb.R"))
  dir.create(file.path(root, "data"))
  .writeModelFixtures(file.path(root, "inst", "modeldb"), files)
  invisible(root)
}

# Spies installed with local_mocked_bindings(). Each records into an
# environment and delegates to the real function stored there by the test.
.modeldbSpy <- new.env(parent = emptyenv())

.spyFileRow <- function(dir, file) {
  .modeldbSpy$files <- c(.modeldbSpy$files, file)
  .modeldbSpy$fileRow(dir = dir, file = file)
}

.spyParseCovariateColumns <- function(path) {
  .modeldbSpy$parses <- .modeldbSpy$parses + 1L
  .modeldbSpy$parseCovariateColumns(path)
}

.spyBuildCovariateAliasMap <- function(conv) {
  .modeldbSpy$aliasBuilds <- .modeldbSpy$aliasBuilds + 1L
  .modeldbSpy$buildCovariateAliasMap(conv)
}

.spyCheckModelConventions <- function(model, verbose = TRUE) {
  .modeldbSpy$checkedNames <- c(.modeldbSpy$checkedNames, model$modelName)
  res <- .modeldbSpy$checkModelConventions(model, verbose = verbose)
  .modeldbSpy$issues <- rbind(.modeldbSpy$issues, res)
  res
}

.noCallr <- function() {
  FALSE
}

.conventionErrorState <- function() {
  env <- nlmixr2lib:::.conventionErrorEnv
  list(byModel = env$byModel, failures = env$failures)
}

test_that("a build that stops on errors caches the models that loaded cleanly", {
  skip_on_cran()
  pkg <- withr::local_tempdir("modeldb-pkg")
  .fixturePackage(pkg)
  withr::local_dir(pkg)
  withr::defer(nlmixr2lib:::.conventionErrorsReset())
  cachePath <- file.path(pkg, "data-raw", "modeldb-cache.rds")

  err <- tryCatch(suppressMessages(nlmixr2lib:::buildModelDb(workers = 1L)), error = identity)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  # Both kinds of error are reported, each naming its model or file.
  expect_match(msg, "buildModelDb(): 1 model file(s) could not be loaded.", fixed = TRUE)
  expect_match(
    msg,
    paste0(
      file.path(normalizePath(pkg, winslash = "/"), "inst/modeldb", "sub/Fixture_2026_mismatch.R"),
      ": Loading model failed due to filename/modelName mismatch"
    ),
    fixed = TRUE
  )
  expect_match(msg, "buildModelDb(): 1 convention error(s) in 1 model(s).", fixed = TRUE)
  expect_match(msg, "Fixture_2026_nodesc: description (file_metadata)", fixed = TRUE)
  # Nothing is published, but the three clean models are cached and the two
  # that failed are not.
  expect_false(file.exists(file.path(pkg, "data", "modeldb.rda")))
  expect_true(file.exists(cachePath))
  cache <- readRDS(cachePath)
  expect_identical(names(cache$entries), c("PK_1cmt.R", "PK_2cmt.R", "sub/Fixture_2026_clean.R"))

  # Fix the two broken models: the rebuild parses only those two.
  .writeModelFixtures(
    file.path(pkg, "inst", "modeldb"),
    list(
      "sub/Fixture_2026_mismatch.R" = .renamedTemplate("Fixture_2026_mismatch"),
      "sub/Fixture_2026_nodesc.R" = .renamedTemplate("Fixture_2026_nodesc")
    )
  )
  .modeldbSpy$files <- character()
  .modeldbSpy$fileRow <- nlmixr2lib:::.modeldbFileRow
  testthat::local_mocked_bindings(.modeldbFileRow = .spyFileRow, .package = "nlmixr2lib")
  doc <- suppressMessages(nlmixr2lib:::buildModelDb(workers = 1L))
  expect_type(doc, "character")
  expect_identical(.modeldbSpy$files, c("sub/Fixture_2026_mismatch.R", "sub/Fixture_2026_nodesc.R"))
  built <- readRDS(file.path(pkg, "inst", "modeldb.rds"))
  expect_identical(
    built$name,
    c("PK_1cmt", "PK_2cmt", "Fixture_2026_clean", "Fixture_2026_mismatch", "Fixture_2026_nodesc")
  )
  expect_identical(
    names(readRDS(cachePath)$entries),
    c(
      "PK_1cmt.R",
      "PK_2cmt.R",
      "sub/Fixture_2026_clean.R",
      "sub/Fixture_2026_mismatch.R",
      "sub/Fixture_2026_nodesc.R"
    )
  )
})

test_that("editing a register makes the next build re-check every model", {
  # A cached model is not re-checked, so its cached convention verdict is only
  # as current as the registers it was checked against.
  skip_on_cran()
  pkg <- withr::local_tempdir("modeldb-pkg")
  .fixturePackage(pkg, files = .buildFixtureFiles()[c("PK_1cmt.R", "PK_2cmt.R")])
  register <- file.path(pkg, "inst", "references", "fixture-register.md")
  dir.create(dirname(register))
  writeLines(c("# Fixture register", "", "### FIXTURE_COV", "- **Type:** continuous"), register)
  withr::local_dir(pkg)
  withr::defer(nlmixr2lib:::.conventionErrorsReset())
  .modeldbSpy$files <- character()
  .modeldbSpy$fileRow <- nlmixr2lib:::.modeldbFileRow
  testthat::local_mocked_bindings(.modeldbFileRow = .spyFileRow, .package = "nlmixr2lib")

  suppressMessages(nlmixr2lib:::buildModelDb(workers = 1L))
  expect_identical(.modeldbSpy$files, c("PK_1cmt.R", "PK_2cmt.R"))
  # Unchanged register: every model is a cache hit.
  suppressMessages(nlmixr2lib:::buildModelDb(workers = 1L))
  expect_identical(.modeldbSpy$files, c("PK_1cmt.R", "PK_2cmt.R"))

  cat("- **Scope:** general\n", file = register, append = TRUE)
  suppressMessages(
    expect_message(
      nlmixr2lib:::buildModelDb(workers = 1L),
      "modeldb cache invalidated (global signature changed)",
      fixed = TRUE
    )
  )
  expect_identical(.modeldbSpy$files, rep(c("PK_1cmt.R", "PK_2cmt.R"), 2L))
})

test_that("a duplicate model name stops the walk with every error gathered so far", {
  skip_on_cran()
  dir <- withr::local_tempdir("modeldb-fixture")
  files <- .buildFixtureFiles()
  files[["sub/PK_1cmt.R"]] <- .templateModelLines("PK_1cmt")
  .writeModelFixtures(dir, files)
  nlmixr2lib:::.conventionErrorsReset()
  withr::defer(nlmixr2lib:::.conventionErrorsReset())
  acc <- new.env(parent = emptyenv())

  err <- tryCatch(
    suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list(), acc = acc)),
    error = identity
  )
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(
    msg,
    paste0(
      "buildModelDb(): 1 duplicated model name(s).\n  PK_1cmt: ",
      file.path(dir, "PK_1cmt.R"),
      ", ",
      file.path(dir, "sub/PK_1cmt.R")
    ),
    fixed = TRUE
  )
  expect_match(msg, "buildModelDb(): 1 model file(s) could not be loaded.", fixed = TRUE)
  expect_match(msg, "buildModelDb(): 1 convention error(s) in 1 model(s).", fixed = TRUE)
  # What a stopped buildModelDb() saves on exit: every model that loaded
  # cleanly, both copies of the duplicate included.
  expect_identical(
    names(acc$entries),
    c("PK_1cmt.R", "PK_2cmt.R", "sub/Fixture_2026_clean.R", "sub/PK_1cmt.R")
  )
})

test_that("a build parses each register once, and not at all when every model is cached", {
  skip_on_cran()
  dir <- withr::local_tempdir("modeldb-fixture")
  files <- .buildFixtureFiles()[c("PK_1cmt.R", "PK_2cmt.R", "sub/Fixture_2026_clean.R")]
  .writeModelFixtures(dir, files)
  .modeldbSpy$parses <- 0L
  .modeldbSpy$aliasBuilds <- 0L
  .modeldbSpy$parseCovariateColumns <- nlmixr2lib:::.parseCovariateColumns
  .modeldbSpy$buildCovariateAliasMap <- nlmixr2lib:::.buildCovariateAliasMap
  testthat::local_mocked_bindings(
    .parseCovariateColumns = .spyParseCovariateColumns,
    .buildCovariateAliasMap = .spyBuildCovariateAliasMap,
    .package = "nlmixr2lib"
  )

  cold <- suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list(), reloadRegisters = TRUE))
  expect_identical(cold$modeldb$name, c("PK_1cmt", "PK_2cmt", "Fixture_2026_clean"))
  # Three models are checked, and every check asks for the alias map.
  expect_identical(.modeldbSpy$parses, 1L)
  expect_identical(.modeldbSpy$aliasBuilds, 1L)

  warm <- suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, cold$entries, reloadRegisters = TRUE))
  expect_identical(warm$modeldb, cold$modeldb)
  expect_identical(.modeldbSpy$parses, 1L)

  # A worker is seeded with the build's parse instead of parsing its own, and
  # rebuilds the alias map from it.
  nlmixr2lib:::.modeldbSeedRegisters(nlmixr2lib:::.modeldbRegisterSnapshot())
  expect_null(nlmixr2lib:::.covariateRegisterCache$aliasMap)
  nlmixr2lib:::.nlmixr2libCovariateAliasMap()
  expect_identical(.modeldbSpy$parses, 1L)
  expect_identical(.modeldbSpy$aliasBuilds, 2L)
})

test_that("the build checks a model under its own name, so its specific-scope covariates are allowed", {
  # CTROUGH is registered `Scope: specific`, and Chen_2021_lorlatinib_teae_grade3
  # is one of the models its entry lists. The model used to be checked under
  # the name rxode2 gave it, "eval(parsedFile)", which no entry lists.
  src <- system.file("modeldb", "specificDrugs", "Chen_2021_lorlatinib_teae_grade3.R", package = "nlmixr2lib")
  skip_if_not(nzchar(src), "Chen_2021_lorlatinib_teae_grade3 is not in the library")
  dir <- withr::local_tempdir("modeldb-fixture")
  lines <- readLines(src)
  writeLines(lines, file.path(dir, "Chen_2021_lorlatinib_teae_grade3.R"))
  # The same model under a name the CTROUGH entry does not list.
  writeLines(
    sub("^Chen_2021_lorlatinib_teae_grade3 <- function", "Fixture_2026_ctrough <- function", lines),
    file.path(dir, "Fixture_2026_ctrough.R")
  )
  .modeldbSpy$checkedNames <- character()
  .modeldbSpy$issues <- NULL
  .modeldbSpy$checkModelConventions <- nlmixr2lib::checkModelConventions
  testthat::local_mocked_bindings(checkModelConventions = .spyCheckModelConventions, .package = "nlmixr2lib")

  suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list()))
  expect_identical(.modeldbSpy$checkedNames, c("Chen_2021_lorlatinib_teae_grade3", "Fixture_2026_ctrough"))
  scope <- .modeldbSpy$issues[grepl("scoped 'specific'", .modeldbSpy$issues$message, fixed = TRUE), , drop = FALSE]
  expect_false("Chen_2021_lorlatinib_teae_grade3" %in% scope$model)
  expect_identical(scope$model[scope$name == "CTROUGH"], "Fixture_2026_ctrough")
})

test_that("a parallel build gives the same rows, cache and errors as a serial one", {
  skip_on_cran()
  skip_if_not_installed("callr")
  dir <- withr::local_tempdir("modeldb-fixture")
  .writeModelFixtures(dir, .buildFixtureFiles())
  withr::defer(nlmixr2lib:::.conventionErrorsReset())

  nlmixr2lib:::.conventionErrorsReset()
  serial <- suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list(), workers = 1L))
  serialErrors <- .conventionErrorState()

  # Chunks of two over five misses: three worker processes, at most two at a
  # time, so a finished worker is replaced while another is still running.
  nlmixr2lib:::.conventionErrorsReset()
  parallelRun <- suppressMessages(
    nlmixr2lib:::.addDirToModelDbCached(dir, list(), workers = 2L, chunkSize = 2L)
  )
  parallelErrors <- .conventionErrorState()

  expect_identical(parallelRun$modeldb, serial$modeldb)
  expect_identical(parallelRun$entries, serial$entries)
  expect_identical(parallelErrors, serialErrors)
  # The comparison covers clean rows, a convention error and a load failure.
  expect_identical(serial$modeldb$name, c("PK_1cmt", "PK_2cmt", "Fixture_2026_clean", "Fixture_2026_nodesc"))
  expect_identical(names(serial$entries), c("PK_1cmt.R", "PK_2cmt.R", "sub/Fixture_2026_clean.R"))
  expect_named(serialErrors$byModel, "Fixture_2026_nodesc")
  expect_named(serialErrors$failures, file.path(dir, "sub/Fixture_2026_mismatch.R"))
})

test_that("a worker that is killed fails every file of its chunk, by name", {
  # The crash path: the process dies (a segfault, the DLL limit, the OOM
  # killer) and leaves no result, only its exit.
  skip_on_cran()
  skip_if_not_installed("callr")
  dir <- withr::local_tempdir("modeldb-fixture")
  files <- c("PK_1cmt.R", "PK_2cmt.R")
  .writeModelFixtures(dir, .buildFixtureFiles()[files])
  w <- nlmixr2lib:::.modeldbStartWorker(
    loader = nlmixr2lib:::.modeldbWorkerLoader(),
    dir = dir,
    files = files,
    hashes = c("h1", "h2"),
    registers = nlmixr2lib:::.modeldbRegisterSnapshot()
  )
  w$proc$kill()
  expect_false(w$proc$is_alive())
  res <- suppressMessages(nlmixr2lib:::.modeldbCollectWorker(w))
  expect_identical(vapply(res, `[[`, "", "file"), files)
  expect_identical(vapply(res, `[[`, "", "hash"), c("h1", "h2"))
  expect_true(all(startsWith(vapply(res, `[[`, "", "failure"), "worker process failed: ")))
  expect_false(file.exists(w$log))
})

test_that("a worker that cannot load the package fails every file of its chunk, by name", {
  skip_on_cran()
  skip_if_not_installed("callr")
  w <- nlmixr2lib:::.modeldbStartWorker(
    loader = list(type = "installed", path = file.path(tempdir(), "no-such-lib", "nlmixr2lib")),
    dir = tempdir(),
    files = c("a/One_2026_x.R", "Two_2026_y.R"),
    hashes = c("h1", "h2"),
    registers = list()
  )
  w$proc$wait(timeout = 60000)
  res <- suppressMessages(nlmixr2lib:::.modeldbCollectWorker(w))
  expect_identical(vapply(res, `[[`, "", "file"), c("a/One_2026_x.R", "Two_2026_y.R"))
  expect_identical(vapply(res, `[[`, "", "hash"), c("h1", "h2"))
  expect_true(all(startsWith(vapply(res, `[[`, "", "failure"), "worker process failed: ")))
  expect_false(file.exists(w$log))
})

test_that("workers load the code this session runs, or the build stays in this session", {
  ns <- asNamespace("nlmixr2lib")
  # Code sourced into a plain environment, as roxygen2::load_source() does: no
  # other process can load it.
  sourced <- nlmixr2lib:::.modeldbWorkerLoader(codeEnv = new.env(), haveCallr = TRUE)
  expect_null(sourced$type)
  expect_match(sourced$serialReason, "not a loaded nlmixr2lib namespace", fixed = TRUE)
  # A pkgload development namespace: workers run load_all() on its tree.
  expect_identical(
    nlmixr2lib:::.modeldbWorkerLoader(codeEnv = ns, isDev = TRUE, path = "/src/nlmixr2lib", haveCallr = TRUE),
    list(type = "source", path = "/src/nlmixr2lib")
  )
  # An installed namespace: workers load the copy this session loaded, only
  # while it is still installed.
  inst <- file.path(withr::local_tempdir("lib"), "nlmixr2lib")
  dir.create(file.path(inst, "Meta"), recursive = TRUE)
  file.create(file.path(inst, "Meta", "package.rds"))
  expect_identical(
    nlmixr2lib:::.modeldbWorkerLoader(codeEnv = ns, isDev = FALSE, path = inst, haveCallr = TRUE),
    list(type = "installed", path = inst)
  )
  unlink(file.path(inst, "Meta"), recursive = TRUE)
  gone <- nlmixr2lib:::.modeldbWorkerLoader(codeEnv = ns, isDev = FALSE, path = inst, haveCallr = TRUE)
  expect_null(gone$type)
  expect_identical(gone$serialReason, paste0("the installed nlmixr2lib this session loaded is no longer at ", inst))
  expect_identical(
    nlmixr2lib:::.modeldbWorkerLoader(codeEnv = ns, haveCallr = FALSE),
    list(serialReason = "the 'callr' package is not installed")
  )
})

test_that("under a development load, workers re-run load_all() on the same tree", {
  skip_if_not_installed("callr")
  skip_if_not(nlmixr2lib:::.modeldbIsDevNamespace(), "nlmixr2lib is installed here, not a pkgload load")
  loader <- nlmixr2lib:::.modeldbWorkerLoader()
  expect_identical(loader, list(type = "source", path = getNamespaceInfo("nlmixr2lib", "path")))
  expect_true(file.exists(file.path(loader$path, "DESCRIPTION")))
})

test_that("under an installed package, workers load the copy this session loaded", {
  # The branch R CMD check takes; the parallel test above then runs its
  # workers from the installed copy.
  skip_if_not_installed("callr")
  skip_if(nlmixr2lib:::.modeldbIsDevNamespace(), "nlmixr2lib is a pkgload load here, not installed")
  loader <- nlmixr2lib:::.modeldbWorkerLoader()
  expect_identical(loader, list(type = "installed", path = getNamespaceInfo("nlmixr2lib", "path")))
  expect_true(file.exists(file.path(loader$path, "Meta", "package.rds")))
})

test_that("a build that cannot start workers parses in this session and says why", {
  skip_on_cran()
  dir <- withr::local_tempdir("modeldb-fixture")
  .writeModelFixtures(dir, .buildFixtureFiles()[c("PK_1cmt.R", "PK_2cmt.R")])
  serial <- suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list()))
  testthat::local_mocked_bindings(.modeldbHaveCallr = .noCallr, .package = "nlmixr2lib")
  suppressMessages(
    expect_message(
      fallback <- nlmixr2lib:::.addDirToModelDbCached(dir, list(), workers = 2L, chunkSize = 1L),
      "modeldb: parsing 2 model(s) in this session: the 'callr' package is not installed",
      fixed = TRUE
    )
  )
  expect_identical(fallback, serial)
})

test_that("misses are split into bounded chunks, and a few stay in this session", {
  plan <- nlmixr2lib:::.modeldbChunkPlan(20L, 8L)
  expect_identical(plan$workers, 1L)
  plan <- nlmixr2lib:::.modeldbChunkPlan(21L, 8L)
  expect_identical(plan$chunk, rep(1:2, c(20L, 1L)))
  expect_identical(plan$workers, 2L)
  plan <- nlmixr2lib:::.modeldbChunkPlan(4662L, 8L)
  expect_identical(max(plan$chunk), 24L)
  expect_identical(max(tabulate(plan$chunk)), 200L)
  expect_identical(plan$workers, 8L)
  plan <- nlmixr2lib:::.modeldbChunkPlan(5L, 2L, chunkSize = 2L)
  expect_identical(plan$chunk, c(1L, 1L, 2L, 2L, 3L))
  expect_identical(plan$workers, 2L)
  plan <- nlmixr2lib:::.modeldbChunkPlan(0L, 8L)
  expect_identical(plan$chunk, integer())
  expect_identical(plan$workers, 0L)
})

test_that("the default worker count is half of rxode2's threads, from 1 to 8", {
  expect_identical(nlmixr2lib:::.modeldbDefaultWorkers(1L), 1L)
  expect_identical(nlmixr2lib:::.modeldbDefaultWorkers(2L), 1L)
  expect_identical(nlmixr2lib:::.modeldbDefaultWorkers(5L), 2L)
  expect_identical(nlmixr2lib:::.modeldbDefaultWorkers(16L), 8L)
  expect_identical(nlmixr2lib:::.modeldbDefaultWorkers(64L), 8L)
  expect_identical(
    nlmixr2lib:::.modeldbDefaultWorkers(),
    min(8L, max(1L, as.integer(rxode2::getRxThreads()) %/% 2L))
  )
})

test_that("a worker count must be one positive whole number", {
  expect_identical(nlmixr2lib:::.modeldbCheckWorkers(3), 3L)
  expect_identical(nlmixr2lib:::.modeldbCheckWorkers(1L), 1L)
  bad <- list(0, -1L, NA, NA_integer_, 2.5, c(2L, 3L), integer(), "4", TRUE, Inf, NULL)
  for (b in bad) {
    expect_error(
      nlmixr2lib:::.modeldbCheckWorkers(b),
      "workers (or option nlmixr2lib.buildModelDb.workers)",
      fixed = TRUE
    )
  }
  # The option is validated the same way, before the build touches any file.
  withr::local_options(nlmixr2lib.buildModelDb.workers = 0)
  expect_error(nlmixr2lib:::buildModelDb(), "nlmixr2lib.buildModelDb.workers", fixed = TRUE)
})
