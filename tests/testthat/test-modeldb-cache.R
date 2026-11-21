test_that(".modeldbCacheRead returns NULL for missing and corrupt files", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  expect_null(nlmixr2lib:::.modeldbCacheRead(tmp))

  writeBin(as.raw(c(0x00, 0x01, 0x02, 0x03, 0x04)), tmp)
  expect_null(nlmixr2lib:::.modeldbCacheRead(tmp))

  saveRDS(list(foo = 1), file = tmp)
  expect_null(nlmixr2lib:::.modeldbCacheRead(tmp))
})

test_that(".modeldbCacheWrite round-trips through .modeldbCacheRead", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  payload <- list(
    globalSig = list(a = 1L, b = "x"),
    entries = list(
      "foo.R" = list(hash = "abc", row = data.frame(name = "foo"))
    )
  )
  nlmixr2lib:::.modeldbCacheWrite(tmp, payload)
  expect_identical(nlmixr2lib:::.modeldbCacheRead(tmp), payload)
})

test_that(".modeldbCacheWrite warns instead of erroring on a write failure", {
  # Pointing at a path whose parent is an existing regular file forces
  # dir.create() to fail silently, then saveRDS() to error on open.
  blocker <- tempfile()
  file.create(blocker)
  on.exit(unlink(blocker), add = TRUE)
  badPath <- file.path(blocker, "cache.rds")
  expect_warning(
    nlmixr2lib:::.modeldbCacheWrite(badPath, list(globalSig = list(), entries = list())),
    "Could not write modeldb cache"
  )
  expect_false(file.exists(badPath))
})

test_that(".modeldbCacheWrite leaves only the cache behind, whether or not it succeeds", {
  dir <- withr::local_tempdir("modeldb-cache")
  path <- file.path(dir, "modeldb-cache.rds")
  payload <- list(globalSig = list(a = 1L), entries = list())
  nlmixr2lib:::.modeldbCacheWrite(path, payload)
  nlmixr2lib:::.modeldbCacheWrite(path, payload)
  expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), "modeldb-cache.rds")

  # The rename fails when the target is a non-empty directory; the temporary
  # file written beside it must not be left over.
  unlink(path)
  dir.create(file.path(path, "occupied"), recursive = TRUE)
  expect_warning(nlmixr2lib:::.modeldbCacheWrite(path, payload), "Could not write modeldb cache")
  expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), "modeldb-cache.rds")
  expect_true(dir.exists(path))
})

test_that(".modeldbGlobalSig has expected structure and reacts to build, checker and register edits", {
  pkgDir <- tempfile("pkg")
  dir.create(file.path(pkgDir, "R"), recursive = TRUE)
  on.exit(unlink(pkgDir, recursive = TRUE), add = TRUE)
  writeLines("fake", file.path(pkgDir, "R/modeldb.R"))

  writeLines("checker", file.path(pkgDir, "R/checkModelConventions.R"))
  writeLines("conventions", file.path(pkgDir, "R/conventions.R"))
  dir.create(file.path(pkgDir, "inst/references"), recursive = TRUE)
  writeLines("b", file.path(pkgDir, "inst/references/b-register.md"))
  writeLines("a", file.path(pkgDir, "inst/references/a-register.md"))

  sig <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
  expect_named(
    sig,
    c(
      "cacheVersion",
      "rVersion",
      "nlmixr2",
      "nlmixr2est",
      "rxode2",
      "modeldbRhash",
      "checkerHash",
      "registerHash"
    ),
    ignore.order = TRUE
  )
  expect_identical(sig$cacheVersion, 3L)
  expect_identical(sig$rVersion, R.version.string)
  expect_equal(nchar(sig$modeldbRhash), 32L)
  expect_named(sig$checkerHash, c("checkModelConventions.R", "conventions.R"))
  expect_identical(nchar(unname(sig$checkerHash)), c(32L, 32L))
  # Registers in C-collation order, whatever order they were written in.
  expect_named(sig$registerHash, c("a-register.md", "b-register.md"))
  expect_identical(nchar(unname(sig$registerHash)), c(32L, 32L))
  expect_identical(nlmixr2lib:::.modeldbGlobalSig(pkgDir), sig)

  writeLines("different", file.path(pkgDir, "R/modeldb.R"))
  sig2 <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
  expect_false(identical(sig$modeldbRhash, sig2$modeldbRhash))

  # Each checker file and each register is part of the signature.
  for (f in c("R/checkModelConventions.R", "R/conventions.R", "inst/references/b-register.md")) {
    before <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
    cat("edit\n", file = file.path(pkgDir, f), append = TRUE)
    expect_false(identical(nlmixr2lib:::.modeldbGlobalSig(pkgDir), before), label = f)
  }
  # A register added or removed changes it too.
  before <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
  writeLines("c", file.path(pkgDir, "inst/references/c-register.md"))
  expect_named(
    nlmixr2lib:::.modeldbGlobalSig(pkgDir)$registerHash,
    c("a-register.md", "b-register.md", "c-register.md")
  )
  unlink(file.path(pkgDir, "inst/references/c-register.md"))
  expect_identical(nlmixr2lib:::.modeldbGlobalSig(pkgDir), before)
})

# Counts calls to the per-file loader, delegating to the real one.
.cacheSpy <- new.env(parent = emptyenv())

.countingFileRow <- function(dir, file) {
  .cacheSpy$count <- .cacheSpy$count + 1L
  .cacheSpy$fileRow(dir = dir, file = file)
}

test_that(".addDirToModelDbCached parses on cold run and reuses cache on warm run", {
  skip_if_not_installed("nlmixr2est")

  fixtureDir <- tempfile("modeldb-fixture")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE), add = TRUE)

  srcFile <- system.file("modeldb/PK_1cmt.R", package = "nlmixr2lib")
  if (!nzchar(srcFile)) {
    skip("PK_1cmt fixture not installed")
  }
  file.copy(srcFile, file.path(fixtureDir, "PK_1cmt.R"))

  # Spy on the per-file loader so we can assert cache hits avoid re-parsing.
  .cacheSpy$count <- 0L
  .cacheSpy$fileRow <- nlmixr2lib:::.modeldbFileRow
  testthat::local_mocked_bindings(.modeldbFileRow = .countingFileRow, .package = "nlmixr2lib")

  cold <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = list())
  expect_equal(.cacheSpy$count, 1L)
  expect_equal(nrow(cold$modeldb), 1L)
  expect_identical(cold$modeldb$name, "PK_1cmt")
  expect_named(cold$entries, "PK_1cmt.R")
  expect_equal(nchar(cold$entries[["PK_1cmt.R"]]$hash), 32L)

  warm <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = cold$entries)
  expect_equal(.cacheSpy$count, 1L) # unchanged — no re-parse on cache hit
  expect_equal(warm$modeldb, cold$modeldb)
  expect_identical(warm$entries, cold$entries)
})

test_that(".addDirToModelDbCached re-parses a file whose content changes", {
  skip_if_not_installed("nlmixr2est")

  fixtureDir <- tempfile("modeldb-fixture")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE), add = TRUE)

  srcFile <- system.file("modeldb/PK_1cmt.R", package = "nlmixr2lib")
  if (!nzchar(srcFile)) {
    skip("PK_1cmt fixture not installed")
  }
  fixturePath <- file.path(fixtureDir, "PK_1cmt.R")
  file.copy(srcFile, fixturePath)

  cold <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = list())
  oldHash <- cold$entries[["PK_1cmt.R"]]$hash

  cat("\n", file = fixturePath, append = TRUE)
  warm <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = cold$entries)
  expect_false(identical(warm$entries[["PK_1cmt.R"]]$hash, oldHash))
  expect_identical(warm$modeldb$name, cold$modeldb$name)
})

test_that(".addDirToModelDbCached drops stale entries for files not on disk", {
  skip_if_not_installed("nlmixr2est")

  fixtureDir <- tempfile("modeldb-fixture")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE), add = TRUE)

  srcFile <- system.file("modeldb/PK_1cmt.R", package = "nlmixr2lib")
  if (!nzchar(srcFile)) {
    skip("PK_1cmt fixture not installed")
  }
  file.copy(srcFile, file.path(fixtureDir, "PK_1cmt.R"))

  staleEntries <- list(
    "ghost.R" = list(hash = "deadbeef", row = data.frame(name = "ghost"))
  )
  result <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = staleEntries)
  expect_false("ghost.R" %in% names(result$entries))
  expect_true("PK_1cmt.R" %in% names(result$entries))
})
