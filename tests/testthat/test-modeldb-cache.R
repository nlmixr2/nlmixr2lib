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

test_that(".modeldbGlobalSig has expected structure and reacts to R/modeldb.R edits", {
  pkgDir <- tempfile("pkg")
  dir.create(file.path(pkgDir, "R"), recursive = TRUE)
  on.exit(unlink(pkgDir, recursive = TRUE), add = TRUE)
  writeLines("fake", file.path(pkgDir, "R/modeldb.R"))

  sig <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
  expect_named(
    sig,
    c("cacheVersion", "rVersion", "nlmixr2", "nlmixr2est", "rxode2", "modeldbRhash"),
    ignore.order = TRUE
  )
  expect_identical(sig$cacheVersion, 1L)
  expect_identical(sig$rVersion, R.version.string)
  expect_equal(nchar(sig$modeldbRhash), 32L)

  writeLines("different", file.path(pkgDir, "R/modeldb.R"))
  sig2 <- nlmixr2lib:::.modeldbGlobalSig(pkgDir)
  expect_false(identical(sig$modeldbRhash, sig2$modeldbRhash))
})

test_that(".addDirToModelDbCached parses on cold run and reuses cache on warm run", {
  skip_if_not_installed("nlmixr2est")

  fixtureDir <- tempfile("modeldb-fixture")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE), add = TRUE)

  srcFile <- system.file("modeldb/PK_1cmt.R", package = "nlmixr2lib")
  if (!nzchar(srcFile)) skip("PK_1cmt fixture not installed")
  file.copy(srcFile, file.path(fixtureDir, "PK_1cmt.R"))

  # Spy on addFileToModelDb so we can assert cache hits avoid re-parsing.
  ns <- asNamespace("nlmixr2lib")
  orig <- get("addFileToModelDb", envir = ns)
  count <- 0L
  replacement <- function(...) {
    count <<- count + 1L
    orig(...)
  }
  unlockBinding("addFileToModelDb", ns)
  assign("addFileToModelDb", replacement, envir = ns)
  on.exit(
    {
      assign("addFileToModelDb", orig, envir = ns)
      lockBinding("addFileToModelDb", ns)
    },
    add = TRUE
  )

  cold <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = list())
  expect_equal(count, 1L)
  expect_equal(nrow(cold$modeldb), 1L)
  expect_identical(cold$modeldb$name, "PK_1cmt")
  expect_named(cold$entries, "PK_1cmt.R")
  expect_equal(nchar(cold$entries[["PK_1cmt.R"]]$hash), 32L)

  warm <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = cold$entries)
  expect_equal(count, 1L) # unchanged — no re-parse on cache hit
  expect_equal(warm$modeldb, cold$modeldb)
  expect_identical(warm$entries, cold$entries)
})

test_that(".addDirToModelDbCached re-parses a file whose content changes", {
  skip_if_not_installed("nlmixr2est")

  fixtureDir <- tempfile("modeldb-fixture")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE), add = TRUE)

  srcFile <- system.file("modeldb/PK_1cmt.R", package = "nlmixr2lib")
  if (!nzchar(srcFile)) skip("PK_1cmt fixture not installed")
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
  if (!nzchar(srcFile)) skip("PK_1cmt fixture not installed")
  file.copy(srcFile, file.path(fixtureDir, "PK_1cmt.R"))

  staleEntries <- list(
    "ghost.R" = list(hash = "deadbeef", row = data.frame(name = "ghost"))
  )
  result <- nlmixr2lib:::.addDirToModelDbCached(fixtureDir, entries = staleEntries)
  expect_false("ghost.R" %in% names(result$entries))
  expect_true("PK_1cmt.R" %in% names(result$entries))
})
