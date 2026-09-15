# Locale-independence gate for the generated artifacts.
#
# `buildModelDb()` writes data/modeldb.rda, inst/modeldb.rds and the
# auto-generated regions of _pkgdown.yml. An ordering that follows LC_COLLATE
# makes those artifacts depend on the locale of the machine that rebuilt them:
# C and en_US.UTF-8 disagree about case, so most of the modeldb rows and most
# of the navbar entries move between the two, and a rebuild reads as a diff of
# thousands of lines that says nothing about the change being made. These
# tests re-run the generators under a collation that disagrees with C and
# require byte-identical output.

test_that("the collating wrappers ignore LC_COLLATE", {
  labels <- c("Abacavir (Archary 2019)", "ABT 102 (Othman 2013)", "azd (x 2001)")
  res <- .collateBothLocales(function() {
    list(
      sorted = nlmixr2lib:::.collateSort(labels),
      ordered = nlmixr2lib:::.collateOrder(labels)
    )
  })
  expect_identical(res$other, res$C, info = res$locale)
  # Pinned, not merely stable: C collation is the order that is committed.
  expect_identical(
    res$C$sorted,
    c("ABT 102 (Othman 2013)", "Abacavir (Archary 2019)", "azd (x 2001)")
  )
  expect_identical(res$C$ordered, c(2L, 1L, 3L))
})

test_that("the collating wrappers re-sort the directory enumerators", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("Abb.R", "ABa.R", "azz.R", "AZa.R")))
  res <- .collateBothLocales(function() {
    list(
      files = nlmixr2lib:::.collateListFiles(dir, pattern = "\\.R$"),
      glob = basename(nlmixr2lib:::.collateGlob(file.path(dir, "*.R")))
    )
  })
  expect_identical(res$other, res$C, info = res$locale)
  expect_identical(res$C$files, c("ABa.R", "AZa.R", "Abb.R", "azz.R"))
  expect_identical(res$C$glob, c("ABa.R", "AZa.R", "Abb.R", "azz.R"))
})

test_that("the collating wrappers accept empty input", {
  expect_identical(nlmixr2lib:::.collateSort(character(0)), character(0))
  expect_identical(nlmixr2lib:::.collateOrder(character(0)), integer(0))
  dir <- withr::local_tempdir()
  expect_identical(nlmixr2lib:::.collateListFiles(dir, pattern = "\\.R$"), character(0))
  expect_identical(nlmixr2lib:::.collateGlob(file.path(dir, "*.R")), character(0))
})

test_that(".writePkgdownNavbar() writes the same bytes in any collation locale", {
  res <- .collateBothLocales(function() {
    dir <- withr::local_tempdir()
    .collateNavbarFixture(dir)
    suppressMessages(
      nlmixr2lib:::.writePkgdownNavbar(.collateNavbarModeldb(), dir)
    )
    readLines(file.path(dir, "_pkgdown.yml"), encoding = "UTF-8", warn = FALSE)
  })
  expect_identical(res$other, res$C, info = res$locale)
  expect_identical(
    res$C,
    c(
      "navbar:",
      "  components:",
      "    specific_drugs:",
      "      menu:",
      "        # AUTOGEN:specific_drugs:BEGIN",
      '        - text: "ABT 102 (Othman 2013)"',
      "          href: articles/ABT_102.html",
      '        - text: "AZD6088 rat (Viberg 2012)"',
      "          href: articles/AZD6088_rat.html",
      '        - text: "Abacavir (Archary 2019)"',
      "          href: articles/abacavir.html",
      '        - text: "Zidovudine (Zhou 2001)"',
      "          href: articles/zidovudine.html",
      "        # AUTOGEN:specific_drugs:END",
      "    ddmore:",
      "      menu:",
      "        # AUTOGEN:ddmore:BEGIN",
      '        - text: "DDMoRe: Zolpidem"',
      "          href: articles/NA_NA_Zolpidem.html",
      '        - text: "DDMoRe: lidocaine"',
      "          href: articles/NA_NA_lidocaine.html",
      "        # AUTOGEN:ddmore:END",
      "articles:",
      "# AUTOGEN:articles:BEGIN",
      "- title: General",
      "  desc: Cross-cutting guides and the list of models.",
      "  contents:",
      "  - Apple2",
      "  - Zebra",
      "  - apple",
      "  - banana",
      "- title: internal",
      paste(
        "  desc: Drug-specific validation vignettes. They are reached from the",
        "Specific drug models and DDMoRe models navbar dropdowns; this group",
        "hides them from the main Articles index."
      ),
      "  contents:",
      "  - articles/ABT_102",
      "  - articles/Zed_2020_x",
      "  - articles/abacavir",
      "  internal: true",
      "# AUTOGEN:articles:END"
    )
  )
})

test_that("the committed modeldb rows are in C-collation order", {
  # The row order is the order inst/modeldb was walked, and it is frozen into
  # data/modeldb.rda and inst/modeldb.rds. If a rebuild on a differently
  # collating machine reordered it, this is what catches it.
  filenames <- nlmixr2lib::modeldb$filename
  expect_identical(filenames, nlmixr2lib:::.collateSort(filenames))
  res <- .collateBothLocales(function() nlmixr2lib:::.collateSort(filenames))
  expect_identical(res$other, filenames, info = res$locale)
})

test_that("the committed _pkgdown.yml regenerates byte-for-byte in any locale", {
  root <- testthat::test_path("..", "..")
  skip_if(
    !file.exists(file.path(root, "_pkgdown.yml")),
    "_pkgdown.yml not present (installed package)"
  )
  skip_if(
    !dir.exists(file.path(root, "vignettes")),
    "vignettes/ not present (installed package)"
  )
  committed <- readLines(file.path(root, "_pkgdown.yml"), encoding = "UTF-8", warn = FALSE)
  res <- .collateBothLocales(function() .collateRegenerateNavbar(root))
  expect_identical(res$C, committed)
  expect_identical(res$other, committed, info = res$locale)
})

test_that("the generator orders nothing with the collation locale", {
  # R/modeldb.R is the only file that writes a committed artifact. Every
  # ordering in it must go through R/collate.R, so a bare sort(), order(),
  # list.files() or Sys.glob() added later fails here instead of surfacing as
  # a mystery diff months afterwards.
  src <- testthat::test_path("..", "..", "R", "modeldb.R")
  skip_if(!file.exists(src), "package source not present (installed package)")
  expect_identical(.collateBareOrderingCalls(src), character(0))
})

test_that("the ordering scanner finds what it is meant to find", {
  # Without this the test above could be green because the scanner is broken
  # rather than because the generator is clean.
  src <- withr::local_tempfile(fileext = ".R")
  writeLines(
    c(
      "f <- function(x) {",
      "  # sort(x) in a comment does not count",
      '  y <- "neither does order(x) in a string"',
      "  z <- .collateSort(x)",
      "  a <- sort(",
      "    x",
      "  )",
      "  b <- base::order(x)",
      '  d <- list.files(x, pattern = "[.]R$")',
      "  e <- Sys.glob(x)",
      "  c(z, a, b, d, e)",
      "}"
    ),
    src
  )
  expect_identical(
    .collateBareOrderingCalls(src),
    c("sort(x)", "base::order(x)", 'list.files(x, pattern = "[.]R$")', "Sys.glob(x)")
  )
})
