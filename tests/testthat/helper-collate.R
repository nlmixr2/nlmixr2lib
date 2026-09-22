# Support for the locale-independence gate (test-collate.R).
#
# The gate compares a generator's output under the C collation against its
# output under some other collation. A locale that happens to collate like C
# would make every one of those comparisons vacuous -- they could never go
# red -- so a disagreeing locale is probed for rather than assumed, and the
# tests skip when none is installed.

# C orders every upper-case initial ahead of every lower-case one, so "ABT"
# sorts before "Abacavir". en_US.UTF-8 and its relatives interleave case and
# put "Abacavir" first.
.collateProbe <- c("Abacavir (Archary 2019)", "ABT 102 (Othman 2013)")

.collateDiffersFromC <- function() {
  !identical(sort(.collateProbe), sort(.collateProbe, method = "radix"))
}

# The first installed locale whose collation actually disagrees with C, or
# NULL. Windows spells these differently from glibc, and a minimal container
# may carry only C/C.UTF-8, so several candidates are tried.
.collateNonCLocale <- function() {
  candidates <- c(
    "en_US.UTF-8",
    "en_US.utf8",
    "English_United States.utf8",
    "English_United States.1252",
    "en_GB.UTF-8",
    "de_DE.UTF-8"
  )
  for (loc in candidates) {
    differs <- tryCatch(
      withr::with_collate(loc, .collateDiffersFromC()),
      warning = function(w) FALSE,
      error = function(e) FALSE
    )
    if (isTRUE(differs)) {
      return(loc)
    }
  }
  NULL
}

# Run `fn` twice, once under each collation, and return both results plus the
# name of the non-C locale used. Skips the calling test if there is no locale
# to contrast C against.
.collateBothLocales <- function(fn) {
  loc <- .collateNonCLocale()
  testthat::skip_if(
    is.null(loc),
    "no installed locale collates differently from C; the comparison would be vacuous"
  )
  list(
    C = withr::with_collate("C", fn()),
    other = withr::with_collate(loc, fn()),
    locale = loc
  )
}

# A package directory that .writePkgdownNavbar() can refresh: a `_pkgdown.yml`
# carrying the three AUTOGEN signposts, and a vignettes tree whose basenames
# collate differently under C than under en_US.
.collateNavbarFixture <- function(dir) {
  vign <- file.path(dir, "vignettes")
  dir.create(file.path(vign, "articles"), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(vign, paste0(c("Zebra", "apple", "Apple2", "banana"), ".Rmd")))
  file.create(file.path(vign, "articles", paste0(c("Zed_2020_x", "abacavir", "ABT_102"), ".Rmd")))
  writeLines(
    c(
      "navbar:",
      "  components:",
      "    specific_drugs:",
      "      menu:",
      "        # AUTOGEN:specific_drugs:BEGIN",
      "        # AUTOGEN:specific_drugs:END",
      "    ddmore:",
      "      menu:",
      "        # AUTOGEN:ddmore:BEGIN",
      "        # AUTOGEN:ddmore:END",
      "articles:",
      "# AUTOGEN:articles:BEGIN",
      "# AUTOGEN:articles:END"
    ),
    file.path(dir, "_pkgdown.yml")
  )
  invisible(dir)
}

# The modeldb columns .writePkgdownNavbar() reads. Covers both menus, the
# `is.na(vignette)` filter and the "other" category that belongs in neither.
.collateNavbarModeldb <- function() {
  data.frame(
    label = c(
      "Abacavir (Archary 2019)",
      "ABT 102 (Othman 2013)",
      "AZD6088 rat (Viberg 2012)",
      "Zidovudine (Zhou 2001)",
      "DDMoRe: lidocaine",
      "DDMoRe: Zolpidem",
      "Aspirin (No Vignette 1900)",
      "Some template"
    ),
    vignette = c(
      "abacavir",
      "ABT_102",
      "AZD6088_rat",
      "zidovudine",
      "NA_NA_lidocaine",
      "NA_NA_Zolpidem",
      NA_character_,
      "template"
    ),
    category = c(
      rep("specificDrugs", 4),
      rep("ddmore", 2),
      "specificDrugs",
      "other"
    ),
    stringsAsFactors = FALSE
  )
}

# Regenerate the package's real `_pkgdown.yml` in a scratch copy. Empty
# stand-ins carry the vignette names, so the navbar rewrite runs against the
# committed inputs without the ~45 minute model parse buildModelDb() needs.
.collateRegenerateNavbar <- function(root) {
  dir <- tempfile("pkgdownRegen")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  dir.create(file.path(dir, "vignettes", "articles"), recursive = TRUE)
  file.copy(file.path(root, "_pkgdown.yml"), file.path(dir, "_pkgdown.yml"))
  .collateStubRmd(file.path(root, "vignettes"), file.path(dir, "vignettes"))
  .collateStubRmd(
    file.path(root, "vignettes", "articles"),
    file.path(dir, "vignettes", "articles")
  )
  suppressMessages(
    nlmixr2lib:::.writePkgdownNavbar(nlmixr2lib::modeldb, dir)
  )
  readLines(file.path(dir, "_pkgdown.yml"), encoding = "UTF-8", warn = FALSE)
}

.collateStubRmd <- function(from, to) {
  files <- list.files(from, pattern = "\\.Rmd$")
  if (length(files) > 0) {
    file.create(file.path(to, files))
  }
  invisible(NULL)
}

# Every call node in a parsed R file, flattened.
.collateCallNodes <- function(expr) {
  if (is.call(expr)) {
    return(c(
      list(expr),
      unlist(lapply(as.list(expr), .collateCallNodes), recursive = FALSE)
    ))
  }
  if (is.pairlist(expr) || is.expression(expr) || is.list(expr)) {
    return(unlist(lapply(as.list(expr), .collateCallNodes), recursive = FALSE))
  }
  list()
}

.collateCalleeName <- function(cl) {
  fn <- cl[[1]]
  if (is.name(fn)) {
    return(as.character(fn))
  }
  # `base::sort(x)` and `base:::sort(x)` name the callee in the third slot.
  if (is.call(fn) && as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  ""
}

# Deparsed calls to the four locale-collating primitives in `path`. Parsing
# rather than grepping means a `#` in a string, or a call split over several
# lines, cannot hide one.
.collateBareOrderingCalls <- function(path) {
  calls <- .collateCallNodes(parse(path, keep.source = FALSE))
  callee <- vapply(calls, .collateCalleeName, character(1))
  bare <- calls[callee %in% c("sort", "order", "list.files", "Sys.glob")]
  # unname(): flattening the AST leaves names on the call list, which would
  # otherwise make an empty result differ from character(0).
  unname(vapply(bare, function(cl) paste(deparse(cl), collapse = " "), character(1)))
}
