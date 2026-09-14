# Tests for the `- **Type:**` field of the inst/references/*.md registers.
#
# Four places in the package read that one field, and they used to disagree
# about a trailing parenthetical qualifier. `.parseRegister()` strips it;
# `.parseCovariateColumns()` did not, so `WHO_PS` -- written `continuous
# (semantically ordinal but treated as continuous in the covariate model)` --
# stored the whole sentence where every other covariate stored a bare token.
# Nothing branched on a covariate's type at the time, which is exactly why it
# went unnoticed: the first code to filter covariates by type would have
# silently dropped that entry. These tests are the mechanical gate that keeps
# the parsers reading the field the same way.

test_that("every covariate entry's type is a bare routing token", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  expect_gt(length(canon), 1000)
  types <- vapply(canon, function(e) e$type %||% "", character(1))
  # A bare token is lower-case letters and internal hyphens, nothing else --
  # no whitespace, no parentheses. An empty type is legitimate only for a
  # DEPRECATED tombstone, which carries no `Type:` line at all.
  offenders <- names(canon)[nzchar(types) &
                              !grepl("^[a-z]+(-[a-z]+)*$", types)]
  expect_equal(offenders, character())
})

test_that("covariate types stay inside the register's documented vocabulary", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  types <- unique(vapply(canon, function(e) e$type %||% "", character(1)))
  # Checked in BOTH directions, the way .caseExemptions() is: an unexpected
  # new value fails, and so does the disappearance of a known deviation, so
  # the allowance below cannot rot into a lie after the deviation is resolved.
  #
  # `""` is the DEPRECATED tombstone `DIAL`, which has no `Type:` line.
  # `"ordinal"` is a real, separate, pre-existing defect: it is absent from
  # checkNamingRegisters.R::.knownTypes, so checkNamingRegisters() reports it
  # today (DIS_COPD_GOLD, DIS_COPD_GOLD_LOW / _HIGH). It is listed here rather
  # than silently tolerated -- either ratify `ordinal` into .knownTypes or
  # retype those three entries, then delete it from this vector.
  expect_setequal(setdiff(types, nlmixr2lib:::.knownTypes),
                  c("", "ordinal"))
})

test_that("the covariate and register parsers agree on every shared type", {
  path <- system.file("references", "covariate-columns.md",
                      package = "nlmixr2lib")
  expect_true(nzchar(path))

  regEntries <- nlmixr2lib:::.parseRegister(path)
  regType <- character()
  for (e in regEntries) {
    for (nm in e$names) {
      if (nzchar(nm)) {
        regType[[nm]] <- if (is.na(e$type)) "" else e$type
      }
    }
  }
  covEntries <- nlmixr2lib:::.parseCovariateColumns(path)
  covType <- vapply(covEntries, function(e) e$type %||% "", character(1))

  shared <- intersect(names(covType), names(regType))
  # An enumerating check, not a hand-picked example: every canonical the two
  # parsers both see must get the same routing token from both. Before the
  # split this reported WHO_PS out of 1736 shared names.
  expect_gt(length(shared), 1000)
  expect_equal(unname(covType[shared]), unname(regType[shared]))
})

test_that("a Type parenthetical is split into type and typeQualifier", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  expect_equal(canon$WHO_PS$type, "continuous")
  expect_equal(
    canon$WHO_PS$typeQualifier,
    "semantically ordinal but treated as continuous in the covariate model")
  # A type with no parenthetical gets an empty qualifier, not NULL.
  expect_equal(canon$WT$type, "continuous")
  expect_equal(canon$WT$typeQualifier, "")
  hasScalarQualifier <- function(e) {
    is.character(e$typeQualifier) && length(e$typeQualifier) == 1L
  }
  expect_true(all(vapply(canon, hasScalarQualifier, logical(1))))
})

test_that("the parenthetical split is driven by the file, not by WHO_PS", {
  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "### PLAIN (**canonical for a plain thing**)",
    "- **Units:** kg",
    "- **Type:** continuous",
    "- **Scope:** general",
    "",
    "### QUALIFIED (**canonical for a qualified thing**)",
    "- **Units:** none",
    "- **Type:** binary (only ever 0 or 1 in the founding model)",
    "- **Scope:** specific",
    "",
    "### NESTED (**canonical for a nested parenthetical**)",
    "- **Type:** count (integer sum of per-domain items (see Notes))",
    "- **Scope:** general"
  ), tmp)
  e <- nlmixr2lib:::.parseCovariateColumns(tmp)

  expect_equal(e$PLAIN$type, "continuous")
  expect_equal(e$PLAIN$typeQualifier, "")
  expect_equal(e$QUALIFIED$type, "binary")
  expect_equal(e$QUALIFIED$typeQualifier,
               "only ever 0 or 1 in the founding model")
  # A nested parenthetical keeps the inner parentheses in the qualifier and
  # still yields a bare leading token.
  expect_equal(e$NESTED$type, "count")
  expect_equal(e$NESTED$typeQualifier,
               "integer sum of per-domain items (see Notes)")
})

# Deprecated metabolite suffixes ---------------------------------------------
#
# .parseTypedNamesMd() deliberately does NOT strip the parenthetical, and the
# difference from the two parsers above is load-bearing rather than cosmetic.
# See the NOTE beside the `Type:` branch in R/conventions.R. These tests pin
# the behaviour AND the mechanism, so that "make the parsers consistent" done
# without reading that note fails here instead of silently re-admitting two
# retired suffixes.

test_that("deprecated metabolite suffixes stay out of registeredMetabolites", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  # `as` and `ag` were retired on 2026-06-19 (R-reserved-word and
  # chemistry-symbol collisions) in favour of `apaps` / `apapg`.
  expect_false("as" %in% conv$registeredMetabolites)
  expect_false("ag" %in% conv$registeredMetabolites)
  # The replacements are live, so this is an exclusion of the retired spelling
  # rather than of the chemical species.
  expect_true("apaps" %in% conv$registeredMetabolites)
  expect_true("apapg" %in% conv$registeredMetabolites)
})

test_that("the deprecated exclusion rests on the unstripped Type qualifier", {
  comp <- nlmixr2lib:::.loadCanonicalCompartments()
  byName <- function(nm) Filter(function(e) identical(e$name, nm), comp)
  for (nm in c("as", "ag")) {
    entries <- byName(nm)
    expect_equal(length(entries), 1L)
    # This is the exact string .namesByType() fails to match. If a future
    # change strips it, this test goes red first and points at the NOTE
    # explaining that an explicit `deprecated` flag is the intended fix.
    expect_equal(entries[[1]]$type, "metabolite-suffix (deprecated)")
  }
  expect_false("metabolite-suffix (deprecated)" %in% nlmixr2lib:::.knownTypes)
})
