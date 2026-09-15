# Tests for the `- **Type:**` field of the inst/references/*.md registers.
#
# Four places in the package read that one field and they used to disagree
# about a trailing parenthetical qualifier. Only
# checkNamingRegisters.R::.parseRegister() split it off; the three others kept
# it, so `WHO_PS` -- written `continuous (semantically ordinal but treated as
# continuous in the covariate model)` -- stored that whole sentence as its
# routing tag where every other covariate stored a bare token. Nothing
# branched on a covariate's type at the time, which is exactly why it went
# unnoticed: the first code to filter covariates by type would have silently
# dropped that entry.
#
# The tag is matched with identical() wherever it IS consumed
# (.namesByType(), and the `name \r type` duplicate key in
# checkModelConventions.R), so a qualifier left on it drops an entry out of
# the canonical lists with no warning anywhere. These tests are the mechanical
# gate that keeps all four readers splitting the field the same way.

test_that(".splitRegisterType splits a tag from its qualifier", {
  f <- nlmixr2lib:::.splitRegisterType
  expect_equal(f("continuous"), list(type = "continuous", qualifier = ""))
  expect_equal(f("  metabolite-suffix  "), list(type = "metabolite-suffix", qualifier = ""))
  expect_equal(f("binary (only ever 0 or 1)"), list(type = "binary", qualifier = "only ever 0 or 1"))
  # Greedy to the LAST close paren, so a nested parenthetical stays whole.
  expect_equal(
    f("count (integer sum of items (see Notes))"),
    list(type = "count", qualifier = "integer sum of items (see Notes)")
  )
  # A missing Type line reaches this as NA and must stay NA rather than
  # becoming the empty string, which .parseRegister() treats as "declared".
  expect_equal(f(NA_character_), list(type = NA_character_, qualifier = ""))
  expect_equal(f(character()), list(type = NA_character_, qualifier = ""))
})

test_that("every covariate entry's type is a bare routing token", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  expect_gt(length(canon), 1000)
  types <- vapply(canon, function(e) e$type %||% "", character(1))
  # A bare token is lower-case letters and internal hyphens, nothing else --
  # no whitespace, no parentheses. An empty type is legitimate only for a
  # DEPRECATED tombstone, which carries no `Type:` line at all.
  offenders <- names(canon)[
    nzchar(types) &
      !grepl("^[a-z]+(-[a-z]+)*$", types)
  ]
  expect_equal(offenders, character())
})

test_that("covariate types stay inside the register's documented vocabulary", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  types <- unique(vapply(canon, function(e) e$type %||% "", character(1)))
  # .knownTypes is a closed set precisely so a new value has to be ratified
  # rather than minted in passing. Every covariate type must now be in it --
  # no allowance list, so there is nothing here to rot.
  #
  # test-checkNamingRegisters.R has a stronger per-file version of this via
  # .schemaTypes(), but it reads the register through .parseRegister() only.
  # This one reads it through .parseCovariateColumns(), which is the parser
  # that actually feeds conventions$canonicalCovariates -- so the two cover
  # the same register through different readers, which is the whole point
  # given those two readers once disagreed.
  expect_equal(setdiff(types, nlmixr2lib:::.knownTypes), character())
})

test_that("the covariate and register parsers agree on every shared type", {
  path <- system.file("references", "covariate-columns.md", package = "nlmixr2lib")
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
    "semantically ordinal but treated as continuous in the covariate model"
  )
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
  writeLines(
    c(
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
    ),
    tmp
  )
  e <- nlmixr2lib:::.parseCovariateColumns(tmp)

  expect_equal(e$PLAIN$type, "continuous")
  expect_equal(e$PLAIN$typeQualifier, "")
  expect_equal(e$QUALIFIED$type, "binary")
  expect_equal(e$QUALIFIED$typeQualifier, "only ever 0 or 1 in the founding model")
  # A nested parenthetical keeps the inner parentheses in the qualifier and
  # still yields a bare leading token.
  expect_equal(e$NESTED$type, "count")
  expect_equal(e$NESTED$typeQualifier, "integer sum of per-domain items (see Notes)")
})

# Every register, every reader ------------------------------------------------

test_that("no parser stores a parenthetical in the routing tag", {
  paths <- list.files(system.file("references", package = "nlmixr2lib"), pattern = "\\.md$", full.names = TRUE)
  # Several files under references/ are follow-up notes rather than registers
  # and carry no `Type:` field; discovering the registers by content rather
  # than by name means a NEW register is covered the day it is added.
  hasType <- function(path) {
    any(grepl("^- \\*\\*Type:\\*\\*", readLines(path, warn = FALSE)))
  }
  registers <- Filter(hasType, paths)
  expect_gte(length(registers), 3L)
  for (path in registers) {
    # .parseRegister() reads every register; .parseTypedNamesMd() reads the
    # compartment and parameter ones; .parseCovariateColumns() the covariate
    # one. Run all three over each -- a parser that owns a different file
    # still parses the shared `### name` / `- **Type:**` shape, so all three
    # readings of every register get checked.
    got <- c(
      lapply(nlmixr2lib:::.parseRegister(path), `[[`, "type"),
      lapply(nlmixr2lib:::.parseTypedNamesMd(path), `[[`, "type"),
      lapply(nlmixr2lib:::.parseCovariateColumns(path), `[[`, "type")
    )
    got <- as.character(unlist(got[!vapply(got, is.null, logical(1))]))
    got <- got[!is.na(got) & nzchar(got)]
    expect_gt(length(got), 0L)
    bad <- as.character(unique(got[!grepl("^[a-z]+(-[a-z]+)*$", got)]))
    expect_equal(bad, character(), info = paste("parenthetical left on a routing tag in", basename(path)))
  }
})

test_that("registered metabolite suffixes are all bare tokens", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  # .namesByType() matches with identical(), so any compartment entry whose
  # Type still carried a qualifier would be missing from these vectors
  # entirely rather than merely misspelled.
  expect_gt(length(conv$registeredMetabolites), 200L)
  expect_gt(length(conv$compartments), 500L)
  expect_true(all(nzchar(conv$registeredMetabolites)))
  # The Zurlinden paracetamol suffixes were migrated to the Cook 2016 forms
  # before 0.3.2 shipped, so only the canonical spellings exist now.
  expect_true(all(c("apaps", "apapg") %in% conv$registeredMetabolites))
  expect_false(any(c("as", "ag") %in% conv$registeredMetabolites))
})

test_that("the removed deprecation tombstones are gone from the registers", {
  # `as`, `ag` and `DIAL` were deprecated and replaced entirely within the
  # 0.3.2.9000 development cycle, so no released version ever carried them and
  # a tombstone pointing at the replacement had no audience. `DIAL` survives
  # where it is actually useful -- as a source alias recording what the source
  # papers called the column.
  expect_false(any(
    c("as", "ag") %in%
      vapply(nlmixr2lib:::.loadCanonicalCompartments(), `[[`, character(1), "name")
  ))
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  expect_false("DIAL" %in% names(canon))
  expect_true("DIAL" %in% canon$RRT_HEMODIAL_ACTIVE$aliases)
  expect_true("DIAL" %in% canon$RRT_HEMODIAL_STATUS$aliases)
})
