test_that("the naming registers are structurally clean", {
  # This is the gate. checkModelConventions() checks the names a MODEL uses;
  # this checks the REGISTERS that define them, which had never been checked
  # and had accumulated (audit 2026-09-01): 8 duplicate canonicals, 12
  # canonicals no model used, and 86 citations of model files that did not
  # exist -- 42 distinct, and 22 of those had never existed at ANY commit,
  # i.e. filenames an extraction invented rather than the file it wrote.
  #
  # If this fails, read the `detail` column: each row names the register, the
  # entry and the line.
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  skip_if(!dir.exists(file.path(root, "inst", "references")),
          "registers not present (installed package without inst/references)")
  issues <- checkNamingRegisters(root)
  expect_equal(nrow(issues), 0L,
               info = paste0(
                 "\n", paste(utils::capture.output(print(issues)), collapse = "\n")))
})

test_that("checkNamingRegisters detects each defect class it claims to", {
  # A gate that cannot go red is worse than no gate, so prove each check fires
  # on a register that contains exactly that defect.
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "inst", "references"), recursive = TRUE)
  dir.create(file.path(tmp, "inst", "modeldb"), recursive = TRUE)
  writeLines("x <- 1  # token: real_token",
             file.path(tmp, "inst", "modeldb", "Real_2020_drug.R"))

  writeLines(c(
    "## Section one",
    "",
    "### real_token (**canonical that is fine**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### real_token (**a duplicate of the entry above**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### orphan_cite (**cites a file that does not exist**)",
    "- **Example models:** `Nonexistent_1999_ghost.R`.",
    "",
    "### unused_canonical (**no model uses this**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### no_examples_here (**has no example line**)",
    "- **Description:** nothing.",
    "",
    "### bad_xref (**links nowhere**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "- **Notes:** see [[does_not_exist]].",
    "",
    "### deprecated_one (**DEPRECATED -- superseded by `real_token`**)",
    "- **Notes:** tombstone; no examples and no use is correct here."
  ), file.path(tmp, "inst", "references", "covariate-columns.md"))

  issues <- checkNamingRegisters(tmp)
  for (chk in c("duplicate-canonical", "orphan-example-model",
                "registered-but-unused", "no-example-model", "broken-xref")) {
    expect_true(chk %in% issues$check, info = chk)
  }
  # The tombstone must NOT be reported: it legitimately has neither an example
  # nor a use, and flagging it would push authors to delete deprecation
  # records or to fabricate examples for them.
  expect_false("deprecated_one" %in% issues$name)
})

test_that("checkNamingRegisters does not flag legitimate register patterns", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "inst", "references"), recursive = TRUE)
  dir.create(file.path(tmp, "inst", "modeldb"), recursive = TRUE)
  writeLines(c("d/dt(central_dox) <- -k * central_dox",
               "d/dt(igg) <- -kel * igg"),
             file.path(tmp, "inst", "modeldb", "Real_2020_drug.R"))

  writeLines(c(
    "## Bare compartments",
    "",
    "### igg (**bare compartment**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "## Metabolite suffixes",
    "",
    "# The same token may be BOTH a bare compartment and a suffix, so the",
    "# duplicate check is per-section outside covariate-columns.md.",
    "### igg (**suffix form**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "# A suffix appears in source only inside a compound token (`central_dox`),",
    "# never bare, so the usage check has to match that form.",
    "### dox (**suffix used only as _dox**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "# `###` is also used for policy notes and patterns; neither is a canonical.",
    "### ROUTE_* family -- section-header policy",
    "",
    "### `<tissue>_slab<n>` (**a pattern, not a canonical**)"
  ), file.path(tmp, "inst", "references", "compartment-names.md"))

  issues <- checkNamingRegisters(tmp)
  expect_equal(nrow(issues), 0L,
               info = paste0(
                 "\n", paste(utils::capture.output(print(issues)), collapse = "\n")))
})
