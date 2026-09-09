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
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### real_token (**a duplicate of the entry above**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### orphan_cite (**cites a file that does not exist**)",
    "- **Type:** binary",
    "- **Example models:** `Nonexistent_1999_ghost.R`.",
    "",
    "### unused_canonical (**no model uses this**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### no_examples_here (**has no example line**)",
    "- **Type:** binary",
    "- **Description:** nothing.",
    "",
    "### bad_xref (**links nowhere**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "- **Notes:** see [[does_not_exist]].",
    "",
    "### no_type_here (**declares no Type at all**)",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### bogus_type_here (**Type outside the known vocabulary**)",
    "- **Type:** wibble",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### deprecated_one (**DEPRECATED -- superseded by `real_token`**)",
    "- **Notes:** tombstone; no examples, no use and no Type is correct here."
  ), file.path(tmp, "inst", "references", "covariate-columns.md"))

  issues <- checkNamingRegisters(tmp)
  for (chk in c("duplicate-canonical", "orphan-example-model",
                "registered-but-unused", "no-example-model", "broken-xref",
                "no-type", "unknown-type")) {
    expect_true(chk %in% issues$check, info = chk)
  }
  # A missing `Type:` is the quietest defect of the set -- it drops the entry
  # out of the canonical list checkModelConventions() builds while every other
  # register check stays green -- so pin which entry each fires on.
  expect_true("no_type_here" %in% issues$name[issues$check == "no-type"])
  expect_true("bogus_type_here" %in% issues$name[issues$check == "unknown-type"])
  # The tombstone must NOT be reported: it legitimately has neither an example,
  # nor a use, nor a Type, and flagging it would push authors to delete
  # deprecation records or to fabricate examples and types for them.
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
    "- **Type:** compartment",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "## Metabolite suffixes",
    "",
    "# The same token may be BOTH a bare compartment and a suffix, so the",
    "# duplicate check is per-section outside covariate-columns.md.",
    "### igg (**suffix form**)",
    "- **Type:** metabolite-suffix",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "# A suffix appears in source only inside a compound token (`central_dox`),",
    "# never bare, so the usage check has to match that form.",
    "### dox (**suffix used only as _dox**)",
    "- **Type:** metabolite-suffix",
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

test_that("checkNamingRegisters enforces the covariate ALL-CAPS convention", {
  # The `## Case convention` section of covariate-columns.md used to name its
  # exceptions in a hand-maintained sentence. It said there were two
  # (`dilution`, `nonECZTRA`) when the file held ten, and had said so through
  # every register audit, because nothing read it. It is now a machine-read
  # list and this is what reads it -- in BOTH directions, since the drift that
  # prompted the check was a documented list falling behind the file rather
  # than a bad name being added.
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "inst", "references"), recursive = TRUE)
  dir.create(file.path(tmp, "inst", "modeldb"), recursive = TRUE)
  writeLines("cl <- exp(lcl) * (1 + e_mixedCase * mixedCase + e * LEGACY_ok)",
             file.path(tmp, "inst", "modeldb", "Real_2020_drug.R"))

  register <- c(
    "## Case convention",
    "",
    "Prose naming `NOT_AN_EXEMPTION` on a parent bullet must not be read as one:",
    "",
    "- **Legacy source-preserved**:",
    "  - `legacy_lower`",
    "  - `gone_missing`",
    "",
    "## Entries",
    "",
    "### LEGACY_OK (**all caps, never flagged**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### legacy_lower (**lower case, but documented above**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### mixedCase (**mixed case and NOT documented above**)",
    "- **Type:** binary",
    "- **Example models:** `Real_2020_drug.R`.",
    "",
    "### old_lower (**DEPRECATED -- superseded by `LEGACY_OK`**)",
    "- **Notes:** a tombstone records history, so its case is not a choice."
  )
  writeLines(register, file.path(tmp, "inst", "references",
                                 "covariate-columns.md"))
  issues <- checkNamingRegisters(tmp)
  caseIss <- issues[issues$check == "case-convention", ]

  # Fires on exactly the undocumented one.
  expect_equal(caseIss$name, "mixedCase")
  # An exemption listed but no longer present in the file is the drift that
  # made the old sentence wrong, so it is reported too.
  expect_equal(issues$name[issues$check == "stale-case-exemption"],
               "gone_missing")
  # Parent bullets carry prose, not names; reading them would silently widen
  # the exemption set to every backticked token in the section.
  expect_false("NOT_AN_EXEMPTION" %in% issues$name)
  # A DEPRECATED tombstone is exempt, as it is from the example/type checks.
  expect_false("old_lower" %in% caseIss$name)
})

test_that("the covariate register's documented case exceptions are complete", {
  # Pins the shipped list itself, not just the mechanism: every non-ALL-CAPS
  # covariate canonical must be one a reader was told about. Enumerates the
  # register rather than sampling it, so a newly-merged lower-case canonical
  # fails here until it is either renamed or documented.
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  path <- file.path(root, "inst", "references", "covariate-columns.md")
  skip_if(!file.exists(path), "register not present in an installed package")

  exempt <- nlmixr2lib:::.caseExemptions(path)
  expect_true(length(exempt) > 0L)

  entries <- nlmixr2lib:::.parseRegister(path)
  real <- Filter(function(e) !isTRUE(e$pseudo) && !isTRUE(e$deprecated), entries)
  canonicals <- unlist(lapply(real, `[[`, "names"))
  expect_setequal(canonicals[canonicals != toupper(canonicals)], exempt)
})
