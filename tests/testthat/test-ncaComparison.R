test_that("ncaParamLabel maps known codes to friendly labels", {
  expect_equal(
    ncaParamLabel(c("cmax", "tmax", "auclast", "aucinf.obs", "half.life")),
    c("Cmax", "Tmax", "AUClast", "AUC0-∞ (obs)", "t½")
  )
})

test_that("ncaParamLabel appends units when supplied", {
  expect_equal(
    ncaParamLabel(c("cmax", "auclast"),
                  units = c(cmax = "ng/mL", auclast = "ng*h/mL")),
    c("Cmax (ng/mL)", "AUClast (ng*h/mL)")
  )
})

test_that("ncaParamLabel passes through unknown codes with a warning", {
  expect_warning(
    out <- ncaParamLabel(c("cmax", "totally_fake_param")),
    "unknown PKNCA code"
  )
  expect_equal(out, c("Cmax", "totally_fake_param"))
})

test_that("ncaParamLabel preserves NA values", {
  expect_equal(suppressWarnings(ncaParamLabel(NA_character_)), NA_character_)
})

test_that("ncaComparisonTable builds a side-by-side frame from long sim + wide ref", {
  simulated <- data.frame(
    treatment = rep(c("50 mg", "100 mg"), each = 3),
    PPTESTCD  = rep(c("cmax", "tmax", "auclast"), 2),
    PPORRES   = c(15.2, 2.0, 96.0, 29.1, 2.1, 191.0)
  )
  reference <- data.frame(
    treatment = c("50 mg", "100 mg"),
    cmax      = c(14.8, 28.5),
    tmax      = c(2.0,  2.1),
    auclast   = c(95.0, 190.0)
  )
  tbl <- ncaComparisonTable(
    simulated, reference,
    by = "treatment",
    units = c(cmax = "ng/mL", auclast = "ng*h/mL")
  )
  expect_named(
    tbl,
    c("NCA parameter", "treatment", "Reference", "Simulated", "% diff")
  )
  expect_equal(nrow(tbl), 6L)
  expect_true(all(c("Cmax (ng/mL)", "Tmax", "AUClast (ng*h/mL)") %in%
                    tbl[["NCA parameter"]]))
  expect_equal(tbl[["NCA parameter"]][1L], "Cmax (ng/mL)")
  expect_equal(tbl[["NCA parameter"]][3L], "Tmax")
})

test_that("ncaComparisonTable handles ungrouped input", {
  simulated <- data.frame(
    PPTESTCD  = c("cmax", "auclast"),
    PPORRES   = c(15.2, 96.0)
  )
  reference <- data.frame(cmax = 14.8, auclast = 95.0)
  tbl <- ncaComparisonTable(simulated, reference)
  expect_named(tbl, c("NCA parameter", "Reference", "Simulated", "% diff"))
  expect_equal(nrow(tbl), 2L)
  expect_equal(tbl[["NCA parameter"]], c("Cmax", "AUClast"))
})

test_that("ncaComparisonTable aggregates per-subject simulated rows via median", {
  simulated <- data.frame(
    id        = rep(1:4, each = 1L),
    PPTESTCD  = "cmax",
    PPORRES   = c(10, 12, 14, 16)
  )
  reference <- data.frame(cmax = 13.0)
  tbl <- ncaComparisonTable(simulated, reference)
  expect_equal(tbl[["Simulated"]], "13")
})

test_that("ncaComparisonTable flags rows exceeding the tolerance", {
  simulated <- data.frame(PPTESTCD = c("cmax", "tmax"),
                          PPORRES  = c(18.0, 2.0))
  reference <- data.frame(cmax = 10.0, tmax = 2.0)
  tbl <- ncaComparisonTable(simulated, reference, tolerance_pct = 20)
  expect_match(tbl[["% diff"]][1L], "\\*$")
  expect_false(grepl("\\*$", tbl[["% diff"]][2L]))
  expect_match(attr(tbl, "footnote"), "exceeds|differs from reference")
})

test_that("ncaComparisonTable: tolerance disabled with Inf", {
  simulated <- data.frame(PPTESTCD = "cmax", PPORRES = 100)
  reference <- data.frame(cmax = 10)
  tbl <- ncaComparisonTable(simulated, reference, tolerance_pct = Inf)
  expect_false(grepl("\\*$", tbl[["% diff"]][1L]))
  expect_null(attr(tbl, "footnote"))
})

test_that("ncaComparisonTable formats NA values with em-dash", {
  simulated <- data.frame(PPTESTCD = c("cmax", "tmax"),
                          PPORRES  = c(NA_real_, 2.0))
  reference <- data.frame(cmax = NA_real_, tmax = 2.0)
  tbl <- ncaComparisonTable(simulated, reference)
  expect_equal(tbl[["Reference"]][match("Cmax", tbl[["NCA parameter"]])], "—")
  expect_equal(tbl[["Simulated"]][match("Cmax", tbl[["NCA parameter"]])], "—")
  expect_equal(tbl[["% diff"]][match("Cmax", tbl[["NCA parameter"]])], "—")
})

test_that("ncaComparisonTable respects params subset", {
  simulated <- data.frame(
    PPTESTCD = c("cmax", "tmax", "auclast"),
    PPORRES  = c(15.0, 2.0, 95.0)
  )
  reference <- data.frame(cmax = 14.0, tmax = 2.0, auclast = 90.0)
  tbl <- ncaComparisonTable(simulated, reference, params = c("cmax", "auclast"))
  expect_equal(nrow(tbl), 2L)
  expect_setequal(tbl[["NCA parameter"]], c("Cmax", "AUClast"))
})

test_that("ncaComparisonTable errors when no parameters overlap", {
  simulated <- data.frame(PPTESTCD = "cmax", PPORRES = 10)
  reference <- data.frame(tmax = 2)
  expect_error(
    ncaComparisonTable(simulated, reference),
    "No NCA parameters overlap"
  )
})

test_that("ncaComparisonTable supports custom label_first_column", {
  simulated <- data.frame(PPTESTCD = "cmax", PPORRES = 10)
  reference <- data.frame(cmax = 10)
  tbl <- ncaComparisonTable(simulated, reference,
                            label_first_column = "Parameter")
  expect_named(tbl, c("Parameter", "Reference", "Simulated", "% diff"))
})

test_that("ncaComparisonTable orders rows by canonical PKNCA parameter order", {
  simulated <- data.frame(
    PPTESTCD = c("half.life", "cmax", "tmax", "auclast"),
    PPORRES  = c(10, 15, 2, 95)
  )
  reference <- data.frame(half.life = 11, cmax = 14, tmax = 2, auclast = 90)
  tbl <- ncaComparisonTable(simulated, reference)
  expect_equal(tbl[["NCA parameter"]],
               c("Cmax", "Tmax", "AUClast", "t½"))
})
