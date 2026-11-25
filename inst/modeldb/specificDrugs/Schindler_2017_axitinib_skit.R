Schindler_2017_axitinib_skit <- function() {
  description <- "Empirical time-course model of soluble stem cell factor receptor (sKIT) in Japanese adults with cytokine-refractory metastatic renal cell carcinoma on oral axitinib (Schindler 2017, Supplementary Material 1 Equation S1). sKIT changes linearly with time from an individual baseline, sKIT(t) = Base + c * t, with no axitinib effect identified; baseline has exponential IIV, the slope has additive IIV, and residual error is additive on the log scale."
  reference <- paste(
    "Schindler E, Amantea MA, Karlsson MO, Friberg LE.",
    "A Pharmacometric Framework for Axitinib Exposure, Efficacy, and Safety",
    "in Metastatic Renal Cell Carcinoma Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(6):373-382.",
    "doi:10.1002/psp4.12193. Supplementary Material 1, 'Final sKIT model'.",
    "Sister model files from the same paper:",
    "modellib('Schindler_2017_axitinib_biomarkers_sld'),",
    "modellib('Schindler_2017_axitinib_dbp'),",
    "modellib('Schindler_2017_axitinib_os').",
    sep = " "
  )
  vignette <- "Schindler_2017_axitinib"
  units <- list(time = "day", dosing = "mg", concentration = "pg/mL (sKIT)")

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 1L,
    age_range = "adults (age not reported in Schindler 2017)",
    weight_range = "not reported in Schindler 2017",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Cytokine-refractory metastatic renal cell carcinoma (mRCC)",
    dose_range = "Axitinib 5 mg PO twice daily starting dose, titrated 2-10 mg twice daily.",
    regions = "Japan",
    biomarkers = "sKIT, n = 436 observations; median follow-up 168 days (range 32-624).",
    notes = "Phase II study NCT00569946 (Tomita et al. Eur J Cancer 2011;47:2592-2602). No axitinib effect on sKIT was identified, so the sKIT model was not included in the joint biomarker model and sKIT was not tested as a predictor of tumour size or survival (Schindler 2017 Results 'Biomarkers models')."
  )

  ini({
    # Schindler 2017 Supplementary Material 1 'Final sKIT model' (Equation S1).
    lrbase_skit <- log(37700); label("Baseline sKIT concentration (pg/mL)") # Supplementary Material 1: Base_sKIT = 37700 pg/mL (RSE 3.0%)
    slope_skit <- 7.71; label("Linear change of sKIT with time, c (pg/mL/day)") # Supplementary Material 1: c = 7.71 pg/mL/day (RSE 62%)

    etalrbase_skit ~ 0.0441 # Supplementary Material 1: IIV on Base_sKIT 21% CV -> 0.21^2 (same %CV = 100 * sqrt(omega^2) convention as Table 2)
    etaslope_skit ~ 361 # Supplementary Material 1: additive IIV on c with SD 19 pg/mL/day -> 19^2

    expSd_skit <- 0.124; label("Additive residual error on the log scale for sKIT") # Supplementary Material 1: 'additive error on log scale ... corresponding to a proportional error of 12.4%' (RSE 5.4%)
  })

  model({
    rbase_skit <- exp(lrbase_skit + etalrbase_skit)
    c_skit <- slope_skit + etaslope_skit
    # Equation S1: sKIT(t) = Base_sKIT + c * t.
    skit <- rbase_skit + c_skit * t
    skit ~ lnorm(expSd_skit)
  })
}
