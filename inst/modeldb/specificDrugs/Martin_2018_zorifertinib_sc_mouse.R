Martin_2018_zorifertinib_sc_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, subcutaneous NSCLC xenograft). Dose-driven tumour",
    "growth inhibition model for zorifertinib (AZD3759), an EGFR tyrosine",
    "kinase inhibitor: exponential tumour-volume growth with a sigmoid Emax",
    "cell-kill term driven by the daily dose in mg/kg (no PK layer; the",
    "effect is constant under daily dosing). Martin 2018, Equation 1 and",
    "Table 1."
  )
  reference <- paste(
    "Martin EC, Aarons L, Yates JWT. Pharmacodynamic modelling of resistance",
    "to epidermal growth factor receptor inhibition in brain metastasis mouse",
    "models. Cancer Chemother Pharmacol. 2018;82(4):669-675.",
    "doi:10.1007/s00280-018-3630-8"
  )
  vignette <- "Martin_2018_zorifertinib"
  units <- list(
    time = "day",
    dosing = "mg/kg (daily dose, supplied as the DOSE_ZORIFERTINIB_MGKGD covariate)",
    concentration = "mm^3 (tumour volume tumor_size)"
  )

  covariateData <- list(
    DOSE_ZORIFERTINIB_MGKGD = list(
      description = "Daily zorifertinib (AZD3759) dose of the animal's treatment group",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Paper symbol D in Equation 1. 0 for vehicle controls; the study tested",
        "3.75, 7.5, 15 and 30 mg/kg once daily. The model has no PK layer: the",
        "drug effect is a function of the daily dose and is assumed constant",
        "throughout the treatment period (Methods)."
      ),
      source_name = "D"
    )
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "tumour volume (calliper, length x width^2 / 2)",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (subcutaneous NSCLC xenograft)",
    n_subjects = 38L,
    n_studies = 1L,
    disease_state = paste(
      "Subcutaneous human non-small cell lung cancer xenografts (about 15 mm^3",
      "of cells implanted; about 180 mm^3 at the start of dosing)"
    ),
    dose_range = "zorifertinib (AZD3759) 3.75, 7.5, 15 or 30 mg/kg once daily, or vehicle, for 15 days",
    n_per_group = "7 controls; 7 or 8 per dose group",
    notes = paste(
      "Tumour volume measured by calliper roughly every 4 days for 15 days from",
      "randomisation (Martin 2018, Methods 'Subcutaneous xenograft mouse",
      "models'). The xenograft cell line and mouse strain are not stated in",
      "the source."
    )
  )

  ini({
    lkg <- log(0.0855)
    label("Exponential tumour proliferation rate constant lambda (1/day)") # Table 1: lambda = 0.0855 1/day (RSE 10.8%)
    lemax <- log(0.519)
    label("Maximum drug-induced cell-kill rate constant Emax (1/day)") # Table 1: Emax = 0.519 (RSE 27.8%); printed unit mm3/day is inconsistent with Eq 1, Fig 5 confirms 1/day
    led50 <- log(21.5)
    label("Daily dose giving half-maximal cell kill ED50 (mg/kg)") # Table 1: ED50 = 21.5 mg/kg (RSE 83.2%)
    lhill <- log(0.617)
    label("Hill coefficient of the dose-response gamma (unitless)") # Table 1: Gamma = 0.617 (RSE 45.2%)
    lrbase <- log(175)
    label("Tumour volume at the start of dosing (mm^3)") # Table 1: Baseline = 175 mm3 (RSE 5.5%)

    # IIV: Methods define CV% as the square root of the variance, so omega^2 = (CV/100)^2
    etalkg ~ 0.121801 # Table 1: lambda IIV CV 34.9% -> 0.349^2
    etalrbase ~ 0.093025 # Table 1: Baseline IIV CV 30.5% -> 0.305^2

    # Residual error: log-normal (Methods), 12.7% read as the log-scale SD
    expSd <- 0.127
    label("Log-scale residual SD of tumour volume (unitless)") # Table 1: Residual error = 12.7%
  })
  model({
    kg <- exp(lkg + etalkg)
    emax <- exp(lemax)
    ed50 <- exp(led50)
    hill <- exp(lhill)
    rbase <- exp(lrbase + etalrbase)

    # Equation 1: sigmoid Emax cell kill proportional to tumour volume
    dose_d <- DOSE_ZORIFERTINIB_MGKGD
    kill <- emax * dose_d^hill / (ed50^hill + dose_d^hill)

    tumor_size(0) <- rbase
    d/dt(tumor_size) <- kg * tumor_size - kill * tumor_size

    tumor_size ~ lnorm(expSd)
  })
}
