Martin_2018_zorifertinib_brainmet_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, PC9_Luc NSCLC brain metastasis model). Tumour growth",
    "and resistance model for zorifertinib (AZD3759), an EGFR tyrosine kinase",
    "inhibitor, fitted to bioluminescence: the tumour is a treatment-",
    "sensitive and a treatment-resistant cell population growing",
    "exponentially at separate rates, a fraction of cells being resistant at",
    "baseline (Hypothesis 1, the final model), with a cell-kill term on the",
    "sensitive cells proportional to the daily dose in mg/kg (no PK layer).",
    "Martin 2018, Equations 2-4 and Table 2."
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
    concentration = "1e6 photons/s (total bioluminescence flux tumor_size)"
  )
  paper_specific_compartments <- c("tumor_sensitive", "tumor_resistant")

  covariateData <- list(
    DOSE_ZORIFERTINIB_MGKGD = list(
      description = "Daily zorifertinib (AZD3759) dose of the animal's treatment group",
      units = "mg/kg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Paper symbol D in Equation 2. 0 for vehicle controls; the study tested",
        "7.5 and 15 mg/kg once daily. The model has no PK layer: the cell-kill",
        "rate kD * D is assumed constant throughout the study (Methods)."
      ),
      source_name = "D"
    )
  )

  compartmentData <- list(
    tumor_sensitive = list(
      analyte = "treatment-sensitive tumour cells (bioluminescence)",
      units = "1e6 photons/s",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_resistant = list(
      analyte = "treatment-resistant tumour cells (bioluminescence)",
      units = "1e6 photons/s",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (nude; PC9_Luc NSCLC brain metastasis model)",
    n_subjects = 51L,
    n_studies = 1L,
    disease_state = paste(
      "Brain metastasis model of human NSCLC: PC9 cells (EGFR exon 19",
      "deletion) transfected with a luciferase vector (PC9_Luc), implanted",
      "into immunodeficient nude mice"
    ),
    dose_range = "zorifertinib (AZD3759) 7.5 or 15 mg/kg once daily, or vehicle, for up to 70 days",
    n_per_group = "25 controls; 13 per dose group",
    notes = paste(
      "Tumour burden quantified by bioluminescence imaging (Xenogen); follow-up",
      "up to 70 days from randomisation; about 10% of animals dropped out early",
      "with large tumour burdens (Martin 2018, Methods 'Brain metastasis mouse",
      "models' and Discussion)."
    )
  )

  ini({
    lkg_sens <- log(0.108)
    label("Proliferation rate constant of sensitive cells lambda_S (1/day)") # Table 2 Hypothesis 1: lambda_S = 0.108 (RSE 8.2%); printed unit '10^6 (p/s)/day' is a misprint, Results compare it with lambda = 0.0855 1/day
    lkg_res <- log(0.0556)
    label("Proliferation rate constant of resistant cells lambda_R (1/day)") # Table 2 Hypothesis 1: lambda_R = 0.0556 (RSE 49.1%)
    lkkill <- log(0.0139)
    label("Cell-kill rate constant of sensitive cells per unit daily dose kD (1/day per mg/kg)") # Table 2 Hypothesis 1: kD = 0.0139 (RSE 11.7%); multiplies D in Eq 2, Fig 5 confirms 0.21 1/day at 15 mg/kg
    lrbase <- log(137)
    label("Total tumour bioluminescence at randomisation (1e6 photons/s)") # Table 2 Hypothesis 1: Baseline = 137 x 10^6 p/s (RSE 19.6%)
    logitfsens <- logit(0.987)
    label("Logit of the fraction of cells sensitive at baseline S (unitless)") # Table 2 Hypothesis 1: Percentage sensitive at baseline = 98.7% (RSE 9.3%); p = exp(S)/(1+exp(S))

    # IIV: Methods define CV% as the square root of the variance, so omega^2 = (CV/100)^2
    etalkg_sens ~ 0.232324 # Table 2 Hypothesis 1: lambda_S IIV CV 48.2% -> 0.482^2
    etalkkill ~ 0.0121 # Table 2 Hypothesis 1: kD IIV CV 11.0% -> 0.110^2
    etalrbase ~ 6.0025 # Table 2 Hypothesis 1: Baseline IIV CV 245% -> 2.45^2
    etalogitfsens ~ 0.234256 # Table 2 Hypothesis 1: Percentage sensitive IIV CV 48.4% -> 0.484^2, log-normal on S

    # Residual error: log-normal (Methods), 19.0% read as the log-scale SD
    expSd <- 0.190
    label("Log-scale residual SD of bioluminescence (unitless)") # Table 2 Hypothesis 1: Residual error = 19.0%
  })
  model({
    kg_sens <- exp(lkg_sens + etalkg_sens)
    kg_res <- exp(lkg_res)
    kkill <- exp(lkkill + etalkkill)
    rbase <- exp(lrbase + etalrbase)
    # IIV on S is log-normal, as the Methods state for every random effect
    fsens <- expit(logitfsens * exp(etalogitfsens))

    dose_d <- DOSE_ZORIFERTINIB_MGKGD

    # Initial conditions: Sensitive(0) = p Total(0), Resistant(0) = (1 - p) Total(0)
    tumor_sensitive(0) <- fsens * rbase
    tumor_resistant(0) <- (1 - fsens) * rbase

    # Equations 2-3 with kR = 0 (Hypothesis 1: no conversion during the study)
    d/dt(tumor_sensitive) <- kg_sens * tumor_sensitive - kkill * dose_d * tumor_sensitive
    d/dt(tumor_resistant) <- kg_res * tumor_resistant

    # Equation 4
    tumor_size <- tumor_sensitive + tumor_resistant
    frac_resistant <- tumor_resistant / tumor_size

    tumor_size ~ lnorm(expSd)
  })
}
