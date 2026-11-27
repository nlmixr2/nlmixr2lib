Martin_2018_zorifertinib_brainmet_h12_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, PC9_Luc NSCLC brain metastasis model). Alternative",
    "(not preferred) tumour growth and resistance model for zorifertinib",
    "(AZD3759) fitted to bioluminescence, combining Hypotheses 1 and 2:",
    "a fraction of cells resistant at baseline plus first-order conversion of",
    "sensitive to resistant cells during the study, with a cell-kill term on",
    "the sensitive cells proportional to the daily dose in mg/kg (no PK",
    "layer). The authors preferred Hypothesis 1 (see",
    "Martin_2018_zorifertinib_brainmet_mouse). Martin 2018, Equations 2-4 and",
    "Table 2."
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
      "Same data as Martin_2018_zorifertinib_brainmet_mouse. AIC 4028.59 for",
      "this model against 4029.19 for Hypothesis 1; the authors preferred",
      "Hypothesis 1 on AIC, model fit and the negligible contribution of",
      "conversion (Results)."
    )
  )

  ini({
    lkg_sens <- log(0.0992)
    label("Proliferation rate constant of sensitive cells lambda_S (1/day)") # Table 2 Hypothesis 1 + 2: lambda_S = 0.0992 (RSE 11.7%); printed unit '10^6 (p/s)/day' is a misprint
    lkg_res <- log(0.0708)
    label("Proliferation rate constant of resistant cells lambda_R (1/day)") # Table 2 Hypothesis 1 + 2: lambda_R = 0.0708 (RSE 49.3%)
    lkr <- log(1.54e-6)
    label("Conversion rate constant of sensitive to resistant cells kR (1/day)") # Table 2 Hypothesis 1 + 2: kR = 1.54 x 10^-6 1/day (RSE 139%)
    lkkill <- log(0.0132)
    label("Cell-kill rate constant of sensitive cells per unit daily dose kD (1/day per mg/kg)") # Table 2 Hypothesis 1 + 2: kD = 0.0132 (RSE 19.0%); multiplies D in Eq 2
    lrbase <- log(177)
    label("Total tumour bioluminescence at randomisation (1e6 photons/s)") # Table 2 Hypothesis 1 + 2: Baseline = 177 x 10^6 p/s (RSE 34.0%)
    logitfsens <- logit(0.996)
    label("Logit of the fraction of cells sensitive at baseline S (unitless)") # Table 2 Hypothesis 1 + 2: Percentage sensitive at baseline = 99.6% (RSE 18.9%)

    # IIV: Methods define CV% as the square root of the variance, so omega^2 = (CV/100)^2
    etalkg_sens ~ 0.239121 # Table 2 Hypothesis 1 + 2: lambda_S IIV CV 48.9% -> 0.489^2
    etalkkill ~ 0.00824464 # Table 2 Hypothesis 1 + 2: kD IIV CV 9.08% -> 0.0908^2
    etalrbase ~ 6.1504 # Table 2 Hypothesis 1 + 2: Baseline IIV CV 248% -> 2.48^2
    etalogitfsens ~ 1.8225 # Table 2 Hypothesis 1 + 2: Percentage sensitive IIV CV 135% -> 1.35^2, log-normal on S

    # Residual error: log-normal (Methods), 16.5% read as the log-scale SD
    expSd <- 0.165
    label("Log-scale residual SD of bioluminescence (unitless)") # Table 2 Hypothesis 1 + 2: Residual error = 16.5%
  })
  model({
    kg_sens <- exp(lkg_sens + etalkg_sens)
    kg_res <- exp(lkg_res)
    kr <- exp(lkr)
    kkill <- exp(lkkill + etalkkill)
    rbase <- exp(lrbase + etalrbase)
    # IIV on S is log-normal, as the Methods state for every random effect
    fsens <- expit(logitfsens * exp(etalogitfsens))

    dose_d <- DOSE_ZORIFERTINIB_MGKGD

    # Initial conditions: Sensitive(0) = p Total(0), Resistant(0) = (1 - p) Total(0)
    tumor_sensitive(0) <- fsens * rbase
    tumor_resistant(0) <- (1 - fsens) * rbase

    # Equations 2-3
    d/dt(tumor_sensitive) <- kg_sens * tumor_sensitive - kr * tumor_sensitive - kkill * dose_d * tumor_sensitive
    d/dt(tumor_resistant) <- kg_res * tumor_resistant + kr * tumor_sensitive

    # Equation 4
    tumor_size <- tumor_sensitive + tumor_resistant
    frac_resistant <- tumor_resistant / tumor_size

    tumor_size ~ lnorm(expSd)
  })
}
