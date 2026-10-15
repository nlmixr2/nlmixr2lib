Zou_2020_leuprorelin <- function() {
  description <- paste(
    "Population disease-progression model of serum prostate-specific antigen",
    "(PSA) in men with hormone-sensitive prostate cancer treated with the",
    "LHRH agonist leuprorelin, developed from a US medical claims database",
    "(Humana, 2007-2011; 264 subjects, 1113 PSA observations). The final",
    "'clonal selection' structure (the paper's Model II) splits the observed",
    "baseline PSA into a drug-resistant fraction R = exp(-RP) that grows",
    "first-order at GR from the baseline draw onwards, and a drug-sensitive",
    "fraction (1 - R) that grows first-order at GS until the first",
    "leuprorelin dose and is then killed first-order at DS. Covariates:",
    "hemoglobin (power, on RP), baseline PSA (power, on DS) and antiandrogen",
    "use within 30 days of leuprorelin initiation (exponential, on DS).",
    "There is no PK input; treatment enters only through the per-subject",
    "time of the first leuprorelin dose, T_SCAN_TO_DOSE. Additive residual",
    "error on log PSA."
  )
  reference <- paste(
    "Zou Y, Tang F, Talbert JC, Ng CM.",
    "Using medical claims database to develop a population disease",
    "progression model for leuprorelin-treated subjects with",
    "hormone-sensitive prostate cancer.",
    "PLoS ONE. 2020;15(3):e0230571.",
    "doi:10.1371/journal.pone.0230571. PMCID: PMC7092991.",
    "Structural equations from Methods Eqs 2-3 (Model II), IIV form from",
    "Eq 4, residual error from Eq 5, covariate forms from Eqs 16-17; all",
    "final parameter values from Table 2.",
    sep = " "
  )
  vignette <- "Zou_2020_leuprorelin"

  units <- list(
    time = "day",
    dosing = "n/a (no PK input; leuprorelin treatment enters only through the per-subject start time T_SCAN_TO_DOSE)",
    concentration = "ng/mL (the observable `PSA` is serum prostate-specific antigen)"
  )

  # Both ODE states are PSA sub-fractions of the measured serum PSA (ng/mL):
  # `growth` is the PSA produced by the drug-resistant clone (paper's PSA_R),
  # `shrink` the PSA produced by the drug-sensitive clone (paper's PSA_S).
  # The split is stated explicitly in Methods Eqs 2-3.
  compartmentData <- list(
    growth = list(analyte = "PSA", units = "ng/mL", specimen = "serum", verified = TRUE),
    shrink = list(analyte = "PSA", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    PSA_BL = list(
      description = "Observed baseline serum PSA, the last PSA measured before the first leuprorelin dose; the model's time origin is the date of this draw.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      source_name = "BAS",
      notes = paste(
        "Two roles, both from the paper. (a) It is a per-subject regressor, NOT an estimated parameter, and scales both sub-states (Methods Eqs 2-3): growth(0) = R * PSA_BL and shrink(0) = (1 - R) * PSA_BL, so PSA(0) = PSA_BL exactly.",
        "(b) Power covariate on the drug kill rate DS, centered at the cohort median 8.5 ng/mL (Results Eq 17; Table 2 'BAS on DS' = 0.174).",
        "Cohort (Table 1): median 8.50 ng/mL, range 0.200-782. Subjects whose baseline PSA was undetectable were excluded (Methods, Study subjects).",
        sep = " "
      )
    ),
    HGB = list(
      description = "Baseline blood hemoglobin concentration.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      source_name = "HGB",
      notes = paste(
        "Power covariate on RP, the log-scale parameter of the resistant fraction R = exp(-RP), centered at the cohort median 13.6 g/dL (Results Eq 16; Table 2 'HGB on RP' = 2.30).",
        "Lower hemoglobin gives a smaller RP and therefore a LARGER drug-resistant fraction: typical R = 9.36%, 1.94% and 0.326% at the cohort 5th percentile (10.9 g/dL), median (13.6) and 95th percentile (16.0) (Results).",
        "Cohort (Table 1): median 13.6 g/dL, range 6.80-17.4. Subjects with hemoglobin < 6 g/dL were excluded as likely acutely ill (Methods).",
        sep = " "
      )
    ),
    CONMED_ANTIANDROGEN = list(
      description = "Antiandrogen use indicator: 1 = an antiandrogen (bicalutamide, enzalutamide, flutamide or nilutamide) was dispensed within 30 days of leuprorelin initiation, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no antiandrogen within 30 days of leuprorelin initiation; 231 of 264 subjects)",
      source_name = "AND (IND_AND in Eq 17)",
      notes = paste(
        "Exponential covariate on DS: DS * exp(0.677 * CONMED_ANTIANDROGEN) (Results Eq 17; Table 2 'AND on DS' = 0.677), i.e. a 96.8% higher kill rate with antiandrogen use (Results).",
        "The paper presumes these short courses were given to prevent the testosterone flare of LHRH-agonist initiation (Methods). PSA observed after the start of CONTINUOUS antiandrogen therapy was excluded from the dataset, so this indicator describes a peri-initiation course only.",
        "The class composition (bicalutamide, enzalutamide, flutamide, nilutamide) is this paper's definition.",
        "Cohort (Table 1): 33 yes, 231 no.",
        sep = " "
      )
    ),
    T_SCAN_TO_DOSE = list(
      description = "Per-subject time from the baseline PSA draw (the model's time origin) to the first leuprorelin dose.",
      units = "day",
      type = "continuous",
      reference_category = NULL,
      source_name = "t1 (time of the first LHRH dose)",
      notes = paste(
        "Enters as the switch point of the drug-sensitive sub-state (Methods Eqs 2-3): t_s = min(t, t1) and t_k = max(0, t - t1). Before t1 the sensitive clone grows at GS; from t1 onwards it is killed at DS. The resistant clone grows at GR from time 0 regardless.",
        "Only the FIRST dose matters: once started, treatment is assumed to act continuously whether the subject was on continuous or intermittent leuprorelin, which the claims data could not distinguish (Discussion, limitations). Dose amount and dose intensity did not enter the model (dose intensity was tested and was not significant).",
        "Per-subject DATA from the claims pharmacy fill dates; the paper does not report its distribution.",
        sep = " "
      )
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Tested in the full covariate model (WAM-BE backward elimination) and not retained. Cohort median 80 years (range 60-100; Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "IU/L",
      type = "continuous",
      notes = "Tested and not retained. Cohort median 20 IU/L (range 9-91; Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "IU/L",
      type = "continuous",
      notes = "Tested and not retained. Cohort median 18 IU/L (range 4-110; Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested and not retained. Cohort median 1.10 mg/dL (range 0.700-9.30; Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase.",
      units = "IU/L",
      type = "continuous",
      notes = "Tested and not retained. Cohort median 76.5 IU/L (range 23.0-3640; Table 1)."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/dL",
      type = "continuous",
      notes = "Tested and not retained. Cohort median 4.13 g/dL (range 2.90-4.80; Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 264,
    n_observations = 1113,
    n_studies = 1,
    age_range = "60-100 years (median 80)",
    sex_female_pct = 0,
    race_ethnicity = c(Caucasian = 189, Black = 59, `Hispanic/other` = 16),
    regions = "United States (Humana commercially insured; South 196, Midwest 42, West 20, Northeast 6)",
    disease_state = "Malignant, hormone-sensitive prostate cancer (ICD-9-CM 185 / ICD-10-CM C61) on leuprorelin as the only LHRH agonist",
    dose_range = "Leuprorelin, continuous or intermittent, at any dose; dose intensity (dose received relative to an expected 7.5 mg per month) was not a significant covariate",
    baseline_psa = "median 8.50 ng/mL (range 0.200-782)",
    baseline_hgb = "median 13.6 g/dL (range 6.80-17.4)",
    antiandrogen_use = "33 of 264 within 30 days of leuprorelin initiation",
    notes = paste(
      "Retrospective medical claims data (Humana, 1 January 2007 to 31 December 2011). Inclusion: a PSA before leuprorelin initiation and at least one during treatment.",
      "Exclusions (86.4% of eligible patients in total): PSA falling before leuprorelin, undetectable baseline PSA or same-day duplicate measurements; undetectable PSA throughout treatment; incomplete demographics; hemoglobin < 6 g/dL.",
      "PSA observations after the start of continuous antiandrogen therapy, surgery, radiotherapy or chemotherapy (whichever first) were dropped.",
      "PSA below 0.1 ng/mL (LLOQ) was retained as censored (M3-type likelihood).",
      "Estimation: NONMEM 7.3, MCPEM; covariate selection by Wald's approximation with backward elimination (WAM-BE).",
      sep = " "
    )
  )

  ini({
    # Structural parameters: typical values at HGB = 13.6 g/dL,
    # PSA_BL = 8.5 ng/mL and no antiandrogen use. All from Table 2.
    lkse <- log(3.78e-2)
    label("Kill rate DS of the drug-sensitive PSA fraction after leuprorelin start (1/day)") # Table 2: DS = 3.78 x 10^-2 day^-1 (%CV 6.19)
    lkge_sens <- log(1.96e-3)
    label("Growth rate GS of the drug-sensitive PSA fraction before leuprorelin start (1/day)") # Table 2: GS = 1.96 x 10^-3 day^-1 (%CV 22.5)
    lrp <- log(3.94)
    label("RP, minus the log of the drug-resistant fraction R = exp(-RP) (unitless)") # Table 2: RP = 3.94 (%CV 7.44); R = exp(-3.94) = 1.94%
    lkge <- log(6.54e-4)
    label("Growth rate GR of the drug-resistant PSA fraction (1/day)") # Table 2: GR = 6.54 x 10^-4 day^-1 (%CV 28.4)

    # Covariate effects (Table 2; functional forms Eqs 16-17).
    e_hgb_rp <- 2.30
    label("Power exponent of hemoglobin (/13.6 g/dL) on RP (unitless)") # Table 2: theta HGB_RP = 2.30 (%CV 24.7)
    e_psa_bl_kse <- 0.174
    label("Power exponent of baseline PSA (/8.5 ng/mL) on DS (unitless)") # Table 2: theta BAS_DS = 0.174 (%CV 24.5)
    e_antiandrogen_kse <- 0.677
    label("Exponential coefficient of antiandrogen use on DS (unitless)") # Table 2: theta AND_DS = 0.677 (%CV 30.0)

    # IIV: exponential (Eq 4); Table 2 reports the variances omega^2.
    etalkse ~ 0.453 # Table 2: omega^2 DS = 0.453
    etalkge_sens ~ 2.59 # Table 2: omega^2 GS = 2.59
    etalrp ~ 0.944 # Table 2: omega^2 RP = 0.944
    etalkge ~ 3.76 # Table 2: omega^2 'DR' = 3.76 (the IIV of GR; the row label is a typo)

    # Residual error: additive on log PSA (Eq 5).
    expSd <- 0.201
    label("Additive residual error SD on log PSA (log ng/mL)") # Table 2: sigma add = 2.01 x 10^-1 (%CV 3.49)
  })
  model({
    # Individual parameters (Eqs 4, 16, 17)
    kse <- exp(lkse + etalkse) *
      (PSA_BL / 8.5)^e_psa_bl_kse *
      exp(e_antiandrogen_kse * CONMED_ANTIANDROGEN)
    kge_sens <- exp(lkge_sens + etalkge_sens)
    rp <- exp(lrp + etalrp) * (HGB / 13.6)^e_hgb_rp
    kge <- exp(lkge + etalkge)

    # Drug-resistant fraction of the baseline PSA (Methods: R = exp(-RP))
    fres <- exp(-rp)

    # Model II (Methods Eqs 2-3):
    #   PSA_R(t) = BAS * R       * exp(GR * t)
    #   PSA_S(t) = BAS * (1 - R) * exp(GS * t_s) * exp(-DS * t_k)
    # with t_s = min(t, t1), t_k = max(0, t - t1), t1 = first leuprorelin
    # dose. Written as two first-order ODEs whose sensitive-clone rate
    # switches from +GS to -DS at t1, which reproduces the closed form
    # exactly.
    dosed <- time >= T_SCAN_TO_DOSE

    d / dt(growth) <- kge * growth
    d / dt(shrink) <- (kge_sens * (1 - dosed) - kse * dosed) * shrink

    growth(0) <- fres * PSA_BL
    shrink(0) <- (1 - fres) * PSA_BL

    PSA <- growth + shrink
    PSA ~ lnorm(expSd)
  })
}
