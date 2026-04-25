Chen_2022_guselkumab <- function() {
  description <- "One-compartment population PK model with first-order SC absorption for guselkumab (anti-IL-23 human IgG1 lambda mAb) in patients with active psoriatic arthritis (DISCOVER-1 and DISCOVER-2 phase 3 trials; Chen 2022)"
  reference <- "Chen Y, Miao X, Hsu CH, Zhuang Y, Kollmeier A, Xu Z, Zhou H, Sharma A. Population pharmacokinetics and exposure-response modeling analyses of guselkumab in patients with psoriatic arthritis. Clin Transl Sci. 2022;15(3):749-760. doi:10.1111/cts.13197"
  vignette <- "Chen_2022_guselkumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F and V/F per Chen 2022 Table 1 footnotes f and g; reference value 84 kg (median of pooled DISCOVER-1 + DISCOVER-2 PK analysis population). Time-fixed at baseline.",
      source_name        = "BWT"
    ),
    DIAB = list(
      description        = "Diabetes mellitus comorbidity at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes mellitus comorbidity)",
      notes              = "Multiplicative effect on CL/F per Chen 2022 Table 1 footnote f: CL/F = 0.596 * (BWT/84)^0.926 * 1.15^DIAB. With DIAB=1, CL/F is 15% higher than the non-diabetic reference (Results section, page 752, and Table 1). Diabetes-mellitus-positive comorbidity prevalence in the analysis population was approximately 9% (Results section).",
      source_name        = "DIAB"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Adults with active psoriatic arthritis (specific age range in Table S1 supplement, not extractable from main-text trim).",
    weight_median  = "84 kg (25th-75th percentile 71.0-97.3 kg, per Results page 752)",
    weight_range   = "25th-75th percentile 71.0-97.3 kg (Results, page 752); full range in Table S1 supplement.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Majority White (Table S1 supplement; specific percentages not in main-text trim).",
    disease_state  = "Active psoriatic arthritis (PsA) with inadequate response to standard nonbiologic disease-modifying therapy and/or anti-tumor necrosis factor alpha therapy. Pooled data from DISCOVER-1 (NCT03162796) and DISCOVER-2 (NCT03158285) phase 3 trials.",
    dose_range     = "Subcutaneous guselkumab 100 mg q4w or 100 mg at weeks 0 and 4 then q8w (the two phase 3 regimens). Sparse PK sampling.",
    regions        = "Global (multicenter phase 3 trials).",
    diabetes_pct   = 9.0,
    prior_anti_tnf_pct = 11.0,
    ada_positive_pct   = 2.0,
    notes          = "Demographics from Chen 2022 Table S1 (supplement). Median baseline weight 84 kg, ~9% with diabetes comorbidity, majority White, ~11% had prior anti-TNF-alpha therapy, 2.0% positive for antidrug antibodies. Median baseline PASI = 5.8 and DAS28 = 5.1 in the exposure-response analysis datasets (Results, page 752). The longitudinal exposure-response analysis used pooled data through week 24 from 1116 patients (Results, page 754); the size of the population PK analysis dataset itself is reported in Table S1 of the supplement, not in the main-text trim. Estimated terminal half-life for the typical reference subject is approximately 18.1 days (Results, page 752)."
  )

  ini({
    # Structural parameters from Chen 2022 Table 1 (final population PK model).
    # Typical values are for the reference subject: 84 kg body weight (median),
    # no diabetes comorbidity (DIAB = 0).
    lcl <- log(0.596); label("Apparent clearance at reference covariates (CL/F, L/day)")               # Chen 2022 Table 1 (CL/F = 0.596 L/day for 84 kg, non-diabetic)
    lvc <- log(15.5);  label("Apparent central volume of distribution at reference covariates (V/F, L)") # Chen 2022 Table 1 (V/F = 15.5 L for 84 kg)
    lka <- log(0.572); label("First-order SC absorption rate constant (Ka, 1/day)")                    # Chen 2022 Table 1 (Ka = 0.572 1/day)

    # Covariate effect exponents and multipliers (Chen 2022 Table 1 / Table 1 footnotes f-g):
    #   CL/F = 0.596 * (BWT/84)^0.926 * 1.15^DIAB
    #   V/F  = 15.5  * (BWT/84)^0.861
    e_wt_cl   <- 0.926; label("Power exponent of body weight on CL/F (unitless)")             # Chen 2022 Table 1 (BWT on CL/F)
    e_wt_vc   <- 0.861; label("Power exponent of body weight on V/F (unitless)")              # Chen 2022 Table 1 (BWT on V/F)
    e_diab_cl <- 1.15;  label("Multiplier on CL/F for diabetes comorbidity (1.15^DIAB)")      # Chen 2022 Table 1 (Diabetes on CL/F = 1.15)

    # IIV. Table 1 reports IIV as %CV (linear scale) for log-normal IIV
    # entries; convert with omega^2 = log(CV^2 + 1).
    #   CL/F 38.9% -> omega^2 = log(1 + 0.389^2) = 0.14092
    #   V/F  33.3% -> omega^2 = log(1 + 0.333^2) = 0.10515
    #   Ka   93.4% -> omega^2 = log(1 + 0.934^2) = 0.62725
    # Correlation between IIV of CL/F and V/F is 0.101 (Chen 2022 Table 1):
    #   covariance = 0.101 * sqrt(0.14092 * 0.10515) = 0.012295
    etalcl + etalvc ~ c(0.14092,
                        0.012295, 0.10515)                                               # Chen 2022 Table 1 (IIV CL/F 38.9%, IIV V/F 33.3%, correlation 0.101)
    etalka ~ 0.62725                                                                     # Chen 2022 Table 1 (IIV Ka 93.4%; shrinkage 61.7%)

    # Residual error: combined additive and proportional (Chen 2022 Methods, Base
    # model section, and Table 1).
    #   Proportional CV% = 19.1  -> propSd = 0.191
    #   Additive (ug/mL) = 0.00289
    propSd <- 0.191;   label("Proportional residual error (SD, fraction)")               # Chen 2022 Table 1 (Proportional residual error CV% = 19.1)
    addSd  <- 0.00289; label("Additive residual error (SD, ug/mL)")                      # Chen 2022 Table 1 (Additive residual error 0.00289 ug/mL)
  })
  model({
    # Individual PK parameters. Reference subject: 84 kg body weight, no
    # diabetes comorbidity. Covariate forms per Chen 2022 Table 1
    # footnotes f (CL/F) and g (V/F):
    #   CL/F_i = 0.596 * (BWT/84)^0.926 * 1.15^DIAB * exp(eta_CL)
    #   V/F_i  = 15.5  * (BWT/84)^0.861             * exp(eta_V)
    cl <- exp(lcl + etalcl) * (WT / 84)^e_wt_cl * e_diab_cl^DIAB
    vc <- exp(lvc + etalvc) * (WT / 84)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    # One-compartment linear PK with first-order SC absorption and first-order
    # elimination (Chen 2022 base model description; Results, page 751-752).
    # Bioavailability is implicit in the apparent clearance and volume
    # (CL/F, V/F); SC dosing is administered into the depot compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
