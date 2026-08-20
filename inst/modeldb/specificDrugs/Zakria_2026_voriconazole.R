Zakria_2026_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order elimination for voriconazole in adult Pakistani cancer patients receiving therapeutic drug monitoring (Zakria 2026); a binary age-group indicator (> 65 years vs <= 65 years) is the only retained covariate and acts on clearance"
  reference <- "Zakria KZ, Usman M, Bilal H, Akbar Z, Hussain T, Ali M, Sattar A, Alvi I, Rasheed H, Zulfiqar S, Khan MR, Mushtaq MZ. Comparative pharmacokinetics of voriconazole between elderly and young adult patients: a population pharmacokinetic study. Journal of Pharmaceutical Policy and Practice. 2026;19(1):2601420. doi:10.1080/20523211.2025.2601420"
  vignette <- "Zakria_2026_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. `central` is the single disposition compartment of the
  # one-compartment model (Zakria 2026 Results, "Base model"); plasma is the
  # matrix because all observations are plasma TDM trough concentrations
  # (Zakria 2026 Methods, "Study design and data collection").
  compartmentData <- list(
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE_GT65 = list(
      description        = "Age greater than 65 years indicator (elderly age group)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (age <= 65 years; the 'young adult' group, n = 31 / 61% of the cohort)",
      notes              = "Additive-fractional effect on CL: CL = 6.46 * (1 + (-0.519) * AGE_GT65), i.e. 6.46 L/h when AGE_GT65 = 0 and 6.46 * (1 - 0.519) = 3.11 L/h when AGE_GT65 = 1 (Zakria 2026 Eqs 1 and 2, Table 2 rows CL-AGE_GROUP 1 / CL-AGE_GROUP 2). The threshold is strictly greater than 65 years: Zakria 2026 assigns subjects aged exactly 65 years to the reference (<= 65 y) group, which is why the canonical name is AGE_GT65 and not AGE_GE65. The underlying continuous AGE was also screened as a covariate on CL and was significant on forward inclusion (Supplementary Table S1, dOFV 7.99, p = 0.0047) but was dropped in favour of the dichotomised age group (dOFV 15.22, p = 0.000096); AGE is therefore recorded in covariatesDataExcluded, not here.",
      source_name        = "AGE_GROUP (1 = age <= 65 y, 2 = age > 65 y)"
    )
  )

  # Covariates that Zakria 2026 screened during stepwise covariate modelling
  # but did NOT retain in the final model (Zakria 2026 Supplementary Table S1).
  # Documentation only -- these are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL. Significant on forward inclusion (base OFV 153.18 -> 145.19, dOFV 7.99, p = 0.0047) but superseded by the dichotomised age group (dOFV 15.22), and not significant when tested on top of the age group (dOFV 0.75, p = 0.386). Not retained. Cohort median 53 years (range 18-77) per Zakria 2026 Table 1."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL; not significant (dOFV 1.08, p = 0.30 on forward inclusion; dOFV 0.20, p = 0.655 on top of the age group). Not retained. Cohort median 59 kg (range 42-92) per Zakria 2026 Table 1. The Discussion notes that no subject in this cohort was obese, so obesity was not tested."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on CL. Nominally significant on forward inclusion (dOFV 4.00, p = 0.045) but superseded by the age group and not significant on top of it (dOFV 1.17, p = 0.279). Not retained. Cohort median 0.65 mg/dL (range 0.18-4.21) per Zakria 2026 Table 1."
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL; not significant (dOFV 1.22, p = 0.271 on forward inclusion; dOFV 0.89, p = 0.345 on top of the age group). Not retained. Derived with the Cockcroft-Gault equation per Zakria 2026 Methods. Cohort median 100.9 mL/min (range 9.5-379.2) per Zakria 2026 Table 1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL. Entered the model during forward inclusion on top of the age group (dOFV 5.33, p = 0.0209 against alpha = 0.05) but was removed during backward elimination (dOFV 5.33, p = 0.021 against the stricter alpha = 0.01) and is absent from the final model (Zakria 2026 Supplementary Table S1). Cohort 24 female (47%) / 27 male (53%) per Zakria 2026 Table 1."
    ),
    TUMTP_LEUK = list(
      description = "Leukemia cancer-type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cancer type was screened as a covariate on CL per Zakria 2026 Methods and Results, and no cancer-type relation was retained. Supplementary Table S1 does not tabulate a CL-cancer-type step, so no dOFV is available for this covariate. Leukemia n = 26 (51%) per Zakria 2026 Table 1."
    ),
    TUMTP_LYMPH = list(
      description = "Lymphoma cancer-type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the cancer-type covariate on CL; not retained. Lymphoma n = 11 (21.5%) per Zakria 2026 Table 1."
    ),
    TUMTP_BREAST = list(
      description = "Breast-cancer indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the cancer-type covariate on CL; not retained. Breast cancer n = 12 (23.5%) per Zakria 2026 Table 1."
    ),
    TUMTP_SARC = list(
      description = "Sarcoma cancer-type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the cancer-type covariate on CL; not retained. Sarcoma n = 2 (4%) per Zakria 2026 Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    n_observations = 56L,
    age_range      = "18-77 years",
    age_median     = "53 years",
    weight_range   = "42-92 kg",
    weight_median  = "59 kg",
    sex_female_pct = 47,
    age_group      = c(le65y_n = 31, le65y_pct = 61, gt65y_n = 20, gt65y_pct = 39),
    cancer_type    = c(Leukemia_pct = 51, Lymphoma_pct = 21.5, BreastCancer_pct = 23.5, Sarcoma_pct = 4),
    renal_function = "Serum creatinine median 0.65 mg/dL (range 0.18-4.21); Cockcroft-Gault creatinine clearance median 100.9 mL/min (range 9.5-379.2).",
    disease_state  = "Adult cancer patients (leukemia, lymphoma, breast cancer, sarcoma) receiving voriconazole for serious systemic fungal infections and undergoing routine therapeutic drug monitoring.",
    dose_range     = "Recorded single doses 56-400 mg (median 190 mg). Zakria 2026 does not state the route or the dosing interval; the same hospital's voriconazole TDM protocol described by Akbar 2025 (a companion study from the same centre, with an identical 56-400 mg recorded dose range) is intravenous 6 mg/kg q12h loading for 24 h then 4 mg/kg q12h maintenance. See the vignette Assumptions section.",
    regions        = "Single center: Shaukat Khanum Memorial Cancer Hospital and Research Centre, Lahore, Pakistan.",
    notes          = "Retrospective, single-centre, open-label, non-interventional analysis of routine therapeutic drug monitoring (TDM) records. 56 samples from 51 patients (Zakria 2026 Table 1; the Results text states '53 data samples from 51 patients' while Table 1 reports 56 samples -- see vignette Errata). All samples were drawn at trough. Observed trough concentrations 0.5-8.7 mg/L (median 3.35). Estimation by FOCE-I in NONMEM with the PsN toolkit and Pirana. Final-model parameter estimates and the 1000-replicate bootstrap are in Zakria 2026 Table 2; the stepwise-covariate-modelling OFV trail is in Supplementary Table S1."
  )

  ini({
    # Structural parameters. The typical-value reference subject is a patient
    # aged <= 65 years (Zakria 2026 Eq 1: "IF (AGE <= 65) CL = 6.46").
    lcl <- log(6.46); label("Clearance in the age <= 65 years reference group (L/h)")  # Zakria 2026 Table 2 (final estimate) and Eq 1
    lvc <- log(97.1); label("Volume of distribution (L)")                              # Zakria 2026 Table 2 (final estimate); Results, "Covariate analysis"

    # Covariate effect on CL. Additive-fractional change relative to the
    # age <= 65 years reference: CL = 6.46 * (1 - 0.519) = 3.11 L/h for
    # subjects older than 65 years, matching the 3.11 L/h quoted in the
    # Abstract, Results and Discussion.
    e_agegt65_cl <- -0.519; label("Additive-fractional effect of age > 65 years vs age <= 65 years on CL (unitless)")  # Zakria 2026 Table 2 (row CL-AGE_GROUP 2) and Eq 2

    # IIV. Zakria 2026 reports interindividual variability on CL only; the
    # Limitations paragraph states "the interindividual variability was
    # evaluated only for the CL of VCZ due to the availability of trough
    # concentrations", so there is no eta on Vd. NONMEM exponential IIV is
    # `CL = TVCL * exp(eta)` with var(eta) = omega^2; the paper reports
    # IIV-CL as a CV% (23.7%). Using the standard NONMEM reporting
    # convention CV% ~= sqrt(omega^2) * 100 gives omega^2 = 0.237^2 =
    # 0.056169. (The exact log-normal variance log(1 + 0.237^2) = 0.054697
    # is an alternative reading; the approximate form is used here to match
    # the convention already used for the companion model
    # Akbar_2025_voriconazole from the same centre.)
    etalcl ~ 0.056169                        # Zakria 2026 Table 2 (IIV-CL 23.7%); see comment above for the omega^2 derivation

    # Residual error. Zakria 2026 selected a proportional error model
    # (Results, "Base model") and reports a single "Proportional error"
    # estimate of 0.519, taken here as the proportional SD (i.e. 51.9%).
    # See the vignette Errata for the scale ambiguity.
    propSd <- 0.519; label("Proportional residual error (fraction)")  # Zakria 2026 Table 2 (row Proportional error)
  })

  model({
    # Age-group effect on CL as an additive-fractional multiplier
    # (Zakria 2026 Eqs 1 and 2). AGE_GT65 = 0 recovers Eq 1 exactly.
    age_factor <- 1 + e_agegt65_cl * AGE_GT65

    # Individual parameters. Exponential IIV on CL only.
    cl <- exp(lcl + etalcl) * age_factor
    vc <- exp(lvc)

    # One-compartment model with first-order elimination (Zakria 2026
    # Results, "Base model"); linCmt() resolves to the analytical solution
    # given cl and vc, equivalent to NONMEM ADVAN1 TRANS2. No absorption
    # stage is parameterised because the paper reports no ka and no
    # bioavailability term.
    Cc <- linCmt()

    # Concentration units: dose mg / volume L = mg/L, matching the
    # "Concentration (mg/l)" units of Zakria 2026 Table 1.
    Cc ~ prop(propSd)
  })
}
