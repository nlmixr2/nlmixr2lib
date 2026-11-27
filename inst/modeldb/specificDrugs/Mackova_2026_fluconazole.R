Mackova_2026_fluconazole <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for",
    "intravenous fluconazole in hospitalised adults (82.8% in intensive care",
    "units, most with Candida infection), developed from routine therapeutic",
    "drug monitoring data at a single Czech centre (183 serum concentrations",
    "from 116 patients, 2022-2025; patients on renal replacement therapy",
    "excluded). Clearance is scaled by a power function of the time-varying",
    "CKD-EPI 2021 creatinine-based estimated glomerular filtration rate",
    "normalised to 90 mL/min, CL = 0.81 * (eGFR/90)^0.43 L/h, and the volume",
    "of distribution by a power function of adjusted body weight normalised",
    "to 75 kg, Vd = 41.58 * (ABW/75)^2.28 L. Inter-individual variability is",
    "log-normal on Vd, on CL and on the eGFR power exponent itself (eGFR",
    "entered the Monolix model as a regressor, which makes the exponent an",
    "individual parameter carrying its own random effect). Residual",
    "variability is proportional. The paper converts the clearance-eGFR",
    "relationship into a once-daily maintenance-dose nomogram,",
    "MD = 769.5 * (eGFR/90)^0.43 mg/day, preceded by a loading dose of 2.8",
    "times the maintenance dose."
  )
  reference <- paste(
    "Mackova T, Vodicka M, Slanar O, Sima M.",
    "Population Pharmacokinetic Modeling of Intravenous Fluconazole to Inform",
    "eGFR-Based Dosing: A Simulation Study Demonstrating Superiority over",
    "Standard Regimens. Pharmaceutics. 2026;18(8):929.",
    "doi:10.3390/pharmaceutics18080929."
  )
  vignette <- "Mackova_2026_fluconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "fluconazole",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate calculated with the CKD-EPI",
        "2021 creatinine equation. Time-varying: entered the Monolix model as",
        "a regressor, so each concentration record carries the eGFR",
        "contemporaneous with that sample rather than a single baseline value."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 90 mL/min ('an individual with normal renal function",
        "(eGFR = 90 mL/min)', Mackova 2026 Results 3.2); the model applies",
        "(CRCL/90)^e_crcl_cl to clearance (Table 2, Fixed effects block).",
        "SIZE NORMALISATION: Table 1, the Table 2 footnote and the Figure 3",
        "axis label all give 'mL/min', but the CKD-EPI 2021 equation natively",
        "returns mL/min/1.73 m^2, Results 3.1 quotes augmented renal clearance",
        "as eGFR >= 130 mL/min/1.73 m2, and the paper never states that it",
        "de-indexed for body surface area. Supply the CKD-EPI value as",
        "reported; the cohort the model was fitted to had median 86 (IQR",
        "55-103, range 8-160) (Table 1). Cystatin C-based CKD-EPI 2012 eGFR",
        "was tested in the 53 patients with cystatin C and performed",
        "equivalently, but runs much lower in this cohort (median 46 vs 90",
        "mL/min in that subgroup) and must not be substituted for this column."
      ),
      source_name = "eGFR"
    ),
    ABW = list(
      description = paste(
        "Adjusted body weight, ABW = IBW + 0.4 * (TBW - IBW), with ideal body",
        "weight from the Devine formula (50 kg + 2.3 kg per inch over 5 ft in",
        "men; 45.5 kg + 2.3 kg per inch over 5 ft in women)."
      ),
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 75 kg, the cohort median (Mackova 2026 Table 1:",
        "median 75, IQR 68-84, range 43-110 kg); enters Vd as",
        "(ABW/75)^2.28 (Table 2, Fixed effects block). The exponent is steep:",
        "across the observed 43-110 kg range it scales Vd by 0.28- to",
        "2.4-fold, so the model should not be extrapolated outside that range.",
        "Formula from Mackova 2026 Methods 2.2."
      ),
      source_name = "ABW"
    )
  )

  # Screened in the covariate analysis (Methods 2.3) but not retained:
  # "After incorporating these covariates into the model, none of the
  # remaining tested variables demonstrated a meaningful influence on
  # fluconazole pharmacokinetics" (Results 3.2). Documentation only; these
  # names are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Median 67 years (IQR 58-74, range 21-94; Table 1). Screened as a continuous covariate, not retained."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Median 175 cm (IQR 169-180, range 149-192; Table 1). Screened, not retained."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Median 85 kg (IQR 72-101, range 43-178; Table 1). Screened against adjusted body weight on Vd; ABW was retained instead."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Median 2.00 m^2 (IQR 1.84-2.20, range 1.36-2.75; Table 1). Screened, not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "38 of 116 patients female (32.8%; Table 1). The only categorical covariate tested; not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Median 70 umol/L (IQR 53-112, range 29-483; Table 1). Tested as a time-varying regressor and rejected in favour of creatinine-based eGFR."
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      notes = "Available in 53 patients only (median 1.41 mg/L, IQR 1.05-1.98, range 0.66-5.20; Table 1). Tested as a time-varying regressor, as was the cystatin C-based CKD-EPI 2012 eGFR; their predictive performance was equivalent to creatinine-based eGFR (Supplementary Figure S1), and creatinine-based eGFR was preferred because it was available for the whole cohort."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 116L,
    n_studies = 1L,
    n_observations = 183L,
    age_range = "21-94 years",
    age_median = "67 years",
    weight_range = "43-178 kg",
    weight_median = "85 kg",
    adjusted_weight = "median 75 kg (IQR 68-84, range 43-110)",
    height_range = "149-192 cm",
    sex_female_pct = 32.8,
    race_ethnicity = "Single-centre Czech cohort; race/ethnicity not reported in the source.",
    disease_state = paste(
      "Hospitalised adults treated with intravenous fluconazole: culture-proven",
      "Candida infection in most (C. albicans n = 75, C. glabrata n = 5,",
      "C. tropicalis n = 5, C. parapsilosis n = 2, mixed n = 6) and empirical",
      "therapy in 14 (12.1%). 96 of 116 (82.8%) were in intensive care units,",
      "mainly anaesthesiology/resuscitation (n = 78)."
    ),
    renal_function = paste(
      "Creatinine-based CKD-EPI 2021 eGFR median 86 mL/min (IQR 55-103, range",
      "8-160); only 4 patients (3.4%) had eGFR >= 130 mL/min/1.73 m^2.",
      "Patients on renal replacement therapy were EXCLUDED."
    ),
    dose_range = paste(
      "Median 400 mg (range 200-800 mg) intravenous fluconazole every 12 or",
      "24 h as a 30-, 60- or 120-min infusion, chosen in routine practice."
    ),
    regions = "Czech Republic (single centre: Tomas Bata Hospital Zlin).",
    notes = paste(
      "Retrospective single-centre TDM study, 1 January 2022 to 31 March 2025.",
      "Demographics in Mackova 2026 Table 1. 183 serum concentrations (1-4 per",
      "patient, mean 1.6): 30 troughs, 9 peaks and 144 mid-interval samples",
      "with exact sampling times. Total (not unbound) fluconazole was",
      "measured. Estimation in Monolix 2021R2 (SAEM); simulations in Simulx",
      "2021R2."
    )
  )

  ini({
    # ---- Structural parameters (Mackova 2026 Table 2, Fixed effects) ----
    # Both are typical values AT the reference covariates (ABW 75 kg, eGFR
    # 90 mL/min). Results 3.2: 'an individual with normal renal function
    # (eGFR = 90 mL/min) and a median ABW of 75 kg had a predicted Vd of
    # 41.58 L and a CL of 0.81 L/h, corresponding to an elimination half-life
    # of 35.6 h' (0.693 * 41.58 / 0.81 = 35.6 h).
    lvc <- log(41.58)
    label("Volume of distribution at ABW 75 kg (L)") # Table 2: Vd_pop = 41.58 L (R.S.E. 6.98%; bootstrap 42.12, 41.55-42.64)
    lcl <- log(0.81)
    label("Clearance at eGFR 90 mL/min (L/h)") # Table 2: CL_pop = 0.81 L/h (R.S.E. 4.77%; bootstrap 0.82, 0.82-0.83)

    # ---- Covariate effects ----
    e_abw_vc <- 2.28
    label("Power exponent on (ABW/75) for Vd (unitless)") # Table 2: beta_Vd_ABW = 2.28 (R.S.E. 16.0%; bootstrap 2.24, 2.21-2.30)

    # Power exponent on (eGFR / 90) for CL, log-transformed so that its
    # log-normal IIV is mu-referenced; the canonical untransformed name
    # e_crcl_cl is recovered in model(). Same encoding as the same group's
    # Subrtova_2026_daptomycin.R.
    le_crcl_cl <- log(0.43)
    label("Log power exponent on (CRCL/90) for CL (unitless)") # Table 2: beta_CL_eGFR = 0.43 (R.S.E. 30.3%; bootstrap 0.43, 0.41-0.45)

    # ---- Inter-individual variability (Table 2, 'Between subject
    # variability') ----
    # Monolix reports omega, the SD of the normally distributed random effect
    # on the log scale; nlmixr2 takes the variance, so each value is squared.
    etalvc ~ 0.1936 # Table 2: BSV Vd = 0.44 (R.S.E. 12.5%); 0.44^2
    etalcl ~ 0.0625 # Table 2: BSV CL = 0.25 (R.S.E. 19.0%); 0.25^2
    # IIV on the eGFR exponent. Methods 2.3 state that the time-varying
    # covariates 'were tested by incorporating them as regressors' (a Monolix
    # regressor's coefficient has to be declared as an individual parameter)
    # and that 'All model parameters were assumed to follow log-normal
    # distributions', which fixes the distribution as log-normal around
    # 0.43. Its eta-shrinkage is reported (0.27%, Results 3.2), confirming a
    # genuine estimated random effect.
    etale_crcl_cl ~ 1.5376 # Table 2: BSV beta_CL_eGFR = 1.24 (R.S.E. 19.8%); 1.24^2

    # ---- Residual variability ----
    # Monolix proportional error y = f * (1 + b * eps), i.e. nlmixr2 prop().
    propSd <- 0.18
    label("Proportional residual error (fraction)") # Table 2: Error model parameter, Proportional = 0.18 (R.S.E. 10.2%)
  })

  model({
    # 1. Covariate-effect exponent, recovered from its log-scale,
    #    mu-referenced counterpart.
    e_crcl_cl <- exp(le_crcl_cl + etale_crcl_cl)

    # 2. Individual parameters. Mackova 2026 Table 2, Fixed effects block:
    #      Vd = Vd_pop * (ABW / 75)^beta_Vd_ABW
    #      CL = CL_pop * (eGFR / 90)^beta_CL_eGFR
    #    CRCL is time-varying, so cl is re-evaluated at every record.
    vc <- exp(lvc + etalvc) * (ABW / 75)^e_abw_vc
    cl <- exp(lcl + etalcl) * (CRCL / 90)^e_crcl_cl

    # 3. Micro-constant. Written as an ODE rather than linCmt() because a
    #    time-varying regressor on CL makes kel time-varying.
    kel <- cl / vc

    # 4. Intravenous infusion only (30-120 min); no absorption compartment.
    d/dt(central) <- -kel * central

    # 5. Total fluconazole serum concentration (mg/L): dose mg / volume L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
