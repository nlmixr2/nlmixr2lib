Wattanakul_2026_paracetamol_alt <- function() {
  description <- "Exposure-response model for the natural logarithm of serum alanine aminotransferase (ALT) over 0-168 h in Malaysian patients with Plasmodium knowlesi malaria randomised to oral paracetamol 1 g every 6 h for 72 h or no paracetamol (PACKNOW trial; Wattanakul 2026). ln(ALT) rises linearly in time from an intercept of 3.49 (scaled linearly by body weight) with a slope of 0.00236 1/h that is scaled linearly by the individual paracetamol AUC0-72h (+0.3% per mg*h/L), age (-1.1% per year) and baseline ALT (-1.2% per U/L), each centred on the printed value. Additive IIV on slope and intercept; additive residual error on ln(ALT). Paracetamol exposure enters only as the per-subject AUC_PARACETAMOL covariate, computed from the companion PK model modellib('Wattanakul_2026_paracetamol')."
  reference <- paste(
    "Wattanakul T, Cooper DJ, Plewes K, Grigg MJ, Rajahram GS, William T,",
    "Dondorp AM, Edstein MD, Birrell GW, Anstey NM, Hoglund RM, Barber BE,",
    "Tarning J. Population Pharmacokinetics and Pharmacodynamics of",
    "Paracetamol in Malaysian Patients With Plasmodium knowlesi Malaria.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70283.",
    "doi:10.1002/psp4.70283.",
    "Exposure driver from the companion PK model; see",
    "modellib('Wattanakul_2026_paracetamol').",
    sep = " "
  )
  vignette <- "Wattanakul_2026_paracetamol"
  units <- list(
    time = "h",
    dosing = "n/a (no dosing events; paracetamol exposure is supplied as the AUC_PARACETAMOL covariate in mg*h/L)",
    concentration = "log(U/L) (lnalt is the natural logarithm of serum ALT in U/L)"
  )

  covariateData <- list(
    AUC_PARACETAMOL = list(
      description = "Individual model-predicted paracetamol AUC from 0 to 72 h after the first dose.",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject. Derived in the source from each patient's empirical-Bayes parameters in the final PK model (modellib('Wattanakul_2026_paracetamol')); 0 for patients who received no paracetamol. Enters linearly, centred at 114.76 mg*h/L (supplement ALT stream 'SLOPEAUC72MG = (1 + THETA(5)*(AUC72MG - 114.76))'), so the untreated arm sits below the centring value. Paracetamol-arm median (IQR) 408 (331-504) mg*h/L (Table 2).",
      source_name = "AUC72MG"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed enrollment value. Linear effect on the ln(ALT) slope centred at 36 years (supplement ALT stream 'SLOPEAGE = (1 + THETA(4)*(AGE - 36))').",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed enrollment value. Linear effect on the ln(ALT) intercept centred at 60 kg (supplement ALT stream 'INTERCEPTWT = (1 + THETA(3)*(WT - 60))').",
      source_name = "WT"
    ),
    ALT = list(
      description = "Baseline (enrollment) serum alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed enrollment value, NOT the time-varying observation (the observation is lnalt). Linear effect on the ln(ALT) slope centred at 34 U/L (supplement ALT stream 'SLOPEBALT = (1 + THETA(6)*(BALT - 34))'). The stream sets the factor to 1 when baseline ALT is missing (BALT = -99); supply ALT = 34 for a patient without a baseline value to reproduce that. Paracetamol-arm median (IQR) 35 (23-55) U/L, control 33 (22-48) U/L (Table 1).",
      source_name = "BALT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 372L,
    n_studies = 1L,
    n_observations = "ALT measured at admission, day 3, day 7 and day 28 (plus clinician-ordered extras); the 0-168 h window was modelled.",
    age_range = "5 years and older (inclusion criterion); IQR 26-49 (paracetamol) and 24-48 (control) years",
    age_median = "35 (paracetamol) and 36 (control) years",
    weight_range = "IQR 53.5-68.4 (paracetamol) and 52.0-66.5 (control) kg",
    weight_median = "61.0 (paracetamol) and 58.0 (control) kg",
    sex_female_pct = 15.6,
    race_ethnicity = "Malaysian (Sabah); not further broken down.",
    disease_state = "Hospitalised patients with PCR-confirmed Plasmodium knowlesi monoinfection; 7.5% WHO-defined severe malaria and 28.5% with acute kidney injury at enrollment. Patients with known cirrhosis or consuming more than six standard alcoholic drinks per day were excluded.",
    dose_range = "Paracetamol 1 g orally every 6 h for 72 h (n = 183) or no paracetamol (n = 189).",
    regions = "Sabah, Malaysia.",
    notes = "PACKNOW trial (NCT03056391), both arms. Demographics from Wattanakul 2026 Table 1; ALT sampling schedule from the supplement S1."
  )

  ini({
    # Wattanakul 2026 Table 4. Covariate effects are printed as percent per
    # unit and are converted to fractions here (0.3% -> 0.003). The supplement
    # ALT stream carries the centring values but its $THETA block holds
    # earlier initial estimates, not the final values.
    slope_lnalt <- 0.00236
    label("Slope of ln(ALT) over time at the centring covariate values (1/h)") # Table 4: Slope = 0.00236 1/h (RSE 11.4%)
    bl_lnalt <- 3.49
    label("Intercept of ln(ALT) at time zero and WT = 60 kg (log U/L)") # Table 4: Intercept = 3.49 (RSE 0.9%)
    e_age_slope_lnalt <- -0.011
    label("Linear effect of age on the ln(ALT) slope (fraction per year from 36 years)") # Table 4: theta_AGE on slope = -1.1% per year (RSE 36.1%)
    e_auc_paracetamol_slope_lnalt <- 0.003
    label("Linear effect of paracetamol AUC0-72h on the ln(ALT) slope (fraction per mg*h/L from 114.76 mg*h/L)") # Table 4: theta_AUC0-72H on slope = 0.3% per mg*h/L (RSE 26.0%)
    e_alt_slope_lnalt <- -0.012
    label("Linear effect of baseline ALT on the ln(ALT) slope (fraction per U/L from 34 U/L)") # Table 4: theta_BASEALT on slope = -1.2% per U/L (RSE 19.2%)
    e_wt_bl_lnalt <- 0.004
    label("Linear effect of body weight on the ln(ALT) intercept (fraction per kg from 60 kg)") # Table 4: theta_WT on intercept = 0.4% per kg (RSE 23.6%)

    # Additive IIV. Table 4 reports these as standard deviations: the
    # stream's $OMEGA initial values (7.17E-06, 0.292) are the squares of
    # numbers within 2% of the Table 4 entries.
    etaslope_lnalt ~ 7.5076e-06 # Table 4: IIV Slope = 0.00274 (SD); 0.00274^2
    etabl_lnalt ~ 0.299209 # Table 4: IIV Intercept = 0.547 (SD); 0.547^2

    # Table 4 sigma = 0.118 is read as a VARIANCE: the stream '$SIGMA 0.12'
    # is within 2% of it, whereas reading it as an SD would put the final
    # variance (0.0139) nine-fold below that initial value while every
    # $OMEGA moved by under 2%.
    addSd <- sqrt(0.118)
    label("Additive residual SD on ln(ALT) (log U/L)") # Table 4: sigma = 0.118 (variance; RSE 10.3%)
  })

  model({
    slope_cov <- (1 + e_age_slope_lnalt * (AGE - 36)) *
      (1 + e_auc_paracetamol_slope_lnalt * (AUC_PARACETAMOL - 114.76)) *
      (1 + e_alt_slope_lnalt * (ALT - 34))
    slope_lnalt_i <- slope_lnalt * slope_cov + etaslope_lnalt
    bl_lnalt_i <- bl_lnalt * (1 + e_wt_bl_lnalt * (WT - 60)) + etabl_lnalt

    # Methods Eq. 10: ln(ALT) = Base + Slope * t
    lnalt <- bl_lnalt_i + slope_lnalt_i * t
    lnalt ~ add(addSd)
  })
}
