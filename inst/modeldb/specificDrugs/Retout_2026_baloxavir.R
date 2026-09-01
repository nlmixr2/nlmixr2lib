Retout_2026_baloxavir <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for baloxavir (the active form of the prodrug baloxavir marboxil) in influenza patients aged 1 year and older, with bodyweight, race (Asian vs non-Asian), sex and age covariates (Retout 2026)"
  reference <- "Retout S, De Buck S, Gaudreault J, Jolivet S, Duval V, Cosson V, Delporte ML. Population Pharmacokinetic and Exposure-Efficacy Analysis of Baloxavir Marboxil for Influenza Treatment and Post-Exposure Prophylaxis in Children. Clin Pharmacol Ther. 2026;119(4):1047-1056. doi:10.1002/cpt.70204"
  vignette <- "Retout_2026_baloxavir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Doses are administered as baloxavir marboxil (the prodrug), which is
  # converted pre-systemically and near-completely to baloxavir acid
  # ("baloxavir"); baloxavir marboxil itself was below the limit of
  # quantification in most samples (Retout 2026 Results, PopPK). CL/F and V/F
  # are therefore apparent parameters for baloxavir referenced to the
  # administered baloxavir marboxil dose, so compartment amounts are expressed
  # in mg of baloxavir marboxil dose equivalents.
  compartmentData <- list(
    depot       = list(analyte = "baloxavir marboxil", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "baloxavir", units = "mg (baloxavir marboxil dose equivalents)", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "baloxavir", units = "mg (baloxavir marboxil dose equivalents)", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric) scaling centered on a 70 kg reference weight; a single exponent is shared by CL/F and Q/F and a second exponent is shared by Vc/F and Vp/F (Retout 2026 Table 1 footnote d). Pooled dataset median 65.7 kg, range 4-217 kg (Table S2).",
      source_name        = "bodyweight"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Bodyweight-independent fractional reduction of CL/F, Vc/F and Q/F in Asian patients relative to non-Asian patients; Vp/F carries no race effect. Retout 2026 Table 1 footnote d defines Asian = 1 for Asian patients and 0 for non-Asian patients. Note that the demographic listing in Table S3 prints the opposite index labels ('0: Asian', '1: Non-Asian'); the model equations, the direction of the covariate estimates and the Discussion (CL/F ~50% lower in Asian patients) all pin the model indicator to Asian = 1.",
      source_name        = "Asian"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Fractional reduction of ka in female patients (Retout 2026 Table 1 footnote d, Sex = 1 for female patients).",
      source_name        = "Sex"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on ka centered on 37 years, the median age of the pooled popPK dataset (Retout 2026 Table 1 footnote d; Table S2 median age 37.0 years).",
      source_name        = "age"
    ),
    FORM_BXM_TAB10 = list(
      description        = "Baloxavir marboxil 10 mg tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other formulations used in the analysis: 20 mg tablet, 2% granules, granules for oral suspension 2 mg/mL)",
      notes              = "Relative bioavailability was set to 1.00 for every formulation in the analysis except the 10 mg tablet used in study T0822, for which F_rel was set to 0.88 (Retout 2026 Methods, PopPK analysis; 'Roche data on file'). Set to 0 for the currently marketed formulations.",
      source_name        = "10 mg tablet formulation (T0822)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1795,
    n_studies      = 6,
    n_observations = 6399,
    age_range      = "0.12-85 years",
    age_median     = "37.0 years",
    weight_range   = "4-217 kg",
    weight_median  = "65.7 kg",
    sex_female_pct = 50.2,
    race_ethnicity = c(Asian = 54.9, `Non-Asian` = 45.1),
    disease_state  = "Uncomplicated influenza; 63.0% otherwise healthy and 37.0% at high risk of influenza-related complications",
    dose_range     = "Single oral dose of baloxavir marboxil: 1 mg/kg or 2 mg/kg bodyweight-based dosing, or fixed doses of 5, 10, 20, 40 or 80 mg depending on study and bodyweight",
    regions        = "Asia 53.2%, North America / Europe 45.0%, other 1.8%",
    notes          = "Pooled from six treatment studies: T0821 (JapicCTI-153090), CAPSTONE-1 (NCT02954354), CAPSTONE-2 (NCT02949011), miniSTONE-2 (NCT03629184), T0822 (JapicCTI-163417) and T0833 (JapicCTI-173811). BLOCKSTONE (JapicCTI-184180, post-exposure prophylaxis) and T0835 (JapicCTI-194577) were not used for estimation; Bayesian post hoc estimates for those two studies were derived from this model. Baseline demographics: Retout 2026 Tables S1-S3."
  )

  ini({
    # Structural parameters - typical values for a 70 kg, 37-year-old, non-Asian male
    lcl   <- log(11.02); label("Apparent clearance (CL/F, L/h)")                                   # Retout 2026 Table 1: CL/F = 11.02 L/h (RSE 1.87%)
    lvc   <- log(735);   label("Apparent central volume of distribution (Vc/F, L)")                # Retout 2026 Table 1: Vc/F = 735 L (RSE 2.12%)
    lq    <- log(2.12);  label("Apparent intercompartmental clearance (Q/F, L/h)")                 # Retout 2026 Table 1: Q/F = 2.12 L/h (RSE 5.2%)
    lvp   <- log(260);   label("Apparent peripheral volume of distribution (Vp/F, L)")             # Retout 2026 Table 1: Vp/F = 260 L (RSE 6.68%)
    lka   <- log(1.39);  label("First-order absorption rate constant (ka, 1/h)")                   # Retout 2026 Table 1: ka = 1.39 1/h (RSE 14.6%)
    ltlag <- log(0.223); label("Absorption lag time (Tlag, h)")                                    # Retout 2026 Table 1: Tlag = 0.223 h (RSE 39.6%)

    # Relative bioavailability was set, not estimated
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the reference formulations (F_rel, unitless)")  # Retout 2026 Methods: F_rel set at 1.00 for all formulations except the T0822 10 mg tablet
    e_form_bxm_tab10_fdepot <- fixed(0.88); label("Relative bioavailability of the 10 mg baloxavir marboxil tablet used in T0822 (unitless)")  # Retout 2026 Methods: F_rel 0.88 (Roche data on file)

    # Covariate effects
    e_wt_cl_q       <- 0.467; label("Power (allometric) exponent on (WT/70) shared by CL/F and Q/F (unitless)")   # Retout 2026 Table 1: effect of BW on CL/F, Q/F = 0.467 (RSE 5.61%)
    e_wt_vc_vp      <- 0.887; label("Power (allometric) exponent on (WT/70) shared by Vc/F and Vp/F (unitless)")  # Retout 2026 Table 1: effect of BW on Vc/F, Vp/F = 0.887 (RSE 4.09%)
    e_race_asian_cl <- 0.504; label("Fractional reduction of CL/F in Asian vs non-Asian patients (unitless)")     # Retout 2026 Table 1: effect of race (Asian) on CL/F = 0.504 (RSE 2.4%)
    e_race_asian_vc <- 0.335; label("Fractional reduction of Vc/F in Asian vs non-Asian patients (unitless)")     # Retout 2026 Table 1: effect of race (Asian) on V/F = 0.335 (RSE 5.28%)
    e_race_asian_q  <- 0.391; label("Fractional reduction of Q/F in Asian vs non-Asian patients (unitless)")      # Retout 2026 Table 1: effect of race (Asian) on Q/F = 0.391 (RSE 7.21%)
    e_sexf_ka       <- 0.205; label("Fractional reduction of ka in female vs male patients (unitless)")           # Retout 2026 Table 1: effect of sex (female) on ka = 0.205 (RSE 36.4%)
    e_age_ka        <- 0.242; label("Power exponent on (AGE/37) for ka (unitless)")                               # Retout 2026 Table 1: effect of age on ka = 0.242 (RSE 21.4%)

    # Inter-individual variability, omega^2 = log(CV^2 + 1) from the CV% in Table 1.
    # CL/F 45.6%, Vc/F 45.7%, Q/F 48.6%, ka 113%, Tlag 62.6%. Off-diagonals are
    # corr * sqrt(omega2_i * omega2_j) using the correlations in Table 1:
    # CL-V 0.86, CL-Q 0.700, CL-ka 0.218, CL-Tlag -0.032, V-Q 0.899, V-ka 0.224,
    # V-Tlag -0.085, Q-ka 0.369, Q-Tlag -0.121, ka-Tlag -0.733.
    etalcl + etalvc + etalq + etalka + etaltlag ~
      c(0.188913,
        0.162790, 0.189669,
        0.140100, 0.180287, 0.212039,
        0.085949, 0.088490, 0.154129, 0.822815,
        -0.007998, -0.021286, -0.032039, -0.382332, 0.330653)  # Retout 2026 Table 1, random effects between-patient variability
    etalvp ~ fixed(0.022251)  # Retout 2026 Table 1: Vp/F between-patient CV 15%; omega^2 = log(1 + 0.15^2)

    # Residual error
    addSd  <- 0.257; label("Additive residual error (ng/mL)")         # Retout 2026 Table 1: sigma1 (additive) = 0.257 ng/mL (RSE 10.5%)
    propSd <- 0.149; label("Proportional residual error (fraction)")  # Retout 2026 Table 1: sigma2 (proportional) = 14.9% (RSE 4.09%)
  })
  model({
    # Race multipliers; non-Asian is the reference (RACE_ASIAN = 0).
    # Retout 2026 Table 1 footnote d:
    #   CL/F = 11.02 * (BW/70)^0.467 * (1 - 0.504 * Asian)
    #   Vc/F = 735   * (BW/70)^0.887 * (1 - 0.335 * Asian)
    #   Q/F  = 2.12  * (BW/70)^0.467 * (1 - 0.391 * Asian)
    #   Vp/F = 260   * (BW/70)^0.887                          (no race effect)
    race_cl <- 1 - e_race_asian_cl * RACE_ASIAN
    race_vc <- 1 - e_race_asian_vc * RACE_ASIAN
    race_q  <- 1 - e_race_asian_q  * RACE_ASIAN

    # ka = 1.39 * (1 - 0.205 * Sex) * (age/37)^0.242 (Table 1 footnote d)
    ka   <- exp(lka + etalka) * (1 - e_sexf_ka * SEXF) * (AGE / 37)^e_age_ka
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * race_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * race_vc
    q    <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q  * race_q
    vp   <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # F_rel is 1 for every formulation in the analysis except the 10 mg tablet
    # used in T0822, for which it was set to 0.88.
    f(depot)    <- exp(lfdepot) * (1 + (e_form_bxm_tab10_fdepot - 1) * FORM_BXM_TAB10)
    alag(depot) <- tlag

    # Amounts are in mg and volumes in L, so central/vc is mg/L; x1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
