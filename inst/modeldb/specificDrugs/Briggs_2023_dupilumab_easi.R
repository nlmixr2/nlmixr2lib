Briggs_2023_dupilumab_easi <- function() {
  description <- "Exposure-response model for the Eczema Area and Severity Index (EASI) under dupilumab in children (6-11 years), adolescents and adults with moderate-to-severe atopic dermatitis (Briggs 2023). PD-only: the time-varying serum functional dupilumab concentration is an input column (CP_DUPILUMAB_MGL, mg/L), as in the paper's sequential PK-then-PD fit to individual predictions from the age-group dupilumab popPK models (see modellib('Kovalenko_2020_dupilumab_covariate'), modellib('Kovalenko_2021_dupilumab_adolescent_covariate') and modellib('Kovalenko_2021_dupilumab_children_covariate')). The model is fitted on the Tukey-transformed scale (EASI + 1)^0.4 as the sum of a baseline, a hyperbolic (Pmax * t / (ET50 + t)) placebo time course and a type I indirect-response drug effect with Imax inhibition of the production rate. Covariates: topical corticosteroids, prior systemic immunotherapy, age, eosinophil count and TARC on the baseline; topical corticosteroids and prior systemic immunotherapy on Pmax; body weight, Asian race and eosinophil count on Imax. Companion IGA model from the same paper: modellib('Briggs_2023_dupilumab_iga')."

  reference <- paste(
    "Briggs E, Kamal MA, Kosloski MP, Linsmeier I, Jusko N, Dolphin N,",
    "Chittenden J, Simpson EL, Paller AS, Siegfried EC, Shumel B, Levit NA,",
    "Bansal A, Davis JD, Chapel S, Smith DE, Huniti N. Integrated",
    "Exposure-Response of Dupilumab in Children, Adolescents, and Adults With",
    "Atopic Dermatitis Using Categorical and Continuous Efficacy Assessments:",
    "A Population Analysis. Pharm Res. 2023;40(11):2653-2666.",
    "doi:10.1007/s11095-023-03616-8.",
    sep = " "
  )
  vignette <- "Briggs_2023_dupilumab"

  units <- list(
    time = "day",
    dosing = "(none; dupilumab exposure enters as the CP_DUPILUMAB_MGL concentration column, not as dosing events)",
    concentration = "(EASI score; the observation easi_tukey is (EASI + 1)^0.4, and easi is the back-transformed score on the 0-72 clinical scale; the drug input CP_DUPILUMAB_MGL is in mg/L)"
  )

  compartmentData <- list(
    latent_easi = list(
      analyte = "latent indirect-response variable driving the EASI drug effect",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_DUPILUMAB_MGL = list(
      description = "Serum functional dupilumab concentration",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input supplied on every record: 0 on placebo. In the source analysis it was the individual prediction of the age-group dupilumab popPK model (Briggs 2023 Online Resource 'Population Pharmacokinetic Model' and Table S1, which reprints Kovalenko 2020 / 2021), evaluated with each patient's actual dosing history. The paper's covariate and exposure-response simulations hold it constant (68 mg/L in Figure 3; a range of trough values in Figure 5).",
      source_name = "C (dupilumab concentration)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE / 30)^-0.0372 on the transformed-scale EASI baseline. 30 years is the Figure 3 reference patient; the paper does not print the functional form or the centring value, which were chosen because they reproduce the Figure 3 forest plot.",
      source_name = "Age"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (WT / 70)^-0.185 on Imax, centred on the 70 kg Figure 3 reference patient (form and centring not printed in the paper; chosen to reproduce Figure 3).",
      source_name = "Baseline body weight"
    ),
    EOS = list(
      description = "Baseline blood eosinophil count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper reports eosinophils in 10^9/L; multiply by 1000 to get cells/uL. Table II names the effects 'Log baseline eosinophil count', i.e. theta * log(EOS / ref), which is the power form (EOS / 481)^theta used here: exponent 0.0173 on the transformed-scale baseline and -0.0570 on Imax. 481 cells/uL (0.481 x 10^9/L) is the median observed value used for the Figure 3 reference patient.",
      source_name = "EOS"
    ),
    TARC = list(
      description = "Baseline serum thymus and activation-regulated chemokine (TARC / CCL17)",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table II 'Log baseline TARC on baseline' = 0.0508, encoded as (TARC / 2033.6)^0.0508 on the transformed-scale baseline. 2033.6 pg/mL is the median observed value used for the Figure 3 reference patient.",
      source_name = "TARC"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = any other race)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian; the Figure 3 reference patient is White)",
      notes = "Linear effect (1 - 0.374 * RACE_ASIAN) on Imax. Black and Other race were screened (Table S2) but not retained.",
      source_name = "Asian race"
    ),
    CONMED_STEROID_TOPICAL = list(
      description = "Concomitant topical corticosteroid (TCS) therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dupilumab or placebo without TCS)",
      notes = "Subject-level grain, set by the study design: all patients in the children study (NCT03345914) and the adult CHRONOS study (NCT02260986) received TCS, no patient in the other four studies did (Briggs 2023 Table I and Table S3). Linear effects (1 - 0.0218 * TCS) on the transformed-scale baseline and (1 + 0.232 * TCS) on Pmax.",
      source_name = "TCS co-administration"
    ),
    PRIOR_SYSTEMIC_IMMUNOTHERAPY = list(
      description = "Prior exposure to systemic immunotherapy for atopic dermatitis (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior systemic immunotherapy)",
      notes = "Linear effects (1 + 0.0233 * PRIOR) on the transformed-scale baseline and (1 - 0.166 * PRIOR) on Pmax. Briggs 2023 Table S3 tabulates the column as 'Prior exposure to systemic immunotherapies' (31.8% yes overall).",
      source_name = "Prior exposure to systemic immunotherapies"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened with Other race on the placebo and drug-effect parameters (Briggs 2023 Table S2) and removed by backward elimination (alpha = 0.001)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "The Methods list gender among the baseline demographics considered, but Table S2 does not list it among the covariates evaluated on EASI and it is not in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2968L,
    n_studies = 6L,
    age_range = "6-75 years (children 6-11 n = 359, adolescents 12-17 n = 243, adults 18-75 n = 2366); mean 32.2 (SD 16.4) years",
    weight_range = "mean 69.8 (SD 23.2) kg overall; adults 76.1, adolescents 65.0, children 31.5 kg",
    sex_female_pct = 41.9,
    race_ethnicity = "White 67.2%, Asian 21.1%, Black 8.1%, Other / not reported 3.7%; 9.0% Japanese (all adults)",
    disease_state = "Moderate-to-severe atopic dermatitis (baseline IGA 3 in 45.8%, IGA 4 in 54.2%); median baseline EASI 30.8 (range 10.7-72)",
    dose_range = "Dupilumab SC with a Day-1 loading dose of twice the maintenance dose: adults 100-300 mg q4w to 300 mg qw; adolescents 200 or 300 mg q2w or 300 mg q4w; children 100 or 200 mg q2w or 300 mg q4w with TCS; or placebo (with or without TCS)",
    regions = "Multinational Phase IIb / III programme (NCT01859988, NCT02277769, NCT02277743, NCT02260986, NCT03054428, NCT03345914)",
    co_medication = "Topical corticosteroids in 35.6% (all children and the adult CHRONOS study)",
    notes = "Demographics from Briggs 2023 Table S3; 29,413 EASI observations. 66.2% received dupilumab, 33.8% placebo."
  )

  ini({
    # Briggs 2023 Table II (final EASI model) and Online Resource equations:
    #   (EASI + 1)^0.4 = BASE + Pmax * t / (ET50 + t) + IDR
    #   dR/dt = kin * (1 - Imax * C / (IC50 + C)) - kout * R, kin = kout * BASE
    # with R(0) = BASE and IDR = R - BASE.
    lrbase <- log(30.4)
    label("Typical baseline EASI score, converted to the transformed scale as (BASE + 1)^0.4 (EASI units)") # Table II 'Baseline EASI score units' = 30.4 (RSE 0.4%)
    pmax <- -1.24
    label("Maximum placebo effect on the transformed (EASI + 1)^0.4 scale (unitless)") # Table II 'Pmax' = -1.24 (RSE 2.6%)
    let50 <- log(29.2)
    label("Time to half the maximum placebo effect (day)") # Table II 'ET50 (days)' = 29.2 (95% CI 27.5-31.0)
    limax <- log(0.266)
    label("Maximum fractional inhibition of the EASI production rate by dupilumab (fraction)") # Table II 'Drug effect (Imax)' = 0.266 (RSE 2.7%)
    lic50 <- log(20.3)
    label("Dupilumab concentration giving half the maximum inhibition (mg/L)") # Table II 'IC50 (mg/L)' = 20.3 (95% CI 16.1-25.5)
    lthalf_onset <- log(13.7)
    label("Half-life of drug effect onset; kout = log(2) / half-life (day)") # Table II 'Half-life for drug effect onset (day)' = 13.7 (95% CI 12.9-14.7)

    e_conmed_steroid_topical_rbase <- -0.0218
    label("Fractional change in the transformed-scale baseline with topical corticosteroids (unitless)") # Table II 'TCS on baseline' = -0.0218 (RSE 27.9%)
    e_prior_systemic_immunotherapy_rbase <- 0.0233
    label("Fractional change in the transformed-scale baseline with prior systemic immunotherapy (unitless)") # Table II 'Prior immunotherapy on baseline' = 0.0233 (RSE 27.5%)
    e_age_rbase <- -0.0372
    label("Power exponent of AGE/30 on the transformed-scale baseline (unitless)") # Table II 'Age on baseline' = -0.0372 (RSE 13.4%)
    e_eos_rbase <- 0.0173
    label("Power exponent of EOS/481 on the transformed-scale baseline (unitless)") # Table II 'Log baseline eosinophil count on baseline' = 0.0173 (RSE 15.6%)
    e_tarc_rbase <- 0.0508
    label("Power exponent of TARC/2033.6 on the transformed-scale baseline (unitless)") # Table II 'Log baseline TARC on baseline' = 0.0508 (RSE 4.4%)
    e_conmed_steroid_topical_pmax <- 0.232
    label("Fractional change in Pmax with topical corticosteroids (unitless)") # Table II 'TCS on Pmax' = 0.232 (RSE 16.9%)
    e_prior_systemic_immunotherapy_pmax <- -0.166
    label("Fractional change in Pmax with prior systemic immunotherapy (unitless)") # Table II 'Prior immunotherapy on Pmax' = -0.166 (RSE 18.2%)
    e_wt_imax <- -0.185
    label("Power exponent of WT/70 on Imax (unitless)") # Table II 'Baseline body weight on drug effect' = -0.185 (RSE 23.3%)
    e_race_asian_imax <- -0.374
    label("Fractional change in Imax for Asian race (unitless)") # Table II 'Asian race on drug effect' = -0.374 (RSE 8.3%)
    e_eos_imax <- -0.0570
    label("Power exponent of EOS/481 on Imax (unitless)") # Table II 'Log baseline EOS on drug effect' = -0.0570 (RSE 29.3%)

    # IIV: exponential on the transformed-scale baseline, additive on Pmax
    etalrbase ~ 0.0177 # Table II 'EASI baseline - exponential (omega2)' = 0.0177 (RSE 3.4%; shrinkage 11.8%)
    etapmax ~ 0.788 # Table II 'EASI Pmax - additive (omega2)' = 0.788 (RSE 3.7%; shrinkage 14.3%)

    # Residual error is additive on the transformed scale. 0.422 is read as a
    # standard deviation: it reproduces the simulated 5th-95th percentile
    # band at time 0 of the Figure 2 VPC, which sqrt(0.422) does not.
    addSd <- 0.422
    label("Additive residual error on (EASI + 1)^0.4 (unitless)") # Table II 'Residual variability - additive' = 0.422 (RSE 0.5%)
  })

  model({
    # Tukey transformation exponent (Online Resource: '(EASI + 1)^lambda
    # where lambda = 0.4').
    lambda <- 0.4

    # Baseline on the transformed scale with exponential IIV; covariate
    # functional forms and centring values reproduce Figure 3.
    rbase_t <- (exp(lrbase) + 1)^lambda *
      (1 + e_conmed_steroid_topical_rbase * CONMED_STEROID_TOPICAL) *
      (1 + e_prior_systemic_immunotherapy_rbase * PRIOR_SYSTEMIC_IMMUNOTHERAPY) *
      (AGE / 30)^e_age_rbase *
      (EOS / 481)^e_eos_rbase *
      (TARC / 2033.6)^e_tarc_rbase *
      exp(etalrbase)

    # Placebo (non-drug) time course with additive IIV on Pmax
    pmax_i <- pmax *
      (1 + e_conmed_steroid_topical_pmax * CONMED_STEROID_TOPICAL) *
      (1 + e_prior_systemic_immunotherapy_pmax * PRIOR_SYSTEMIC_IMMUNOTHERAPY) +
      etapmax
    et50 <- exp(let50)
    placebo <- pmax_i * time / (et50 + time)

    # Indirect response with Imax inhibition of kin. The state is the
    # paper's R divided by BASE, so it starts at 1 and kin = kout.
    imax <- exp(limax) *
      (WT / 70)^e_wt_imax *
      (1 + e_race_asian_imax * RACE_ASIAN) *
      (EOS / 481)^e_eos_imax
    ic50 <- exp(lic50)
    thalf_onset <- exp(lthalf_onset)
    kout <- log(2) / thalf_onset
    kin <- kout
    inh <- imax * CP_DUPILUMAB_MGL / (ic50 + CP_DUPILUMAB_MGL)
    latent_easi(0) <- 1
    d / dt(latent_easi) <- kin * (1 - inh) - kout * latent_easi
    idr <- rbase_t * (latent_easi - 1)

    # Observation on the transformed scale and back-transformed EASI score.
    # A transformed value below 1 would back-transform to a negative EASI,
    # so the score is floored at 0.
    easi_tukey <- rbase_t + placebo + idr
    easi <- max(easi_tukey, 1)^(1 / lambda) - 1
    easi_tukey ~ add(addSd)
  })
}
