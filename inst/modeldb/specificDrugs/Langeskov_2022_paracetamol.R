Langeskov_2022_paracetamol <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for a single 1500 mg oral dose of paracetamol (solubilised in yoghurt) in healthy obese adults with and without steady-state once-weekly subcutaneous semaglutide 1.0 mg (Langeskov 2022). Semaglutide co-administration reduces the absorption rate constant by a fraction of 0.525 and body weight enters the apparent central volume linearly; proportional residual error."
  reference <- paste(
    "Langeskov EK, Kristensen K. Population pharmacokinetic of paracetamol",
    "and atorvastatin with co-administration of semaglutide.",
    "Pharmacol Res Perspect. 2022;10(4):e00962. doi:10.1002/prp2.962"
  )
  vignette <- "Langeskov_2022_semaglutide_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "paracetamol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the apparent central volume linearly, V1/F = tvV1/F * (1 + 0.0312 * (WT - MedianBW)) (Langeskov 2022 Section 3.2 PML equation and Table 2). The paper does not print MedianBW; the cohort mean of 102 kg (Table 1; 102.3 kg in Section 2) is used as the centring value. The linear form reaches zero at WT = 102 - 1/0.0312 = 70 kg, so the model is only meaningful inside the studied 81.5-121 kg range.",
      source_name = "BW"
    ),
    CONMED_SEMAGLUTIDE = list(
      description = "Semaglutide co-administration indicator: 1 = paracetamol dosed at semaglutide steady state (once-weekly SC semaglutide escalated 0.25 -> 0.5 -> 1.0 mg over 12 weeks), 0 = paracetamol dosed in the placebo period.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo period, no semaglutide)",
      notes = "Per-occasion indicator in a two-period crossover. The source data-set column was named `placebo` but was coded 1 for semaglutide co-administration and 0 without (Section 3.2), so CONMED_SEMAGLUTIDE equals the source column without transformation. Enters ka as ka = tvka * (1 - 0.525 * CONMED_SEMAGLUTIDE).",
      source_name = "placebo"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    n_observations = 503L,
    age_range = "21-65 years",
    age_mean = "42 years",
    weight_range = "81.5-121 kg",
    weight_mean = "102 kg",
    bmi_range = "30.5-42.8 kg/m^2 (mean 33.2)",
    sex_female_pct = 31.0,
    race_ethnicity = "Not reported in source paper.",
    disease_state = "Healthy obese adults (BMI 30-45 kg/m^2) without diabetes.",
    dose_range = "Single oral 1500 mg paracetamol solubilised in yoghurt, given at the end of each 12-week period of once-weekly SC semaglutide (escalated 0.25 mg x 4 weeks, 0.5 mg x 4 weeks, 1.0 mg x 4 weeks) or placebo.",
    regions = "Single centre (trial NCT02079870, Hjerpsted 2018).",
    notes = "Randomised, double-blind, placebo-controlled two-period crossover; 29 subjects received paracetamol with placebo and 27 with semaglutide (two withdrew after period 1). Samples at 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4 and 5 h post-dose; all pre-dose samples were below LLOQ and excluded. Demographics from Langeskov 2022 Table 1 (20 male / 9 female). Fitted in Phoenix NLME 8.1 with FOCE-ELS."
  )

  ini({
    # Final-model estimates, Langeskov 2022 Table 2 (paracetamol).
    lka <- log(9.4); label("Absorption rate constant without semaglutide (1/h)") # Table 2 'ka (h-1) (Placebo)' = 9.4
    e_conmed_semaglutide_ka <- 0.525; label("Fractional reduction of ka with semaglutide co-administration (unitless)") # Table 2 'Theta kacovariate' = 0.525; Section 3.2 dKadplacebo1
    lvc <- log(48.5); label("Apparent central volume V1/F at the centring weight (L)") # Table 2 'V1/F (L)' = 48.5
    e_wt_vc <- 0.0312; label("Linear body-weight coefficient on V1/F (1/kg)") # Table 2 'Theta V1/Fcovariate' = 0.0312; Section 3.2 dVBW
    lcl <- log(25.9); label("Apparent oral clearance CL/F (L/h)") # Table 2 'Cl/F (L/h)' = 25.9
    lvp <- log(55.4); label("Apparent peripheral volume V2/F (L)") # Table 2 'V2/F (L)' = 55.4
    lq <- log(199); label("Apparent intercompartmental clearance Cl2/F (L/h)") # Table 2 'Cl2/F (L/h)' = 199
    ltlag <- log(0.16); label("Absorption lag time (h)") # Table 2 'Tlag (h)' = 0.16

    # Between-subject variability: Table 2 reports omega^2 (variances of
    # exponential etas, Section 2.1). No BSV on Cl2/F (Section 3.1).
    etalka ~ 0.514 # Table 2 'omega2 Ka' = 0.514
    etalvc ~ 0.270 # Table 2 'omega2 V1/F' = 0.270
    etalcl ~ 0.0606 # Table 2 'omega2 Cl/f' = 0.0606
    etalvp ~ 0.0703 # Table 2 'omega2 V2/F' = 0.0703
    etaltlag ~ 0.0586 # Table 2 'omega2 Tlag' = 0.0586

    # Proportional residual error. Table 2 prints 'Ceps' = 0.094; Phoenix
    # NLME reports CEps as the standard deviation of epsilon, and the tight
    # observed-vs-IPRED scatter of Figure 3A (about +/-10%) matches an SD of
    # 9.4% rather than the 31% that a variance reading (sqrt(0.094)) implies.
    propSd <- 0.094; label("Proportional residual error (fraction)") # Table 2 'Residual unexplained variability (Ceps)' = 0.094
  })

  model({
    # Individual parameters
    ka <- exp(lka + etalka) * (1 - e_conmed_semaglutide_ka * CONMED_SEMAGLUTIDE)
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 102))
    cl <- exp(lcl + etalcl)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L. The paper plots concentrations in
    # umol/L (paracetamol molecular weight 151.16 g/mol).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
