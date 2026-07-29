Luo_2016_efavirenz <- function() {
  description <- "Two-compartment population PK model with first-order absorption and elimination for oral efavirenz in pediatric HIV-1-infected patients (Luo 2016). Capsule / capsule-sprinkle formulation; body weight is a power covariate on CL/F, Vc/F, and Ka with reference 20 kg. The adult cohort (n = 24 healthy adults) and oral-solution formulation (study-specific Frel) reported in the same paper are documented in the validation vignette but not encoded as separate sub-models."
  reference <- paste(
    "Luo M, Chapel S, Sevinsky H, Savant I, Cirincione B, Bertz R, Roy A.",
    "Population pharmacokinetics analysis to inform efavirenz dosing",
    "recommendations in pediatric HIV patients aged 3 months to 3 years.",
    "Antimicrob Agents Chemother. 2016;60(6):3676-3686.",
    "doi:10.1128/AAC.02678-15.",
    sep = " "
  )
  vignette <- "Luo_2016_efavirenz"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Power covariate on CL/F (exponent 0.57), Vc/F (exponent 1.35), and Ka (exponent 0.768) with reference weight 20 kg (the approximate mean pediatric weight in the source studies). Luo 2016 Table 3 final model (updated data set).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    PART = list(
      description = "Prior antiretroviral therapy indicator (1 = yes; PACTG1021 study). Retained in Luo 2016 Table 3 final model on pediatric CL with point estimate 0.381 +/- 0.401 (SE > estimate). The paper Discussion explicitly cautions that this effect is likely confounded with the PACTG1021 study identity rather than a true PART effect.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not encoded in this packaged model. The pediatric dosing-recommendation simulations in Luo 2016 Figure 3 use the typical-value pediatric model (PART = 0); the user can re-introduce PART externally if needed. See vignette Assumptions and deviations."
    ),
    FORM_SOLUTION = list(
      description = "Oral-solution formulation indicator (1 = solution, 0 = capsule or capsule sprinkle). Luo 2016 Table 3 reports study-specific Frel for solution: PACTG382 -0.346 (Frel 0.654), PACTG1021 -0.0509 (Frel 0.949), AI266922 -0.754 (Frel 0.246).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not encoded. Capsule sprinkles were shown bioequivalent to capsules in adults (Luo 2016 Methods, study AI266059); the packaged model assumes capsule / capsule-sprinkle formulation with F = 1. The solution Frel adjustments were estimated to fit historical pediatric solution data and are study-specific."
    ),
    ADULT = list(
      description = "Adult-vs-pediatric cohort indicator. Luo 2016 estimated separate reference values for adult subjects (n = 24 healthy adults, study AI266059): CL_ref,adult = 3.66 L/h, Vc_ref,adult = 188 L, Ka same as pediatric reference 0.414 1/h, T_lag,adult = 0.633 h; Q (6.01 L/h) and Vp (287 L) were shared across cohorts.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not encoded. This file packages the pediatric structural model that drives the Luo 2016 pediatric dosing recommendations (the focus of the paper). Adult typical values and adult IIV (CL var 0.158, Vc var 0.132) are documented in the vignette."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 168L,
    n_studies      = 3L,
    age_range      = "0.2-24.7 years (pediatric studies; updated data set)",
    age_median     = "6.73 years (pediatric)",
    weight_range   = "3.3-117 kg",
    weight_median  = "25.3 kg (pediatric)",
    sex_female_pct = 52.4,
    race_ethnicity = c(White = 33.3, Black_AA = 52.4, Other = 14.3, Asian = 0),
    disease_state  = "HIV-1 infection on combination antiretroviral therapy. Pediatric patients 3 months to 21 years of age. Patients in PACTG382, PACTG1021, and AI266922 contributed PK data after at least 2 weeks of daily dosing (efavirenz autoinduction assumed to be at steady state).",
    dose_range     = "Pediatric efavirenz doses ranging from 100 mg QD (children >= 3.5 to < 5 kg) up to a maximum of 1000 mg QD per dosing algorithm. Formulations: oral solution, capsule, or capsule sprinkle.",
    regions        = "Multi-center pediatric AIDS clinical trials (PACTG382 and PACTG1021 Phase I/II open-label; AI266922 Phase II open-label).",
    notes          = "Updated data set (Luo 2016 Methods 'Study data' paragraph 2): 4521 plasma efavirenz concentrations from 192 subjects = 3289 concentrations from 168 pediatric patients (sparse) + 1232 concentrations from 24 healthy adults (intensive, study AI266059). The packaged model focuses on the pediatric sub-population (n = 168). Adult cohort demographics: mean age 32.8 y (20-45), mean weight 80.6 kg (59.6-98.1), 95.8% male, 50% White / 50% Black-AA. The model was fitted with NONMEM VI FOCE-I. Pharmacogenomic data (CYP2B6 15631GT and 21563CT) were available for 102 pediatric subjects and analyzed in two ad hoc analyses (covariate-on-CL and mixture model); the paper concluded CYP2B6 genotype information is not informative for guiding pediatric dosing and is not encoded in this packaged model."
  )

  ini({
    # =========================================================================
    # Structural disposition parameters -- pediatric reference values from
    # Luo 2016 Table 3 (final model, updated data set, right column). The
    # reference pediatric subject in Table 3 is: male, body weight 20 kg,
    # age 6 years, race non-African-American and not Other, no prior
    # antiretroviral therapy (PART = 0), no concomitant protease inhibitor
    # (PINT = 0). Capsule / capsule-sprinkle formulation (F = 1).
    # =========================================================================
    lcl <- log(4.8);   label("Apparent oral clearance CL/F (L/h) for a 20 kg reference pediatric subject")              # Luo 2016 Table 3 final 'CL_ref,ped = 4.8 +/- 0.33 L/h'
    lvc <- log(84.9);  label("Apparent oral central volume Vc/F (L) for a 20 kg reference pediatric subject")           # Luo 2016 Table 3 final 'V_C ref,ped = 84.9 +/- 8.13 L'
    lvp <- log(287);   label("Apparent oral peripheral volume Vp/F (L)")                                                # Luo 2016 Table 3 final 'V_P ref = 287 +/- 34.4 L' (shared between adult and pediatric per Luo 2016 Results 'Full-model development' final paragraph)
    lq  <- log(6.01);  label("Apparent oral intercompartmental clearance Q/F (L/h)")                                    # Luo 2016 Table 3 final 'Q_ref = 6.01 +/- 0.839 L/h' (shared between adult and pediatric per Luo 2016 Results 'Full-model development' final paragraph)
    lka <- log(0.414); label("Absorption rate constant Ka (1/h) for a 20 kg reference pediatric subject")               # Luo 2016 Table 3 final 'K_a ref = 0.414 +/- 0.0387 1/h'

    # =========================================================================
    # Body-weight power exponents on pediatric typical values (final model
    # covariate set: WT on CL, WT on Vc, WT on Ka). The covariate-parameter
    # functional form per Luo 2016 Methods 'Full-model development':
    #   P_TV = P1 * (R / R_ref)^P2
    # with R = WT and R_ref = 20 kg.
    # =========================================================================
    e_wt_cl <- 0.57;  label("Power exponent of (WT / 20) on pediatric CL/F (unitless)")  # Luo 2016 Table 3 final 'CL_ped,WT = 0.57 +/- 0.107'
    e_wt_vc <- 1.35;  label("Power exponent of (WT / 20) on pediatric Vc/F (unitless)")  # Luo 2016 Table 3 final 'V_C ped,WT = 1.35 +/- 0.152'
    e_wt_ka <- 0.768; label("Power exponent of (WT / 20) on pediatric Ka (unitless)")    # Luo 2016 Table 3 final 'K_a ped,WT = 0.768 +/- 0.0844'

    # =========================================================================
    # Inter-individual variability (variances, pediatric values from Luo 2016
    # Table 3 final model, updated data set). The paper specifies the IIV
    # model as log-normal: P_i = P_TV * exp(eta_i), eta_i ~ N(0, omega^2),
    # so omega^2 is directly the variance on the log scale (no CV-to-variance
    # conversion required).
    # =========================================================================
    etalcl ~ 0.602  # Luo 2016 Table 3 final 'IIV_CL ped = 0.602 +/- 0.231'
    etalvc ~ 0.234  # Luo 2016 Table 3 final 'IIV_V_C, ped = 0.234 +/- 0.0651'
    etalq  ~ 0.695  # Luo 2016 Table 3 final 'IIV_Q = 0.695 +/- 0.164' (shared between cohorts)
    etalvp ~ 0.296  # Luo 2016 Table 3 final 'IIV_V_P = 0.296 +/- 0.0877' (shared between cohorts)
    etalka ~ 0.202  # Luo 2016 Table 3 final 'IIV_K_a = 0.202 +/- 0.0570' (shared between cohorts)

    # =========================================================================
    # Residual variability -- Luo 2016 Methods 'Development of the base model'
    # paragraph 2: "The residual variability was also assumed to follow a
    # log-normal distribution, characterized by a log-transformed normal
    # distribution with a zero mean and variance sigma^2." That is, the
    # paper explicitly used a log-normal (multiplicative-exp) error
    # structure: Y_obs = Y_pred * exp(eps), eps ~ N(0, sigma^2). In nlmixr2
    # this is the lnorm error model with expSd = sqrt(sigma^2). The packaged
    # model uses the pediatric capsule / capsule-sprinkle residual (variance
    # 0.461 -> expSd = sqrt(0.461) = 0.679). The pediatric oral-solution
    # residual (variance 0.784) and the adult residual (variance 0.212) are
    # documented in the vignette but not encoded here (capsule /
    # capsule-sprinkle is the formulation used in the Luo 2016 pediatric
    # dosing-recommendation simulations).
    # =========================================================================
    expSd <- sqrt(0.461); label("Log-normal residual error SD (unitless; capsule / capsule-sprinkle pediatric)")  # Luo 2016 Table 3 final 'Capsule, pediatric studies = 0.461 +/- 0.0286' (log-scale variance; expSd = sqrt(0.461))
  })

  model({
    # 1. Individual PK parameters with body-weight power covariates
    #    (reference WT = 20 kg per the Luo 2016 pediatric reference subject)
    cl <- exp(lcl + etalcl) * (WT / 20)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 20)^e_wt_vc
    ka <- exp(lka + etalka) * (WT / 20)^e_wt_ka
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # 2. Micro-constants for the explicit 2-compartment ODE system
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. Two-compartment ODE system with first-order oral absorption
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    # 4. Concentration in mg/L (= ug/mL): dose in mg, volume in L
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
