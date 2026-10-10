Liu_2022_dabigatran_healthy <- function() {
  description <- paste(
    "Population PK/PD model for a single oral dose of dabigatran etexilate in",
    "healthy Chinese adults (Liu 2022, 'Healthy volunteers model' columns of",
    "Table 2). PK: two-compartment model with first-order absorption and an",
    "absorption lag for total dabigatran, with age and HDL cholesterol power",
    "effects on apparent clearance. PD: direct linear effects of the",
    "dabigatran plasma concentration on activated partial thromboplastin time",
    "(aPTT), prothrombin time (PT) and anti-factor IIa activity, without",
    "covariates. The sibling Liu_2022_dabigatran is the paper's final",
    "integrated model refitted to healthy volunteers and patients with",
    "non-valvular atrial fibrillation."
  )
  reference <- paste(
    "Liu Y-O, Xie Q-F, Liu Z-Y, Wang Z, Mu G-Y, Zhang Y-T, Zhao Z-N, Yuan D-D,",
    "Guo L-P, Wang N, Xiang J, Song H-T, Jiang J, Xiang Q, Cui Y-M. Population",
    "pharmacokinetic analysis for dabigatran etexilate in Chinese patients with",
    "non-valvular atrial fibrillation. Front Cardiovasc Med. 2022;9:998751.",
    "doi:10.3389/fcvm.2022.998751"
  )
  vignette <- "Liu_2022_dabigatran"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on apparent clearance, (22 / AGE)^e_age_cl, with reference 22 years,",
        "the healthy-volunteer median (Liu 2022 Table 1). The printed",
        "Equation 3 form theta_cov^(COV / COVpop) is not used; see the",
        "vignette 'Assumptions and deviations' for the reasoning."
      ),
      source_name = "AGE"
    ),
    HDLC = list(
      description = "Serum high-density lipoprotein cholesterol",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on apparent clearance, (1.42 / HDLC)^e_hdlc_cl, with",
        "reference 1.42 mmol/L, the healthy-volunteer median (Liu 2022",
        "Table 1). Supply in mmol/L (divide mg/dL by 38.67)."
      ),
      source_name = "HDL-C"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dabigatran etexilate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "dabigatran (total, after alkaline hydrolysis of the acyl glucuronide)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "dabigatran (total)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 118L,
    n_studies = 1L,
    age_range = "18-41 years",
    age_median = "22 years",
    weight_range = "50.3-73.8 kg",
    weight_median = "60.3 kg",
    sex_female_pct = 27.97,
    race_ethnicity = "Chinese (all subjects)",
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Single oral 150 mg dabigatran etexilate (Pradaxa reference period of a",
      "four-period bioequivalence crossover), fasted (n = 59) or within 30 min",
      "after a standard high-fat meal (n = 59) (Liu 2022 Table 1)."
    ),
    regions = "China (Peking University First Hospital)",
    renal_function = "CrCL median 120.18 (86.96-171.02) mL/min (Liu 2022 Table 1)",
    notes = paste(
      "Baseline characteristics from Liu 2022 Table 1. 1,849 dabigatran",
      "concentrations (17-18 samples per subject) and 472 aPTT, 468 PT and",
      "447 anti-factor IIa observations (0, 2, 8 and 12 h post-dose). The PD",
      "model was fitted sequentially on the individual empirical-Bayes PK",
      "predictions. Study registered as NCT03161496."
    )
  )

  ini({
    # ---- PK (Liu 2022 Table 2, 'Healthy volunteers model' column) --------
    lka <- log(0.261)
    label("Absorption rate constant (1/h)") # Table 2 Ka = 0.261 (RSE 3.3%)
    ltlag <- fixed(log(0.634))
    label("Absorption lag time (h)") # Table 2 ALAG1 = 0.634, FIX (taken from the literature, Results 'PopPK model')
    lfdepot <- fixed(log(1))
    label("Relative bioavailability (unitless)") # Table 2 F1 = 1, FIX
    lcl <- log(166)
    label("Apparent clearance CL/F at AGE 22 y and HDLC 1.42 mmol/L (L/h)") # Table 2 CL/F = 166 (RSE 3.7%)
    lvc <- log(235)
    label("Apparent central volume V2/F (L)") # Table 2 V2/F = 235 (RSE 11.4%)
    lq <- fixed(log(35.5))
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 Q/F = 35.5, FIX
    lvp <- fixed(log(345))
    label("Apparent peripheral volume V3/F (L)") # Table 2 V3/F = 345, FIX

    e_age_cl <- 0.506
    label("Power exponent of (22 / AGE) on CL/F (unitless)") # Table 2 'Impact of AGE on CL/F' = 0.506 (RSE 37.4%)
    e_hdlc_cl <- 0.694
    label("Power exponent of (1.42 / HDLC) on CL/F (unitless)") # Table 2 'Impact of HDL on CL/F' = 0.694 (RSE 33.0%)

    # ---- aPTT (s) ---------------------------------------------------------
    lrbase_aPTT <- log(30.4)
    label("Baseline aPTT E0 (s)") # Table 2 APTT E0 = 30.4 (RSE 1.3%)
    lslope_aPTT <- log(0.0911)
    label("Linear aPTT slope K1 (s per ng/mL)") # Table 2 K1 = 0.0911 (RSE 7.4%)

    # ---- PT (s) -----------------------------------------------------------
    lrbase_PT <- log(11.7)
    label("Baseline PT E1 (s)") # Table 2 PT E1 = 11.7 (RSE 0.5%)
    lslope_PT <- log(0.00836)
    label("Linear PT slope K2 (s per ng/mL)") # Table 2 K2 = 0.00836 (RSE 9.2%)

    # ---- Anti-factor IIa activity (ng/mL) -----------------------------------
    lrbase_antiFIIa <- log(14.9)
    label("Baseline anti-factor IIa activity E2 (ng/mL)") # Table 2 Anti-IIa E2 = 14.9 (RSE 10.7%)
    lslope_antiFIIa <- log(0.628)
    label("Linear anti-factor IIa slope K3 (ng/mL per ng/mL)") # Table 2 K3 = 0.628 (RSE 8.5%)

    # ---- IIV (Table 2 'IIV' columns are NONMEM omega^2 variances) ----------
    etalcl ~ 0.142 # Table 2 healthy CL/F IIV 0.142 (RSE 15.4%)
    etalvc ~ 1.14 # Table 2 healthy V2/F IIV 1.14 (RSE 17.3%)
    etalrbase_aPTT ~ 0.0144 # Table 2 healthy APTT E0 IIV 0.0144 (RSE 22.5%)
    etalslope_PT ~ 0.0215 # Table 2 healthy K2 IIV 0.0215 (RSE 34.7%)
    etalslope_antiFIIa ~ 0.257 # Table 2 healthy K3 IIV 0.257 (RSE 24.9%)

    # ---- Residual error (Table 2 values are NONMEM sigma variances; SD = sqrt) ----
    addSd <- sqrt(221)
    label("Additive residual SD for dabigatran (ng/mL)") # Table 2 additive residual error 221 (variance, RSE 11.4%) -> SD 14.87
    addSd_aPTT <- sqrt(30.1)
    label("Additive residual SD for aPTT (s)") # Table 2 aPTT residual 30.1 (variance, RSE 9.7%) -> SD 5.49
    addSd_PT <- sqrt(0.592)
    label("Additive residual SD for PT (s)") # Table 2 PT residual 0.592 (variance, RSE 7.8%) -> SD 0.769
    addSd_antiFIIa <- sqrt(643)
    label("Additive residual SD for anti-factor IIa activity (ng/mL)") # Table 2 anti-IIa residual 643 (variance, RSE 11.5%) -> SD 25.4
  })
  model({
    # PK: Liu 2022 Table 2. Covariates enter as (reference / covariate)^exponent
    # with the healthy-volunteer medians as references (see vignette).
    ka <- exp(lka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl) * (22 / AGE)^e_age_cl * (1.42 / HDLC)^e_hdlc_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc

    # PD: Liu 2022 Equation 5, PD = E0 + K * X1, with X1 the dabigatran
    # plasma concentration (ng/mL).
    rbase_aPTT <- exp(lrbase_aPTT + etalrbase_aPTT)
    slope_aPTT <- exp(lslope_aPTT)
    rbase_PT <- exp(lrbase_PT)
    slope_PT <- exp(lslope_PT + etalslope_PT)
    rbase_antiFIIa <- exp(lrbase_antiFIIa)
    slope_antiFIIa <- exp(lslope_antiFIIa + etalslope_antiFIIa)

    aPTT <- rbase_aPTT + slope_aPTT * Cc
    PT <- rbase_PT + slope_PT * Cc
    antiFIIa <- rbase_antiFIIa + slope_antiFIIa * Cc

    Cc ~ add(addSd)
    aPTT ~ add(addSd_aPTT)
    PT ~ add(addSd_PT)
    antiFIIa ~ add(addSd_antiFIIa)
  })
}
