Zhang_2026_tebipenem <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for tebipenem, the active moiety of the oral",
    "carbapenem pro-drug tebipenem pivoxil, in children. The structural",
    "model and every parameter estimate are carried unchanged from the",
    "Sato 2008 Japanese pediatric popPK analysis (217 patients aged 0.5-16",
    "years with otolaryngological infection or pneumonia). Zhang 2026",
    "re-expressed that model so it could be driven by tebipenem dose rather",
    "than tebipenem-pivoxil dose, applying two rescalings documented in its",
    "Supporting Information: both apparent clearance and apparent volume are",
    "multiplied by the tebipenem : tebipenem-pivoxil molecular-weight ratio",
    "383.5 / 497.63 = 0.7707, and the creatinine-clearance slope is",
    "unit-converted. Weight-normalized apparent clearance is a linear",
    "function of weight-normalized creatinine clearance; weight-normalized",
    "apparent volume is a power function of age. Zhang 2026 used the model to",
    "simulate tebipenem exposure in Bangladeshi children aged 24-59 months",
    "with shigellosis and to compute the fraction of the 72-hour treatment",
    "period during which the free plasma concentration exceeds the Shigella",
    "MIC (40% fT > MIC). The plasma unbound fraction 0.58 is carried in the",
    "model as fu so the free-concentration driver Ccfree is available",
    "directly. Distinct from Ganesan_2023_tebipenem.R, which is a",
    "two-compartment transit-absorption model in adults with complicated",
    "urinary tract infection."
  )
  reference <- paste(
    "Zhang CX, Nuzhat S, Islam MR, Bashar SJ, Das S, Amin R, Qadri F,",
    "Khanam F, Ahmed D, Pavlinac PB, Chisti MJ, Ahmed T, Arnold SLM (2026).",
    "Population pharmacokinetic and pharmacodynamic prediction for",
    "tebipenem pivoxil treatment of pediatric shigellosis.",
    "Clinical and Translational Science 19(1): e70453.",
    "doi:10.1111/cts.70453.",
    "Structural model and all parameter estimates from",
    "Sato N, Kijima K, Koresawa T, Kanazu T, Itoh Y, Ito Y, Yamaguchi Y,",
    "Sunakawa K, Totsuka K, Miyazaki S, Nakayama I (2008).",
    "Population pharmacokinetics of tebipenem pivoxil (ME1211), a novel oral",
    "carbapenem antibiotic, in pediatric patients with otolaryngological",
    "infection or pneumonia.",
    "Drug Metabolism and Pharmacokinetics 23(6): 434-446.",
    "doi:10.2133/dmpk.23.434.",
    sep = " "
  )
  vignette <- "Zhang_2026_tebipenem"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot   = list(analyte = "tebipenem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tebipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; converts the weight-normalized apparent clearance (L/h/kg) and weight-normalized apparent volume (L/kg) of the Sato 2008 parameterization into whole-body CL/F (L/h) and V/F (L)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters model() as a plain linear multiplier on both CL/F and V/F, not",
        "as an allometric power term: Sato 2008 fitted 'Basic Model 2', whose",
        "primary parameters are weight-normalized apparent clearance CL/F",
        "(L/h/kg) and weight-normalized apparent volume Vd/F (L/kg) with dosage",
        "expressed as amount per body weight (Sato 2008 Methods, 'Population",
        "pharmacokinetic analysis ...'). There is therefore no reference weight",
        "and no exponent to estimate.",
        "Sato 2008 estimation cohort median (range) weight 15.4 (7.06-49.5) kg",
        "(Sato 2008 Results, 'Assay of plasma tebipenem concentrations'",
        "paragraph, and Table 2).",
        "Zhang 2026 simulation target population weight mean (SD) 11.89 (2.64)",
        "kg in girls and 12.42 (2.58) kg in boys aged 24-59 months",
        "(Zhang 2026 Table 1). Weight also sets the mg/kg dose."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age in years; power covariate on weight-normalized apparent volume of distribution",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as the uncentred power term AGE^-0.132 (Sato 2008 Table 4:",
        "Vd/F (L/kg) = 1.18 x Age^-0.132), so the multiplier is 1 at age 1 year",
        "and the coefficient 1.18 L/kg is the 1-year-old value. No reference age",
        "is used and none is reported.",
        "Sato 2008 estimation cohort median (range) age 4.50 (0.67-15.23) years.",
        "Zhang 2026 simulation target population 24-59 months (2.0-4.9 years),",
        "mean 36.27 months in girls and 35.58 months in boys (Zhang 2026",
        "Table 1); the model is therefore extrapolated only mildly within the",
        "Sato 2008 age range."
      ),
      source_name        = "Age"
    ),
    CRCL = list(
      description        = "Weight-normalized creatinine clearance; linear covariate on weight-normalized apparent clearance",
      units              = "mL/min/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WEIGHT-normalized (mL/min/kg), NOT the BSA-normalized mL/min/1.73 m^2",
        "form that is the CRCL register default. Sato 2008 derived it from",
        "serum creatinine with the pediatric Schwartz-type predictive equations",
        "Ccr (mL/min/kg) = 0.00187 x HT^1.725 / (WT^0.575 x Scr) for infants",
        "under one year and Ccr (mL/min/kg) = 0.00199 x HT^1.725 / (WT^0.575 x",
        "Scr) for children aged 1-18 years (Sato 2008 Equations 5 and 6), with",
        "HT in cm, WT in kg and Scr in mg/dL, using the value on the first day",
        "of administration.",
        "Sato 2008 cohort mean (SD) 3.72 (0.97) mL/min/kg, median 3.73, range",
        "1.46-6.97 (Sato 2008 Table 2).",
        "Creatinine clearance was NOT collected in the Zhang 2026 screening",
        "data set; Zhang 2026 sampled it for every virtual participant from a",
        "normal distribution with the Sato 2008 mean and SD (Zhang 2026 Methods",
        "2.3, final paragraph, and the deposited virtual-population script in",
        "Appendix S1, 'Randomly generate Ccr from normal distribution').",
        "The effect on CL/F is additive linear with no centering:",
        "cl = mwRatio * (exp(lcl) + e_crcl_cl * CRCL) * WT.",
        "Zhang 2026 Equation (6) prints the equivalent form with CrCl supplied",
        "in L/h/kg and the slope multiplied by 1000/60; this file uses Sato's",
        "native mL/min/kg so the published slope 0.104 applies directly."
      ),
      source_name        = "Ccr"
    )
  )

  # Documented in the source but not used by model(): sex shapes the Zhang 2026
  # virtual population (it selects the sex-specific height-from-age and
  # weight-from-height regressions and the 60% male sampling probability) but
  # is not a covariate on any PK parameter in Sato 2008's final model.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Not a PK covariate. Used only to construct the Zhang 2026 virtual",
        "population: virtual participants were assigned male with probability",
        "0.6 (Zhang 2026 Methods 2.3), and the age-to-height and",
        "height-to-weight linear regressions in Zhang 2026 Table S1 are",
        "sex-specific. Sato 2008 did not test sex as a covariate; the only",
        "categorical covariate it screened was FAC1 (0 = otolaryngological",
        "infection, 1 = bacterial or mycoplasmal pneumonia) on CL/F and Vd/F",
        "(Sato 2008 Equations 8 and 9), and both FAC1 terms were dropped from",
        "the Full Model by the likelihood-ratio test (Sato 2008 Table 3 Final",
        "Model, which retains only theta2 + theta3 x Ccr on CL/F)."
      ),
      source_name        = "Sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 217L,
    n_studies      = 5L,
    age_range      = "0.67-15.23 years",
    age_median     = "4.50 years",
    weight_range   = "7.06-49.5 kg",
    weight_median  = "15.4 kg",
    sex_female_pct = 48.4,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "pediatric otolaryngological infection or bacterial / mycoplasmal pneumonia",
    dose_range     = "Oral tebipenem pivoxil 4 or 6 mg/kg twice daily",
    regions        = "Japan",
    renal_function = "Weight-normalized creatinine clearance mean (SD) 3.72 (0.97) mL/min/kg, median 3.73, range 1.46-6.97 (Sato 2008 Table 2)",
    notes          = paste(
      "The PARAMETER-ESTIMATION population is Sato 2008's: 217 Japanese",
      "children (112 boys, 105 girls) pooled from five clinical studies,",
      "median age 4.50 years, median weight 15.4 kg, median serum creatinine",
      "0.30 mg/dL (Sato 2008 Results and Table 2). NONMEM V, ADVAN2 / TRANS2,",
      "one-compartment model with first-order absorption, relative (proportional)",
      "inter-subject error model and additive residual error.",
      "The SIMULATION TARGET population of Zhang 2026 is different: 2249",
      "Bangladeshi children aged 24-59 months screened for the shigellosis",
      "trial NCT05121974 between May and September 2022 (59.5% male; girls",
      "mean age 36.27 months, height 91.68 cm, weight 11.89 kg; boys 35.58",
      "months, 92.66 cm, 12.42 kg; Zhang 2026 Table 1). Zhang 2026 built a",
      "33,000-participant virtual population (500 trials x 66 participants)",
      "from those demographics and simulated 4 mg/kg tebipenem pivoxil three",
      "times daily for 3 days and 3 mg/kg four times daily for 3 days.",
      "Zhang 2026 compared the Sato 2008 model against 15 Bangladeshi children",
      "from the trial's pilot study and reported a median prediction error of",
      "21.6% and a median absolute prediction error of 40.9% (Zhang 2026",
      "Methods 2.4 and Figure S1); the model was judged adequate and no",
      "re-estimation was performed, so every value in ini() is fixed."
    ),
    target_population = paste(
      "2249 Bangladeshi children aged 24-59 months with suspected Shigella",
      "infection, screened May-September 2022 for trial NCT05121974",
      "(Zhang 2026 Table 1)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # EVERY parameter below is FIXED. Zhang 2026 did not estimate anything:
    # it adopted Sato 2008's final estimates unchanged (Zhang 2026 Methods
    # 2.4: "Hence, there was insufficient data for building a new population
    # PK model. Therefore, Sato et al.'s model was deemed adequate and used
    # for subsequent simulations."). The two rescalings Zhang 2026 applied
    # (molecular-weight ratio, creatinine-clearance unit conversion) live in
    # model(), not here, so the values below are Sato 2008's as published.
    # ---------------------------------------------------------------------

    lka <- fixed(log(5.85))
    label("First-order absorption rate constant ka (1/h)")
    # Sato 2008 Table 4 population mean ka = 5.85 1/h; Table 5 final estimate
    # theta1 = 5.85 +/- 1.02 (S.E.). Deposited Zhang 2026 mrgsolve stream
    # (Appendix S1 file CTS-2025-0490-s02.docx): $PARAM TVKA = 5.85.

    lcl <- fixed(log(0.363))
    label("Creatinine-clearance-independent intercept of weight-normalized CL/F (L/h/kg)")
    # Sato 2008 Table 4: CL/F (L/h/kg) = 0.363 + 0.104 x Ccr; Table 5 final
    # estimate theta2 = 0.363 +/- 0.138. Deposited stream: $PARAM TVCL = 0.363.
    # Reproduced as the bracketed intercept of Zhang 2026 Equation (6).

    e_crcl_cl <- fixed(0.104)
    label("Slope of weight-normalized creatinine clearance on weight-normalized CL/F ((L/h/kg) per (mL/min/kg))")
    # Sato 2008 Table 4 (same equation as lcl); Table 5 final estimate
    # theta3 = 0.104 +/- 0.0354. Deposited stream: the literal 0.104 in
    # "double CL = 0.7707*(TVCL + 0.104*1000/60*CCR) * exp(ECL) * WT;", where
    # the 1000/60 converts the stream's L/h/kg CCR input back to the
    # mL/min/kg in which this slope is defined. This file supplies CRCL in
    # mL/min/kg directly, so no conversion factor appears in model().

    lvc <- fixed(log(1.18))
    label("Weight-normalized V/F at age 1 year (L/kg)")
    # Sato 2008 Table 4: Vd/F (L/kg) = 1.18 x Age^-0.132; Table 5 final
    # estimate theta5 = 1.18 +/- 0.141. Deposited stream: $PARAM TVV = 1.18.

    e_age_vc <- fixed(-0.132)
    label("Power exponent of age (years) on weight-normalized V/F (unitless)")
    # Sato 2008 Table 4 prints the exponent with the minus sign inside the
    # equation (Age^-0.132) and Table 5 reports the magnitude as
    # theta6 = 0.132 +/- 0.0885. Stored here with the sign carried by the
    # coefficient so model() can write AGE^e_age_vc. Deposited stream:
    # "-0.132 * log(AGE)" inside the exp() for V.

    ltlag <- fixed(log(0.239))
    label("Absorption lag time (h)")
    # Sato 2008 Table 4: Tlag (hr) = 0.239; Table 5 final estimate
    # theta8 = 0.239 +/- 0.0377. Deposited stream: $PARAM TVALAG1 = 0.239.

    fu <- fixed(0.58)
    label("Plasma unbound fraction of tebipenem")
    # Zhang 2026 Methods 2.5: "The free tebipenem concentrations were
    # calculated from the total concentrations using a mean unbound plasma
    # fraction of 0.58", citing Rodvold 2022 (doi:10.1128/aac.00590-22).
    # Deposited stream: "sims.total.df$freeTBPM <- sims.total.df$IPRED * 0.58".

    # ---------------------------------------------------------------------
    # Interindividual variability. Sato 2008 reports these both as variances
    # (Table 5) and as the corresponding CV% (Table 4); the two agree exactly
    # under Sato's RELATIVE (proportional) error model, CV = sqrt(omega^2):
    # sqrt(1.28) = 113%, sqrt(0.0407) = 20.2%, sqrt(0.558) = 74.7%,
    # sqrt(2.10) = 145%. The deposited Zhang 2026 stream applies these same
    # variances as LOG-NORMAL (exponential) IIV, which is what is encoded
    # here so that this file reproduces Zhang 2026's published results; see
    # the vignette "Assumptions and deviations" section.
    # ---------------------------------------------------------------------

    etalcl ~ fixed(0.0407)
    # Sato 2008 Table 5 omega^2(CL/F) = 0.0407 +/- 0.0183; Table 4 omega(CL/F)
    # = 20.2%. Deposited stream $OMEGA first element (label ECL) = 0.0407.

    etalvc ~ fixed(0.5558)
    # Deposited stream $OMEGA second element (label EV) = 0.5558. Sato 2008
    # Table 5 reports omega^2(Vd/F) = 0.558 +/- 0.0968 and Table 4 reports
    # omega(Vd/F) = 74.7% (sqrt(0.558) = 0.747). The deposited value carries
    # one extra digit; the 0.4% difference is immaterial and the deposited
    # value is used because it is what generated the published results.

    etalka ~ fixed(1.28)
    # Sato 2008 Table 5 omega^2(ka) = 1.28 +/- 0.981; Table 4 omega(ka) =
    # 113%. Deposited stream $OMEGA third element (label EKA) = 1.28.

    etaltlag ~ fixed(2.1)
    # Sato 2008 Table 5 omega^2(Tlag) = 2.10 +/- 1.28; Table 4 omega(Tlag) =
    # 145%. Deposited stream $OMEGA fourth element (label EALAG1) = 2.1.

    addSd <- fixed(0.453)
    label("Additive residual error (ug/mL)")
    # Sato 2008 Table 4 "Intra individual variance, sigma (ug/mL) 0.453" is
    # the residual standard deviation, corroborated by Table 5, which reports
    # the variance sigma^2 = 0.205 +/- 0.0584 (sqrt(0.205) = 0.4528). Sato
    # 2008 Methods states "an additive error model was used for residual
    # variation". The deposited Zhang 2026 stream writes "$SIGMA 0.453",
    # which mrgsolve interprets as a VARIANCE, i.e. an SD of 0.673; that slip
    # does not affect any published result because every Zhang 2026 analysis
    # is computed from IPRED, never from DV. nlmixr2's addSd is a standard
    # deviation, so the Sato 2008 published SD is used here.
  })

  model({
    # -------------------------------------------------------------------
    # Zhang 2026 correction (Supporting Information Appendix S1, file
    # CTS-2025-0490-s01.docx): Sato 2008 fitted CL/F and Vd/F against the
    # tebipenem-PIVOXIL (pro-drug) dose. To drive the model with a
    # tebipenem dose instead, both are multiplied by the molecular-weight
    # ratio tebipenem / tebipenem-pivoxil = 383.5 / 497.63 = 0.7707. This
    # reproduces the leading factor in Zhang 2026 Equations (6) and (7).
    # Because the dose supplied to the model must be rescaled by the same
    # factor, the rescaling is a re-parameterization: it leaves the
    # predicted concentration-time profile identical to Sato 2008's.
    # -------------------------------------------------------------------
    mwRatio <- 383.5 / 497.63

    ka   <- exp(lka + etalka)
    # CL/F (L/h). Zhang 2026 Equation (6). CRCL is weight-normalized and in
    # mL/min/kg, the unit in which Sato 2008's slope 0.104 is defined, so
    # the 1000/60 factor printed in Equation (6) (which converts a L/h/kg
    # input) is not required here.
    cl   <- mwRatio * (exp(lcl) + e_crcl_cl * CRCL) * WT * exp(etalcl)
    # V/F (L). Zhang 2026 Equation (7); AGE in years.
    vc   <- mwRatio * exp(lvc + etalvc) * AGE^e_age_vc * WT
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Sato 2008 Table 4 Tlag; NM-TRAN ALAG1 on the absorption compartment.
    alag(depot) <- tlag

    Cc <- central / vc
    # Free (unbound) plasma tebipenem, the driver of the 40% fT > MIC and
    # fAUC0-24/MIC/tau PK/PD targets in Zhang 2026 Methods 2.5 and 2.8.
    Ccfree <- Cc * fu
    Cc ~ add(addSd)
  })
}
