Wallender_2021_piperaquine_qtc <- function() {
  description <- "Concentration-QTc model for piperaquine in Ugandan children, from Wallender 2021. The Bazett-corrected QT interval is a linear function of the time-varying plasma piperaquine concentration on top of an individual pre-dose baseline: QTcB = baseline * exp(eta) + 0.0463 * Cc, so each 100 ng/mL of piperaquine adds 4.6 msec. Fit simultaneously with the pharmacokinetics on the intensive-PK substudy of 32 children with paired electrocardiograms before dosing and 2 to 3 h after the third daily dihydroartemisinin-piperaquine dose at 32 and 104 weeks of age; age, sex and weight were screened as covariates on the QTc parameters and none was retained. The three-compartment piperaquine population PK model of the same paper is embedded so the QTc prediction is driven by simulated concentration. The authors caution that peak concentrations under the WHO 2015 and proposed age-based regimens exceed those used to build this model, so extrapolation to those regimens overpredicts the risk. Sister model files from the same paper: modellib('Wallender_2021_piperaquine') (population PK) and modellib('Wallender_2021_piperaquine_malaria') (incident-malaria hazard)."
  reference <- paste(
    "Wallender E, Ali AM, Hughes E, Kakuru A, Jagannathan P, Muhindo MK,",
    "Opira B, Whalen M, Huang L, Duvalsaint M, Legac J, Kajubi R, Aweeka F,",
    "Dorsey G, Kamya MR, Rosenthal PJ, Savic RM.",
    "Identifying an optimal dihydroartemisinin-piperaquine dosing regimen",
    "for malaria prevention in young Ugandan children.",
    "Nat Commun. 2021;12(1):6714.",
    "doi:10.1038/s41467-021-27051-8. PMC8602248.",
    "Open Access under CC BY 4.0.",
    "The QTc parameters are in Supplementary Table 2 of",
    "41467_2021_27051_MOESM1_ESM.pdf and the equation is in the Results",
    "section 'PK-QTc model'; the embedded PK parameters are in Table 2.",
    "Sister model files from the same paper:",
    "modellib('Wallender_2021_piperaquine') and",
    "modellib('Wallender_2021_piperaquine_malaria').",
    sep = " "
  )
  vignette <- "Wallender_2021_piperaquine"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Acts on the EMBEDDED PK model only (allometric exponent 0.75 on all clearances and 1 on all volumes, referenced to the study-median 8.6 kg), and thereby on QTcB indirectly through piperaquine concentration. Weight was tested as a direct covariate on the QTc parameters and was not retained (Methods, 'PK-QTc analysis').",
      source_name        = "WEIGHT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age).",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months. Acts on the EMBEDDED PK model only, through the Emax maturation term PMA / (PMA + 96) on clearance. Age was tested as a direct covariate on the QTc parameters and was not retained.",
      source_name        = "Post menstrual age"
    ),
    WAZ = list(
      description        = "Weight-for-age z-score against the WHO Child Growth Standards.",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Acts on the EMBEDDED PK model only, through relative oral bioavailability, and therefore reaches QTcB only through piperaquine exposure. This is the mechanism behind Supplementary Figure 6, where predicted maximum QTcB differs between malnourished and better-nourished children.",
      source_name        = "WAZ"
    ),
    SELFADMIN = list(
      description        = "Self-administered dosing-occasion indicator: 1 = the DP course was taken at home without direct observation, 0 = every daily dose of the course was directly observed.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (directly observed therapy)",
      notes              = "Per dosing occasion. Acts on the EMBEDDED PK model only. Every course contributing to this sub-model was directly observed, because the electrocardiograms were paired with the intensive-PK substudy in which all three daily doses were given in clinic, so SELFADMIN = 0 is the setting that reproduces the data this model was fit to.",
      source_name        = "Self-administered DP"
    ),
    OCC = list(
      description        = "Integer-valued occasion index identifying the DP treatment course, for between-occasion variability on relative oral bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "One occasion per three-day DP course; values 1 to 25 are multiplexed inside model() onto etaiov_fdepot_1 .. etaiov_fdepot_25 sharing one variance. OCC = 0, or any value outside 1..25, yields the occasion-free typical value. See modellib('Wallender_2021_piperaquine') for the full annotation.",
      source_name        = "OCC"
    )
  )

  # Covariates the source screened on the QTc parameters but did not retain.
  # Documentation only.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Tested as a covariate on the PK-QTc parameters and not retained: 'Age, sex and weight were tested as covariates for the PK-QTc model' (Methods, added in revision at the reviewer's request; Peer Review File response to the reviewer question on line 178). No coefficient is reported, so none can be encoded."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "32 and 104 weeks of age (the two intensive-PK visits)",
    weight_median  = "8.6 kg (the parent study median, used as the allometric reference)",
    disease_state  = "Healthy Ugandan infants and toddlers enrolled in the intensive pharmacokinetic substudy of an intermittent preventive treatment trial; 22 from the every-12-weeks arm and 10 from the every-4-weeks arm.",
    dose_range     = "Weight-band dosed dihydroartemisinin-piperaquine once daily for three consecutive days, directly observed at both intensive-PK visits.",
    regions        = "Tororo District, Uganda",
    biomarkers     = "QT interval corrected by Bazett's formula (QT divided by the square root of the RR interval). Median pre-drug QTcB 413 msec (range 347-472) and post-drug 424 msec (range 388-482).",
    notes          = "Substudy of NCT02163447. An electrocardiogram was obtained before the first dose of DP at 32 and 104 weeks of age and again 2 to 3 h after the third daily dose. Bazett's correction was chosen because it best corrected for heart rate in this cohort. The PK sampling schedule was venipuncture at 0.5, 1, 2, 3, 4, 6, 8 and 24 h after the third daily dose, then finger-prick at 24 h and at 4, 7, 14 and 21 days. Linear and Emax models were tested for the concentration-QTcB relationship and the linear model was retained. Estimation in NONMEM 7.4."
  )

  ini({
    # ==================================================================
    # EMBEDDED PIPERAQUINE PK -- identical to
    # modellib('Wallender_2021_piperaquine'), Table 2 "Pharmacokinetic
    # parameters". See that file for the full annotation, including the
    # note that the Peer Review File carries a SUPERSEDED earlier
    # revision of Table 2 that must not be used. The capillary-to-venous
    # measurement conversion is deliberately absent: the QTc model is
    # driven by plasma concentration.
    # ==================================================================
    lcl  <- log(867)   ; label("Apparent oral clearance at the reference weight, extrapolated to full enzyme maturation (CL/F, L/day)")     # Table 2 Clearance 867 L/d
    lvc  <- log(592)   ; label("Apparent central volume of distribution at the reference weight (Vc/F, L)")                                 # Table 2 Volume of central compartment 592 L
    lq   <- log(511)   ; label("Apparent intercompartmental clearance to the first peripheral compartment at the reference weight (Q1/F, L/day)")   # Table 2 Intercompartmental clearance 1 = 511 L/d
    lvp  <- log(7240)  ; label("Apparent first peripheral volume of distribution at the reference weight (Vp1/F, L)")                       # Table 2 Volume of peripheral compartment 1 = 7240 L
    lq2  <- log(671)   ; label("Apparent intercompartmental clearance to the second peripheral compartment at the reference weight (Q2/F, L/day)")  # Table 2 Intercompartmental clearance 2 = 671 L/d
    lvp2 <- log(1060)  ; label("Apparent second peripheral volume of distribution at the reference weight (Vp2/F, L)")                      # Table 2 Volume of peripheral compartment 2 = 1060 L
    lmtt <- log(0.045) ; label("Mean transit time of the absorption chain (MTT, day)")                                                      # Table 2 Absorption transit time 0.045 d; ktr = (n+1)/MTT with n = 2

    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on all clearance parameters, referenced to 8.6 kg (unitless)")       # Methods, "Population PK model"
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of body weight on all volume parameters, referenced to 8.6 kg (unitless)")          # Methods, "Population PK model"

    ltm50_cl <- log(96) ; label("Postmenstrual age at which the maturation of piperaquine clearance is half-maximal (weeks)")               # Table 2 theta Post menstrual age EC50 = 96 weeks

    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability at the reference covariates and on a directly observed dosing occasion (F, unitless)")  # Table 2 Relative bioavailability (F) = 1, pre-specified anchor

    e_waz_fdepot       <- 0.113 ; label("Fractional change in relative oral bioavailability per unit of weight-for-age z-score above the cohort median of -0.5 (unitless)")  # Table 2 theta Weight for age z-score = 0.113
    e_selfadmin_fdepot <- 0.397 ; label("Multiplicative effect of a self-administered rather than directly observed dosing occasion on relative oral bioavailability (unitless)")  # Table 2 theta Self-administered DP = 0.397

    # Table 2 Clearance IIV 27.1%; log(0.271^2 + 1) = 0.0708694
    etalcl  ~ 0.0708694
    # Table 2 Volume of central compartment IIV 32.8%; log(0.328^2 + 1) = 0.1021811
    etalvc  ~ 0.1021811
    # Table 2 Absorption transit time IIV 43.2%; log(0.432^2 + 1) = 0.1711123
    etalmtt ~ 0.1711123

    # Table 2 theta Between occasion variability 66.9%; log(0.669^2 + 1) = 0.3698801.
    etaiov_fdepot_1  ~ 0.3698801
    etaiov_fdepot_2  ~ fixed(0.3698801)
    etaiov_fdepot_3  ~ fixed(0.3698801)
    etaiov_fdepot_4  ~ fixed(0.3698801)
    etaiov_fdepot_5  ~ fixed(0.3698801)
    etaiov_fdepot_6  ~ fixed(0.3698801)
    etaiov_fdepot_7  ~ fixed(0.3698801)
    etaiov_fdepot_8  ~ fixed(0.3698801)
    etaiov_fdepot_9  ~ fixed(0.3698801)
    etaiov_fdepot_10 ~ fixed(0.3698801)
    etaiov_fdepot_11 ~ fixed(0.3698801)
    etaiov_fdepot_12 ~ fixed(0.3698801)
    etaiov_fdepot_13 ~ fixed(0.3698801)
    etaiov_fdepot_14 ~ fixed(0.3698801)
    etaiov_fdepot_15 ~ fixed(0.3698801)
    etaiov_fdepot_16 ~ fixed(0.3698801)
    etaiov_fdepot_17 ~ fixed(0.3698801)
    etaiov_fdepot_18 ~ fixed(0.3698801)
    etaiov_fdepot_19 ~ fixed(0.3698801)
    etaiov_fdepot_20 ~ fixed(0.3698801)
    etaiov_fdepot_21 ~ fixed(0.3698801)
    etaiov_fdepot_22 ~ fixed(0.3698801)
    etaiov_fdepot_23 ~ fixed(0.3698801)
    etaiov_fdepot_24 ~ fixed(0.3698801)
    etaiov_fdepot_25 ~ fixed(0.3698801)

    # ==================================================================
    # CONCENTRATION-QTcB MODEL -- Supplementary Table 2 and the Results
    # section "PK-QTc model":
    #   QTcB = modeled baseline QTcB + [PPQ] * slope
    #
    # UNITS OF THE SLOPE. Supplementary Table 2 labels the row
    # "theta Slope /1000 (msec/ng/mL)" with a value of .0463, and the
    # body text writes the equation as "+ [PQ] x 0.046/1000". Both are
    # NONMEM estimation-scaling notation: the model's native
    # concentration unit is mg/L, so 46.3 msec per mg/L is 0.0463 msec
    # per ng/mL. The body text pins the effective scale unambiguously --
    # "each 100 ng/mL increase in PPQ concentration was associated with a
    # 4.6 msec increase in the QTcB" -- and 100 * 0.0463 = 4.63 msec.
    # The tabulated .0463 is therefore used directly against a
    # concentration in ng/mL.
    # ==================================================================
    le0_qtc <- log(410) ; label("Pre-dose Bazett-corrected QT interval (msec)")                                        # Supplementary Table 2 Pre-drug QTcB = 410 msec (0.5% RSE)

    e_ppq_qtc <- 0.0463 ; label("Increase in the Bazett-corrected QT interval per ng/mL of plasma piperaquine (msec per ng/mL)")  # Supplementary Table 2 theta Slope = .0463 (45% RSE); 100 ng/mL gives 4.63 msec

    # Supplementary Table 2 interindividual variability on the pre-drug
    # QTcB of 1.4% (41% RSE); log(0.014^2 + 1) = 0.000196
    etale0_qtc ~ 0.000196

    # Supplementary Table 2 reports no residual error for the QTcB
    # observation. This placeholder additive residual is attached so the
    # nlmixr2 likelihood machinery accepts the model for forward
    # simulation; it is NOT from the source and must not be read as a
    # published estimate. See the vignette Assumptions and deviations
    # section.
    addSd <- fixed(1) ; label("Additive residual standard deviation on QTcB (msec); not paper-derived")
  })

  model({
    # ---- embedded PK ------------------------------------------------
    allom_cl <- (WT / 8.6)^e_wt_cl
    allom_v  <- (WT / 8.6)^e_wt_vc
    mat_cl  <- PAGE / (PAGE + exp(ltm50_cl))

    iov_fdepot <-
      (OCC ==  1) * etaiov_fdepot_1  + (OCC ==  2) * etaiov_fdepot_2  +
      (OCC ==  3) * etaiov_fdepot_3  + (OCC ==  4) * etaiov_fdepot_4  +
      (OCC ==  5) * etaiov_fdepot_5  + (OCC ==  6) * etaiov_fdepot_6  +
      (OCC ==  7) * etaiov_fdepot_7  + (OCC ==  8) * etaiov_fdepot_8  +
      (OCC ==  9) * etaiov_fdepot_9  + (OCC == 10) * etaiov_fdepot_10 +
      (OCC == 11) * etaiov_fdepot_11 + (OCC == 12) * etaiov_fdepot_12 +
      (OCC == 13) * etaiov_fdepot_13 + (OCC == 14) * etaiov_fdepot_14 +
      (OCC == 15) * etaiov_fdepot_15 + (OCC == 16) * etaiov_fdepot_16 +
      (OCC == 17) * etaiov_fdepot_17 + (OCC == 18) * etaiov_fdepot_18 +
      (OCC == 19) * etaiov_fdepot_19 + (OCC == 20) * etaiov_fdepot_20 +
      (OCC == 21) * etaiov_fdepot_21 + (OCC == 22) * etaiov_fdepot_22 +
      (OCC == 23) * etaiov_fdepot_23 + (OCC == 24) * etaiov_fdepot_24 +
      (OCC == 25) * etaiov_fdepot_25

    cl  <- exp(lcl + etalcl) * allom_cl * mat_cl
    vc  <- exp(lvc + etalvc) * allom_v
    q   <- exp(lq)   * allom_cl
    vp  <- exp(lvp)  * allom_v
    q2  <- exp(lq2)  * allom_cl
    vp2 <- exp(lvp2) * allom_v

    mtt <- exp(lmtt + etalmtt)

    fdepot <- exp(lfdepot + iov_fdepot) *
      (1 + e_waz_fdepot * (WAZ - (-0.5))) *
      e_selfadmin_fdepot^SELFADMIN

    ktr <- 3 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # Plasma piperaquine concentration in ng/mL.
    Cc <- 1000 * central / vc

    # ---- concentration-QTcB relationship ----------------------------
    # Linear, no threshold and no Emax: an Emax form was tested and the
    # linear model was retained (Methods, 'PK-QTc analysis'). The max()
    # is a numerical guard against a solver undershoot only.
    e0_qtc <- exp(le0_qtc + etale0_qtc)

    QTcB <- e0_qtc + e_ppq_qtc * max(Cc, 0)

    # Change from the individual's own pre-dose baseline, the quantity
    # the paper reports as the drug effect (4.6 msec per 100 ng/mL).
    dQTcB <- QTcB - e0_qtc

    QTcB ~ add(addSd)
  })
}
