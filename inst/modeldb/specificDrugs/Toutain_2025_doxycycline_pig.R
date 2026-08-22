Toutain_2025_doxycycline_pig <- function() {
  description <- "Preclinical (pig). Three-compartment population PK meta-analysis of doxycycline in pigs, parameterised per kg body weight, with a body-weight power model on the three clearances and the deep peripheral volume and four oral administration modalities (medicated feed under field or laboratory conditions, drinking water, stomach tube) each carrying its own absorption rate constant, bioavailability and residual error (Toutain 2025)"
  reference <- "Toutain PL, Bousquet-Melou A, Ferran AA, Roques BB, del Castillo JRE, Lees P, Croubels S, Bousquet E, Pelligand L. Pharmacokinetic-pharmacodynamic cutoff values for doxycycline in pigs to support the establishment of clinical breakpoints for antimicrobial susceptibility testing. J Vet Pharmacol Ther. 2025;48(4):300-317. doi:10.1111/jvp.13511"
  paper_specific_residual_sds <- c(
    "propSdIv", "addSdIv",
    "propSdFeedField", "addSdFeedField",
    "propSdFeedLab", "addSdFeedLab",
    "propSdWater", "addSdWater",
    "propSdTube", "addSdTube"
  )
  vignette <- "Toutain_2025_doxycycline_pig"
  # Every structural parameter is normalised to body weight: doses are given in
  # mg/kg, compartment amounts are carried in mg/kg, volumes in L/kg and
  # clearances in L/kg/h, exactly as published. central/vc is therefore already
  # mg/L = ug/mL and no unit conversion is applied.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "doxycycline", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power model (WT/50)^theta on the three clearances and on the deep peripheral volume, per Toutain 2025 Equation 1; the 50 kg scaling factor is the rounded observed median body weight (actual median 44.15 kg, range 8.5-100.6 kg). Because every parameter is already expressed per kg body weight, the absence of a WT term on Vc and V2 means the whole-animal Vc and V2 are exactly proportional to body weight, while whole-animal CL scales as WT^1.299.",
      source_name        = "BW"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral administration: medicated feed, drinking water or stomach tube)",
      notes              = "Per-dose-record indicator (1 = single IV dose through an indwelling catheter, 0 = one of the four oral modalities). Selects the IV residual-error pair (propSdIv, addSdIv) and switches off oral absorption (ka and fdepot both collapse to 0). IV doses must additionally be placed in the central compartment via cmt = 'central' on the event record; oral doses go to cmt = 'depot'. Toutain 2025 Methods 2.3.2 block A of the Appendix S3 Phoenix script.",
      source_name        = "route of administration (IV)"
    ),
    FORM_DOX_FEED = list(
      description        = "Indicator for doxycycline administered mixed in medicated feed",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (doxycycline given as an aqueous solution, i.e. in drinking water or by stomach tube)",
      notes              = "Per-dose-record indicator distinguishing the two feed sub-models (submodel_1 and submodel_2 of Toutain 2025 Methods 2.3.2) from the two aqueous-solution sub-models (submodel_3 and submodel_4). Ignored when ROUTE_IV = 1. Combined with STUDY_TLS to pick between the field-condition and the laboratory-condition feed parameters.",
      source_name        = "route of administration (oral feed)"
    ),
    STUDY_TLS = list(
      description        = "Indicator for the TLS field trial (doxycycline in feed under on-farm field conditions)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the laboratory-condition feed trials AFSSA, BIOEQ, PARADOX, Company 9203 and Company 9204)",
      notes              = "Per-dose-record indicator for the 215-pig TLS field trial (del Castillo 2006), the single largest cohort in the meta-analysis and the only one dosed in feed under field conditions with group access to medicated meal. Only meaningful when FORM_DOX_FEED = 1. Toutain 2025 Table 1 and Table 6 ('Trial TLS (feed, field conditions)' vs 'Trials AFSSA, BIOEQ, PARADOX, Company 9203 and Company 9204 (feed, laboratory conditions)'). The paper attributes the much larger bioavailability variability under field conditions to competition between animals for access to feed.",
      source_name        = "Trial ID = TLS"
    ),
    ROUTE_NGT = list(
      description        = "Indicator for an oral solution delivered by gastric tube",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral solution taken spontaneously in drinking water)",
      notes              = "Per-dose-record indicator separating the stomach-tubing sub-model (submodel_4: trials 104NL, 3205NL, GHENT) from the spontaneous drinking-water sub-model (submodel_3: trials KING_NL, Bea). Only meaningful when ROUTE_IV = 0 and FORM_DOX_FEED = 0. The canonical ROUTE_NGT column covers tube-delivered vs voluntarily-ingested enteral dosing; Toutain 2025 uses an oro-gastric ('stomach tubing') rather than a nasogastric tube.",
      source_name        = "SOL (ST)"
    )
  )

  population <- list(
    species          = "pig (Sus scrofa domesticus; piglets to finishing pigs)",
    n_subjects       = 300L,
    n_datasets       = 380L,
    n_observations   = 3295L,
    n_studies        = 11L,
    weight_range     = "8.5-100.6 kg",
    weight_median    = "44.15 kg",
    sex_female_pct   = 49.3,
    disease_state    = "Mostly healthy production pigs; in the 215-pig TLS field trial 146 pigs were assessed as healthy and 66 as sick. Health status was tested as a covariate on bioavailability and was not significant, so it is absent from the final model.",
    dose_range       = "IV 5-10.5 mg/kg single dose (57 data sets); oral in medicated feed 2.4-13.3 mg/kg per administration (265 data sets, single dose to 15 administrations at 12 h intervals); oral solution 8.68-10.59 mg/kg by stomach tube or in drinking water (58 data sets)",
    regions          = "France (6 trials), the Netherlands (4 trials), Belgium (1 trial)",
    sampling         = "3295 plasma samples analysed by HPLC-UV or HPLC-MS/MS with LLOQ 0.025-0.2 ug/mL; 2% of samples were below LLOQ and were discarded for the final analysis. Rich sampling (8-15 samples per pig over 12-48 h) in the laboratory trials; sparse sampling (one pre-dose plus 5-6 post-dose samples) in the 215-pig TLS field trial.",
    reference_subject = "50 kg body weight (rounded observed median), the scaling factor of the body-weight power model.",
    notes            = "VetCAST meta-analysis of raw individual plasma concentration-time data from 11 trials (5 published, 6 unpublished marketing-authorisation trials); the same pigs were dosed on 2 or 3 occasions in the laboratory trials and each data set was treated as a separate individual, giving 380 analyzable data sets from 300 pigs. Doses are expressed as doxycycline base. Estimation was FOCE ELS in Phoenix NLME 8.3 with standard errors from a QRPEM re-run; the model was used to compute PK/PD cutoffs from fAUC/MIC with an unbound fraction of 0.31 (Portugal 2023)."
  )

  ini({
    # ---- Systemic disposition (shared by every route) ---------------------
    # Toutain 2025 Table 4 'Bootstraps (With BW covariate)' mean column; the
    # full-precision values below are the frozen thetas of the final Phoenix
    # run reproduced in Appendix S3 block A. They were held at their IV-only
    # estimates when the merged IV + oral data were fitted (Toutain 2025
    # Methods 2.3.2: "the fixed effects (i.e., the thetas) of the IV route
    # were set at their optimal values obtained with the adjustment of the IV
    # data alone"), hence fixed() here.
    lvc  <- fixed(log(0.191853096058459)); label("Central volume of distribution Vc (L/kg)")                            # Toutain 2025 Table 4 (0.192) / Appendix S3 tvVc(freeze)
    lvp  <- fixed(log(0.594964536703773)); label("Superficial peripheral volume of distribution V2 (L/kg)")             # Toutain 2025 Table 4 (0.595) / Appendix S3 tvV2(freeze)
    lvp2 <- fixed(log(0.535786748090738)); label("Deep peripheral volume of distribution V3 at 50 kg (L/kg)")           # Toutain 2025 Table 4 (0.536) / Appendix S3 tvV3(freeze)
    lcl  <- fixed(log(0.259372338153985)); label("Plasma clearance at 50 kg (L/kg/h)")                                  # Toutain 2025 Table 4 (0.259) / Appendix S3 tvClearance(freeze)
    lq   <- fixed(log(1.17871023400196));  label("Distributional clearance to peripheral1 Cl2 at 50 kg (L/kg/h)")       # Toutain 2025 Table 4 (1.179) / Appendix S3 tvCld2(freeze)
    lq2  <- fixed(log(0.0723295319908054)); label("Distributional clearance to peripheral2 Cl3 at 50 kg (L/kg/h)")      # Toutain 2025 Table 4 (0.072) / Appendix S3 tvCld3(freeze)

    # ---- Body-weight power model, Toutain 2025 Equation 1 ------------------
    # parameter = theta_pop * (BW / 50)^theta_BW, applied to the per-kg values
    # above. Only Cl, Cl2, Cl3 and V3 carried a significant BW effect
    # (Toutain 2025 Results 3.2); Vc and V2 did not.
    e_wt_cl  <- fixed(0.299259388973256);  label("Power exponent on (WT/50) for per-kg plasma clearance (unitless)")            # Toutain 2025 Table 4 Cov BWCl (0.299) / Appendix S3 dClearancedBW(freeze)
    e_wt_q   <- fixed(-0.224035667026942); label("Power exponent on (WT/50) for per-kg distributional clearance Cl2 (unitless)") # Toutain 2025 Table 4 Cov BWCl2 (-0.224) / Appendix S3 dCld2dBW(freeze)
    e_wt_q2  <- fixed(-0.544270965800396); label("Power exponent on (WT/50) for per-kg distributional clearance Cl3 (unitless)") # Toutain 2025 Table 4 Cov BWCl3 (-0.544) / Appendix S3 dCld3dBW(freeze)
    e_wt_vp2 <- fixed(0.375752281395802);  label("Power exponent on (WT/50) for per-kg deep peripheral volume V3 (unitless)")    # Toutain 2025 Table 4 CovBWV3 (0.376) / Appendix S3 dV3BW(freeze)

    # ---- Oral absorption, one stratum per administration modality ---------
    # Toutain 2025 Table 6 and Appendix S3 blocks B-E. Every stratum carries
    # its own suffix; there is no bare lka / lfdepot.
    lka_feedfield <- log(0.0719834777737194); label("Absorption rate constant, doxycycline in feed under field conditions (1/h)")   # Toutain 2025 Table 6 KaOR_FEED_TLS (0.072)
    lka_feedlab   <- log(0.143923093064807);  label("Absorption rate constant, doxycycline in feed under laboratory conditions (1/h)") # Toutain 2025 Table 6 KaOR_FEED_OTHERS (0.144)
    lka_water     <- log(0.688749514620522);  label("Absorption rate constant, doxycycline solution in drinking water (1/h)")       # Toutain 2025 Table 6 KaOR_SOL_DW (0.689)
    lka_tube      <- log(0.724820162038053);  label("Absorption rate constant, doxycycline solution by stomach tube (1/h)")         # Toutain 2025 Table 6 KaOR_SOLTUBING (0.725)

    lfdepot_feedfield <- log(0.500572676610548); label("Absolute bioavailability, doxycycline in feed under field conditions (fraction)")   # Toutain 2025 Table 6 F_FEED_TLS (0.501)
    lfdepot_feedlab   <- log(0.339942298566399); label("Absolute bioavailability, doxycycline in feed under laboratory conditions (fraction)") # Toutain 2025 Table 6 F_FEED_OTHERS (0.340)
    lfdepot_water     <- log(0.307406622217815); label("Absolute bioavailability, doxycycline solution in drinking water (fraction)")       # Toutain 2025 Table 6 F_SOL_DW (0.307)
    lfdepot_tube      <- log(0.258107364278074); label("Absolute bioavailability, doxycycline solution by stomach tube (fraction)")         # Toutain 2025 Table 6 F_SOLTUBING (0.258)

    # ---- Between-subject variability --------------------------------------
    # Full 6x6 covariance block on the disposition parameters, re-estimated in
    # the merged IV + oral run (Toutain 2025 Methods 2.3.2: "all random
    # effects (OMEGA) were re-estimated in this final run"). Reproduced from
    # the Appendix S3 ranef(block(nVc, nClearance, nV2, nV3, nCld2, nCld3));
    # the diagonal reproduces the Table 6 BSV% via CV = sqrt(exp(omega^2) - 1)
    # (107.8, 27.1, 63.4, 47.5, 80.0 and 246.5% respectively).
    etalvc + etalcl + etalvp + etalvp2 + etalq + etalq2 ~ c(
      0.77092389,
      0.087723279, 0.070991386,
      -0.29277707, 0.04116736, 0.33771914,
      -0.082136838, 0.038806937, 0.071118333, 0.20311218,
      0.38993596, 0.070765694, -0.026903024, 0.1295456, 0.49517412,
      -0.71515683, -0.014934437, 0.18413379, 0.4896846, -0.13000376, 1.9564546
    )  # Toutain 2025 Appendix S3 block A ranef(block(...)); order nVc, nClearance, nV2, nV3, nCld2, nCld3

    # Per-modality 2x2 blocks on absorption rate and bioavailability
    # (Appendix S3 blocks B-E). The diagonals reproduce the Table 6 BSV%.
    # For the field-condition feed stratum the variance of the bioavailability
    # eta is taken from the published BSV of 84.8% (Toutain 2025 Abstract,
    # Table 6 and Discussion) as log(1 + 0.848^2) = 0.541803 rather than from
    # the deposited script, whose 0.671391 implies 97.8%; the covariance is
    # rescaled to preserve the script's correlation of 0.2268. See the
    # vignette Errata.
    etalka_feedfield + etalfdepot_feedfield ~ c(
      0.028277318,
      0.028075540, 0.541803200
    )  # Toutain 2025 Table 6: BSV Ka 16.9%, BSV F 84.8%; correlation 0.2268 from Appendix S3
    etalka_feedlab + etalfdepot_feedlab ~ c(
      0.23420791,
      -0.14655547, 0.12474181
    )  # Toutain 2025 Appendix S3 block C; Table 6 BSV Ka 51.4%, BSV F 36.4%
    etalka_water + etalfdepot_water ~ c(
      0.073474253,
      0.050621961, 0.111363700
    )  # Toutain 2025 Appendix S3 block E; Table 6 BSV Ka 27.6%, BSV F 34.3%
    etalka_tube + etalfdepot_tube ~ c(
      0.25547940,
      -0.05638097, 0.08228053
    )  # Toutain 2025 Appendix S3 block D; Table 6 BSV Ka 54.0%, BSV F 29.3%

    # ---- Residual error, one additive + proportional pair per stratum -----
    # Phoenix writes the combined model as
    #   CEps * sqrt(1 + C^2 * (CMultStdev / sigma())^2)
    # whose standard deviation is sqrt(sigma^2 + C^2 * CMultStdev^2), i.e. an
    # additive SD of sigma (the stdev rows of Table 6) plus a proportional SD
    # of CMultStdev. The IV proportional term was frozen with the other IV
    # thetas; every other residual parameter was estimated in the final run.
    propSdIv        <- fixed(0.139034481021228); label("Proportional residual error, IV route (fraction)")                          # Toutain 2025 Table 6 IV block (0.139) / Appendix S3 tvCMultStdevIV(freeze)
    addSdIv         <- 0.0128320457812293;       label("Additive residual error, IV route (ug/mL)")                                 # Toutain 2025 Table 6 stdev0 IV (0.013)
    propSdFeedField <- 0.228343287097325;        label("Proportional residual error, feed under field conditions (fraction)")       # Toutain 2025 Table 6 CMultStdevOR_FEEDTLS (0.228)
    addSdFeedField  <- 0.113176534645855;        label("Additive residual error, feed under field conditions (ug/mL)")              # Toutain 2025 Table 6 stdev1 Trial TLS (0.113)
    propSdFeedLab   <- 0.184140074875283;        label("Proportional residual error, feed under laboratory conditions (fraction)")  # Toutain 2025 Table 6 CMultStdevOR_FEEDOTHERS (0.184)
    addSdFeedLab    <- 0.0189106141804147;       label("Additive residual error, feed under laboratory conditions (ug/mL)")         # Toutain 2025 Table 6 stdev2 trials FEED_OTHERS (0.019)
    propSdWater     <- 0.29328492186821;         label("Proportional residual error, solution in drinking water (fraction)")        # Toutain 2025 Table 6 CMultStdevOR_SOL_DW (0.293)
    addSdWater      <- 0.00571571651736011;      label("Additive residual error, solution in drinking water (ug/mL)")               # Toutain 2025 Table 6 stdev4 Trials drinking water (0.006)
    propSdTube      <- 0.274744873698196;        label("Proportional residual error, solution by stomach tube (fraction)")          # Toutain 2025 Table 6 CMultStdevOR_SOL_TUBING (0.275)
    addSdTube       <- 0.00255318847588493;      label("Additive residual error, solution by stomach tube (ug/mL)")                 # Toutain 2025 Table 6 stdev3 trials Stomach tubing (0.003)
  })
  model({
    # Mutually exclusive administration-modality selectors. Exactly one of the
    # five is 1 for any dose record; they map onto the five equation blocks of
    # the Appendix S3 Phoenix script (A = IV, B = feed field, C = feed
    # laboratory, D = stomach tube, E = drinking water).
    sOral      <- 1 - ROUTE_IV
    sFeedField <- sOral * FORM_DOX_FEED * STUDY_TLS
    sFeedLab   <- sOral * FORM_DOX_FEED * (1 - STUDY_TLS)
    sTube      <- sOral * (1 - FORM_DOX_FEED) * ROUTE_NGT
    sWater     <- sOral * (1 - FORM_DOX_FEED) * (1 - ROUTE_NGT)

    # Body-weight power model on the per-kg parameters (Equation 1).
    cl  <- exp(lcl  + etalcl)  * (WT / 50)^e_wt_cl
    q   <- exp(lq   + etalq)   * (WT / 50)^e_wt_q
    q2  <- exp(lq2  + etalq2)  * (WT / 50)^e_wt_q2
    vp2 <- exp(lvp2 + etalvp2) * (WT / 50)^e_wt_vp2
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)

    # Modality-specific absorption. Each stratum is given its own simple line
    # so every theta sits in a mu-referenced position, then the active one is
    # selected. Both ka and fdepot collapse to 0 for an IV record, so a dose
    # placed in central is unaffected by the (unused) depot equation.
    kaFeedField <- exp(lka_feedfield + etalka_feedfield)
    kaFeedLab   <- exp(lka_feedlab   + etalka_feedlab)
    kaWater     <- exp(lka_water     + etalka_water)
    kaTube      <- exp(lka_tube      + etalka_tube)
    fFeedField  <- exp(lfdepot_feedfield + etalfdepot_feedfield)
    fFeedLab    <- exp(lfdepot_feedlab   + etalfdepot_feedlab)
    fWater      <- exp(lfdepot_water     + etalfdepot_water)
    fTube       <- exp(lfdepot_tube      + etalfdepot_tube)

    ka <- sFeedField * kaFeedField + sFeedLab * kaFeedLab +
      sWater * kaWater + sTube * kaTube
    fdepot <- sFeedField * fFeedField + sFeedLab * fFeedLab +
      sWater * fWater + sTube * fTube

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # Doses are mg/kg and vc is L/kg, so central / vc is already mg/L = ug/mL.
    Cc <- central / vc

    # Modality-conditional additive-plus-proportional residual error.
    propSd <- ROUTE_IV * propSdIv + sFeedField * propSdFeedField +
      sFeedLab * propSdFeedLab + sWater * propSdWater + sTube * propSdTube
    addSd <- ROUTE_IV * addSdIv + sFeedField * addSdFeedField +
      sFeedLab * addSdFeedLab + sWater * addSdWater + sTube * addSdTube
    Cc ~ add(addSd) + prop(propSd)
  })
}
