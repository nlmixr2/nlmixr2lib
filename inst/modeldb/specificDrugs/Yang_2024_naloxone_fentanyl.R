Yang_2024_naloxone_fentanyl <- function() {
  description <- paste(
    "QSP. Mechanistic PK-PD model of fentanyl-induced respiratory",
    "depression and its reversal by intramuscular naloxone auto-injector",
    "(NAI). Couples a two-compartment fentanyl PK model (the parent,",
    "unsuffixed) to the Yang 2024 naloxone PK model (depot plus three",
    "transit absorption compartments and a two-compartment disposition,",
    "all _naloxone-suffixed). Both drugs equilibrate with their own",
    "biophase (effect) compartment. Unlike the buprenorphine, morphine",
    "and carfentanil models in this family, fentanyl carries no explicit",
    "receptor-occupancy state: because fentanyl binds and dissociates",
    "rapidly (dissociation half-life about 6.8 min, complete within an",
    "hour) the receptor layer collapses to an instantaneous fractional",
    "Emax function of effect-site concentration, and naloxone acts as a",
    "classical competitive antagonist that right-shifts the fentanyl",
    "EC50 by the factor (1 + Ce_naloxone / EC50_naloxone). Fentanyl",
    "PK-PD parameters are inherited unchanged from Yassen 2007; the",
    "naloxone PK layer is Yang 2024's own fit."
  )
  reference <- paste(
    "Yang TE, Del Bene F, Lavezzi SM, Iavarone L, Zhang J, Kim J,",
    "Gjurich B, Kessler C. Mechanistic pharmacokinetic-pharmacodynamic",
    "modeling and simulations of naloxone auto-injector 10 mg reversal",
    "of opioid-induced respiratory depression.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1722-1733.",
    "doi:10.1002/psp4.13215. PMCID PMC11494827.",
    "Parameter source: Appendix S1 Table S7 (Fentanyl Population PK-PD",
    "Parameter Estimates) and Appendix S1 example NONMEM control stream",
    "(C) Fentanyl, whose $THETA / $OMEGA blocks and $DES / $ERROR",
    "equations this file reproduces one-for-one.",
    "Fentanyl PK-PD parameters originate from Yassen A, Olofsen E,",
    "Romberg R, Sarton E, Teppema L, Danhof M, Dahan A. Mechanism-based",
    "PK/PD modeling of the respiratory depressant effect of",
    "buprenorphine and fentanyl in healthy volunteers.",
    "Clin Pharmacol Ther. 2007;81(1):50-58. doi:10.1038/sj.clpt.6100025.",
    "The naloxone EC50 of 0.6021768 ng/mL is taken from Olofsen E, et",
    "al. Anesthesiology. 2010;112(6):1417-1427",
    "(doi:10.1097/ALN.0b013e3181d5e29d), as stated in Appendix S1 and",
    "hardcoded in the $ERROR block of control stream (C).",
    "Naloxone PK parameters are Yang 2024 Appendix S1 Table S4.",
    "The competitive-antagonism form follows Salahudeen MS, Nishtala PS.",
    "An overview of pharmacodynamic modelling, ligand-binding approach",
    "and its application in clinical practice. Saudi Pharm J.",
    "2017;25(2):165-175. doi:10.1016/j.jsps.2016.07.002."
  )
  vignette <- "Yang_2024_naloxone_opioid_reversal"
  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ng/mL"
  )
  dosing <- c("central", "depot_naloxone")

  compartmentData <- list(
    central              = list(analyte = "fentanyl", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "fentanyl", units = "ug", specimen = "plasma", verified = TRUE),
    effect               = list(analyte = "fentanyl biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    depot_naloxone       = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit1_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit2_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit3_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    central_naloxone     = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_naloxone = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    effect_naloxone      = list(analyte = "naloxone biophase (effect-site) concentration", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Acts only on naloxone apparent clearance:",
        "CL/F = 3.26 * (WT/70)^0.538 (Yang 2024 Table S4). The reference",
        "weight of 70 kg is explicit in control stream (C)",
        "('WT= 70 ; kg'). The fentanyl PK parameters carry no weight",
        "scaling in Yassen 2007."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Healthy adults (naloxone layer 23-54 years)",
    weight_range   = "Simulations fix WT = 70 kg; naloxone PK data spanned 57.2-100.2 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy adult volunteers. The fentanyl PK-PD layer comes from",
      "the Yassen 2007 healthy-volunteer study of fentanyl-induced",
      "respiratory depression; the naloxone PK layer is Yang 2024's own",
      "auto-injector population PK in 48 healthy adults."
    ),
    dose_range     = paste(
      "Yang 2024 simulated IV fentanyl at 2.1, 23.1 and 44.1 ug/kg (the",
      "lowest corresponding to 50 percent ventilation suppression, the",
      "highest 21-fold higher and representing a lethal dose), each",
      "with naloxone auto-injector 2 mg or 10 mg IM/SC given once",
      "ventilation had fallen 60 percent from baseline; and a",
      "prophylactic arm giving NAI 10 mg at 5, 15, 30 or 60 min BEFORE",
      "the 44.1 ug/kg fentanyl dose."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Yang 2024 validated this constructed model by digitising the",
      "ventilation-time profiles for fentanyl 0.15 and 0.3 mg/70 kg",
      "from Figure 1 of Yassen 2007 and confirming that the observed",
      "points fell inside the 90 percent prediction interval of 1000",
      "simulations (Appendix S1 Figure S9).",
      "Baseline ventilation V0 was carried as a per-subject data column",
      "(E0) in the control stream rather than as a THETA; it is encoded",
      "here as the canonical le0 parameter with the Table S7 typical",
      "value and IIV so that the model is self-contained.",
      "The fentanyl PK parameters reconcile with the Appendix S1",
      "carfentanil-scaling text, which states that human fentanyl",
      "clearance and total volume are 58.8 L/h and 170 L: CL = 0.98",
      "L/min = 58.8 L/h and V1 + V2 = 19.5 + 150 = 169.5 L."
    )
  )

  ini({
    # Every value below is held FIX in Appendix S1 control stream (C):
    # the fentanyl PK-PD block is inherited unchanged from Yassen 2007
    # and the naloxone PK block from Yang 2024's own Table S4, and
    # neither was re-fitted for this analysis.

    # --- Naloxone PK (Yang 2024 Table S4; control stream (C) THETA 1-6)
    lktr_naloxone <- fixed(log(0.696))
    label("Naloxone transit absorption rate constant (KTR, 1/min)")            # (C) THETA 1 = 0.696; Table S4
    lcl_naloxone <- fixed(log(3.26))
    label("Naloxone apparent clearance at 70 kg (CL/F, L/min)")                # (C) THETA 2 = 3.26; Table S4
    e_wt_cl_naloxone <- fixed(0.538)
    label("Allometric exponent of body weight on naloxone CL/F (unitless)")    # (C) THETA 3 = 0.538; Table S4
    lvc_naloxone <- fixed(log(404))
    label("Naloxone apparent central volume (V2/F, L)")                        # (C) THETA 4 = 404; Table S4
    lq_naloxone <- fixed(log(0.0847))
    label("Naloxone apparent intercompartmental clearance (Q/F, L/min)")       # (C) THETA 5 = 0.0847; Table S4
    lvp_naloxone <- fixed(log(81.8))
    label("Naloxone apparent peripheral volume (V3/F, L)")                     # (C) THETA 6 = 81.8; Table S4

    # --- Fentanyl PK (Table S7; control stream (C) THETA 7-10)
    lcl <- fixed(log(0.98))
    label("Fentanyl clearance (CL, L/min)")                                    # (C) THETA 7 = 0.98; Table S7
    lvc <- fixed(log(19.5))
    label("Fentanyl central volume (V1, L)")                                   # (C) THETA 8 = 19.5; Table S7
    lq <- fixed(log(3.51))
    label("Fentanyl intercompartmental clearance (Q2, L/min)")                 # (C) THETA 9 = 3.51; Table S7
    lvp <- fixed(log(150))
    label("Fentanyl peripheral volume (V2, L)")                                # (C) THETA 10 = 150; Table S7

    # --- Biophase equilibration (equation 1)
    lke0_naloxone <- fixed(log(0.106))
    label("Naloxone biophase equilibration rate constant (ke0, 1/min)")        # (C) THETA 11 = 0.106; Table S7
    lke0 <- fixed(log(0.0422))
    label("Fentanyl biophase equilibration rate constant (ke0, 1/min)")        # (C) THETA 12 = 0.0422; Table S7

    # --- Fractional Emax PD with competitive naloxone antagonism
    # (equations 11-12)
    lec50 <- fixed(log(1.140))
    label("Fentanyl half-maximal effect-site concentration (EC50, ng/mL)")     # (C) THETA 13 = 1.140; Table S7
    lalpha <- fixed(log(0.91))
    label("Fentanyl intrinsic activity (alpha, unitless 0-1)")                 # (C) THETA 14 = 0.91; Table S7
    lhill <- fixed(log(2.68))
    label("Fentanyl sigmoidicity slope parameter (n, unitless)")               # (C) THETA 15 = 2.68; Table S7
    lec50_naloxone <- fixed(log(0.6021768))
    label("Naloxone half-maximal effect-site concentration, from Olofsen 2010 (EC50, ng/mL)")  # (C) $ERROR EC50N = 0.6021768; Table S7 rounds to 0.602
    le0 <- fixed(log(20.2))
    label("Baseline minute ventilation (V0, L/min)")                           # Table S7: V0 20.2 L/min (carried as the E0 data column in control stream (C))

    # --- IIV. Control stream (C) $OMEGA entries 1-8, in order. Each
    # equals the square of the Table S7 "IIV (%)" column (e.g. fentanyl
    # CL 67.7% -> 0.458329 = 0.677^2). Q2 and V2 have no IIV, matching
    # the "-" entries in Table S7.
    etalktr_naloxone ~ fixed(0.111)
    etalcl_naloxone ~ fixed(0.0129)
    etalvc_naloxone ~ fixed(0.0658)
    etalcl ~ fixed(0.458329)
    etalvc ~ fixed(0.329476)
    etalke0 ~ fixed(0.142884)
    etalec50 ~ fixed(0.035721)
    etalalpha ~ fixed(0.064516)
    # Baseline ventilation IIV. Control stream (C) has no $OMEGA entry
    # for it because E0 arrived as a per-subject data column already
    # carrying the variability; Table S7 reports V0 IIV = 26.7%, giving
    # the variance 0.267^2 = 0.071289.
    etale0 ~ fixed(0.071289)

    propSd <- fixed(0.217)
    label("Fentanyl proportional residual error SD (fraction)")                # Table S7: proportional error 21.7%
    addSd <- fixed(6.1)
    label("Ventilation additive residual error SD (L/min)")                    # (C) THETA 16 = 6.1; Table S7
  })

  model({
    # ---- Naloxone individual parameters
    ktr_naloxone <- exp(lktr_naloxone + etalktr_naloxone)
    cl_naloxone <- exp(lcl_naloxone + etalcl_naloxone) * (WT / 70)^e_wt_cl_naloxone
    vc_naloxone <- exp(lvc_naloxone + etalvc_naloxone)
    q_naloxone <- exp(lq_naloxone)
    vp_naloxone <- exp(lvp_naloxone)
    ke0_naloxone <- exp(lke0_naloxone)
    ec50_naloxone <- exp(lec50_naloxone)

    # ---- Fentanyl individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0 + etalke0)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)
    e0 <- exp(le0 + etale0)
    # Control stream (C) $PK caps the intrinsic activity at 1
    # ("IF(ALPHA.GT.1) ALPHA = 1 ;; ALPHA is between 0 and 1"), which
    # matters because the 25.4% IIV on a typical value of 0.91 puts a
    # sizeable share of the population above 1. Without the cap those
    # subjects would be driven to negative ventilation at high
    # effect-site concentrations.
    alpha_uncapped <- exp(lalpha + etalalpha)
    alpha <- alpha_uncapped - (alpha_uncapped > 1) * (alpha_uncapped - 1)

    # ---- Naloxone PK, control stream (C) $DES DADT(1)-DADT(6).
    d/dt(depot_naloxone) <- -ktr_naloxone * depot_naloxone
    d/dt(transit1_naloxone) <- ktr_naloxone * depot_naloxone - ktr_naloxone * transit1_naloxone
    d/dt(transit2_naloxone) <- ktr_naloxone * transit1_naloxone - ktr_naloxone * transit2_naloxone
    d/dt(transit3_naloxone) <- ktr_naloxone * transit2_naloxone - ktr_naloxone * transit3_naloxone
    d/dt(central_naloxone) <- ktr_naloxone * transit3_naloxone -
      (cl_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vc_naloxone) * central_naloxone +
      (q_naloxone / vp_naloxone) * peripheral1_naloxone
    d/dt(peripheral1_naloxone) <- (q_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vp_naloxone) * peripheral1_naloxone
    d/dt(effect_naloxone) <- ke0_naloxone * (central_naloxone / vc_naloxone - effect_naloxone)

    # ---- Fentanyl PK, control stream (C) $DES DADT(8)-DADT(9).
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(effect) <- ke0 * (central / vc - effect)

    # ---- Outputs
    Cc <- central / vc
    Cnal <- central_naloxone / vc_naloxone
    # Fractional Emax with competitive naloxone antagonism, equation 12 /
    # control stream (C) $ERROR:
    #   NALEFF = 1 + Ce_nal/EC50_nal
    #   E = E0*(1 - alpha*Ce_f^n / (Ce_f^n + (EC50_f*NALEFF)^n))
    # Naloxone shifts the fentanyl EC50 to the right without changing
    # the attainable maximum, the signature of competitive antagonism.
    naleff <- 1 + effect_naloxone / ec50_naloxone
    VE <- e0 * (1 - alpha * effect^hill /
      (effect^hill + (ec50 * naleff)^hill))
    ERATIO <- VE / e0

    Cc ~ prop(propSd)
    VE ~ add(addSd)
  })
}
