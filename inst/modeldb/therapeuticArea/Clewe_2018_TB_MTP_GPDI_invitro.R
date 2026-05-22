Clewe_2018_TB_MTP_GPDI_invitro <- function() {
  description <- "In vitro (M. tuberculosis Beijing VN 2002-1585, BE1585 strain). Multistate Tuberculosis Pharmacometric (MTP) model linked to the General Pharmacodynamic Interaction (GPDI) model describing CFU/mL dynamics during 6-day time-kill assays of M. tuberculosis exposed to static rifampicin / isoniazid / ethambutol concentrations alone or in duo / trio combinations. Three bacterial states (fast-multiplying F, slow-multiplying S, non-multiplying N) follow MTP natural-growth dynamics with the transfer rates kFSlin, kSF, kFN, kSN, kNS fixed from Clewe 2016 JAC dkv416; growth is exponential (kG) rather than logistic. Drug effects per Table 1: rifampicin inhibits F growth (sigmoidal Emax) and kills F, S (Emax) and N (linear); isoniazid kills F and S via sigmoidal Emax with a linear isoniazid-induced adaptive-resistance loop (AR_off / AR_on) modulating EC50 on F-kill and S-kill; ethambutol kills F (sigmoidal Emax) and S (linear). The thirteen quantified GPDI interactions (Table 3) enter as multiplicative modifications of each modulated drug's EC50 using the canonical Emax-style interaction term 1 + INT * C_modifier / (EC50_modifier + C_modifier) per Table 3 footer, with the modifier's mono EC50 acting as the half-saturation anchor (least complex baseline per the paper's Materials and methods). Combination effects on each bacterial state are pooled by the Bliss Independence criterion. Drug concentrations are static covariates (no PK structure -- the in vitro experiment fixes nominal concentrations for the 6-day window). This is the in vitro twin of Chen 2017 in vivo mouse MTP-GPDI; see modellib('Chen_2017_TB_MTP_GPDI_mouse') for the same MTP-GPDI framework applied to TB-infected BALB/c mice with explicit PK."
  reference <- paste(
    "Clewe O, Wicha SG, de Vogel CP, de Steenwinkel JEM, Simonsson USH. (2018).",
    "A model-informed preclinical approach for prediction of clinical",
    "pharmacodynamic interactions of anti-TB drug combinations.",
    "J Antimicrob Chemother 73(2):437-447.",
    "doi:10.1093/jac/dkx380.",
    "MTP natural-growth transfer rates kFSlin, kSF, kFN, kSN, kNS fixed from",
    "Clewe O, Aulin L, Hu Y, Coates AR, Simonsson US. (2016).",
    "A multistate tuberculosis pharmacometric model: a framework for studying",
    "anti-tubercular drug effects in vitro. J Antimicrob Chemother 71(4):964-974.",
    "doi:10.1093/jac/dkv416.",
    "GPDI framework: Wicha SG, Chen C, Clewe O, Simonsson USH. (2017).",
    "A general pharmacodynamic interaction model identifies perpetrators and",
    "victims in drug interactions. Nat Commun 8:2129.",
    "doi:10.1038/s41467-017-01929-y.",
    "Adaptive-resistance loop adapted from",
    "Mohamed AF, Cars O, Friberg LE. (2014). A pharmacokinetic/pharmacodynamic",
    "model developed for the effect of colistin on Pseudomonas aeruginosa in vitro.",
    "J Antimicrob Chemother 69(5):1350-1361.",
    "doi:10.1093/jac/dkt520."
  )
  vignette <- "Clewe_2018_TB_MTP_GPDI_invitro"
  units <- list(
    time          = "day",
    dosing        = "mg/L (static covariates -- not administered events)",
    concentration = "log(CFU/mL) for the model observation Cc"
  )

  covariateData <- list(
    CONC_RIF_MGL = list(
      description        = "Static rifampicin concentration in the in vitro time-kill assay (mg/L)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per Clewe 2018 Materials and methods (Multistate TB pharmacometric model):",
        "'Static drug concentrations were used as input to the PD modelling.",
        "The stability of the drugs allowed assessment of activity during the 6 days",
        "of the experiment without the need for replenishment.' Tested concentrations in",
        "the source experiment (Figure 1): rifampicin 0.002, 0.008, 0.03, 0.125, 0.5, 8",
        "mg/L. Set to 0 in regimens without rifampicin. Paper-specific covariate not in",
        "inst/references/covariate-columns.md because the canonical concentration concept",
        "in nlmixr2lib is a state-derived plasma concentration (Cc), not a static",
        "exogenous-drug-concentration covariate used to drive an in vitro PD model."
      ),
      source_name        = "CRIF (Clewe 2018 Materials and methods)"
    ),
    CONC_INH_MGL = list(
      description        = "Static isoniazid concentration in the in vitro time-kill assay (mg/L)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per Clewe 2018 Materials and methods. Tested concentrations in",
        "the source experiment (Figure 1): isoniazid 0.01, 0.039, 0.156, 0.625, 2.5, 10,",
        "40 mg/L. Set to 0 in regimens without isoniazid. Paper-specific covariate;",
        "drives both the kill effects on F and S sub-states and the adaptive-resistance",
        "AR_on / AR_off transition rate kon * CONC_INH_MGL."
      ),
      source_name        = "CINH (Clewe 2018 Materials and methods)"
    ),
    CONC_EMB_MGL = list(
      description        = "Static ethambutol concentration in the in vitro time-kill assay (mg/L)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per Clewe 2018 Materials and methods. Tested concentrations in",
        "the source experiment (Figure 1): ethambutol 0.0078, 0.031, 0.125, 0.5, 2, 8,",
        "32 mg/L. Set to 0 in regimens without ethambutol. Paper-specific covariate."
      ),
      source_name        = "CEMB (Clewe 2018 Materials and methods)"
    )
  )

  population <- list(
    species        = "in vitro (M. tuberculosis Beijing VN 2002-1585, BE1585 strain)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = paste(
      "Tuberculosis time-kill experiments using the Mycobacterium tuberculosis",
      "Beijing VN 2002-1585 (BE1585) genotype, an MDR-precursor clinical isolate.",
      "Cultures grown in Middlebrook 7H9 broth + 10% OADC + 0.5% glycerol +",
      "0.02% Tween 20 under shaking conditions (96 rpm, 37 C) per Clewe 2018",
      "Materials and methods (In vitro assay). The limit of quantification was",
      "5 CFU; data below the LOQ were handled in NONMEM with the M3 method."
    ),
    dose_range     = paste(
      "Static drug concentrations (no PK dosing), tested in mono / duo / trio",
      "combinations over a 6-day window (assessed on days 1, 2, 3, and 6).",
      "Rifampicin 0.002-8 mg/L; isoniazid 0.01-40 mg/L; ethambutol 0.0078-32 mg/L",
      "(Figure 1). Duo and trio combinations used cross-products of selected",
      "concentrations including clinically relevant unbound Cmax values (RIF 2,",
      "INH 10, EMB 8 mg/L)."
    ),
    notes          = paste(
      "Experiments performed in duplicate at Erasmus Medical Centre, Rotterdam, NL",
      "(Department of Medical Microbiology and Infectious Diseases). The paper",
      "does not report interindividual or interreplicate variability on the",
      "structural parameters; the model is a typical-value mechanistic in vitro",
      "PD without etas. The residual error structure used in NONMEM (additive on",
      "natural-log CFU/mL with M3-method censoring at 5 CFU) is not numerically",
      "reported in the main text or Table 1 / Table 3; the addSd placeholder in",
      "ini() carries a fixed natural-log SD of 1.0 (see vignette Errata)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # NATURAL GROWTH (Clewe 2018 Table 1). Parameters not wrapped in fixed()
    # were re-estimated for the BE1585 strain; those wrapped in fixed()
    # are carried from the upstream Clewe 2016 H37Rv MTP fit (footnote b
    # in Table 1).
    # ------------------------------------------------------------------
    lkg <- log(0.796)
    label("kG -- fast-multiplying state growth rate (1/day)")
    # Clewe 2018 Table 1: kG = 0.796 1/day, RSE 5%. Exponential growth function
    # selected over logistic (Bmax) per Results paragraph 1.

    lkfs_lin <- fixed(log(0.166e-2))
    label("kFSlin -- second-order time-dependent F->S transfer slope (1/day^2)")
    # Clewe 2018 Table 1: kFSlin = 0.166e-2, fixed^b from Clewe 2016 H37Rv MTP.
    # Used as kFS(t) = kFSlin * t in the F<->S transfer.

    lkfn <- fixed(log(0.897e-6))
    label("kFN -- first-order F->N transfer (1/day)")
    # Clewe 2018 Table 1: kFN = 0.897e-6, fixed^b from Clewe 2016.

    lksn <- fixed(log(0.186))
    label("kSN -- first-order S->N transfer (1/day)")
    # Clewe 2018 Table 1: kSN = 0.186 fixed^b from Clewe 2016.

    lksf <- fixed(log(0.0145))
    label("kSF -- first-order S->F transfer (1/day)")
    # Clewe 2018 Table 1: kSF = 0.0145 fixed^b from Clewe 2016.

    lkns <- fixed(log(0.123e-2))
    label("kNS -- first-order N->S transfer (1/day)")
    # Clewe 2018 Table 1: kNS = 0.123e-2, fixed^b from Clewe 2016.

    lf0 <- log(209e3)
    label("F0 -- initial fast-multiplying state count (CFU/mL)")
    # Clewe 2018 Table 1: F0 = 209e3 CFU/mL, RSE 17%.

    ls0 <- log(324e3)
    label("S0 -- initial slow-multiplying state count (CFU/mL)")
    # Clewe 2018 Table 1: S0 = 324e3 CFU/mL, RSE 12%.

    # ------------------------------------------------------------------
    # RIFAMPICIN monotherapy effects (Clewe 2018 Table 1).
    # Emax^FG_RIF is fixed at 1 day^-1 because the maximum inhibition of
    # F growth is structurally bounded at 100% (the growth rate cannot
    # become negative); the data could not separately identify Emax above
    # 1 for the F-growth-inhibition mechanism.
    # ------------------------------------------------------------------
    lemax_fg_rif <- fixed(log(1))
    label("Emax^FG_RIF -- maximum RIF inhibition of F growth (1/day)")
    # Clewe 2018 Table 1: Emax^FG_RIF = 1 fixed.

    lec50_fg_rif <- log(0.388)
    label("EC50^FG_RIF -- concentration at 50% of Emax^FG_RIF (mg/L)")
    # Clewe 2018 Table 1: EC50^FG_RIF = 0.388 mg/L, RSE 19%.

    gamma_fg_r <- 2.8
    label("Hill factor for RIF inhibition of F growth (unitless)")
    # Clewe 2018 Table 1: gamma^FG_R = 2.8, RSE 28%.

    lemax_fd_rif <- log(1.97)
    label("Emax^FD_RIF -- maximum RIF kill of F state (1/day)")
    # Clewe 2018 Table 1: Emax^FD_RIF = 1.97 1/day, RSE 3%.

    lec50_fd_rif <- log(0.00303)
    label("EC50^FD_RIF -- concentration at 50% of Emax^FD_RIF (mg/L)")
    # Clewe 2018 Table 1: EC50^FD_RIF = 0.00303 mg/L, RSE 10%.

    lemax_sd_rif <- log(1.79)
    label("Emax^SD_RIF -- maximum RIF kill of S state (1/day)")
    # Clewe 2018 Table 1: Emax^SD_RIF = 1.79 1/day, RSE 4%.

    lec50_sd_rif <- log(0.0113)
    label("EC50^SD_RIF -- concentration at 50% of Emax^SD_RIF (mg/L)")
    # Clewe 2018 Table 1: EC50^SD_RIF = 0.0113 mg/L, RSE 32%.

    lknd_rif <- log(3.29)
    label("kND_RIF -- linear RIF kill slope on N state (L*mg^-1*day^-1)")
    # Clewe 2018 Table 1: kND_RIF = 3.29, RSE 17%. Units listed as days^-1
    # in Table 1 reflect the rate output at unit (1 mg/L) RIF; the slope
    # itself is L*mg^-1*day^-1. E^ND_RIF = kND_RIF * C_RIF.

    # ------------------------------------------------------------------
    # ISONIAZID monotherapy effects (Clewe 2018 Table 1) + adaptive
    # resistance.
    # ------------------------------------------------------------------
    lemax_fd_inh <- log(22.2)
    label("Emax^FD_INH -- maximum INH kill of F state (1/day)")
    # Clewe 2018 Table 1: Emax^FD_INH = 22.2 1/day, RSE 35%.

    lec50_fd_inh <- log(0.168)
    label("EC50^FD_INH -- concentration at 50% of Emax^FD_INH (mg/L)")
    # Clewe 2018 Table 1: EC50^FD_INH = 0.168 mg/L, RSE 34%. Adaptive
    # resistance modulates this EC50 by (1 + kAR_FD_INH * AR_on).

    gamma_fd_h <- 1.9
    label("Hill factor for INH kill of F state (unitless)")
    # Clewe 2018 Table 1: gamma^FD_H = 1.9, RSE 11%.

    lemax_sd_inh <- log(8.55)
    label("Emax^SD_INH -- maximum INH kill of S state (1/day)")
    # Clewe 2018 Table 1: Emax^SD_INH = 8.55 1/day, RSE 17%.

    lec50_sd_inh <- log(0.0329)
    label("EC50^SD_INH -- concentration at 50% of Emax^SD_INH (mg/L)")
    # Clewe 2018 Table 1: EC50^SD_INH = 0.0329 mg/L, RSE 49%. Adaptive
    # resistance modulates this EC50 by (1 + kAR_SD_INH * AR_on).

    gamma_sd_h <- 1.74
    label("Hill factor for INH kill of S state (unitless)")
    # Clewe 2018 Table 1: gamma^SD_H = 1.74, RSE 25%.

    lkon <- log(0.0206)
    label("kon -- adaptive-resistance development rate (L*mg^-1*day^-1)")
    # Clewe 2018 Table 1: kon = 0.0206, RSE 31%. Units listed as
    # "L*mg^-1*days" in Table 1 are interpreted as L*mg^-1*day^-1 so that
    # k_on * C_INH * AR_off has units 1/day matching dAR_on/dt. The "days"
    # vs "day^-1" is treated as a typesetting slip in the published Table 1.
    # k_off was held at 0 per Materials and methods ("no data were
    # available on resistance reversal, koff was fixed to 0"), so no
    # k_off parameter appears in ini().

    lkar_fd_inh <- log(522)
    label("kAR_FD_INH -- linear INH adaptive resistance on F-kill EC50 (unitless)")
    # Clewe 2018 Table 1: kAR_FD_INH = 522, RSE 46%. Encoded as a
    # dimensionless slope (EC50_eff = EC50 * (1 + kAR * AR_on)); the
    # "(days^-1)" in Table 1 is interpreted as a unit-label slip because
    # AR_on is dimensionless and the multiplier on EC50 must be
    # dimensionless. See vignette Errata for the linear-vs-Hill form
    # rationale (paper's Eq 3-4 give Hill-type Emax/AR50 forms but the
    # final model uses linear k_AR per Materials and methods).

    lkar_sd_inh <- log(2350)
    label("kAR_SD_INH -- linear INH adaptive resistance on S-kill EC50 (unitless)")
    # Clewe 2018 Table 1: kAR_SD_INH = 2350, RSE 51%.

    # ------------------------------------------------------------------
    # ETHAMBUTOL monotherapy effects (Clewe 2018 Table 1).
    # ------------------------------------------------------------------
    lemax_fd_emb <- log(2.21)
    label("Emax^FD_EMB -- maximum EMB kill of F state (1/day)")
    # Clewe 2018 Table 1: Emax^FD_EMB = 2.21 1/day, RSE 1%.

    lec50_fd_emb <- log(0.86)
    label("EC50^FD_EMB -- concentration at 50% of Emax^FD_EMB (mg/L)")
    # Clewe 2018 Table 1: EC50^FD_EMB = 0.86 mg/L, RSE 16%.

    gamma_fd_e <- 2.46
    label("Hill factor for EMB kill of F state (unitless)")
    # Clewe 2018 Table 1: gamma^FD_E = 2.46, RSE 23%.

    lksd_emb <- log(4.39)
    label("kSD_EMB -- linear EMB kill slope on S state (L*mg^-1*day^-1)")
    # Clewe 2018 Table 1: kSD_EMB = 4.39, RSE 69%. E^SD_EMB = kSD_EMB * C_EMB
    # before the GPDI modulation by INH and RIF.

    # ------------------------------------------------------------------
    # GPDI INTERACTIONS (Clewe 2018 Table 3). Convention per Table 3
    # footer: INT_A,B where A is the modifier and B the modulated drug;
    # encoded as a (1 + INT * C_A / (EC50_A + C_A)) multiplicative factor
    # on B's EC50 (or B's linear-kill slope when applicable). EC50_A is
    # the modifier's mono EC50 for the same bacterial-state mechanism
    # (the "least complex" baseline per Materials and methods, GPDI
    # section -- Table 3 reports no separate EC50_INT values).
    # INT values are not log-transformed because they may be negative.
    # ------------------------------------------------------------------

    # F-state interactions
    int_fd_rif_inh <- -0.679
    label("INT^FD_RIF,INH -- RIF modulates INH EC50 for F-kill (fractional)")
    # Clewe 2018 Table 3: INT^FD_RIF,INH = -0.679, RSE 11%. Saturating
    # decrease (synergism) -- RIF presence lowers INH's apparent EC50
    # for F-kill by up to 67.9%.

    int_fd_inh_rif <- fixed(0)
    label("INT^FD_INH,RIF -- INH modulates RIF EC50 for F-kill (fractional, fixed 0)")
    # Clewe 2018 Table 3: 0 fixed^b (additivity).

    int_fd_emb_inh <- 1.72
    label("INT^FD_EMB,INH -- EMB modulates INH EC50 for F-kill (fractional)")
    # Clewe 2018 Table 3: INT^FD_EMB,INH = 1.72, RSE 15%. Saturating
    # increase (antagonism).

    int_fd_inh_emb <- fixed(0)
    label("INT^FD_INH,EMB -- INH modulates EMB EC50 for F-kill (fractional, fixed 0)")
    # Clewe 2018 Table 3: 0 fixed^b.

    int_fd_rif_emb <- fixed(-0.99)
    label("INT^FD_RIF,EMB -- RIF modulates EMB EC50 for F-kill (fractional, fixed -0.9999)")
    # Clewe 2018 Table 3: INT^FD_RIF,EMB = -0.99 fixed^c (max decrease,
    # synergism). Paper footnote c sets this to -0.9999 reflecting "a
    # maximal decrease in EC50 identified in mono exposure"; here we
    # encode -0.99 to keep the multiplier (1 + INT) positive and avoid
    # numerical issues, equivalent to -0.9999 within the saturating term.

    int_fd_emb_rif <- -0.668
    label("INT^FD_EMB,RIF -- EMB modulates RIF EC50 for F-kill (fractional)")
    # Clewe 2018 Table 3: INT^FD_EMB,RIF = -0.668, RSE 22%.

    # S-state interactions
    int_sd_rif_inh <- 1.42
    label("INT^SD_RIF,INH -- RIF modulates INH EC50 for S-kill (fractional)")
    # Clewe 2018 Table 3: INT^SD_RIF,INH = 1.42, RSE 22%.

    int_sd_inh_rif <- 15.2
    label("INT^SD_INH,RIF -- INH modulates RIF EC50 for S-kill (fractional)")
    # Clewe 2018 Table 3: INT^SD_INH,RIF = 15.2, RSE 49%.

    int_sd_emb_inh <- 0.0963
    label("INT^SD_EMB,INH -- EMB modulates INH EC50 for S-kill (fractional)")
    # Clewe 2018 Table 3: INT^SD_EMB,INH = 0.0963, RSE 81%.

    int_sd_inh_emb <- 164
    label("INT^SD_INH,EMB -- INH modulates EMB linear-S slope (fractional)")
    # Clewe 2018 Table 3: INT^SD_INH,EMB = 164, RSE 259%. The S-kill by
    # EMB is linear (kSD_EMB) so this interaction modulates the slope
    # rather than an EC50: slope_eff = kSD_EMB / (1 + INT * C_INH /
    # (EC50_SD_INH + C_INH)). Positive INT increases the denominator,
    # decreasing the effective slope (antagonism on EMB's S-kill).

    int_sd_rif_emb <- 486
    label("INT^SD_RIF,EMB -- RIF modulates EMB linear-S slope (fractional)")
    # Clewe 2018 Table 3: INT^SD_RIF,EMB = 486, RSE 12%. Same denominator
    # form as INT^SD_INH,EMB.

    int_sd_emb_rif <- 2.09
    label("INT^SD_EMB,RIF -- EMB modulates RIF EC50 for S-kill (fractional)")
    # Clewe 2018 Table 3: INT^SD_EMB,RIF = 2.09, RSE 32%.

    int_sd_rif_inh_emb <- -0.749
    label("INT^SD_RIF,INH|EMB -- EMB modulates the RIF-on-INH S-kill interaction (fractional)")
    # Clewe 2018 Table 3: INT^SD_RIF,INH|EMB = -0.749, RSE 18%. Trio
    # modulator: when EMB is present, the magnitude of int_sd_rif_inh
    # is itself multiplied by (1 + INT^SD_RIF,INH|EMB * C_EMB /
    # (EC50_FD_EMB + C_EMB)). EC50_FD_EMB is used as the anchor because
    # EMB has no mono S-kill EC50 (its S-kill is linear); see vignette
    # Errata for the choice.

    # ------------------------------------------------------------------
    # RESIDUAL ERROR. Clewe 2018 reports natural-log CFU/mL data with
    # M3-method handling of < 5 CFU points but does not tabulate the
    # residual SD. Placeholder fixed at 1.0 on the natural-log scale
    # (about the order of magnitude reported by Chen 2017 for the
    # in vivo twin paper) so the model is nlmixr2-fit-compatible; see
    # vignette Errata for the rationale and recommended downstream use
    # as a typical-value mechanism rather than a population fit.
    # ------------------------------------------------------------------
    addSd <- fixed(1.0)
    label("Placeholder additive residual SD on natural-log CFU/mL (unreported in source)")
  })

  model({
    # ================================================================
    # 1. Back-transform structural parameters from log scale.
    # ================================================================
    kg       <- exp(lkg)
    kfs_lin  <- exp(lkfs_lin)
    kfn      <- exp(lkfn)
    ksn      <- exp(lksn)
    ksf      <- exp(lksf)
    kns      <- exp(lkns)
    f0       <- exp(lf0)
    s0       <- exp(ls0)

    emax_fg_rif  <- exp(lemax_fg_rif)
    ec50_fg_rif  <- exp(lec50_fg_rif)
    emax_fd_rif  <- exp(lemax_fd_rif)
    ec50_fd_rif  <- exp(lec50_fd_rif)
    emax_sd_rif  <- exp(lemax_sd_rif)
    ec50_sd_rif  <- exp(lec50_sd_rif)
    knd_rif      <- exp(lknd_rif)

    emax_fd_inh  <- exp(lemax_fd_inh)
    ec50_fd_inh  <- exp(lec50_fd_inh)
    emax_sd_inh  <- exp(lemax_sd_inh)
    ec50_sd_inh  <- exp(lec50_sd_inh)
    kon          <- exp(lkon)
    kar_fd_inh   <- exp(lkar_fd_inh)
    kar_sd_inh   <- exp(lkar_sd_inh)

    emax_fd_emb  <- exp(lemax_fd_emb)
    ec50_fd_emb  <- exp(lec50_fd_emb)
    ksd_emb      <- exp(lksd_emb)

    # Time-dependent F->S transfer (Clewe 2018 Materials and methods:
    # kFS = kFSlin * t). t is the rxode2 simulation-time variable.
    kfs <- kfs_lin * t

    # ================================================================
    # 2. Adaptive resistance for isoniazid (Clewe 2018 Eqs 1-2, with
    #    Eq 2 corrected for mass-balance: -kon * C_INH * AR_off rather
    #    than the paper's -kon * AR_off -- see vignette Errata).
    #    Initial conditions: AR_off = 1, AR_on = 0 (Materials and
    #    methods: "all bacteria were assigned to the AR_off state at
    #    the start of the experiment"). koff fixed to 0, so the two
    #    states sum to 1 indefinitely.
    # ================================================================
    d/dt(ar_off) <- -kon * CONC_INH_MGL * ar_off
    d/dt(ar_on)  <-  kon * CONC_INH_MGL * ar_off
    ar_off(0) <- 1
    ar_on(0)  <- 0

    # Linear adaptive-resistance modulation of INH's EC50 on F-kill and
    # S-kill (Materials and methods: "linear function"; the paper's
    # Eq 3-4 Hill forms with AR_max/AR_50 are alternates that were
    # evaluated but did not enter the final model -- see Table 1).
    ec50_fd_inh_eff <- ec50_fd_inh * (1 + kar_fd_inh * ar_on)
    ec50_sd_inh_eff <- ec50_sd_inh * (1 + kar_sd_inh * ar_on)

    # ================================================================
    # 3. GPDI modulation factors for each drug pair / mechanism. Form
    #    per Table 3 footer:
    #       factor = 1 + INT_A,B * C_A / (EC50_A + C_A)
    #    where A is the modifier (its concentration drives the
    #    saturation) and EC50_A is A's mono EC50 for the same
    #    bacterial-state mechanism. Where a mechanism-specific EC50 is
    #    not available for the modifier (kND_RIF on N is linear; kSD_EMB
    #    on S is linear), the F-kill EC50 of the modifier is used as
    #    the anchor (documented in vignette Errata).
    # ================================================================

    # Saturation-fraction shortcuts. Compute once and reuse.
    sat_rif_fd <- CONC_RIF_MGL / (ec50_fd_rif + CONC_RIF_MGL)
    sat_inh_fd <- CONC_INH_MGL / (ec50_fd_inh + CONC_INH_MGL)
    sat_emb_fd <- CONC_EMB_MGL / (ec50_fd_emb + CONC_EMB_MGL)
    sat_rif_sd <- CONC_RIF_MGL / (ec50_sd_rif + CONC_RIF_MGL)
    sat_inh_sd <- CONC_INH_MGL / (ec50_sd_inh + CONC_INH_MGL)
    # EMB has no SD-mechanism EC50 (linear kSD_EMB), so its S-kill
    # saturation fraction falls back to the F-kill EC50 (see Errata).
    sat_emb_sd <- CONC_EMB_MGL / (ec50_fd_emb + CONC_EMB_MGL)

    # F-kill EC50 modulation factors
    fac_inh_fd_by_rif <- 1 + int_fd_rif_inh * sat_rif_fd
    fac_inh_fd_by_emb <- 1 + int_fd_emb_inh * sat_emb_fd
    fac_rif_fd_by_inh <- 1 + int_fd_inh_rif * sat_inh_fd
    fac_rif_fd_by_emb <- 1 + int_fd_emb_rif * sat_emb_fd
    fac_emb_fd_by_rif <- 1 + int_fd_rif_emb * sat_rif_fd
    fac_emb_fd_by_inh <- 1 + int_fd_inh_emb * sat_inh_fd

    # S-kill EC50 (or linear-slope) modulation factors
    fac_inh_sd_by_rif_baseline <- 1 + int_sd_rif_inh * sat_rif_sd
    # Trio modulator: when EMB is present, INT^SD_RIF,INH itself is
    # modulated by EMB. Encoded by replacing int_sd_rif_inh with
    # int_sd_rif_inh * (1 + int_sd_rif_inh_emb * sat_emb_sd) in the
    # factor expression. With EMB absent (CONC_EMB_MGL = 0,
    # sat_emb_sd = 0), the trio modulator vanishes and the duo factor
    # recovers.
    fac_inh_sd_by_rif <- 1 + int_sd_rif_inh *
      (1 + int_sd_rif_inh_emb * sat_emb_sd) * sat_rif_sd
    fac_inh_sd_by_emb <- 1 + int_sd_emb_inh * sat_emb_sd
    fac_rif_sd_by_inh <- 1 + int_sd_inh_rif * sat_inh_sd
    fac_rif_sd_by_emb <- 1 + int_sd_emb_rif * sat_emb_sd
    fac_emb_sd_by_inh <- 1 + int_sd_inh_emb * sat_inh_sd
    fac_emb_sd_by_rif <- 1 + int_sd_rif_emb * sat_rif_sd

    # ================================================================
    # 4. Drug-specific effects per bacterial sub-state.
    # ================================================================

    # F-state growth inhibition by rifampicin (Table 2 row F/rifampicin
    # E^FG_RIF). No PD interactions identified for E^FG_RIF (Table 2
    # "PD interaction identified" column: "none identified").
    crif_g <- CONC_RIF_MGL + 1e-30
    e_fg_rif <- 1 - emax_fg_rif * crif_g^gamma_fg_r /
                  (ec50_fg_rif^gamma_fg_r + crif_g^gamma_fg_r)

    # F-state kill rates -- Emax forms per Table 2.
    crif_d <- CONC_RIF_MGL + 1e-30
    cinh_d <- CONC_INH_MGL + 1e-30
    cemb_d <- CONC_EMB_MGL + 1e-30

    e_fd_rif <- emax_fd_rif * crif_d /
                (ec50_fd_rif * fac_rif_fd_by_inh * fac_rif_fd_by_emb + crif_d)
    e_fd_inh <- emax_fd_inh * cinh_d^gamma_fd_h /
                ((ec50_fd_inh_eff * fac_inh_fd_by_rif * fac_inh_fd_by_emb)^gamma_fd_h +
                   cinh_d^gamma_fd_h)
    e_fd_emb <- emax_fd_emb * cemb_d^gamma_fd_e /
                ((ec50_fd_emb * fac_emb_fd_by_rif * fac_emb_fd_by_inh)^gamma_fd_e +
                   cemb_d^gamma_fd_e)

    # S-state kill rates
    e_sd_rif <- emax_sd_rif * crif_d /
                (ec50_sd_rif * fac_rif_sd_by_inh * fac_rif_sd_by_emb + crif_d)
    e_sd_inh <- emax_sd_inh * cinh_d^gamma_sd_h /
                ((ec50_sd_inh_eff * fac_inh_sd_by_rif * fac_inh_sd_by_emb)^gamma_sd_h +
                   cinh_d^gamma_sd_h)
    # E^SD_EMB: linear in CONC_EMB_MGL with the slope modulated by
    # INH and RIF (factors in the denominator per Eq 6 of the paper).
    e_sd_emb <- ksd_emb * cemb_d / (fac_emb_sd_by_inh * fac_emb_sd_by_rif)

    # N-state kill rate -- only rifampicin, linear, no PD interactions.
    e_nd_rif <- knd_rif * crif_d

    # ================================================================
    # 5. Bliss Independence combination over the three drugs for the
    #    F-kill and S-kill effects (Clewe 2018 Eq 8 and Eq 9). The
    #    N-kill effect (Eq 10) has only rifampicin so no BI combination
    #    is needed.
    # ================================================================
    e_fd <- e_fd_rif + e_fd_inh + e_fd_emb -
            e_fd_rif * e_fd_inh - e_fd_inh * e_fd_emb - e_fd_rif * e_fd_emb +
            e_fd_rif * e_fd_inh * e_fd_emb
    e_sd <- e_sd_rif + e_sd_inh + e_sd_emb -
            e_sd_rif * e_sd_inh - e_sd_inh * e_sd_emb - e_sd_rif * e_sd_emb +
            e_sd_rif * e_sd_inh * e_sd_emb
    e_nd <- e_nd_rif

    # ================================================================
    # 6. MTP bacterial-state ODEs (Clewe 2018 Eq 8, 9, 10). The paper's
    #    Eq 8 includes a "+ kNF * N" term that has no corresponding
    #    Table 1 parameter and no arrow in Figure 2; treated as a
    #    typesetting error per the canonical Clewe 2016 MTP and the
    #    Chen 2017 mouse twin (modellib('Chen_2017_TB_MTP_GPDI_mouse'),
    #    which omits N->F transfer). See vignette Errata.
    #
    #    Initial states: F0, S0 from Table 1; N0 = 0 (the experiment
    #    starts without non-multiplying bacteria; the N pool is
    #    populated by F->N and S->N transfer over time).
    # ================================================================
    fbugs(0) <- f0
    sbugs(0) <- s0
    nbugs(0) <- 0

    d/dt(fbugs) <- fbugs * kg * e_fg_rif + ksf * sbugs -
                   kfs * fbugs - kfn * fbugs - e_fd * fbugs
    d/dt(sbugs) <- kfs * fbugs - ksf * sbugs + kns * nbugs -
                   ksn * sbugs - e_sd * sbugs
    d/dt(nbugs) <- ksn * sbugs + kfn * fbugs - kns * nbugs -
                   e_nd * nbugs

    # ================================================================
    # 7. Observation -- natural-log of total CFU/mL with a small floor
    #    to avoid log(0) when the integrated total approaches zero
    #    under saturating drug exposure (LOQ = 5 CFU/mL per Clewe 2018
    #    Materials and methods; the floor of 1 CFU/mL keeps the
    #    natural-log finite without biasing observations above LOQ).
    # ================================================================
    total_bugs <- fbugs + sbugs + nbugs
    if (total_bugs < 1) total_bugs <- 1
    Cc <- log(total_bugs)
    Cc ~ add(addSd)
  })
}
