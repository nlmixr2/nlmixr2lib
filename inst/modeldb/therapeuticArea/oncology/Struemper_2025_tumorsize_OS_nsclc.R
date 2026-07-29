Struemper_2025_tumorsize_OS_nsclc <- function() {
  description <- paste(
    "Joint tumor-size (TS) / overall-survival (OS) framework model for",
    "non-small cell lung cancer (NSCLC), developed by Struemper et al.",
    "(GSK) on pooled individual-level data from 786 participants across",
    "seven GSK-sponsored clinical trials (INDUCE-1, INDUCE-2, Entree",
    "Lung Part 2, GARNET, AMBER, INTR@PID LUNG 037, PERLA) spanning",
    "immunotherapy, chemotherapy, and combinations thereof. The TS",
    "sub-model is the bi-exponential Stein model with per-treatment-arm",
    "typical tumor-growth (kge) and tumor-shrinkage (kse) rates",
    "selected at simulation / fit time via the TRT categorical",
    "covariate (12 levels; see covariateData[[TRT]] for the integer",
    "coding). The OS sub-model is an accelerated failure time (AFT)",
    "log-normal survival with a treatment-agnostic link to individual",
    "TS parameters: tumor growth rate kge enters via a fixed Emax",
    "function, baseline TS and time-to-tumor-growth (TTG) enter",
    "linearly, and three baseline laboratory covariates (albumin, total",
    "protein, neutrophil-to-lymphocyte ratio) are additive on the",
    "log-time scale. TS is observed in mm (sum of longest diameters of",
    "target lesions per RECIST 1.1); survival is reported on the day",
    "scale in the source paper and converted to weeks inside model() so",
    "the entire model runs on a single weeks time axis."
  )
  reference <- paste(
    "Struemper H, Rathi C, Muliaditan M, Goulooze SC, Franzese RC,",
    "Mantero A, Melhem M, Post TM, Visser SAG.",
    "Development of a Joint Tumor Size-Overall Survival Modeling and",
    "Simulation Framework Supporting Oncology Development",
    "Decision-Making.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(6):1006-1017.",
    "doi:10.1002/psp4.70002. PMID: 39985158.",
    sep = " "
  )
  vignette <- "Struemper_2025_tumorsize_OS_nsclc"

  units <- list(
    time          = "weeks (TS dynamics and OS hazard; OS mu_OS published in log-days is converted to log-weeks inside model() so a single weeks time axis carries both sub-models)",
    dosing        = "n/a (no PK input; treatment effect encoded entirely via per-arm typical kge and kse selected by the TRT categorical covariate)",
    concentration = "mm (TS observable = sum of longest diameters of target lesions per RECIST 1.1; the canonical covariate column corresponding to the time-varying TS observable is TUM_SLD)"
  )

  covariateData <- list(
    TRT = list(
      description = "Per-subject treatment-arm integer indicator selecting the per-arm typical TS parameters (kge, kse) at simulation / fit time.",
      units       = "(categorical / integer-coded)",
      type        = "categorical",
      source_name = "TRT",
      notes       = paste(
        "Integer coding (Struemper 2025 Table 1 + Table 2):",
        "  1 = PEMBRO (pembrolizumab; INTR@PID LUNG 037; n = 152) -- canonical default, largest single-agent PD-1 cohort",
        "  2 = FELAD (feladilimab; INDUCE-1; n = 52)",
        "  3 = CHEMO (docetaxel; Entree Lung Part 2; n = 34)",
        "  4 = FELAD+CHEMO (feladilimab + docetaxel; pooled INDUCE-1 + Entree Lung Part 2; n = 78)",
        "  5 = IO-COMBO (feladilimab + IO [cobolimab / tremelimumab / dostarlimab / GSK3174998]; pooled INDUCE-1 + INDUCE-2; n = 23)",
        "  6 = DOSTAR (dostarlimab; GARNET; n = 67)",
        "  7 = DOSTAR+CHEMO (PERLA; n = 121)",
        "  8 = PEMBRO+CHEMO (PERLA; n = 122)",
        "  9 = COBO100MG+DOSTAR cohort B (AMBER; n = 14)",
        " 10 = COBO300MG+DOSTAR cohort B (AMBER; n = 41)",
        " 11 = COBO900MG+DOSTAR cohort B (AMBER; n = 29)",
        " 12 = COBO300MG+DOSTAR cohort D (AMBER; n = 53)",
        "Arms 1 and 6-12 contain pembrolizumab or dostarlimab and therefore receive the PDL1_TUM effect on kse (gated by the derived has_pd1 indicator inside model()); arms 2-5 do not.",
        sep = "\n"
      )
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase activity (paper alias LDHBL).",
      units       = "IU/L (equivalently U/L; the canonical register treats IU/L and U/L interchangeably)",
      type        = "continuous",
      source_name = "LDHBL",
      notes       = "Power scaling on kge: kge = exp(lkge_arm + eta) * (LDH / 225.5) ^ e_ldh_kge with reference 225.5 IU/L (Struemper 2025 Figure 3 caption typical-subject value)."
    ),
    PDL1_TUM = list(
      description = "Baseline tumor PD-L1 expression by immunohistochemistry, Tumor Proportion Score (TPS): percent of viable tumor cells with complete or partial membrane staining at any intensity.",
      units       = "percent (0-100)",
      type        = "continuous",
      source_name = "PD-L1",
      notes       = "Exponential effect on kse: kse = exp(lkse_arm + eta) * exp(PDL1_TUM * e_pdl1_tum_ks * has_pd1) where has_pd1 = 1 if TRT is 1 or in {6, 7, 8, 9, 10, 11, 12} (pembrolizumab- or dostarlimab-containing arms) and 0 otherwise. PD-L1 assay was 22C3 in all studies except INTR@PID (TRT = 1) which used 73-10 and had no individual values available -- the source paper imputed 70 percent for all INTR@PID subjects, the population median among high (>=50 percent) expressers; users replicating the paper's PEMBRO arm should set PDL1_TUM = 70 for all TRT = 1 subjects."
    ),
    NTARGET_GE3 = list(
      description = "Binary indicator: 1 if baseline number of target lesions (per RECIST 1.1) is three or more, 0 if one or two.",
      units       = "(binary)",
      type        = "binary",
      source_name = "NTARGET",
      notes       = "Multiplicative effect on the typical-value baseline TS: TSb = exp(lrbase + eta) * (1 + e_ntarget_ge3_tsb * NTARGET_GE3). The source paper reports the linear-continuous form TVTSb * (1 + (NTARGET - 3) * 0.288) centred at NTARGET = 3 (Table 3 footnote a). The library binarises at the paper's reference value of 3 per the count-covariate policy (cf. MET_GE4 in Bruno 2005); the per-lesion linear coefficient 0.288 is reused as the single-step binary coefficient. Deviation documented in vignette Errata. Reference category 0 = one or two target lesions. Derive from a raw count column via NTARGET_GE3 = as.integer(NTARGET >= 3)."
    ),
    ALB = list(
      description = "Baseline serum albumin (paper alias ALBBL).",
      units       = "g/L",
      type        = "continuous",
      source_name = "ALBBL",
      notes       = "Linear centred deviation on mu_OS (log-days): CV4 = (ALB - 39.4) * 0.0452. Reference 39.4 g/L = population median (Struemper 2025 Figure 3 caption)."
    ),
    TPRO = list(
      description = "Baseline serum total protein (paper alias TPROBL).",
      units       = "g/L",
      type        = "continuous",
      source_name = "TPROBL",
      notes       = "Linear centred deviation on mu_OS (log-days): CV5 = (TPRO - 71) * 0.0194. Reference 71 g/L = population median (Struemper 2025 Figure 3 caption)."
    ),
    NLR = list(
      description = "Baseline neutrophil-to-lymphocyte ratio (paper alias NLRBL); capped at 100 in the source analysis per Methods Section 2.4 to avoid stability issues from four outlying values (>3000).",
      units       = "ratio (unitless)",
      type        = "continuous",
      source_name = "NLRBL",
      notes       = "Linear centred deviation on mu_OS (log-days): CV6 = (NLR - 4.1) * (-0.0141). Reference 4.1 = population median (Struemper 2025 Figure 3 caption). Cap NLR at 100 on data ingestion to match the source paper."
    )
  )

  population <- list(
    species         = "human (adults with advanced/metastatic NSCLC)",
    n_subjects      = 786L,
    n_studies       = 7L,
    age_range       = "not reported in this paper; advanced/metastatic NSCLC trial cohorts",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "not reported in this paper at the pooled level (per-study demographics are in the underlying trial publications)",
    disease_state   = "advanced/metastatic NSCLC (locally advanced, Stage IIIb/IIIc, Stage IV, recurrent, or metastatic per each study's inclusion criteria)",
    dose_range      = "n/a (no PK input; per-arm dosing was the protocol-defined dose per study)",
    regions         = "multiregional across seven clinical trials (per-study geographic mix not pooled in this paper)",
    notes           = paste(
      "Pooled-cohort baseline covariate medians (Struemper 2025 Figure 3 caption; population-typical reference values used in the model's centring constants):",
      "  Baseline lactate dehydrogenase (LDH) = 225.5 IU/L",
      "  Baseline tumor PD-L1 expression (PDL1_TUM) = 15 percent",
      "  Baseline number of target lesions = 2 (population median); model equation centres at NTARGET = 3 (linear-continuous form); library binarises at NTARGET_GE3",
      "  Baseline tumor size (TSb) = 74 mm (population median; model typical-value parameter TVTSb = 86.07 mm corresponds to NTARGET = 3 reference)",
      "  Tumor growth rate (kge) = 0.013 1/week (illustrative typical subject)",
      "  Time to tumor growth (TTG) = 14 weeks (illustrative typical subject)",
      "  Baseline albumin (ALB) = 39.4 g/L",
      "  Baseline total protein (TPRO) = 71 g/L",
      "  Baseline neutrophil-to-lymphocyte ratio (NLR) = 4.1",
      "",
      "Per-arm subject counts (Struemper 2025 Table 1) and TS-model rate constants (Struemper 2025 Table 2) are encoded inline as per-arm `lkge_<arm>` and `lkse_<arm>` parameters selected by the TRT categorical covariate; see covariateData[[TRT]]$notes for the integer-to-cohort mapping.",
      "",
      "PEMBRO (TRT = 1, n = 152) is the model's canonical default arm: the largest single-agent PD-1-inhibitor cohort. Users simulating other arms supply the corresponding TRT integer in the event table; no parameter override is required because the per-arm typical values are already in ini().",
      sep = "\n"
    )
  )

  ini({
    # ----- Baseline TS (TUM_SLD) typical value, Stein TVTSb at reference NTARGET = 3 -----
    lrbase <- log(86.07)             ; label("log baseline tumor size (TVTSb, mm) at reference NTARGET = 3 -- Stein lrbase")  # Struemper 2025 Table 2 row TVTSb = 86.07 (95% CI 82.35, 89.95)

    # ----- Per-arm typical tumor growth rate constants kge (1/week), log-transformed -----
    # All from Struemper 2025 Table 2 (per-arm TVKG estimates from the single joint NONMEM fit).
    lkge_pembro      <- log(0.00651)  ; label("log TVKG PEMBRO (TRT = 1; 1/week)")                                # Table 2: TVKG PEMBRO = 0.00651
    lkge_felad       <- log(0.01135)  ; label("log TVKG FELAD (TRT = 2; 1/week)")                                 # Table 2: TVKG FELAD = 0.01135
    lkge_chemo       <- log(0.01822)  ; label("log TVKG CHEMO (TRT = 3; 1/week)")                                 # Table 2: TVKG CHEMO = 0.01822
    lkge_feladchemo  <- log(0.01672)  ; label("log TVKG FELAD+CHEMO (TRT = 4; 1/week)")                           # Table 2: TVKG CHEMO+FELAD = 0.01672
    lkge_iocombo     <- log(0.01442)  ; label("log TVKG IO-COMBO (TRT = 5; 1/week)")                              # Table 2: TVKG IO-COMBO = 0.01442
    lkge_dostar      <- log(0.008105) ; label("log TVKG DOSTAR (TRT = 6; 1/week)")                                # Table 2: TVKG DOSTAR = 0.008105
    lkge_dostarchemo <- log(0.008845) ; label("log TVKG DOSTAR+CHEMO (TRT = 7; 1/week)")                          # Table 2: TVKG DOSTAR+CHEMO = 0.008845
    lkge_pembrochemo <- log(0.01153)  ; label("log TVKG PEMBRO+CHEMO (TRT = 8; 1/week)")                          # Table 2: TVKG PEMBRO+CHEMO = 0.01153
    lkge_cobo100b    <- log(0.01372)  ; label("log TVKG COBO 100 mg + DOSTAR cohort B (TRT = 9; 1/week)")         # Table 2: TVKG COBO 100mg+DOSTAR (cohort B) = 0.01372
    lkge_cobo300b    <- log(0.01257)  ; label("log TVKG COBO 300 mg + DOSTAR cohort B (TRT = 10; 1/week)")        # Table 2: TVKG COBO 300mg+DOSTAR (cohort B) = 0.01257
    lkge_cobo900b    <- log(0.01603)  ; label("log TVKG COBO 900 mg + DOSTAR cohort B (TRT = 11; 1/week)")        # Table 2: TVKG COBO 900mg+DOSTAR (cohort B) = 0.01603
    lkge_cobo300d    <- log(0.01267)  ; label("log TVKG COBO 300 mg + DOSTAR cohort D (TRT = 12; 1/week)")        # Table 2: TVKG COBO 300mg+DOSTAR (cohort D) = 0.01267

    # ----- Per-arm typical tumor shrinkage rate constants kse (1/week), log-transformed -----
    # All from Struemper 2025 Table 2 (per-arm TVKS estimates from the single joint NONMEM fit).
    lkse_pembro      <- log(0.01700)   ; label("log TVKS PEMBRO (TRT = 1; 1/week)")                               # Table 2: TVKS PEMBRO = 0.01700
    lkse_felad       <- log(3.396e-06) ; label("log TVKS FELAD (TRT = 2; 1/week); poorly identified (95% CI [1.941e-15, 5943])") # Table 2: TVKS FELAD = 3.396e-06 (negligible shrinkage; see vignette Errata)
    lkse_chemo       <- log(0.03566)   ; label("log TVKS CHEMO (TRT = 3; 1/week)")                                # Table 2: TVKS CHEMO = 0.03566
    lkse_feladchemo  <- log(0.02805)   ; label("log TVKS FELAD+CHEMO (TRT = 4; 1/week)")                          # Table 2: TVKS CHEMO+FELAD = 0.02805
    lkse_iocombo     <- log(0.007894)  ; label("log TVKS IO-COMBO (TRT = 5; 1/week)")                             # Table 2: TVKS IO-COMBO = 0.007894
    lkse_dostar      <- log(0.01615)   ; label("log TVKS DOSTAR (TRT = 6; 1/week)")                               # Table 2: TVKS DOSTAR = 0.01615
    lkse_dostarchemo <- log(0.03504)   ; label("log TVKS DOSTAR+CHEMO (TRT = 7; 1/week)")                         # Table 2: TVKS DOSTAR+CHEMO = 0.03504
    lkse_pembrochemo <- log(0.03499)   ; label("log TVKS PEMBRO+CHEMO (TRT = 8; 1/week)")                         # Table 2: TVKS PEMBRO+CHEMO = 0.03499
    lkse_cobo100b    <- log(0.005006)  ; label("log TVKS COBO 100 mg + DOSTAR cohort B (TRT = 9; 1/week)")        # Table 2: TVKS COBO 100mg+DOSTAR (cohort B) = 0.005006
    lkse_cobo300b    <- log(0.01125)   ; label("log TVKS COBO 300 mg + DOSTAR cohort B (TRT = 10; 1/week)")       # Table 2: TVKS COBO 300mg+DOSTAR (cohort B) = 0.01125
    lkse_cobo900b    <- log(0.006273)  ; label("log TVKS COBO 900 mg + DOSTAR cohort B (TRT = 11; 1/week)")       # Table 2: TVKS COBO 900mg+DOSTAR (cohort B) = 0.006273
    lkse_cobo300d    <- log(0.00448)   ; label("log TVKS COBO 300 mg + DOSTAR cohort D (TRT = 12; 1/week)")       # Table 2: TVKS COBO 300mg+DOSTAR (cohort D) = 0.00448

    # ----- TS covariate effects (Struemper 2025 Table 3 footnotes a, b, c) -----
    e_ntarget_ge3_tsb <- fixed(0.288)   ; label("Multiplicative-additive effect of NTARGET_GE3 on TSb (single-step binary; per-lesion linear coefficient from paper Table 3 footnote a reused on the binarised indicator)")  # Table 3: NTARGET on TSb = 0.288; paper form (1 + (NTARGET - 3) * 0.288), library form (1 + 0.288 * NTARGET_GE3)
    e_pdl1_tum_ks     <- fixed(0.00902) ; label("Exponential effect of PDL1_TUM (%) on kse (gated by has_pd1)")                                                                                                                # Table 3 footnote b: KS = TVKS * EXP(PDL1 * 0.00902) * EXP(ETA(KS))
    e_ldh_kge         <- fixed(0.474)   ; label("Power exponent of (LDH / 225.5) on kge")                                                                                                                                       # Table 3 footnote c: KG = TVKG * (LDHBL / 225.5)^0.474 * EXP(ETA(KG))

    # ----- OS sub-model parameters (Struemper 2025 Table 3) -----
    mu_os_int_logdays <- 6.94    ; label("Intercept of mu_OS in log(days) at reference subject")                # Table 3: Intercept mu_OS,int = 6.94 (95% CI 6.75, 7.13)
    sigma_os_log      <- -0.331  ; label("log(sigma_OS) -- log standard deviation of the log(t) AFT log-normal survival distribution")  # Table 3: Log SD of log(mu_OS), sigma_OS = -0.331

    emax_kge <- fixed(-6.91) ; label("Emax of kge on mu_OS (FIX; paper-equivalent 1000-fold reduction in mu_OS)")  # Table 3 footnote d: Maximum effect of kg = -6.91 [FIX]
    ec50_kge <- 0.109        ; label("EC50 (kge associated with 50% of Emax on mu_OS) (1/week)")                   # Table 3 footnote d: kg associated with 50% of maximum effect = 0.109 (95% CI 0.0764, 0.142)

    e_tsb_mu_os  <- fixed(-0.00366) ; label("Slope of TSb on mu_OS (centred at TSb = 100 mm)")                  # Table 3 footnote e: CV2 = (TSb - 100) * -0.00366
    e_ttg_mu_os  <- fixed(0.0124)   ; label("Slope of TTG on mu_OS (centred at TTG = 39.9 weeks)")              # Table 3 footnote f: CV3 = (TTG - 39.9) * 0.0124
    e_alb_mu_os  <- fixed(0.0452)   ; label("Slope of ALB on mu_OS (centred at ALB = 39.4 g/L)")                # Table 3 footnote g: CV4 = (ALBBL - 39.4) * 0.0452
    e_tpro_mu_os <- fixed(0.0194)   ; label("Slope of TPRO on mu_OS (centred at TPRO = 71 g/L)")                # Table 3 footnote h: CV5 = (TPROBL - 71) * 0.0194
    e_nlr_mu_os  <- fixed(-0.0141)  ; label("Slope of NLR on mu_OS (centred at NLR = 4.1)")                     # Table 3 footnote i: CV6 = (NLRBL - 4.1) * -0.0141

    # ----- Inter-individual variability (Struemper 2025 Table 2; correlated 3x3 block on TSb, kge, kse) -----
    # Lower-triangular packing for `c(...)`:
    #   omega^2 TSb           = 0.2645     (Table 2)
    #   omega_xy TSb x KG     = 0.01911    (Table 2)
    #   omega^2 KG            = 0.8715     (Table 2)
    #   omega_xy TSb x KS     = -0.06367   (Table 2)
    #   omega_xy KG x KS      = 0.1504     (Table 2; the source PDF typesets this entry
    #                                       as 'omega xy KG x KG' -- a typesetting error
    #                                       since omega_xy(KG,KG) would equal omega^2 KG;
    #                                       we treat it as omega_xy(KG,KS) per the
    #                                       correlated 3x3 block context.)
    #   omega^2 KS            = 0.9174     (Table 2)
    # Box-Cox transformation on eta_kge (shape -0.3744, all studies except INTR@PID) is
    # intentionally OMITTED here in favour of a standard log-normal IIV. Deviation
    # documented in vignette Errata.
    etalrbase + etalkge + etalkse ~ c(0.2645,
                                      0.01911, 0.8715,
                                     -0.06367, 0.1504, 0.9174)

    # ----- Residual error on TS (Struemper 2025 Table 2: combined additive + proportional) -----
    addSd  <- sqrt(6.324)   ; label("Additive residual SD on TS (mm)")           # Table 2: sigma^2 additive = 6.324 -> SD = sqrt(6.324)
    propSd <- sqrt(0.01276) ; label("Proportional residual SD on TS (fraction)")  # Table 2: sigma^2 proportional = 0.01276 -> SD = sqrt(0.01276)
  })

  model({
    # ----- Per-arm typical TS rate constants selected by TRT (arithmetic indicator form) -----
    # Each (TRT == k) returns 1 or 0; the sum picks the matching lkge_/lkse_ value.
    lkge_arm <- (TRT == 1)  * lkge_pembro      +
                (TRT == 2)  * lkge_felad       +
                (TRT == 3)  * lkge_chemo       +
                (TRT == 4)  * lkge_feladchemo  +
                (TRT == 5)  * lkge_iocombo     +
                (TRT == 6)  * lkge_dostar      +
                (TRT == 7)  * lkge_dostarchemo +
                (TRT == 8)  * lkge_pembrochemo +
                (TRT == 9)  * lkge_cobo100b    +
                (TRT == 10) * lkge_cobo300b    +
                (TRT == 11) * lkge_cobo900b    +
                (TRT == 12) * lkge_cobo300d

    lkse_arm <- (TRT == 1)  * lkse_pembro      +
                (TRT == 2)  * lkse_felad       +
                (TRT == 3)  * lkse_chemo       +
                (TRT == 4)  * lkse_feladchemo  +
                (TRT == 5)  * lkse_iocombo     +
                (TRT == 6)  * lkse_dostar      +
                (TRT == 7)  * lkse_dostarchemo +
                (TRT == 8)  * lkse_pembrochemo +
                (TRT == 9)  * lkse_cobo100b    +
                (TRT == 10) * lkse_cobo300b    +
                (TRT == 11) * lkse_cobo900b    +
                (TRT == 12) * lkse_cobo300d

    # has_pd1 = 1 if the treatment arm contains pembrolizumab or dostarlimab (TRT 1, 6, 7, 8, 9, 10, 11, 12).
    # The PDL1_TUM effect on kse is gated by has_pd1 per Struemper 2025 Methods: "PD-L1 tumor expression
    # was tested only as a covariate for participants who received a PD-1 inhibitor".
    has_pd1 <- (TRT == 1) + (TRT >= 6)

    # ----- Individual TS parameters -----
    TSb <- exp(lrbase + etalrbase) * (1 + e_ntarget_ge3_tsb * NTARGET_GE3)   # Table 3 footnote a (binarised at NTARGET_GE3; deviation per count-covariate policy)
    kge <- exp(lkge_arm + etalkge) * (LDH / 225.5) ^ e_ldh_kge               # Table 3 footnote c
    kse <- exp(lkse_arm + etalkse) * exp(PDL1_TUM * e_pdl1_tum_ks * has_pd1) # Table 3 footnote b (gated by has_pd1)

    # ----- TS sub-model: Stein bi-exponential, TS(t) = TSb * (exp(kge*t) + exp(-kse*t) - 1) -----
    # Encoded as two ODE compartments with per-subject initials:
    #   growth(t) = TSb * exp(kge * t),  shrink(t) = TSb * exp(-kse * t),  TS = growth + shrink - TSb.
    # At t = 0: TS = TSb + TSb - TSb = TSb (per Stein 2008 Eq. 1 boundary).
    d/dt(growth) <- kge * growth
    d/dt(shrink) <- -kse * shrink
    growth(0) <- TSb
    shrink(0) <- TSb
    TS <- growth + shrink - TSb

    # ----- OS sub-model: AFT log-normal -----
    # TTG = (log(kse) - log(kge)) / (kge + kse), clipped at 0 per Table S2.
    ttg_raw <- (log(kse) - log(kge)) / (kge + kse)
    ttg     <- (ttg_raw > 0) * ttg_raw

    # Covariate-and-link contributions to mu_OS (Table 3 footnotes d-i, summed per footnote j).
    cv1 <- emax_kge * kge / (kge + ec50_kge)   # Table 3 footnote d
    cv2 <- (TSb - 100)   * e_tsb_mu_os         # Table 3 footnote e
    cv3 <- (ttg - 39.9)  * e_ttg_mu_os         # Table 3 footnote f
    cv4 <- (ALB - 39.4)  * e_alb_mu_os         # Table 3 footnote g
    cv5 <- (TPRO - 71)   * e_tpro_mu_os        # Table 3 footnote h
    cv6 <- (NLR - 4.1)   * e_nlr_mu_os         # Table 3 footnote i

    # Paper mu_OS is in log(days); convert to log(weeks) so the entire model runs on weeks.
    mu_os_logdays <- mu_os_int_logdays + cv1 + cv2 + cv3 + cv4 + cv5 + cv6
    mu_os_logwks  <- mu_os_logdays - log(7)
    sigma_os      <- exp(sigma_os_log)

    # Log-normal AFT survival probability S(t) and hazard h(t).
    # del_t avoids log(0) at t = 0; sur is bounded below by a tiny epsilon to keep hazard finite.
    del_t  <- 1e-6
    t_pos  <- t + del_t
    z_os   <- (log(t_pos) - mu_os_logwks) / sigma_os
    sur    <- 1 - pnorm(z_os)
    pdf_lt <- exp(-0.5 * z_os * z_os) / (t_pos * sigma_os * sqrt(2 * pi))
    hazard <- pdf_lt / (sur + 1e-30)

    # ----- TS observable + combined residual error -----
    TS ~ add(addSd) + prop(propSd)
  })
}
attr(Struemper_2025_tumorsize_OS_nsclc, "message") <-
  "Joint TS-OS framework model for advanced/metastatic NSCLC (Struemper 2025, n=786, 7 GSK trials). TS = Stein bi-exponential with per-treatment-arm kge/kse (selected by TRT integer 1-12); OS = AFT log-normal with treatment-agnostic link (kge Emax + TSb + TTG linear) plus baseline ALB, TPRO, NLR. Outputs: TS (observable, mm), mu_os_logwks, sigma_os, sur (S(t)), hazard. PDL1_TUM effect on kse is gated by a derived has_pd1 indicator (1 if TRT == 1 or TRT >= 6, else 0). NTARGET (count) is binarised to NTARGET_GE3 per the library's count-covariate policy; Box-Cox IIV on kge (paper shape -0.3744) is omitted in favour of standard log-normal IIV. Both deviations are documented in the vignette Errata."
Struemper_2025_tumorsize_OS_nsclc
