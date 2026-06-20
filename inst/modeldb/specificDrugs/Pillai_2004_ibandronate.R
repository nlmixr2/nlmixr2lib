Pillai_2004_ibandronate <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) model for ibandronate (a nitrogen-containing bisphosphonate) suppression of urinary C-telopeptide of type-I collagen (uCTX) in postmenopausal women with osteoporosis. A virtual K-PD effect compartment receives the administered dose and decays at rate KDE, producing a dose-driving rate (DODR = KDE * central * F) that inhibits the uCTX synthesis rate KS via a sigmoid Emax (inhibition-fraction form). uCTX follows an indirect-response synthesis-degradation turnover (KS, KD). A multiplicative placebo disease-progression drift (1 + SLOPE * t) and a calcium + vitamin D supplementation suppression term (1 - VIT * [1 - exp(-KVIT * t)]) modify the observed uCTX. The supplementation indicator CONMED_CAVITD also switches KDE between the no-supplement (0.112 /day) and with-supplement (0.014 /day) typical values. The model handles intravenous and oral routes via the bioavailability factor F (default F=1 for IV; oral users override lfdepot to log(0.008) for 2.5 mg oral or log(0.007) for 20 mg oral, both fixed from a separate absolute-bioavailability study)."
  reference <- paste(
    "Pillai G, Gieschke R, Goggin T, Jacqmin P, Schimmer RC, Steimer JL (2004).",
    "A semimechanistic and mechanistic population PK-PD model for biomarker",
    "response to ibandronate, a new bisphosphonate for the treatment of osteoporosis.",
    "Br J Clin Pharmacol 58(6):618-631.",
    "doi:10.1111/j.1365-2125.2004.02224.x.",
    sep = " "
  )
  vignette <- "Pillai_2004_ibandronate"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mmolCR",
    outcome       = "uCTX (ug/mmolCR)"
    # units$concentration intentionally documents the urinary CTX (uCTX)
    # biomarker output rather than a plasma drug concentration: this is a
    # K-PD model with no observed plasma drug concentration. The Cc
    # observation in model() is the uCTX biomarker concentration (in
    # ug/mmolCR, creatinine-normalised) after the multiplicative
    # placebo-drift (PBO) and calcium + vitamin D supplementation (VITD)
    # adjustments. The classical 4-compartment PK-PD model (Table 3 of
    # the paper) is the accompanying-PK reference for plasma + urine
    # ibandronate concentrations; see vignette narrative.
  )

  covariateData <- list(
    CONMED_CAVITD = list(
      description        = "Concomitant calcium + vitamin D supplementation indicator (1 = on daily oral calcium + vitamin D coadministration during the observation period, 0 = not on supplementation)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no supplementation)",
      notes              = paste(
        "Per-subject baseline indicator. Drives the multiplicative covariate effect",
        "on KDE (the K-PD virtual-compartment elimination rate constant: 0.112 /day",
        "without supplements, 0.014 /day with supplements) AND gates the multiplicative",
        "VITD = 1 - VIT * (1 - exp(-KVIT * t)) suppression term so that VITD = 1 when",
        "CONMED_CAVITD = 0 (the term contributes no extra uCTX suppression) and adds",
        "up to a 31.4% additional uCTX suppression at long times when CONMED_CAVITD = 1.",
        "Paper studies MF4361 and MF4411 (Thiebaud et al. 1997; Delmas et al. 1996)",
        "used daily oral calcium + vitamin D supplementation; study MF9853 did not.",
        "Source label: 'supplemental therapy' (Pillai 2004 narrative + Table 4 legend)."
      ),
      source_name        = "supplemental therapy"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 174,
    n_studies      = 3,
    age_range      = "Postmenopausal women",
    weight_range   = NA_character_,
    sex_female_pct = 100,
    race_ethnicity = c(White = NA_real_, Asian_Japanese = NA_real_, Other = NA_real_),
    disease_state  = "Postmenopausal osteoporosis (PMO) and postmenopausal osteopenia. Two cohorts of i.v. dose-finding data (study MF9853 in 50 Japanese women with osteopenia; Thiebaud et al. 1997 in 124 Caucasian women with osteoporosis) and one cohort of oral data (Delmas et al. 1996 in 676 Caucasian women with osteoporosis); validation data drawn from Adami 2002 (n=520), Recker 2000 (n=2860), and Ravn 1996 (n=180) post-hoc.",
    dose_range     = "I.v. ibandronate 0.25, 0.5, 1, or 2 mg every 3 months (development studies MF9853 and Thiebaud 1997); oral ibandronate 2.5 mg once daily (development study Delmas 1996); validation regimens 0.25-5 mg oral daily and 0.5-2 mg i.v. q3mo.",
    regions        = "Japan (MF9853); Europe (Thiebaud 1997, Delmas 1996, Adami 2002, Recker 2000, Ravn 1996).",
    notes          = paste(
      "Development pool of 174 women (50 + 124) for i.v. K-PD and 676 women for the",
      "extended K-PD with oral; demographics not tabulated by cohort in the paper.",
      "Race / ethnicity inferred from study geography (MF9853 = Japanese; Thiebaud /",
      "Delmas / Adami / Recker / Ravn = Caucasian). See Pillai 2004 Tables 1 and 2",
      "for the full study list."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters (Pillai 2004 Table 4, p. 627)
    # All point estimates from Table 4 (K-PD final model with oral +
    # supplementation extensions). BSV (CV %) translated to log-normal
    # variance via omega^2 = log(CV^2 + 1) where the IIV is log-normal.
    # ---------------------------------------------------------------

    # uCTX indirect-response (synthesis-degradation) turnover
    lksyn <- log(255.00)  ; label("uCTX synthesis rate KS (ug mmolCR-1 day-1)")        # Pillai 2004 Table 4: KS = 255.00 ug mmolCR-1 day-1
    lkdeg <- log(1.06)    ; label("uCTX degradation rate KD (1/day)")                  # Pillai 2004 Table 4: KD = 1.06 /day

    # K-PD virtual-compartment kinetics (lkel typical value = no-supplement KDE)
    lkel   <- log(0.112)  ; label("K-PD virtual compartment elimination rate KDE, no-supplement reference (1/day)")  # Pillai 2004 Table 4: KDE (MF9853; no supplements) = 0.112 /day
    lekd50 <- log(17.20)  ; label("Dose-driving rate at half-Emax EKD50 (ug/day)")     # Pillai 2004 Table 4: EKD50 = 17.20 ug/day
    lhill  <- log(0.913)  ; label("Hill coefficient HILL (unitless)")                  # Pillai 2004 Table 4: HILL = 0.913

    # Calcium + vitamin D supplementation suppression term parameters
    lvit  <- log(0.314)   ; label("Maximum supplement-induced uCTX suppression fraction VIT (unitless)")  # Pillai 2004 Table 4: VIT = 0.314
    lkvit <- log(0.012)   ; label("Supplement-effect attainment rate constant KVIT (1/day)")              # Pillai 2004 Table 4: KVIT = 0.012 /day

    # Placebo / disease-progression drift (linear-scale because IIV is additive)
    slope_pbo <- 2.32e-4  ; label("Placebo / disease-progression slope SLOPE for uCTX, linear scale (1/day; ~0.7% per month)")  # Pillai 2004 Table 4: SLOPE = 2.32E-04 /day

    # Bioavailability (fixed from upstream Phase I/II absolute-bioavailability
    # study, Pillai 2004 ref [31]). Default value log(1) corresponds to F = 1
    # for intravenous administration; for oral simulation the user overrides
    # lfdepot to log(0.008) for 2.5 mg oral or log(0.007) for 20 mg oral (the
    # two F values Table 4 lists from the upstream popPK fit).
    lfdepot <- fixed(log(1))  ; label("Bioavailability F (default 1 = IV; override to log(0.008) for 2.5 mg oral or log(0.007) for 20 mg oral per Pillai 2004 Table 4 / ref [31])")  # Pillai 2004 Table 4: F(2.5 mg) = 0.008, F(20 mg) = 0.007 (fixed from upstream bioavailability study [31])

    # ---------------------------------------------------------------
    # Covariate effect: calcium + vitamin D supplementation on KDE
    # ---------------------------------------------------------------
    # kel = exp(lkel + e_cavitd_kel * CONMED_CAVITD)
    # When CONMED_CAVITD = 0: kel = 0.112 /day (no-supplement)
    # When CONMED_CAVITD = 1: kel = 0.014 /day (with supplements)
    e_cavitd_kel <- log(0.014 / 0.112)  ; label("Effect of CONMED_CAVITD on KDE (log-scale ratio; supplements give 8-fold slower KDE)")  # Pillai 2004 Table 4: KDE (MF4361/MF4411; with supplements) = 0.014 /day vs KDE (MF9853; no supplements) = 0.112 /day; ratio = 0.014/0.112 = 0.125

    # ---------------------------------------------------------------
    # Inter-individual variability (Pillai 2004 Table 4 BSV CV %)
    # omega^2 = log(CV^2 + 1) for log-normal IIV
    # ---------------------------------------------------------------
    etalksyn   ~ log(0.26^2 + 1)        # Pillai 2004 Table 4: BSV(KS) = 26% CV -> log-normal variance ~ 0.0654
    etalkdeg   ~ log(0.34^2 + 1)        # Pillai 2004 Table 4: BSV(KD) = 34% CV -> log-normal variance ~ 0.1094
    etalkel    ~ log(0.30^2 + 1)        # Pillai 2004 Table 4: BSV(KDE, MF9853 / no-supp) = 30% CV -> log-normal variance ~ 0.0862. NB: Table 4 also reports BSV 139% CV in the supplemented studies; the higher with-supplement variability is documented in the vignette but the single etalkel variance is taken from the better-controlled no-supplement cohort per the paper's narrative.
    etalekd50  ~ log(0.83^2 + 1)        # Pillai 2004 Table 4: BSV(EKD50) = 83% CV -> log-normal variance ~ 0.5240
    etalvit    ~ log(0.50^2 + 1)        # Pillai 2004 Table 4: BSV(VIT) = 50% CV -> log-normal variance ~ 0.2231
    etalkvit   ~ log(1.92^2 + 1)        # Pillai 2004 Table 4: BSV(KVIT) = 192% CV -> log-normal variance ~ 1.5446
    etalfdepot ~ log(0.52^2 + 1)        # Pillai 2004 Table 4: BSV(F, 2.5 mg) = 52% CV -> log-normal variance ~ 0.2393

    # Additive IIV on placebo slope (paper Table 4 footnote: "BSV for slope is
    # modelled as additive and reported as variance"). Variance is on the
    # linear scale of slope_pbo; the eta is added (not multiplied) to slope_pbo
    # in model(). Hill coefficient has no BSV in Table 4 (entry is "-").
    etaslope_pbo ~ 3.16e-7              # Pillai 2004 Table 4: BSV(SLOPE) = 3.16E-07 (additive variance, linear scale)

    # ---------------------------------------------------------------
    # Residual error (Pillai 2004 Table 4)
    # ---------------------------------------------------------------
    propSd <- 0.33        ; label("Proportional residual SD on uCTX (CV fraction)")  # Pillai 2004 Table 4: residual error 33% CV on uCTX concentration
  })

  model({
    # ---------------------------------------------------------------
    # 1. Individual parameters
    # ---------------------------------------------------------------
    ksyn  <- exp(lksyn  + etalksyn)
    kdeg  <- exp(lkdeg  + etalkdeg)
    kel   <- exp(lkel   + etalkel + e_cavitd_kel * CONMED_CAVITD)
    ekd50 <- exp(lekd50 + etalekd50)
    hill  <- exp(lhill)
    vit   <- exp(lvit   + etalvit)
    kvit  <- exp(lkvit  + etalkvit)
    slope_ind <- slope_pbo + etaslope_pbo
    fdepot <- exp(lfdepot + etalfdepot)

    # ---------------------------------------------------------------
    # 2. Bioavailability applied at dose load. Default lfdepot = log(1)
    #    corresponds to F = 1 for IV; the user overrides lfdepot to
    #    log(0.008) (2.5 mg oral) or log(0.007) (20 mg oral) when
    #    simulating oral dosing.
    # ---------------------------------------------------------------
    f(central) <- fdepot

    # ---------------------------------------------------------------
    # 3. K-PD virtual compartment and dose-driving rate.
    #    central holds the K-PD virtual amount in mg (paper's A1);
    #    EKD50 in Table 4 is in ug/day, so DODR is scaled by 1000 to
    #    align units (mg/day * 1000 = ug/day). The inhibition factor
    #    INH is the KS-fraction-remaining form per the paper's
    #    appendix: INH = 1 - DODR^HILL / [EKD50^HILL + DODR^HILL],
    #    so INH = 1 when DODR = 0 (no drug, full KS) and INH -> 0 as
    #    DODR -> infinity (full KS inhibition).
    # ---------------------------------------------------------------
    dodr_ug <- kel * central * 1000
    inh_factor <- 1 - dodr_ug^hill / (ekd50^hill + dodr_ug^hill)

    # ---------------------------------------------------------------
    # 4. ODE system. central decays at KDE (rate kel); the uCTX pool
    #    is a turnover state with inhibited synthesis (ksyn * INH)
    #    and first-order degradation (kdeg). uCTX baseline is the
    #    no-drug steady state ksyn / kdeg (Pillai 2004 Table 3
    #    reports uCTXss = 339 = 231.43 / 0.68 for the PK-PD model;
    #    Table 4 K-PD: ksyn / kdeg = 255.00 / 1.06 = 240.6 ug/mmolCR).
    # ---------------------------------------------------------------
    d/dt(central) <- -kel * central
    d/dt(effect)  <-  ksyn * inh_factor - kdeg * effect
    effect(0)     <-  ksyn / kdeg

    # ---------------------------------------------------------------
    # 5. Multiplicative placebo drift and supplement-driven suppression
    #    on observed uCTX. PBO = 1 + slope_ind * t is the time-linear
    #    placebo drift (~0.7% / month per Pillai 2004 Discussion).
    #    VITD = 1 - VIT * [1 - exp(-KVIT * t)] is the supplement
    #    suppression term, gated by CONMED_CAVITD so that VITD = 1
    #    (no extra suppression) when CONMED_CAVITD = 0. Combined
    #    observation: uCTX(composite) = uCTX * PBO * VITD per the
    #    paper's appendix.
    # ---------------------------------------------------------------
    pbo  <- 1 + slope_ind * t
    vitd <- 1 - vit * (1 - exp(-kvit * t)) * CONMED_CAVITD

    # ---------------------------------------------------------------
    # 6. Observation and residual error. Cc is named per the
    #    nlmixr2lib single-output canonical and represents the
    #    composite uCTX biomarker concentration (ug/mmolCR,
    #    creatinine-normalised) after PBO and VITD adjustments, NOT a
    #    plasma drug concentration. The K-PD model has no plasma-PK
    #    observation; the classical 4-cmt PK-PD model in the paper
    #    (Table 3) is the accompanying-PK reference. Pillai 2004
    #    Table 4 reports residual error as 33% CV on uCTX.
    # ---------------------------------------------------------------
    Cc <- effect * pbo * vitd
    Cc ~ prop(propSd)
  })
}
