Thorsted_2019_endotoxin <- function() {
  description <- paste(
    "Preclinical (piglet). Semi-mechanistic model of intravenously infused",
    "E. coli O111:B4 endotoxin (ETX) and the innate cytokine response it",
    "induces in anaesthetised piglets (Thorsted 2019). ETX follows",
    "one-compartment disposition with Michaelis-Menten elimination plus a",
    "constant assay baseline. ETX stimulates TNF-alpha production through a",
    "sigmoidal Emax whose EC50 is itself raised by ETX through a three-stage",
    "moderator chain (endotoxin tolerance). IL-6 production is stimulated",
    "linearly by the relative increase in TNF-alpha and by ETX, with",
    "appearance in plasma delayed by four transit compartments. A mixture",
    "sub-population (MIX_ETX_CONTAM) carries pre-study contaminating ETX",
    "that is measured but does not drive the cytokines, so the model keeps",
    "a second, identically eliminated ETX pool (central_contam) for the",
    "assay read-out. Every infusion must be given to both central and",
    "central_contam."
  )
  reference <- paste(
    "Thorsted A, Bouchene S, Tano E, Castegren M, Lipcsey M, Sjolin J,",
    "Karlsson MO, Friberg LE, Nielsen EI. (2019). A non-linear mixed effect",
    "model for innate immune response: In vivo kinetics of endotoxin and its",
    "induction of the cytokines tumor necrosis factor alpha and interleukin-6.",
    "PLoS ONE 14(2):e0211981. doi:10.1371/journal.pone.0211981."
  )
  vignette <- "Thorsted_2019_endotoxin"
  units <- list(
    time = "h",
    dosing = "EU",
    concentration = "EU/L (endotoxin); ng/L (TNF-alpha, IL-6)"
  )

  covariateData <- list(
    MIX_ETX_CONTAM = list(
      description = paste(
        "Per-subject binary mixture-model class indicator. 1 = subject in",
        "the sub-population whose ETX compartment starts at a contaminating",
        "amount BASE_ETX * Propcont (S1 Text $MIX sub-population 2,",
        "P = 0.399); 0 = subject whose ETX compartment starts at zero."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no pre-study contamination)",
      notes = paste(
        "Not a measured covariate: the latent class of the NONMEM $MIX block",
        "(NSPOP = 2, P(2) = THETA(7) = 0.399, fixed in the final run). The",
        "contamination enters only the measured ETX (Cc); it does not",
        "stimulate TNF-alpha or IL-6. For population simulation draw",
        "MIX_ETX_CONTAM ~ Bernoulli(0.399) per subject."
      ),
      source_name = "MIXNUM"
    ),
    STUDY_ETX_C = list(
      description = "1 = piglet from study C of Thorsted 2019 Table 1 (6 h ETX infusion dose-response, ETX batch 3000 EU/ug); 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0. All four STUDY_ETX_* = 0 means study A or B.",
      notes = "Selects the BioSource TNF-alpha baseline (3.06 ng/L) and the BioSource/Quantikine IL-6 baseline (8.19 ng/L). S1 Text STUDY == 3.",
      source_name = "STUDY"
    ),
    STUDY_ETX_D = list(
      description = "1 = piglet from study D of Thorsted 2019 Table 1 (3 h accelerating / decelerating ETX infusion, batch 3000 EU/ug); 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0. All four STUDY_ETX_* = 0 means study A or B.",
      notes = "Selects the study-D BioSource TNF-alpha baseline (68.2 ng/L) and the Quantikine IL-6 baseline (8.19 ng/L). S1 Text STUDY == 4.",
      source_name = "STUDY"
    ),
    STUDY_ETX_E = list(
      description = "1 = piglet from study E of Thorsted 2019 Table 1 (1-3 h ETX infusion, batch 1000 EU/ug); 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0. All four STUDY_ETX_* = 0 means study A or B.",
      notes = "Selects the study-E DuoSet TNF-alpha baseline (25.4 ng/L); IL-6 uses the DuoSet baseline shared with studies A and B. S1 Text STUDY == 5.",
      source_name = "STUDY"
    ),
    STUDY_ETX_F = list(
      description = "1 = piglet from study F of Thorsted 2019 Table 1 (6 h ETX infusion with a higher first 30 min, batch 3000 EU/ug); 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0. All four STUDY_ETX_* = 0 means study A or B.",
      notes = "Selects the Coatest ETX assay baseline (1810 EU/L), the BioSource TNF-alpha baseline (3.06 ng/L) and the Quantikine IL-6 baseline (8.19 ng/L). S1 Text STUDY == 6.",
      source_name = "STUDY"
    )
  )

  compartmentData <- list(
    central = list(analyte = "endotoxin (infused)", units = "EU", specimen = "plasma", verified = TRUE),
    central_contam = list(
      analyte = "endotoxin (infused + pre-study contamination)",
      units = "EU",
      specimen = "plasma",
      verified = TRUE
    ),
    tnf = list(analyte = "TNF-alpha", units = "ng/L", specimen = "plasma", verified = TRUE),
    moderator1 = list(
      analyte = "EC50 of the ETX-TNF-alpha relationship, stage 1",
      units = "EU/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator2 = list(
      analyte = "EC50 of the ETX-TNF-alpha relationship, stage 2",
      units = "EU/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    moderator3 = list(
      analyte = "EC50 of the ETX-TNF-alpha relationship (operative)",
      units = "EU/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "IL-6 appearance delay, stage 1",
      units = "ng/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "IL-6 appearance delay, stage 2",
      units = "ng/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "IL-6 appearance delay, stage 3",
      units = "ng/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "IL-6 appearance delay, stage 4",
      units = "ng/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    il6 = list(analyte = "IL-6", units = "ng/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "pig (anaesthetised piglets)",
    n_subjects = 116L,
    n_studies = 6L,
    age_range = "approximately 12 weeks",
    weight_median = "mean approximately 27 kg",
    disease_state = "experimental endotoxaemia (porcine endotoxin model); saline-infused controls included",
    dose_range = paste(
      "Continuous IV infusion of E. coli O111:B4 endotoxin for 1-30 h at",
      "0 (saline) to 16.0 ug/kg/h; batch potency 1000 EU/ug (studies A, B, E)",
      "or 3000 EU/ug (studies C, D, F)"
    ),
    regions = "Sweden (Uppsala University Hospital)",
    notes = paste(
      "Six studies (A-F, Thorsted 2019 Table 1) pooled: 278 ETX, 1068",
      "TNF-alpha and 1045 IL-6 observations. ETX was measured in studies C,",
      "E and F only. A seventh study (G, n = 6) was held out for external",
      "evaluation and is not part of the fit. Sex not reported."
    )
  )

  ini({
    # Endotoxin disposition (Thorsted 2019 Table 2; structure from S1 Text)
    lvmax <- log(442000); label("Maximum ETX elimination capacity Vmax (EU/h)") # Table 2 'Vmax' 442000 EU/h (RSE 1%)
    lkm <- log(12600); label("ETX concentration at half of Vmax, KM (EU/L)") # Table 2 'KM' 12600 EU/L (RSE 17%)
    lvc <- log(36.1); label("Volume of distribution of ETX, Vc (L)") # Table 2 'Vc' 36.1 L (RSE 20%)
    lrbase_etx <- log(155); label("Constant ETX assay baseline, Endochrome assay (EU/L)") # Table 2 'BASE ETX (Endochrome, study C+E)' 155 EU/L (RSE 20%)
    lrbase_etx_f <- log(1810); label("Constant ETX assay baseline, Coatest assay, study F (EU/L)") # Table 2 'BASE ETX (Coatest, study F)' 1810 EU/L (RSE 38%)
    lprop_cont <- log(132); label("Proportionality of BASE_ETX to the contaminating initial ETX amount, Propcont (unitless)") # Table 2 'Propcont' 132 (RSE 15%); Eq 1

    # TNF-alpha turnover with ETX-stimulated production
    lmtt_tnf <- log(1.04); label("Mean transit time of TNF-alpha in plasma, 1/kdeg (h)") # Table 2 'MTT TNF-alpha' 1.04 h (RSE 13%)
    lrbase_tnf <- log(281); label("Baseline TNF-alpha, DuoSet assay, studies A+B (ng/L)") # Table 2 'S0,TNF-alpha (DuoSet, study A+B)' 281 ng/L (RSE 9%)
    lrbase_tnf_e <- log(25.4); label("Baseline TNF-alpha, DuoSet assay, study E (ng/L)") # Table 2 'S0,TNF-alpha (DuoSet, study E)' 25.4 ng/L (RSE 12%)
    lrbase_tnf_cf <- log(3.06); label("Baseline TNF-alpha, BioSource assay, studies C+F (ng/L)") # Table 2 'S0,TNF-alpha (BioSource, study C+F)' 3.06 ng/L (RSE 29%)
    lrbase_tnf_d <- log(68.2); label("Baseline TNF-alpha, BioSource assay, study D (ng/L)") # Table 2 'S0,TNF-alpha (BioSource, study D)' 68.2 ng/L (RSE 15%)
    lemax_tnf <- log(2540); label("Maximum fractional increase of TNF-alpha production, Emax (unitless)") # Table 2 'Emax' 2540 (RSE 9%)
    lec50_tnf <- log(286); label("ETX concentration at half of Emax before tolerance, EC0 (EU/L)") # Table 2 'EC0' 286 EU/L (RSE 5%)
    lhill_tnf <- log(2.10); label("Sigmoidicity coefficient of the ETX-TNF-alpha Emax, gamma (unitless)") # Table 2 'gamma' 2.10 (RSE 26%)

    # Endotoxin tolerance: ETX-driven rise of EC50 through a three-stage chain
    lmtt_tol <- log(6.33); label("Mean transit time of the EC50 moderator chain, MTT_EC50 (h)") # Table 2 'MTT EC50' 6.33 h (RSE 5%)
    lemax_tol <- log(45100); label("Maximum fractional increase of EC0, Tmax (unitless)") # Table 2 'Tmax' 45100 (RSE 2%)
    lec50_tol <- log(29300); label("ETX concentration at half of Tmax, TC50 (EU/L)") # Table 2 'TC50' 29300 EU/L (RSE 3%)

    # IL-6 turnover with TNF-alpha- and ETX-stimulated production
    lmtt_il6 <- log(1.45); label("Mean transit time of IL-6 through the five-stage chain, MTT_IL-6 (h)") # Table 2 'MTT IL-6' 1.45 h (RSE 11%)
    lrbase_il6 <- log(49.2); label("Baseline IL-6, DuoSet assay, studies A+B+E (ng/L)") # Table 2 'S0,IL-6 (DuoSet, study A+B+E)' 49.2 ng/L (RSE 13%)
    lrbase_il6_cdf <- log(8.19); label("Baseline IL-6, BioSource/Quantikine assay, studies C+D+F (ng/L)") # Table 2 'S0,IL-6 (BioSource, study C+D+F)' 8.19 ng/L (RSE 10%)
    lslope_tnf <- log(0.961); label("Slope of the relative TNF-alpha increase on IL-6 production, ESLP,TNF (unitless)") # Table 2 'ESLP,TNF-alpha' 0.961 (RSE 9%); Eq 2
    lslope_etx <- log(0.00937); label("Slope of ETX concentration on IL-6 production, ESLP,ETX (L/EU)") # Table 2 'ESLP,ETX' 0.00937 (RSE 28%)

    # IIV: Table 2 CV% = 100 * sqrt(omega^2), so omega^2 = (CV/100)^2.
    # KM-Vc correlation -0.42 (Results 'Infused endotoxin'):
    # cov = -0.42 * 1.20 * 0.753 = -0.379512.
    etalkm + etalvc ~ c(1.44, -0.379512, 0.567009) # Table 2 'KM' 120 CV%, 'Vc' 75.3 CV%; r = -0.42
    etalrbase_etx ~ 0.927369 # Table 2 'BASE ETX' 96.3 CV% (shared by both assay baselines)
    etalmtt_tnf ~ 0.664225 # Table 2 'MTT TNF-alpha' 81.5 CV%
    etalrbase_tnf ~ 0.292681 # Table 2 'S0,TNF-alpha' 54.1 CV% (shared by all four baselines)
    etalmtt_il6 ~ 0.224676 # Table 2 'MTT IL-6' 47.4 CV%
    etalrbase_il6 ~ 0.638401 # Table 2 'S0,IL-6' 79.9 CV% (shared by both baselines)

    # Residual error: additive on log-transformed data (Methods 'Modelling process')
    expSd <- 0.321; label("Log-scale residual SD, ETX (unitless)") # Table 2 'sigma ETX' 32.1%
    expSd_tnf <- 0.485; label("Log-scale residual SD, TNF-alpha (unitless)") # Table 2 'sigma TNF-alpha' 48.5%
    expSd_il6 <- 0.505; label("Log-scale residual SD, IL-6 (unitless)") # Table 2 'sigma IL-6' 50.5%
  })
  model({
    # Study / assay indicators. All four 0 = study A or B (S1 Text $PK IF-branches).
    study_ab <- 1 - STUDY_ETX_C - STUDY_ETX_D - STUDY_ETX_E - STUDY_ETX_F
    study_cf <- STUDY_ETX_C + STUDY_ETX_F
    study_cdf <- STUDY_ETX_C + STUDY_ETX_D + STUDY_ETX_F

    # Endotoxin
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm)
    vc <- exp(lvc + etalvc)
    rbase_etx <- exp(lrbase_etx * (1 - STUDY_ETX_F) + lrbase_etx_f * STUDY_ETX_F + etalrbase_etx)
    prop_cont <- exp(lprop_cont)

    # TNF-alpha
    mtt_tnf <- exp(lmtt_tnf + etalmtt_tnf)
    rbase_tnf <- exp(
      lrbase_tnf * study_ab + lrbase_tnf_cf * study_cf +
        lrbase_tnf_d * STUDY_ETX_D + lrbase_tnf_e * STUDY_ETX_E + etalrbase_tnf
    )
    emax_tnf <- exp(lemax_tnf)
    ec0_tnf <- exp(lec50_tnf)
    hill_tnf <- exp(lhill_tnf)
    kout_tnf <- 1 / mtt_tnf
    kin_tnf <- rbase_tnf * kout_tnf

    # Tolerance chain: three stages, kout = 3 / MTT_EC50, started at EC0
    mtt_tol <- exp(lmtt_tol)
    emax_tol <- exp(lemax_tol)
    ec50_tol <- exp(lec50_tol)
    kout_tol <- 3 / mtt_tol
    kin_tol <- ec0_tnf * kout_tol

    # IL-6: five stages (four transits + plasma IL-6), kout = 5 / MTT_IL-6
    mtt_il6 <- exp(lmtt_il6 + etalmtt_il6)
    rbase_il6 <- exp(lrbase_il6 * (1 - study_cdf) + lrbase_il6_cdf * study_cdf + etalrbase_il6)
    slope_tnf <- exp(lslope_tnf)
    slope_etx <- exp(lslope_etx)
    kout_il6 <- 5 / mtt_il6
    kin_il6 <- rbase_il6 * kout_il6

    # Only the infused ETX (central, started at zero) drives the cytokines;
    # the assay baseline and the contamination do not (Results 'Host response').
    cp_etx <- central / vc
    ec50_tnf <- moderator3
    tol <- emax_tol * cp_etx / (ec50_tol + cp_etx)
    eff_tnf <- emax_tnf * cp_etx^hill_tnf / (ec50_tnf^hill_tnf + cp_etx^hill_tnf)
    eff_il6 <- slope_tnf * (tnf / rbase_tnf - 1) + slope_etx * cp_etx

    d/dt(central) <- -vmax * (central / vc) / (km + central / vc)
    d/dt(central_contam) <- -vmax * (central_contam / vc) / (km + central_contam / vc)
    d/dt(tnf) <- kin_tnf * (1 + eff_tnf) - kout_tnf * tnf
    d/dt(moderator1) <- kin_tol * (1 + tol) - kout_tol * moderator1
    d/dt(moderator2) <- kout_tol * moderator1 - kout_tol * moderator2
    d/dt(moderator3) <- kout_tol * moderator2 - kout_tol * moderator3
    d/dt(transit1) <- kin_il6 * (1 + eff_il6) - kout_il6 * transit1
    d/dt(transit2) <- kout_il6 * transit1 - kout_il6 * transit2
    d/dt(transit3) <- kout_il6 * transit2 - kout_il6 * transit3
    d/dt(transit4) <- kout_il6 * transit3 - kout_il6 * transit4
    d/dt(il6) <- kout_il6 * transit4 - kout_il6 * il6

    # S1 Text A_0(1) = EBASE * PROP: the baseline in EU/L times a unitless
    # factor is used directly as an amount (EU), as printed in Eq 1.
    central_contam(0) <- MIX_ETX_CONTAM * rbase_etx * prop_cont
    tnf(0) <- rbase_tnf
    moderator1(0) <- ec0_tnf
    moderator2(0) <- ec0_tnf
    moderator3(0) <- ec0_tnf
    transit1(0) <- rbase_il6
    transit2(0) <- rbase_il6
    transit3(0) <- rbase_il6
    transit4(0) <- rbase_il6
    il6(0) <- rbase_il6

    # Measured ETX = assay baseline + all ETX in the contamination-carrying pool
    Cc <- rbase_etx + central_contam / vc
    Cc ~ lnorm(expSd)
    tnf ~ lnorm(expSd_tnf)
    il6 ~ lnorm(expSd_il6)
  })
}
