Chigutsa_2026_abemaciclib <- function() {
  description <- paste(
    "Joint pharmacokinetic - tumor size - progression-free survival (PFS)",
    "model for oral abemaciclib plus fulvestrant (or placebo plus",
    "fulvestrant) in 663 women with HR+/HER2- metastatic breast cancer",
    "from MONARCH 2 (Chigutsa 2026). PK layer: the semi-mechanistic",
    "parent model of Chigutsa 2020 (parallel lagged zero-order and",
    "four-transit saturable absorption with gut metabolism, well-stirred",
    "liver compartment, two-compartment distribution; metabolites not",
    "carried). Tumor layer: zero-order growth of proliferating cells,",
    "first-order damage driven by a fulvestrant rate plus a linear",
    "abemaciclib-concentration effect that wanes exponentially after a",
    "delay (arm-specific resistance half-life), and three damaged-cell",
    "transit compartments; tumor size is the sum of all four states.",
    "PFS layer: constant baseline hazard (separate for measurable and",
    "non-measurable disease) multiplied by exp of a change-from-baseline",
    "tumor-size term (two-class mixture on its coefficient) minus a",
    "waning linear concentration effect, with liver-metastasis and ECOG",
    "effects on the baseline hazard."
  )
  reference <- paste(
    "Chigutsa E, Chapman SC, Turner PK.",
    "Longitudinal Tumor Size and Survival Modeling for Exposure-Response",
    "Analysis of Drugs with Frequent Dose Reductions: Dose Justification of",
    "Abemaciclib in Patients with Metastatic Breast Cancer.",
    "Clin Pharmacol Ther. 2026;119(4):1070.",
    "doi:10.1002/cpt.70212.",
    "PK layer: Chigutsa E, Kambhampati SRP, Karen Sykes A, Posada MM,",
    "van der Walt JS, Turner PK. Development and application of a",
    "mechanistic population modeling approach to describe abemaciclib",
    "pharmacokinetics. CPT Pharmacometrics Syst Pharmacol.",
    "2020;9(9):523-533. doi:10.1002/psp4.12544."
  )
  vignette <- "Chigutsa_2026_abemaciclib"
  units <- list(time = "week", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling (reference 70 kg): exponent 1 on the central and",
        "peripheral volumes, 0.75 on intercompartmental clearance, hepatic",
        "blood flow and intrinsic clearance; also enters body-surface area",
        "for the liver volume (Chigutsa 2020 supplementary control stream)."
      ),
      source_name = "WTE"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the DuBois body-surface area that sets the liver volume",
        "(Urata 1995). The Chigutsa 2026 control stream imputes a missing",
        "height as 160 cm."
      ),
      source_name = "HTE"
    ),
    AE_DIARRHEA = list(
      description = "Diarrhea indicator (onset-and-persist)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Additive shift of -1.28 on the logit of the fraction absorbed (Fa).",
        "In the source model the switch turns on at the first diarrhea event",
        "(or at baseline for pre-existing diarrhea) and stays on for the rest",
        "of the record, so carry AE_DIARRHEA = 1 from first onset onwards."
      ),
      source_name = "SWITCH (from MOTIL / MDSFL)"
    ),
    FORM_ABEMACICLIB_C2 = list(
      description = "Abemaciclib 50 mg 50% w/w capsule (C2) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = C2 capsule; 0 = otherwise. Multiplies the absorption lag time",
        "by (1 + 0.822). Reference formulations are the drug-in-capsule (C1)",
        "and the oral solution."
      ),
      source_name = "FORMU == 3"
    ),
    FORM_ABEMACICLIB_C3 = list(
      description = "Abemaciclib 50 mg 25% w/w capsule (C3, commercial-equivalent) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = C3 capsule; 0 = otherwise. Multiplies the absorption lag time by",
        "(1 + 1.06) and the maximum saturable absorption rate by (1 - 0.16).",
        "MONARCH 2 used the C3 capsule (Chigutsa 2020 Table S1)."
      ),
      source_name = "FORMU == 4"
    ),
    TRT_ABEMACICLIB = list(
      description = "Abemaciclib treatment-arm indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = randomized to abemaciclib + fulvestrant; 0 = placebo +",
        "fulvestrant. Selects the arm-specific half-life of the waning",
        "tumor drug effect (8.44 vs 2.51 weeks). Abemaciclib exposure itself",
        "comes from the dosing records."
      ),
      source_name = "TRT"
    ),
    RECIST_NONMEASURABLE = list(
      description = "Non-measurable (bone-only) disease indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = no RECIST-measurable target lesions (bone-only disease; 186 of",
        "663 patients); 0 = measurable disease with tumor-size data. For",
        "RECIST_NONMEASURABLE = 1 the baseline hazard is 0.0103/week, the",
        "tumor-size term is removed and the drug effect on the hazard does",
        "not wane."
      ),
      source_name = "NO_MEASR"
    ),
    MIX_HIGH_TS_HAZARD = list(
      description = "Mixture-model class: stronger tumor-size-change effect on the PFS hazard",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Latent $MIXTURE class index (MIXNUM == 1). 1 = subpopulation in",
        "which the change-from-baseline tumor-size coefficient is",
        "multiplied by (1 + 13.6); estimated population probability 0.479",
        "(Table 2). For population simulation draw",
        "MIX_HIGH_TS_HAZARD ~ Bernoulli(0.479) per measurable-disease",
        "patient; it has no effect when RECIST_NONMEASURABLE = 1."
      ),
      source_name = "MIXNUM"
    ),
    LMET = list(
      description = "Baseline liver (hepatic) metastases",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies the baseline PFS hazard by (1 + 1.23).",
      source_name = "HEPM_BL (2 = present)"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Multiplies the baseline PFS hazard by (1 + 0.479). MONARCH 2",
        "enrolled ECOG 0 or 1 only."
      ),
      source_name = "ECOG_BL"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "abemaciclib", units = "umol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "abemaciclib", units = "umol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "abemaciclib", units = "umol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "abemaciclib", units = "umol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "abemaciclib", units = "umol", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "abemaciclib", units = "umol", specimen = "tissue", verified = TRUE),
    central = list(analyte = "abemaciclib", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "abemaciclib", units = "umol", specimen = "plasma", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor", units = "mm", specimen = "not applicable", verified = TRUE),
    damaged_cells1 = list(analyte = "damaged tumor", units = "mm", specimen = "not applicable", verified = TRUE),
    damaged_cells2 = list(analyte = "damaged tumor", units = "mm", specimen = "not applicable", verified = TRUE),
    damaged_cells3 = list(analyte = "damaged tumor", units = "mm", specimen = "not applicable", verified = TRUE),
    cumhaz = list(analyte = "progression", units = "(unitless)", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 663,
    n_studies = 1,
    sex_female_pct = 100,
    disease_state = paste(
      "HR+/HER2- advanced or metastatic breast cancer that progressed on",
      "endocrine therapy; 477 with measurable disease (tumor-size data) and",
      "186 with non-measurable bone-only disease."
    ),
    dose_range = paste(
      "Abemaciclib 200 mg (pre-amendment) or 150 mg (post-amendment) orally",
      "every 12 h with 50 mg dose reductions (average 160 and 136 mg), plus",
      "fulvestrant; control arm placebo plus fulvestrant."
    ),
    regions = "Multinational (MONARCH 2, NCT02107703).",
    notes = paste(
      "440 patients randomized to abemaciclib + fulvestrant and 223 to",
      "placebo + fulvestrant (Results). The paper does not tabulate age,",
      "weight or race for this analysis set. The PK layer was developed in",
      "486 healthy subjects and cancer patients and re-estimated with",
      "MONARCH 1 and MONARCH 2 data (Chigutsa 2020 Table S2)."
    )
  )

  ini({
    # ---- PK layer: Chigutsa 2020 Table S2 (model re-estimated with MONARCH 1 + 2) ----
    logitfa <- logit(0.92); label("Fraction of dose entering the gut, Fa (logit scale)") # Chigutsa 2020 Table S2: Fa THETA1 = 0.92
    logitfdepot <- logit(0.493); label("Fraction of Fa absorbed via the transit route (logit scale)") # Chigutsa 2020 Table S2: THETA12 = 0.493
    ltmax_abs <- log(94.1); label("Maximum saturable absorption rate from the last transit compartment, ka,max (umol/h)") # Chigutsa 2020 Table S2: THETA2 = 94.1
    lkt_abs <- log(53.6); label("Amount giving half-maximal absorption rate, Km (umol)") # Chigutsa 2020 Table S2: THETA18 = 53.6
    lmtt <- log(1.54); label("Mean transit time of the transit route (h)") # Chigutsa 2020 Table S2: THETA3 = 1.54
    lvc <- log(549); label("Central volume for 70 kg (L)") # Chigutsa 2020 Table S2: THETA4 = 549
    lq <- log(4.37); label("Intercompartmental clearance for 70 kg (L/h)") # Chigutsa 2020 Table S2: THETA5 = 4.37
    lvp <- log(174); label("Peripheral volume for 70 kg (L)") # Chigutsa 2020 Table S2: THETA6 = 174
    ld1 <- log(2.84); label("Duration of the zero-order absorption route into the liver (h)") # Chigutsa 2020 Table S2: THETA7 = 2.84
    ltlag <- log(1.71); label("Lag time of the zero-order absorption route (h)") # Chigutsa 2020 Table S2: THETA9 = 1.71
    lkgutex_m2 <- log(0.347); label("Rate constant of gut formation of M2 from the last transit compartment (1/h)") # Chigutsa 2020 Table S2: THETA8 = 0.347
    lkgutex_m20 <- log(0.259); label("Rate constant of gut formation of M20 from the last transit compartment (1/h)") # Chigutsa 2020 Table S2: THETA15 = 0.259
    lclint_m2 <- log(285); label("Unbound hepatic intrinsic clearance of abemaciclib forming M2, 70 kg (L/h)") # Chigutsa 2020 Table S2: THETA10 = 285
    lclint_m20 <- log(819); label("Unbound hepatic intrinsic clearance of abemaciclib forming M20, 70 kg (L/h)") # Chigutsa 2020 Table S2: THETA11 = 819
    boxcox_ld1 <- 0.292; label("Box-Cox shape of the zero-order-duration random effect (unitless)") # Chigutsa 2020 Table S2: THETA16 = 0.292
    e_ae_diarrhea_fa <- -1.28; label("Additive effect of diarrhea on logit(Fa) (unitless)") # Chigutsa 2020 Table S2: THETA17 = -1.28
    e_form_c2_tlag <- 0.822; label("Fractional increase in lag time for the C2 capsule (unitless)") # Chigutsa 2020 Table S2: THETA19 = 0.822
    e_form_c3_tlag <- 1.06; label("Fractional increase in lag time for the C3 capsule (unitless)") # Chigutsa 2020 Table S2: THETA20 = 1.06
    e_form_c3_tmax_abs <- -0.16; label("Fractional change in ka,max for the C3 capsule (unitless)") # Chigutsa 2020 Table S2: THETA21 = -0.16
    fu <- fixed(0.027); label("Fraction unbound of abemaciclib in plasma (unitless)") # Chigutsa 2020 control stream: FU = 0.027 (measured, ADME study)
    bp <- fixed(0.84); label("Blood-to-plasma concentration ratio of abemaciclib (unitless)") # Chigutsa 2020 / 2026 control streams: EH uses FU/0.84

    # PK IIV: omega^2 = log(1 + CV^2) from Chigutsa 2020 Table S2 %CV.
    etaltmax_abs ~ 0.316342 # Table S2: ka,max BSV 61.0%
    etald1 + etalmtt ~ c(2.78995, 0.381398, 0.46459) # Table S2: D BSV 390.9%, MTT BSV 76.9%, correlation 0.335
    etalvc ~ 0.401597 # Table S2: central volume BSV 70.3%
    etalclint ~ 0.180674 # Table S2: total CLint BSV in patients 44.5%
    etalogitfa ~ fixed(0.005) # Table S2: Fa 7.1% (OMEGA 0.005, not estimated)
    etalogitfdepot ~ fixed(0.005) # Table S2: transit fraction 7.1% (OMEGA 0.005, not estimated)
    etalkt_abs ~ fixed(0.005) # Table S2: Km 7.1% (OMEGA 0.005, not estimated)
    etalq ~ fixed(0.005) # Table S2: Q 7.1% (OMEGA 0.005, not estimated)
    etalvp ~ fixed(0.005) # Table S2: peripheral volume 7.1% (OMEGA 0.005, not estimated)
    etaltlag ~ fixed(0.005) # Table S2: lag time 7.1% (OMEGA 0.005, not estimated)
    etalkgutex_m2 ~ fixed(0.005) # Table S2: gut M2 formation 7.1% (OMEGA 0.005, not estimated)
    etalkgutex_m20 ~ fixed(0.005) # Table S2: gut M20 formation 7.1% (OMEGA 0.005, not estimated)
    etalclint_m2 ~ fixed(0.005) # Chigutsa 2020 control stream: ETA5 on CLint forming M2, OMEGA 0.005, not estimated
    etalclint_m20 ~ fixed(0.005) # Chigutsa 2020 control stream: ETA8 on CLint forming M20, OMEGA 0.005, not estimated

    expSd <- 0.170; label("Residual error of log abemaciclib concentration (log-scale SD)") # Chigutsa 2020 Table S2: proportional error parent (oral) 17.0%

    # ---- Tumor-size layer: Chigutsa 2026 Table 1 ----
    lrbase <- log(36.7); label("Baseline tumor size, sum of longest diameters (mm)") # Table 1: 36.7 mm
    lktr <- log(0.934); label("Damaged-cell transit rate constant (1/week)") # Table 1: 0.934 /week
    lkgl <- log(0.00715); label("Zero-order tumor growth rate (mm/week)") # Table 1: 0.00715 mm/week
    lkkill <- log(0.00647); label("Tumor damage rate constant due to fulvestrant (1/week)") # Table 1: 0.00647 /week
    lkkill_slope <- log(0.0399); label("Abemaciclib tumor damage rate per unit plasma concentration (1/week per umol/L)") # Table 1: 0.0399 /week
    lthalf_res <- log(2.51); label("Half-life of the waning tumor drug effect, control arm (week)") # Table 1: 2.51 weeks
    lthalf_res_abemaciclib <- log(8.44); label("Half-life of the waning tumor drug effect, abemaciclib arm (week)") # Table 1: 8.44 weeks
    ltdelay_res <- log(7.11); label("Delay before the tumor drug effect starts to wane (week)") # Table 1: 7.11 weeks

    # Tumor IIV: omega^2 = log(1 + CV^2) (Table 1 footnote a); correlations from Table 1.
    etalrbase + etalkgl ~ c(0.484954, 0.567973, 9.47247) # Table 1: BSV 79.0% and 11400%, correlation 0.265
    etalktr ~ fixed(0.0225) # Table 1: 15.1% (OMEGA 0.0225, not estimated; footnote b)
    etalkkill + etalkkill_slope + etalthalf_res + etalthalf_res_abemaciclib + etaltdelay_res ~ c(
      2.10199,
      -0.305151, 2.26019,
      -0.107502, -0.186573, 0.0380749,
      0.394668, -1.69080, 0.0586036, 1.63333,
      1.06042, -1.35314, -0.0298600, 1.43390, 3.15158
    ) # Table 1: BSV 268, 293, 19.7, 203, 473 %CV; 10 correlations from Table 1

    addSd_tumor_size <- 1.39; label("Additive residual error on tumor size (mm)") # Table 1: 1.39 mm
    propSd_tumor_size <- 0.0904; label("Proportional residual error on tumor size (fraction)") # Table 1: 9.04%

    # ---- PFS hazard layer: Chigutsa 2026 Table 2 ----
    llam_haz <- log(0.0157); label("Baseline PFS hazard, measurable disease (1/week)") # Table 2: 0.0157
    llam_haz_nonmeasurable <- log(0.0103); label("Baseline PFS hazard, non-measurable disease (1/week)") # Table 2: 0.0103
    e_cfb_haz <- 0.0156; label("Coefficient of change-from-baseline tumor size on the log hazard (1/mm)") # Table 2: 0.0156
    e_mix_high_ts_hazard_cfb <- 13.6; label("Fractional increase in the tumor-size coefficient in the mixture subpopulation (unitless)") # Table 2: 13.6
    slope_haz <- 1.42; label("Coefficient of the non-tumor-size drug effect on the log hazard (L/umol)") # Table 2: 1.42
    lthalf_haz <- log(12.2); label("Half-life of the waning drug effect on the hazard, measurable disease (week)") # Table 2: 12.2 weeks
    e_lmet_haz <- 1.23; label("Fractional increase in baseline hazard with hepatic metastases (unitless)") # Table 2: 1.23
    e_ecog_ge1_haz <- 0.479; label("Fractional increase in baseline hazard for ECOG 1 (unitless)") # Table 2: 0.479
  })

  model({
    # Units: model time is weeks; Chigutsa 2020 PK rates are per hour and
    # are multiplied by 168 h/week (Chigutsa 2026 control stream).
    hweek <- 168
    mw <- 0.5066 # abemaciclib molecular weight (mg/umol); dose mg -> umol

    # Body-surface area (DuBois, exponents as printed in both control
    # streams) and liver volume (Urata 1995), L.
    bsa <- WT^0.4245 * HT^0.725 * 0.007184
    vliv <- (706.2 * bsa + 2.4) / 1000
    # Hepatic blood flow 1.35 L/min at 70 kg, allometric exponent 0.75.
    qh <- 1.35 * 60 * hweek * (WT / 70)^0.75

    # ---- PK individual parameters ----
    fa <- expit(logitfa + etalogitfa + e_ae_diarrhea_fa * AE_DIARRHEA)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    tmax_abs <- exp(ltmax_abs + etaltmax_abs) * (1 + e_form_c3_tmax_abs * FORM_ABEMACICLIB_C3) * hweek
    kt_abs <- exp(lkt_abs + etalkt_abs)
    mtt <- exp(lmtt + etalmtt) / hweek
    ktr_abs <- (4 + 1) / mtt
    vc <- exp(lvc + etalvc) * (WT / 70)
    vp <- exp(lvp + etalvp) * (WT / 70)
    q <- exp(lq + etalq) * (WT / 70)^0.75 * hweek
    bc_d1 <- (exp(etald1)^boxcox_ld1 - 1) / boxcox_ld1
    d1 <- exp(ld1 + bc_d1) / hweek
    tlag <- exp(ltlag + etaltlag) * (1 + e_form_c2_tlag * FORM_ABEMACICLIB_C2) *
      (1 + e_form_c3_tlag * FORM_ABEMACICLIB_C3) / hweek
    kgutex_m2 <- exp(lkgutex_m2 + etalkgutex_m2)
    kgutex_m20 <- exp(lkgutex_m20 + etalkgutex_m20)
    kgutex <- (kgutex_m2 + kgutex_m20) * hweek

    # Well-stirred liver: total intrinsic clearance (sum of the M2 and M20
    # formation pathways) with patient BSV, unbound in blood via fu / bp.
    clint_m2 <- exp(lclint_m2 + etalclint_m2)
    clint_m20 <- exp(lclint_m20 + etalclint_m20)
    clint <- (clint_m2 + clint_m20) * exp(etalclint) * (WT / 70)^0.75 * hweek
    eh <- clint * (fu / bp) / (qh + clint * fu / bp)
    k_liver_central <- qh * (1 - eh) / vliv
    k_liver_elim <- qh * eh / vliv
    k_central_liver <- qh / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PK ODEs (amounts in umol) ----
    abs_sat <- tmax_abs * transit4 / (kt_abs + transit4)
    d/dt(depot) <- -ktr_abs * depot
    d/dt(transit1) <- ktr_abs * depot - ktr_abs * transit1
    d/dt(transit2) <- ktr_abs * transit1 - ktr_abs * transit2
    d/dt(transit3) <- ktr_abs * transit2 - ktr_abs * transit3
    d/dt(transit4) <- ktr_abs * transit3 - abs_sat - kgutex * transit4
    d/dt(liver) <- abs_sat - k_liver_central * liver + k_central_liver * central - k_liver_elim * liver
    d/dt(central) <- k_liver_central * liver - k_central_liver * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Each oral dose is entered twice (same amt): once into depot (transit
    # route) and once into liver with rate = -2 (lagged zero-order route).
    f(depot) <- fa * fdepot / mw
    f(liver) <- fa * (1 - fdepot) / mw
    alag(liver) <- tlag
    dur(liver) <- d1

    cmolar <- central / vc # umol/L, the concentration that drives both effects
    if (cmolar < 0) cmolar <- 0

    # ---- Tumor-size layer ----
    rbase <- exp(lrbase + etalrbase)
    ktr <- exp(lktr + etalktr)
    kgl <- exp(lkgl + etalkgl)
    kkill <- exp(lkkill + etalkkill)
    kkill_slope <- exp(lkkill_slope + etalkkill_slope)
    thalf_res <- exp(lthalf_res + etalthalf_res)
    thalf_res_abemaciclib <- exp(lthalf_res_abemaciclib + etalthalf_res_abemaciclib)
    tdelay_res <- exp(ltdelay_res + etaltdelay_res)
    thalf_arm <- thalf_res * (1 - TRT_ABEMACICLIB) + thalf_res_abemaciclib * TRT_ABEMACICLIB
    kres <- log(2) / thalf_arm
    t_res <- t - tdelay_res
    if (t_res < 0) t_res <- 0
    kdamage <- (kkill + kkill_slope * cmolar) * exp(-kres * t_res)

    cycling_cells(0) <- rbase
    d/dt(cycling_cells) <- kgl - kdamage * cycling_cells
    d/dt(damaged_cells1) <- kdamage * cycling_cells - ktr * damaged_cells1
    d/dt(damaged_cells2) <- ktr * damaged_cells1 - ktr * damaged_cells2
    d/dt(damaged_cells3) <- ktr * damaged_cells2 - ktr * damaged_cells3
    tumor_size <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    # ---- PFS hazard layer ----
    measurable <- 1 - RECIST_NONMEASURABLE
    hbase <- (exp(llam_haz) * measurable + exp(llam_haz_nonmeasurable) * RECIST_NONMEASURABLE) *
      (1 + e_lmet_haz * LMET) * (1 + e_ecog_ge1_haz * ECOG_GE1)
    thaz <- e_cfb_haz * (1 + e_mix_high_ts_hazard_cfb * MIX_HIGH_TS_HAZARD) * measurable
    cfb <- (tumor_size - rbase) * measurable
    kdecay_haz <- log(2) / exp(lthalf_haz) * measurable
    hazard <- hbase * exp(thaz * cfb - slope_haz * cmolar * exp(-kdecay_haz * t))
    # Numerical guard: stop integrating once survival < exp(-50). The
    # 11400% CV on the growth rate gives a few subjects tumor growth of tens
    # of mm/week, whose hazard then overflows and fails the solve.
    hazard_ode <- hazard
    if (cumhaz > 50) hazard_ode <- 0
    d/dt(cumhaz) <- hazard_ode
    surv <- exp(-cumhaz)

    Cc <- cmolar * mw * 1000 # ng/mL
    Cc ~ lnorm(expSd)
    tumor_size ~ add(addSd_tumor_size) + prop(propSd_tumor_size)
  })
}
