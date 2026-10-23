Kudryashova_2024_nsclc_survival_mbma <- function() {
  description <- paste(
    "MBMA. Integrative overall-survival model of a mixed non-small cell lung",
    "cancer (NSCLC) cohort across all TNM stages, built from twelve",
    "Royston-Parmar flexible parametric survival sub-models (restricted cubic",
    "spline on the log cumulative hazard in log time). Three sub-models",
    "describe early-stage disease (stages I, II and IIIa after surgery with or",
    "without adjuvant therapy, from IASLC 8th-edition staging epidemiology);",
    "nine describe first-line treatment arms for stage IIIb/IV disease from",
    "five randomised trials (ALEX, FLAURA, KEYNOTE-407, KEYNOTE-189,",
    "KEYNOTE-024). Cohort survival is a prevalence-weighted mixture of the",
    "sub-model survival curves (Eq. 7), evaluated simultaneously for the",
    "'conventional' and 'optimal' stage IIIb/IV treatment schemes. Two ODE",
    "accumulators per scheme integrate the cohort survival curve (life-years,",
    "Eqs. 8-10) and its stage-utility-weighted counterpart (quality-adjusted",
    "life-years, Eq. 11); life-years gained and QALY gained are differences of",
    "these areas between two scenarios at the same time horizon. Stage",
    "distribution defaults to the paper's 'Standard' scenario and is",
    "overridable through the frac_stage* parameters. Deterministic: no",
    "random effects, no residual error; the paper's confidence bands come",
    "from parameter-uncertainty sampling, not from between-subject",
    "variability. Simulates cohort-level survival curves, not individual",
    "event times. Time is in months."
  )
  reference <- paste(
    "Kudryashova N, Shulgin B, Katuninks N, Kulesh V, Helmlinger G,",
    "Zhudenkov K, Peskov K. Assessment of NSCLC disease burden: A survival",
    "model-based meta-analysis study. Comput Struct Biotechnol J.",
    "2024;24:611-621. doi:10.1016/j.csbj.2024.09.012.",
    "Equations: main text Eqs. 2-11. Spline coefficients and knots:",
    "Supplementary Table S2a. Stage IIIb/IV sub-group prevalence: Figure 2.",
    "Stage distributions and utility scores: Table 1."
  )
  vignette <- "Kudryashova_2024_nsclc_survival_mbma"

  # Cohort-level survival model with no drug input: nothing is dosed and the
  # outputs are dimensionless survival probabilities and areas under survival
  # curves in months. The dosing / concentration strings are n/a in the
  # Zierhut_2016_hcc_antiangiogenic_os_mbma style.
  units <- list(
    time = "month",
    dosing = "n/a (no dosing; treatment enters as the choice of sub-model)",
    concentration = "n/a (no drug concentration is modelled)",
    response = "survival probability (unitless); auc_* states in months"
  )

  covariateData <- list()

  # The four integrator states are bookkeeping accumulators for the paper's
  # pharmacoeconomic areas (Eqs. 8-11). They are not drug pools, so they sit
  # outside the compartment register.
  paper_specific_compartments <- c(
    "auc_sur_conv",
    "auc_sur_opt",
    "auc_qaly_conv",
    "auc_qaly_opt"
  )

  compartmentData <- list(
    auc_sur_conv = list(
      analyte = "area under the cohort overall-survival curve, conventional stage IIIb/IV treatment",
      units = "month",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_sur_opt = list(
      analyte = "area under the cohort overall-survival curve, optimal stage IIIb/IV treatment",
      units = "month",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_qaly_conv = list(
      analyte = "area under the stage-utility-weighted cohort survival curve, conventional stage IIIb/IV treatment",
      units = "month",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_qaly_opt = list(
      analyte = "area under the stage-utility-weighted cohort survival curve, optimal stage IIIb/IV treatment",
      units = "month",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15808L,
    n_studies = 6L,
    age_range = "Not reported at the cohort level; summary-level data from one staging-epidemiology source and five randomised trials.",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Described as representative of a Caucasian population; trials enrolling only Asian patients were excluded (Kudryashova 2024 Section 3.1).",
    disease_state = paste(
      "Non-small cell lung cancer, all TNM stages.",
      "Early stages (I, II, IIIa): 13469 patients treated with surgery with or without adjuvant therapy,",
      "digitised from the IASLC 8th-edition staging project survival curves (Goldstraw 2016):",
      "stage I 8231 (IA1 781, IA2 3105, IA3 2417, IB 1928), stage II 2038 (IIA 585, IIB 1453), stage IIIa 3200.",
      "Advanced stages (IIIb/IV): 2339 patients across ten first-line Kaplan-Meier curves:",
      "ALEX (alectinib 152, crizotinib 151; ALK-positive), FLAURA (osimertinib 279, gefitinib/erlotinib 277; EGFR-mutated),",
      "KEYNOTE-407 (pembrolizumab + chemotherapy 278, chemotherapy 281; squamous),",
      "KEYNOTE-189 (pembrolizumab + chemotherapy 410, chemotherapy 206; non-squamous),",
      "KEYNOTE-024 (pembrolizumab 154, chemotherapy 151; PD-L1 tumour proportion score of 50 or more percent)",
      "(Kudryashova 2024 Supplementary Table S6).",
      sep = " "
    ),
    dose_range = "Not modelled. Each treatment arm is a separate survival sub-model; no dose or exposure enters the model.",
    regions = "Multinational trials; Asian-only, elderly-only, poor-performance-status and heavy-smoker studies were excluded.",
    notes = paste(
      "Summary-level model-based meta-analysis: individual patient data were reconstructed from digitised",
      "Kaplan-Meier curves with the Guyot algorithm, then one Royston-Parmar model was fitted per curve",
      "(the two ALEX arms pooled into one alectinib/crizotinib sub-model), with 0, 1 or 2 internal knots chosen by AIC",
      "and an RSE below 50 percent (Supplementary Table S1). n_studies counts the staging-epidemiology source plus the five trials.",
      sep = " "
    )
  )

  ini({
    # Royston-Parmar coefficients (rpg*) and knot positions on the log-time
    # scale (rpk*), Supplementary Table S2a. The knot list in Table S2a
    # includes both boundary knots: every sub-model lists (internal knots + 2)
    # knots and (internal knots + 2) gamma coefficients. The boundary knots of
    # the zero-internal-knot (Weibull) sub-models only bound the fitted data
    # window and do not enter the survival function, so they are not carried
    # as parameters; they are recorded in the trailing comments.

    # Model 1: stage I, surgery with or without adjuvant therapy (2 internal knots)
    rpg0_stage1 <- fixed(-6.19); label("Royston-Parmar intercept gamma0, stage I (log cumulative hazard)") # Table S2a model 1 gamma0 = -6.19 (RSE 3.2)
    rpg1_stage1 <- fixed(1.26); label("Royston-Parmar log-time coefficient gamma1, stage I (unitless)") # Table S2a model 1 gamma1 = 1.26 (RSE 7.5)
    rpg2_stage1 <- fixed(0.24); label("Royston-Parmar first spline-basis coefficient gamma2, stage I (unitless)") # Table S2a model 1 gamma2 = 0.24 (RSE 22.2)
    rpg3_stage1 <- fixed(-0.62); label("Royston-Parmar second spline-basis coefficient gamma3, stage I (unitless)") # Table S2a model 1 gamma3 = -0.62 (RSE 20.3)
    rpkmin_stage1 <- fixed(-0.43); label("Lower boundary knot, stage I (log month)") # Table S2a model 1 knots -0.43, 3.24, 3.90, 4.26
    rpk1_stage1 <- fixed(3.24); label("First internal knot, stage I (log month)") # Table S2a model 1
    rpk2_stage1 <- fixed(3.90); label("Second internal knot, stage I (log month)") # Table S2a model 1
    rpkmax_stage1 <- fixed(4.26); label("Upper boundary knot, stage I (log month)") # Table S2a model 1

    # Model 2: stage II, surgery with or without adjuvant therapy (2 internal knots)
    rpg0_stage2 <- fixed(-4.89); label("Royston-Parmar intercept gamma0, stage II (log cumulative hazard)") # Table S2a model 2 gamma0 = -4.89 (RSE 4.1)
    rpg1_stage2 <- fixed(1.24); label("Royston-Parmar log-time coefficient gamma1, stage II (unitless)") # Table S2a model 2 gamma1 = 1.24 (RSE 11.9)
    rpg2_stage2 <- fixed(-0.09); label("Royston-Parmar first spline-basis coefficient gamma2, stage II (unitless)") # Table S2a model 2 gamma2 = -0.09 (RSE 42.1)
    rpg3_stage2 <- fixed(0.21); label("Royston-Parmar second spline-basis coefficient gamma3, stage II (unitless)") # Table S2a model 2 gamma3 = 0.21 (RSE 28.7)
    rpkmin_stage2 <- fixed(-1.50); label("Lower boundary knot, stage II (log month)") # Table S2a model 2 knots -1.50, 2.69, 3.52, 4.26
    rpk1_stage2 <- fixed(2.69); label("First internal knot, stage II (log month)") # Table S2a model 2
    rpk2_stage2 <- fixed(3.52); label("Second internal knot, stage II (log month)") # Table S2a model 2
    rpkmax_stage2 <- fixed(4.26); label("Upper boundary knot, stage II (log month)") # Table S2a model 2

    # Model 3: stage IIIa, surgery with or without adjuvant therapy (2 internal knots)
    rpg0_stage3a <- fixed(-3.94); label("Royston-Parmar intercept gamma0, stage IIIa (log cumulative hazard)") # Table S2a model 3 gamma0 = -3.94 (RSE 2.5)
    rpg1_stage3a <- fixed(1.11); label("Royston-Parmar log-time coefficient gamma1, stage IIIa (unitless)") # Table S2a model 3 gamma1 = 1.11 (RSE 8.0)
    rpg2_stage3a <- fixed(-0.13); label("Royston-Parmar first spline-basis coefficient gamma2, stage IIIa (unitless)") # Table S2a model 3 gamma2 = -0.13 (RSE 14.9)
    rpg3_stage3a <- fixed(0.23); label("Royston-Parmar second spline-basis coefficient gamma3, stage IIIa (unitless)") # Table S2a model 3 gamma3 = 0.23 (RSE 11.1)
    rpkmin_stage3a <- fixed(-1.98); label("Lower boundary knot, stage IIIa (log month)") # Table S2a model 3 knots -1.98, 2.29, 3.17, 4.23
    rpk1_stage3a <- fixed(2.29); label("First internal knot, stage IIIa (log month)") # Table S2a model 3
    rpk2_stage3a <- fixed(3.17); label("Second internal knot, stage IIIa (log month)") # Table S2a model 3
    rpkmax_stage3a <- fixed(4.23); label("Upper boundary knot, stage IIIa (log month)") # Table S2a model 3

    # Model 4: stage IIIb/IV squamous, pembrolizumab + chemotherapy (Weibull)
    rpg0_sq_immuno <- fixed(-3.66); label("Royston-Parmar intercept gamma0, squamous pembrolizumab + chemotherapy (log cumulative hazard)") # Table S2a model 4 gamma0 = -3.66 (RSE 6.7); boundary knots -1.15, 3.28
    rpg1_sq_immuno <- fixed(1.16); label("Royston-Parmar log-time coefficient gamma1, squamous pembrolizumab + chemotherapy (unitless)") # Table S2a model 4 gamma1 = 1.16 (RSE 6.5)

    # Model 5: stage IIIb/IV squamous, chemotherapy (Weibull)
    rpg0_sq_chemo <- fixed(-3.16); label("Royston-Parmar intercept gamma0, squamous chemotherapy (log cumulative hazard)") # Table S2a model 5 gamma0 = -3.16 (RSE 6.3); boundary knots -0.95, 3.38
    rpg1_sq_chemo <- fixed(1.10); label("Royston-Parmar log-time coefficient gamma1, squamous chemotherapy (unitless)") # Table S2a model 5 gamma1 = 1.10 (RSE 5.6)

    # Model 6: stage IIIb/IV EGFR-mutated, gefitinib/erlotinib (Weibull)
    rpg0_egfr_conv <- fixed(-4.83); label("Royston-Parmar intercept gamma0, EGFR-mutated gefitinib/erlotinib (log cumulative hazard)") # Table S2a model 6 gamma0 = -4.83 (RSE 6.5); boundary knots -0.80, 3.84
    rpg1_egfr_conv <- fixed(1.29); label("Royston-Parmar log-time coefficient gamma1, EGFR-mutated gefitinib/erlotinib (unitless)") # Table S2a model 6 gamma1 = 1.29 (RSE 6.5)

    # Model 7: stage IIIb/IV EGFR-mutated, osimertinib (Weibull)
    rpg0_egfr_osi <- fixed(-6.39); label("Royston-Parmar intercept gamma0, EGFR-mutated osimertinib (log cumulative hazard)") # Table S2a model 7 gamma0 = -6.39 (RSE 6.6); boundary knots 0.53, 3.87
    rpg1_egfr_osi <- fixed(1.66); label("Royston-Parmar log-time coefficient gamma1, EGFR-mutated osimertinib (unitless)") # Table S2a model 7 gamma1 = 1.66 (RSE 6.8)

    # Model 8: stage IIIb/IV ALK-positive, crizotinib/alectinib (2 internal knots)
    rpg0_alk <- fixed(-4.24); label("Royston-Parmar intercept gamma0, ALK-positive crizotinib/alectinib (log cumulative hazard)") # Table S2a model 8 gamma0 = -4.24 (RSE 8.2)
    rpg1_alk <- fixed(0.71); label("Royston-Parmar log-time coefficient gamma1, ALK-positive crizotinib/alectinib (unitless)") # Table S2a model 8 gamma1 = 0.71 (RSE 42.0)
    rpg2_alk <- fixed(-0.23); label("Royston-Parmar first spline-basis coefficient gamma2, ALK-positive crizotinib/alectinib (unitless)") # Table S2a model 8 gamma2 = -0.23 (RSE 40.1)
    rpg3_alk <- fixed(0.39); label("Royston-Parmar second spline-basis coefficient gamma3, ALK-positive crizotinib/alectinib (unitless)") # Table S2a model 8 gamma3 = 0.39 (RSE 0.4)
    rpkmin_alk <- fixed(-2.74); label("Lower boundary knot, ALK-positive crizotinib/alectinib (log month)") # Table S2a model 8 knots -2.74, 2.34, 3.08, 3.90
    rpk1_alk <- fixed(2.34); label("First internal knot, ALK-positive crizotinib/alectinib (log month)") # Table S2a model 8
    rpk2_alk <- fixed(3.08); label("Second internal knot, ALK-positive crizotinib/alectinib (log month)") # Table S2a model 8
    rpkmax_alk <- fixed(3.90); label("Upper boundary knot, ALK-positive crizotinib/alectinib (log month)") # Table S2a model 8

    # Model 9: stage IIIb/IV non-squamous, pembrolizumab + chemotherapy (Weibull)
    rpg0_nsq_immuno <- fixed(-3.89); label("Royston-Parmar intercept gamma0, non-squamous pembrolizumab + chemotherapy (log cumulative hazard)") # Table S2a model 9 gamma0 = -3.89 (RSE 5.7); boundary knots -0.50, 3.29
    rpg1_nsq_immuno <- fixed(1.16); label("Royston-Parmar log-time coefficient gamma1, non-squamous pembrolizumab + chemotherapy (unitless)") # Table S2a model 9 gamma1 = 1.16 (RSE 5.7)

    # Model 10: stage IIIb/IV non-squamous, chemotherapy (1 internal knot)
    rpg0_nsq_chemo <- fixed(-3.43); label("Royston-Parmar intercept gamma0, non-squamous chemotherapy (log cumulative hazard)") # Table S2a model 10 gamma0 = -3.43 (RSE 9.2)
    rpg1_nsq_chemo <- fixed(1.63); label("Royston-Parmar log-time coefficient gamma1, non-squamous chemotherapy (unitless)") # Table S2a model 10 gamma1 = 1.63 (RSE 15.9)
    rpg2_nsq_chemo <- fixed(0.06); label("Royston-Parmar first spline-basis coefficient gamma2, non-squamous chemotherapy (unitless)") # Table S2a model 10 gamma2 = 0.06 (RSE 38.3)
    rpkmin_nsq_chemo <- fixed(-1.21); label("Lower boundary knot, non-squamous chemotherapy (log month)") # Table S2a model 10 knots -1.21, 2.01, 3.36
    rpk1_nsq_chemo <- fixed(2.01); label("Internal knot, non-squamous chemotherapy (log month)") # Table S2a model 10
    rpkmax_nsq_chemo <- fixed(3.36); label("Upper boundary knot, non-squamous chemotherapy (log month)") # Table S2a model 10

    # Model 11: stage IIIb/IV PD-L1 high, pembrolizumab monotherapy (Weibull)
    rpg0_pdl1_immuno <- fixed(-3.40); label("Royston-Parmar intercept gamma0, PD-L1-high pembrolizumab monotherapy (log cumulative hazard)") # Table S2a model 11 gamma0 = -3.40 (RSE 9.8); boundary knots -3.08, 3.17
    rpg1_pdl1_immuno <- fixed(0.94); label("Royston-Parmar log-time coefficient gamma1, PD-L1-high pembrolizumab monotherapy (unitless)") # Table S2a model 11 gamma1 = 0.94 (RSE 10.6)

    # Model 12: stage IIIb/IV PD-L1 high, chemotherapy (1 internal knot)
    rpg0_pdl1_chemo <- fixed(-3.29); label("Royston-Parmar intercept gamma0, PD-L1-high chemotherapy (log cumulative hazard)") # Table S2a model 12 gamma0 = -3.29 (RSE 10.5)
    rpg1_pdl1_chemo <- fixed(1.37); label("Royston-Parmar log-time coefficient gamma1, PD-L1-high chemotherapy (unitless)") # Table S2a model 12 gamma1 = 1.37 (RSE 19.3)
    rpg2_pdl1_chemo <- fixed(0.06); label("Royston-Parmar first spline-basis coefficient gamma2, PD-L1-high chemotherapy (unitless)") # Table S2a model 12 gamma2 = 0.06 (RSE 0.0 as printed)
    rpkmin_pdl1_chemo <- fixed(-0.84); label("Lower boundary knot, PD-L1-high chemotherapy (log month)") # Table S2a model 12 knots -0.84, 2.04, 3.21
    rpk1_pdl1_chemo <- fixed(2.04); label("Internal knot, PD-L1-high chemotherapy (log month)") # Table S2a model 12
    rpkmax_pdl1_chemo <- fixed(3.21); label("Upper boundary knot, PD-L1-high chemotherapy (log month)") # Table S2a model 12

    # Prevalence of each sub-group within stage IIIb/IV (Figure 2 ovals; sum to 1)
    frac_alk <- fixed(0.05); label("Proportion of stage IIIb/IV patients who are ALK-positive (fraction)") # Figure 2 oval '5%'
    frac_egfr <- fixed(0.15); label("Proportion of stage IIIb/IV patients who are EGFR-mutated (fraction)") # Figure 2 oval '15%'
    frac_sq <- fixed(0.25); label("Proportion of stage IIIb/IV patients with squamous histology and no driver mutation (fraction)") # Figure 2 oval '25%'
    frac_nsq <- fixed(0.45); label("Proportion of stage IIIb/IV patients with non-squamous histology and no driver mutation (fraction)") # Figure 2 oval '45%'
    frac_pdl1 <- fixed(0.10); label("Proportion of stage IIIb/IV patients with PD-L1 tumour proportion score above 50 (fraction)") # Figure 2 oval '10%'

    # Stage distribution at diagnosis; defaults are the 'Standard' scenario.
    # Scenario 1: 0.29 / 0.08 / 0.13 / 0.50. Scenario 2: 0.35 / 0.15 / 0.10 / 0.40.
    frac_stage1 <- fixed(0.15); label("Proportion of the cohort diagnosed at stage I (fraction)") # Table 1 'Standard' stage I = 15
    frac_stage2 <- fixed(0.14); label("Proportion of the cohort diagnosed at stage II (fraction)") # Table 1 'Standard' stage II = 14
    frac_stage3a <- fixed(0.14); label("Proportion of the cohort diagnosed at stage IIIa (fraction)") # Table 1 'Standard' stage IIIa = 14
    frac_stage3b4 <- fixed(0.57); label("Proportion of the cohort diagnosed at stage IIIb/IV (fraction)") # Table 1 'Standard' stage IIIb/IV = 57

    # Health-state utility scores by stage (Table 1 'Utility Score' row, from the paper's ref. 4)
    util_stage1 <- fixed(0.87); label("Health-state utility score, stage I (unitless)") # Table 1 utility stage I = 0.87
    util_stage2 <- fixed(0.87); label("Health-state utility score, stage II (unitless)") # Table 1 utility stage II = 0.87
    util_stage3a <- fixed(0.77); label("Health-state utility score, stage IIIa (unitless)") # Table 1 utility stage IIIa = 0.77
    util_stage3b4 <- fixed(0.57); label("Health-state utility score, stage IIIb/IV (unitless)") # Table 1 utility stage IIIb/IV = 0.57
  })
  model({
    # Log time in months. The floor keeps log() finite at t = 0, where every
    # sub-model survival is 1; at 1e-6 months the largest cumulative hazard
    # across the twelve sub-models is below 1e-5.
    lt <- log(max(t, 1e-6))

    # Restricted cubic spline basis (Eqs. 4-6):
    #   v_j = (x - k_j)^3_+ - lam_j (x - kmin)^3_+ - (1 - lam_j) (x - kmax)^3_+
    #   lam_j = (kmax - k_j) / (kmax - kmin)
    # The (x - k)^3_+ truncated power is written (x - k)^3 * (x > k).
    # log H(t) = gamma0 + gamma1 x + sum_j gamma_{j+1} v_j (Eqs. 2-3), and
    # S(t) = exp(-H(t)).

    # Stage I
    rplam1_stage1 <- (rpkmax_stage1 - rpk1_stage1) / (rpkmax_stage1 - rpkmin_stage1)
    rplam2_stage1 <- (rpkmax_stage1 - rpk2_stage1) / (rpkmax_stage1 - rpkmin_stage1)
    rpcmin_stage1 <- (lt - rpkmin_stage1)^3 * (lt > rpkmin_stage1)
    rpcmax_stage1 <- (lt - rpkmax_stage1)^3 * (lt > rpkmax_stage1)
    rpv1_stage1 <- (lt - rpk1_stage1)^3 * (lt > rpk1_stage1) - rplam1_stage1 * rpcmin_stage1 - (1 - rplam1_stage1) * rpcmax_stage1
    rpv2_stage1 <- (lt - rpk2_stage1)^3 * (lt > rpk2_stage1) - rplam2_stage1 * rpcmin_stage1 - (1 - rplam2_stage1) * rpcmax_stage1
    lcumhaz_stage1 <- rpg0_stage1 + rpg1_stage1 * lt + rpg2_stage1 * rpv1_stage1 + rpg3_stage1 * rpv2_stage1
    sur_stage1 <- exp(-exp(lcumhaz_stage1))

    # Stage II
    rplam1_stage2 <- (rpkmax_stage2 - rpk1_stage2) / (rpkmax_stage2 - rpkmin_stage2)
    rplam2_stage2 <- (rpkmax_stage2 - rpk2_stage2) / (rpkmax_stage2 - rpkmin_stage2)
    rpcmin_stage2 <- (lt - rpkmin_stage2)^3 * (lt > rpkmin_stage2)
    rpcmax_stage2 <- (lt - rpkmax_stage2)^3 * (lt > rpkmax_stage2)
    rpv1_stage2 <- (lt - rpk1_stage2)^3 * (lt > rpk1_stage2) - rplam1_stage2 * rpcmin_stage2 - (1 - rplam1_stage2) * rpcmax_stage2
    rpv2_stage2 <- (lt - rpk2_stage2)^3 * (lt > rpk2_stage2) - rplam2_stage2 * rpcmin_stage2 - (1 - rplam2_stage2) * rpcmax_stage2
    lcumhaz_stage2 <- rpg0_stage2 + rpg1_stage2 * lt + rpg2_stage2 * rpv1_stage2 + rpg3_stage2 * rpv2_stage2
    sur_stage2 <- exp(-exp(lcumhaz_stage2))

    # Stage IIIa
    rplam1_stage3a <- (rpkmax_stage3a - rpk1_stage3a) / (rpkmax_stage3a - rpkmin_stage3a)
    rplam2_stage3a <- (rpkmax_stage3a - rpk2_stage3a) / (rpkmax_stage3a - rpkmin_stage3a)
    rpcmin_stage3a <- (lt - rpkmin_stage3a)^3 * (lt > rpkmin_stage3a)
    rpcmax_stage3a <- (lt - rpkmax_stage3a)^3 * (lt > rpkmax_stage3a)
    rpv1_stage3a <- (lt - rpk1_stage3a)^3 * (lt > rpk1_stage3a) - rplam1_stage3a * rpcmin_stage3a - (1 - rplam1_stage3a) * rpcmax_stage3a
    rpv2_stage3a <- (lt - rpk2_stage3a)^3 * (lt > rpk2_stage3a) - rplam2_stage3a * rpcmin_stage3a - (1 - rplam2_stage3a) * rpcmax_stage3a
    lcumhaz_stage3a <- rpg0_stage3a + rpg1_stage3a * lt + rpg2_stage3a * rpv1_stage3a + rpg3_stage3a * rpv2_stage3a
    sur_stage3a <- exp(-exp(lcumhaz_stage3a))

    # ALK-positive, crizotinib/alectinib (same arm in both treatment schemes)
    rplam1_alk <- (rpkmax_alk - rpk1_alk) / (rpkmax_alk - rpkmin_alk)
    rplam2_alk <- (rpkmax_alk - rpk2_alk) / (rpkmax_alk - rpkmin_alk)
    rpcmin_alk <- (lt - rpkmin_alk)^3 * (lt > rpkmin_alk)
    rpcmax_alk <- (lt - rpkmax_alk)^3 * (lt > rpkmax_alk)
    rpv1_alk <- (lt - rpk1_alk)^3 * (lt > rpk1_alk) - rplam1_alk * rpcmin_alk - (1 - rplam1_alk) * rpcmax_alk
    rpv2_alk <- (lt - rpk2_alk)^3 * (lt > rpk2_alk) - rplam2_alk * rpcmin_alk - (1 - rplam2_alk) * rpcmax_alk
    lcumhaz_alk <- rpg0_alk + rpg1_alk * lt + rpg2_alk * rpv1_alk + rpg3_alk * rpv2_alk
    sur_alk <- exp(-exp(lcumhaz_alk))

    # Non-squamous chemotherapy (1 internal knot)
    rplam1_nsq_chemo <- (rpkmax_nsq_chemo - rpk1_nsq_chemo) / (rpkmax_nsq_chemo - rpkmin_nsq_chemo)
    rpv1_nsq_chemo <- (lt - rpk1_nsq_chemo)^3 * (lt > rpk1_nsq_chemo) -
      rplam1_nsq_chemo * (lt - rpkmin_nsq_chemo)^3 * (lt > rpkmin_nsq_chemo) -
      (1 - rplam1_nsq_chemo) * (lt - rpkmax_nsq_chemo)^3 * (lt > rpkmax_nsq_chemo)
    lcumhaz_nsq_chemo <- rpg0_nsq_chemo + rpg1_nsq_chemo * lt + rpg2_nsq_chemo * rpv1_nsq_chemo
    sur_nsq_chemo <- exp(-exp(lcumhaz_nsq_chemo))

    # PD-L1-high chemotherapy (1 internal knot)
    rplam1_pdl1_chemo <- (rpkmax_pdl1_chemo - rpk1_pdl1_chemo) / (rpkmax_pdl1_chemo - rpkmin_pdl1_chemo)
    rpv1_pdl1_chemo <- (lt - rpk1_pdl1_chemo)^3 * (lt > rpk1_pdl1_chemo) -
      rplam1_pdl1_chemo * (lt - rpkmin_pdl1_chemo)^3 * (lt > rpkmin_pdl1_chemo) -
      (1 - rplam1_pdl1_chemo) * (lt - rpkmax_pdl1_chemo)^3 * (lt > rpkmax_pdl1_chemo)
    lcumhaz_pdl1_chemo <- rpg0_pdl1_chemo + rpg1_pdl1_chemo * lt + rpg2_pdl1_chemo * rpv1_pdl1_chemo
    sur_pdl1_chemo <- exp(-exp(lcumhaz_pdl1_chemo))

    # Zero internal knots: Weibull, log H(t) = gamma0 + gamma1 log(t)
    sur_sq_immuno <- exp(-exp(rpg0_sq_immuno + rpg1_sq_immuno * lt))
    sur_sq_chemo <- exp(-exp(rpg0_sq_chemo + rpg1_sq_chemo * lt))
    sur_egfr_conv <- exp(-exp(rpg0_egfr_conv + rpg1_egfr_conv * lt))
    sur_egfr_osi <- exp(-exp(rpg0_egfr_osi + rpg1_egfr_osi * lt))
    sur_nsq_immuno <- exp(-exp(rpg0_nsq_immuno + rpg1_nsq_immuno * lt))
    sur_pdl1_immuno <- exp(-exp(rpg0_pdl1_immuno + rpg1_pdl1_immuno * lt))

    # Stage IIIb/IV mixtures (Eq. 7). Conventional: alectinib/crizotinib,
    # gefitinib/erlotinib and chemotherapy. Optimal: alectinib/crizotinib,
    # osimertinib, pembrolizumab + chemotherapy (squamous, non-squamous) and
    # pembrolizumab monotherapy (PD-L1 high). Section 3.2.
    sur_stage3b4_conv <- frac_alk * sur_alk + frac_egfr * sur_egfr_conv +
      frac_sq * sur_sq_chemo + frac_nsq * sur_nsq_chemo + frac_pdl1 * sur_pdl1_chemo
    sur_stage3b4_opt <- frac_alk * sur_alk + frac_egfr * sur_egfr_osi +
      frac_sq * sur_sq_immuno + frac_nsq * sur_nsq_immuno + frac_pdl1 * sur_pdl1_immuno

    # Whole-cohort overall survival (Eq. 7 over stages)
    sur_early <- frac_stage1 * sur_stage1 + frac_stage2 * sur_stage2 + frac_stage3a * sur_stage3a
    sur_conv <- sur_early + frac_stage3b4 * sur_stage3b4_conv
    sur_opt <- sur_early + frac_stage3b4 * sur_stage3b4_opt

    # Stage-utility-weighted survival, the integrand of Eq. 11
    surq_early <- util_stage1 * frac_stage1 * sur_stage1 +
      util_stage2 * frac_stage2 * sur_stage2 +
      util_stage3a * frac_stage3a * sur_stage3a
    surq_conv <- surq_early + util_stage3b4 * frac_stage3b4 * sur_stage3b4_conv
    surq_opt <- surq_early + util_stage3b4 * frac_stage3b4 * sur_stage3b4_opt

    # Areas under the curves from 0 to t (Eqs. 8-9 and the integrals of Eq. 11)
    d/dt(auc_sur_conv) <- sur_conv
    d/dt(auc_sur_opt) <- sur_opt
    d/dt(auc_qaly_conv) <- surq_conv
    d/dt(auc_qaly_opt) <- surq_opt
  })
}
