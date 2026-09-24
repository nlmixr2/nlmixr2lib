Bi_2018_testosteroneCypionate <- function() {
  description <- "Coupled hypothalamic-pituitary-gonadal (HPG) axis PK/PD model for weekly intramuscular depot testosterone cypionate in 31 healthy men. One-compartment first-order-absorption PK for total testosterone whose endogenous secretion rate is up-regulated by luteinizing hormone through a sigmoid Emax term, linked to an effect-compartment-driven indirect response model for LHRH-stimulated luteinizing hormone measured 30 minutes after stimulation (LH30). Clearance and volume carry Wahlby-style baseline and change-from-baseline body-weight and serum-albumin covariates; LH30 potency carries baseline weight and LH30 synthesis carries baseline thyroxine."
  reference <- paste(
    "Bi Y, Perry PJ, Ellerby M, Murry DJ.",
    "Population Pharmacokinetic/Pharmacodynamic Modeling of Depot Testosterone Cypionate in Healthy Male Subjects.",
    "CPT Pharmacometrics Syst Pharmacol 2018;7(4):259-268.",
    "doi:10.1002/psp4.12287.",
    "Structural equations and covariate reference/centring constants taken from the",
    "NONMEM control streams in Supplementary Material S6; final parameter values from Tables 2 and 3.",
    sep = " "
  )
  vignette <- "Bi_2018_testosteroneCypionate"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # LHRH-stimulated luteinizing hormone measured 30 minutes after the stimulation
  # test is the paper's own PD endpoint and has no canonical compartment name.
  paper_specific_compartments <- c("lh30")

  # Volumes are in kL and amounts in mg, so central / vc is mg/kL = ug/L = ng/mL,
  # which is the unit total testosterone is reported in throughout the paper.
  compartmentData <- list(
    depot = list(analyte = "testosterone cypionate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "testosterone (total)", units = "mg", specimen = "serum", verified = TRUE),
    effect = list(analyte = "testosterone (total)", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    lh30 = list(analyte = "luteinizing hormone", units = "IU/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT_BASE = list(
      description = "Per-subject baseline body weight, time-fixed.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Bi 2018 follows Wahlby 2004 (Supplementary Text, reference 7) in splitting body weight into a time-fixed baseline column and a change-from-baseline column. Enters CL and V as a power function referenced to 85 kg (Supplementary S6 PK stream: CLWT = (WT/85)**THETA(11), VCWT = (WT/85)**THETA(12)) and LH30 potency as a power function referenced to 84.70 kg (S6 LH stream: IC50WT = ((WT/84.70)**THETA(7))). Table 1 group medians 82.2 / 88.8 / 84.7 kg, pooled range 60.7-115 kg.",
      source_name = "WT"
    ),
    WT = list(
      description = "Body weight, time-varying within a subject.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only through the within-subject deviation (WT - WT_BASE), which is the source paper's CWT column. Enters CL exponentially, centred at +2.90 kg (Supplementary S6 PK stream: CLCWT = EXP(THETA(13)*(CWT - 2.90))). The centring constant is the cohort median weight gain during dosing and appears nowhere in the main article.",
      source_name = "CWT (as WT - WT_BASE)"
    ),
    ALB_BASE = list(
      description = "Per-subject baseline serum albumin, time-fixed.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Canonical SI units are g/L; Bi 2018 reports and calibrates in US-convention g/dL, so model() applies the register-mandated inline conversion ALB_BASE * 0.1 before the power term. Reference 4.60 g/dL = 46.0 g/L (Supplementary S6 PK stream: VALBUMIN = (ALBUMIN/4.60)**THETA(15)). Table 1 group medians 4.55 / 4.55 / 4.6 g/dL.",
      source_name = "ALBUMIN"
    ),
    ALB = list(
      description = "Serum albumin, time-varying within a subject.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only through the within-subject deviation (ALB - ALB_BASE), the source paper's CALBUMIN column, converted to g/dL inside model(). Enters V linearly, centred at -0.2 g/dL (Supplementary S6 PK stream: VCALBUMIN = (1 + THETA(14)*(CALBUMIN + 0.2))). Bi 2018 Results: 9 of 496 albumin covariate records (1.81%) were missing and carried forward from the previous measurement.",
      source_name = "CALBUMIN (as ALB - ALB_BASE)"
    ),
    T4 = list(
      description = "Baseline serum total thyroxine (T4).",
      units = "ug/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LH30 zero-order synthesis as a power function referenced to 7.40 ug/dL (Supplementary S6 LH stream: KINTHYROXIN = ((THYROXIN/7.40)**THETA(8))). Bi 2018 Table 1 does NOT tabulate this analyte - the table lists thyroxine-BINDING GLOBULIN (group medians 19.5 / 18.5 / 21 ug/mL), which is a different analyte on a different scale. The LH dataset shipped as Supplementary S8 carries thyroxin = 7.4 for its example subject, matching the reference constant and the conventional total-T4 adult range of roughly 4.5-12 ug/dL rather than the TBG range; the baseline distribution of this covariate is therefore not published. Time-fixed per subject.",
      source_name = "THYROXIN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "21-39 years",
    age_median = "26.5 / 25.5 / 30 years in the 100 / 250 / 500 mg groups",
    weight_range = "60.7-115 kg",
    weight_median = "82.2 / 88.8 / 84.7 kg in the 100 / 250 / 500 mg groups",
    sex_female_pct = 0,
    disease_state = "healthy men",
    dose_range = "100, 250 or 500 mg testosterone cypionate intramuscularly once weekly for 14 consecutive weeks (study weeks 2-15), preceded by 2 and followed by 12 weekly placebo injections",
    n_observations = 729L,
    follow_up = "40 weeks",
    regions = "United States",
    notes = "Randomised, double-blind trial originally reported by MacIndoe JH et al. J Investig Med 1997;45:441-447 (Bi 2018 reference 11); baseline demographics and laboratory values in Bi 2018 Table 1. 729 total-testosterone and 379 luteinizing-hormone serum samples were available; 299 testosterone samples fell within 26 days of the last dose and carried the disposition information. The 250 and 500 mg/week arms deliberately mimic supratherapeutic doses used illicitly rather than replacement-therapy doses."
  )

  ini({
    # --- Testosterone disposition (Bi 2018 Table 2, Final model estimate column) ---
    lka <- log(1.22); label("First-order absorption rate constant from the intramuscular depot (1/day)") # Table 2, Ka 1.22 /day
    lcl <- log(2.6); label("Apparent clearance of total testosterone (CL/F, kL/day)") # Table 2, CL/F 2.6 kL/day
    lvc <- log(14.4); label("Apparent central volume of distribution (V/F, kL)") # Table 2, V/F 14.4 kL
    lrbase_te <- log(6.24); label("Baseline total testosterone concentration setting the central-compartment initial condition (ng/mL)") # Table 2, 'A_0, TC' 6.24; S6 PK stream INT=THETA(10)*EXP(ETA(6)), A_0(2)=INT*V, so the estimate is a concentration

    # --- Endogenous testosterone secretion, up-regulated by LH (Bi 2018 Table 2) ---
    lkin_te <- log(6.24); label("Basal LH-independent endogenous testosterone secretion rate (mg/day)") # Table 2, 'b base' 6.24 mg/day
    lemax_te <- log(12.5); label("Maximum LH-driven increment to endogenous testosterone secretion (mg/day)") # Table 2, Emax 12.5 mg/day
    lec50_te <- log(12.9); label("LH30 concentration giving half-maximal up-regulation of testosterone secretion (IU/L)") # Table 2, LH50 12.9 IU/L
    # Held on the natural scale, not log-transformed, because Supplementary S6
    # estimates it as an UNBOUNDED $THETA (1.17938, no parentheses) and gives it a
    # proportional rather than log-normal random effect.
    hill_te <- 1.18; label("Hill coefficient for LH up-regulation of testosterone secretion (unitless)") # Table 2, 'Lambda (k)' 1.18

    # --- Testosterone PK covariate effects (Bi 2018 Table 2; forms and centring from Supplementary S6) ---
    e_bwt_cl <- 0.785; label("Power exponent of (WT_BASE / 85) on CL (unitless)") # Table 2, 'h BWT-CL' 0.785
    e_dwt_cl <- 0.016; label("Coefficient on (WT - WT_BASE - 2.90) in the exponential CL term (1/kg)") # Table 2, 'h CWT-CL' 0.016
    e_bwt_vc <- 1.71; label("Power exponent of (WT_BASE / 85) on V (unitless)") # Table 2, 'h BWT-V' 1.71
    e_balb_vc <- -1.55; label("Power exponent of (ALB_BASE / 4.60 g/dL) on V (unitless)") # Table 2, 'h Balbumin-V' -1.55
    e_dalb_vc <- -0.27; label("Coefficient on (ALB - ALB_BASE + 0.2 g/dL) in the linear V term (dL/g)") # Table 2, 'h Calbumin-V' -0.27

    # --- LH30 indirect response (Bi 2018 Table 3, LH block) ---
    lkin_lh <- log(1.7); label("Zero-order LH30 synthesis rate (IU/L/day)") # Table 3, Kin 1.7
    lkout_lh <- log(0.11); label("First-order LH30 loss rate constant (1/day)") # Table 3, Kout 0.11
    limax_lh <- fixed(log(1)); label("Maximum fractional inhibition of LH30 synthesis by testosterone (unitless)") # Bi 2018 Results: 'The Emax is fixed to 1 as LH synthesis is fully inhibited during the study'
    lic50_lh <- log(9.33); label("Effect-compartment testosterone concentration inhibiting LH30 synthesis by half (ng/mL)") # Table 3, 'Test50' 9.33 ng/mL (called TC50 in the Results text and IC50 in Figure 1 and the S6 LH stream)
    lhill_lh <- log(18.3); label("Hill coefficient for testosterone inhibition of LH30 synthesis (unitless)") # Table 3, 'Lambda; k' 18.3
    lke0 <- log(0.0321); label("Effect-compartment equilibration rate constant (1/day)") # Table 3, ke0 0.0321

    # --- LH30 covariate effects (Bi 2018 Table 3; forms and references from Supplementary S6) ---
    e_bwt_ic50_lh <- -1.14; label("Power exponent of (WT_BASE / 84.70) on the LH30 inhibitory potency (unitless)") # Table 3, 'h Bwt-IC50' -1.14
    e_t4_kin_lh <- 1.19; label("Power exponent of (T4 / 7.40 ug/dL) on LH30 synthesis (unitless)") # Table 3, 'h Bthyroxin-Kin' 1.19

    # --- Inter-individual variability ---
    # Variances are the Supplementary S6 $OMEGA entries, which reproduce the
    # percentages in Tables 2 and 3 as 100 * sqrt(omega).
    etalcl ~ 0.00932598 # S6 PK stream $OMEGA 1; Table 2 'IIV_CL' 9.66%
    etalvc ~ 0.137296 # S6 PK stream $OMEGA 2; Table 2 'IIV_V' 37%
    etalka ~ 0.289713 # S6 PK stream $OMEGA 3; Table 2 'IIV_Ka' 53.9%
    etalemax_te ~ 0.205164 # S6 PK stream $OMEGA 4; Table 2 'IIV_Emax' 45.3%
    etahill_te ~ 0.575982 # S6 PK stream $OMEGA 5; Table 2 'IIV_k' 75.9%. Enters PROPORTIONALLY, not log-normally: S6 codes LAM=TVLAM*(1+ETA(5))
    etalrbase_te ~ 0.018456 # S6 PK stream $OMEGA 6; Table 2 'IIV_A_0, TC' 13.6%
    etalkin_lh ~ 0.153716 # S6 LH stream $OMEGA 1; Table 3 'IIV_Kin' 39.2%
    etalhill_lh ~ 0.177732 # S6 LH stream $OMEGA 2; Table 3 'IIV_k' 42.2%
    etalke0 ~ 0.217495 # S6 LH stream $OMEGA 3; Table 3 'IIV_ke0' 46.6%
    etalic50_lh ~ 0.0767912 # S6 LH stream $OMEGA 4; Table 3 'IIV_IC50' 27.7%

    # --- Residual error ---
    # The S6 PK stream builds W = SQRT(IPRED*IPRED*THETA(8) + THETA(9)) with $SIGMA 1 FIX,
    # so Table 2's two 'sigma^2' rows really are variances and the SDs are their square roots.
    propSd <- 0.239; label("Proportional residual error for total testosterone (fraction)") # sqrt(0.0572167) from S6 PK $THETA 8; Table 2 'sigma^2 proportional' 0.057
    addSd <- 0.511; label("Additive residual error for total testosterone (ng/mL)") # sqrt(0.261394) from S6 PK $THETA 9; Table 2 'sigma^2 additive' 0.261
    # The S6 LH stream instead sets W = THETA(6) directly on log(LH30 + 1), so Table 3's
    # 'sigma^2 additive' row is a standard deviation despite its label.
    addSd_logLH30 <- 0.397; label("Additive residual error on log(LH30 + 1) (log IU/L)") # S6 LH $THETA 6 magnitude 0.396587; Table 3 'sigma^2 additive' 0.397
  })

  model({
    # --- Covariate transforms -------------------------------------------------
    # The register holds albumin in SI g/L; Bi 2018 calibrated in US-convention g/dL.
    alb_base_gdL <- ALB_BASE * 0.1
    dalb_gdL <- (ALB - ALB_BASE) * 0.1
    dwt <- WT - WT_BASE

    # --- Individual disposition parameters ------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT_BASE / 85)^e_bwt_cl *
      exp(e_dwt_cl * (dwt - 2.90))
    vc <- exp(lvc + etalvc) *
      (WT_BASE / 85)^e_bwt_vc *
      (alb_base_gdL / 4.60)^e_balb_vc *
      (1 + e_dalb_vc * (dalb_gdL + 0.2))
    kel <- cl / vc

    # --- Endogenous secretion parameters --------------------------------------
    kin_te <- exp(lkin_te)
    emax_te <- exp(lemax_te + etalemax_te)
    ec50_te <- exp(lec50_te)
    # Proportional (not log-normal) random effect, replicating S6's LAM = TVLAM * (1 + ETA(5)).
    # Contrast the LH and sperm layers, where the same Hill coefficient carries a
    # log-normal effect (LAM = TVLAM * EXP(ETA)).
    hill_tei <- hill_te * (1 + etahill_te)
    rbase_te <- exp(lrbase_te + etalrbase_te)

    # --- LH30 parameters ------------------------------------------------------
    kin_lh <- exp(lkin_lh + etalkin_lh) * (T4 / 7.40)^e_t4_kin_lh
    kout_lh <- exp(lkout_lh)
    imax_lh <- exp(limax_lh)
    ic50_lh <- exp(lic50_lh + etalic50_lh) * (WT_BASE / 84.70)^e_bwt_ic50_lh
    hill_lh <- exp(lhill_lh + etalhill_lh)
    ke0 <- exp(lke0 + etalke0)

    # --- Initial conditions ---------------------------------------------------
    # S6 PK stream: A_0(2) = INT * V. S6 LH stream: A_0(4) = kin / kout, and the
    # effect compartment starts in equilibrium with the baseline concentration.
    central(0) <- rbase_te * vc
    effect(0) <- rbase_te
    lh30(0) <- kin_lh / kout_lh

    # --- Testosterone PK with LH-driven endogenous secretion ------------------
    # The floor reproduces S6's IF (LHORI .EQ. 0) k1 = BASE branch and keeps the
    # power term finite for the subjects whose proportional eta drives hill_te
    # below zero.
    lh_drive <- max(lh30, 1e-8)
    ksec_te <- kin_te + emax_te * lh_drive^hill_tei / (lh_drive^hill_tei + ec50_te^hill_tei)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ksec_te + ka * depot - kel * central

    Cc <- central / vc

    # --- Effect compartment and LH30 indirect response ------------------------
    d/dt(effect) <- ke0 * (Cc - effect)
    # S6 LH stream clamps a negative effect-site concentration to zero before the
    # power term: IF (Ce .LE. 0) Ce = 0.
    ce <- max(effect, 0)
    inh_lh <- imax_lh * ce^hill_lh / (ce^hill_lh + ic50_lh^hill_lh)
    d/dt(lh30) <- kin_lh * (1 - inh_lh) - kout_lh * lh30

    # --- Observations ---------------------------------------------------------
    # Bi 2018 Supplementary Text: 'natural log-transformation was applied to the
    # data before the PD analysis'; S6 sets IPRED = LOG(A(4) + 1), so the +1
    # offset is part of the observation model and matters while LH30 is near zero.
    logLH30 <- log(lh30 + 1)

    Cc ~ add(addSd) + prop(propSd)
    logLH30 ~ add(addSd_logLH30)
  })
}
