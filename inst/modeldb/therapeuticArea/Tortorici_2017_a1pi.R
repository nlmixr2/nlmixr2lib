Tortorici_2017_a1pi <- function() {
  description <- "Empirical disease-progression model of intravenous alpha-1 proteinase inhibitor (A1-PI) augmentation therapy in alpha-1 antitrypsin deficiency (Tortorici 2017). Combines a sequential dose-exposure regression that predicts the per-subject median trough serum A1-PI from average dose rate (mg/day, supplied as DOSE covariate), baseline body weight, and baseline endogenous A1-PI, with a piecewise-linear-in-time exposure-response model that predicts CT lung density at total lung capacity from the dose-derived A1-PI exposure and baseline FEV1, with a slope transition at 720 days separating the RAPID-RCT and RAPID-OLE study phases."

  reference <- paste(
    "Tortorici MA, Rogers JA, Vit O, Bexon M, Sandhaus RA, Burdon J,",
    "Chorostowska-Wynimko J, Thompson P, Stocks J, McElvaney NG, Chapman KR,",
    "Edelman JM. Quantitative disease progression model of alpha-1 proteinase",
    "inhibitor therapy on computed tomography lung density in patients with",
    "alpha-1 antitrypsin deficiency. Br J Clin Pharmacol 2017;83(11):2386-2397.",
    "doi:10.1111/bcp.13358.",
    sep = " "
  )
  vignette <- "Tortorici_2017_a1pi"
  units <- list(
    time          = "day",
    dosing        = "mg/day (per-subject average dose rate, supplied as the DOSE covariate column; the model is an empirical regression and does not consume rxode2 dose events)",
    concentration = "umol/L (serum A1-PI, observation Cc); g/L (CT lung density at total lung capacity, observation lungDens)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject time-fixed average A1-PI dose rate during the modelled study phase (mg/day; not weight-normalised). Computed from the prescribed weekly weight-based regimen as DOSE = (mg_per_kg_per_week * WT) / 7. The reference 60 mg/kg/wk regimen for a 77 kg subject corresponds to DOSE ~ 660 mg/day; the placebo arm has DOSE = 0.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tortorici 2017 dose-exposure model (DOSE in equation 6 is the average daily mg dose, denoted Dij). Time-fixed within a study phase; the user can switch DOSE between simulation phases by stratifying the dataset, but the algebraic dose-exposure relationship has no PK lag, so transitions appear as instantaneous step changes in Cc. Use case (a) of the canonical DOSE entry: per-subject assigned dose level used as a regressor in a PD / dose-exposure model that does not instantiate a PK ODE.",
      source_name        = "Dij"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tortorici 2017 used baseline weight; mean (SD) reported in Table 1 was 75.9 (16.2) kg in the A1-PI arm and 79.5 (13.9) kg in placebo. Power-form effect on the dose-exposure slope: (WT/77)^theta3 with reference 77 kg (the median-weight reference individual the paper uses for covariate predictions).",
      source_name        = "WT"
    ),
    A1PI = list(
      description        = "Baseline (pre-treatment) endogenous serum alpha-1 proteinase inhibitor concentration. Each subject's screening A1-PI value, used as the Cbase covariate in equation 6.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Two power-form effects: (A1PI/5.5)^theta5 with theta5 = 0.73 on the placebo-arm post-treatment exposure intercept (treats subjects with higher endogenous A1-PI as having proportionally higher placebo-arm levels), and (A1PI/5.5)^theta4 with theta4 = -0.12 on the dose-rate slope. The reference 5.5 umol/L corresponds approximately to the median pre-treatment A1-PI in the RAPID-RCT placebo arm (paper Methods, dose-exposure section). RAPID enrolment criteria restricted A1PI to <= 11 umol/L (severe deficiency).",
      source_name        = "Cbase"
    ),
    FEV1 = list(
      description        = "Baseline forced expiratory volume in 1 second; absolute volume in litres.",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on the lung-density decline rate: theta5 * (FEV1 - 1.6), where 1.6 L is the cohort median (Table 1 shows mean (SD) of 1.6 (0.5) L in both arms of RAPID-RCT). theta5 = +0.56 (g/L/year per L FEV1); patients with lower FEV1 have steeper natural decline rates independent of A1-PI exposure. The covariate effect modifies the natural disease-progression rate, not the exposure response itself.",
      source_name        = "FEV1"
    )
  )

  population <- list(
    n_subjects       = 134L,
    n_studies        = 2L,
    age_range        = "approximately 30-65 years",
    age_median       = "53 years",
    weight_range     = "approximately 50-150 kg",
    weight_median    = "77 kg",
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Adults with severe alpha-1 antitrypsin deficiency and clinical emphysema (RAPID enrolment criteria: PI*ZZ or rare deficient genotype, A1-PI level <= 11 umol/L, FEV1 35-70% predicted, post-bronchodilator FEV1/FVC < 70%).",
    dose_range       = "60 mg/kg weekly intravenous infusion (single dose level studied in RAPID-RCT and RAPID-OLE). Bootstrap dose-exposure simulations in the paper extrapolated to 90 and 120 mg/kg/week.",
    regions          = "International multi-center: Australia, Canada, multiple European countries, USA in RAPID-RCT (n=180); RAPID-OLE continued non-USA enrolment (n=140) because A1-PI augmentation therapy was already approved in the USA.",
    fev1_median      = "1.6 L",
    a1pi_bl_median   = "approximately 6 umol/L (Table 1: 6.38 (4.62) umol/L active arm; 5.94 (2.42) umol/L placebo arm)",
    notes            = "Population is the exposure-response analysis set: 134 subjects (61 placebo and 73 A1-PI active) who had at least one post-baseline CT lung-density measurement. The dose-exposure analysis included a slightly larger 308-subject set (170 from RAPID-RCT plus 138 from the RAPID-OLE continuation, with overlap) but the structural-model parameters reported here come from the exposure-response analysis set's contribution to the joint sequential model. Trial NCT identifiers: NCT00261833 (RAPID-RCT), NCT00670007 (RAPID-OLE). Population characteristics are reproduced from Table 1 of Tortorici 2017 (RAPID-RCT ITT n=180 column)."
  )

  ini({
    # ============================================================
    # Dose-exposure layer (Equation 6, parameter values from Table 2).
    # The published equation is
    #   Cij = ( theta1 * exp(eta1) * (Cbase/5.5)^theta5
    #         + theta2 * (WT/77)^theta3 * (Cbase/5.5)^theta4 * Dij )
    #         * exp(epsilon_ij).
    # theta1 is reported as 5.42 with the units "(umol l-1)" but the column
    # description in Table 2 reads "Log A1-PI exposure for placebo"; reading
    # 5.42 as a back-transformed (linear) value is the only interpretation
    # consistent with the equation, with the placebo-arm baseline in Table 1
    # (~6 umol/L), and with the paper's bootstrap predictions in Figure 2A
    # (~16 umol/L typical post-treatment A1-PI at 60 mg/kg/wk for a 77 kg
    # subject). The omega1 of 0.07 is the log-scale SD for IIV on this
    # back-transformed value; we encode it as a log-normal eta on log(theta1).
    # ============================================================

    # Linear-scale typical placebo-arm post-treatment A1-PI exposure
    lc_pbo    <- log(5.42)         ; label("Typical placebo-arm A1-PI exposure (umol/L)")  # Table 2 theta1 = 5.42
    # Linear-scale slope of A1-PI exposure with respect to average dose rate
    lcslope   <- log(0.02)         ; label("A1-PI exposure slope vs dose rate ((umol/L)/(mg/day))")  # Table 2 theta2 = 0.02
    # Allometric weight exponent on the dose-rate slope; -0.85 is consistent
    # with Kleiber's law (weight-based clearance) per the paper Discussion
    e_wt_cslope    <- -0.85        ; label("WT power exponent on dose-rate slope")  # Table 2 theta3 = -0.85
    # Power exponent on baseline endogenous A1-PI for the dose-rate slope
    e_a1pi_cslope  <- -0.12        ; label("Baseline A1PI power exponent on dose-rate slope")  # Table 2 theta4 = -0.12
    # Power exponent on baseline endogenous A1-PI for the placebo intercept
    e_a1pi_cpbo    <-  0.73        ; label("Baseline A1PI power exponent on placebo intercept")  # Table 2 theta5 = 0.73

    # IIV on placebo-arm typical exposure (log-scale SD = 0.07; variance 0.0049)
    etalc_pbo ~ 0.0049             # Table 2 omega1 = 0.07 (log-scale SD); 0.07^2 = 0.0049

    # Residual: multiplicative on the linear-scale prediction (Eq 6 has
    # exp(epsilon_ij) wrapping the entire dose-exposure expression, so this
    # is log-additive / log-transform-both-sides ~ proportional in nlmixr2).
    propSd <- 0.15              ; label("Proportional residual SD on Cc (log-scale; ~CV)")  # Table 2 sigma = 0.15

    # ============================================================
    # Exposure-response layer (Equations 7-10, parameter values from Table 3).
    # The published structure is
    #   Y_ijk = Int_i + DP_i1 * min(720_days, t_ijk)
    #                  + DP_i2 * max(0, t_ijk - 720_days) + epsilon_ijk
    # with intercept and slopes:
    #   Int_i  = theta1 + eta1_i
    #   DP_i1  = (theta2 + eta2_i) + (theta3 + eta3_i) * Ci_centered
    #            + theta5 * (FEV1_i - FEV1_median)
    #   DP_i2  = DP_i1 + (theta4 + eta4_i)
    # Time t in equations 7-10 is days since RAPID-RCT randomisation;
    # disease-progression slopes are reported in g/L/year, so a /365.25
    # conversion is applied inside model() (see "t_yr" below).
    # ============================================================

    # Pre-treatment lung density intercept (g/L)
    bld          <-  46.89         ; label("Pre-treatment lung density at TLC (g/L)")  # Table 3 theta1 = 46.89
    # Placebo-arm lung-density decline rate (g/L/year); negative = density loss
    dpr_pbo      <-  -2.18         ; label("Placebo-arm lung-density decline rate (g/L/year)")  # Table 3 theta2 = -2.18
    # Slope of decline rate vs A1-PI exposure (centered at typical placebo level)
    e_cc_dpr     <-  0.06          ; label("A1-PI exposure effect on decline rate ((g/L/year)/(umol/L))")  # Table 3 theta3 = 0.06
    # Phase-2 (RAPID-OLE) increment on decline rate (g/L/year)
    dpr_phase2   <-  0.20          ; label("RAPID-OLE phase increment on decline rate (g/L/year)")  # Table 3 theta4 = 0.20
    # Baseline FEV1 effect on decline rate ((g/L/year)/L deviation from median)
    e_fev1_dpr   <-  0.56          ; label("Baseline FEV1 effect on decline rate ((g/L/year)/L)")  # Table 3 theta5 = 0.56

    # 3x3 correlated IIV block on (intercept, decline-rate, exposure-effect):
    # variances and covariances reconstructed from Table 3 SDs and correlations.
    # SDs:   omega1 = 15.28 (Int), omega2 = 1.32 (DPR), omega3 = 0.09 (Cc effect)
    # Corrs: omega12 = -0.270, omega13 = 0.26, omega23 = -0.75
    # Lower-triangle (row-wise):
    #   var1                         = 15.28^2          = 233.4784
    #   cov12 = -0.270*15.28*1.32    = -5.4476
    #   var2                         = 1.32^2           = 1.7424
    #   cov13 = 0.26*15.28*0.09      = 0.3576
    #   cov23 = -0.75*1.32*0.09      = -0.0891
    #   var3                         = 0.09^2           = 0.0081
    etabld + etadpr_pbo + etae_cc_dpr ~ c(
      233.4784,
      -5.4476, 1.7424,
      0.3576, -0.0891, 0.0081
    )                              # Table 3 omega1..3 SDs and corr matrix

    # Phase-2 IIV (poorly estimated per Table 3 footnote; SD = 0.23 with very
    # wide CI 0-6644, retained on conservative-prediction grounds per the
    # paper's Model validity and limitations section).
    etadpr_phase2 ~ 0.0529         # Table 3 omega4 = 0.23 (SD); 0.23^2 = 0.0529

    # Additive residual error on lung density (g/L)
    addSd_lungDens <- 2.60         ; label("Additive residual SD on lung density (g/L)")  # Table 3 sigma = 2.60
  })

  model({
    # ----- Reference values used inside derived quantities -----
    # Reference body weight for the dose-exposure allometric term: median 77 kg
    # (Methods, dose-exposure section: "a reference individual with a weight of
    # 77.0 kg"). Reference baseline A1-PI for the power-form effects: median
    # 5.5 umol/L (Methods, dose-exposure section: "5.5, equates approximately
    # to the median pre-treatment A1-PI concentration among those who were
    # randomized to placebo in the RAPID-RCT study"). Reference FEV1 for the
    # decline-rate covariate effect: median 1.6 L (Table 1 mean = 1.6 L).
    ref_wt_kg    <- 77
    ref_a1pi     <- 5.5
    ref_fev1_l   <- 1.6
    knot_days    <- 720
    days_per_yr  <- 365.25

    # ----- Dose-exposure layer (Equation 6) -----
    # Subject-specific placebo-arm exposure: theta1 * exp(eta1) * (Cbase/5.5)^theta5
    cc_pbo   <- exp(lc_pbo + etalc_pbo) * (A1PI / ref_a1pi)^e_a1pi_cpbo
    # Subject-specific dose-rate slope: theta2 * (WT/77)^theta3 * (Cbase/5.5)^theta4
    cc_slope <- exp(lcslope) * (WT / ref_wt_kg)^e_wt_cslope * (A1PI / ref_a1pi)^e_a1pi_cslope
    # Predicted A1-PI exposure (Cc) given subject-specific dose rate (DOSE)
    Cc <- cc_pbo + cc_slope * DOSE

    # ----- Exposure-response layer (Equations 7-10) -----
    # Convert time to years to match the g/L/year units of the decline-rate
    # parameters in Table 3. The published equation 7 writes t in days but the
    # decline rates DP are in g/L/year, so t/days_per_yr is the consistent
    # dimensional choice.
    t_yr      <- t / days_per_yr
    knot_yr   <- knot_days / days_per_yr  # ~1.971 years

    # Center A1-PI exposure at the typical placebo-arm level (~ theta1 = 5.42).
    # The paper centers Ci* at the median placebo level observed in RAPID-RCT
    # so that theta2 + eta2 represents the decline rate that would be seen at
    # placebo-equivalent exposure ("the natural rate of lung density decline").
    Cc_centered <- Cc - exp(lc_pbo)

    # Subject-specific intercept and slopes (each eta placed on its own
    # mu-referenced line for rxode2 / nlmixr2 mu-reference parsing).
    Int_i      <- bld + etabld
    dpr_pbo_i  <- dpr_pbo + etadpr_pbo
    ecc_dpr_i  <- e_cc_dpr + etae_cc_dpr
    dpr_p2_i   <- dpr_phase2 + etadpr_phase2

    DP1_i <- dpr_pbo_i + ecc_dpr_i * Cc_centered + e_fev1_dpr * (FEV1 - ref_fev1_l)
    DP2_i <- DP1_i + dpr_p2_i

    # Piecewise-linear lung density (g/L) over time
    seg1 <- min(knot_yr, t_yr)
    seg2 <- max(0, t_yr - knot_yr)
    lungDens <- Int_i + DP1_i * seg1 + DP2_i * seg2

    # ----- Residual error -----
    Cc       ~ prop(propSd)
    lungDens ~ add(addSd_lungDens)
  })
}
