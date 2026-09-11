Bender_2024_mosunetuzumab <- function() {
  description <- "Two-compartment population PK model of mosunetuzumab (CD20xCD3 T-cell engaging bispecific antibody) in adults with relapsed/refractory B-cell non-Hodgkin lymphoma, with time-dependent clearance transitioning from a baseline clearance CLbase to a steady-state clearance CLss with a transition half-life HLtrans. Body weight, sex and tumor SPD act on CLss; albumin and the composite baseline anti-CD20 drug concentration act on CLbase; body weight, albumin and sex act on V1. Residual predose rituximab and obinutuzumab from prior therapy are carried as states decaying at fixed literature terminal half-lives and drive a competitive equilibrium-binding CD20 receptor-occupancy percentage (RO%) observable (Bender 2024)."
  reference <- "Bender B, Li C-C, Marchand M, Turner DC, Li F, Vadhavkar S, Wang B, Deng R, Lu J, Jin J, Li C, Yin S, Wei M, Chanu P. Population pharmacokinetics and CD20 binding dynamics for mosunetuzumab in relapsed/refractory B-cell non-Hodgkin lymphoma. Clin Transl Sci. 2024;17(5):e13825. doi:10.1111/cts.13825"
  vignette <- "Bender_2024_mosunetuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # `ritux` and `obin` hold the residual plasma CONCENTRATION (ug/mL) of a
  # competing anti-CD20 antibody left over from the patient's prior lines of
  # therapy. They are not amounts, are never dosed, and exist only to supply
  # the competing-ligand terms of the receptor-occupancy observable. Bender
  # 2024 Table S3 declares them as NONMEM compartments 3 and 5
  # (COMP=(RITUXIMAB), COMP=(OBINUTUZUMAB)) for exactly the same reason.
  #
  # `auc` and `auc_ro` are the paper's own cumulative-endpoint integrators,
  # Table S3 compartments 4 and 6 (COMP=(PK_AUC), COMP=(RO_AUC)); they carry
  # no drug and follow the established auc_<scope> bookkeeping-state pattern.
  paper_specific_compartments <- c("ritux", "obin", "auc", "auc_ro")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power model normalised to the cohort median of 78 kg, on CLss, V1 and V2 only. Bender 2024 Table S3 $PK: CLBWT=(BBWT/78)**THETA(8), VBWT=(BBWT/78)**THETA(9), V2BWT=(BBWT/78)**THETA(11). The weight effect on Q, THETA(10), was fixed to 0 and is therefore absent here. Note the $THETA comment block mislabels THETA(8) as 'WT_CLbase', but the code applies it to TVCLSS; Table 2 and the Results text both name the affected parameter as CLss.",
      source_name        = "BWT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power model normalised to the cohort median of 39 g/L, acting on both CLbase and V1. Bender 2024 Table S3 $PK: CL0ALBUM=(ALBUMT/39.00)**THETA(12) and V1ALBUM=(ALBUMT/39.00)**THETA(16). The control stream guards a unit-error record with IF(ALBUM.GT.200) ALBUMT=39.00, visible in Table 1 as an implausible 480 g/L maximum; that guard is a data-cleaning step and is not reproduced in model().",
      source_name        = "ALBUM"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = "Bender 2024 Table S3 codes SEX = 2 for male (the reference, 64.7% of the cohort) and SEX = 1 for female, and applies the effect as the linear multiplier (1 + THETA), not as a power or exponential term. SEXF = 1 - (SEX == 2) recovers the canonical coding with no change of reference level, so the printed coefficients carry over unchanged.",
      source_name        = "SEX"
    ),
    TUMSZ = list(
      description        = "Baseline tumor burden, sum of the products of perpendicular diameters (SPD)",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters on CLss as a power model in the SQUARE ROOT of SPD, normalised to 54.5 mm: Bender 2024 Table S3 $PK LTS2=SQRT(LTS1); CLSSBSPD=(LTS2/54.5)**THETA(14). The reference is exactly sqrt(2970 mm^2) = 54.5 mm, the cohort median SPD of Table 1, and Figure 2 plots the square root of tumor size on its x-axis for this reason. Units are mm^2 SPD, matching the anti-CD20 sibling model Gibiansky_2014_obinutuzumab.R for the same disease.",
      source_name        = "BSPD"
    ),
    CP_RITUXIMAB_UGML = list(
      description        = "Observed predose (baseline) plasma rituximab concentration remaining from prior therapy",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "BASELINE usage: read once at t = 0 to seed the initial condition of the `ritux` state (Bender 2024 Table S3 $PK: A_0(3) = BLRITUX), which then decays at a fixed 24-day terminal half-life. It does NOT carry the clearance covariate effect -- that acts on the composite CP_ACD20_UGML. 195 of 439 patients had detectable residual rituximab; the reported values are floored at the 0.5 ug/mL assay LOQ (Table 1 median 0.500, maximum 151).",
      source_name        = "BLRITUX"
    ),
    CP_OBINUTUZUMAB_UGML = list(
      description        = "Observed predose (baseline) plasma obinutuzumab concentration remaining from prior therapy",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "BASELINE usage: read once at t = 0 to seed the initial condition of the `obin` state (Bender 2024 Table S3 $PK: A_0(5) = BLOBIN), which then decays at a fixed 28-day terminal half-life. 35 of 439 patients had detectable residual obinutuzumab; Table 1 reports a median of 0 with a maximum of 305 ug/mL. Set to 0 for a patient with no prior obinutuzumab exposure.",
      source_name        = "BLOBIN"
    ),
    CP_ACD20_UGML = list(
      description        = "Composite baseline anti-CD20 drug concentration: the maximum of the predose rituximab and obinutuzumab concentrations",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Bender 2024 Table 1 footnote e: 'aCD20 is the maximum concentration between rituximab and obinutuzumab'; Table S3 derives it as IF(BLOBIN.GT.BLRITUX2) ACD20=BLOBIN / IF(BLRITUX2.GT.BLOBIN) ACD20=BLRITUX2. Supplied as a column rather than computed with max() inside model() because the control stream feeds the covariate an NHL-type-dependent IMPUTED rituximab value (aggressive or unknown NHL -> 2105 ng/mL, indolent NHL -> 500 ng/mL) for the 4.6% of patients with a missing measurement, while seeding the ODE from the raw value. The effect enters as a ratio of LOGARITHMS on the ng/mL scale and is therefore not scale-invariant; see the conversion comment in model().",
      source_name        = "ACD20"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "mosunetuzumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mosunetuzumab", units = "mg", specimen = "plasma", verified = TRUE),
    ritux       = list(analyte = "rituximab", units = "ug/mL", specimen = "plasma", verified = TRUE),
    obin        = list(analyte = "obinutuzumab", units = "ug/mL", specimen = "plasma", verified = TRUE),
    auc         = list(analyte = "mosunetuzumab", units = "ug/mL*day", specimen = "not applicable", verified = TRUE),
    auc_ro      = list(analyte = "mosunetuzumab", units = "%*day", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 439,
    n_studies      = 2,
    age_range      = "19-96 years",
    age_median     = "63 years",
    weight_range   = "37.1-163 kg",
    weight_median  = "77.9 kg",
    sex_female_pct = 35.3,
    race_ethnicity = c(White = 75.9, Asian = 17.5, Black = 2.7,
                       `American Indian/Alaskan Native` = 0.5, Multiple = 0.5,
                       Unknown = 3.0),
    disease_state  = "relapsed/refractory B-cell non-Hodgkin lymphoma (61.5% aggressive, 38.3% indolent; DLBCL 35.8%, FL 37.1%, MCL 8.9%, transformed FL 13.0%)",
    dose_range     = "0.05-2.8 mg IV q3w fixed dosing (Group A, n = 32) and 0.4/1/2.8 up to 1/2/60/30 mg IV q3w Cycle 1 step-up dosing (Group B, n = 407); 19 dose levels; approved regimen 1/2/60/30 mg IV q3w",
    albumin_median = "39 g/L (range 19-480; the 480 g/L maximum is a unit-error record guarded in the control stream)",
    tumor_median   = "SPD 2970 mm^2 (range 96.0-70,900)",
    prior_therapy  = "median 3 prior lines; ~50% of patients carried residual anti-CD20 drug at baseline (rituximab n = 195, obinutuzumab n = 35, both n = 7)",
    notes          = "Study GO29781 (phase I/II), 7250 PK observations from 439 patients after exclusions. Baseline characteristics from Bender 2024 Table 1. Fitted in NONMEM 7.4.3 with ADVAN13/TRANS1 and FOCE-I (Table S3)."
  )

  ini({
    # ---- Structural disposition (Bender 2024 Table 2, final model estimates) ----
    # NOTE: Table S3's $THETA block lists INITIAL estimates (1, 5.4, 0.57, 18,
    # 6.1, 1.47, 0.25, ...). Every value below is the FINAL estimate from
    # Table 2.
    lcl <- log(1.08); label("Baseline clearance CLbase at t = 0 (L/day)") # Bender 2024 Table 2 CLbase = 1.08 L/day (%RSE 5.6; 95% CI 0.962, 1.20)
    lcl_exp_inf <- log(0.584); label("Steady-state (asymptotic) clearance CLss (L/day)") # Bender 2024 Table 2 CLss = 0.584 L/day (%RSE 2.0; 95% CI 0.561, 0.607)
    lcl_exp_thalf <- log(16.3); label("Half-life of the CLbase to CLss transition, HLtrans (day)") # Bender 2024 Table 2 HLtrans = 16.3 day (%RSE 7.1; 95% CI 14.026, 18.6)
    lvc <- log(5.49); label("Central volume of distribution V1 (L)") # Bender 2024 Table 2 V1 = 5.49 L (%RSE 2.5; 95% CI 5.221, 5.76)
    lvp <- log(6.17); label("Peripheral volume of distribution V2 (L)") # Bender 2024 Table 2 V2 = 6.17 L (%RSE 3.6; 95% CI 5.729, 6.61)
    lq <- log(1.46); label("Intercompartmental clearance Q (L/day)") # Bender 2024 Table 2 Q = 1.46 L/day (%RSE 3.7; 95% CI 1.354, 1.57)

    # ---- Covariate effects (Bender 2024 Table 2) ----
    # Body weight enters as a power model on (WT / 78 kg). The weight effect on
    # Q is THETA(10) = 0 FIX in Table S3 and is not reported in Table 2, so no
    # e_wt_q term exists.
    e_wt_cl_exp_inf <- 0.549; label("Power exponent for body weight on CLss (unitless)") # Bender 2024 Table 2 WT_CLss = 0.549 (%RSE 10.5; 95% CI 0.436, 0.662). Table S3 applies THETA(8) to TVCLSS despite its stale 'WT_CLbase' $THETA comment.
    e_wt_vc <- 0.433; label("Power exponent for body weight on V1 (unitless)") # Bender 2024 Table 2 WT_V1 = 0.433 (%RSE 13.4; 95% CI 0.319, 0.547)
    e_wt_vp <- 0.737; label("Power exponent for body weight on V2 (unitless)") # Bender 2024 Table 2 WT_V2 = 0.737 (%RSE 15.9; 95% CI 0.508, 0.966)

    # Albumin enters as a power model on (ALB / 39 g/L).
    e_alb_cl <- -1.51; label("Power exponent for serum albumin on CLbase (unitless)") # Bender 2024 Table 2 ALB_CLbase = -1.51 (%RSE 19.1; 95% CI -2.074, -0.946)
    e_alb_vc <- -0.481; label("Power exponent for serum albumin on V1 (unitless)") # Bender 2024 Table 2 ALB_V1 = -0.481 (%RSE 23.3; 95% CI -0.701, -0.261)

    # Composite baseline anti-CD20 drug concentration on CLbase. See the
    # ratio-of-logarithms comment in model() -- this exponent is NOT applied to
    # a concentration ratio.
    e_acd20_cl <- -0.573; label("Power exponent for the composite baseline anti-CD20 concentration on CLbase (unitless)") # Bender 2024 Table 2 aCD20_CLbase = -0.573 (%RSE 20.2; 95% CI -0.800, -0.346)

    # Tumor SPD enters as a power model on sqrt(SPD) / 54.5 mm.
    e_tumsz_cl_exp_inf <- 0.0935; label("Power exponent for sqrt(tumor SPD) on CLss (unitless)") # Bender 2024 Table 2 SPD_CLss = 0.0935 (%RSE 26.6; 95% CI 0.045, 0.142)

    # Sex enters as the LINEAR multiplier (1 + theta) for female relative to
    # the male reference, not as a power or exponential term (Table S3
    # V1SEX / CLSSSEX definition blocks).
    e_sexf_cl_exp_inf <- -0.128; label("Fractional change in CLss for female vs male (unitless)") # Bender 2024 Table 2 Sex_CLss = -0.128 (%RSE 18.8; 95% CI -0.175, -0.081); Results: '12.8% slower in female subjects'
    e_sexf_vc <- -0.126; label("Fractional change in V1 for female vs male (unitless)") # Bender 2024 Table 2 Sex_V1 = -0.126 (%RSE 18.9; 95% CI -0.173, -0.079); Results: '12.6% lower in female subjects'

    # ---- CD20 equilibrium-binding constants (Bender 2024 Table S3 $PK) ----
    # Scatchard-derived dissociation constants, fixed (they are not $THETAs and
    # carry no uncertainty). Converted from the control stream's ng/mL to the
    # model's declared ug/mL by dividing by 1000; the receptor-occupancy
    # expression is a ratio of concentrations, so the conversion is exact
    # provided every term shares one scale (see model()).
    lkd_mosun <- fixed(log(10.2)); label("Mosunetuzumab CD20 dissociation constant KD (ug/mL)") # Bender 2024 Table S3 KD_TDB = 10200 ng/mL
    lkd_ritux <- fixed(log(0.675)); label("Rituximab CD20 dissociation constant KD (ug/mL)") # Bender 2024 Table S3 KD_R = 675 ng/mL
    lkd_obin <- fixed(log(0.600)); label("Obinutuzumab CD20 dissociation constant KD (ug/mL)") # Bender 2024 Table S3 KD_G = 600 ng/mL

    # ---- Competing-drug elimination half-lives (Bender 2024 Table S3 $DES) ----
    # Fixed to published terminal half-lives, not estimated. Methods: 'Initial
    # values for rituximab (Ritux) and obinutuzumab (Obin) model compartments
    # were set to the observed baseline value, and elimination rates fixed to
    # the respective terminal half-life value: HL_Ritux = 24 days and
    # HL_Obin = 28 days.'
    lthalf_ritux <- fixed(log(24)); label("Rituximab terminal half-life (day)") # Bender 2024 Table S3 DADT(3) = (-0.693/24)*A(3)
    lthalf_obin <- fixed(log(28)); label("Obinutuzumab terminal half-life (day)") # Bender 2024 Table S3 DADT(5) = (-0.693/28)*A(5)

    # ---- Interindividual variability (Bender 2024 Table 2, variances) ----
    # $OMEGA BLOCK(2) on (CLbase, V1). Table 2 labels the covariance row with
    # the symbol 'omega_CLbase,CLss', but the $OMEGA BLOCK comments in Table S3
    # and the reported correlation both identify it as CLbase-V1:
    # 0.180 / sqrt(0.426 * 0.0981) = 0.881, matching footnote c's 0.882.
    etalcl + etalvc ~ c(0.426,
                        0.180, 0.0981) # Bender 2024 Table 2: omega^2 CLbase 0.426 (%RSE 8.4, shrinkage 4.8%), covariance 0.180 (%RSE 7.3), omega^2 V1 0.0981 (%RSE 5.8, shrinkage 4.6%)

    # $OMEGA BLOCK(2) on (CLss, HLtrans); correlation
    # -0.0892 / sqrt(0.0343 * 0.739) = -0.560, matching footnote d.
    etalcl_exp_inf + etalcl_exp_thalf ~ c(0.0343,
                                          -0.0892, 0.739) # Bender 2024 Table 2: omega^2 CLss 0.0343 (%RSE 11.5, shrinkage 33.8%), covariance -0.0892 (%RSE 24.8), omega^2 HLtrans 0.739 (%RSE 15.8, shrinkage 40.9%)

    etalvp ~ 0.0621 # Bender 2024 Table 2 omega^2 V2 = 0.0621 (%RSE 16.7, shrinkage 49.9%)

    # No IIV on Q: Table S3 declares $OMEGA 0 FIX for it and Table 2 reports no
    # omega^2 Q row.

    # ---- Residual unexplained variability ----
    # Table S3 uses the log-transformed-both-sides pattern: IPRED = LOG(F),
    # Y = IPRED + ERR(1)*W with W = THETA(7) and $SIGMA 1 FIX. Additive on the
    # log scale is exponential (log-normal) on the linear scale, which Table 2
    # footnote e states outright ('Corresponds to proportional on normal
    # scale').
    expSd <- 0.259; label("Exponential, log-scale additive, residual error SD (unitless)") # Bender 2024 Table 2 residual variability = 0.259 (%RSE 0.257; 95% CI 0.258, 0.260)
  })

  model({
    # ---- 1. Derived covariate multipliers ----
    # All power models are normalised to the cohort median of the covariate
    # (Bender 2024 Table 1): WT 78 kg, ALB 39 g/L, SPD 2970 mm^2. Sex enters as
    # a linear (1 + theta) multiplier with male as the reference.
    cov_cl_alb <- (ALB / 39)^e_alb_cl
    cov_vc_alb <- (ALB / 39)^e_alb_vc

    # The composite anti-CD20 covariate is a power model on the RATIO OF
    # LOGARITHMS, not on the concentration ratio:
    #   Bender 2024 Table S3: CL0BLRITUX = (LOG(ACD20)/LOG(500))**THETA(13)
    # with ACD20 in ng/mL and a reference of 500 ng/mL (= 0.5 ug/mL, the
    # rituximab LOQ and the value Figure 2 assigns the typical patient). This
    # form is NOT scale-invariant, so the ug/mL column must be converted back
    # to ng/mL INSIDE the logarithm -- rewriting it as log(CP_ACD20_UGML/0.5)
    # gives a different, and at high concentrations undefined, function.
    #
    # The naive concentration-ratio reading (ACD20/500)^-0.573 is falsified by
    # the paper's own text: at the 95th percentile (55.91 ug/mL) it predicts a
    # 93% fall in CLbase, whereas the Results state that 'with the exception of
    # albumin, all covariate effects resulted in <=31% change from the typical
    # parameter values of CLbase, CLss, and V1 when evaluated at the extremes'.
    # The ratio-of-logs form gives 28%, inside that bound. Methods likewise
    # says the anti-CD20 concentrations 'were log-transformed'.
    cov_cl_acd20 <- (log(CP_ACD20_UGML * 1000) / log(500))^e_acd20_cl

    # Tumor burden enters through the SQUARE ROOT of SPD, referenced to
    # sqrt(2970 mm^2) = 54.5 mm (Bender 2024 Table S3 CLSSBSPD block).
    cov_cl_exp_inf_tumsz <- (sqrt(TUMSZ) / 54.5)^e_tumsz_cl_exp_inf

    # ---- 2. Individual parameters ----
    # CLbase carries albumin and aCD20; CLss carries body weight, sex and
    # tumor SPD; V1 carries body weight, albumin and sex; V2 carries body
    # weight only; Q carries no covariate (its weight exponent was fixed to 0).
    clbase <- exp(lcl + etalcl) * cov_cl_alb * cov_cl_acd20
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) * (WT / 78)^e_wt_cl_exp_inf *
      cov_cl_exp_inf_tumsz * (1 + e_sexf_cl_exp_inf * SEXF)
    cl_exp_thalf <- exp(lcl_exp_thalf + etalcl_exp_thalf)
    vc <- exp(lvc + etalvc) * (WT / 78)^e_wt_vc * cov_vc_alb * (1 + e_sexf_vc * SEXF)
    vp <- exp(lvp + etalvp) * (WT / 78)^e_wt_vp
    q <- exp(lq)

    kd_mosun <- exp(lkd_mosun)
    kd_ritux <- exp(lkd_ritux)
    kd_obin <- exp(lkd_obin)
    kel_ritux <- log(2) / exp(lthalf_ritux)
    kel_obin <- log(2) / exp(lthalf_obin)

    # ---- 3. Time-dependent clearance and micro-constants ----
    # Bender 2024 supplement Table S1, Model 2 (the selected model):
    #   CL = CLbase + (CLss - CLbase) * {1 - exp[-(ln(2)/HLtrans) * t]}
    # so CL starts at CLbase, decays exponentially and approaches CLss. Table
    # S3's $DES writes the same expression with 0.693 substituted for ln(2);
    # log(2) is used here because it is the form the published equation
    # prints.
    #
    # `time` is elapsed simulation time from t = 0, which for this analysis is
    # the first mosunetuzumab dose (Table S3 uses NONMEM's $DES variable T).
    cl <- clbase + (cl_exp_inf - clbase) * (1 - exp(-log(2) / cl_exp_thalf * time))

    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ----
    d/dt(central) <- -cl / vc * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Residual competing anti-CD20 antibody from prior therapy. These states
    # hold CONCENTRATIONS (ug/mL), are seeded from the patient's observed
    # predose value, and decay first-order at a fixed literature terminal
    # half-life. They are never dosed.
    ritux(0) <- CP_RITUXIMAB_UGML
    obin(0) <- CP_OBINUTUZUMAB_UGML
    d/dt(ritux) <- -kel_ritux * ritux
    d/dt(obin) <- -kel_obin * obin

    # ---- 5. Observation, receptor occupancy and cumulative endpoints ----
    Cc <- central / vc

    # Bender 2024 Equation 1 / Table S3 DADT(6). Competitive equilibrium
    # binding of mosunetuzumab, rituximab and obinutuzumab to CD20:
    #   RO% = 100*Cmosun / (Cmosun + KD_M + (KD_M/KD_R)*Critux + (KD_M/KD_G)*Cobin)
    # The published expression multiplies each concentration by 1000 to work in
    # ng/mL against ng/mL dissociation constants. Because RO% is a ratio in
    # which every term carries the same concentration units, the factor cancels
    # exactly; here all four concentrations and all three KDs are in ug/mL.
    RO <- 100 * Cc / (Cc + kd_mosun + (kd_mosun / kd_ritux) * ritux +
                        (kd_mosun / kd_obin) * obin)

    # Cumulative exposure endpoints, carried as states so that the paper's
    # AUC0-42 (Table S2, Figure 2) and average-RO metrics are available
    # directly from a solve. Bender 2024 Table S3 DADT(4) and DADT(6).
    d/dt(auc) <- Cc
    d/dt(auc_ro) <- RO

    Cc ~ lnorm(expSd)
  })
}
