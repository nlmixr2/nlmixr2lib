vanLieshout_2026_bedaquiline <- function() {
  description <- "Translational PKPD model predicting the clinical early bactericidal activity (EBA) of oral bedaquiline in adults with pulmonary tuberculosis from a zebrafish-larva exposure-response relationship. The Svensson 2016 population PK model (two-transit absorption, three-compartment bedaquiline plus one-compartment metabolite M2, allometric weight, age and Black-race effects on clearance, and a logistic serum-albumin trajectory that scales clearance, fraction metabolised and the volumes) drives first-order killing of sputum Mycobacterium tuberculosis. The Emax kill rate is 109% (first 72 h) then 137% of the zebrafish M. marinum growth rate (0.0853 1/h), the zebrafish whole-body EC50 (26.5 mg/L) is scaled by the M. tuberculosis / M. marinum MIC ratio (0.0415 / 0.1 mg/L), and net M. tuberculosis growth is zero. All parameters are fixed; the bacterial load starts at the subject's baseline log10 CFU/mL."
  reference <- paste(
    "van Lieshout B, Hanuun A, van Hasselt JGC, van Wijk RC. (2026).",
    "Successful Translational Prediction of Bedaquiline's Phase 2a Clinical",
    "Trial Results against Tuberculosis: The Zebrafish as a New Approach",
    "Methodology. Clin Pharmacol Ther 120(5):1281. doi:10.1002/cpt.70430.",
    "PK structure and values from Svensson EM, Dosne AG, Karlsson MO (2016)",
    "CPT Pharmacometrics Syst Pharmacol 5(12):682-691, doi:10.1002/psp4.12147;",
    "see modellib('Svensson_2016_bedaquiline'). Zebrafish PD from Forn-Cuni",
    "et al. (2026) bioRxiv doi:10.64898/2026.05.29.728504.",
    sep = " "
  )
  vignette <- "vanLieshout_2026_bedaquiline"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L (Cc, Cc_m2); log10 CFU/mL (log_cfu)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling around 70 kg: (WT/70)^0.180912 on CL, Q1, Q2 and CLM2 and (WT/70)^1 on V, VP1, VP2 and VM2 (Supplementary Materials model code, ALL_CL / ALL_V; THETA(19), THETA(20)). The model uses each subject's own weight in place of the Svensson 2016 weight-trajectory model; the study mean weight was used where individual data were not available (Methods, 'Simulations of the clinical pharmacokinetics').",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL and CLM2 centred on 32 years: 1 + (32 - AGE) * 0.00880756 (model code AGECL; THETA(24)).",
      source_name = "AGE"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Multiplicative effect on CL and CLM2: 1 + 0.838743 * RACE_BLACK (model code BLACKCL; THETA(23)). The source column RACE is coded 2 = Black; RACE_BLACK = 1 when RACE == 2.",
      source_name = "RACE"
    ),
    SPUTUM_BL_LOG10CFU = list(
      description = "Baseline (pre-treatment) sputum bacterial load on the log10 CFU/mL scale",
      units = "log10 CFU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sets the initial bacterial state, cfu(0) = 10^SPUTUM_BL_LOG10CFU (model code A_0(9) = 10**BASELINE). For the 14-day study the authors supplied individual empirical Bayes estimates from a mixed-effects intercept model on pre-treatment counts (regimen means 5.96, 5.69, 5.82 and 6.20 log10 CFU/mL for the 100, 200, 300 and 400 mg groups, CV 98.6%); for the 7-day study, values simulated from the reported regimen means (6.66, 6.32 and 6.82 log10 CFU/mL for 25, 100 and 400 mg; Results, 'Simulated exposures and baselines'). Because killing is first order in the bacterial load, the baseline shifts the log10 CFU/mL curve but does not change the EBA.",
      source_name = "BASELINE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bedaquiline", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "bedaquiline", units = "mg", specimen = "tissue", verified = TRUE),
    central_m2 = list(analyte = "M2 (N-desmethyl bedaquiline)", units = "mg", specimen = "plasma", verified = TRUE),
    cfu = list(analyte = "Mycobacterium tuberculosis", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 105L,
    n_studies = 2L,
    disease_state = "Adults with sputum smear-positive pulmonary tuberculosis not previously treated with bedaquiline, enrolled in two randomized, double-blind, placebo-controlled phase 2a EBA trials of bedaquiline monotherapy.",
    dose_range = "Rustomjee 2008 (7 days, aggregate data, n = 45): 25, 100 or 400 mg once daily. Diacon 2013 (14 days, individual data from TB-PACTS, n = 60): loading doses on days 1 and 2 then daily maintenance, 200/100/100, 400/300/200, 500/400/300 or 700/500/400 mg (labelled by the maintenance dose, 100-400 mg).",
    notes = "The model was not fitted to these trials: it is a forward prediction that combines the Svensson 2016 popPK model (335 patients with MDR-TB in the C208 and C209 trials) with a zebrafish-larva (Danio rerio, M. marinum infection, up to 5 days post-fertilization, 0-50 uM waterborne bedaquiline) exposure-response relationship. Only the baseline bacterial loads were estimated from the clinical data. The paper does not tabulate demographics of either trial; it used the study-mean covariates (7-day study) or individual covariates (14-day study)."
  )

  ini({
    # Every value is a $THETA / $OMEGA / $SIGMA entry of the Supplementary
    # Materials NONMEM model code. All $THETA are FIX in that code, which is a
    # simulation-only translation of Svensson 2016 (THETA(9)-THETA(24) equal
    # Svensson 2016 Table 3).

    # --- Absorption: two transit compartments then first-order absorption.
    lmat <- fixed(log(6 * 0.662045)); label("Mean absorption time, delay plus absorption (h)") # THETA(9) TVMAT = 0.662045, a fraction of 6 h: MAT = 6*expit(logit(TVMAT)) with IOV fixed to 0
    fr <- fixed(0.466443); label("Fraction of the mean absorption time spent in the transit delay (unitless)") # THETA(10) FR

    # --- Bedaquiline and M2 disposition (apparent values, 70 kg, mg/L units).
    lcl <- fixed(log(2.61592)); label("Bedaquiline clearance CL/F (L/h)") # THETA(11)
    lvc <- fixed(log(198.34)); label("Bedaquiline central volume V/F (L)") # THETA(12)
    lq <- fixed(log(3.658)); label("Bedaquiline intercompartmental clearance to peripheral 1, Q1/F (L/h)") # THETA(13)
    lvp <- fixed(log(8549.06)); label("Bedaquiline peripheral volume 1, VP1/F (L)") # THETA(14)
    lq2 <- fixed(log(7.33504)); label("Bedaquiline intercompartmental clearance to peripheral 2, Q2/F (L/h)") # THETA(15)
    lvp2 <- fixed(log(2690.91)); label("Bedaquiline peripheral volume 2, VP2/F (L)") # THETA(16)
    lcl_m2 <- fixed(log(10.0496)); label("M2 clearance CLM2/(F*fm) (L/h)") # THETA(17)
    lvc_m2 <- fixed(log(2203.71)); label("M2 volume VM2/(F*fm) (L)") # THETA(18)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability (unitless)") # model code F1 = 1*EXP(IOV_F+IIV_F)

    # --- Covariate effects.
    e_wt_cl <- fixed(0.180912); label("Allometric weight exponent on CL, Q1, Q2 and CLM2 (unitless)") # THETA(19)
    e_wt_vc <- fixed(1); label("Allometric weight exponent on V, VP1, VP2 and VM2 (unitless)") # THETA(20)
    e_alb_vc <- fixed(1); label("Exponent of the albumin-driven relative unbound fraction on the observation volumes (unitless)") # THETA(21)
    e_alb_cl <- fixed(1.64021); label("Exponent of relative albumin on CL, CLM2 and fm (unitless)") # THETA(22)
    e_race_black_cl <- fixed(0.838743); label("Fractional increase in CL and CLM2 for Black race (unitless)") # THETA(23)
    e_age_cl <- fixed(0.00880756); label("Linear age effect on CL and CLM2 per year below 32 years (1/year)") # THETA(24)

    # --- Serum albumin trajectory (logistic return to steady state).
    lalb0 <- fixed(log(3.64865)); label("Typical serum albumin at the start of treatment (g/dL)") # THETA(1) X0
    lalb_ss <- fixed(log(4.04068)); label("Typical steady-state serum albumin (g/dL)") # THETA(2) Xss
    lkrep_alb <- fixed(log(0.03399109)); label("Rate constant of albumin return to steady state (1/week)") # THETA(3)
    bc_alb0 <- fixed(-2.43936); label("Box-Cox shape of the IIV on starting albumin (unitless)") # THETA(4)
    bc_alb_ss <- fixed(-5.37749); label("Box-Cox shape of the IIV on steady-state albumin (unitless)") # THETA(5)

    # --- Zebrafish-derived pharmacodynamics (model code, 'Translation PD').
    lkg_mm <- fixed(log(0.0853)); label("M. marinum growth rate constant in zebrafish larvae (1/h)") # GROWTH_Mm = 0.0853; Methods '109-137% of 0.0853 h-1'
    kg <- fixed(0); label("Net M. tuberculosis growth rate constant in chronic infection (1/h)") # GROWTH_MTB = 0
    f_emax_kg <- fixed(1.09); label("Maximum kill rate as a multiple of the M. marinum growth rate, before the switch time (unitless)") # EMAX_BASE = 1.09; Methods '109% of the bacterial growth rate constant'
    f_emax_kg_late <- fixed(1.37); label("Maximum kill rate as a multiple of the M. marinum growth rate, after the switch time (unitless)") # IF(T.GT.72) EMAX_BASE = 1.37; Methods 'increased to 137%'
    t_emax_switch <- fixed(72); label("Time after which the higher maximum kill rate applies (h)") # IF(T.GT.72)
    lec50_zf <- fixed(log(26.5)); label("Bedaquiline EC50 in zebrafish larvae, whole-body concentration (mg/L)") # BDQ_EC50 = 26.5 * ...; Methods 'EC50 of 26.5 mg/L'
    mic_mtb <- fixed(0.04145085); label("Bedaquiline MIC for M. tuberculosis, geometric mean of the EUCAST distribution (mg/L)") # MIC_MTB_TARGET; Methods '~0.04 mg/L'
    mic_mm <- fixed(0.1); label("Bedaquiline MIC for M. marinum (mg/L)") # MIC_Mm_ORIGIN; Methods '0.1 mg/L was used'
    hill <- fixed(1); label("Hill coefficient of the kill-rate Emax function (unitless)") # BDQ_GAMMA = 1

    # --- Between-subject variability (variances). Between-occasion
    # variability on F and MAT is 0 FIX in the model code ('inter-occasion
    # variability was omitted', Methods).
    etalalb0 + etalalb_ss + etalkrep_alb ~ c(
      0.0253964,
      0.00787165, 0.00965573,
      -0.0831085, 0.00887115, 0.979314
    ) # $OMEGA BLOCK(5) elements 1-3 (IIV_X0, IIV_XSS, IIV_REP); elements 4-5 belong to the unused weight model
    etalfdepot ~ 0.0803271 # $OMEGA 10 IIV_F
    etalcl + etalcl_m2 ~ c(
      0.152776,
      0.13486, 0.212124
    ) # $OMEGA BLOCK(2) 11-12 IIV_CL, IIV_CLM2
    etalvc ~ 0.171909 # $OMEGA 13 IIV_V
    etalq ~ 0.181182 # $OMEGA 14 IIV_Q1
    etalvc_m2 ~ 0.150223 # $OMEGA 15 IIV_VM2

    # --- Residual error: additive on log concentration (Y = log(C) + EPS).
    # The bacterial load carries no residual error ('0 RUV for bacteria').
    expSd <- sqrt(0.0518161); label("Bedaquiline residual error, additive on the log scale (SD)") # $SIGMA BLOCK(2) element 3 'Prop error BDQ' variance 0.0518161
    expSd_m2 <- sqrt(0.0366836); label("M2 residual error, additive on the log scale (SD)") # $SIGMA BLOCK(2) element 4 'Prop error M2' variance 0.0366836
  })

  model({
    # Serum albumin (g/dL). The model code integrates the logistic
    # dALB/dt = krep * ALB * (1 - ALB/Xss) from ALB(0) = X0, with krep in
    # 1/week converted to 1/h by 7*24; the closed form below is the exact
    # solution, with t = 0 at the start of treatment. X0 and Xss carry
    # Box-Cox-transformed log-normal IIV.
    bcx_alb0 <- ((exp(etalalb0))^bc_alb0 - 1) / bc_alb0
    bcx_alb_ss <- ((exp(etalalb_ss))^bc_alb_ss - 1) / bc_alb_ss
    alb0 <- exp(lalb0 + bcx_alb0)
    alb_ss <- exp(lalb_ss + bcx_alb_ss)
    krep_alb <- exp(lkrep_alb + etalkrep_alb) / (7 * 24)
    alb <- alb_ss / (1 + (alb_ss / alb0 - 1) * exp(-krep_alb * t))

    # Albumin effects. alb_rel is relative to the subject's own steady state
    # (ALB_RELI); fu_rel is relative to the typical population steady state
    # (FURELIP = THETA(2)/A(8)) and enters only the observation volumes.
    alb_rel <- alb / alb_ss
    cov_alb <- alb_rel^e_alb_cl
    alb_fm_fac <- alb_rel^(-e_alb_cl)
    cov_fu <- (exp(lalb_ss) / alb)^e_alb_vc

    blackcl <- 1 + e_race_black_cl * RACE_BLACK
    agecl <- 1 + (32 - AGE) * e_age_cl
    allcl <- (WT / 70)^e_wt_cl
    allv <- (WT / 70)^e_wt_vc

    # Absorption: MAT = delay (fraction fr, two transit compartments with rate
    # ktr) plus first-order absorption whose half-life is 1/3.3 of the rest.
    mat <- exp(lmat)
    mtt <- mat * fr
    ka <- log(2) / (mat * (1 - fr) / 3.3)
    ktr <- 2 / mtt

    cl <- exp(lcl + etalcl) * blackcl * agecl * allcl * cov_alb
    vc <- exp(lvc + etalvc) * allv
    q <- exp(lq + etalq) * allcl
    vp <- exp(lvp) * allv
    q2 <- exp(lq2) * allcl
    vp2 <- exp(lvp2) * allv
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * blackcl * agecl * allcl / alb_fm_fac * cov_alb
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * allv / alb_fm_fac

    # Zebrafish-to-human translation (model code BDQ_EMAX, BDQ_EC50): the
    # maximum kill rate is a multiple of the M. marinum growth rate and the
    # EC50 is scaled by the MIC ratio M. tuberculosis / M. marinum.
    kg_mm <- exp(lkg_mm)
    f_emax <- f_emax_kg
    if (t > t_emax_switch) f_emax <- f_emax_kg_late
    emax <- f_emax * kg_mm
    ec50 <- exp(lec50_zf) * mic_mtb / mic_mm

    # The kill rate is driven by central / V without the unbound-fraction
    # volume correction (CONC = A(2)/V in $DES).
    conc_pd <- central / vc
    kill <- emax * conc_pd^hill / (conc_pd^hill + ec50^hill)

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ka * transit2
    d/dt(central) <- ka * transit2 - cl / vc * central -
      q / vc * central + q / vp * peripheral1 -
      q2 / vc * central + q2 / vp2 * peripheral2
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(peripheral2) <- q2 / vc * central - q2 / vp2 * peripheral2
    d/dt(central_m2) <- cl / vc * central - cl_m2 / vc_m2 * central_m2
    d/dt(cfu) <- kg * cfu - kill * cfu

    f(depot) <- exp(lfdepot + etalfdepot)
    cfu(0) <- 10^SPUTUM_BL_LOG10CFU

    # Observed concentrations use the albumin-corrected volumes of $ERROR
    # (VE = VB*ALL_VE*COVFUIP, VM2E = VM2B*ALL_VE*COVFUIP/FME).
    Cc <- central / (vc * cov_fu)
    Cc_m2 <- central_m2 / (vc_m2 * cov_fu)
    log_cfu <- log10(cfu)

    Cc ~ lnorm(expSd)
    Cc_m2 ~ lnorm(expSd_m2)
  })
}
