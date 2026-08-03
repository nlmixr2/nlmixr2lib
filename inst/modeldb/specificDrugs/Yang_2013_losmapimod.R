Yang_2013_losmapimod <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination",
    "and time-dependent first-order absorption (change-point MTIME",
    "between two absorption rate constants XKA1 and XKA2) plus an",
    "absorption lag time ALAG1, for oral losmapimod (a p38 alpha/beta",
    "mitogen-activated protein kinase inhibitor). Data pooled from 30",
    "healthy volunteers, 23 patients with rheumatoid arthritis, and 24",
    "patients with chronic obstructive pulmonary disease across four",
    "GSK studies (Yang 2013). Retained covariates on inter-individual",
    "variability are sex on CL/F (males reference; females get a lower",
    "CL/F), age on V1 (power exponent centred at median age 54 years),",
    "and body weight on both absorption rate constants (shared power",
    "exponent centred at median 74 kg). Residual error is proportional",
    "and larger in COPD patients than in the pooled healthy + RA",
    "reference. Doses were 5, 7, 7.5, 10, and 20 mg oral single- or",
    "repeat-dose; 60 mg data were excluded from the analysis due to",
    "lack of dose-proportionality (possible saturable absorption).")
  reference   <- "Yang S, Lukey P, Beerahee M, Hoke F. Population pharmacokinetics of losmapimod in healthy subjects and patients with rheumatoid arthritis and chronic obstructive pulmonary diseases. Clin Pharmacokinet. 2013;52(3):187-198. doi:10.1007/s40262-012-0025-6"
  vignette    <- "Yang_2013_losmapimod"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on both absorption rate constants XKA1 and XKA2 with reference value 74 kg (median in the pooled 77-subject analysis dataset; Yang 2013 Table 3 footnote e). Applied as (WT/74)^-0.787. Body weight was tested but not retained on CL/F, V1, V2, or Q. Time-fixed at baseline.",
      source_name        = "Bodyweight"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on V1 (apparent central volume of distribution) with reference value 54 years (median in the pooled 77-subject analysis dataset; Yang 2013 Table 3 footnote d). Applied as (AGE/54)^-0.881. Age was tested but not retained on CL/F, Q, V2, XKA1, or XKA2. Time-fixed at baseline.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; Yang 2013 Table 3 footnote c 'Sex = male as reference').",
      notes              = "Log-linear effect on CL/F: log(CL/F_female) - log(CL/F_male) = -0.238 (Yang 2013 Table 3). Encoded via lcl + e_sex_cl * SEXF so that SEXF = 0 recovers exp(3.44) = 31.19 L/h for males and SEXF = 1 recovers exp(3.44 - 0.238) = 24.60 L/h for females (paper reports 31.2 L/h for males and 24.6 L/h for females in Section 3). Sex was tested but not retained on V1, Q, V2, XKA1, or XKA2. Time-fixed at baseline.",
      source_name        = "Sex"
    ),
    DIS_COPD = list(
      description        = "Chronic obstructive pulmonary disease patient indicator, 1 = COPD patient, 0 = non-COPD subject (healthy volunteer or rheumatoid arthritis patient).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-COPD; the reference group is the pooled healthy volunteers plus rheumatoid arthritis patients that share the smaller proportional residual variance).",
      notes              = "Not retained on any structural PK parameter (CL/F, V1, Q, V2, XKA1, XKA2, MTIME, ALAG1). Retained only on the proportional residual error magnitude: sigma^2_prop = 0.061 for non-COPD subjects vs sigma^2_prop_COPD = 0.268 for COPD patients (Yang 2013 Table 3 Residual variability rows). Encoded in model() by making the proportional SD a linear switch on DIS_COPD: propSdEff = (1 - DIS_COPD) * propSd + DIS_COPD * propSd_copd. This reproduces the NONMEM $ERROR pattern IF (COPD.EQ.1) Y = IPRED*(1+EPS(2)) ELSE Y = IPRED*(1+EPS(1)). Time-fixed at baseline.",
      source_name        = "Population indicator (Yang 2013 Section 2.3 dataset stratification)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 4L,
    age_range      = "22-74 years across the pooled analysis dataset (healthy: 22-53 y, median 43; RA: 34-73 y, median 57; COPD: 42-74 y, median 56)",
    age_median     = "54 years (median across the pooled analysis dataset; Yang 2013 Table 3 footnote d and reference age for the V1 power covariate)",
    weight_range   = "38.5-116 kg across the pooled analysis dataset (healthy: 56.7-100.8 kg, median 79; RA: 38.5-98 kg, median 69; COPD: 52-116 kg, median 80)",
    weight_median  = "74 kg (median across the pooled analysis dataset; Yang 2013 Table 3 footnote e and reference weight for the ka power covariate)",
    sex_female_pct = 45,
    race_ethnicity = NA_character_,
    disease_state  = "Pooled adult healthy volunteers (n = 30, 20 M / 10 F; studies MKI101678 and MKI102422) and patients with active rheumatoid arthritis on stable anti-rheumatic therapy (n = 23, 3 M / 20 F; study RA3103730 / NCT00256919) and patients with chronic obstructive pulmonary disease (n = 24, 19 M / 5 F; study MKI106209 / NCT00392587).",
    dose_range     = "5-20 mg oral losmapimod, single dose or repeated dosing: 5 mg or 10 mg once daily for 14 days followed by 7.5 mg twice daily for 14 days (healthy); 7.5, 20 mg single doses (RA); 7.5 mg once or twice daily for 14 days (COPD). 60 mg data (from healthy and RA cohorts) were excluded due to lack of dose-proportionality.",
    regions        = "GlaxoSmithKline multi-region program (United Kingdom, United States and other sites; specific country composition not tabulated in the paper).",
    n_observations = "1545 losmapimod plasma concentrations from 77 subjects in the analysis dataset (1247 from 30 healthy volunteers, 182 from 24 RA patients, 216 from 23 COPD patients; Yang 2013 Section 2.3).",
    notes          = "Analysis dataset used for model building. Losmapimod plasma concentrations were quantified by validated HPLC-MS/MS assay with LLQ = 0.2 ng/mL (Yang 2013 Section 2.2). An additional 204 concentrations from 34 RA subjects (7.5 mg BID x 28 days; study RA3103718 / NCT00393146) formed the external evaluation dataset (not part of model building). Fitting used NONMEM VI with FOCEI and the MTIME feature to estimate the absorption change point (Yang 2013 Sections 2.4-2.5)."
  )

  ini({
    # ============================================================
    # Final model parameter estimates from Yang 2013 Table 3.
    # Structural theta_pop values are reported on the natural-log
    # scale; typical PK parameters are recovered by exp(theta_pop).
    # See paper Section 2.4:
    #   theta_ij = exp(theta_pop,j + eta_ij)
    # so IIVs enter as additive shifts on log-scale parameters
    # (canonical mu-referenced encoding for FOCEI). Covariate
    # effects on the log scale reproduce the paper's Sect. 2.6
    # Eqs. 2-4 (continuous power covariates; categorical additive
    # log-shifts). The reference values (median age 54 y, median
    # weight 74 kg, male sex) come from Table 3 footnotes c-e.
    # ============================================================

    # ---------------- Structural (typical-value) parameters ----------------
    lcl    <- 3.44
    label("Log apparent total oral clearance CL/F at male reference (log(L/h))")     # Yang 2013 Table 3: theta_pop,1 = 3.44 (CL/F male reference = exp(3.44) = 31.19 L/h)
    lvc    <- 3.53
    label("Log apparent central volume of distribution V1 at median age 54 y (log(L))")  # Yang 2013 Table 3: theta_pop,2 = 3.53 (V1 at reference = exp(3.53) = 34.12 L)
    lq     <- 2.42
    label("Log apparent inter-compartmental clearance Q (log(L/h))")                  # Yang 2013 Table 3: theta_pop,3 = 2.42 (Q = exp(2.42) = 11.25 L/h)
    lvp    <- 4.80
    label("Log apparent peripheral volume of distribution V2 (log(L))")               # Yang 2013 Table 3: theta_pop,4 = 4.80 (V2 = exp(4.80) = 121.51 L)
    lka1   <- -2.02
    label("Log slow-phase absorption rate constant XKA1 at median weight 74 kg (log(1/h))")  # Yang 2013 Table 3: theta_pop,5 = -2.02 (XKA1 at reference = exp(-2.02) = 0.133 1/h; slow phase before MTIME)
    lka2   <- -1.62
    label("Log fast-phase absorption rate constant XKA2 at median weight 74 kg (log(1/h))")  # Yang 2013 Table 3: theta_pop,6 = -1.62 (XKA2 at reference = exp(-1.62) = 0.198 1/h; fast phase after MTIME)
    lmtime <- log(0.811)
    label("Log time of absorption ka change-point relative to dose (log(h))")         # Yang 2013 Table 3: MTIME = 0.811 h (fixed effect only, no IIV reported)
    ltlag  <- log(0.233)
    label("Log absorption lag time ALAG1 (log(h))")                                    # Yang 2013 Table 3: ALAG1 = 0.233 h (fixed effect only, no IIV reported)

    # ---------------- Covariate effects (Yang 2013 Table 3 covariate rows) ----------------
    # Sex on CL/F: log-linear categorical (Yang 2013 Section 2.6 Eq. 3;
    # male = reference per Table 3 footnote c). log(CL/F_female) -
    # log(CL/F_male) = -0.238, i.e. female CL/F = exp(-0.238) = 0.788 x
    # male CL/F. Paper reports 31.2 L/h for males and 24.6 L/h for
    # females (24.6 / 31.2 = 0.789).
    e_sex_cl <- -0.238
    label("Log-linear coefficient of SEXF on CL/F (unitless; female shift on log(CL/F))")  # Yang 2013 Table 3: Sex coefficient on CL/F row = -0.238 (%RSE 34)

    # Age on V1: continuous power (Yang 2013 Section 2.6 Eq. 2) with
    # reference median age = 54 years (Table 3 footnote d).
    # V1 = exp(theta_pop,2) * (AGE/54)^-0.881.
    e_age_vc <- -0.881
    label("Power exponent of (AGE/54) on V1 (unitless)")                                   # Yang 2013 Table 3: Age coefficient on V1 row = -0.881 (%RSE 40.3)

    # Bodyweight on ka: continuous power (Yang 2013 Section 2.6 Eq. 2)
    # with reference median weight = 74 kg (Table 3 footnote e). The
    # single "Bodyweight" row in Table 3 sits within the block that
    # defines both XKA1 = exp(theta_pop,5) and XKA2 = exp(theta_pop,6),
    # so the same exponent applies to both absorption rate constants.
    e_wt_ka <- -0.787
    label("Power exponent of (WT/74) shared by XKA1 and XKA2 (unitless)")                  # Yang 2013 Table 3: Bodyweight coefficient in ka block = -0.787 (%RSE 28.6)

    # ---------------- Between-subject variability (Table 3 IIV rows) ----------------
    # Yang 2013 reports omega^2 directly ("Inter-individual variability
    # omega^2" column; the %CV column back-calculates as
    # 100 * sqrt(omega^2), e.g. sqrt(0.131) = 0.362 -> 36.2% for CL/F).
    # MTIME and ALAG1 have no IIV row -- they are typical-value only.
    etalcl  ~ 0.131                                                                        # Yang 2013 Table 3: omega^2(CL/F) = 0.131 (36.2%)
    etalvc  ~ 0.167                                                                        # Yang 2013 Table 3: omega^2(V1) = 0.167 (40.9%)
    etalq   ~ 0.257                                                                        # Yang 2013 Table 3: omega^2(Q) = 0.257 (50.7%)
    etalvp  ~ 0.230                                                                        # Yang 2013 Table 3: omega^2(V2) = 0.230 (48%)
    etalka1 ~ 1.44                                                                         # Yang 2013 Table 3: omega^2(XKA1) = 1.44 (120% -- large by design for the slow-phase absorption)
    etalka2 ~ 0.0268                                                                       # Yang 2013 Table 3: omega^2(XKA2) = 0.0268 (16.4%)

    # ---------------- Residual variability (Table 3 residual rows) ----------------
    # Two proportional residual variances estimated jointly, mutually
    # exclusive on subject population per Yang 2013 Section 3 first
    # paragraph ("A slightly higher residual error was associated with
    # COPD patients compared with the healthy and RA subjects"). NONMEM
    # $ERROR pattern:
    #   IF (COPD.EQ.1) Y = IPRED * (1 + EPS(2))
    #   ELSE           Y = IPRED * (1 + EPS(1))
    # with SIGMA(1,1) = 0.061 (%CV 24.7) and SIGMA(2,2) = 0.268
    # (%CV 51.8). Converted to standard-deviation form for nlmixr2's
    # prop() error syntax by sqrt(sigma^2).
    propSd      <- sqrt(0.061)
    label("Proportional residual SD on Cc for non-COPD subjects (healthy or RA; fraction)")   # Yang 2013 Table 3: sigma^2_prop = 0.061 (24.7%; %RSE 10.6)
    propSd_copd <- sqrt(0.268)
    label("Proportional residual SD on Cc for COPD subjects (fraction)")                       # Yang 2013 Table 3: sigma^2_prop,COPD = 0.268 (51.8%; %RSE 17.2)
  })

  model({
    # ------------------------------------------------------------------
    # Reference values for covariate scaling (Yang 2013 Table 3 footnotes
    # c-e). These are constants -- do not add IIV.
    # ------------------------------------------------------------------
    ref_age <- 54    # Yang 2013 Table 3 footnote d: "Median age = 54 years"
    ref_wt  <- 74    # Yang 2013 Table 3 footnote e: "Median bodyweight = 74 kg"

    # ------------------------------------------------------------------
    # Individual PK parameters. Covariate effects follow Yang 2013
    # Section 2.6 Eqs. 2-4:
    #   * Continuous power: log(TVPK) = theta_pop + theta_cov * log(cov/median_cov)
    #   * Categorical:       log(TVPK) = theta_pop + theta_cov * indicator
    # Reference conditions when all covariates are at reference:
    #   * Male (SEXF = 0), age 54 y, weight 74 kg -> CL/F = 31.19 L/h,
    #     V1 = 34.12 L, Q = 11.25 L/h, V2 = 121.51 L,
    #     XKA1 = 0.133 h^-1, XKA2 = 0.198 h^-1.
    # ------------------------------------------------------------------
    cl  <- exp(lcl + etalcl + e_sex_cl * SEXF)
    vc  <- exp(lvc + etalvc) * (AGE / ref_age)^e_age_vc
    q   <- exp(lq  + etalq)
    vp  <- exp(lvp + etalvp)
    ka1 <- exp(lka1 + etalka1) * (WT / ref_wt)^e_wt_ka
    ka2 <- exp(lka2 + etalka2) * (WT / ref_wt)^e_wt_ka

    # Change-point time and absorption lag (typical-value only per Table 3).
    # The published parameter name MTIME collides with an rxode2 internal
    # symbol (same shadowing as Hajjar 2018 DMD 6MWT), so the change-point
    # time is carried here as `t_change`.
    t_change <- exp(lmtime)
    tlag     <- exp(ltlag)

    # ------------------------------------------------------------------
    # Time-dependent absorption rate constant (Yang 2013 Section 2.4):
    #   MTDIFF = 1   ; single change point
    #   MTIME  = theta ; time of change (estimated)
    #   ka = XKA1 * (1 - MPAST) + XKA2 * MPAST
    # where MPAST = 0 for TIME < MTIME and 1 for TIME >= MTIME. NONMEM's
    # MTIME feature is subject-level (fires once per subject at absolute
    # time MTIME). In rxode2 the equivalent indicator is (time >= t_change)
    # on the simulation clock.
    # ------------------------------------------------------------------
    ka <- ka1 * (time < t_change) + ka2 * (time >= t_change)

    # ------------------------------------------------------------------
    # Two-compartment ODE system parameterised in CL/V terms (paper
    # uses NONMEM ADVAN4 TRANS1: CL/F, V1, Q, V2 with first-order
    # absorption from a depot).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag on the oral depot (Yang 2013 Table 3 ALAG1 = 0.233 h).
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # Observation and population-dependent proportional residual error.
    # Concentration units: dose in mg, volumes in L, so central / vc is
    # mg/L. Multiply by 1000 to convert to ng/mL (the paper's reporting
    # unit; Yang 2013 Figs. 1-2).
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc

    # Population-dependent proportional SD: SD reference (non-COPD) for
    # DIS_COPD = 0; SD_copd for DIS_COPD = 1. Reproduces the NONMEM
    # IF (COPD.EQ.1) EPS(2) ELSE EPS(1) pattern by making the SD a
    # linear switch on the covariate.
    propSdEff <- (1 - DIS_COPD) * propSd + DIS_COPD * propSd_copd

    Cc ~ prop(propSdEff)
  })
}
