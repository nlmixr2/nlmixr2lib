Diao_2016_daclizumab_cd56bright <- function() {
  description <- "Indirect-response PK/PD model of CD56 bright natural killer (NK) cell expansion following subcutaneous daclizumab high-yield process (HYP) in adults with relapsing-remitting multiple sclerosis (Diao 2016). Daclizumab HYP serum concentration stimulates the zero-order production rate (Kin) of CD56 bright NK cells (% of all lymphocytes) via a saturable Smax function; first-order elimination rate Kout is fixed by the median baseline. The PK backbone is the two-compartment, first-order SC absorption + lag model from Othman 2014 (file inst/modeldb/specificDrugs/Othman_2014_daclizumab.R), copied verbatim with weight-based allometric scaling."
  reference <- "Diao L, Hang Y, Othman AA, Nestorov I, Tran JQ, Mehta D, Amaravadi L. Population PK/PD analyses of CD25 occupancy, CD56 bright NK cell expansion and regulatory T cell reduction by daclizumab HYP in subjects with multiple sclerosis. Br J Clin Pharmacol. 2016;82(5):1333-1342. doi:10.1111/bcp.13051 (PMID 27333593). PK backbone: Othman AA, Tran JQ, Tang MT, Dutta S. Population Pharmacokinetics of Daclizumab High-Yield Process in Healthy Volunteers. Clin Pharmacokinet. 2014;53(10):907-918. doi:10.1007/s40262-014-0159-9."
  vignette <- "Diao_2016_daclizumab_cd56bright"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", response = "% of total lymphocytes (T + B + NK)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of the inherited Othman 2014 PK parameters (CL, Q, Vc, Vp) with reference 70 kg; exponents 0.54 on CL/Q and 0.64 on Vc/Vp. CD56 bright NK PD parameters do not carry weight covariates in Diao 2016.",
      source_name        = "WT"
    ),
    DOSE_50MG = list(
      description        = "Record-level indicator for the 50 mg SC dose (1 = 50 mg SC, 0 = any other SC dose or any IV dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (100, 150, 200, or 300 mg SC dose, or any IV dose)",
      notes              = "Inherited from the Othman 2014 PK backbone. Diao 2016 dosing is 150 or 300 mg SC every 4 weeks, so leave DOSE_50MG = 0 in clinical simulations.",
      source_name        = "(derived from AMT)"
    )
  )

  population <- list(
    n_subjects     = 1405L,
    n_records      = 9630L,
    n_studies      = 4L,
    study_names    = c("205MS201 / SELECT (Phase 2, RRMS)",
                       "205MS202 / SELECTION (Phase 2 extension with washout cohort)",
                       "205MS302 / OBSERVE (immunogenicity / PK / PD with intensive substudy)",
                       "205MS301 / DECIDE (Phase 3 vs IFN beta-1a)"),
    disease_state  = "Relapsing-remitting multiple sclerosis (RRMS)",
    dose_range     = "Daclizumab HYP 150 or 300 mg SC every 4 weeks",
    notes          = paste0(
      "Pooled PK/PD dataset of 1405 RRMS subjects with 9630 CD56 bright NK ",
      "cell records from four daclizumab HYP clinical studies (Diao 2016 ",
      "Table 2). Median (mean) baseline CD56 bright NK is 0.6% (0.75%) of ",
      "all lymphocytes (Diao 2016 Results, NK section)."
    ),
    pd_subgroups = list(
      `205MS201/202 (SELECT/SELECTION)` = list(subjects = 561L, records = 5071L),
      `205MS302 (OBSERVE)`              = list(subjects = 107L, records =  922L),
      `205MS301 (DECIDE)`                = list(subjects = 737L, records = 3637L)
    ),
    baseline_cd56bright_nk_pct = list(median = 0.6, mean = 0.75)
  )

  ini({
    # ----------------------------------------------------------------------
    # PK backbone (copied verbatim from Othman_2014_daclizumab.R).
    # ----------------------------------------------------------------------
    lka      <- log(0.009 * 24); label("Absorption rate constant (Ka, 1/day)")                  # Othman 2014 Table 2
    lcl      <- log(0.010 * 24); label("Clearance for a 70 kg adult (CL, L/day)")               # Othman 2014 Table 2
    lvc      <- log(3.89);       label("Central volume of distribution for a 70 kg adult (Vc, L)") # Othman 2014 Table 2
    lvp      <- log(2.52);       label("Peripheral volume of distribution for a 70 kg adult (Vp, L)") # Othman 2014 Table 2
    lq       <- log(0.044 * 24); label("Inter-compartmental clearance for a 70 kg adult (Q, L/day)") # Othman 2014 Table 2
    lfdepot  <- log(0.84);       label("SC bioavailability for 100-300 mg doses (F, fraction)") # Othman 2014 Table 2
    lalag    <- log(2 / 24);     label("Absorption lag time for SC doses (Tlag, day; 2 h)")     # Othman 2014 Table 2

    allo_cl <- 0.54; label("Allometric exponent on CL and Q (unitless)") # Othman 2014 Table 2
    allo_v  <- 0.64; label("Allometric exponent on Vc and Vp (unitless)") # Othman 2014 Table 2

    e_dose_50mg_f <- -0.32143; label("Relative change in F for 50 mg SC vs 100-300 mg SC (fraction)") # Othman 2014 Table 2

    # Othman 2014 Table 2 (ka 58% CV, cl 27% CV, corr -0.72): ka -> omega^2 = log(1 + 0.58^2) = 0.29003;
    # cl -> omega^2 = log(1 + 0.27^2) = 0.07038; cov = -0.72 * sqrt(0.29003 * 0.07038) = -0.10290.
    etalka + etalcl ~ c(0.29003,
                        -0.10290, 0.07038)
    etalvc ~ 0.09175  # Othman 2014 Table 2 (Vc 31% CV)

    propSd <- 0.22; label("Proportional residual error on daclizumab HYP serum concentration (fraction)") # Othman 2014 Table 2
    addSd  <- 0.33; label("Additive residual error on daclizumab HYP serum concentration (ug/mL)")        # Othman 2014 Table 2

    # ----------------------------------------------------------------------
    # CD56 bright NK cell PD parameters (Diao 2016 Table 4, indirect-response).
    # Equation 2: dNK/dt = Kin * (1 + Smax * Cc / (EC50 + Cc)) - Kout * NK
    # Kout was constrained to Kin / baseline_NK in the source fitting (Diao 2016
    # Methods: "Kout was expressed as a function of Kin using baseline information
    # ... and therefore was not independently estimated").
    # The library implementation reproduces that constraint by computing kout
    # from kin and the median baseline of 0.6% inside model().
    # ----------------------------------------------------------------------
    # Kin reported as 4.12e-04 %/h. Convert to %/day = 4.12e-04 * 24 = 0.009888 %/day.
    lcd56Kin  <- log(4.12e-04 * 24); label("Zero-order CD56 bright NK production rate (Kin, %/day)") # Diao 2016 Table 4 (Kin = 4.12e-04 %/h)
    lcd56Smax <- log(7.89);          label("Maximum stimulatory factor on Kin (Smax, unitless)")     # Diao 2016 Table 4 (Smax = 7.89)
    lcd56EC50 <- log(18.0);          label("Daclizumab Cc giving 50% of Smax (EC50, mg/L = ug/mL)")  # Diao 2016 Table 4 (EC50 = 18.0 mg/L)

    cd56baseline <- fixed(0.6); label("Median baseline CD56 bright NK (% of total lymphocytes)")  # Diao 2016 Results (NK section): "median ... baseline ... is 0.6%"

    # IIV (Diao 2016 Table 4): Kin 97% CV, Smax 67% CV. EC50 has no IIV reported.
    # omega^2 = log(1 + CV^2)
    # Kin  CV 97% -> omega^2 = log(1 + 0.97^2) = 0.66022
    # Smax CV 67% -> omega^2 = log(1 + 0.67^2) = 0.36398
    etalcd56Kin  ~ 0.66022    # Diao 2016 Table 4 (Kin IIV 97% CV)
    etalcd56Smax ~ 0.36398    # Diao 2016 Table 4 (Smax IIV 67% CV)

    # Residual error: proportional on CD56 bright NK %, 29.1%.
    propSd_cd56bright <- 0.291; label("Proportional residual error on CD56 bright NK (fraction)") # Diao 2016 Table 4
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (Othman 2014 PK backbone).
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_v
    vp <- exp(lvp)          * (WT / 70)^allo_v
    q  <- exp(lq)           * (WT / 70)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 2. Two-compartment SC PK ODE system.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- exp(lalag)
    f(depot)    <- exp(lfdepot) * (1 + e_dose_50mg_f * DOSE_50MG)

    Cc <- central / vc

    # ------------------------------------------------------------------
    # 3. Individual CD56 bright NK PD parameters.
    # ------------------------------------------------------------------
    cd56KinI  <- exp(lcd56Kin  + etalcd56Kin)
    cd56SmaxI <- exp(lcd56Smax + etalcd56Smax)
    cd56EC50I <- exp(lcd56EC50)

    # Kout derived from Kin and the baseline (Diao 2016 Methods).
    cd56KoutI <- cd56KinI / cd56baseline

    # ------------------------------------------------------------------
    # 4. Indirect-response ODE for CD56 bright NK cell percentage
    #    (Diao 2016 Equation 2). Uses the canonical `effect` compartment;
    #    `cd56bright` is the named observation.
    # ------------------------------------------------------------------
    d/dt(effect) <- cd56KinI * (1 + cd56SmaxI * Cc / (cd56EC50I + Cc)) -
                    cd56KoutI * effect
    effect(0)    <- cd56baseline

    cd56bright <- effect

    # ------------------------------------------------------------------
    # 5. Observation and error model.
    # ------------------------------------------------------------------
    Cc         ~ add(addSd) + prop(propSd)
    cd56bright ~ prop(propSd_cd56bright)
  })
}
