Langdon_2010_PF00821385_human <- function() {
  description <- paste(
    "First-in-human popPK-PD model for PF-00821385, a Pfizer HIV-1 gp120",
    "cell-fusion inhibitor candidate (molecular weight 440.49 g/mol),",
    "developed in 24 healthy male volunteers from a single-ascending-dose",
    "study (Langdon 2010 Tables 2 and Figure 3). PK is a two-compartment",
    "model with first-order oral absorption and an additive residual error",
    "on the log-transformed plasma concentrations (i.e., a log-normal",
    "residual). PD describes supine pulse rate as the sum of (a) a typical-",
    "value baseline rate with log-normal inter-subject variability, (b) a",
    "24-h cosine circadian rhythm with typical-value amplitude and log-",
    "normal inter-subject variable peak time, and (c) a linear drug effect",
    "on free plasma concentration with no IIV. The PD-slope SLOPE = 0.76",
    "bpm per micromolar free drug is from Langdon 2010 Table 2; plasma",
    "unbound fraction fu = 0.64 is FIXED via back-calculation from the",
    "canine toxicology Cmax data (see vignette Errata). The PD layer was",
    "fit sequentially to individual Bayesian post hoc PK estimates from",
    "the population PK fit (NONMEM VI with FOCE INTER; 500-iteration",
    "nonparametric bootstrap for SE / CI).",
    sep = " "
  )
  reference <- paste(
    "Langdon G, Davis JD, McFadyen LM, Dewhurst M, Brunton NS, Rawal JK,",
    "Van der Graaf PH, Benson N.",
    "Translational pharmacokinetic-pharmacodynamic modelling; application",
    "to cardiovascular safety data for PF-00821385, a novel HIV agent.",
    "Br J Clin Pharmacol. 2010 Apr;69(4):335-345.",
    "doi:10.1111/j.1365-2125.2009.03594.x.",
    sep = " "
  )
  vignette <- "Langdon_2010_PF00821385"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "21 - 55 years",
    weight_range   = "> 50 kg",
    sex_female_pct = 0.0,
    race_ethnicity = NULL,
    disease_state  = "Healthy adult male volunteers (no diagnosed condition).",
    dose_range     = paste(
      "Single ascending oral suspension doses, two cohorts of 12 volunteers",
      "each. Cohort 1 (fasted): 3, 10, 30, 100 mg in dose escalation, one",
      "placebo. Cohort 2 (fasted): 250, 500, 1000, 1300 mg in dose escalation,",
      "one placebo; the cohort then received a repeat 500 mg dose in the fed",
      "state. 7-day washout between treatments.",
      sep = " "
    ),
    regions        = "Singapore (Singapore General Hospital).",
    notes          = paste(
      "Body mass index approximately 18 - 30 kg/m^2 per inclusion criteria.",
      "Plasma PK samples at predose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16,",
      "24, 36 h post dose (and 48 h for PD). Supine pulse rate was recorded",
      "predose and at 1, 2, 4, 6, 8, 12, 24, 48 h. Dosing was at approximately",
      "08:00. PK-data set: 969 plasma concentrations contributed to the popPK",
      "fit. Dose escalation followed exposure- and safety-driven stopping",
      "rules; the PD data set was correspondingly truncated.",
      "NONMEM VI with FOCE INTER; nonparametric bootstrap (500 iterations)",
      "provided SE and 95% CI on parameter estimates. PK and PD fit",
      "sequentially (individual Bayesian post hoc PK as input to PD).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # PK structural parameters -- Langdon 2010 Table 2 (human final NONMEM
    # estimate; SE and bootstrap 95% CI as reported).
    # =====================================================================
    lfdepot <- fixed(log(1)); label("Oral bioavailability (FIXED at 1; absolute F unknown)")  # Table 2: F = 1 (FIXED)
    lka     <- log(0.599);    label("First-order oral absorption rate constant ka (1/h)")     # Table 2: ka = 0.599 1/h (SE 0.0101; bootstrap 95% CI 0.582 - 0.616)
    lcl     <- log(36.7);     label("Apparent clearance CL/F (L/h)")                          # Table 2: CL/F = 36.7 L/h (SE 2.55; bootstrap 95% CI 32.8 - 40.9)
    lvc     <- log(18.4);     label("Apparent central volume Vc/F (L)")                       # Table 2: Vc/F = 18.4 L (SE 3.08; bootstrap 95% CI 14.2 - 24.3)
    lvp     <- log(6.88);     label("Apparent peripheral volume Vp/F (L)")                    # Table 2: Vp/F = 6.88 L (SE 1.03; bootstrap 95% CI 5.51 - 8.89)
    lq      <- log(0.704);    label("Apparent inter-compartmental clearance Q/F (L/h)")       # Table 2: Q = 0.704 L/h (SE 0.117; bootstrap 95% CI 0.546 - 0.945)

    # =====================================================================
    # PD structural parameters -- Langdon 2010 Table 2 (human final NONMEM
    # estimate). Equation 1 (placebo pulse rate cosine) and Equation 2
    # (linear free-drug effect) on page 338.
    # =====================================================================
    lbase  <- log(61.5); label("Baseline supine pulse rate (bpm)")                            # Table 2: BASE = 61.5 bpm (SE 1.36; bootstrap 95% CI 59.3 - 63.6)
    lamp   <- log(4.01); label("Circadian-rhythm amplitude on pulse rate (bpm)")              # Table 2: AMP = 4.01 bpm (SE 0.416; bootstrap 95% CI 3.33 - 4.82)
    lpeak  <- log(14.3); label("Time of circadian peak post morning dose (h)")                # Table 2: PEAK = 14.3 h (SE 0.55; bootstrap 95% CI 13.3 - 15.2)
    lslope <- log(0.76); label("Linear slope of pulse rate on free plasma concentration (bpm per uM free)")  # Table 2: SLOPE = 0.76 bpm/uM-free (SE 0.16; bootstrap 95% CI 0.54 - 1.14)

    # Reference unbound fraction. FIXED via back-calculation from the
    # canine Introduction text ('unbound Cmax = 56 uM' at 20 mg/kg oral
    # divided by the typical-value total Cmax = 87.8 uM computed from the
    # Langdon 2010 Table 1 canine PK parameters). The paper does not
    # quantify fu directly for either species; the operator-chosen value
    # fu = 0.64 is applied to both dog and human models for consistency
    # (no species-difference in fu is reported). See vignette Errata.
    lfu <- fixed(log(0.64)); label("Plasma unbound fraction (FIXED; back-calculated)")        # vignette Errata (not in source); back-calculation from paper Introduction

    # =====================================================================
    # IIV. Per Methods paragraph 'Human pharmacokinetic-pharmacodynamic
    # analysis' ('IIV in the pharmacokinetic parameters was modelled using
    # multiplicative exponential random effects'), every IIV term below is
    # encoded as a log-normal random effect on the log-transformed
    # structural parameter; the reported variances are omega^2 directly.
    #
    # PK IIV. Per Methods + Table 2: 'IIV was included on the following
    # structural parameters: relative bioavailability (F), V/F, CL/F'.
    # IIV on F is encoded around the FIXED population value of 1.
    # =====================================================================
    etalfdepot ~ 0.075   # Table 2: Variance for F = 0.075 (SE 0.030; bootstrap 95% CI 0.017 - 0.123); CV approx 28.0%
    etalcl     ~ 0.039   # Table 2: Variance for CL = 0.039 (point estimate); SE 0.028; bootstrap 95% CI 8.73e-9 - 9.35e-4 (CI is degenerate, see vignette Errata); CV approx 19.9%
    etalvc     ~ 0.332   # Table 2: Variance for Vc = 0.332 (SE 0.116; bootstrap 95% CI 0.164 - 0.505); CV approx 62.5%

    # PD IIV (Methods + Table 2: 'IIV was included on the following
    # structural parameters: baseline HR and peak time of the circadian
    # effect'). The text and Table 2 agree.
    etalbase ~ 0.0094  # Table 2: Variance for BASE = 0.0094 (SE 0.0026; bootstrap 95% CI 0.0050 - 0.0133); CV approx 9.7% (matches text 'approximately 10%')
    etalpeak ~ 0.0206  # Table 2: Variance for PEAK = 0.0206 (SE 0.0055; bootstrap 95% CI 0.0099 - 0.0299); CV approx 14.5%

    # =====================================================================
    # Residual error.
    # PK: 'additive error on the log-transformed concentrations' per
    # Methods. Encoded as a log-normal (lnorm) residual; the reported
    # value is the log-scale SD. The Table 2 unit annotation 'ng/mL' is an
    # editorial inconsistency -- a log-scale SD is dimensionless.
    # PD: additive on supine pulse rate per Methods + Table 2.
    # =====================================================================
    expSd    <- 0.658; label("PK residual error -- log-scale additive SD (dimensionless)")    # Table 2: PK residual error = 0.658 (SE 0.0615; bootstrap 95% CI 0.546 - 0.759); Methods 'additive error on the log-transformed concentrations'
    addSd_HR <- 4.77;  label("Additive PD residual error on supine pulse rate (bpm)")         # Table 2: PD additive residual error = 4.77 bpm (SE 0.309; bootstrap 95% CI 4.28 - 5.23)
  })

  model({
    # Constant: PF-00821385 molecular weight (paper Introduction,
    # 'PF-00821385, molecular weight = 440.49 Da'). Used to convert plasma
    # total ng/mL into micromolar for the free-drug PD argument.
    mw <- 440.49  # g/mol

    # Individual structural parameters.
    fdepot <- exp(lfdepot + etalfdepot)
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    fu     <- exp(lfu)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual PD parameters.
    base_i  <- exp(lbase + etalbase)
    amp_i   <- exp(lamp)
    peak_i  <- exp(lpeak + etalpeak)
    slope_i <- exp(lslope)

    # PK ODE -- 2-compartment disposition with first-order absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability anchor (FIXED at 1 population; eta carries IIV).
    f(depot) <- fdepot

    # Plasma concentration: central is in mg, vc in L so central/vc is in
    # mg/L; multiply by 1000 to express Cc in ng/mL (paper unit).
    Cc <- 1000 * central / vc

    # Free plasma concentration in micromolar (the argument of the SLOPE
    # in Equation 2 of Langdon 2010). Cc_free_uM = fu * Cc_ng_per_mL / MW
    # (g/mol).
    Cc_free_uM <- fu * Cc / mw

    # Equation 1: placebo pulse-rate cosine circadian rhythm.
    HR_placebo <- base_i + amp_i * cos(2 * pi / 24 * (t - peak_i))

    # Equation 2: HR = HR_placebo + SLOPE * Cc_free.
    HR <- HR_placebo + slope_i * Cc_free_uM

    # Multi-output residual error.
    Cc ~ lnorm(expSd)
    HR ~ add(addSd_HR)
  })
}
