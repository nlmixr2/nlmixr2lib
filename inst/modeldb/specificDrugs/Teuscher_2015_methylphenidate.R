Teuscher_2015_methylphenidate <- function() {
  description <- "Pediatric population PK model for methylphenidate hydrochloride extended-release multilayer beads (MPH-MLR, Aptensio XR) after a single oral dose, parameterized as a two-input, one-compartment, first-order-elimination structure: a fast-release (IR) depot delivers a fraction F1 of dose with first-order absorption rate Ka1, a slow-release (ER) depot delivers the remaining 1 - F1 with first-order rate Ka2 after an absorption lag tlag, and the central compartment eliminates linearly via clearance CL and apparent volume V. Body weight enters CL via a power covariate CL = CL_TV * WT^theta (Eq 4). Between-individual variability is retained on CL and V; IIV on Ka1, Ka2, F1, and tlag was not in the final pediatric fit (Table 1). The companion exposure-response analysis maps simulated Cmax to change-from-baseline ADHD-RS-IV total score via the Emax model E = Emax * Cmax / (EC50 + Cmax) with Emax = -34.96 and EC50 = 5.77 ng/mL (Table 2); the PD step lives outside the ODE system because the published mapping uses a per-period Cmax, not the instantaneous central concentration. The vignette reproduces the full PK simulation, NCA, and Cmax-to-ADHD-RS-IV exposure-response."
  reference <- paste(
    "Teuscher NS, Adjei A, Findling RL, Greenhill LL, Kupper RJ, Wigal S.",
    "Population pharmacokinetics of methylphenidate hydrochloride extended-release",
    "multiple-layer beads in pediatric subjects with attention deficit hyperactivity disorder.",
    "Drug Des Devel Ther. 2015;9:2767-2775.",
    "doi:10.2147/DDDT.S83234.",
    sep = " "
  )
  vignette <- "Teuscher_2015_methylphenidate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight; used in the unnormalized power model CL = CL_TV * WT^theta (Eq 4). No reference weight in the published parameterization.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 17L,
    n_studies      = 1L,
    n_observations = "154 plasma MPH concentrations supplied the pediatric population PK fit (Results).",
    age_range      = "6-12 years",
    weight_range   = "Pediatric ADHD cohort; simulations in the source paper spanned 70-150 lb (~32-68 kg).",
    sex_female_pct = NA_real_,
    disease_state  = "Children with attention deficit hyperactivity disorder (ADHD).",
    dose_range     = "Single oral doses of MPH-MLR (Aptensio XR) 10, 15, 20, 30, or 40 mg strengths in a single-center, randomized, open-label, two-way crossover study comparing MPH-MLR with IR methylphenidate. Only the MPH-MLR arm was used for the pediatric PK fit. Sparse plasma sampling at predose and 1, 2, 3, 4, 5, 6, 8, 10, 12, and 24 h postdose.",
    regions        = "United States.",
    notes          = "Pediatric PK structure inherited from a separately developed adult PK model (Methods, Population PK modeling); the pediatric fit added body weight as a covariate on CL (Eq 4). Source dataset described in Methods, Data sources."
  )

  ini({
    # Structural parameters - Table 1, final pediatric population PK model.
    # CL is reported as the typical value CL_TV in the unnormalized power
    # parameterization CL_i = CL_TV * WT^theta * exp(eta_CL) (Eq 4); the
    # apparent units of CL_TV are L/h per (kg)^theta, but the table labels
    # CL in L/h to match the resulting whole-subject CL.
    lcl     <- log(1.3);      label("Typical-value apparent clearance CL_TV (L/h) in CL = CL_TV * WT^theta") # Teuscher 2015 Table 1: CL = 1.3 (RSE 68.8 %)
    lvc     <- log(64.7);     label("Apparent central volume V/F (L)")                                        # Teuscher 2015 Table 1: V = 64.7 L (RSE 25.8 %)
    lka1    <- log(0.25);     label("First-order absorption rate from the IR depot Ka1 (1/h)")                # Teuscher 2015 Table 1: Ka1 = 0.25 1/h (RSE 26.6 %)
    lka2    <- log(0.16);     label("First-order absorption rate from the ER depot Ka2 (1/h)")                # Teuscher 2015 Table 1: Ka2 = 0.16 1/h (RSE 31.8 %)
    logitf1 <- qlogis(0.65);  label("Fraction of total dose released through the IR depot F1 (logit-scale)")  # Teuscher 2015 Table 1: F1 = 0.65 (RSE 10.5 %); ER fraction = 1 - F1
    ltlag   <- log(5.75);     label("Absorption lag time tlag on the ER depot (h)")                           # Teuscher 2015 Table 1: tlag = 5.75 h (RSE 3.6 %)
    e_wt_cl <- 1.53;          label("Body weight exponent on CL (unitless) in CL = CL_TV * WT^theta")         # Teuscher 2015 Table 1, Eq 4: theta = 1.53 (RSE 12.9 %)

    # Inter-individual variability - Table 1. The paper reports omega values
    # of 0.057 (CL) and 0.096 (V) as the lognormal variances on the log-scale
    # individual parameters, paired with the lognormal IIV model in Eq 1.
    # Shrinkage 18 % (CL) and 29 % (V) is recorded in the table footnote;
    # see vignette Assumptions and deviations for the shrinkage commentary.
    etalcl  ~ 0.057   # Teuscher 2015 Table 1: omega_CL = 0.057 (shrinkage 18 %)
    etalvc  ~ 0.096   # Teuscher 2015 Table 1: omega_V  = 0.096 (shrinkage 29 %)

    # Residual error - Table 1. The paper selected the proportional residual
    # error model `C_ij = C_hat_ij * (1 + epsilon_2ij)` with epsilon_2ij ~
    # N(0, omega_1^2) (Eq 3). The reported point estimate epsilon = 1.94
    # (RSE 6.2 %) is encoded literally as the SD of the proportional error
    # term. The value is anomalously large relative to the goodness-of-fit
    # plots; see the vignette Assumptions and deviations section for the
    # flagged caveat and recommended sensitivity check.
    propSd  <- 1.94;          label("Proportional residual error SD (fraction)")  # Teuscher 2015 Table 1: epsilon = 1.94 (RSE 6.2 %)
  })

  model({
    # Individual structural parameters. Body weight enters CL via the
    # unnormalized power model from Eq 4 (no reference weight scaling).
    cl   <- exp(lcl + etalcl) * WT^e_wt_cl
    vc   <- exp(lvc + etalvc)
    ka1  <- exp(lka1)
    ka2  <- exp(lka2)
    f1   <- expit(logitf1)
    tlag <- exp(ltlag)

    kel  <- cl / vc

    # Two-input, one-compartment, first-order-elimination ODE system
    # (Figure 1): the IR depot (`depot`) delivers fraction F1 of dose with
    # rate Ka1, the ER depot (`depot2`) delivers 1 - F1 with rate Ka2 after
    # absorption lag tlag, and the central compartment eliminates linearly
    # with first-order rate kel = CL/V.
    d/dt(depot)   <- -ka1 * depot
    d/dt(depot2)  <- -ka2 * depot2
    d/dt(central) <-  ka1 * depot + ka2 * depot2 - kel * central

    # Dose splitting and lag time. A single oral dose event records targets
    # depot and depot2 in turn (cmt = depot, depot2) with the same amt;
    # the f() multipliers split the dose into the IR (F1) and ER (1 - F1)
    # fractions, and lag() imposes the ER absorption lag.
    f(depot)    <- f1
    f(depot2)   <- 1 - f1
    lag(depot2) <- tlag

    # Plasma observation. central is in mg and vc in L, so central/vc is
    # in mg/L; the 1000 multiplier converts to ng/mL (the unit declared in
    # `units$concentration` and the unit the paper reports throughout).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
