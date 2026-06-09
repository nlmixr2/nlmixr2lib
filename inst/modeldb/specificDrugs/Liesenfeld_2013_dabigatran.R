Liesenfeld_2013_dabigatran <- function() {
  description <- "Two-compartment population PK model for oral dabigatran (after dabigatran etexilate prodrug) in seven end-stage renal disease (ESRD) subjects undergoing intermittent hemodialysis, with first-order absorption, absorption lag, an apparent total body clearance (renal + non-renal), and an apparent dialysis clearance described by the Michaels equation as a function of blood and dialysate flow rates and a hemodialyzer mass transfer-area coefficient (Liesenfeld 2013)."
  reference <- "Liesenfeld KH, Staab A, Haertter S, Formella S, Clemens A, Lehr T. Pharmacometric Characterization of Dabigatran Hemodialysis. Clin Pharmacokinet. 2013;52(6):453-462. doi:10.1007/s40262-013-0049-6"
  vignette <- "Liesenfeld_2013_dabigatran"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HEMODIALYSIS = list(
      description        = "Hemodialysis-active indicator (1 during a dialysis session, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = "Time-varying within subject. Gates the Michaels-equation apparent dialysis clearance: CLdialysis/F is added to CL/F only when HEMODIALYSIS = 1. The source paper reports four-hour intermittent hemodialysis sessions on Days 1, 3, and 5 of each study period (Methods, Study Design). For non-ESRD patients with no dialysis, set HEMODIALYSIS = 0 throughout. The source paper's data column was named `DIAL`; renamed to the canonical `HEMODIALYSIS` per inst/references/covariate-columns.md.",
      source_name        = "DIAL"
    ),
    BFR = list(
      description        = "Blood flow rate during the active hemodialysis session",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. Enters the Michaels equation together with DFR and KoA to compute the apparent dialysis clearance. Values investigated in the source study were 200, 300, and 400 mL/min (Methods, Study Design; Table 1). Should be 0 (or any sentinel) when HEMODIALYSIS = 0 since the Michaels term is gated off.",
      source_name        = "BFR"
    ),
    DFR = list(
      description        = "Dialysate flow rate during the active hemodialysis session",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. Enters the Michaels equation together with BFR and KoA. The source study fixed DFR at 700 mL/min throughout (Methods, Study Design); the simulation scenarios also examined 500 mL/min (Methods, Simulations). Should be 0 (or any sentinel) when HEMODIALYSIS = 0.",
      source_name        = "DFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 7L,
    n_studies      = 1L,
    age_range      = "27-53 years (mean 38.3)",
    weight_range   = "60-87 kg (mean 74.0)",
    sex_female_pct = 0,
    race_ethnicity = "All seven subjects were white males (Results, Population Pharmacokinetic Model).",
    disease_state  = "Dialysis-dependent end-stage renal disease (ESRD) without atrial fibrillation. Source paper notes that renal elimination of dabigatran in this population is negligible (< 0.04 L/h); the estimated total body clearance therefore reflects predominantly non-renal clearance.",
    dose_range     = "Three oral doses of dabigatran etexilate per study period, separated by 21 h: 150 mg (post-dialysis on Day 1), 110 mg (Day 2), and 75 mg (Day 3, 8 h before dialysis). Repeated across two periods with different Day-3 blood flow rates (200 and 400 mL/min). Hemodialysis on Days 1, 3, and 5 of each period; DFR fixed at 700 mL/min; high-flux Polyflux PF-210H filter (Gambro).",
    regions        = "Single-centre phase I dialysis study (Khadzhynov et al. 2013, ref [12] of the source paper).",
    n_observations = 308L,
    notes          = "Baseline demographics from Liesenfeld 2013 Results, Population Pharmacokinetic Model. Geometric mean trough concentrations after the second dose were 140 ng/mL (CV 54.2%) in period 1 and 128 ng/mL (CV 44.5%) in period 2. The structural PK parameters reported here characterise dabigatran disposition in ESRD; the paper explicitly used parameter estimates from the RE-LY trial (Reilly et al. 2013) when simulating typical AF patients, so the ESRD-cohort parameters in ini() are not appropriate for simulating non-ESRD AF dosing scenarios -- they describe the cohort that informed the dialysis-clearance component."
  )

  ini({
    # Structural PK parameters for the ESRD population (Liesenfeld 2013 Table 2).
    # All apparent (CL/F, V/F) because relative bioavailability F is fixed at 1.
    lcl     <- log(12.4);  label("Apparent total body clearance in ESRD subjects (CL/F, L/h)")            # Liesenfeld 2013 Table 2
    lvc     <- log(531);   label("Apparent central volume of distribution (V2/F, L)")                    # Liesenfeld 2013 Table 2
    lq      <- log(152);   label("Apparent inter-compartmental clearance (Q/F, L/h)")                    # Liesenfeld 2013 Table 2
    lvp     <- log(499);   label("Apparent peripheral volume of distribution (V3/F, L)")                 # Liesenfeld 2013 Table 2
    lka     <- log(0.821); label("First-order absorption rate constant (ka, 1/h)")                       # Liesenfeld 2013 Table 2
    ltlag   <- log(1.67);  label("Absorption lag time for the fed condition (ALAG, h)")                  # Liesenfeld 2013 Table 2
    lkoa    <- log(313);   label("Hemodialyzer mass transfer-area coefficient (KoA, mL/min)")            # Liesenfeld 2013 Table 2

    # Relative bioavailability fixed at 1.0 (Liesenfeld 2013 Table 2 footnote b).
    lfdepot <- fixed(log(1.00)); label("Relative bioavailability (F, fraction)")                          # Liesenfeld 2013 Table 2 (fixed)

    # Inter-individual variability. Liesenfeld 2013 Table 2 reports CV%; the
    # exponential random-effect model implies omega^2 = log(CV^2 + 1).
    #   CL/F   : 40.4% CV -> omega^2 = log(0.404^2 + 1) = 0.15113
    #   V2/F   : 14.3% CV -> omega^2 = log(0.143^2 + 1) = 0.02024
    etalcl ~ 0.15113   # Liesenfeld 2013 Table 2 (IIV CL/F, 40.4% CV)
    etalvc ~ 0.02024   # Liesenfeld 2013 Table 2 (IIV V2/F, 14.3% CV)

    # The source paper additionally reports inter-occasion variability (IOV) on
    # ka (64.0% CV) and on F (48.0% CV) with one occasion per 21-h dosing
    # interval. nlmixr2lib encodes these as conventional IIV on ka and on
    # bioavailability so that simulated typical-population profiles carry the
    # same total stochastic spread; the "per-occasion" structure of the source
    # is not preserved. See vignette Assumptions and deviations.
    #   ka : 64.0% CV -> omega^2 = log(0.64^2 + 1) = 0.34336
    #   F  : 48.0% CV -> omega^2 = log(0.48^2 + 1) = 0.20733
    etalka      ~ 0.34336   # Liesenfeld 2013 Table 2 (IOV ka, 64.0% CV; recast as IIV)
    etalfdepot  ~ 0.20733   # Liesenfeld 2013 Table 2 (IOV F, 48.0% CV; recast as IIV)

    # Residual error: proportional (Liesenfeld 2013 Table 2, PRV 8.5% CV).
    propSd <- 0.085; label("Proportional residual error (fraction)")  # Liesenfeld 2013 Table 2 (PRV)
  })
  model({
    # Individual PK parameters at the ESRD-cohort typical level (no allometric
    # or covariate scaling applied -- the source paper retained no continuous
    # covariates after testing serum creatinine).
    ka       <- exp(lka  + etalka)
    cl       <- exp(lcl  + etalcl)
    vc       <- exp(lvc  + etalvc)
    vp       <- exp(lvp)
    q        <- exp(lq)
    alag_d   <- exp(ltlag)
    fdepot   <- exp(lfdepot + etalfdepot)
    koa      <- exp(lkoa)

    # Apparent dialysis clearance via the Michaels equation (Liesenfeld 2013
    # Equation 1; mass-transfer / dialyzer-clearance form attributed to
    # Michaels [17] in the paper). The numerical output is interpreted as the
    # apparent dialysis clearance CLdialysis/F in L/h when BFR, DFR, and KoA
    # are all supplied in mL/min -- the equivalence holds because the
    # mL/min -> L/h conversion factor (60/1000) cancels the absolute
    # bioavailability F_abs = 0.06 used by the paper for the apparent /
    # actual reporting bridge (60/1000/0.06 = 1).
    #
    # bfr_safe / dfr_safe guard against the BFR = 0 / DFR = 0 condition that
    # holds in the interdialytic period: when HEMODIALYSIS = 0 the safe
    # values fall back to 1 mL/min so the rational expression evaluates
    # without divide-by-zero, then the entire term is gated to zero by
    # the leading HEMODIALYSIS multiplier.
    bfr_safe         <- BFR * HEMODIALYSIS + (1 - HEMODIALYSIS)
    dfr_safe         <- DFR * HEMODIALYSIS + (1 - HEMODIALYSIS)
    michaels_arg     <- -koa * (dfr_safe - bfr_safe) / (bfr_safe * dfr_safe)
    michaels_exp     <- exp(michaels_arg)
    cl_dialysis_raw  <- bfr_safe * (1 - michaels_exp) / (1 - (bfr_safe / dfr_safe) * michaels_exp)
    cl_dialysis      <- HEMODIALYSIS * cl_dialysis_raw

    # Total apparent clearance during the active session is CL/F + CLdialysis/F
    # (Liesenfeld 2013 Results, Table 1 reports the sum as CLtotal/F).
    cl_total <- cl + cl_dialysis

    kel <- cl_total / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- alag_d

    # Concentration: dose in mg, vc in L -> central / vc in mg/L = ug/mL.
    # The source paper reports dabigatran plasma concentrations in ng/mL, so
    # the model output is rescaled by 1000 to ng/mL to match the source
    # tables and figures.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
