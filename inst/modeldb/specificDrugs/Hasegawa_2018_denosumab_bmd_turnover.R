Hasegawa_2018_denosumab_bmd_turnover <- function() {
  description <- paste(
    "Empirical turnover (indirect-response) model of lumbar spine BMD after subcutaneous denosumab:",
    "free serum denosumab inhibits the first-order BMD loss rate through an Imax function, with",
    "R_in = k_out * 100 so BMD starts at 100% of baseline. Free denosumab comes from the same fixed",
    "two-compartment quasi-steady-state TMDD model used for the reduced QSP model",
    "(Hasegawa_2018_denosumab_qsp). Fitted by Hasegawa and Duffull 2018 to the first 12 months of",
    "mean BMD data in postmenopausal women as an empirical comparator; it does not capture the",
    "continued BMD gain beyond 12 months."
  )
  reference <- paste(
    "Hasegawa C, Duffull SB. Automated Scale Reduction of Nonlinear QSP Models With an Illustrative",
    "Application to a Bone Biology System. CPT Pharmacometrics Syst Pharmacol. 2018;7(9):562-572.",
    "doi:10.1002/psp4.12324"
  )
  vignette <- "Hasegawa_2018_denosumab_bone_scale_reduction"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "denosumab", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "denosumab (total, free + bound)", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "denosumab", units = "ug", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "total RANKL (free + denosumab-bound)",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    BMD_LS = list(
      analyte = "lumbar spine bone mineral density",
      units = "% of baseline (100 at time 0)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 190L,
    n_studies = 1L,
    sex_female_pct = 100,
    weight_median = "66 kg (fixed in the denosumab absorption rate)",
    disease_state = "Postmenopausal women with low bone mineral density",
    dose_range = "Denosumab 14, 60, 100 or 210 mg subcutaneously every 6 months",
    regions = "United States (phase 2 denosumab dose-ranging study)",
    notes = paste(
      "Fitted to digitized arm-mean lumbar-spine BMD (percent change from baseline) for months 0-12 of",
      "the every-6-month arms (54, 47, 42 and 47 women at 14, 60, 100 and 210 mg; Hasegawa 2018",
      "Results and Table S2). See Hasegawa_2018_denosumab_qsp for the shared PK provenance."
    )
  )

  ini({
    # ---- Denosumab PK (fixed; identical to Hasegawa_2018_denosumab_qsp) ----
    lcl <- fixed(log(0.00306))
    label("Linear clearance of denosumab (L/h)") # nonlinear_lumped.mod THETA(22) = 3.06 FIX (mL/h)
    lvc <- fixed(log(2.49))
    label("Central volume of distribution (L)") # nonlinear_lumped.mod THETA(23) = 2490 FIX (mL)
    lq <- fixed(log(0.0379))
    label("Intercompartmental clearance (L/h)") # nonlinear_lumped.mod THETA(24) = 37.9 FIX (mL/h)
    lvp <- fixed(log(1.36))
    label("Peripheral volume of distribution (L)") # nonlinear_lumped.mod THETA(25) = 1360 FIX (mL)
    lka <- fixed(log(0.212 / 24))
    label("Subcutaneous absorption rate constant at 71.5 kg (1/h)") # nonlinear_lumped.mod THETA(26) = 0.212 FIX (1/day)
    e_wt_ka <- fixed(-0.577)
    label("Power exponent of (WT/71.5) on ka, evaluated at a body weight of 66 kg (unitless)") # nonlinear_lumped.mod: KA = TVKA*(66/71.5)**(-0.577)
    lfdepot <- fixed(log(0.638))
    label("Subcutaneous bioavailability (fraction)") # nonlinear_lumped.mod THETA(27) = 0.638 FIX
    lkint <- fixed(log(0.00795))
    label("Internalization rate constant of the denosumab-RANKL complex (1/h)") # Table 1 kint = 0.00795
    lkss <- fixed(log(1380))
    label("Quasi-steady-state constant of the denosumab PK model (ng/mL)") # MATLAB Model_parameter_values0.m Kss = 138*10 (see vignette)
    lrbase_target <- fixed(log(614))
    label("Baseline total RANKL in the PK model (ng/mL)") # MATLAB Model_parameter_values0.m R0 = 614
    lkdeg <- fixed(log(0.00148))
    label("Degradation rate constant of free RANKL in the PK model (1/h)") # MATLAB Model_parameter_values0.m kdeg = 0.00148

    # ---- Turnover model (Eq 12; Table 2) ----
    limax <- log(0.0577)
    label("Maximal fractional inhibition of BMD loss by denosumab, Imax (fraction)") # Table 2 Imax = 0.0577 (15.2% RSE)
    lkout <- log(0.000336)
    label("First-order loss rate constant of BMD, k_out (1/h)") # Table 2 k_out = 0.000336 /h (26.5% RSE)
    lec50 <- log(52.5)
    label("Free denosumab concentration giving half the maximal inhibition, C50 (ng/mL)") # Table 2 C50 = 52.5 ng/mL (67.6% RSE)

    # ---- Residual error ----
    addSd <- 0.6220932
    label("Additive residual error on BMD (% of baseline)") # Table 2 sigma^2 = 0.387 (26.4% RSE); sqrt(0.387)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka) * (66 / 71.5)^e_wt_ka
    kint <- exp(lkint)
    kss <- exp(lkss)
    rbase_target <- exp(lrbase_target)
    kdeg <- exp(lkdeg)
    ksyn <- rbase_target * kdeg

    ctot <- central / vc
    qssdisc <- ctot - total_target - kss
    Cc <- 0.5 * (qssdisc + sqrt(qssdisc * qssdisc + 4 * kss * ctot))
    cltot <- cl + kint * vc * total_target / (kss + Cc)

    d/dt(depot) <- -ka * depot
    f(depot) <- exp(lfdepot)
    d/dt(central) <- ka * depot - cltot * Cc - q * (Cc - peripheral1 / vp)
    d/dt(peripheral1) <- q * (Cc - peripheral1 / vp)
    d/dt(total_target) <- ksyn - (kdeg + (kint - kdeg) * Cc / (kss + Cc)) * total_target
    total_target(0) <- rbase_target

    # Eq 12 with R_in = k_out * 100 (Methods: 'Rin is derived as kout times 100').
    imax <- exp(limax)
    kout <- exp(lkout)
    ec50 <- exp(lec50)
    kin <- kout * 100
    d/dt(BMD_LS) <- kin - kout * (1 - imax * Cc / (ec50 + Cc)) * BMD_LS
    BMD_LS(0) <- 100

    bmd_pcfb <- BMD_LS - 100
    BMD_LS ~ add(addSd)
  })
}
