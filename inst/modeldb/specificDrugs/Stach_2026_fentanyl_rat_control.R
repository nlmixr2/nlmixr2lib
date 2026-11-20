Stach_2026_fentanyl_rat_control <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment pharmacokinetic model of intravenous",
    "fentanyl in anaesthetised (isoflurane) normovolaemic male Wistar rats, the",
    "control group of Stach 2026, after a 10 ug/kg dose infused over 5 min.",
    "Stach 2026 fitted each rat individually in Phoenix WinNonlin 8.2 (two-",
    "compartment model, zero-order 5-min input, first-order elimination from",
    "the central compartment, 1/Yhat^2 weighting) and reported only the mean or",
    "median of the six individual estimates, so no between-subject variability",
    "or residual-error magnitude is available: every parameter is fixed at the",
    "published group summary and the proportional residual SD is fixed at zero.",
    "Clearance and both volumes are the Table 3 group means; the",
    "intercompartmental clearance, which Stach 2026 estimated (as k12/k21) but",
    "did not tabulate, is back-solved from the Table 3 median elimination",
    "half-life. All states are per kg body weight (amounts in ng/kg, volumes in",
    "L/kg), so concentration is in pg/mL. Norfentanyl was analysed only by",
    "non-compartmental methods and is not modelled. The companion",
    "hemorrhagic-shock group is",
    "modellib('Stach_2026_fentanyl_rat_hemorrhagicShock')."
  )
  reference <- paste(
    "Stach P, Skowron K, Rojek S, Cios A, Wesolowska A, Huestis MA, Gil K",
    "(2026). Pharmacokinetics of Fentanyl and Norfentanyl in a Rat Model of",
    "Fixed-Volume Hemorrhagic Shock. Eur J Drug Metab Pharmacokinet",
    "51(5):637-647. doi:10.1007/s13318-026-01023-z"
  )
  vignette <- "Stach_2026_fentanyl_rat"
  units <- list(time = "min", dosing = "ng/kg", concentration = "pg/mL")

  compartmentData <- list(
    central = list(analyte = "fentanyl", units = "ng/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fentanyl", units = "ng/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 6,
    n_studies = 1,
    age_range = "18-20 weeks",
    weight_range = "441.3 +/- 26.4 g (mean +/- SD); median 433.0 g (Q1-Q3 425.3-453.5 g)",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy, normovolaemic control animals under isoflurane/oxygen",
      "anaesthesia with jugular-vein and femoral-artery catheters; no blood",
      "withdrawal before dosing."
    ),
    dose_range = "fentanyl 10 ug/kg as a single 5-min intravenous infusion",
    regions = "Poland",
    notes = paste(
      "Stach 2026 Methods (Experimental Design): twelve male Wistar rats aged",
      "18-20 weeks (average 450 g) were randomised to a control group (n = 6)",
      "and a fixed-volume hemorrhagic-shock group (n = 6); this file is the",
      "control group. Core temperature was held at 37.0 +/- 0.5 C. Fentanyl",
      "was infused through the jugular-vein catheter over 5 min and arterial",
      "plasma sampled at 5, 10, 15, 20, 30, 40, 50 and 60 min; 100 uL of saline",
      "replaced each sample. Body weight is from the Results (General Animal",
      "Parameters); end-of-study biochemistry (albumin 33.7 +/- 1.5 g/L,",
      "creatinine 27.1 +/- 6.8 umol/L) is from Table 1. No race/ethnicity data",
      "apply."
    )
  )

  ini({
    # Table 3 (compartmental analysis, control group, n = 6): CL 5.9 +/- 1.2
    # L/h/kg, Vc 0.9 +/- 0.4 L/kg, Vp 3.5 +/- 0.8 L/kg. The model runs in
    # minutes (the Table 3 rate constants and half-lives are per minute), so
    # CL = 5.9 / 60 = 0.098333 L/min/kg. Vc + Vp = 4.4 L/kg reproduces the
    # Table 3 Vss of 4.4 +/- 0.9 L/kg.
    lcl <- fixed(log(5.9 / 60)); label("Clearance CL (L/min/kg)") # Table 3, control 'CL (L/h/kg)' = 5.9, divided by 60
    lvc <- fixed(log(0.9)); label("Central volume of distribution Vc (L/kg)") # Table 3, control 'Vc (L/kg)' = 0.9
    lvp <- fixed(log(3.5)); label("Peripheral volume of distribution Vp (L/kg)") # Table 3, control 'Vp (L/kg)' = 3.5

    # Intercompartmental clearance. Stach 2026 estimated k12 and k21 (Methods,
    # Pharmacokinetic Analysis) but Table 3 does not print them. With CL, Vc
    # and Vp fixed above, Q is the single remaining unknown; it is solved so
    # the beta (terminal) half-life of the two-compartment system equals the
    # Table 3 control median t1/2beta = 44.1 min, giving Q = 0.1587 L/min/kg.
    # The implied distribution half-life is then t1/2alpha = 2.20 min against
    # the independently reported 1.9 +/- 0.7 min (Table 3).
    lq <- fixed(log(0.1587)); label("Intercompartmental clearance Q (L/min/kg)") # back-solved from Table 3 control t1/2beta = 44.1 min

    # Section 'Compartmental Pharmacokinetics of Fentanyl': WinNonlin fits used
    # a 1/Yhat^2 weighting scheme, i.e. a proportional error structure; its
    # magnitude is not reported.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported)") # 1/Yhat^2 weighting, magnitude not reported
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
