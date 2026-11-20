Stach_2026_fentanyl_rat_hemorrhagicShock <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment pharmacokinetic model of intravenous",
    "fentanyl in anaesthetised (isoflurane) male Wistar rats in fixed-volume",
    "hemorrhagic shock (30 percent of estimated blood volume withdrawn over",
    "15 min before dosing), the HS group of Stach 2026, after a 10 ug/kg dose",
    "infused over 5 min.",
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
    "normovolaemic control group is",
    "modellib('Stach_2026_fentanyl_rat_control')."
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
    weight_range = "471.8 +/- 41.1 g (mean +/- SD); median 482.0 g (Q1-Q3 446.5-497.3 g)",
    sex_female_pct = 0,
    disease_state = paste(
      "Acute fixed-volume hemorrhagic shock under isoflurane/oxygen",
      "anaesthesia: 30 percent of the estimated blood volume (EBV, mL =",
      "0.06 x body weight in g + 0.77) withdrawn through the femoral-artery",
      "catheter within 15 min, no fluid resuscitation; fentanyl given as the",
      "next step after the withdrawal."
    ),
    dose_range = "fentanyl 10 ug/kg as a single 5-min intravenous infusion",
    regions = "Poland",
    notes = paste(
      "Stach 2026 Methods (Experimental Design): twelve male Wistar rats aged",
      "18-20 weeks (average 450 g) were randomised to a control group (n = 6)",
      "and a fixed-volume hemorrhagic-shock group (n = 6); this file is the",
      "hemorrhagic-shock (HS) group. Core temperature was held at 37.0 +/- 0.5 C. Fentanyl",
      "was infused through the jugular-vein catheter over 5 min and arterial",
      "plasma sampled at 5, 10, 15, 20, 30, 40, 50 and 60 min; 100 uL of saline",
      "replaced each sample. Body weight is from the Results (General Animal",
      "Parameters); end-of-study biochemistry (albumin 27.8 +/- 1.6 g/L,",
      "creatinine 61.4 +/- 12.1 umol/L, lactate median 2.3 mmol/L) is from",
      "Table 1. No race/ethnicity data",
      "apply."
    )
  )

  ini({
    # Table 3 (compartmental analysis, HS group, n = 6): CL 2.6 +/- 0.5
    # L/h/kg, Vc 0.6 +/- 0.4 L/kg, Vp 3.5 +/- 1.6 L/kg. The model runs in
    # minutes (the Table 3 rate constants and half-lives are per minute), so
    # CL = 2.6 / 60 = 0.043333 L/min/kg. Vc + Vp = 4.1 L/kg reproduces the
    # Table 3 Vss of 4.1 +/- 1.9 L/kg.
    lcl <- fixed(log(2.6 / 60)); label("Clearance CL (L/min/kg)") # Table 3, HS 'CL (L/h/kg)' = 2.6, divided by 60
    lvc <- fixed(log(0.6)); label("Central volume of distribution Vc (L/kg)") # Table 3, HS 'Vc (L/kg)' = 0.6
    lvp <- fixed(log(3.5)); label("Peripheral volume of distribution Vp (L/kg)") # Table 3, HS 'Vp (L/kg)' = 3.5

    # Intercompartmental clearance. Stach 2026 estimated k12 and k21 (Methods,
    # Pharmacokinetic Analysis) but Table 3 does not print them. With CL, Vc
    # and Vp fixed above, Q is the single remaining unknown; it is solved so
    # the beta (terminal) half-life of the two-compartment system equals the
    # Table 3 HS median t1/2beta = 78.1 min, giving Q = 0.1700 L/min/kg.
    # The implied distribution half-life is then t1/2alpha = 1.75 min against
    # the independently reported 2.4 +/- 0.9 min (Table 3).
    lq <- fixed(log(0.1700)); label("Intercompartmental clearance Q (L/min/kg)") # back-solved from Table 3 HS t1/2beta = 78.1 min

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
