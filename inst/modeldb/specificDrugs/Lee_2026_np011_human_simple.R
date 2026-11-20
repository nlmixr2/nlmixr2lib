Lee_2026_np011_human_simple <- function() {
  description <- paste(
    "Two-compartment PK model with first-order elimination for intravenous",
    "NP-011, a truncated recombinant human milk fat globule-EGF factor 8",
    "(MFG-E8) protein for MASH-associated liver fibrosis, in a 70 kg human,",
    "projected by simple allometry from the mouse, rat and monkey population",
    "PK fits of Lee 2026 (the paper's primary scaling method). Structural",
    "values are fixed; between-subject variability on CL and Vc and the",
    "residual error are carried from the monkey model, as in the authors'",
    "simulation control stream."
  )
  reference <- paste(
    "Lee H, Tran QT, Tham BT, Kim MK, Yun HY, Chae JW, Park JH, Lee S.",
    "Population Pharmacokinetics and Extrapolation of NP-011 from Animal to",
    "Humans: A Novel Therapeutic Candidate Protein for Metabolic",
    "Dysfunction-Associated Steatohepatitis.",
    "Pharmaceutics. 2026;18(9):1194. doi:10.3390/pharmaceutics18091194",
    sep = " "
  )
  vignette <- "Lee_2026_np011"

  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "NP-011", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "NP-011", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human (projected by allometry from mouse, rat and monkey)",
    n_subjects = 0L,
    n_studies = 0L,
    weight_range = "70 kg (reference body weight of the projection)",
    disease_state = paste(
      "not applicable -- forward projection; retrospectively compared with",
      "the single ascending dose part of the Phase I study in healthy adults",
      "(NCT05387499; 30 NP-011-treated participants, 0.25-4 mg IV bolus,",
      "mean body weight 72.81 kg)"
    ),
    dose_range = "0.25, 0.5, 1, 2 and 4 mg single IV bolus (250-4000 ug)",
    regions = "Republic of Korea (Phase I comparison data)",
    notes = paste(
      "Simple allometry Y = a * BW^b fitted by unweighted least squares on",
      "log total parameters of the three species (Equation 1; Table S10:",
      "exponents CL 0.703, VC 0.914, Q 0.867, VP 1.374). Only the 70 kg",
      "prediction is encoded, so there is no weight covariate. The model",
      "over-predicted the observed Phase I AUCinf by 1.62- to 1.96-fold and",
      "Cmax by 1.33- to 1.48-fold (Table 4). The paper's simulation intervals",
      "exclude residual variability."
    )
  )

  ini({
    # Supplementary SIM_SA.mod $THETA (the authors' simulation control
    # stream); rounded in Table 2 to 3.80 L/h, 2.88 L, 0.25 L/h and 24.61 L.
    # Fixed: the paper states the human parameters were not estimated from
    # clinical data but obtained by allometry and fixed for simulation
    # (Methods 2.4).
    lcl <- fixed(log(3.79590)); label("Clearance CL (L/h)") # SIM_SA.mod THETA(1) = 3.79590; Table 2 CL = 3.80 L/h
    lvc <- fixed(log(2.88054)); label("Central volume of distribution VC (L)") # SIM_SA.mod THETA(2) = 2.88054; Table 2 VC = 2.88 L
    lq <- fixed(log(0.25010)); label("Intercompartmental clearance Q (L/h)") # SIM_SA.mod THETA(4) = 0.25010; Table 2 Q = 0.25 L/h
    lvp <- fixed(log(24.60960491)); label("Peripheral volume of distribution VP (L)") # SIM_SA.mod THETA(3) = 24.60960491; Table 2 VP = 24.61 L

    # Human IIV assumed equal to the monkey estimates (Methods 2.3).
    etalcl ~ fixed(0.00216) # SIM_SA.mod $OMEGA IIV CL = 0.00216 (Table 1 Monkeys 0.0022)
    etalvc ~ fixed(0.0371) # SIM_SA.mod $OMEGA IIV Vc = 0.0371 (Table 1 Monkeys 0.0371)

    # SIM_SA.mod carries the monkey log-scale residual SD; the paper's
    # prediction intervals were built from the error-free prediction.
    expSd <- fixed(0.457); label("Residual error SD on the log-concentration scale") # SIM_SA.mod THETA(5) = 0.457 (Table 1 Monkeys proportional error)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses in ug and vc in L give ug/L, which equals ng/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
