Lee_2026_np011_human_mlp <- function() {
  description <- paste(
    "Two-compartment PK model with first-order elimination for intravenous",
    "NP-011, a truncated recombinant human milk fat globule-EGF factor 8",
    "(MFG-E8) protein for MASH-associated liver fibrosis, in a 70 kg human,",
    "projected from the mouse, rat and monkey population PK fits of Lee 2026",
    "with maximum-lifespan-potential (MLP) corrected allometry for CL and Q",
    "and simple allometry for Vc and Vp (the paper's supplementary",
    "sensitivity analysis). Structural values are fixed; between-subject",
    "variability on CL and Vc and the residual error are carried from the",
    "monkey model, as for the simple-allometry projection."
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
      "(NCT05387499; 30 NP-011-treated participants, 0.25-4 mg IV bolus)"
    ),
    dose_range = "0.25, 0.5, 1, 2 and 4 mg single IV bolus (250-4000 ug)",
    regions = "Republic of Korea (Phase I comparison data)",
    notes = paste(
      "Rule of exponents: the CL and Q exponents (0.703, 0.867) fall in",
      "0.7 <= b < 1, which calls for the MLP correction (Table S8);",
      "MLP = 185.4 * BrW^0.636 * BW^-0.225 with the brain weights of Table S9.",
      "The paper retained simple allometry (Lee_2026_np011_human_simple) as",
      "its primary projection because the MLP correction over-predicted the",
      "observed Phase I AUCinf by 2.38- to 2.89-fold (Table S13)."
    )
  )

  ini({
    # Supplementary Results 'Human PK prediction using MLP-corrected
    # allometric scaling': CL 2.58 L/h and Q 0.17 L/h, with VC and VP
    # unchanged from simple allometry (2.88 L and 24.61 L). VC and VP use the
    # SIM_SA.mod values that Table 2 rounds. Fixed: allometric predictions,
    # not estimates from human data (Methods 2.4).
    lcl <- fixed(log(2.58)); label("Clearance CL (L/h)") # Supplementary Results MLP: CL = 2.58 L/h
    lvc <- fixed(log(2.88054)); label("Central volume of distribution VC (L)") # SIM_SA.mod THETA(2) = 2.88054; Supplementary Results MLP: Vc unchanged at 2.88 L
    lq <- fixed(log(0.17)); label("Intercompartmental clearance Q (L/h)") # Supplementary Results MLP: Q = 0.17 L/h
    lvp <- fixed(log(24.60960491)); label("Peripheral volume of distribution VP (L)") # SIM_SA.mod THETA(3) = 24.60960491; Supplementary Results MLP: Vp unchanged at 24.61 L

    # Variability as for the simple-allometry projection ('human profiles
    # were simulated using the same procedure', Supplementary Methods).
    etalcl ~ fixed(0.00216) # SIM_SA.mod $OMEGA IIV CL = 0.00216 (Table 1 Monkeys 0.0022)
    etalvc ~ fixed(0.0371) # SIM_SA.mod $OMEGA IIV Vc = 0.0371 (Table 1 Monkeys 0.0371)

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
