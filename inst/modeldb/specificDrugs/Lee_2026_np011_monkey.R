Lee_2026_np011_monkey <- function() {
  description <- paste(
    "Preclinical (monkey). Two-compartment population PK model with",
    "first-order elimination for NP-011, a truncated recombinant human milk",
    "fat globule-EGF factor 8 (MFG-E8) protein for MASH-associated liver",
    "fibrosis, after single and repeated (once daily for 14 days) intravenous",
    "bolus doses. Parameters are per kg body weight, so doses are entered in",
    "ug/kg. Between-animal variability on CL and Vc; log-normal (proportional",
    "on the log scale) residual error."
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

  units <- list(time = "h", dosing = "ug/kg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "NP-011", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "NP-011", units = "ug/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "monkey (species not stated)",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "6-10 months",
    weight_range = "2.55-4.75 kg",
    sex_female_pct = 25,
    disease_state = "healthy",
    dose_range = paste(
      "0.00534 and 0.0854 mg/kg single IV bolus; 0.0213 mg/kg IV bolus once",
      "daily for 14 days (5.34-85.4 ug/kg)"
    ),
    regions = "preclinical (Republic of Korea)",
    notes = paste(
      "Table S1: 3 male monkeys per dose group plus 3 female monkeys in the",
      "0.0213 mg/kg repeated-dose group; serial sampling (predose, 0.083-24 h;",
      "intensive on days 1 and 14 of the repeated-dose study, predose on days",
      "3-13). Plasma NP-011 by ELISA; BLQ observations excluded (M1). NONMEM",
      "7.5 FOCE-I, ADVAN6/TRANS1, natural-log-transformed DV. The monkey IIV",
      "and residual error were carried into the human projection",
      "(Lee_2026_np011_human_simple)."
    )
  )

  ini({
    # Lee 2026 Table 1, 'Monkeys' column (final estimates). Structural
    # parameters are normalized to body weight (Table 1 notes; Supplementary
    # Methods).
    lcl <- log(0.0904); label("Clearance CL (L/h/kg)") # Table 1 Monkeys: CL = 0.0904 L/h/kg (RSE 3.8)
    lvc <- log(0.0475); label("Central volume of distribution VC (L/kg)") # Table 1 Monkeys: VC = 0.0475 L/kg (RSE 6.0)
    lq <- log(0.00342); label("Intercompartmental clearance Q (L/h/kg)") # Table 1 Monkeys: Q = 0.00342 L/h/kg (RSE 9.3)
    lvp <- log(0.103); label("Peripheral volume of distribution VP (L/kg)") # Table 1 Monkeys: VP = 0.103 L/kg (RSE 22.5)

    # IIV reported as omega^2 (variance of the exponential random effect).
    # The supplementary SIM_SA.mod carries 0.00216 for the CL term.
    etalcl ~ 0.0022 # Table 1 Monkeys: IIV on CL omega2 = 0.0022 (RSE 131.5)
    etalvc ~ 0.0371 # Table 1 Monkeys: IIV on VC omega2 = 0.0371 (RSE 41.8)

    # Residual error on the natural-log concentration scale (Supplementary
    # Methods; Table S6 'Proportional (log scale)'), i.e. lnorm in nlmixr2.
    # The additive component was fixed to a negligible value and is omitted.
    expSd <- 0.457; label("Residual error SD on the log-concentration scale") # Table 1 Monkeys: proportional error = 0.457 (RSE 10.0)
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

    # Doses in ug/kg and vc in L/kg give ug/L, which equals ng/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
