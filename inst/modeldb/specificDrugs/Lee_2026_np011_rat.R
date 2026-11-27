Lee_2026_np011_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model with first-order",
    "elimination for NP-011, a truncated recombinant human milk fat globule-EGF",
    "factor 8 (MFG-E8) protein for MASH-associated liver fibrosis, after",
    "intravenous bolus doses. Parameters are per kg body weight, so doses are",
    "entered in ug/kg. Between-animal variability on CL only; log-normal",
    "(proportional on the log scale) residual error."
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
    species = "rat (strain not stated)",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "6-9 weeks",
    weight_range = "228.3-368.3 g",
    sex_female_pct = 25,
    disease_state = "healthy",
    dose_range = "0.00798, 0.0426 and 0.171 mg/kg single IV bolus (7.98-171 ug/kg)",
    regions = "preclinical (Republic of Korea)",
    notes = paste(
      "Table S1: 3 male rats per dose group plus 3 female rats at 0.0426",
      "mg/kg; serial sampling (predose, 0.083-24 h) in the same animals.",
      "Plasma NP-011 by ELISA; BLQ observations excluded (M1); one anomalous",
      "concentration excluded. NONMEM 7.5 FOCE-I, ADVAN6/TRANS1,",
      "natural-log-transformed DV. The rat fit showed strong parameter",
      "correlations (CL-VC 0.958; Q-IIV CL 0.968; Results 3.2)."
    )
  )

  ini({
    # Lee 2026 Table 1, 'Rats' column (final estimates). Structural parameters
    # are normalized to body weight (Table 1 notes; Supplementary Methods).
    lcl <- log(0.552); label("Clearance CL (L/h/kg)") # Table 1 Rats: CL = 0.552 L/h/kg (RSE 20.1)
    lvc <- log(0.0814); label("Central volume of distribution VC (L/kg)") # Table 1 Rats: VC = 0.0814 L/kg (RSE 22.2)
    lq <- log(0.0176); label("Intercompartmental clearance Q (L/h/kg)") # Table 1 Rats: Q = 0.0176 L/h/kg (RSE 44.0)
    lvp <- log(0.0627); label("Peripheral volume of distribution VP (L/kg)") # Table 1 Rats: VP = 0.0627 L/kg (RSE 22.6)

    # IIV reported as omega^2 (variance of the exponential random effect).
    etalcl ~ 0.0198 # Table 1 Rats: IIV on CL omega2 = 0.0198 (RSE 166.2)

    # Residual error on the natural-log concentration scale (Supplementary
    # Methods; Table S6 'Proportional (log scale)'), i.e. lnorm in nlmixr2.
    # The additive component was fixed to a negligible value and is omitted.
    expSd <- 0.985; label("Residual error SD on the log-concentration scale") # Table 1 Rats: proportional error = 0.985 (RSE 30.8)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
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
