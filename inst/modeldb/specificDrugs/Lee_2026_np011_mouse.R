Lee_2026_np011_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment population PK model with first-order",
    "elimination for NP-011, a truncated recombinant human milk fat globule-EGF",
    "factor 8 (MFG-E8) protein for MASH-associated liver fibrosis, after",
    "intravenous bolus doses. Parameters are per kg body weight, so doses are",
    "entered in ug/kg. Between-animal variability on CL and Vp; log-normal",
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
    species = "mouse (strain not stated)",
    n_subjects = 120L,
    n_studies = 1L,
    age_range = "6-10 weeks",
    weight_range = "23.60-38.12 g",
    sex_female_pct = 25,
    disease_state = "healthy",
    dose_range = "0.0214, 0.0854 and 0.342 mg/kg single IV bolus (21.4-342 ug/kg)",
    regions = "preclinical (Republic of Korea)",
    notes = paste(
      "Table S1: 30 male mice per dose group plus 30 female mice at 0.0854",
      "mg/kg. Cross-sectional (terminal) sampling, three animals per nominal",
      "time point (predose, 0.083-24 h), so each mouse contributed one",
      "concentration. Plasma NP-011 by ELISA; BLQ observations excluded (M1).",
      "NONMEM 7.5 FOCE-I, ADVAN6/TRANS1, natural-log-transformed DV."
    )
  )

  ini({
    # Lee 2026 Table 1, 'Mice' column (final estimates). Structural parameters
    # are normalized to body weight (Table 1 notes; Supplementary Methods).
    lcl <- log(0.363); label("Clearance CL (L/h/kg)") # Table 1 Mice: CL = 0.363 L/h/kg (RSE 8.9)
    lvc <- log(0.071); label("Central volume of distribution VC (L/kg)") # Table 1 Mice: VC = 0.071 L/kg (RSE 8.2)
    lq <- log(0.0063); label("Intercompartmental clearance Q (L/h/kg)") # Table 1 Mice: Q = 0.0063 L/h/kg (RSE 18.9)
    lvp <- log(0.0174); label("Peripheral volume of distribution VP (L/kg)") # Table 1 Mice: VP = 0.0174 L/kg (RSE 44.8)

    # IIV reported as omega^2 (variance of the exponential random effect).
    etalcl ~ 0.0133 # Table 1 Mice: IIV on CL omega2 = 0.0133 (RSE 56.9)
    etalvp ~ 0.227 # Table 1 Mice: IIV on VP omega2 = 0.227 (RSE 45.8)

    # Residual error on the natural-log concentration scale (Supplementary
    # Methods; Table S6 'Proportional (log scale)'), i.e. lnorm in nlmixr2.
    # The additive component was fixed to a negligible value and is omitted.
    expSd <- 0.361; label("Residual error SD on the log-concentration scale") # Table 1 Mice: proportional error = 0.361 (RSE 8.0)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

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
