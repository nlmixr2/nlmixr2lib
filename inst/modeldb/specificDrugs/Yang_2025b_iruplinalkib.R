Yang_2025b_iruplinalkib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "first-order elimination for oral iruplinalkib (WX-0593, a selective",
    "ALK/ROS1 tyrosine kinase inhibitor), fitted to the early clinical data",
    "used to select the recommended phase 2 dose: the phase 1 dose-escalation",
    "and expansion trial WX-0593-001 in patients with ALK- or ROS1-positive",
    "advanced non-small-cell lung cancer and the single-dose 120 mg",
    "food-effect study WX-0593-002 in healthy volunteers. Apparent oral",
    "clearance carries power effects of baseline body weight and serum",
    "albumin, and the apparent central volume a power effect of baseline body",
    "weight. Dosing in the fed state multiplies the absorption rate constant",
    "by 0.796 and adds a 0.479 h absorption lag. Residual error is",
    "proportional, with a separate magnitude for each of the two studies."
  )

  reference <- paste(
    "Yang G, Wang Y, Zhao H, Zheng Q, Wang X, Jia L, Xin Q, Ma C, Zhang Y,",
    "Zheng S, Kang X.",
    "Translational Model-Informed Dose Selection for Iruplinalkib, a",
    "Selective Oral ALK/ROS1 Tyrosine Kinase Inhibitor.",
    "Clin Transl Sci. 2025;18(7):e70287.",
    "doi:10.1111/cts.70287"
  )

  vignette <- "Yang_2025b_iruplinalkib"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "iruplinalkib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "iruplinalkib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iruplinalkib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT_BASE = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed BASELINE body weight: Table 3 names the rows 'Baseline",
        "body weight on CL/F' and 'Baseline body weight on V1/F', and the",
        "Table 3 footnote equations write the term as (BBWT/63). Power",
        "effects on CL/F (exponent 0.713) and V1/F (exponent 2.02), both",
        "normalised to 63 kg. The paper prints no demographic table, so the",
        "cohort weight distribution behind the 63 kg reference is not",
        "reported."
      ),
      source_name = "BBWT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TIME-VARYING serum albumin (Results: 'time-varying albumin was a",
        "significant covariate of CL/F'). Power effect on CL/F with exponent",
        "1.23, normalised to 42 (Table 3 footnote equation (ALB/42)^1.23).",
        "The paper does not print the units; 42 is read as g/L, the only",
        "plausible scale for a human serum albumin reference value (42 g/dL",
        "would be ten times the physiological range). The exponent is",
        "positive, so lower albumin lowers CL/F and raises exposure."
      ),
      source_name = "ALB"
    ),
    FED = list(
      description = "Fed-versus-fasted dose-record indicator, 1 = fed",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Two effects on absorption only (Table 3): the absorption rate",
        "constant is multiplied by 0.796 when fed ('Food on Ka', footnote:",
        "'Effect of food on Ka with reference to fasting'), and a 0.479 h",
        "absorption lag is added when fed (row 'ALAG for fed subjects';",
        "there is no lag in the fasted state). The only fed data come from",
        "the healthy-volunteer food-effect study WX-0593-002, in which",
        "subjects received a single 120 mg dose. The meal composition is not",
        "stated, so the general FED canonical is used rather than",
        "FED_HIGHFAT."
      ),
      source_name = "Food"
    ),
    STUDY_WX0593_002 = list(
      description = "WX-0593-002 study indicator: 1 = the healthy-volunteer food-effect study, 0 = the phase 1 patient trial WX-0593-001",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (WX-0593-001, the phase 1 trial in patients with non-small-cell lung cancer)",
      notes = paste(
        "Switches the proportional residual-error magnitude only: Table 3",
        "reports 'Residual proportional errors for WX-0593-001' = 32.6% and",
        "'... for WX-0593-002' = 54.7%. This analysis pooled only these two",
        "studies, so the reference group is WX-0593-001 alone. Subject type",
        "(healthy subject or patient) was screened as a structural covariate",
        "and not retained, so this indicator carries no structural effect."
      ),
      source_name = "Study"
    )
  )

  # Covariates screened by stepwise covariate modelling (Methods 'Modeling')
  # and not retained in the final model (Table 3).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened by SCM and not retained."
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Screened by SCM and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened by SCM and not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened by SCM and not retained (it was retained in the later four-study analysis, Yang_2025_iruplinalkib)."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator, 1 = healthy subject, 0 = patient",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'subject type' by SCM and not retained on any structural parameter."
    ),
    WHO_PS = list(
      description = "Eastern Cooperative Oncology Group performance status",
      units = "(score 0-5)",
      type = "continuous",
      notes = "Screened by SCM and not retained."
    )
  )

  population <- list(
    species = "human",
    n_studies = 2L,
    disease_state = paste(
      "WX-0593-001: Chinese patients with ALK- or ROS1-positive advanced",
      "non-small-cell lung cancer in the first-in-human dose-escalation and",
      "expansion trial (ChiCTR20170871 / NCT03389815). WX-0593-002: healthy",
      "Chinese volunteers in a food-effect study (CTR20190737)"
    ),
    dose_range = paste(
      "WX-0593-001: once-daily oral doses across the escalation range",
      "(54 patients dosed at 120, 180 or 240 mg QD were the cohort for the",
      "steady-state simulations). WX-0593-002: a single oral 120 mg dose"
    ),
    regions = "China",
    notes = paste(
      "The paper prints no demographic table and does not give the number",
      "of subjects or observations in the clinical PopPK dataset; the only",
      "count stated is the 54 WX-0593-001 patients used as the simulation",
      "cohort (Methods 'Establishing the Clinical PopPK Model and RP2D').",
      "Estimation was FOCEI in NONMEM 7.5 with PsN 4.8.1. Covariates screened",
      "by SCM (forward p < 0.05, backward p < 0.001): food, body weight, sex,",
      "age, albumin, AST, creatinine clearance, subject type and ECOG",
      "performance status. The same programme's later four-study pooled",
      "analysis (Yang 2025, doi:10.1111/cts.70099) is Yang_2025_iruplinalkib."
    )
  )

  ini({
    # Structural parameters. Typical values for the reference subject of the
    # Table 3 footnote equations: baseline body weight 63 kg and albumin
    # 42 g/L, dosed FASTED (FED = 0: no lag, no slowing of Ka). All
    # disposition parameters are apparent (per unit bioavailability);
    # iruplinalkib was given orally only.
    lka <- log(1.11); label("First-order absorption rate constant Ka, fasted (1/h)")    # Table 3: Ka = 1.11 1/h (RSE 12%)
    lcl <- log(20.5); label("Apparent oral clearance CL/F (L/h)")                       # Table 3: CL/F = 20.5 L/h (RSE 4%); footnote equation CL/F = 20.5*(BBWT/63)^0.713*(ALB/42)^1.23*e^eta1
    lvc <- log(342); label("Apparent central volume of distribution V1/F (L)")          # Table 3: V1/F = 342 L (RSE 7%); the footnote equation prints 348 (see vignette Errata)
    lq <- log(14.8); label("Apparent intercompartmental clearance Q/F (L/h)")           # Table 3: Q/F = 14.8 L/h (RSE 14%)
    lvp <- log(336); label("Apparent peripheral volume of distribution V2/F (L)")       # Table 3: V2/F = 336 L (RSE 12%)

    # Absorption lag in the fed state only; gated by FED in model().
    ltlag <- log(0.479); label("Absorption lag time in the fed state (h)")              # Table 3: ALAG for fed subjects = 0.479 h (RSE 2%)

    # Food effect on Ka: a multiplicative factor with fasting as reference
    # (Table 3 footnote), applied as e_fed_ka^FED.
    e_fed_ka <- 0.796; label("Multiplicative factor of the fed state on Ka, fasted reference (unitless)")  # Table 3: Food on Ka = 0.796 (RSE 27%); Results text says 58.2% (see vignette Errata)

    # Covariate effects (Table 3 footnote equations).
    e_alb_cl <- 1.23; label("Power exponent of serum albumin (/42 g/L) on CL/F (unitless)")              # Table 3: Albumin on CL/F = 1.23 (RSE 29%)
    e_wt_base_cl <- 0.713; label("Power exponent of baseline body weight (/63 kg) on CL/F (unitless)")   # Table 3: Baseline body weight on CL/F = 0.713 (RSE 37%)
    e_wt_base_vc <- 2.02; label("Power exponent of baseline body weight (/63 kg) on V1/F (unitless)")    # Table 3: Baseline body weight on V1/F = 2.02 (RSE 21%)

    # Inter-individual variability. Methods: 'The interindividual
    # variabilities (IIVs) for fixed effect parameters were modeled using the
    # exponential error model'; Table 3 reports IIV as a percentage and its
    # abbreviation list defines CV, so omega^2 = log(1 + CV^2). No IIV on
    # Q/F, V2/F or the lag ('-' in Table 3); no correlations reported.
    etalka ~ 0.513279 # Table 3: IIV Ka = 81.9% (IIV_RSE 13%, shrinkage 7%); log(1 + 0.819^2)
    etalcl ~ 0.094033 # Table 3: IIV CL/F = 31.4% (IIV_RSE 12%, shrinkage 3%); log(1 + 0.314^2)
    etalvc ~ 0.163174 # Table 3: IIV V1/F = 42.1% (IIV_RSE 24%, shrinkage 0.1%); log(1 + 0.421^2)

    # Residual error: proportional (Methods), one magnitude per study.
    propSdPatient <- 0.326; label("Proportional residual error SD for study WX-0593-001, patients (fraction)")       # Table 3: Residual proportional errors for WX-0593-001 = 32.6% (RSE 10%)
    propSdHv <- 0.547; label("Proportional residual error SD for study WX-0593-002, healthy volunteers (fraction)")  # Table 3: Residual proportional errors for WX-0593-002 = 54.7% (RSE 6%)
  })

  model({
    # Individual parameters (Table 3 footnote equations).
    ka <- exp(lka + etalka) * e_fed_ka^FED
    cl <- exp(lcl + etalcl) * (WT_BASE / 63)^e_wt_base_cl * (ALB / 42)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT_BASE / 63)^e_wt_base_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # Lag only when fed; zero in the fasted reference state.
    tlag <- exp(ltlag) * FED

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    propSd <- propSdHv * STUDY_WX0593_002 + propSdPatient * (1 - STUDY_WX0593_002)

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
