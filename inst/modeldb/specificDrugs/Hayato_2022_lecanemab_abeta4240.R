Hayato_2022_lecanemab_abeta4240 <- function() {
  description <- paste(
    "Indirect-response PK/PD model for the plasma amyloid-beta 42/40 ratio",
    "in subjects with early Alzheimer's disease receiving intravenous",
    "lecanemab: dR/dt = Kin * (1 + Slope * Cc) - Kout * R, with Kin = Kout",
    "* baseline, so serum lecanemab linearly stimulates production and",
    "raises the ratio. Exponential IIV on baseline and slope, proportional",
    "residual error, no covariates retained. Fit to 1254 measurements",
    "(C2N PrecivityAD IP/LC-MS/MS assay) from 284 subjects in the phase 2",
    "study 201 Core and open-label extension. The serum lecanemab driver",
    "is the companion two-compartment population PK model",
    "Hayato_2022_lecanemab, whose parameters are held fixed here (the",
    "paper used its post hoc individual PK parameters).",
    sep = " "
  )
  reference <- paste(
    "Hayato S, Takenaka O, Sreerama Reddy SH, Landry I, Reyderman L,",
    "Koyama A, Swanson C, Yasuda S, Hussein Z (2022).",
    "Population pharmacokinetic-pharmacodynamic analyses of amyloid",
    "positron emission tomography and plasma biomarkers for lecanemab in",
    "subjects with early Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol 11(12):1578-1591.",
    "doi:10.1002/psp4.12862.",
    "PK driver: modellib('Hayato_2022_lecanemab') (same paper).",
    sep = " "
  )
  vignette <- "Hayato_2022_lecanemab"

  # Plasma Abeta42/40 ratio turnover state; paper-specific PD readout (same
  # state name as Bhagunde_2026_lecanemab_abeta4240).
  paper_specific_compartments <- c("abeta4240")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "lecanemab", units = "mg", specimen = "tissue", verified = TRUE),
    abeta4240 = list(
      analyte = "amyloid-beta 42/40 ratio",
      units = "(ratio)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effects on CL and V1 normalised to 73.7 kg (Hayato 2022 Table 1). No covariate was retained on the Abeta42/40 parameters.",
      source_name = "WGT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effect on CL normalised to 42.9 g/L (Hayato 2022 Table 1).",
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "PK driver only: ratios on CL (0.792) and V1 (0.893).",
      source_name = "SEXN"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "PK driver only: time-varying ratio 1.09 on CL.",
      source_name = "ADA"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator (1 = Japanese, 0 = non-Japanese)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese)",
      notes = "PK driver only: ratio 0.455 on V2.",
      source_name = "RACEN == 3.1"
    ),
    FORM_LEC_PROCESSB = list(
      description = "Lecanemab manufacturing Process B drug-product indicator (1 = Process B, 0 = Process A)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Process A; F fixed at 1)",
      notes = "PK driver only: relative bioavailability 0.998 with 34.2% IIV on Process B doses (study 201 OLE). Set to 0 for study 201 Core dosing.",
      source_name = "FORM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 284L,
    n_studies = 1L,
    n_observations = "1254 plasma Abeta42/40 ratio measurements (study 201 Core and OLE)",
    age_range = "50-88 years (mean 71.1, SD 8.2; Table S2)",
    weight_range = "35.9-111.8 kg (mean 70.9, SD 14.5; Table S2)",
    sex_female_pct = 51.4,
    apoe4_carrier_pct = 71.1,
    baseline_abeta4240 = "mean 0.0845, SD 0.0082, range 0.0591-0.145 (Table S2)",
    disease_state = "Early Alzheimer's disease: MCI due to AD (66.5%) or mild AD dementia (33.5%)",
    dose_range = paste(
      "Study 201 Core (18 months): placebo (n = 88), 2.5 mg/kg biweekly",
      "(13), 5 mg/kg monthly (16), 5 mg/kg biweekly (29), 10 mg/kg monthly",
      "(95), 10 mg/kg biweekly (43); then a 9-59 month off-treatment gap",
      "and 10 mg/kg biweekly in the OLE (Table S2, Figure S1)",
      sep = " "
    ),
    regions = "Multinational (North America, Europe, Japan, South Korea)",
    notes = paste(
      "Abeta42/40 analysis column of supplement Table S2. ADA-positive at",
      "subject level 40.8%. Eta shrinkage: baseline 11.9%, slope 63.1%",
      "(Table 3 footnote).",
      sep = " "
    )
  )

  ini({
    # Lecanemab PK driver -- fixed at the Hayato_2022_lecanemab estimates
    # (Table 1); the Abeta42/40 analysis used post hoc individual PK
    # parameters (Model S3 $PK: 'CL = ICL', 'V1 = IV1', ...).
    lcl <- fixed(log(0.0181)); label("Lecanemab clearance for the reference subject (L/h)") # Table 1 CL = 0.0181 L/h
    lvc <- fixed(log(3.22)); label("Lecanemab central volume for the reference subject (L)") # Table 1 V1 = 3.22 L
    lq <- fixed(log(0.0349)); label("Lecanemab intercompartmental clearance (L/h)") # Table 1 Q = 0.0349 L/h
    lvp <- fixed(log(2.19)); label("Lecanemab peripheral volume for the reference subject (L)") # Table 1 V2 = 2.19 L
    lfcentral <- fixed(log(1)); label("Relative bioavailability for Process A (unitless)") # Model S1 $PK 'F1=1'
    e_processb_f <- fixed(0.998); label("Relative bioavailability for Process B vs Process A (unitless)") # Table 1 'F for process B' = 0.998
    e_wt_cl <- fixed(0.403); label("Power exponent on (WT/73.7 kg) for CL (unitless)") # Table 1
    e_alb_cl <- fixed(-0.243); label("Power exponent on (ALB/42.9 g/L) for CL (unitless)") # Table 1
    e_female_cl <- fixed(0.792); label("CL ratio for females vs males (unitless)") # Table 1
    e_ada_cl <- fixed(1.09); label("CL ratio for ADA-positive vs ADA-negative (unitless)") # Table 1
    e_wt_vc <- fixed(0.606); label("Power exponent on (WT/73.7 kg) for V1 (unitless)") # Table 1
    e_female_vc <- fixed(0.893); label("V1 ratio for females vs males (unitless)") # Table 1
    e_japanese_vp <- fixed(0.455); label("V2 ratio for Japanese vs non-Japanese (unitless)") # Table 1

    # Plasma Abeta42/40 indirect-response model (Table 3; Model S3)
    lrbase <- log(0.0842); label("Baseline plasma Abeta42/40 ratio (ratio)") # Table 3 'Baseline plasma Abeta42/40 ratio' = 0.0842 (%RSE 4.28)
    lkout <- log(0.367); label("First-order degradation rate constant Kout (1/year)") # Table 3 'Kout (1/year)' = 0.367 (%RSE 1.97); Model S3 'ABOUT = THETA(2)/8760'
    lslope <- log(0.00155); label("Linear stimulation of production per unit serum lecanemab (1/(ug/mL))") # Table 3 'Slope (1/ug/ml)' = 0.00155 (%RSE 9.32)

    # Lecanemab PK IIV, fixed at the Hayato_2022_lecanemab estimates (Table 1;
    # CV% = sqrt(variance) x 100, so omega^2 = (CV/100)^2).
    etalcl ~ fixed(0.151321) # Table 1 IIV CL 38.9%
    etalvc ~ fixed(0.0196) # Table 1 IIV V1 14.0%
    etalvp ~ fixed(0.990025) # Table 1 IIV V2 99.5%
    etalfcentral ~ fixed(0.116964) # Table 1 IIV F 34.2%; Process B doses only

    # Model S3 $PK: 'BSL = THETA(1)*EXP(ETA(1))', 'SLOPE = THETA(3)*EXP(ETA(2))';
    # diagonal $OMEGA. CV% = sqrt(variance) x 100 (Table 3 footnote).
    etalrbase ~ 0.00459684 # Table 3 IIV baseline 6.78% -> 0.0678^2
    etalslope ~ 0.194481 # Table 3 IIV slope 44.1% -> 0.441^2

    # Model S3 $ERROR: 'W = IPRED*THETA(4); Y = IPRED + W*ERR(1)' with
    # $SIGMA 1 FIX, so THETA(4) is the proportional SD.
    propSd <- 0.0641; label("Proportional residual error on the Abeta42/40 ratio (fraction)") # Table 3 'Proportional (CV%)' = 6.41 (%RSE 2.82)
  })

  model({
    hrs_per_year <- 8760 # Model S3 'ABOUT = THETA(2)/8760': 365 d x 24 h

    # Lecanemab PK driver (Hayato 2022 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / 73.7)^e_wt_cl * (ALB / 42.9)^e_alb_cl *
      e_female_cl^SEXF * e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 73.7)^e_wt_vc * e_female_vc^SEXF
    vp <- exp(lvp + etalvp) * e_japanese_vp^RACE_JAPANESE
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual Abeta42/40 parameters (Model S3 $PK)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout) / hrs_per_year
    kin <- kout * rbase
    slope <- exp(lslope + etalslope)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Model S1 $PK: 'F1=1; IF (FORM.EQ.1) F1=THETA(11)*EXP(ETA(4))'
    fprocessb <- exp(lfcentral + etalfcentral) * e_processb_f
    f(central) <- (1 - FORM_LEC_PROCESSB) + FORM_LEC_PROCESSB * fprocessb

    Cc <- central / vc

    # Figure 1 / Model S3 $DES: dR/dt = Kin * (1 + Slope * C) - Kout * R
    d/dt(abeta4240) <- kin * (1 + slope * Cc) - kout * abeta4240
    abeta4240(0) <- rbase

    abeta4240 ~ prop(propSd)
  })
}
