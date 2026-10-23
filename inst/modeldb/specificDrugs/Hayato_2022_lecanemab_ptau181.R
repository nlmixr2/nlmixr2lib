Hayato_2022_lecanemab_ptau181 <- function() {
  description <- paste(
    "Indirect-response PK/PD model for plasma tau phosphorylated at",
    "threonine 181 (p-tau181) in subjects with early Alzheimer's disease",
    "receiving intravenous lecanemab: dR/dt = Kin * (1 - Slope * (Cc +",
    "0.01)) - Kout * R, with Kin = Kout * baseline, so serum lecanemab",
    "linearly inhibits production. Baseline falls with body weight (power,",
    "reference 72.2 kg); IIV is exponential on baseline and additive on",
    "slope; residual error is proportional. Fit to 2021 measurements",
    "(Quanterix Simoa assay) from 562 subjects in the phase 2 study 201",
    "Core and open-label extension. The serum lecanemab driver is the",
    "companion two-compartment population PK model Hayato_2022_lecanemab,",
    "whose parameters are held fixed here (the paper used its post hoc",
    "individual PK parameters).",
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

  # Plasma p-tau181 turnover state; paper-specific PD readout (same state
  # name as Bhagunde_2026_lecanemab_ptau181).
  paper_specific_compartments <- c("ptau181")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "lecanemab", units = "mg", specimen = "tissue", verified = TRUE),
    ptau181 = list(
      analyte = "tau phosphorylated at threonine 181 (p-tau181)",
      units = "pg/mL",
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
      notes = paste(
        "Enters twice. (1) PK driver: power effects on CL and V1 normalised",
        "to 73.7 kg (Hayato 2022 Table 1). (2) Baseline p-tau181: power",
        "effect (WT/72.2)^-0.300 (Model S4 $PK '(BWGTC/72.2)**THETA(5)';",
        "Table 3). A 50 kg subject has an 11.7% higher and a 96 kg subject",
        "an 8.2% lower baseline than the 72 kg reference (Results).",
        sep = " "
      ),
      source_name = "WGT (PK); BWGTC (Model S4)"
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
      notes = "PK driver only: ratios on CL (0.792) and V1 (0.893). Sex was tested on the p-tau181 parameters and not retained (Table S9).",
      source_name = "SEXN"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = paste(
        "PK driver only: time-varying ratio 1.09 on CL. ADA status on the",
        "p-tau181 baseline was significant in univariate screening (Table",
        "S9) but removed in backward elimination (Table S10); Model S4 keeps",
        "its THETA(6) as '1 FIX', so it has no effect and is not encoded.",
        sep = " "
      ),
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

  covariatesDataExcluded <- list(
    DIS_AD_MILD = list(
      description = "Baseline diagnosis (mild AD dementia vs MCI due to AD)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significant on baseline p-tau181 in univariate screening (Table S9,",
        "dOFV 9.16) and removed in backward elimination (Table S10). Model",
        "S4 keeps its THETA(4) as '1 FIX', so it has no effect.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 562L,
    n_studies = 1L,
    n_observations = "2021 plasma p-tau181 measurements (study 201 Core and OLE)",
    age_range = "50-89 years (mean 71.0, SD 8.3; Table S2)",
    weight_range = "29.2-118.7 kg (mean 72.5, SD 14.6; Table S2)",
    sex_female_pct = 50.9,
    race_ethnicity = c(
      White = 89.0,
      `Black/African American` = 2.1,
      `Asian/Other Asian (excluding Chinese and Japanese)` = 2.7,
      Japanese = 5.7,
      `Chinese/other/missing` = 0.5
    ),
    apoe4_carrier_pct = 69.9,
    baseline_ptau181 = "mean 4.43 pg/mL, SD 1.85, range 0.84-17.4 (Table S2)",
    disease_state = "Early Alzheimer's disease: MCI due to AD (66.5%) or mild AD dementia (33.5%)",
    dose_range = paste(
      "Study 201 Core (18 months): placebo (n = 179), 2.5 mg/kg biweekly",
      "(36), 5 mg/kg monthly (38), 5 mg/kg biweekly (70), 10 mg/kg monthly",
      "(155), 10 mg/kg biweekly (84); then a 9-59 month off-treatment gap",
      "and 10 mg/kg biweekly in the OLE (Table S2, Figure S1)",
      sep = " "
    ),
    regions = "Multinational (North America, Europe, Japan, South Korea)",
    notes = paste(
      "p-tau181 analysis column of supplement Table S2. ADA-positive at",
      "subject level 41.8%. Eta shrinkage: baseline 4.65%, slope 68.4%",
      "(Table 3 footnote).",
      sep = " "
    )
  )

  ini({
    # Lecanemab PK driver -- fixed at the Hayato_2022_lecanemab estimates
    # (Table 1); the p-tau181 analysis used post hoc individual PK
    # parameters (Model S4 $PK: 'CL = ICL', 'V1 = IV1', ...).
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

    # Plasma p-tau181 indirect-response model (Table 3; Model S4)
    lrbase <- log(4.06); label("Baseline plasma p-tau181 at 72.2 kg (pg/mL)") # Table 3 'Baseline plasma p-tau181' = 4.06 (%RSE 1.61)
    lkout <- log(0.468); label("First-order degradation rate constant Kout (1/year)") # Table 3 'Kout (1/year)' = 0.468 (%RSE 20.7); Model S4 'KOUT = (THETA(2)/8760)'
    slope <- 0.00313; label("Linear inhibition of production per unit serum lecanemab (1/(ug/mL))") # Table 3 'Slope (1/ug/ml)' = 0.00313 (%RSE 15.6); Model S4 'SLOPE = THETA(3)+ETA(2)' (additive IIV, so not log-transformed)
    e_wt_rbase <- -0.300; label("Power exponent on (WT/72.2 kg) for baseline p-tau181 (unitless)") # Table 3 'Weight ~ baseline (exponent)' = -0.300 (%RSE 24.2)

    # Lecanemab PK IIV, fixed at the Hayato_2022_lecanemab estimates (Table 1;
    # CV% = sqrt(variance) x 100, so omega^2 = (CV/100)^2).
    etalcl ~ fixed(0.151321) # Table 1 IIV CL 38.9%
    etalvc ~ fixed(0.0196) # Table 1 IIV V1 14.0%
    etalvp ~ fixed(0.990025) # Table 1 IIV V2 99.5%
    etalfcentral ~ fixed(0.116964) # Table 1 IIV F 34.2%; Process B doses only

    # Model S4 diagonal $OMEGA: exponential on baseline, additive on slope.
    etalrbase ~ 0.123201 # Table 3 IIV baseline 35.1% -> 0.351^2
    etaslope ~ 2.2801e-6 # Table 3 IIV slope SD 0.00151 (1/(ug/mL)) -> 0.00151^2

    propSd <- 0.194; label("Proportional residual error on plasma p-tau181 (fraction)") # Table 3 'Proportional (CV%)' = 19.4 (%RSE 2.39)
  })

  model({
    hrs_per_year <- 8760 # Model S4 'KOUT = (THETA(2)/8760)': 365 d x 24 h

    # Lecanemab PK driver (Hayato 2022 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / 73.7)^e_wt_cl * (ALB / 42.9)^e_alb_cl *
      e_female_cl^SEXF * e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 73.7)^e_wt_vc * e_female_vc^SEXF
    vp <- exp(lvp + etalvp) * e_japanese_vp^RACE_JAPANESE
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual p-tau181 parameters (Model S4 $PK; the diagnosis and ADA
    # thetas are '1 FIX' there and drop out)
    rbase <- exp(lrbase + etalrbase) * (WT / 72.2)^e_wt_rbase
    kout <- exp(lkout) / hrs_per_year
    kin <- kout * rbase
    slope_i <- slope + etaslope

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Model S1 $PK: 'F1=1; IF (FORM.EQ.1) F1=THETA(11)*EXP(ETA(4))'
    fprocessb <- exp(lfcentral + etalfcentral) * e_processb_f
    f(central) <- (1 - FORM_LEC_PROCESSB) + FORM_LEC_PROCESSB * fprocessb

    Cc <- central / vc

    # Model S4 $DES: 'DADT(3) = KIN*(1 - SLOPE*(C1A+0.01)) - A(3)*KOUT'.
    # The 0.01 ug/mL offset is in the published control stream (Figure 1
    # prints the equation without it); it lowers the untreated steady state
    # by a factor of 1 - 0.00313 * 0.01, i.e. 0.003%.
    d/dt(ptau181) <- kin * (1 - slope_i * (Cc + 0.01)) - kout * ptau181
    ptau181(0) <- rbase

    ptau181 ~ prop(propSd)
  })
}
