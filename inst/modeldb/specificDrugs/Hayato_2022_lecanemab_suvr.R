Hayato_2022_lecanemab_suvr <- function() {
  description <- paste(
    "Indirect-response PK/PD model for the amyloid PET standard uptake",
    "value ratio (SUVr, florbetapir, whole-cerebellum reference) in",
    "subjects with early Alzheimer's disease receiving intravenous",
    "lecanemab: dSUVr/dt = Kin - Kout * SUVr * (1 + Emax * Cc / (EC50 +",
    "Cc)), with Kin estimated and Kout = Kin / baseline. APOE4 carriers",
    "have a 4% higher baseline and Emax rises with age (power, reference",
    "72 years); baseline and Emax IIV are correlated; residual error is",
    "proportional. Fit to 1213 SUVr measurements from 374 subjects in the",
    "phase 2 study 201 Core and open-label extension. The serum lecanemab",
    "driver is the companion two-compartment population PK model",
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

  # Amyloid PET SUVr turnover state; paper-specific PD readout (same state
  # name as KandadiMuralidharan_2022_aducanumab).
  paper_specific_compartments <- c("suvr")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "lecanemab", units = "mg", specimen = "tissue", verified = TRUE),
    suvr = list(
      analyte = "amyloid-beta plaque (florbetapir PET SUVr, whole-cerebellum reference)",
      units = "(ratio)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effects on CL and V1 normalised to 73.7 kg (Hayato 2022 Table 1). Not retained on the SUVr parameters (Table S5).",
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
      notes = "PK driver only: ratios on CL (0.792) and V1 (0.893). Sex was tested on the SUVr parameters and not retained (Table S5).",
      source_name = "SEXN"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "PK driver only: time-varying ratio 1.09 on CL. ADA and NAb status were tested on the SUVr parameters and not retained (Table S5).",
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
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Emax, (AGE/72)^1.58 (Model S2 $PK: 'EMAX = THETA(3)*(AGE/72)**THETA(5)'). Older subjects have a larger maximum plaque-removal effect.",
      source_name = "AGE"
    ),
    APOE4_CARRIER = list(
      description = "APOE-epsilon4 carrier status (1 = carrier, 0 = non-carrier)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-carrier)",
      notes = paste(
        "Ratio 1.04 on baseline SUVr, applied as ratio^APOE4_CARRIER. The",
        "Table 2 row label ('APOE4 carrier ~ baseline'), the Results ('APOE4",
        "carrier subjects have higher baseline SUVr') and the Discussion",
        "('APOE4 carriers have slightly higher baseline SUVr than",
        "non-carriers') all make carriers the 1-level; the Table 2 footnote",
        "'APOE = 0 (APOE4 carrier) or 1 (non-carrier)' contradicts all three",
        "and is treated as a typographical inversion. See the vignette",
        "Assumptions and deviations section.",
        sep = " "
      ),
      source_name = "APOE (Model S2 $INPUT; APO = 1 if APOE == 1)"
    )
  )

  covariatesDataExcluded <- list(
    DIS_AD_MILD = list(
      description = "Baseline diagnosis (mild AD dementia vs MCI due to AD)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on baseline, Kin, Emax and EC50 and not significant in univariate screening (Table S5, models 18-21)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 374L,
    n_studies = 1L,
    n_observations = "1213 amyloid PET SUVr measurements (study 201 Core and OLE)",
    age_range = "50-89 years (mean 71.5, SD 8.1; Table S2)",
    weight_range = "29.2-118.7 kg (mean 74.3, SD 15.1; Table S2)",
    sex_female_pct = 47.9,
    apoe4_carrier_pct = 69.8,
    baseline_suvr = "mean 1.38, SD 0.17, range 0.76-1.84 (Table S2)",
    disease_state = "Early Alzheimer's disease: MCI due to AD (70.1%) or mild AD dementia (29.9%)",
    dose_range = paste(
      "Study 201 Core (18 months): placebo (n = 115), 2.5 mg/kg biweekly",
      "(30), 5 mg/kg monthly (30), 5 mg/kg biweekly (36), 10 mg/kg monthly",
      "(105), 10 mg/kg biweekly (58); then a 9-59 month off-treatment gap",
      "and 10 mg/kg biweekly in the OLE (Table S2, Figure S1)",
      sep = " "
    ),
    regions = "Multinational (North America, Europe, Japan, South Korea)",
    notes = paste(
      "SUVr analysis column of supplement Table S2. ADA-positive at subject",
      "level 43.0%. Eta shrinkage: baseline 6.94%, Emax 27.7% (Table 2",
      "footnote). Amyloid negativity is SUVr < 1.17.",
      sep = " "
    )
  )

  ini({
    # Lecanemab PK driver -- fixed at the Hayato_2022_lecanemab estimates
    # (Table 1); the SUVr analysis used post hoc individual PK parameters
    # (Model S2 $PK: 'CL = ICL', 'V1 = IV1', ...).
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

    # Amyloid PET SUVr indirect-response model (Table 2; Model S2)
    lrbase <- log(1.34); label("Baseline SUVr for an APOE4 non-carrier (ratio)") # Table 2 'Baseline SUVr' = 1.34 (%RSE 0.873)
    lkin <- log(0.232); label("Zero-order SUVr production rate Kin (SUVr units/year)") # Table 2 'Kin (1/year)' = 0.232 (%RSE 11.1); Model S2 'KIN = (THETA(2)/8760)'
    lemax <- log(1.54); label("Maximum fold-stimulation of SUVr loss Emax at age 72 (unitless)") # Table 2 Emax = 1.54 (%RSE 11.8)
    lec50 <- log(75.0); label("Serum lecanemab concentration giving half of Emax (ug/mL)") # Table 2 EC50 = 75.0 ug/mL (%RSE 19.6)
    e_apoe4_carrier_rbase <- 1.04; label("Baseline SUVr ratio for APOE4 carriers vs non-carriers, applied as ratio^APOE4_CARRIER (unitless)") # Table 2 'APOE4 carrier ~ baseline' = 1.04 (%RSE 1.04)
    e_age_emax <- 1.58; label("Power exponent on (AGE/72) for Emax (unitless)") # Table 2 'Age ~ Emax (exponent)' = 1.58 (%RSE 20.9)

    # Lecanemab PK IIV, fixed at the Hayato_2022_lecanemab estimates (Table 1;
    # CV% = sqrt(variance) x 100, so omega^2 = (CV/100)^2).
    etalcl ~ fixed(0.151321) # Table 1 IIV CL 38.9%
    etalvc ~ fixed(0.0196) # Table 1 IIV V1 14.0%
    etalvp ~ fixed(0.990025) # Table 1 IIV V2 99.5%
    etalfcentral ~ fixed(0.116964) # Table 1 IIV F 34.2%; Process B doses only

    # CV% = sqrt(variance) x 100 (Table 2 footnote), so omega^2 = (CV/100)^2.
    # Covariance = R * omega_baseline * omega_Emax = 0.669 * 0.109 * 0.503.
    etalrbase + etalemax ~ c(
      0.011881,
      0.036679, 0.253009
    ) # Table 2: IIV baseline 10.9%, IIV Emax 50.3%, correlation baseline ~ Emax R = 0.669 (Model S2 $OMEGA BLOCK(2))

    propSd <- 0.0501; label("Proportional residual error on SUVr (fraction)") # Table 2 'Proportional' = 5.01% (%RSE 2.75)
  })

  model({
    hrs_per_year <- 8760 # Model S2 'KIN = (THETA(2)/8760)': 365 d x 24 h

    # Lecanemab PK driver (Hayato 2022 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / 73.7)^e_wt_cl * (ALB / 42.9)^e_alb_cl *
      e_female_cl^SEXF * e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 73.7)^e_wt_vc * e_female_vc^SEXF
    vp <- exp(lvp + etalvp) * e_japanese_vp^RACE_JAPANESE
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual SUVr parameters (Table 2 equations; Model S2 $PK)
    rbase <- exp(lrbase + etalrbase) * e_apoe4_carrier_rbase^APOE4_CARRIER
    kin <- exp(lkin) / hrs_per_year
    kout <- kin / rbase
    emax <- exp(lemax + etalemax) * (AGE / 72)^e_age_emax
    ec50 <- exp(lec50)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Model S1 $PK: 'F1=1; IF (FORM.EQ.1) F1=THETA(11)*EXP(ETA(4))'
    fprocessb <- exp(lfcentral + etalfcentral) * e_processb_f
    f(central) <- (1 - FORM_LEC_PROCESSB) + FORM_LEC_PROCESSB * fprocessb

    Cc <- central / vc

    # Figure 1 / Model S2 $DES:
    # dSUVr/dt = Kin - SUVr * Kout * (1 + Emax * C / (EC50 + C))
    d/dt(suvr) <- kin - suvr * kout * (1 + emax * Cc / (ec50 + Cc))
    suvr(0) <- rbase

    suvr ~ prop(propSd)
  })
}
