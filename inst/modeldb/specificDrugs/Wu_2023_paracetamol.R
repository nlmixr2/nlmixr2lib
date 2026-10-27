Wu_2023_paracetamol <- function() {
  description <- paste(
    "Parent-and-metabolites population PK model for intravenous paracetamol",
    "(PCM) and its glucuronide, sulfate and combined oxidative (cysteine +",
    "mercapturate) metabolites in plasma and urine of preterm neonates born",
    "below 32 weeks of gestation in their first 30 days of life (Wu 2023).",
    "Two-compartment plasma disposition for parent PCM; one-compartment",
    "disposition for each metabolite with volumes set as fractions of the PCM",
    "central volume; four cumulative urine compartments. Total PCM clearance",
    "scales with birthweight (power) and postnatal age (linear) and is",
    "partitioned into glucuronidation, sulfation, oxidation and unchanged",
    "renal excretion by three estimated clearance ratios against the",
    "unchanged-renal route (glucuronidation ratio rising with birthweight,",
    "sulfation ratio declining hyperbolically with postnatal age). The three",
    "metabolite renal clearances share one birthweight/postnatal-age",
    "relationship and one random effect. Urinary concentrations divide the",
    "excreted amount by the measured diaper urine volume times an estimated",
    "volume correction factor. Metabolite amounts and concentrations are in",
    "PCM equivalents."
  )
  reference <- paste(
    "Wu Y, Voller S, Krekels EHJ, Roofthooft DWE, Simons SHP, Tibboel D,",
    "Flint RB, Knibbe CAJ (2023).",
    "Maturation of Paracetamol Elimination Routes in Preterm Neonates Born",
    "Below 32 Weeks of Gestation.",
    "Pharm Res 40(9):2155-2166.",
    "doi:10.1007/s11095-023-03580-3.",
    sep = " "
  )
  vignette <- "Wu_2023_paracetamol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Methods 'Population Pharmacokinetics Model Development': metabolite
  # concentrations were converted to PCM equivalents with the molecular
  # weights (151.16, 327.29, 231.23, 312.24 and 270.30 g/mol), so every state
  # holds mg of PCM or mg PCM-equivalent. Supplement $MODEL names the plasma
  # and urine compartments.
  compartmentData <- list(
    central = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE),
    central_gluc = list(
      analyte = "paracetamol-glucuronide (PCM equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_sulf = list(
      analyte = "paracetamol-sulfate (PCM equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_cysmer = list(
      analyte = "paracetamol-cysteine + paracetamol-mercapturate (PCM equivalents)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    urine = list(analyte = "paracetamol", units = "mg", specimen = "urine", verified = TRUE),
    urine_gluc = list(
      analyte = "paracetamol-glucuronide (PCM equivalents)",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    ),
    urine_sulf = list(
      analyte = "paracetamol-sulfate (PCM equivalents)",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    ),
    urine_cysmer = list(
      analyte = "paracetamol-cysteine + paracetamol-mercapturate (PCM equivalents)",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT_BIRTH = list(
      description = "Body weight at birth",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling of total PCM clearance (exponent 0.902), of the",
        "glucuronidation clearance ratio (exponent 0.479) and of the three",
        "metabolite renal clearances (shared exponent 1.62), all centred on",
        "the 985 g cohort median (Table II). The paper works in grams;",
        "model() divides the kg canonical by 0.985."
      ),
      source_name = "BW"
    ),
    WT = list(
      description = "Current body weight (time-varying)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling of the PCM central volume (exponent 0.836) centred",
        "on the 1,012 g cohort median (Table II). Methods: current weight was",
        "set equal to birthweight in dataset 1 and linearly interpolated",
        "between measurements in dataset 2."
      ),
      source_name = "WT"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. The paper expresses PNA in days (data limited to",
        "0-30 days); model() recovers days with PNA * 30.4375. Linear",
        "effects on total PCM clearance (0.0469 per day) and on the",
        "metabolite renal clearances (0.0854 per day), both centred on",
        "5 days; hyperbolic decline of the sulfation clearance ratio with",
        "PNA50 = 9.18 days (Table II)."
      ),
      source_name = "PNA"
    ),
    STUDY_DINO = list(
      description = "DINO study (dataset 2) record indicator; 1 = DINO, 0 = Flint et al. dataset 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Flint et al. single-dose dataset 1)",
      notes = paste(
        "Selects the dataset-specific plasma residual-error magnitudes only.",
        "Direct copy of the supplement control-stream column ASY ('0 for",
        "dataset1 and 1 for dataset2'), used as",
        "ADD1 = (1-ASY)*THETA(22) + ASY*THETA(26) and PRO3 = (1-ASY)*THETA(31)."
      ),
      source_name = "ASY"
    ),
    URINE_VOL_INTERVAL = list(
      description = "Urine volume recovered from the diaper over the collection interval of a urinary observation",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Observation-equation divisor only (supplement $PK",
        "S5..S8 = UVOL*THETA(18)*EXP(ETA(7))). The urine concentrations are",
        "in mg/L, so the NONMEM UVOL column is in litres; model() converts",
        "the mL canonical with / 1000. Only dataset 1 collected urine (one",
        "sample per neonate over at most 12 h after the single dose)."
      ),
      source_name = "UVOL"
    )
  )

  paper_specific_etas <- c("etalcle")

  population <- list(
    species = "human",
    n_subjects = 143L,
    n_studies = 2L,
    age_range = "Postnatal age 0-30 days at start of treatment (median 5); gestational age 24.0-31.9 weeks (median 27.7); postmenstrual age 24.7-35.1 weeks (median 28.6)",
    weight_range = "Birthweight 462-1,925 g (median 985); current weight 462-1,959 g (median 1,012)",
    sex_female_pct = 47.6,
    race_ethnicity = NA_character_,
    disease_state = "Preterm neonates born below 32 weeks of gestation in Dutch level III NICUs receiving intravenous paracetamol for pain; 35% small for gestational age",
    dose_range = "Dataset 1: single 15-min IV infusion of 10, 15 or 20 mg/kg. Dataset 2 (DINO): IV loading dose median 19.0 mg/kg in 28 of 83 neonates, then median 9.5 mg/kg every ~6 h (median 10 maintenance doses)",
    regions = "Netherlands",
    ga_range = "24.0-31.9 weeks",
    notes = paste(
      "Table I. Dataset 1 (Flint et al., n = 60, rich sampling after one",
      "dose; plasma and one urine sample per neonate in 44) pooled with the",
      "DINO study (n = 83, sparse scavenged/strategic plasma sampling during",
      "multiple dosing). 534 plasma and 44 urine samples of PCM,",
      "PCM-glucuronide, PCM-sulfate, PCM-cysteine and PCM-mercapturate."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Parent PCM disposition (Table II 'PCM' block; supplement $PK TVCL1,
    # TVV1, TVV9, TVQ). Reference subject: birthweight 985 g, current
    # weight 1,012 g, PNA 5 days.
    # ---------------------------------------------------------------------
    lcl <- log(0.137); label("Total PCM clearance CLT at birthweight 985 g and PNA 5 days (L/h)") # Table II TVCLT = 0.137 (RSE 3%)
    e_wt_birth_cl <- 0.902; label("Power exponent of birthweight on total PCM clearance (unitless)") # Table II theta_BWb = 0.902 (RSE 11%)
    e_pna_cl <- 0.0469; label("Linear postnatal-age effect on total PCM clearance (1/day)") # Table II theta_PNA = 0.0469 (RSE 12%)
    lvc <- log(0.955); label("PCM central volume VP at current weight 1,012 g (L)") # Table II TVVP = 0.955 (RSE 3%)
    e_wt_vc <- 0.836; label("Power exponent of current weight on PCM central volume (unitless)") # Table II theta_BWc = 0.836 (RSE 11%)
    lq <- log(0.0494); label("PCM inter-compartmental clearance Qp (L/h)") # Table II Qp = 0.0494 (RSE 19%)
    lvp <- log(0.244); label("PCM peripheral volume VPP (L)") # Table II VPP = 0.244 (RSE 12%)

    # ---------------------------------------------------------------------
    # Metabolite volumes as fractions of the PCM central volume (Table II
    # 'Volume of metabolites'; supplement FRAC1..FRAC3).
    # ---------------------------------------------------------------------
    lf_vol_gluc <- log(0.6); label("PCM-glucuronide volume as a fraction of PCM central volume, fVG (unitless)") # Table II fVG = 0.6 (RSE 9%)
    lf_vol_sulf <- log(0.326); label("PCM-sulfate volume as a fraction of PCM central volume, fVS (unitless)") # Table II fVS = 0.326 (RSE 4%)
    lf_vol_cysmer <- log(0.807); label("Oxidative-metabolite volume as a fraction of PCM central volume, fVO (unitless)") # Table II fVO = 0.807 (RSE 6%)

    # ---------------------------------------------------------------------
    # Metabolite renal clearances (Table II 'Renal clearance of
    # metabolites'). Separate typical values, one shared birthweight /
    # PNA relationship (supplement THETA(11), THETA(12)) and one shared
    # random effect (supplement ETA(3) on CL2, CL4 and CL6).
    # ---------------------------------------------------------------------
    lcle_gluc <- log(0.0425); label("PCM-glucuronide renal clearance CLrG at birthweight 985 g and PNA 5 days (L/h)") # Table II TVCLrG = 0.0425 (RSE 5%)
    lcle_sulf <- log(0.0256); label("PCM-sulfate renal clearance CLrS at birthweight 985 g and PNA 5 days (L/h)") # Table II TVCLrS = 0.0256 (RSE 4%)
    lcle_cysmer <- log(0.0498); label("Oxidative-metabolite renal clearance CLrO at birthweight 985 g and PNA 5 days (L/h)") # Table II TVCLrO = 0.0498 (RSE 5%)
    e_wt_birth_cle <- 1.62; label("Power exponent of birthweight on the metabolite renal clearances (unitless)") # Table II renal block theta_BWb = 1.62 (RSE 9%)
    e_pna_cle <- 0.0854; label("Linear postnatal-age effect on the metabolite renal clearances (1/day)") # Table II renal block theta_PNA = 0.0854 (RSE 12%)

    # ---------------------------------------------------------------------
    # Pathway partition (Methods Eqs. 1-5; Table II 'Formation clearance
    # of metabolites'). Each ratio is that route's clearance divided by
    # the unchanged-renal PCM clearance (the paper's R_FG / R_FS / R_FO);
    # the four fractions R/(R_FG + R_FS + R_FO + 1) sum to one.
    # ---------------------------------------------------------------------
    lclrat_gluc <- log(1.16); label("Glucuronidation-to-unchanged-renal clearance ratio TVR_FG at birthweight 985 g (unitless)") # Table II TVR_FG = 1.16 (RSE 9%)
    e_wt_birth_clrat_gluc <- 0.479; label("Power exponent of birthweight on the glucuronidation clearance ratio (unitless)") # Table II R_FG theta_BWb = 0.479 (RSE 45%)
    lclrat_sulf <- log(24.5); label("Sulfation-to-unchanged-renal clearance ratio at birth TVR_FS (unitless)") # Table II TVR_FS = 24.5 (RSE 11%)
    lpna50_clrat_sulf <- log(9.18); label("Postnatal age at which the sulfation ratio has halved, PNA50 (days)") # Table II theta_PNA50 = 9.18 (RSE 23%)
    lclrat_cysmer <- log(3.51); label("Oxidation-to-unchanged-renal clearance ratio R_FO (unitless)") # Table II R_FO = 3.51 (RSE 8%)

    # ---------------------------------------------------------------------
    # Urine volume correction factor (Results; Table II). The measured
    # diaper volume is multiplied by this factor to recover the excreted
    # amounts.
    # ---------------------------------------------------------------------
    lcorr_urine_vol <- log(3.62); label("Correction factor on the measured urine collection volume (unitless)") # Table II urine volume correction factor = 3.62 (RSE 16%)

    # ---------------------------------------------------------------------
    # IIV. Table II prints CV% = sqrt(exp(omega^2) - 1) * 100 (footnote),
    # so omega^2 = log(1 + CV^2).
    # ---------------------------------------------------------------------
    etalcl ~ 0.09118 # Table II IIV CLT 30.9% -> log(1 + 0.309^2)
    etalvc ~ 0.07704 # Table II IIV VP 28.3% -> log(1 + 0.283^2)
    etalcle ~ 0.13207 # Table II IIV CLrG = CLrS = CLrO 37.6% (one shared eta) -> log(1 + 0.376^2)
    etalclrat_gluc ~ 0.29092 # Table II IIV TVR_FG 58.1% -> log(1 + 0.581^2)
    etalclrat_sulf ~ 0.18069 # Table II IIV TVR_FS 44.5% -> log(1 + 0.445^2)
    etalclrat_cysmer ~ 0.33137 # Table II IIV R_FO 62.7% -> log(1 + 0.627^2)
    etalcorr_urine_vol ~ 1.07003 # Table II inter-sample variability of the urine correction factor 138.4% -> log(1 + 1.384^2)

    # ---------------------------------------------------------------------
    # Residual error on log-transformed concentrations (Methods Eq. 7;
    # supplement $ERROR): log(C) = log(IPRED) + SD * eps with
    # SD = sqrt(SD1^2 + (SD2 / IPRED)^2). SD1 is the log-scale (lnorm)
    # SD; SD2 is in mg/L and was retained only for the dataset-1
    # sulfate and oxidative-metabolite plasma data. Dataset 1 is the
    # reference (STUDY_DINO = 0); the _dino terms apply to dataset 2.
    # ---------------------------------------------------------------------
    expSd <- 0.166; label("Plasma PCM log-scale residual SD, dataset 1") # Table II Plasma Dataset 1 SD1 PCM = 0.166 (RSE 6%)
    expSd_gluc <- 0.325; label("Plasma PCM-glucuronide log-scale residual SD, dataset 1") # Table II Plasma Dataset 1 SD1 PCM-glucuronide = 0.325 (RSE 6%)
    expSd_sulf <- 0.168; label("Plasma PCM-sulfate log-scale residual SD, dataset 1") # Table II Plasma Dataset 1 SD1 PCM-sulfate = 0.168 (RSE 10%)
    expSd_cysmer <- 0.208; label("Plasma oxidative-metabolite log-scale residual SD, dataset 1") # Table II Plasma Dataset 1 SD1 PCM-oxidative pathway metabolites = 0.208 (RSE 11%)
    expSd_dino <- 0.345; label("Plasma PCM log-scale residual SD, dataset 2 (DINO)") # Table II Plasma Dataset 2 SD1 PCM = 0.345 (RSE 5%)
    expSd_gluc_dino <- 0.858; label("Plasma PCM-glucuronide log-scale residual SD, dataset 2 (DINO)") # Table II Plasma Dataset 2 SD1 PCM-glucuronide = 0.858 (RSE 5%)
    expSd_sulf_dino <- 0.325; label("Plasma PCM-sulfate log-scale residual SD, dataset 2 (DINO)") # Table II Plasma Dataset 2 SD1 PCM-sulfate = 0.325 (RSE 6%)
    expSd_cysmer_dino <- 0.369; label("Plasma oxidative-metabolite log-scale residual SD, dataset 2 (DINO)") # Table II Plasma Dataset 2 SD1 PCM-oxidative metabolites = 0.369 (RSE 5%)
    addSd_sulf <- 0.464; label("Plasma PCM-sulfate additive residual SD2, dataset 1 (mg/L)") # Table II Plasma Dataset 1 SD2 PCM-sulfate = 0.464 (RSE 16%)
    addSd_cysmer <- 0.0489; label("Plasma oxidative-metabolite additive residual SD2, dataset 1 (mg/L)") # Table II Plasma Dataset 1 SD2 PCM-oxidative metabolites = 0.0489 (RSE 22%)
    expSd_urine <- 0.115; label("Urine log-scale residual SD shared by PCM and all metabolites, dataset 1") # Table II Urine Dataset 1 SD1 combined error = 0.115 (RSE 26%)
  })

  model({
    # Covariates in the paper's units: grams and days.
    pna_days <- PNA * 30.4375
    bw_ratio <- WT_BIRTH / 0.985
    wt_ratio <- WT / 1.012

    # Total PCM clearance and PCM disposition (Table II).
    cl <- exp(lcl + etalcl) * bw_ratio^e_wt_birth_cl * (1 + (pna_days - 5) * e_pna_cl)
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # Pathway partition (Methods Eqs. 1-5). The sulfation ratio declines
    # hyperbolically with PNA: R_FS = TVR_FS * (1 - PNA / (PNA50 + PNA)).
    clrat_gluc <- exp(lclrat_gluc + etalclrat_gluc) * bw_ratio^e_wt_birth_clrat_gluc
    pna50_clrat_sulf <- exp(lpna50_clrat_sulf)
    clrat_sulf <- exp(lclrat_sulf + etalclrat_sulf) * (1 - pna_days / (pna50_clrat_sulf + pna_days))
    clrat_cysmer <- exp(lclrat_cysmer + etalclrat_cysmer)
    clrat_total <- clrat_gluc + clrat_sulf + clrat_cysmer + 1
    fm_gluc <- clrat_gluc / clrat_total
    fm_sulf <- clrat_sulf / clrat_total
    fm_cysmer <- clrat_cysmer / clrat_total
    fe <- 1 / clrat_total

    cl_gluc <- fm_gluc * cl
    cl_sulf <- fm_sulf * cl
    cl_cysmer <- fm_cysmer * cl
    cl_renal <- fe * cl

    # Metabolite volumes and renal clearances; one random effect is shared
    # by the three renal clearances (supplement ETA(3)).
    vc_gluc <- exp(lf_vol_gluc) * vc
    vc_sulf <- exp(lf_vol_sulf) * vc
    vc_cysmer <- exp(lf_vol_cysmer) * vc
    cle_mat <- exp(etalcle) * bw_ratio^e_wt_birth_cle * (1 + (pna_days - 5) * e_pna_cle)
    cle_gluc <- exp(lcle_gluc) * cle_mat
    cle_sulf <- exp(lcle_sulf) * cle_mat
    cle_cysmer <- exp(lcle_cysmer) * cle_mat

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_gluc) <- cl_gluc / vc * central - cle_gluc / vc_gluc * central_gluc
    d/dt(central_sulf) <- cl_sulf / vc * central - cle_sulf / vc_sulf * central_sulf
    d/dt(central_cysmer) <- cl_cysmer / vc * central - cle_cysmer / vc_cysmer * central_cysmer
    d/dt(urine) <- cl_renal / vc * central
    d/dt(urine_gluc) <- cle_gluc / vc_gluc * central_gluc
    d/dt(urine_sulf) <- cle_sulf / vc_sulf * central_sulf
    d/dt(urine_cysmer) <- cle_cysmer / vc_cysmer * central_cysmer

    # Plasma concentrations (mg/L; metabolites in PCM equivalents).
    Cc <- central / vc
    Cc_gluc <- central_gluc / vc_gluc
    Cc_sulf <- central_sulf / vc_sulf
    Cc_cysmer <- central_cysmer / vc_cysmer

    # Urine concentrations: excreted amount over the corrected collection
    # volume (supplement S5..S8 = UVOL * THETA(18) * EXP(ETA(7))). The
    # urine states are cumulative; reset them at the start of each
    # collection interval.
    corr_urine_vol <- exp(lcorr_urine_vol + etalcorr_urine_vol)
    urine_volume <- max(URINE_VOL_INTERVAL / 1000, 0.001) * corr_urine_vol
    Curine <- urine / urine_volume
    Curine_gluc <- urine_gluc / urine_volume
    Curine_sulf <- urine_sulf / urine_volume
    Curine_cysmer <- urine_cysmer / urine_volume

    # Residual error (Methods Eq. 7). The SD2 / IPRED term is floored at
    # the supplement's 1e-6 mg/L IPRED floor so a zero metabolite
    # concentration at time zero does not divide by zero.
    expSd_i <- expSd * (1 - STUDY_DINO) + expSd_dino * STUDY_DINO
    expSd_gluc_i <- expSd_gluc * (1 - STUDY_DINO) + expSd_gluc_dino * STUDY_DINO
    expSd_sulf_i <- sqrt(
      (expSd_sulf * (1 - STUDY_DINO) + expSd_sulf_dino * STUDY_DINO)^2 +
        (addSd_sulf * (1 - STUDY_DINO) / max(Cc_sulf, 1e-6))^2
    )
    expSd_cysmer_i <- sqrt(
      (expSd_cysmer * (1 - STUDY_DINO) + expSd_cysmer_dino * STUDY_DINO)^2 +
        (addSd_cysmer * (1 - STUDY_DINO) / max(Cc_cysmer, 1e-6))^2
    )

    # One urine SD is shared by the four urinary endpoints (supplement
    # ADD5..ADD8 = THETA(30)); rxode2 needs a distinct name per endpoint.
    expSd_urine_i <- expSd_urine
    expSd_urine_gluc_i <- expSd_urine
    expSd_urine_sulf_i <- expSd_urine
    expSd_urine_cysmer_i <- expSd_urine

    Cc ~ lnorm(expSd_i)
    Cc_gluc ~ lnorm(expSd_gluc_i)
    Cc_sulf ~ lnorm(expSd_sulf_i)
    Cc_cysmer ~ lnorm(expSd_cysmer_i)
    Curine ~ lnorm(expSd_urine_i)
    Curine_gluc ~ lnorm(expSd_urine_gluc_i)
    Curine_sulf ~ lnorm(expSd_urine_sulf_i)
    Curine_cysmer ~ lnorm(expSd_urine_cysmer_i)
  })
}
