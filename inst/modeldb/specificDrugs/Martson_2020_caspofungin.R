Martson_2020_caspofungin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous caspofungin in critically ill adults",
    "in the intensive care unit with (suspected) invasive candidiasis (Groningen / Enschede, the",
    "Netherlands). Fitted non-parametrically with the NPAG algorithm in Pmetrics 1.5.2 on",
    "primary parameters Ke (elimination rate constant), V (central volume) and the",
    "intercompartmental rate constants kcp and kpc. Body weight is the only retained covariate",
    "and scales the central volume linearly, V = V0 * WT / 78 with 78 kg the cohort median; since",
    "Ke carries no covariate, clearance CL = Ke * V is also proportional to weight, which is the",
    "basis of the paper's weight-based (2 mg/kg loading, 1.25 mg/kg maintenance) dosing",
    "recommendation. Typical values are the means of the NPAG marginal distributions (Table 2);",
    "inter-individual variability is a log-normal approximation built from the Table 2 CV%",
    "column. Table 2 prints V0 in 'liters/kg', but the equation V = V0 * weight/78 and the",
    "Discussion's CL ~ 0.7 L/h (= 0.09 x 7.71) show V0 is the volume in litres of a 78 kg",
    "patient; see the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Martson A-G, van der Elst KCM, Veringa A, Zijlstra JG, Beishuizen A, van der Werf TS,",
    "Kosterink JGW, Neely M, Alffenaar J-W. Caspofungin weight-based dosing supported by a",
    "population pharmacokinetic model in critically ill patients. Antimicrob Agents Chemother.",
    "2020;64(9):e00905-20. doi:10.1128/AAC.00905-20. PMCID: PMC7449215.",
    sep = " "
  )
  vignette <- "Martson_2020_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The ONLY covariate retained in the final model. Enters the central volume as the linear",
        "ratio V = V0 * weight/78 (Results, 'Population pharmacokinetic model'), with 78 kg the",
        "cohort median weight (Table 1, median 78 kg, range 48-139 kg); the exponent is",
        "structurally 1, not estimated. Time-fixed in the source (one weight per patient).",
        sep = " "
      ),
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Tested with forward addition (Methods, 'Population pharmacokinetic modeling') but not retained in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Tested ('gender') but not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Tested but not retained; the Discussion attributes this partly to albumin being",
        "measured infrequently. Cohort median 20 g/L (range 14-28; Table 1).",
        sep = " "
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested but not retained. Table 1 prints the unit as 'mmol/liter' (median 7.5, range 3-376), which can only be umol/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained (Discussion: did not improve the population goodness of fit)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained (Discussion: did not improve the population goodness of fit)."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained (Discussion: did not improve the population goodness of fit)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Tested but not retained."
    ),
    WBC = list(
      description = "Leukocyte count",
      units = "10^9/L",
      type = "continuous",
      notes = "Tested ('leukocyte count') but not retained; no cohort summary reported."
    ),
    SAPS3 = list(
      description = "Simplified Acute Physiology Score 3",
      units = "points",
      type = "continuous",
      notes = "Tested but not retained. Cohort median 59 (range 31-104; Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal-replacement-therapy treatment-status indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Tested ('dialysis' / 'hemodialysis') but not retained; 8 of 20 patients (40%) were on continuous venovenous hemofiltration (Table 1)."
    ),
    CONMED_CORTICOSTEROID = list(
      description = "Concomitant prednisolone or hydrocortisone",
      units = "(binary)",
      type = "categorical",
      notes = "Tested but not retained; 11 of 20 patients (55%) received it (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_occasions = 25L,
    n_concentrations = 219L,
    age_range = "25-83 years",
    age_median = "56 years",
    weight_range = "48-139 kg",
    weight_median = "78 kg",
    sex_female_pct = 45,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Adult critically ill ICU patients treated with caspofungin for (suspected) invasive",
      "candidiasis. Median SAPS 3 59 (31-104); 40% on continuous venovenous hemofiltration; 55%",
      "co-administered prednisolone or hydrocortisone; two patients had Child-Pugh C liver",
      "impairment. Median serum albumin 20 g/L (14-28).",
      sep = " "
    ),
    dose_range = paste(
      "70 mg loading dose on day 1, then 50 mg daily (<= 80 kg) or 70 mg daily (> 80 kg);",
      "35 mg / 50 mg with moderate hepatic impairment (Child-Pugh 7-9). Each dose a 1 h",
      "intravenous infusion. Doses were changed (TDM) when AUC0-24 was below 98 mg*h/L.",
      sep = " "
    ),
    regions = "The Netherlands (University Medical Center Groningen; Medisch Spectrum Twente, Enschede)",
    notes = paste(
      "Demographics from Martson 2020 Table 1 (reproduced from the parent TDM study). Nine",
      "samples per dosing occasion (pre-dose and 1, 2, 3, 4, 6, 8, 12 and 24 h after the start",
      "of infusion) on day 3 (range 2-4) of treatment; 15 patients contributed one occasion and",
      "5 patients two occasions, 219 concentrations in total. Total plasma caspofungin by",
      "LC-MS/MS.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural values are the MEANS of the NPAG non-parametric marginal
    # distributions in Martson 2020 Table 2 (medians are lower: Ke 0.08, V0 7.20,
    # kcp 0.28, kpc 0.34). The primary parameters are Ke and V (not CL).
    # ------------------------------------------------------------------------
    lkel <- log(0.09); label("Elimination rate constant Ke (1/h)")
    # Table 2, Ke mean = 0.09 1/h (SD 0.04, median 0.08, CV 42.38%); also Abstract and Results text

    lvc <- log(7.71); label("Central volume V0 for a 78 kg patient (L)")
    # Table 2, V0 mean = 7.71 (SD 2.70, median 7.20, CV 34.98%). Table 2 labels the unit 'liters/kg',
    # but the model equation is V = V0 * weight/78 (Results), so V0 is the volume (L) at the 78 kg
    # median; the Discussion confirms 'with a Ke of 0.09, our CL is approximately 0.7' = 0.09 * 7.71.

    lk12 <- log(0.44); label("Transfer rate constant central -> peripheral1, kcp (1/h)")
    # Table 2, kcp mean = 0.44 1/h (SD 0.38, median 0.28, CV 88.02%); Abstract prints SD 0.39

    lk21 <- log(0.46); label("Transfer rate constant peripheral1 -> central, kpc (1/h)")
    # Table 2, kpc mean = 0.46 1/h (SD 0.35, median 0.34, CV 75.98%)

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete distribution, not
    # an omega; the Table 2 CV% column is carried as a LOG-NORMAL approximation,
    # omega^2 = log(CV^2 + 1). Correlations are not reported and are absent. The
    # support-point table was not published. See vignette Assumptions.
    # ------------------------------------------------------------------------
    etalkel ~ 0.165181 # Table 2 Ke CV 42.38% -> log(0.4238^2 + 1)
    etalvc ~ 0.115434 # Table 2 V0 CV 34.98% -> log(0.3498^2 + 1)
    etalk12 ~ 0.573661 # Table 2 kcp CV 88.02% -> log(0.8802^2 + 1)
    etalk21 ~ 0.455712 # Table 2 kpc CV 75.98% -> log(0.7598^2 + 1)

    # ------------------------------------------------------------------------
    # Residual error. Pmetrics multiplicative (gamma) model on the assay-error
    # polynomial SD = C0 + C1 * C with C0 = 0.05 and C1 = 0.08 (Methods, 'Model
    # diagnostics') and converged gamma = 0.654 (Results). Total SD =
    # gamma * C0 + gamma * C1 * C, a LINEAR sum -> combined1().
    # ------------------------------------------------------------------------
    addSd <- 0.0327; label("Additive residual SD (mg/L)")
    # gamma * C0 = 0.654 * 0.05 = 0.0327 (Methods C0 = 0.05; Results gamma = 0.654)
    propSd <- 0.05232; label("Proportional residual SD (fraction)")
    # gamma * C1 = 0.654 * 0.08 = 0.05232 (Methods C1 = 0.08; Results gamma = 0.654)
  })

  model({
    # 1. Individual primary parameters. Weight scales the central volume only,
    #    as the linear ratio V = V0 * weight/78 (Results).
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc) * (WT / 78)
    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # 2. Derived clearances and peripheral volume (for reporting). The ODEs are
    #    written with the micro-constants below.
    cl <- kel * vc
    q <- kcp * vc
    vp <- q / kpc

    # 3. Micro-constants.
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment ODE system; caspofungin is given as a 1 h IV infusion
    #    into central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation (total plasma caspofungin, mg/L) and Pmetrics error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
