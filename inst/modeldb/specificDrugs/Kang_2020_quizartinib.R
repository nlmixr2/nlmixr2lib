Kang_2020_quizartinib <- function() {
  description <- paste(
    "Parent-metabolite population PK model for oral quizartinib and its",
    "active metabolite AC886 in adult healthy volunteers and adults with",
    "relapsed/refractory FLT3-ITD acute myeloid leukemia (AML), pooled",
    "across 8 phase 1-3 studies including QuANTUM-R (Kang 2020).",
    "Quizartinib is described by a three-compartment model with sequential",
    "zero-order (duration D1) then first-order (ka) absorption from a depot",
    "and an absorption lag time; AC886 is a two-compartment model fed by a",
    "fixed fraction (fMET = 0.5) of quizartinib clearance. Covariates:",
    "AML patient status on ka, CL and relative bioavailability; fed status",
    "on ka and bioavailability (food-effect study AC220-019 only); strong",
    "CYP3A inhibitor use on quizartinib CL and bioavailability and on AC886",
    "CL and central volume; albumin and body surface area on quizartinib",
    "Vc; body surface area and age on Vp1; body surface area on Q1; body",
    "surface area and Black race on AC886 CL. Interindividual variability",
    "on quizartinib CL, Vc, Q1, Vp1, ka, D1, F1 and lag time, and on AC886",
    "Vc and CL (separate CL variances for AML patients and healthy",
    "volunteers); interoccasion variability on quizartinib F1, CL and Vc.",
    "Residual error is combined additive + proportional for quizartinib and",
    "additive on the log scale for AC886, each with separate magnitudes for",
    "healthy volunteers and AML patients."
  )
  reference <- paste(
    "Kang D, Ludwig E, Jaworowicz D, Huang H, Fiedler-Kelly J, Cortes J,",
    "Ganguly S, Khaled S, Kramer A, Levis M, Martinelli G, Perl A,",
    "Russell N, Abutarif M, Choi Y, Mendell J, Yin O. Population",
    "pharmacokinetic analysis of quizartinib in healthy volunteers and",
    "patients with relapsed/refractory acute myeloid leukemia.",
    "J Clin Pharmacol. 2020;60(12):1629-1641. doi:10.1002/jcph.1680"
  )
  vignette <- "Kang_2020_quizartinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  paper_specific_residual_sds <- c("propSd_aml", "expSd_ac886_aml")

  # Amounts are carried in mg of the dosed quizartinib dihydrochloride
  # (Kang 2020 states every dose in salt form: "20, 30, 60, and 90 mg in
  # salt form (equivalent to 17.7, 26.5, 53.0, and 79.5 mg free base)").
  # AC886 amounts are the mass-equivalent formed from fMET * CL * Cc (Kang
  # 2020 Figure 2); the paper applies no molecular-weight conversion.
  compartmentData <- list(
    depot = list(analyte = "quizartinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_ac886 = list(analyte = "AC886", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ac886 = list(analyte = "AC886", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_AML = list(
      description = "Patient-status indicator, 1 = patient with relapsed/refractory AML, 0 = healthy volunteer.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer; 325/649 = 50.1% per Table 2). Every AML subject in this analysis is a patient and every non-AML subject is a healthy volunteer, so this is the paper's flagPATIENT directly.",
      notes = "Kang 2020 equations 1, 2 and 4: ka is 1.68 1/h in patients instead of 0.874 1/h (a level, not a fractional change), relative F1 is 0.599 in patients vs 1 in healthy volunteers, and CL is 8.2% higher in patients. The indicator also selects the population-specific AC886 CL IIV variance and the population-specific residual-error magnitudes of both analytes. Time-fixed per subject.",
      source_name = "flagPATIENT"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A inhibitor indicator, 1 = a strong CYP3A inhibitor is coadministered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (weak or no CYP3A inhibitor, or moderate inhibitor; strong inhibitor use in 122/649 = 18.8% per Table 2).",
      notes = "Time-varying (Kang 2020 Methods: concomitant medications 'were evaluated as time-varying covariates'). Equations 2, 4, 8 and 9: +13.6% on quizartinib F1, -44.1% on quizartinib CL (additive with the patient term inside one bracket, equation 4), +10.6% on AC886 CL (additive with the Black-race term, equation 9) and +192% on AC886 Vc (equation 8). Moderate CYP3A inhibitors were not retained.",
      source_name = "flagSTRCYP"
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Kang 2020 equation 5 is calibrated in g/dL: Vc = 194 * (ALB/4.1)^-0.725 * (BSA/1.9)^1.46, centred on the pooled median 4.1 g/dL (Table 2). The canonical column is g/L, so model() converts with ALB * 0.1 before applying the published centring value.",
      source_name = "ALB"
    ),
    BSA = list(
      description = "Baseline body surface area.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects centred on the pooled median 1.9 m^2 (Table 2): exponent 1.46 on quizartinib Vc (eq 5), 1.50 on Vp1 (eq 6), 0.970 on Q1 (eq 7) and 1.60 on AC886 CL (eq 9). BSA was chosen over body weight because it was the more significant body-size measure in univariate screening; replacing BSA with weight raised the objective function by 12 points (Discussion). The BSA formula is not stated.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Baseline age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on quizartinib Vp1 only, centred on the pooled median 44 years (eq 6: Vp1 = 170 * (BSA/1.9)^1.50 * (AGE/44)^0.453).",
      source_name = "Age"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator, 1 = Black or African American, 0 = any other race.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other races; Black or African American 117/649 = 18.0% per Table 2).",
      notes = "Fractional change +0.586 on AC886 CL (eq 9, additive with the strong-CYP3A-inhibitor term inside one bracket). No race effect on quizartinib.",
      source_name = "flagRACB"
    ),
    FED = list(
      description = "Fed-state indicator for the dose, 1 = dosed after a meal, 0 = fasted.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted).",
      notes = "Only defined for the food-effect study AC220-019; Kang 2020 Results: 'food effect was estimated only for individuals from the AC220-019 study' because food data were not collected comprehensively elsewhere. Set FED = 0 for every record outside AC220-019. Effects: -51.2% on ka in healthy volunteers (eq 1) and +5.09% on the study-019 F1 (eq 3).",
      source_name = "flagFED"
    ),
    STUDY_AC220019 = list(
      description = "Study AC220-019 indicator (phase 1 food-effect study, 64 healthy volunteers, single 30 mg dose), 1 = record from AC220-019, 0 = any other study.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the other 7 studies).",
      notes = "Selects the study-specific relative bioavailability of equation 3, F1 = 0.913 * (1 + 0.0509 * FED), in place of equation 2. Kang 2020 defines equation 2 as 'the F1 in the ith participant (not including participants in AC220-019)'.",
      source_name = "study == AC220-019"
    ),
    OCC = list(
      description = "Occasion index for the interoccasion variability on quizartinib F1, CL and Vc.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Kang 2020 Table 3 reports IOV on occasions 1-4 for CL and Vc and on occasions 1 and 3 for F1, but does not define the occasions. Pass OCC = 1-4 to select the occasion; OCC = 0 (or any other value) switches all IOV terms off. A single-occasion simulation may use OCC = 1 throughout.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not significant (Kang 2020 Discussion)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Tested and not significant, including a sensitivity analysis with eGFR capped at 120 mL/min/1.73 m^2 (Discussion)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Liver-function tests (AST, ALT, ALP, total bilirubin) were tested and not significant (Discussion)."
    ),
    CONMED_PPI = list(
      description = "Gastric acid-reducing agent use (proton pump inhibitors, H2 antagonists, antacids).",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not significant (Discussion), consistent with the lansoprazole interaction study AC220-018."
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate CYP3A inhibitor indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Tabulated in Table 2 (136/649 = 21.0%) but no moderate-inhibitor effect was retained in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 649L,
    n_studies = 8L,
    age_range = "18-81 years (median 44; healthy volunteers median 33, AML patients median 55)",
    age_median = "44 years",
    weight_range = "39.5-153 kg (median 74.4)",
    weight_median = "74.4 kg",
    bsa_median = "1.9 m^2 (range 1.3-2.8)",
    sex_female_pct = 42.1,
    race_ethnicity = c(
      White = 68.7,
      Black = 18.0,
      Asian = 4.8,
      `American Indian or Alaska Native` = 1.1,
      `Native Hawaiian or Pacific Islander` = 0.2,
      Other = 2.8,
      Unknown = 4.5
    ),
    disease_state = "325 healthy volunteers (5 phase 1 single-dose studies) and 324 patients with relapsed/refractory AML (FLT3-ITD-positive in QuANTUM-R), including 239 from the phase 3 QuANTUM-R study.",
    dose_range = "Single oral doses of 30-90 mg (healthy volunteers) or once-daily doses of 20-90 mg (AML patients) of quizartinib dihydrochloride, as tablets or oral solution.",
    regions = "Multinational.",
    albumin = "Median 4.4 g/dL in healthy volunteers and 3.7 g/dL in AML patients (pooled 4.1, range 2.1-5.2).",
    co_medication = "Strong CYP3A inhibitors 122/649 (18.8%), moderate 136/649 (21.0%), weak or none 391/649 (60.2%).",
    n_observations = "11,770 quizartinib and 10,888 AC886 plasma concentrations; LLOQ 2 ng/mL (5 studies) or 0.5 ng/mL (3 studies).",
    notes = "Demographics from Kang 2020 Table 2; study inventory from Table 1. Fitted in NONMEM 7.3 with FOCE. The AC886 model was fitted sequentially with individual quizartinib parameters fixed to their post hoc estimates; this file expresses both moieties as one joint model so they can be simulated together."
  )

  ini({
    # ------------------------------------------------------------------
    # QUIZARTINIB structural parameters -- Kang 2020 Table 3 and eqs 1-7
    # Reference: healthy volunteer, fasted, no strong CYP3A inhibitor,
    # ALB 4.1 g/dL, BSA 1.9 m^2, age 44 years.
    # ------------------------------------------------------------------
    lcl <- log(2.77)
    label("Quizartinib apparent clearance CL/F, healthy volunteer (L/h)") # Table 3: 2.77 L/h, RSE 2.04%; eq 4
    lvc <- log(194)
    label("Quizartinib apparent central volume Vc/F at ALB 4.1 g/dL, BSA 1.9 m^2 (L)") # Table 3: 194 L, RSE 1.89%; eq 5
    lq <- log(27.9)
    label("Quizartinib apparent intercompartmental clearance Q1/F at BSA 1.9 m^2 (L/h)") # Table 3: 27.9 L/h, RSE 3.14%; eq 7
    lvp <- log(170)
    label("Quizartinib apparent peripheral volume Vp1/F at BSA 1.9 m^2, age 44 y (L)") # Table 3: 170 L, RSE 2.53%; eq 6
    lq2 <- log(0.567)
    label("Quizartinib apparent intercompartmental clearance Q2/F (L/h)") # Table 3: 0.567 L/h, RSE 4.87%
    lvp2 <- log(39.3)
    label("Quizartinib apparent second peripheral volume Vp2/F (L)") # Table 3: 39.3 L, RSE 1.86%
    lka <- log(0.874)
    label("Quizartinib first-order absorption rate constant ka, fasted healthy volunteer (1/h)") # Table 3: 0.874 1/h, RSE 3.00%; eq 1
    lka_aml <- log(1.68)
    label("Quizartinib first-order absorption rate constant ka, AML patient (1/h)") # Table 3: 1.68 1/h, RSE 0.484%; eq 1
    ld1 <- log(0.708)
    label("Quizartinib duration of zero-order input to the depot D1 (h)") # Table 3: 0.708 h, RSE 3.86%
    ltlag <- log(0.205)
    label("Quizartinib absorption lag time ALAG1 (h)") # Table 3: 0.205 h, RSE 4.09%
    lfdepot <- fixed(log(1))
    label("Quizartinib relative bioavailability F1, healthy volunteer reference (unitless)") # eq 2: F1 = 1 * (1 - flagPATIENT) + ...; reference level fixed at 1 (no IV data)

    # Quizartinib covariate effects
    e_dis_aml_fdepot <- 0.599
    label("Quizartinib relative F1 for AML patients vs healthy volunteers (unitless multiplier)") # Table 3 'F1 (relative F1 for patients)' = 0.599, RSE 3.26%; eq 2
    e_cyp3a4_inh_strong_fdepot <- 0.136
    label("Fractional change in quizartinib F1 with strong CYP3A inhibitors (unitless)") # Table 3: 0.136, RSE 14.6%; eq 2
    lfdepot_ac220019 <- log(0.913)
    label("Quizartinib relative F1 in study AC220-019, fasted (unitless)") # Table 3 'F1 (relative F1 for study 019 [fasted])' = 0.913, RSE 1.86%; eq 3
    e_fed_fdepot <- 0.0509
    label("Fractional change in study AC220-019 F1 when fed (unitless)") # Table 3: 0.0509, RSE 52.0%; eq 3
    e_fed_ka <- -0.512
    label("Fractional change in healthy-volunteer ka when fed (unitless)") # Table 3: -0.512, RSE 8.75%; eq 1
    e_dis_aml_cl <- 0.0820
    label("Fractional change in quizartinib CL for AML patients (unitless)") # Table 3: 0.0820, RSE 29.1%; eq 4 prints 0.082
    e_cyp3a4_inh_strong_cl <- -0.441
    label("Fractional change in quizartinib CL with strong CYP3A inhibitors (unitless)") # Table 3: -0.441, RSE 3.86%; eq 4
    e_alb_vc <- -0.725
    label("Power exponent of (ALB/4.1 g/dL) on quizartinib Vc (unitless)") # Table 3: -0.725, RSE 12.6%; eq 5
    e_bsa_vc <- 1.46
    label("Power exponent of (BSA/1.9 m^2) on quizartinib Vc (unitless)") # Table 3: 1.46, RSE 7.16%; eq 5
    e_bsa_q <- 0.970
    label("Power exponent of (BSA/1.9 m^2) on quizartinib Q1 (unitless)") # Table 3: 0.970, RSE 23.3%; eq 7
    e_bsa_vp <- 1.50
    label("Power exponent of (BSA/1.9 m^2) on quizartinib Vp1 (unitless)") # Table 3: 1.50, RSE 10.9%; eq 6
    e_age_vp <- 0.453
    label("Power exponent of (AGE/44 years) on quizartinib Vp1 (unitless)") # Table 3: 0.453, RSE 10.7%; eq 6

    # ------------------------------------------------------------------
    # AC886 structural parameters -- Kang 2020 Table 4 and eqs 8-9
    # ------------------------------------------------------------------
    lcl_ac886 <- log(4.09)
    label("AC886 apparent clearance CLm/F at BSA 1.9 m^2 (L/h)") # Table 4: 4.09 L/h, RSE 3.07%; eq 9
    lvc_ac886 <- log(4.95)
    label("AC886 apparent central volume Vcm/F (L)") # Table 4: 4.95 L, RSE 8.27%; eq 8
    lvp_ac886 <- log(70.6)
    label("AC886 apparent peripheral volume Vpm/F (L)") # Table 4: 70.6 L, RSE 1.06%
    lq_ac886 <- log(3.29)
    label("AC886 apparent intercompartmental clearance Qm/F (L/h)") # Table 4: 3.29 L/h, RSE 1.28%
    fmet <- fixed(0.5)
    label("Fraction of quizartinib clearance converted to AC886 fMET (unitless)") # Table 4: 0.500 FIXED; Results: geometric mean AC886/quizartinib AUC ratio in study 2689-CL-2004

    e_bsa_cl_ac886 <- 1.60
    label("Power exponent of (BSA/1.9 m^2) on AC886 CL (unitless)") # Table 4: 1.60, RSE 12.5%; eq 9
    e_race_black_cl_ac886 <- 0.586
    label("Fractional change in AC886 CL for Black or African American race (unitless)") # Table 4: 0.586, RSE 14.1%; eq 9
    e_cyp3a4_inh_strong_cl_ac886 <- 0.106
    label("Fractional change in AC886 CL with strong CYP3A inhibitors (unitless)") # Table 4: 0.106, RSE 27.7%; eq 9
    e_cyp3a4_inh_strong_vc_ac886 <- 1.92
    label("Fractional change in AC886 Vc with strong CYP3A inhibitors (unitless)") # Table 4: 1.92, RSE 19.3%; eq 8

    # ------------------------------------------------------------------
    # INTERINDIVIDUAL VARIABILITY
    # Tables 3 and 4 print each IIV as '% CV'. The paper converts its
    # log-scale variances to %CV as 100 * sqrt(variance): Table 4 footnote a
    # gives variance 0.0690 -> 26.3% CV and 0.169 -> 41.2% CV. The same
    # convention is applied here, so omega^2 = (CV/100)^2.
    # ------------------------------------------------------------------
    etalcl ~ 0.3036 # Table 3: CL 55.1% CV, RSE 6.65% -> 0.551^2
    etalvc ~ 0.07618 # Table 3: Vc 27.6% CV, RSE 8.39% -> 0.276^2
    etalq ~ 0.06150 # Table 3: Q1 24.8% CV, RSE 23.3% -> 0.248^2
    etalvp ~ 0.1576 # Table 3: Vp1 39.7% CV, RSE 10.6% -> 0.397^2
    etalka ~ 0.1482 # Table 3: ka 38.5% CV, RSE 9.14% -> 0.385^2
    etald1 ~ 0.4802 # Table 3: D1 69.3% CV, RSE 8.69% -> 0.693^2
    etalfdepot ~ 0.1211 # Table 3: F1 34.8% CV, RSE 13.2% -> 0.348^2
    etaltlag ~ 0.3844 # Table 3: ALAG1 62.0% CV, RSE 7.36% -> 0.620^2
    etalvc_ac886 ~ 1.369 # Table 4: Vcm 117% CV, RSE 8.30% -> 1.17^2
    etalcl_ac886_aml ~ 0.4109 # Table 4: IIV in CLm in patients 64.1% CV, RSE 8.62% -> 0.641^2
    etalcl_ac886_nonaml ~ 0.2107 # Table 4: IIV in CLm in healthy volunteers 45.9% CV, RSE 8.50% -> 0.459^2

    # ------------------------------------------------------------------
    # INTEROCCASION VARIABILITY -- Table 3; one variance per parameter
    # shared across occasions (later occasions fixed to the first).
    # ------------------------------------------------------------------
    etaiov_fdepot_1 ~ 0.05108 # Table 3: IOV in F1 22.6% CV, RSE 20.6% -> 0.226^2 (occasion 1)
    etaiov_fdepot_3 ~ fix(0.05108) # Table 3: IOV in F1 on occasion 3 (same variance)
    etaiov_cl_1 ~ 0.1673 # Table 3: IOV in CL 40.9% CV, RSE 7.85% -> 0.409^2 (occasion 1)
    etaiov_cl_2 ~ fix(0.1673) # Table 3: IOV in CL on occasion 2 (same variance)
    etaiov_cl_3 ~ fix(0.1673) # Table 3: IOV in CL on occasion 3 (same variance)
    etaiov_cl_4 ~ fix(0.1673) # Table 3: IOV in CL on occasion 4 (same variance)
    etaiov_vc_1 ~ 0.04121 # Table 3: IOV in Vc 20.3% CV, RSE 17.0% -> 0.203^2 (occasion 1)
    etaiov_vc_2 ~ fix(0.04121) # Table 3: IOV in Vc on occasion 2 (same variance)
    etaiov_vc_3 ~ fix(0.04121) # Table 3: IOV in Vc on occasion 3 (same variance)
    etaiov_vc_4 ~ fix(0.04121) # Table 3: IOV in Vc on occasion 4 (same variance)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # Quizartinib: Var = F^2 * sigma_prop + sigma_add (Table 3 footnotes b
    # and c), i.e. nlmixr2 combined2 with SDs sqrt(sigma).
    # AC886: additive on log-transformed concentrations (Table 4).
    # ------------------------------------------------------------------
    propSd <- 0.07503
    label("Quizartinib proportional residual SD, healthy volunteers (fraction)") # Table 3: HV CCV variance 0.00563 -> sqrt = 0.0750
    propSd_aml <- 0.1939
    label("Quizartinib proportional residual SD, AML patients (fraction)") # Table 3: patient CCV variance 0.0376 -> sqrt = 0.1939
    addSd <- 0.9778
    label("Quizartinib additive residual SD, shared by both populations (ng/mL)") # Table 3: additive variance 0.956 (patient value SAME) -> sqrt = 0.9778
    expSd_ac886 <- 0.2627
    label("AC886 log-scale residual SD, healthy volunteers") # Table 4: variance 0.0690 = 0.263 SD
    expSd_ac886_aml <- 0.4111
    label("AC886 log-scale residual SD, AML patients") # Table 4: variance 0.169 = 0.412 SD
  })

  model({
    # Centring constants (Kang 2020 eqs 5-7 and 9: pooled medians, Table 2)
    ref_alb <- 4.1
    ref_bsa <- 1.9
    ref_age <- 44
    alb_gdl <- ALB * 0.1 # canonical g/L -> g/dL used by the published fit

    # Occasion indicators for the interoccasion variability
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 + oc4 * etaiov_vc_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc3 * etaiov_fdepot_3

    # ------------------------------------------------------------------
    # QUIZARTINIB individual parameters -- eqs 1-7
    # ------------------------------------------------------------------
    # eq 1: ka = 0.874 * (1 - PAT) * (1 + FED * -0.512) + 1.68 * PAT
    ka <- (exp(lka) * (1 - DIS_AML) * (1 + e_fed_ka * FED) + exp(lka_aml) * DIS_AML) * exp(etalka)

    # eq 4: both fractional changes sit inside one additive bracket
    cl <- exp(lcl + etalcl + iov_cl) *
      (1 + e_dis_aml_cl * DIS_AML + e_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG)
    vc <- exp(lvc + etalvc + iov_vc) * (alb_gdl / ref_alb)^e_alb_vc * (BSA / ref_bsa)^e_bsa_vc
    vp <- exp(lvp + etalvp) * (BSA / ref_bsa)^e_bsa_vp * (AGE / ref_age)^e_age_vp
    q <- exp(lq + etalq) * (BSA / ref_bsa)^e_bsa_q
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)

    # eqs 2 and 3: study AC220-019 has its own F1 in place of eq 2
    f_main <- exp(lfdepot) * ((1 - DIS_AML) + e_dis_aml_fdepot * DIS_AML) *
      (1 + e_cyp3a4_inh_strong_fdepot * CONMED_CYP3A4_INH_STRONG)
    f_ac220019 <- exp(lfdepot_ac220019) * (1 + e_fed_fdepot * FED)
    fdepot <- (f_main * (1 - STUDY_AC220019) + f_ac220019 * STUDY_AC220019) *
      exp(etalfdepot + iov_fdepot)

    # ------------------------------------------------------------------
    # AC886 individual parameters -- eqs 8 and 9
    # ------------------------------------------------------------------
    iiv_cl_ac886 <- DIS_AML * etalcl_ac886_aml + (1 - DIS_AML) * etalcl_ac886_nonaml
    cl_ac886 <- exp(lcl_ac886 + iiv_cl_ac886) * (BSA / ref_bsa)^e_bsa_cl_ac886 *
      (1 + e_race_black_cl_ac886 * RACE_BLACK + e_cyp3a4_inh_strong_cl_ac886 * CONMED_CYP3A4_INH_STRONG)
    vc_ac886 <- exp(lvc_ac886 + etalvc_ac886) * (1 + e_cyp3a4_inh_strong_vc_ac886 * CONMED_CYP3A4_INH_STRONG)
    vp_ac886 <- exp(lvp_ac886)
    q_ac886 <- exp(lq_ac886)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_ac886 <- cl_ac886 / vc_ac886
    k45 <- q_ac886 / vc_ac886
    k54 <- q_ac886 / vp_ac886

    # Figure 2: fMET * CL of the parent forms AC886; (1 - fMET) * CL leaves
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_ac886) <- fmet * kel * central - kel_ac886 * central_ac886 -
      k45 * central_ac886 + k54 * peripheral1_ac886
    d/dt(peripheral1_ac886) <- k45 * central_ac886 - k54 * peripheral1_ac886

    # Sequential zero-order input over D1 into the depot after the lag,
    # then first-order ka to central (Table 3 'duration of zero-order input
    # to depot compartment'; Figure 2). Dose records need rate = -2 for
    # the modelled dur() to apply.
    alag(depot) <- tlag
    dur(depot) <- d1
    f(depot) <- fdepot

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_ac886 <- 1000 * central_ac886 / vc_ac886

    propSd_i <- propSd * (1 - DIS_AML) + propSd_aml * DIS_AML
    expSd_ac886_i <- expSd_ac886 * (1 - DIS_AML) + expSd_ac886_aml * DIS_AML
    Cc ~ add(addSd) + prop(propSd_i) + combined2()
    Cc_ac886 ~ lnorm(expSd_ac886_i)
  })
}
