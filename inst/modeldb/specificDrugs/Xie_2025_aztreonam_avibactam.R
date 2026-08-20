Xie_2025_aztreonam_avibactam <- function() {
  description <- paste(
    "Simultaneous four-compartment (two compartments per drug) population PK",
    "model for the fixed-ratio (3:1) aztreonam-avibactam combination, fitted",
    "jointly to 4,914 aztreonam plasma samples from 431 subjects and 18,222",
    "avibactam plasma samples from 2,635 subjects pooled across Phase 1,",
    "Phase 2a and Phase 3 aztreonam-avibactam studies plus the",
    "ceftazidime-avibactam development program (Xie 2025). Both drugs are",
    "given as zero-order intravenous infusions with first-order elimination.",
    "Body weight is an allometric covariate on CL, Vc, Q and Vp of both drugs",
    "with fixed exponents. Time-varying BSA-normalized creatinine clearance",
    "acts on the clearance of both drugs through a hinged relationship: a",
    "power function below 80 mL/min/1.73 m^2 and a linear function at or",
    "above it. Infection type shifts clearance and central volume, and",
    "avibactam additionally carries end-stage-renal-disease, hemodialysis,",
    "Chinese-race and APACHE II severity terms. The two drugs are linked by a",
    "four-way OMEGA block across aztreonam and avibactam CL and Vc, whose",
    "cross-drug correlations are 0.98 and 0.99. Aztreonam is the unsuffixed",
    "parent; avibactam carries the sibling-drug suffix _avi throughout."
  )
  reference <- paste(
    "Xie R, Rogers H, Chow JW, Soto E, Raber SR. Population",
    "pharmacokinetic/pharmacodynamic modeling to optimize aztreonam-avibactam",
    "dose regimens for adult patients. Antimicrob Agents Chemother.",
    "2025;69(8):e01950-24. doi:10.1128/aac.01950-24.",
    "All fixed-effect, random-effect and residual-error estimates are the",
    "'Final model (run 46)' column of supplementary Table S3. The functional",
    "FORMS of the covariate relationships are not printed in Xie 2025; they",
    "are taken from the immediate structural predecessor by the same group,",
    "Das S, Riccobene T, Carrothers TJ, Wright JG, MacPherson M, Cristinacce",
    "A, McFadyen L, Xie R, Luckey A, Raber S. Dose selection for",
    "aztreonam-avibactam, including adjustments for renal impairment, for",
    "Phase IIa and Phase III evaluation. Eur J Clin Pharmacol.",
    "2024;80(4):529-543. doi:10.1007/s00228-023-03609-x, supplementary",
    "Tables 4 and 5.",
    sep = " "
  )
  vignette <- "Xie_2025_aztreonam_avibactam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Aztreonam residual variability is stratified by study phase, so the
  # canonical propSd / addSd consumed by the error model are derived inside
  # model() from four phase-specific ini() magnitudes. Same construction as
  # Cammarata_2024_sulbactam_durlobactam and Valenzuela_2025_nipocalimab.
  paper_specific_residual_sds <- c(
    "propSdPhase1", "propSdPhase2", "propSdPhase3", "propSdPhase23",
    "addSdPhase1", "addSdPhase1_avi"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Xie 2025 (doses in mg; plasma
  # concentrations reported in mg/L).
  compartmentData <- list(
    central         = list(analyte = "aztreonam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "aztreonam", units = "mg", specimen = "plasma", verified = TRUE),
    central_avi     = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_avi = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Xie 2025 Results, 'Final population PK model': 'Body weight",
        "(normalized to 70 kg) was a structural covariate on both aztreonam",
        "and avibactam CL, Vc, Vp, and Q using allometric scaling with fixed",
        "exponents of 0.75 on CL and Q, and 1 on Vc and Vp.' The exponents",
        "are therefore shared by the two drugs and by the two members of each",
        "pair, which is why a single e_wt_cl_q and a single e_wt_vc_vp enter",
        "the whole model. Reference weight 70 kg is stated in the same",
        "sentence. Pooled cohort medians 74 kg (aztreonam data set, range",
        "33-130; Table S1) and 70 kg (avibactam data set, range 28-190;",
        "Table S2). Baseline (time-fixed) per subject."
      ),
      source_name        = "body weight"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance normalized to body surface area (nCrCL);",
        "time-varying"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Xie 2025 Results: 'Time-varying nCrCL effects on aztreonam and",
        "avibactam CL were described by both linear (nCrCL >= 80",
        "mL/min/1.73 m^2) and power (nCrCL <80 mL/min/1.73 m^2)",
        "relationships for subjects not on dialysis.' The hinge value 80 is",
        "both the breakpoint and the reference: the power arm is",
        "(nCrCL/80)^exponent and the linear arm is (1 + slope * (nCrCL - 80)),",
        "so the two arms agree exactly at nCrCL = 80 and the typical CL",
        "values in ini() are the values AT nCrCL = 80. Those forms are not",
        "printed in Xie 2025; they are the forms printed verbatim in the",
        "predecessor Das 2024 supplementary Table 5 rows 'theta7: CL,",
        "(CrCL/80)**theta7, CrCL <80 mL/min' and 'theta8: CL, (1 + theta8 *",
        "(CrCL-80)), CrCL >=80 mL/min'.",
        "TIME-VARYING within subject: Xie 2025 calls this out explicitly and",
        "notes that for the simulations 'CrCL values for each subject were",
        "kept constant throughout the time course of the simulations'.",
        "This is nCrCL (BSA-normalized). The renal-function GROUP boundaries",
        "quoted throughout the paper (ARC >150; normal >80 to <=150; mild >50",
        "to <=80; moderate >30 to <=50; severe >15 to <=30; ESRD <15) are",
        "stated in raw CrCL (mL/min), NOT BSA-normalized - do not use them",
        "interchangeably with this column.",
        "Overrides: subjects flagged RENALIMP_ESRD or RRT_HEMODIAL_STATUS",
        "bypass this covariate on avibactam CL entirely (see those entries);",
        "aztreonam CL always uses it."
      ),
      source_name        = "nCrCL"
    ),
    DIS_CIAI = list(
      description        = "Complicated intra-abdominal infection cohort indicator (1 = cIAI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Phase 1 subject; all infection-type indicators 0)",
      notes              = paste(
        "Xie 2025 Table S3 rows 'cIAI on CL_AVI' (theta27 = 0.115), 'cIAI on",
        "CL_ATM' (theta31 = 0.279) and the shared 'cIAI/NP on Vc'",
        "(theta25 = 0.931). cIAI raises the clearance of BOTH drugs, which is",
        "why Xie 2025 selected cIAI as the most conservative population for",
        "the dose-selection simulations: 'exposures of both aztreonam and",
        "avibactam were lowest in cIAI relative to other infection types'.",
        "The four infection-type indicators (DIS_CIAI, DIS_HABP, DIS_VABP,",
        "DIS_CUTI) are mutually exclusive and share the all-zero reference.",
        "226 of 431 subjects in the aztreonam data set and 1,012 of 2,635 in",
        "the avibactam data set had cIAI (Tables S1 and S2)."
      ),
      source_name        = "cIAI"
    ),
    DIS_HABP = list(
      description        = "Hospital-acquired bacterial pneumonia cohort indicator (1 = HAP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Phase 1 subject; all infection-type indicators 0)",
      notes              = paste(
        "Xie 2025 Results: 'same cIAI and NP (HAP/VAP) effect on aztreonam",
        "and avibactam Vc', i.e. the single Table S3 coefficient",
        "'cIAI/NP on Vc' (theta25 = 0.931) is shared by cIAI, HAP and VAP and",
        "by both drugs. Following the DIS_VABP register convention, the HAP",
        "and VAP COLUMNS are kept distinct and the shared coefficient is",
        "applied to their sum inside model(). Nosocomial pneumonia carries no",
        "separate CL term for either drug. 72 subjects in the aztreonam data",
        "set and 482 in the avibactam data set had HAP or VAP (Tables S1",
        "and S2, which pool the two as 'HAP/VAP')."
      ),
      source_name        = "NP (HAP/VAP)"
    ),
    DIS_VABP = list(
      description        = "Ventilator-associated bacterial pneumonia cohort indicator (1 = VAP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Phase 1 subject; all infection-type indicators 0)",
      notes              = paste(
        "Shares the 'cIAI/NP on Vc' coefficient (theta25 = 0.931) with",
        "DIS_CIAI and DIS_HABP on the central volume of both drugs; see the",
        "DIS_HABP notes. Xie 2025 reports that steady-state exposures were",
        "similar between VAP and cIAI patients while HAP patients had the",
        "highest exposures, but the final model carries no VAP-specific",
        "coefficient distinct from HAP."
      ),
      source_name        = "NP (HAP/VAP)"
    ),
    DIS_CUTI = list(
      description        = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Phase 1 subject; all infection-type indicators 0)",
      notes              = paste(
        "Xie 2025 Table S3 rows 'cUTI on Vc_AVI' (theta24 = 1.5) and",
        "'cUTI on CL_AVI' (theta26 = 0.222). Both terms are AVIBACTAM-only:",
        "aztreonam carries no cUTI coefficient, because only 3 of the 431",
        "subjects in the aztreonam data set had cUTI, against 707 of 2,635 in",
        "the avibactam data set (Tables S1 and S2). Do not borrow the",
        "avibactam coefficients for aztreonam."
      ),
      source_name        = "cUTI"
    ),
    RENALIMP_ESRD = list(
      description        = paste(
        "End-stage renal disease indicator, not receiving hemodialysis",
        "(1 = ESRD)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any renal function above ESRD)",
      notes              = paste(
        "Xie 2025 Table S3 row 'ESRD on CL_AVI' (theta11 = -0.923), applied",
        "as a proportional shift so ESRD subjects retain 1 - 0.923 = 7.7% of",
        "the reference avibactam clearance. AVIBACTAM ONLY - aztreonam",
        "carries no ESRD term.",
        "The ESRD term REPLACES the nCrCL relationship rather than",
        "multiplying it: Das 2024 supplementary Table 4 prints the power arm",
        "as 'theta8: CL, (nCrCL/80)**theta8, nCrCL <80 mL/min AND NOT ESRD',",
        "and the resulting 7.7% is consistent with the 6.71% direct ESRD",
        "multiplier of Das 2024 supplementary Table 5. Multiplying the two",
        "instead would give a physiologically impossible ~0.7% of reference",
        "clearance.",
        "Xie 2025 defines ESRD as CrCL <15 mL/min (raw Cockcroft-Gault, NOT",
        "BSA-normalized), so this column is NOT derivable from the CRCL",
        "column by a simple threshold; supply it from the raw CrCL. 12 of the",
        "431 aztreonam-data-set subjects were ESRD, with median raw CrCL",
        "6.6 mL/min (range 4.6-14.8; Table S1).",
        "Mutually exclusive with RRT_HEMODIAL_STATUS: a subject on",
        "hemodialysis takes the separate dialysis clearance instead."
      ),
      source_name        = "ESRD"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = paste(
        "Intermittent-hemodialysis treatment-status indicator",
        "(1 = subject is a dialysis patient)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not receiving hemodialysis)",
      notes              = paste(
        "Xie 2025 Table S3 row 'CL_AVI_DIAL' (theta12 = 17.9 L/h) and",
        "Results: 'a separate avibactam CL for dialysis patients'. This is an",
        "ABSOLUTE clearance that REPLACES the whole reference-plus-nCrCL",
        "construction for avibactam, not a multiplier on it - hence it is",
        "carried as its own lcl_dial_avi ini() parameter rather than as an",
        "e_* covariate effect. AVIBACTAM ONLY.",
        "Subject-level (time-fixed) rather than per-session: neither Xie 2025",
        "nor Das 2024 resolves clearance during versus between dialysis",
        "sessions, so RRT_HEMODIAL_STATUS (the subject-level canonical) is",
        "used rather than the per-session gate RRT_HEMODIAL_ACTIVE. A",
        "consequence is that this model cannot reproduce the intradialytic /",
        "interdialytic sawtooth; see the vignette Errata.",
        "Mutually exclusive with RENALIMP_ESRD in this encoding.",
        "Data source: the dialysis subjects entered through the",
        "ceftazidime-avibactam program, not the aztreonam-avibactam Phase 3",
        "trials, which 'did not enroll patients with end-stage renal disease",
        "(ESRD) requiring dialysis'."
      ),
      source_name        = "DIAL"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage race indicator (1 = Chinese)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese)",
      notes              = paste(
        "Xie 2025 Table S3 row 'China on Vc_AVI' (theta29 = -0.145), a 14.5%",
        "reduction in avibactam central volume. AVIBACTAM ONLY; aztreonam",
        "Vc carries no race term. Introduced by the addition of the Phase 1",
        "healthy-Chinese-subject study C3601007 (NCT04973826). 84 of 431",
        "subjects in the aztreonam data set were recorded as Chinese",
        "(Table S1).",
        "Note that Xie 2025 separately reports, from the post hoc exposure",
        "analysis, 'no obvious differences in aztreonam and avibactam",
        "exposures correlating with ... China and non-China (rest of world)'",
        "- that statement is about total exposure (AUC, driven by CL) and is",
        "not in conflict with a volume-only covariate."
      ),
      source_name        = "China"
    ),
    APACHE_II_SEV = list(
      description        = paste(
        "Elevated-APACHE-II severity stratum indicator",
        "(1 = subject falls in the elevated-severity stratum)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in the elevated-APACHE-II severity stratum)",
      notes              = paste(
        "Xie 2025 Table S3 row 'APACHE II score on CL_AVI CL'",
        "(theta30 = -0.118): an 11.8% reduction in avibactam clearance.",
        "AVIBACTAM ONLY.",
        "ENCODED AS BINARY, and this is an inference rather than a printed",
        "statement. Neither Xie 2025 nor Das 2024 writes the equation, but",
        "Das 2024 supplementary Table 5 prints the row as 'theta16: APACHE",
        "effect on CL, CL*(1 + theta16)' - the same grammar Das uses for",
        "every BINARY population effect in that table ('Population effect on",
        "Vc (cUTI), Vc*(1 + theta11)'), and pointedly NOT the grammar Das",
        "uses for continuous covariates, which always carry the covariate",
        "symbol in the expression ('AGE/35**theta14 on CL',",
        "'(CrCL/80)**theta7'). A raw continuous score cannot be intended: at",
        "a typical ICU APACHE II of 10, CL*(1 + (-0.118)*10) is negative.",
        "THE THRESHOLD THAT SETS THIS FLAG IS NOT STATED IN ANY SOURCE ON",
        "DISK. It is therefore left to the user's data rather than invented",
        "here, and the validation vignette holds it at the 0 reference for",
        "every simulation so that no reproduced published value depends on",
        "it. Supply it yourself if you wish to exercise the term. See the",
        "vignette Errata.",
        "Distinct from the continuous canonical APACHE_II (score in points):",
        "this column is the derived severity stratum, not the score."
      ),
      source_name        = "APACHE II score"
    ),
    STUDY_CIAI_PH2 = list(
      description        = paste(
        "Phase 2 cIAI study cohort indicator from the ceftazidime-avibactam",
        "program (1 = Study 2002)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other study)",
      notes              = paste(
        "Xie 2025 Table S3 rows 'Study2002 on CL_AVI' (theta13 = 0.89) and",
        "'Study2002 on Vc_AVI' (theta14 = 1.64), retained from prior",
        "modelling: Xie 2025 Results lists 'a separate avibactam CL and Vc",
        "for phase 2 patients with cIAI identified in prior modeling' among",
        "the effects carried into the base model. Das 2024 supplementary",
        "Table 5 labels the same pair 'Population effect on CL (cIAI, phase",
        "II), CL*(1 + theta10)' and 'Population effect on Vc (cIAI, phase",
        "II), Vc*(1 + theta9)', which is where the proportional form comes",
        "from. AVIBACTAM ONLY.",
        "This is an ADDITIONAL study-cohort effect that stacks on top of the",
        "DIS_CIAI infection-type effect for those subjects; it is not a",
        "substitute for it. Zero for every aztreonam-avibactam subject and",
        "for the whole simulated population of this paper, which is why it",
        "does not enter the vignette's reproductions."
      ),
      source_name        = "Study2002"
    ),
    STUDY_AZTAVI_PHASE2 = list(
      description        = "Phase 2 study-stratum indicator for the aztreonam residual-error model",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with the other two phase indicators also 0 (Phase 1 stratum)",
      notes              = paste(
        "Selects the aztreonam Phase 2 proportional residual magnitude,",
        "Xie 2025 Table S3 'Prop RSV_ATM phase 2' (theta21 = 22.4%).",
        "Xie 2025 Results lists 'separate residual error models for aztreonam",
        "and avibactam and study phases' among the retained effects. Table S3",
        "reports four aztreonam proportional magnitudes (unqualified,",
        "phase 2, phase 3, phase 2/3) and one Phase 1 additive term; the",
        "unqualified row is taken as the Phase 1 stratum, matching the",
        "predecessor Das 2024 supplementary Table 5, which reports",
        "'Proportional error, phase I / II / III' plus an additive term for",
        "phase I only. Avibactam residual variability is NOT phase-stratified",
        "(a single proportional term), so none of the three phase indicators",
        "enters the avibactam proportional error - only its Phase 1 additive",
        "term. Paired with STUDY_AZTAVI_PHASE3 and STUDY_AZTAVI_PHASE23.",
        "Residual error does not affect typical-value or individual-prediction",
        "simulation, only the simulated observation."
      ),
      source_name        = "study phase"
    ),
    STUDY_AZTAVI_PHASE3 = list(
      description        = "Phase 3 study-stratum indicator for the aztreonam residual-error model",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with the other two phase indicators also 0 (Phase 1 stratum)",
      notes              = paste(
        "Selects the aztreonam Phase 3 proportional residual magnitude,",
        "Xie 2025 Table S3 'Prop RSV_ATM phase 3' (theta22 = 40.3%). See the",
        "STUDY_AZTAVI_PHASE2 notes for the full residual-error rationale."
      ),
      source_name        = "study phase"
    ),
    STUDY_AZTAVI_PHASE23 = list(
      description        = "Pooled phase-2/3 study-stratum indicator for the aztreonam residual-error model",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with the other two phase indicators also 0 (Phase 1 stratum)",
      notes              = paste(
        "Selects the aztreonam pooled phase-2/3 proportional residual",
        "magnitude, Xie 2025 Table S3 'Prop RSV_ATM phase 2/3'",
        "(theta23 = 53.3%). Xie 2025 does not define which records make up",
        "this fourth stratum as distinct from the separate phase 2",
        "(theta21) and phase 3 (theta22) strata, and no supplement or",
        "predecessor on disk resolves it; the coefficient is carried here so",
        "that the published parameter set is complete, with the stratum",
        "membership left to the user. See the vignette Errata. See the",
        "STUDY_AZTAVI_PHASE2 notes for the full residual-error rationale."
      ),
      source_name        = "study phase"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 2635L,
    n_studies        = 5L,
    n_concentrations = 23136L,
    age_range        = "18-89 years",
    age_median       = "50 years (aztreonam data set); 52 years (avibactam data set)",
    weight_range     = "28-190 kg",
    weight_median    = "74 kg (aztreonam data set); 70 kg (avibactam data set)",
    sex_female_pct   = 38.6,
    race_ethnicity   = c(
      White = 55.5, Chinese = 19.5, Asian = 8.4, Black = 3.0,
      `Native American` = 1.9, Other = 1.4, Missing = 10.4
    ),
    disease_state    = paste(
      "Pooled healthy adults and hospitalized adults with complicated",
      "intra-abdominal infection, nosocomial pneumonia (hospital-acquired or",
      "ventilator-associated), complicated urinary tract infection, or",
      "bloodstream infection. Infection-type counts in the aztreonam data set",
      "(Table S1): healthy 126, cIAI 226, HAP/VAP 72, cUTI 3, BSI 4. In the",
      "avibactam data set (Table S2): healthy 430, cIAI 1,012, HAP/VAP 482,",
      "cUTI 707, BSI 4."
    ),
    dose_range       = paste(
      "All doses given as intravenous infusions of the fixed 3:1",
      "aztreonam:avibactam combination. Phase 3 regimen: 500/167 mg loading",
      "dose over 30 min immediately followed by a 1,500/500 mg extended",
      "loading dose over 3 h, then 1,500/500 mg maintenance doses over 3 h",
      "every 6 h, with reductions for renal impairment. Approved regimens",
      "(Table 5): simplified loading dose 2,000/667 mg over 3 h with",
      "1,500/500 mg q6h maintenance for normal renal function and augmented",
      "renal clearance; 2,000/667 mg then 750/250 mg q6h for moderate",
      "impairment; 1,350/450 mg then 675/225 mg q8h for severe impairment;",
      "1,000/334 mg then 675/225 mg q12h for ESRD."
    ),
    regions          = "Multinational; includes a dedicated Phase 1 study in healthy Chinese subjects (C3601007) and Chinese sites within the Phase 3 program",
    renal_function   = paste(
      "Spans augmented renal clearance to end-stage renal disease. Baseline",
      "CrCL group counts in the aztreonam data set (Table S1): ARC (>150",
      "mL/min) 46, normal (>80 to <=150) 244, mild (>50 to <=80) 74,",
      "moderate (>30 to <=50) 38, severe (>15 to <=30) 17, ESRD (<=15) 12.",
      "Group medians of raw CrCL: 176.7, 117.2, 62.9, 38.6, 20.0 and 6.6",
      "mL/min respectively. The Phase 3 trials themselves enrolled no ESRD",
      "or continuous-renal-replacement-therapy patients; those data came",
      "from the ceftazidime-avibactam program and a dedicated Phase 1",
      "severe-renal-impairment study (C3601006, NCT04486625)."
    ),
    bmi_range        = "10.4-62.5 kg/m^2 (medians 25.1 aztreonam / 24.7 avibactam data set)",
    unbound_fraction = paste(
      "Xie 2025 Methods: total plasma concentrations were converted to free",
      "concentrations for the PK/PD target calculations 'using unbound",
      "fractions of 0.616 for aztreonam and 0.92 for avibactam'. These are",
      "fixed conversion constants used downstream of the PK model, not",
      "fitted parameters, so they are not carried in ini(); the validation",
      "vignette applies them."
    ),
    pkpd_target      = paste(
      "Joint target, both achieved simultaneously: free aztreonam above an",
      "aztreonam-avibactam MIC of 8 mg/L for 60% of the dosing interval",
      "(60% fT>MIC) and free avibactam above a threshold concentration of",
      "2.5 mg/L for 50% of the dosing interval (50% fT>CT)."
    ),
    notes            = paste(
      "n_subjects and n_concentrations are the union across analytes: the",
      "simultaneous fit used 4,914 aztreonam concentrations from 431",
      "subjects and 18,222 avibactam concentrations from 2,635 subjects, and",
      "every aztreonam subject also contributed avibactam data. n_studies",
      "counts the five clinical studies newly ADDED in this analysis (two",
      "Phase 1, one Phase 2a REJUVENATE, two Phase 3 REVISIT and ASSEMBLE;",
      "Table 1); these were pooled with the data of the previous",
      "non-simultaneous aztreonam-avibactam analysis and with avibactam data",
      "from the ceftazidime-avibactam program, so the true study count",
      "underlying the fit is larger and is not enumerated in Xie 2025.",
      "Registered on ClinicalTrials.gov as NCT03329092 and NCT03580044.",
      "Estimation was FOCE with interaction in NONMEM 7.4.3."
    )
  )

  ini({
    # =====================================================================
    # AZTREONAM structural parameters (Xie 2025 Table S3, 'Final model
    # (run 46)' column). Values refer to the reference subject: WT = 70 kg,
    # nCrCL = 80 mL/min/1.73 m^2, no infection (healthy Phase 1 subject).
    # Aztreonam is the unsuffixed parent.
    # =====================================================================
    lcl <- log(5.00); label("Aztreonam clearance at nCrCL = 80 mL/min/1.73 m^2, WT = 70 kg (L/h)") # Table S3: CL_ATM (theta1) = 5.00 L/h (RSE 2.37%)
    lvc <- log(7.01); label("Aztreonam central volume of distribution at WT = 70 kg (L)")          # Table S3: Vc_ATM (theta2) = 7.01 L (RSE 3.34%)
    lq  <- log(9.41); label("Aztreonam inter-compartmental clearance at WT = 70 kg (L/h)")         # Table S3: Q_ATM (theta3) = 9.41 L/h (RSE 7.65%)
    lvp <- log(6.12); label("Aztreonam peripheral volume of distribution at WT = 70 kg (L)")       # Table S3: Vp_ATM (theta4) = 6.12 L (RSE 3.31%)

    # =====================================================================
    # AVIBACTAM structural parameters (Xie 2025 Table S3, same column and
    # same reference subject).
    # =====================================================================
    lcl_avi <- log(10.3); label("Avibactam clearance at nCrCL = 80 mL/min/1.73 m^2, WT = 70 kg (L/h)") # Table S3: CL_AVI (theta5) = 10.3 L/h (RSE 1.61%)
    lvc_avi <- log(12.6); label("Avibactam central volume of distribution at WT = 70 kg (L)")          # Table S3: Vc_AVI (theta6) = 12.6 L (RSE 1.8%)
    lq_avi  <- log(4.82); label("Avibactam inter-compartmental clearance at WT = 70 kg (L/h)")         # Table S3: Q_AVI (theta7) = 4.82 L/h (RSE 3.54%)
    lvp_avi <- log(6.95); label("Avibactam peripheral volume of distribution at WT = 70 kg (L)")       # Table S3: Vp_AVI (theta8) = 6.95 L (RSE 1.75%)

    lcl_dial_avi <- log(17.9); label("Avibactam clearance in hemodialysis patients, at WT = 70 kg (L/h)")
    # Table S3: CL_AVI_DIAL (theta12) = 17.9 L/h (RSE 11.2%). An ABSOLUTE
    # clearance that replaces the reference-plus-nCrCL construction for
    # dialysis patients, not a multiplier - hence an lcl_* parameter rather
    # than an e_* covariate effect. See covariateData$RRT_HEMODIAL_STATUS.

    # =====================================================================
    # Body-weight allometry (Xie 2025 Results). The exponents are FIXED and
    # SHARED across the two drugs and across each pair of parameters, so one
    # exponent serves CL and Q and one serves Vc and Vp, for both analytes.
    # =====================================================================
    e_wt_cl_q   <- fixed(0.75); label("Body-weight allometric exponent on CL and Q of both drugs, WT/70 (unitless)") # Results: 'allometric scaling with fixed exponents of 0.75 on CL and Q'
    e_wt_vc_vp  <- fixed(1);    label("Body-weight allometric exponent on Vc and Vp of both drugs, WT/70 (unitless)") # Results: '... and 1 on Vc and Vp'

    # =====================================================================
    # Renal-function effects on clearance. Hinged at nCrCL = 80
    # mL/min/1.73 m^2: a power arm below the hinge and a linear arm at or
    # above it, continuous at the hinge. Strata carry the _lt80 / _ge80
    # suffixes so neither arm silently claims the bare canonical name.
    # Functional forms from Das 2024 supplementary Table 5; values from
    # Xie 2025 Table S3.
    # =====================================================================
    e_crcl_cl_lt80     <- 0.502;   label("Aztreonam nCrCL power exponent on CL for nCrCL < 80 mL/min/1.73 m^2 (unitless)")        # Table S3: nCrCL on CL_ATM (theta9) = 0.502 (RSE 7.11%)
    e_crcl_cl_ge80     <- 0.00383; label("Aztreonam nCrCL linear slope on CL for nCrCL >= 80 (per mL/min/1.73 m^2)")              # Table S3: Slope nCrCL on CL_ATM (theta15) = 0.00383 (RSE 14.6%)
    e_crcl_cl_lt80_avi <- 1.06;    label("Avibactam nCrCL power exponent on CL for nCrCL < 80 mL/min/1.73 m^2 (unitless)")        # Table S3: nCrCL on CL_AVI (theta10) = 1.06 (RSE 4.61%)
    e_crcl_cl_ge80_avi <- 0.00313; label("Avibactam nCrCL linear slope on CL for nCrCL >= 80 (per mL/min/1.73 m^2)")              # Table S3: Slope nCrCL on CL_AVI (theta16) = 0.00313 (RSE 10.0%)

    e_renalimp_esrd_cl_avi <- -0.923; label("Proportional shift in avibactam CL for end-stage renal disease (fraction)")
    # Table S3: ESRD on CL_AVI (theta11) = -0.923 (RSE 1.95%), i.e. 7.7% of
    # reference CL. REPLACES the nCrCL arms rather than multiplying them -
    # see covariateData$RENALIMP_ESRD for the Das 2024 'and not ESRD'
    # evidence.

    # =====================================================================
    # Infection-type and population effects, all proportional shifts of the
    # form parameter * (1 + theta), per Das 2024 supplementary Table 5
    # ('Population effect on Vc (cUTI), Vc*(1 + theta11)').
    # =====================================================================
    e_ciai_cl     <- 0.279; label("Proportional shift in aztreonam CL for cIAI (fraction)")                     # Table S3: cIAI on CL_ATM (theta31) = 0.279 (RSE 12.4%)
    e_ciai_cl_avi <- 0.115; label("Proportional shift in avibactam CL for cIAI (fraction)")                     # Table S3: cIAI on CL_AVI (theta27) = 0.115 (RSE 17.7%)
    e_cuti_cl_avi <- 0.222; label("Proportional shift in avibactam CL for cUTI (fraction)")                     # Table S3: cUTI on CL_AVI (theta26) = 0.222 (RSE 14.7%)

    e_ciai_habp_vabp_vc <- 0.931; label("Proportional shift in the central volume of BOTH drugs for cIAI, HAP or VAP (fraction)")
    # Table S3: cIAI/NP on Vc (theta25) = 0.931 (RSE 6.37%). Results: 'same
    # cIAI and NP (HAP/VAP) effect on aztreonam and avibactam Vc' - a single
    # estimated coefficient shared by three infection types and two drugs.

    e_cuti_vc_avi          <- 1.5;    label("Proportional shift in avibactam Vc for cUTI (fraction)")            # Table S3: cUTI on Vc_AVI (theta24) = 1.5 (RSE 8.58%)
    e_race_chinese_vc_avi  <- -0.145; label("Proportional shift in avibactam Vc for Chinese race (fraction)")    # Table S3: China on Vc_AVI (theta29) = -0.145 (RSE 15.0%)
    e_apache_ii_sev_cl_avi <- -0.118; label("Proportional shift in avibactam CL for the elevated-APACHE-II stratum (fraction)") # Table S3: APACHE II score on CL_AVI (theta30) = -0.118 (RSE 12.3%)

    e_ciai_ph2_cl_avi <- 0.89; label("Proportional shift in avibactam CL for the Phase 2 cIAI study cohort (fraction)") # Table S3: Study2002 on CL_AVI (theta13) = 0.89 (RSE 16.7%)
    e_ciai_ph2_vc_avi <- 1.64; label("Proportional shift in avibactam Vc for the Phase 2 cIAI study cohort (fraction)") # Table S3: Study2002 on Vc_AVI (theta14) = 1.64 (RSE 22.5%)

    # =====================================================================
    # Inter-individual variability (Xie 2025 Table S3, final-model column).
    #
    # SCALE OF THE %CV COLUMN. Table S3 reports the diagonal elements as
    # '[CV%]' and the off-diagonals as bare covariances, so the two are only
    # mutually consistent under one reading of the %CV column. The paper
    # supplies its own answer key: the Discussion states that the
    # correlations between the aztreonam and avibactam variances 'were high
    # at 0.976 and 0.986' for CL and Vc. Taking omega^2 = (CV%/100)^2
    # reproduces those exactly - 0.171/(0.402*0.435) = 0.978 for CL and
    # 0.391/(0.608*0.652) = 0.986 for Vc - whereas the log-normal reading
    # omega^2 = log(CV^2 + 1) does not. The variances below are therefore
    # (CV%/100)^2 and the covariances are the printed values verbatim.
    #
    # The 4x4 block spans aztreonam CL and Vc and avibactam CL and Vc:
    # Results, 'an OMEGA block for variances of aztreonam CL and Vc and
    # avibactam CL and Vc'. All six off-diagonals are reported, so the block
    # is complete. Aztreonam carries NO IIV on Q or Vp - Discussion: 'it was
    # difficult to estimate IIVs for aztreonam on Q and Vp, possibly because
    # the majority of aztreonam data were collected at steady state'.
    # =====================================================================
    etalcl + etalvc + etalcl_avi + etalvc_avi ~
      c(0.161604,
        0.151,    0.369664,
        0.171,    0.173,    0.189225,
        0.182,    0.391,    0.212,    0.425104)
    # Diagonals: IIV CL_ATM 40.2 %CV -> 0.402^2 = 0.161604 (shrinkage 15.52%);
    #            IIV Vc_ATM 60.8 %CV -> 0.608^2 = 0.369664 (shrinkage 28.1%);
    #            IIV CL_AVI 43.5 %CV -> 0.435^2 = 0.189225 (shrinkage 14.47%);
    #            IIV Vc_AVI 65.2 %CV -> 0.652^2 = 0.425104 (shrinkage 25.68%).
    # Off-diagonals verbatim from Table S3: cov CL_ATM-Vc_ATM = 0.151;
    #            cov CL_ATM-CL_AVI = 0.171; cov CL_AVI-Vc_ATM = 0.173;
    #            cov CL_ATM-Vc_AVI = 0.182; cov Vc_ATM-Vc_AVI = 0.391;
    #            cov CL_AVI-Vc_AVI = 0.212.

    etalq_avi + etalvp_avi ~ c(0.098596,
                               0.0309,   0.031329)
    # Results: 'a covariance term between the avibactam IIVs in Q and Vp'.
    # IIV Q_AVI 31.4 %CV -> 0.314^2 = 0.098596 (shrinkage 72.12%);
    # IIV Vp_AVI 17.7 %CV -> 0.177^2 = 0.031329 (shrinkage 71.56%);
    # cov Q_AVI-Vp_AVI = 0.0309 (r = 0.56). Both shrinkages exceed 70%.

    # =====================================================================
    # Residual variability (Xie 2025 Table S3, final-model column).
    # Results: 'separate residual error models for aztreonam and avibactam
    # and study phases'. Aztreonam has four proportional magnitudes selected
    # by study stratum plus a Phase 1 additive term; avibactam has a single
    # proportional magnitude plus a Phase 1 additive term. The table reports
    # the proportional terms as percentages, so the fractions below are the
    # printed percentages divided by 100.
    # =====================================================================
    propSdPhase1  <- 0.125; label("Aztreonam proportional residual SD, Phase 1 stratum (fraction)")   # Table S3: Prop RSV_ATM (theta18) = 12.5% (RSE 7.18%)
    propSdPhase2  <- 0.224; label("Aztreonam proportional residual SD, Phase 2 stratum (fraction)")   # Table S3: Prop RSV_ATM phase 2 (theta21) = 22.4% (RSE 16.3%)
    propSdPhase3  <- 0.403; label("Aztreonam proportional residual SD, Phase 3 stratum (fraction)")   # Table S3: Prop RSV_ATM phase 3 (theta22) = 40.3% (RSE 7.48%)
    propSdPhase23 <- 0.533; label("Aztreonam proportional residual SD, pooled phase-2/3 stratum (fraction)") # Table S3: Prop RSV_ATM phase 2/3 (theta23) = 53.3% (RSE 3.05%)
    addSdPhase1   <- 0.197; label("Aztreonam additive residual SD, Phase 1 stratum (mg/L)")           # Table S3: Additive RSV_ATM phase 1 (theta17) = 0.197 mg/L (RSE 31.8%)

    propSd_avi      <- 0.20;    label("Avibactam proportional residual SD (fraction)")                # Table S3: Prop RSV_AVI (theta20) = 20% (RSE 6.21%)
    addSdPhase1_avi <- 0.00621; label("Avibactam additive residual SD, Phase 1 stratum (mg/L)")       # Table S3: Additive RSV_AVI phase 1 (theta19) = 0.00621 mg/L (RSE 11.1%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Renal-function factor on clearance.
    #
    #   nCrCL <  80:  (nCrCL / 80)^power
    #   nCrCL >= 80:  1 + slope * (nCrCL - 80)
    #
    # Both arms equal 1 at the hinge, so the ini() clearances are the
    # typical values at nCrCL = 80 mL/min/1.73 m^2. Written with 0/1
    # indicator arithmetic rather than a branch so both arms stay finite
    # for every positive nCrCL.
    # ------------------------------------------------------------------
    renal_cl <- (CRCL / 80)^e_crcl_cl_lt80 * (CRCL < 80) +
                (1 + e_crcl_cl_ge80 * (CRCL - 80)) * (CRCL >= 80)

    renal_cl_avi <- (CRCL / 80)^e_crcl_cl_lt80_avi * (CRCL < 80) +
                    (1 + e_crcl_cl_ge80_avi * (CRCL - 80)) * (CRCL >= 80)

    # Avibactam renal handling has three mutually exclusive regimes. A
    # dialysis patient takes the separate absolute clearance; an ESRD
    # patient not on dialysis takes the proportional ESRD shift; everyone
    # else takes the hinged nCrCL relationship above.
    cl_avi_renalarm <-
      exp(lcl_avi) * renal_cl_avi * (1 - RENALIMP_ESRD) * (1 - RRT_HEMODIAL_STATUS) +
      exp(lcl_avi) * (1 + e_renalimp_esrd_cl_avi) * RENALIMP_ESRD * (1 - RRT_HEMODIAL_STATUS) +
      exp(lcl_dial_avi) * RRT_HEMODIAL_STATUS

    # ------------------------------------------------------------------
    # 2. Infection-type factors. The four infection-type indicators are
    #    mutually exclusive, so their shifts are summed inside a single
    #    (1 + ...) bracket; the shared reference is the uninfected Phase 1
    #    subject. Non-exclusive covariates (race, study cohort, APACHE II
    #    stratum) get their own multiplicative brackets.
    # ------------------------------------------------------------------
    infect_vc <- 1 + e_ciai_habp_vabp_vc * (DIS_CIAI + DIS_HABP + DIS_VABP)

    # ------------------------------------------------------------------
    # 3. Aztreonam individual PK parameters.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * renal_cl *
      (WT / 70)^e_wt_cl_q *
      (1 + e_ciai_cl * DIS_CIAI)

    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc_vp *
      infect_vc

    q  <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    # ------------------------------------------------------------------
    # 4. Avibactam individual PK parameters.
    # ------------------------------------------------------------------
    cl_avi <- cl_avi_renalarm * exp(etalcl_avi) *
      (WT / 70)^e_wt_cl_q *
      (1 + e_ciai_cl_avi * DIS_CIAI + e_cuti_cl_avi * DIS_CUTI) *
      (1 + e_apache_ii_sev_cl_avi * APACHE_II_SEV) *
      (1 + e_ciai_ph2_cl_avi * STUDY_CIAI_PH2)

    vc_avi <- exp(lvc_avi + etalvc_avi) *
      (WT / 70)^e_wt_vc_vp *
      (infect_vc + e_cuti_vc_avi * DIS_CUTI) *
      (1 + e_race_chinese_vc_avi * RACE_CHINESE) *
      (1 + e_ciai_ph2_vc_avi * STUDY_CIAI_PH2)

    q_avi  <- exp(lq_avi  + etalq_avi)  * (WT / 70)^e_wt_cl_q
    vp_avi <- exp(lvp_avi + etalvp_avi) * (WT / 70)^e_wt_vc_vp

    # ------------------------------------------------------------------
    # 5. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_avi <- cl_avi / vc_avi
    k12_avi <- q_avi  / vc_avi
    k21_avi <- q_avi  / vp_avi

    # ------------------------------------------------------------------
    # 6. Four-compartment system: two independent two-compartment IV
    #    dispositions. The two drugs are co-administered in one infusion
    #    bag at a fixed 3:1 aztreonam:avibactam ratio but do not
    #    interconvert, so each takes its dose into its own central
    #    compartment as a zero-order infusion.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    d/dt(central_avi)     <- -(kel_avi + k12_avi) * central_avi + k21_avi * peripheral1_avi
    d/dt(peripheral1_avi) <-   k12_avi * central_avi            - k21_avi * peripheral1_avi

    # ------------------------------------------------------------------
    # 7. Observations. Dose in mg, volumes in L -> concentrations in mg/L.
    #    These are TOTAL plasma concentrations; multiply by the unbound
    #    fractions 0.616 (aztreonam) and 0.92 (avibactam) to obtain the
    #    free concentrations that the paper's PK/PD targets are defined on
    #    (see population$unbound_fraction).
    #
    #    Aztreonam residual variability is switched by study stratum; all
    #    three phase indicators 0 selects the Phase 1 stratum, which is the
    #    only stratum carrying an additive component for either drug.
    # ------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_avi <- central_avi / vc_avi

    phase1 <- 1 - STUDY_AZTAVI_PHASE2 - STUDY_AZTAVI_PHASE3 - STUDY_AZTAVI_PHASE23

    propSd <- propSdPhase1  * phase1 +
              propSdPhase2  * STUDY_AZTAVI_PHASE2 +
              propSdPhase3  * STUDY_AZTAVI_PHASE3 +
              propSdPhase23 * STUDY_AZTAVI_PHASE23
    addSd  <- addSdPhase1 * phase1

    addSd_avi <- addSdPhase1_avi * phase1

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_avi ~ add(addSd_avi) + prop(propSd_avi)
  })
}
