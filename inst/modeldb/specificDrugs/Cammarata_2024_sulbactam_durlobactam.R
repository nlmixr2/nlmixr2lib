Cammarata_2024_sulbactam_durlobactam <- function() {
  description <- paste(
    "Combined four-compartment (two compartments per drug) population PK model",
    "for the sulbactam-durlobactam combination, fitted simultaneously to 373",
    "subjects and 5,188 plasma concentrations pooled from six Phase 1 studies,",
    "one Phase 2 cUTI study, and one Phase 3 study in Acinetobacter",
    "baumannii-calcoaceticus complex infections (Cammarata 2024). Total",
    "clearance of each drug is the sum of a renal arm scaled by baseline",
    "BSA-normalized creatinine clearance and a non-renal arm, split by a fixed",
    "fraction excreted renally; total CL carries a further proportional",
    "downward shift in subjects with CLcr < 30 mL/min/1.73 m^2. Body weight",
    "acts as a power function on CL and central volume of both drugs;",
    "infection type and East Asian region act as proportional shifts. The",
    "paper's two sub-models are integrated into the same file: a time-varying",
    "hemodialysis gate multiplies total CL 6.24-fold (durlobactam) and",
    "8.19-fold (sulbactam) while a session is running, and epithelial lining",
    "fluid concentrations are returned as instantaneous ratios of the plasma",
    "concentrations (37.2% and 53.3% of total drug). Durlobactam uses the",
    "unsuffixed canonical compartment / parameter set; sulbactam carries the",
    "sibling-drug suffix _sbt throughout."
  )
  reference <- paste(
    "Cammarata AP, Safir MC, Trang M, Larson KB, O'Donnell JP, Bhavnani SM,",
    "Rubino CM. Population pharmacokinetic analyses for sulbactam-durlobactam",
    "using Phase 1, 2, and 3 data. Antimicrob Agents Chemother.",
    "2024;69(1):e00485-24. doi:10.1128/aac.00485-24.",
    "Plasma model parameters from Table 1, hemodialysis sub-model parameters",
    "from Table S4, and epithelial lining fluid sub-model parameters from",
    "Table S6 of the supplement.",
    sep = " "
  )
  vignette <- "Cammarata_2024_sulbactam_durlobactam"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  # Durlobactam plasma residual variability is stratified by study phase
  # (Phase 1 / 2 / 3), so the canonical propSd / addSd used by the error model
  # are derived inside model() from three phase-specific ini() magnitudes. Same
  # construction as Valenzuela_2025_nipocalimab and vanIersel_2018_posaconazole.
  paper_specific_residual_sds <- c(
    "propSdPhase1", "propSdPhase2", "propSdPhase3", "addSdPhase1"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power function on total CL and on central volume Vc of BOTH drugs",
        "(Cammarata 2024 Results, 'Population pharmacokinetic model",
        "development': 'All parameter-covariate relationships were described",
        "with a proportional shift with the exception of CLcr and weight,",
        "which were described using power functions'). Exponents are the",
        "Table 1 rows CLWTKG1 / V1WTKG1 (durlobactam) and CLWTKG1 / V3WTKG1",
        "(sulbactam). Reference weight is 75 kg, the reference-population",
        "body weight stated in Cammarata 2024 Results ('subjects with HABP",
        "infection who are not from an East Asian region, have a body weight",
        "of 75 kg, and have normal renal function') and repeated in the",
        "Discussion ('relative to those with normal body weight (75 kg)').",
        "Pooled cohort median 75 kg (range 35.8-150 kg; Table S3).",
        "Baseline (time-fixed) per subject."
      ),
      source_name        = "WTKG"
    ),
    CRCL = list(
      description        = paste(
        "Baseline Cockcroft-Gault creatinine clearance, normalized to body",
        "surface area"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cammarata 2024 Results, 'Data': CLcr 'was calculated using the",
        "Cockcroft-Gault equation (19) and normalized to the body surface",
        "area'. BSA normalization was applied deliberately 'in order to avoid",
        "confounding between body size and renal function' (Results,",
        "'Population pharmacokinetic model development'). Drives a power",
        "function on the RENAL clearance arm only (Table 1 row",
        "'CL R, CLcr power'); the non-renal arm is CLcr-independent. Baseline",
        "CLcr was retained over time-varying CLcr because time-varying CLcr",
        "gave no improvement in the objective function (same subsection).",
        "Reference value 100 mL/min/1.73 m^2, the reference-population normal",
        "renal function stated in Cammarata 2024 Results ('have normal renal",
        "function (CLcr of 100 mL/min/1.73 m^2)'). Pooled cohort median",
        "91.5 mL/min/1.73 m^2 (range 5.61-364; Table S3)."
      ),
      source_name        = "CLcr"
    ),
    RENALIMP_SEV = list(
      description        = paste(
        "Severe renal impairment indicator",
        "(1 = baseline CLcr < 30 mL/min/1.73 m^2; 0 otherwise)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (baseline CLcr >= 30 mL/min/1.73 m^2)",
      notes              = paste(
        "Cammarata 2024 Results: 'For both the durlobactam and sulbactam",
        "models, patients with CLcr < 30 mL/min/1.73 m^2 were found to",
        "produce significantly different results from other patients with",
        "higher CLcr values. Total CL was, therefore, adjusted for patients",
        "with CLcr < 30 mL/min/1.73 m^2 using a proportional shift.' The",
        "shift is applied to TOTAL CL (renal + non-renal), on top of the",
        "continuous CRCL power function on the renal arm, and is the Table 1",
        "row BCLCRNLT30. The severity threshold here is the paper's own",
        "CLcr < 30 mL/min/1.73 m^2 cut-off, matching the FDA / EMA",
        "renal-impairment convention noted in the canonical register.",
        "Derived from CRCL; supply as 1 * (CRCL < 30)."
      ),
      source_name        = "BCLCRNLT30"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = paste(
        "Hemodialysis-session gate",
        "(1 while an intermittent hemodialysis session is running; 0 otherwise)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic interval, or a non-dialysed subject)",
      notes              = paste(
        "Time-varying WITHIN subject. Cammarata 2024 Results, 'Hemodialysis",
        "sub-models': 'an HD effect (HDEFFECT) term with IIV was applied to",
        "CL for each of the two analyte sub-models ... The CL-HDEFFECT",
        "parameters were simply multiplicative terms with IIV applied to them",
        "so they varied by subject.' Encoded here as a log-scale",
        "multiplicative factor on TOTAL CL that is switched on only while the",
        "gate is 1, so the interdialytic clearance reduces exactly to the",
        "plasma-model clearance and non-dialysed subjects are unaffected.",
        "The illustrative simulation reported in Table S5 uses a 4-hour",
        "session starting 1 hour after the end of the morning 3-hour infusion",
        "(i.e. RRT_HEMODIAL_ACTIVE = 1 from t = 4 h to t = 8 h after the start",
        "of the morning infusion).",
        "Note on the column choice: the sidecar option ratified for this",
        "extraction described 'a time-varying HD flag' but named the",
        "subject-level canonical RRT_HEMODIAL_STATUS. The canonical register",
        "reserves RRT_HEMODIAL_STATUS for the subject-level treatment-status",
        "indicator and RRT_HEMODIAL_ACTIVE for exactly this per-session gate,",
        "so RRT_HEMODIAL_ACTIVE is used here to match the ratified semantics."
      ),
      source_name        = "HD (on/off during a session)"
    ),
    REGION_EASTASIA = list(
      description        = paste(
        "East Asian region of origin indicator",
        "(1 = enrolled in China, Taiwan, or South Korea; 0 otherwise)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-East-Asian region)",
      notes              = paste(
        "Cammarata 2024 Discussion / Results: 'the East Asian region was",
        "identified as a statistically significant predictor of the",
        "variability in CL and Vc for durlobactam (but not sulbactam)', with",
        "the region 'defined as patients from China. Taiwan, or South Korea'.",
        "Proportional shifts on durlobactam CL (Table 1 CLEASIAFL1) and",
        "durlobactam Vc (Table 1 V1EASIAFL1); no sulbactam parameter carries",
        "a region effect. Note the paper distinguishes region from race and",
        "country of origin: 'Neither race nor country of origin was",
        "identified as a statistically significant predictor'. 45 of 373",
        "pooled subjects (12.1%) were East Asian (Table S3, 'Region')."
      ),
      source_name        = "EASIAFL"
    ),
    DIS_HABP = list(
      description        = "Hospital-acquired bacterial pneumonia cohort indicator (1 = HABP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Cammarata 2024 Table 1 level INFTYPN1. Proportional shift on",
        "sulbactam CL (CLINFTYPN1 = -0.424) and sulbactam Vc",
        "(V3INFTYPN1 = 0.836). For durlobactam the HABP and VABP Vc terms",
        "were merged into a single coefficient because 'they did not yield",
        "significantly different effects on PK from one another' (Results,",
        "covariate-analysis paragraph), so durlobactam uses e_habp_vabp_vc",
        "applied to DIS_HABP + DIS_VABP. The five infection-type indicators",
        "are mutually exclusive; the shared reference category (all five 0)",
        "is the uninfected Phase 1 healthy-subject / renal-impairment cohort.",
        "38 subjects in the pooled data set had HABP (Table S12)."
      ),
      source_name        = "INFTYPN = 1"
    ),
    DIS_VABP = list(
      description        = "Ventilator-associated bacterial pneumonia cohort indicator (1 = VABP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Cammarata 2024 Table 1 level INFTYPN2. Proportional shift on",
        "sulbactam CL (CLINFTYPN2 = -0.298) and sulbactam Vc",
        "(V3INFTYPN2 = 1.43). For durlobactam Vc it shares the merged",
        "HABP-and-VABP coefficient V1INFTYPN1&2 = 1.52. 56 subjects in the",
        "pooled data set had VABP (Table S12)."
      ),
      source_name        = "INFTYPN = 2"
    ),
    DIS_CUTI = list(
      description        = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Cammarata 2024 Table 1 level INFTYPN3. Proportional shift on",
        "durlobactam Vc (V1INFTYPN3 = 0.343), sulbactam CL",
        "(CLINFTYPN3 = -0.157), and sulbactam Vc (V3INFTYPN3 = 0.17).",
        "cUTI patients came from the Phase 2 study CS2514-2017-0003, in",
        "which acute pyelonephritis (DIS_AP) was a separate reported",
        "infection-type level. 35 subjects in the pooled data set had cUTI",
        "(Table S12)."
      ),
      source_name        = "INFTYPN = 3"
    ),
    DIS_BACTEREMIA = list(
      description        = "Bacteremia / bloodstream-infection cohort indicator (1 = bacteremia)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Cammarata 2024 Table 1 level INFTYPN4. Proportional shift on",
        "durlobactam Vc (V1INFTYPN4 = 3.32, the largest single covariate",
        "effect in the model), sulbactam CL (CLINFTYPN4 = -0.444), and",
        "sulbactam Vc (V3INFTYPN4 = 1.85). Bacteremia due to Acinetobacter",
        "baumannii-calcoaceticus complex was one of the Phase 3",
        "(CS2514-2017-0004) enrollment categories; 16 subjects in the pooled",
        "data set had bacteremia (Table S12)."
      ),
      source_name        = "INFTYPN = 4"
    ),
    DIS_AP = list(
      description        = "Acute pyelonephritis cohort indicator (1 = AP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Cammarata 2024 Table 1 level INFTYPN5 (table footnote 'AP, acute",
        "pyelonephritis'). Proportional shift on sulbactam CL",
        "(CLINFTYPN5 = -0.382) and sulbactam Vc (V3INFTYPN5 = -0.704).",
        "Durlobactam carries no AP term - no V1INFTYPN5 row appears in the",
        "durlobactam half of Table 1, and the missing coefficient must NOT be",
        "borrowed from the cUTI term. AP patients came from the Phase 2",
        "cUTI/AP study CS2514-2017-0003; 17 subjects in the pooled data set",
        "had AP (Table S12)."
      ),
      source_name        = "INFTYPN = 5"
    ),
    STUDY_SULDUR_PHASE2 = list(
      description        = "Phase 2 study cohort indicator (1 = Study CS2514-2017-0003)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 studies, when STUDY_SULDUR_PHASE3 is also 0)",
      notes              = paste(
        "Selects the durlobactam Phase 2 proportional residual magnitude.",
        "Cammarata 2024 Results: 'Residual variability was described for",
        "durlobactam as a series of proportional and additive error models",
        "for plasma concentrations obtained from subjects enrolled in Phase",
        "1, 2, and 3 studies, respectively.' Table 1 reports three separate",
        "durlobactam proportional sigma^2 terms (Phase 1 / 2 / 3) and a",
        "single additive term for Phase 1 only; the Phase 3 additive",
        "component 'was determined to not be significant and was consequently",
        "removed'. Sulbactam residual variability is NOT phase-stratified, so",
        "this indicator does not enter the sulbactam error model, and neither",
        "phase indicator enters the ELF error model (a single proportional",
        "term per analyte; Table S6). Paired with STUDY_SULDUR_PHASE3; both 0",
        "selects Phase 1."
      ),
      source_name        = "study phase"
    ),
    STUDY_SULDUR_PHASE3 = list(
      description        = "Phase 3 study cohort indicator (1 = Study CS2514-2017-0004)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 studies, when STUDY_SULDUR_PHASE2 is also 0)",
      notes              = paste(
        "Selects the durlobactam Phase 3 proportional residual magnitude",
        "(Cammarata 2024 Table 1, 'sigma^2 plasma, Proportional Phase 3' =",
        "0.203, 45.0 %CV). No additive term applies to Phase 3 data. See the",
        "STUDY_SULDUR_PHASE2 notes for the full residual-error rationale."
      ),
      source_name        = "study phase"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 373L,
    n_studies        = 8L,
    n_concentrations = 5188L,
    age_range        = "18-91 years",
    age_median       = "46 years",
    weight_range     = "35.8-150 kg",
    weight_median    = "75 kg",
    sex_female_pct   = 37.5,
    race_ethnicity   = c(
      White = 65.4, Black = 13.1, Asian = 16.1,
      `American Indian/Alaska Native` = 1.88, Other = 3.49
    ),
    disease_state    = paste(
      "Pooled healthy adults (six Phase 1 studies including a dedicated",
      "renal-impairment study spanning normal renal function to end-stage",
      "renal disease on hemodialysis), adults with complicated urinary tract",
      "infection or acute pyelonephritis (Phase 2 Study CS2514-2017-0003),",
      "and adults with hospital-acquired or ventilator-associated bacterial",
      "pneumonia or bacteremia caused by Acinetobacter",
      "baumannii-calcoaceticus complex (Phase 3 Study CS2514-2017-0004).",
      "Infection-type counts in the pooled data set (Table S12): HABP 38,",
      "VABP 56, cUTI 35, AP 17, bacteremia 16."
    ),
    dose_range       = paste(
      "Durlobactam single doses 0.25-8 g and 0.25-2 g q6h; combination",
      "dosing 0.5 g or 1 g sulbactam with 0.5 g or 1 g durlobactam, given as",
      "3-hour intravenous infusions every 6 hours (one Phase 1 cohort",
      "received durlobactam over 2 hours). Subjects with CLcr < 30",
      "mL/min/1.73 m^2 received extended 8-12 hour dosing intervals;",
      "augmented renal function (CLcr > 130 mL/min/1.73 m^2) received",
      "1.5 g/1.5 g q6h. Study details in Table S1."
    ),
    regions          = "United States and other non-East-Asian sites; China (Studies ZL-2402-001 and part of CS2514-2017-0004), Taiwan, and South Korea",
    renal_function   = paste(
      "Baseline BSA-normalized Cockcroft-Gault CLcr 5.61-364 mL/min/1.73 m^2",
      "(pooled median 91.5). The renal-impairment study CS2514-2017-0002",
      "enrolled cohorts spanning normal renal function, mild, moderate, and",
      "severe impairment, plus an end-stage renal disease hemodialysis",
      "cohort (Table S1 footnote d)."
    ),
    bmi_range        = "11-52.1 kg/m^2 (pooled median 25.7)",
    notes            = paste(
      "Baseline demographics from Cammarata 2024 Table S3 (pooled dataset,",
      "n = 373). 432 subjects and 8,100 plasma concentrations entered the",
      "analysis data set; after exclusion of records with missing dose or",
      "sample times or missing covariates and of 48 Phase 3 outliers, the",
      "final plasma model data set held 373 subjects and 5,188 concentration",
      "records (3,494 durlobactam from 373 subjects; 1,944 sulbactam from",
      "264 subjects). The hemodialysis sub-model was fitted separately to 202",
      "Period 2 records from six end-stage renal disease subjects in Cohort 5",
      "of Study CS2514-2017-0002, and the epithelial lining fluid sub-model",
      "to 60 ELF records from 30 healthy subjects in Study CS2514-2017-0001;",
      "in both sub-models the plasma terms were fixed to the Table 1",
      "population means, so their estimated terms compose directly onto the",
      "plasma model as encoded here."
    )
  )

  ini({
    # =====================================================================
    # DURLOBACTAM structural parameters (Cammarata 2024 Table 1,
    # 'Final model / Estimate' column, Durlobactam block). Values refer to
    # the reference subject: CRCL = 100 mL/min/1.73 m^2, WT = 75 kg,
    # non-East-Asian region, no infection (Phase 1 healthy subject),
    # CLcr >= 30 mL/min/1.73 m^2, no hemodialysis running.
    #
    # Total CL is the sum of a renal and a non-renal arm, split by the
    # fraction excreted renally (FE) and with only the renal arm scaled by
    # baseline CLcr (Cammarata 2024 Abstract and Results: 'Both renal
    # clearance and nonrenal clearance were estimated, and total clearance
    # was calculated as the sum of renal and nonrenal clearance. Individual
    # renal clearances were scaled by baseline creatinine clearance.').
    # =====================================================================
    lcl <- log(9.33); label("Durlobactam total CL at CRCL = 100 mL/min/1.73 m^2 (L/h)")   # Table 1 Durlobactam: CL = 9.33 L/h (%SEM 3.24; SIR 90% CI [9.07, 9.65])
    lvc <- log(12.5); label("Durlobactam central volume Vc at WT = 75 kg, no infection (L)") # Table 1 Durlobactam: Vc = 12.5 L (%SEM 2.93; SIR 90% CI [12.1, 12.9])
    lq  <- log(4.43); label("Durlobactam inter-compartmental clearance Q (L/h)")           # Table 1 Durlobactam: Q = 4.43 L/h (%SEM 3.77; SIR 90% CI [4.17, 4.71])
    lvp <- log(5.83); label("Durlobactam peripheral volume Vp (L)")                        # Table 1 Durlobactam: Vp = 5.83 L (%SEM 3.44; SIR 90% CI [5.54, 6.10])

    f_renal <- fixed(0.479); label("Durlobactam fraction of total CL excreted renally (FE, unitless; FIXED)")
    # Table 1 Durlobactam: FE (%) = 0.479, reported in the shaded 'fixed or
    # not estimated' style with no %SEM. Results: 'The fraction of durlobactam
    # and sulbactam assumed to be excreted renally (FE) was fixed to the values
    # that had been used during the previous population PK analyses conducted
    # for sulbactam-durlobactam (0.479 for both drugs) (11).' See the vignette
    # Errata for the later FE revision (0.660 durlobactam / 0.648 sulbactam),
    # whose refitted parameter table is not reported in the paper.

    e_crcl_cl_renal <- 0.875; label("Durlobactam CLcr power exponent on the renal CL arm (unitless)") # Table 1 Durlobactam: 'CL R, CLcr power' = 0.875 (%SEM 12; SIR 90% CI [0.766, 0.988])
    e_wt_cl         <- 0.646; label("Durlobactam body-weight power exponent on total CL (unitless)")  # Table 1 Durlobactam: CLWTKG1 = 0.646 (%SEM 19.4; SIR 90% CI [0.527, 0.783])
    e_wt_vc         <- 0.521; label("Durlobactam body-weight power exponent on Vc (unitless)")        # Table 1 Durlobactam: V1WTKG1 = 0.521 (%SEM 27.0; SIR 90% CI [0.359, 0.704])

    e_renalimp_sev_cl     <- -0.58;  label("Durlobactam proportional shift in total CL for CLcr < 30 mL/min/1.73 m^2 (fraction)")  # Table 1 Durlobactam: BCLCRNLT30 = -0.58 (%SEM 8.15; SIR 90% CI [-0.617, -0.547])
    e_region_eastasia_cl  <- -0.199; label("Durlobactam proportional shift in total CL for East Asian region (fraction)")          # Table 1 Durlobactam: CLEASIAFL1 = -0.199 (%SEM 31.7; SIR 90% CI [-0.282, -0.117])
    e_region_eastasia_vc  <- -0.263; label("Durlobactam proportional shift in Vc for East Asian region (fraction)")                # Table 1 Durlobactam: V1EASIAFL1 = -0.263 (%SEM 49.2; SIR 90% CI [-0.369, -0.147])

    e_habp_vabp_vc   <- 1.52;  label("Durlobactam proportional shift in Vc for HABP or VABP (fraction)")  # Table 1 Durlobactam: V1INFTYPN1&2 (HABP and VABP) = 1.52 (%SEM 26.5; SIR 90% CI [1.12, 2.01]); the HABP and VABP terms were merged during model refinement
    e_cuti_vc        <- 0.343; label("Durlobactam proportional shift in Vc for cUTI (fraction)")          # Table 1 Durlobactam: V1INFTYPN3 (cUTI) = 0.343 (%SEM 62.3; SIR 90% CI [0.165, 0.532])
    e_bacteremia_vc  <- 3.32;  label("Durlobactam proportional shift in Vc for bacteremia (fraction)")    # Table 1 Durlobactam: V1INFTYPN4 (bacteremia) = 3.32 (%SEM 53.7; SIR 90% CI [2.05, 4.52])

    # =====================================================================
    # SULBACTAM structural parameters (Cammarata 2024 Table 1, Sulbactam
    # block). Same reference subject as durlobactam.
    # =====================================================================
    lcl_sbt <- log(13.5); label("Sulbactam total CL at CRCL = 100 mL/min/1.73 m^2 (L/h)")    # Table 1 Sulbactam: CL = 13.5 L/h (%SEM 14.0; SIR 90% CI [12.3, 14.6])
    lvc_sbt <- log(12.0); label("Sulbactam central volume Vc at WT = 75 kg, no infection (L)") # Table 1 Sulbactam: Vc = 12 L (%SEM 8.90; SIR 90% CI [10.9, 13.0])
    lq_sbt  <- log(7.88); label("Sulbactam inter-compartmental clearance Q (L/h)")            # Table 1 Sulbactam: Q = 7.88 L/h (%SEM 18.9; SIR 90% CI [6.56, 9.72])
    lvp_sbt <- log(6.99); label("Sulbactam peripheral volume Vp (L)")                         # Table 1 Sulbactam: Vp = 6.99 L (%SEM 9.29; SIR 90% CI [6.24, 7.80])

    f_renal_sbt <- fixed(0.479); label("Sulbactam fraction of total CL excreted renally (FE, unitless; FIXED)")
    # Table 1 Sulbactam: FE (%) = 0.479, shaded 'fixed or not estimated'.

    e_crcl_cl_renal_sbt <- 1.14; label("Sulbactam CLcr power exponent on the renal CL arm (unitless)") # Table 1 Sulbactam: 'CL R, CLcr power' = 1.14 (%SEM 20.8; SIR 90% CI [0.920, 1.37])
    e_wt_cl_sbt         <- 1.01; label("Sulbactam body-weight power exponent on total CL (unitless)")  # Table 1 Sulbactam: CLWTKG1 = 1.01 (%SEM 26.9; SIR 90% CI [0.795, 1.24])
    e_wt_vc_sbt         <- 0.831; label("Sulbactam body-weight power exponent on Vc (unitless)")       # Table 1 Sulbactam: V3WTKG1 = 0.831 (%SEM 26.7; SIR 90% CI [0.573, 1.11])

    e_renalimp_sev_cl_sbt <- -0.635; label("Sulbactam proportional shift in total CL for CLcr < 30 mL/min/1.73 m^2 (fraction)") # Table 1 Sulbactam: BCLCRNLT30 = -0.635 (%SEM 11; SIR 90% CI [-0.682, -0.578])

    e_habp_cl_sbt       <- -0.424; label("Sulbactam proportional shift in total CL for HABP (fraction)")       # Table 1 Sulbactam: CLINFTYPN1 (HABP) = -0.424 (%SEM 19.4; SIR 90% CI [-0.494, -0.346])
    e_vabp_cl_sbt       <- -0.298; label("Sulbactam proportional shift in total CL for VABP (fraction)")       # Table 1 Sulbactam: CLINFTYPN2 (VABP) = -0.298 (%SEM 36.5; SIR 90% CI [-0.393, -0.195])
    e_cuti_cl_sbt       <- -0.157; label("Sulbactam proportional shift in total CL for cUTI (fraction)")       # Table 1 Sulbactam: CLINFTYPN3 (cUTI) = -0.157 (%SEM 96.8; SIR 90% CI [-0.274, -0.00759])
    e_bacteremia_cl_sbt <- -0.444; label("Sulbactam proportional shift in total CL for bacteremia (fraction)") # Table 1 Sulbactam: CLINFTYPN4 (bacteremia) = -0.444 (%SEM 20.9; SIR 90% CI [-0.537, -0.333])
    e_ap_cl_sbt         <- -0.382; label("Sulbactam proportional shift in total CL for acute pyelonephritis (fraction)") # Table 1 Sulbactam: CLINFTYPN5 (AP) = -0.382 (%SEM 28.8; SIR 90% CI [-0.489, -0.253])

    e_habp_vc_sbt       <-  0.836; label("Sulbactam proportional shift in Vc for HABP (fraction)")       # Table 1 Sulbactam: V3INFTYPN1 (HABP) = 0.836 (%SEM 25.3; SIR 90% CI [0.554, 1.15])
    e_vabp_vc_sbt       <-  1.43;  label("Sulbactam proportional shift in Vc for VABP (fraction)")       # Table 1 Sulbactam: V3INFTYPN2 (VABP) = 1.43 (%SEM 19.5; SIR 90% CI [1.04, 1.84])
    e_cuti_vc_sbt       <-  0.17;  label("Sulbactam proportional shift in Vc for cUTI (fraction)")       # Table 1 Sulbactam: V3INFTYPN3 (cUTI) = 0.17 (%SEM 89.8; SIR 90% CI [-0.00793, 0.362])
    e_bacteremia_vc_sbt <-  1.85;  label("Sulbactam proportional shift in Vc for bacteremia (fraction)") # Table 1 Sulbactam: V3INFTYPN4 (bacteremia) = 1.85 (%SEM 44.7; SIR 90% CI [0.911, 2.96])
    e_ap_vc_sbt         <- -0.704; label("Sulbactam proportional shift in Vc for acute pyelonephritis (fraction)") # Table 1 Sulbactam: V3INFTYPN5 (AP) = -0.704 (%SEM 13.2; SIR 90% CI [-0.824, -0.579])

    # =====================================================================
    # HEMODIALYSIS SUB-MODEL (Cammarata 2024 Table S4). The sub-model was
    # fitted separately with every plasma term fixed to its Table 1
    # population mean, so its two estimated terms compose directly onto the
    # plasma model above.
    #
    # HDEFFECT is a multiplicative fold-change on TOTAL CL with its own IIV,
    # i.e. NONMEM HDEFFECT = THETA * EXP(ETA). It is stored here on the log
    # scale so that the same THETA * EXP(ETA) form is reproduced exactly by
    # exp(e_hemodial_active_cl + eta), gated by the time-varying
    # RRT_HEMODIAL_ACTIVE so that clearance reduces exactly to the plasma
    # value between sessions and in non-dialysed subjects.
    # =====================================================================
    e_hemodial_active_cl     <- log(6.24); label("Log fold-increase in durlobactam total CL while hemodialysis is running (unitless)") # Table S4: CL-HDEFFECT (durlobactam) = 6.24 (%SEM 14.5); Results: 'the effects of HD with durlobactam and sulbactam caused an increase in CL by 6.24-fold and 8.19-fold, respectively'
    e_hemodial_active_cl_sbt <- log(8.19); label("Log fold-increase in sulbactam total CL while hemodialysis is running (unitless)")   # Table S4: CL-HDEFFECT (sulbactam) = 8.19 (%SEM 23.2)

    # =====================================================================
    # EPITHELIAL LINING FLUID SUB-MODEL (Cammarata 2024 Table S6). Also
    # fitted with the plasma terms fixed. ELF concentration is an
    # instantaneous ratio of the plasma concentration: 'the ELF sub-models
    # assumed that ELF concentrations were simply a ratio of the plasma
    # concentrations (i.e., with no time lag)' (Results).
    #
    # These ratios are relative to TOTAL-drug plasma concentrations. The
    # paper also reports the corresponding free-drug ratios, obtained by
    # dividing by the unbound fraction implied by protein binding of 10%
    # (durlobactam) and 38% (sulbactam): 0.372 / 0.90 = 0.413 and
    # 0.533 / 0.62 = 0.860. Protein binding itself is not a fitted parameter
    # of this model, so only the total-drug ratios are encoded.
    # =====================================================================
    lrelf     <- log(0.372); label("Durlobactam ELF / total-drug plasma concentration ratio (unitless)") # Table S6: PLASMA-ELF ratio (durlobactam) = 0.372 (%SEM 3.6); Results: 'ELF penetration relative to total-drug plasma concentrations for durlobactam was 37.2%'
    lrelf_sbt <- log(0.533); label("Sulbactam ELF / total-drug plasma concentration ratio (unitless)")   # Table S6: PLASMA-ELF ratio (sulbactam) = 0.533 (%SEM 5.41); Results: 'the ELF penetration for sulbactam was 53.3%'

    # =====================================================================
    # Inter-individual variability (Cammarata 2024 Table 1 for the plasma
    # terms, Table S4 for the hemodialysis terms). The paper reports omega^2
    # directly (the parenthetical %CV is sqrt(omega^2) x 100, not the
    # log-normal sqrt(exp(omega^2) - 1)), so the variances below are the
    # published values verbatim. IIV was estimated on CL, Vc, and Vp of both
    # drugs, 'with off-diagonal relationships on CL and Vc for both
    # sulbactam and durlobactam' (Results). Off-diagonals between the two
    # drugs were not reported and are therefore not encoded, and Table S6
    # reports no IIV on the ELF ratios.
    # =====================================================================
    etalcl + etalvc ~ c(0.0778,
                        0.0494, 0.0757)
    # Table 1 Durlobactam: omega^2 CL = 0.0778 (27.9 %CV, %SEM 9.56, shrinkage 9.13);
    # omega^2 Vc = 0.0757 (27.5 %CV, %SEM 13.6, shrinkage 25.7);
    # Covariance(omega^2 CL, omega^2 Vc) = 0.0494 (r^2 = 0.415, %SEM 17.8)
    etalvp ~ 0.0773
    # Table 1 Durlobactam: omega^2 Vp = 0.0773 (27.8 %CV, %SEM 9.83, shrinkage 37.6)

    etalcl_sbt + etalvc_sbt ~ c(0.221,
                                0.0727, 0.0967)
    # Table 1 Sulbactam: omega^2 CL = 0.221 (47.0 %CV, %SEM 10.6, shrinkage 19.2);
    # omega^2 Vc = 0.0967 (31.1 %CV, %SEM 24.5, shrinkage 44.7);
    # Covariance(omega^2 CL, omega^2 Vc) = 0.0727 (r^2 = 0.2448, %SEM 32.6)
    etalvp_sbt ~ 0.196
    # Table 1 Sulbactam: omega^2 Vp = 0.196 (44.3 %CV, %SEM 25.7, shrinkage 50.0)

    etae_hemodial_active_cl     ~ 0.124
    # Table S4: omega^2 CL-HDEFFECT (durlobactam) = 0.124 (35.2 %CV, %SEM 40.9, shrinkage 57.8)
    etae_hemodial_active_cl_sbt ~ 0.316
    # Table S4: omega^2 CL-HDEFFECT (sulbactam) = 0.316 (56.2 %CV, %SEM 52.3, shrinkage 58.1)

    # =====================================================================
    # Residual variability (Cammarata 2024 Table 1 for plasma, Table S6 for
    # ELF). Both tables report sigma^2; the SDs below are sqrt(sigma^2) and
    # reproduce the parenthetical %CV / mg/L values printed in the tables.
    #
    # Durlobactam plasma: separate proportional terms per study phase, with
    # an additive term for Phase 1 only. Sulbactam plasma: a single
    # proportional plus additive model for all data. ELF: one proportional
    # term per analyte, not phase-stratified.
    # =====================================================================
    propSdPhase1 <- sqrt(0.019);   label("Durlobactam proportional residual SD, Phase 1 plasma (fraction)")   # Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 1 = 0.019 (13.8 %CV, %SEM 1.90, shrinkage 8.67)
    addSdPhase1  <- sqrt(0.00136); label("Durlobactam additive residual SD, Phase 1 plasma (mg/L)")           # Table 1 Durlobactam: sigma^2 plasma, Additive Phase 1 = 0.00136 (0.0369 mg/L, %SEM 3.69, shrinkage 8.67)
    propSdPhase2 <- sqrt(0.0794);  label("Durlobactam proportional residual SD, Phase 2 plasma (fraction)")   # Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 2 = 0.0794 (28.2 %CV, %SEM 13.3, shrinkage 12.1)
    propSdPhase3 <- sqrt(0.203);   label("Durlobactam proportional residual SD, Phase 3 plasma (fraction)")   # Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 3 = 0.203 (45.0 %CV, %SEM 9.85, shrinkage 8.61); the Phase 3 additive component was tested, found not significant, and removed

    propSd_sbt <- sqrt(0.0433); label("Sulbactam proportional residual SD, plasma (fraction)") # Table 1 Sulbactam: sigma^2 plasma, Proportional = 0.0433 (20.8 %CV, %SEM 2.54, shrinkage 11.6)
    addSd_sbt  <- sqrt(0.0054); label("Sulbactam additive residual SD, plasma (mg/L)")         # Table 1 Sulbactam: sigma^2 plasma, Additive = 0.0054 (0.0735 mg/L, %SEM 6.38, shrinkage 11.6)

    propSd_Celf     <- sqrt(0.0322); label("Durlobactam proportional residual SD, ELF (fraction)") # Table S6: sigma^2 ELF, Proportional (durlobactam) = 0.0322 (17.9 %CV, %SEM 31.8, shrinkage 3.91)
    propSd_Celf_sbt <- sqrt(0.0628); label("Sulbactam proportional residual SD, ELF (fraction)")   # Table S6: sigma^2 ELF, Proportional (sulbactam) = 0.0628 (25.1 %CV, %SEM 32.6, shrinkage 2.56)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariate terms.
    #
    # The five infection-type indicators are mutually exclusive, so the
    # proportional shifts are summed inside a single (1 + ...) bracket; the
    # shared reference category is the uninfected Phase 1 subject.
    # ------------------------------------------------------------------
    infect_vc     <- 1 + e_habp_vabp_vc  * (DIS_HABP + DIS_VABP) +
                         e_cuti_vc       * DIS_CUTI +
                         e_bacteremia_vc * DIS_BACTEREMIA

    infect_cl_sbt <- 1 + e_habp_cl_sbt       * DIS_HABP +
                         e_vabp_cl_sbt       * DIS_VABP +
                         e_cuti_cl_sbt       * DIS_CUTI +
                         e_bacteremia_cl_sbt * DIS_BACTEREMIA +
                         e_ap_cl_sbt         * DIS_AP

    infect_vc_sbt <- 1 + e_habp_vc_sbt       * DIS_HABP +
                         e_vabp_vc_sbt       * DIS_VABP +
                         e_cuti_vc_sbt       * DIS_CUTI +
                         e_bacteremia_vc_sbt * DIS_BACTEREMIA +
                         e_ap_vc_sbt         * DIS_AP

    # Hemodialysis fold-change on total CL. The gate is time-varying, so the
    # factor collapses to exp(0) = 1 between sessions and in every subject
    # who is not on dialysis.
    hdeffect     <- exp((e_hemodial_active_cl     + etae_hemodial_active_cl)     * RRT_HEMODIAL_ACTIVE)
    hdeffect_sbt <- exp((e_hemodial_active_cl_sbt + etae_hemodial_active_cl_sbt) * RRT_HEMODIAL_ACTIVE)

    # ------------------------------------------------------------------
    # 2. Durlobactam individual PK parameters.
    #
    #   CL_R  = CL * FE       * (CLcr / 100)^theta_CLcr
    #   CL_NR = CL * (1 - FE)
    #   CL_T  = (CL_R + CL_NR) * (WT / 75)^theta_WT
    #                          * (1 + theta_EASIA * EASIAFL)
    #                          * (1 + theta_CLCRLT30 * (CLcr < 30))
    #                          * exp(eta_CL) * HDEFFECT^HD
    # ------------------------------------------------------------------
    cl_renal  <- exp(lcl) * f_renal * (CRCL / 100)^e_crcl_cl_renal
    cl_nonren <- exp(lcl) * (1 - f_renal)

    cl <- (cl_renal + cl_nonren) *
      (WT / 75)^e_wt_cl *
      (1 + e_region_eastasia_cl * REGION_EASTASIA) *
      (1 + e_renalimp_sev_cl * RENALIMP_SEV) *
      exp(etalcl) *
      hdeffect

    vc <- exp(lvc + etalvc) *
      (WT / 75)^e_wt_vc *
      (1 + e_region_eastasia_vc * REGION_EASTASIA) *
      infect_vc

    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # ------------------------------------------------------------------
    # 3. Sulbactam individual PK parameters (same construction; sulbactam
    #    carries no region effect but does carry infection-type effects on
    #    both CL and Vc).
    # ------------------------------------------------------------------
    cl_renal_sbt  <- exp(lcl_sbt) * f_renal_sbt * (CRCL / 100)^e_crcl_cl_renal_sbt
    cl_nonren_sbt <- exp(lcl_sbt) * (1 - f_renal_sbt)

    cl_sbt <- (cl_renal_sbt + cl_nonren_sbt) *
      (WT / 75)^e_wt_cl_sbt *
      infect_cl_sbt *
      (1 + e_renalimp_sev_cl_sbt * RENALIMP_SEV) *
      exp(etalcl_sbt) *
      hdeffect_sbt

    vc_sbt <- exp(lvc_sbt + etalvc_sbt) *
      (WT / 75)^e_wt_vc_sbt *
      infect_vc_sbt

    q_sbt  <- exp(lq_sbt)
    vp_sbt <- exp(lvp_sbt + etalvp_sbt)

    # ------------------------------------------------------------------
    # 4. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_sbt <- cl_sbt / vc_sbt
    k12_sbt <- q_sbt  / vc_sbt
    k21_sbt <- q_sbt  / vp_sbt

    # ------------------------------------------------------------------
    # 5. Four-compartment system: two independent two-compartment IV
    #    dispositions (the two drugs do not interconvert). Doses are given
    #    as intravenous infusions into the respective central compartments.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    d/dt(central_sbt)     <- -(kel_sbt + k12_sbt) * central_sbt + k21_sbt * peripheral1_sbt
    d/dt(peripheral1_sbt) <-   k12_sbt * central_sbt            - k21_sbt * peripheral1_sbt

    # ------------------------------------------------------------------
    # 6. Observations. Dose in mg, volumes in L -> concentrations in mg/L.
    #    ELF tracks plasma instantaneously (no time lag).
    #
    #    Durlobactam plasma residual variability is switched by study phase;
    #    both phase indicators 0 selects Phase 1 (proportional + additive).
    # ------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_sbt <- central_sbt / vc_sbt

    Celf     <- Cc     * exp(lrelf)
    Celf_sbt <- Cc_sbt * exp(lrelf_sbt)

    phase1 <- 1 - STUDY_SULDUR_PHASE2 - STUDY_SULDUR_PHASE3
    propSd <- propSdPhase1 * phase1 +
              propSdPhase2 * STUDY_SULDUR_PHASE2 +
              propSdPhase3 * STUDY_SULDUR_PHASE3
    addSd  <- addSdPhase1 * phase1

    Cc       ~ add(addSd)     + prop(propSd)
    Cc_sbt   ~ add(addSd_sbt) + prop(propSd_sbt)
    Celf     ~ prop(propSd_Celf)
    Celf_sbt ~ prop(propSd_Celf_sbt)
  })
}
