Larson_2026_sulbactam_durlobactam_pediatric_allometry <- function() {
  description <- paste(
    "Pediatric 'Adult with Allometry' adaptation of the adult",
    "sulbactam-durlobactam population PK model, used to select the dosing",
    "regimens evaluated in the Phase 1b pediatric trial NCT06801223 (Larson",
    "2026 / Cammarata IDWeek 2025 poster P-444). The adult four-compartment",
    "(two compartments per drug) model of Cammarata 2024 is carried over in",
    "full - renal plus non-renal clearance arms, baseline CLcr power function",
    "on the renal arm, severe-renal-impairment shift, infection-type and East",
    "Asian region proportional shifts, hemodialysis gate, and epithelial",
    "lining fluid ratios - with one change: the estimated adult body-weight",
    "power exponents are replaced by fixed allometric exponents of 0.75 on",
    "every clearance term and 1.0 on every volume term, and body weight is",
    "additionally applied to Q and Vp of both drugs, which carried no weight",
    "effect in the adult model. This is the more conservative of the two",
    "approaches the poster explored and is the one that guided dose selection.",
    "Durlobactam uses the unsuffixed canonical compartment / parameter set;",
    "sulbactam carries the sibling-drug suffix _sbt throughout. See",
    "modellib('Larson_2026_sulbactam_durlobactam_pediatric_allometry_crcl')",
    "for the alternative 'Allometry + CLCR' approach and",
    "modellib('Cammarata_2024_sulbactam_durlobactam') for the adult parent",
    "model."
  )
  reference <- paste(
    "Larson KB, O'Donnell J, Tanudra A, Cammarata AP, Rubino CM. P-444.",
    "Population Pharmacokinetics (PPK) Analysis of Sulbactam-Durlobactam",
    "(SUD) to Support Dose Selection for Evaluation in a Clinical Trial in",
    "Pediatric Patients with Acinetobacter Baumannii-Calcoaceticus Complex",
    "(ABC) Infections. Open Forum Infect Dis. 2026;13(Suppl 1):S405.",
    "doi:10.1093/ofid/ofaf695.659. PMCID: PMC12792777.",
    "The model specification is taken from the corresponding IDWeek 2025",
    "poster, which carries a different author order: Cammarata A, Larson KB,",
    "Tanudra A, O'Donnell JP, Bhavnani SM, Rubino CM. 'Population",
    "pharmacokinetics analysis of sulbactam-durlobactam to support the dose",
    "selection for evaluation in a clinical trial in pediatric patients with",
    "Acinetobacter baumannii-calcoaceticus complex infections.' Poster P-444,",
    "IDWeek 2025, Atlanta, GA.",
    "All structural, covariate, IIV and residual-error values are inherited",
    "unchanged from the adult parent model: Cammarata AP, Safir MC, Trang M,",
    "Larson KB, O'Donnell JP, Bhavnani SM, Rubino CM. Population",
    "pharmacokinetic analyses for sulbactam-durlobactam using Phase 1, 2, and",
    "3 data. Antimicrob Agents Chemother. 2025;69(1):e00485-24.",
    "doi:10.1128/aac.00485-24; see",
    "modellib('Cammarata_2024_sulbactam_durlobactam').",
    sep = " "
  )
  vignette <- "Larson_2026_sulbactam_durlobactam_pediatric"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Durlobactam plasma residual variability is stratified by study phase in the
  # adult parent model, so the canonical propSd / addSd used by the error model
  # are derived inside model() from three phase-specific ini() magnitudes. The
  # poster changed only the body-size terms, so the residual structure is
  # carried over verbatim.
  paper_specific_residual_sds <- c(
    "propSdPhase1", "propSdPhase2", "propSdPhase3", "addSdPhase1"
  )

  compartmentData <- list(
    central         = list(analyte = "durlobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "durlobactam", units = "mg", specimen = "plasma", verified = TRUE),
    central_sbt     = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_sbt = list(analyte = "sulbactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The single covariate the pediatric adaptation changed. Larson 2026",
        "Methods / poster P-444 Methods: 'The first approach utilized the full",
        "covariate model from adults with the only change being the use of",
        "allometric scaling with fixed exponents of 0.75 for clearance and 1.0",
        "for volume of distribution.' The poster's Figure 1 notes add: 'WTKG",
        "effects added to Q and Vp for both drugs to account for likely changes",
        "in these parameters in children.' Weight therefore enters ALL FOUR",
        "disposition parameters of BOTH drugs - CL and Q at 0.75, Vc and Vp at",
        "1.0 - whereas the adult parent model carried estimated exponents on CL",
        "and Vc only (0.646 / 0.521 durlobactam, 1.01 / 0.831 sulbactam) and no",
        "weight effect at all on Q or Vp. Reference weight 75 kg is inherited",
        "from the adult model, whose reference subject is defined in Cammarata",
        "2024 Results as having 'a body weight of 75 kg'; the poster restates",
        "neither the reference weight nor the covariate equation, and says only",
        "that the exponents changed. Baseline (time-fixed) per subject.",
        "Pediatric applicability range: birth (28 weeks of gestation) to",
        "< 18 years of age."
      ),
      source_name        = "WTKG"
    ),
    CRCL = list(
      description        = paste(
        "Baseline creatinine clearance, normalized to body surface area"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried over unchanged from the adult parent model, where it drives a",
        "power function on the RENAL clearance arm only; the non-renal arm is",
        "CLcr-independent. Poster P-444 Figure 1 notes confirm the units for",
        "the pediatric simulations: 'CLCR represents GFR normalized to BSA in",
        "units of mL/min/1.73m2', and the poster Methods restate the adult",
        "structure: 'Both renal clearance and nonrenal clearance were",
        "estimated, and total clearance was calculated as the sum of renal and",
        "nonrenal clearance. Individual renal clearances were scaled by",
        "baseline creatinine clearance (CLCR).' Reference value 100",
        "mL/min/1.73 m^2, from Cammarata 2024 Results ('normal renal function",
        "(CLcr of 100 mL/min/1.73 m^2)'). The pediatric simulation cohort",
        "assumed normal renal function; the Larson 2026 abstract adds that",
        "renal function was estimated using the Rhodin formula from fat-free",
        "mass, but neither source on disk reports that formula's constants, so",
        "the covariate VALUES are a simulation input rather than part of this",
        "model - see the vignette Errata."
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
        "Retained unchanged from the adult parent model because the poster",
        "describes this approach as using 'the full covariate model from",
        "adults'. Applies a proportional shift to TOTAL CL on top of the",
        "continuous CRCL power function on the renal arm (Cammarata 2024",
        "Table 1 row BCLCRNLT30). The pediatric simulation cohort assumed",
        "normal renal function, so the term is inert there, but it is kept so",
        "the packaged model is the full adult covariate model the poster says",
        "it is. Derived from CRCL; supply as 1 * (CRCL < 30). Contrast with",
        "the sibling 'Allometry + CLCR' model, which drops it."
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
        "Time-varying WITHIN subject. Retained unchanged from the adult parent",
        "model (Cammarata 2024 Table S4) as part of 'the full covariate model",
        "from adults'. Encoded as a log-scale multiplicative factor on TOTAL",
        "CL that is switched on only while the gate is 1, so clearance reduces",
        "exactly to the plasma-model value between sessions and non-dialysed",
        "subjects are unaffected. The pediatric simulations of poster P-444",
        "used a cohort with normal renal function and no dialysis, so this",
        "column is 0 throughout the published pediatric analysis."
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
        "Retained unchanged from the adult parent model (proportional shifts",
        "on durlobactam CL and Vc; sulbactam carries no region effect). Poster",
        "P-444 Figure 1 notes state the value used for the pediatric",
        "simulations: 'East Asian flag was set to 'not East Asian'', i.e.",
        "REGION_EASTASIA = 0 for the whole simulated pediatric cohort."
      ),
      source_name        = "EASIAFL"
    ),
    DIS_HABP = list(
      description        = "Hospital-acquired bacterial pneumonia cohort indicator (1 = HABP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Retained unchanged from the adult parent model (Cammarata 2024",
        "Table 1 level INFTYPN1). The five infection-type indicators are",
        "mutually exclusive. Poster P-444 set infection type to 'Bacteremia'",
        "for the pediatric simulations, so this indicator is 0 there."
      ),
      source_name        = "INFTYPN = 1"
    ),
    DIS_VABP = list(
      description        = "Ventilator-associated bacterial pneumonia cohort indicator (1 = VABP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Retained unchanged from the adult parent model (Cammarata 2024",
        "Table 1 level INFTYPN2). Shares the merged durlobactam HABP-and-VABP",
        "Vc coefficient. 0 in the poster's pediatric simulations."
      ),
      source_name        = "INFTYPN = 2"
    ),
    DIS_CUTI = list(
      description        = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Retained unchanged from the adult parent model (Cammarata 2024",
        "Table 1 level INFTYPN3). 0 in the poster's pediatric simulations."
      ),
      source_name        = "INFTYPN = 3"
    ),
    DIS_BACTEREMIA = list(
      description        = "Bacteremia / bloodstream-infection cohort indicator (1 = bacteremia)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Retained unchanged from the adult parent model (Cammarata 2024",
        "Table 1 level INFTYPN4). This is the level the pediatric simulations",
        "used: poster P-444 Figure 1 notes state 'Infection type was set to",
        "'Bacteremia''. It carries the largest single covariate effect in the",
        "adult model (durlobactam Vc +3.32) together with a -0.444",
        "proportional shift on sulbactam CL, so this choice raises simulated",
        "exposures relative to an uninfected reference subject and is part of",
        "why this approach is the more conservative of the two."
      ),
      source_name        = "INFTYPN = 4"
    ),
    DIS_AP = list(
      description        = "Acute pyelonephritis cohort indicator (1 = AP)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled healthy-volunteer / non-infected Phase 1 cohort)",
      notes              = paste(
        "Retained unchanged from the adult parent model (Cammarata 2024",
        "Table 1 level INFTYPN5). Durlobactam carries no AP term. 0 in the",
        "poster's pediatric simulations."
      ),
      source_name        = "INFTYPN = 5"
    ),
    STUDY_SULDUR_PHASE2 = list(
      description        = "Phase 2 study cohort indicator (1 = Study CS2514-2017-0003)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 studies, when STUDY_SULDUR_PHASE3 is also 0)",
      notes              = paste(
        "Selects the durlobactam Phase 2 proportional residual magnitude,",
        "carried over unchanged from the adult parent model. The poster",
        "changed only the body-size terms and reports no residual-error",
        "re-estimation, so the adult phase-stratified residual structure is",
        "retained. Paired with STUDY_SULDUR_PHASE3; both 0 selects Phase 1."
      ),
      source_name        = "study phase"
    ),
    STUDY_SULDUR_PHASE3 = list(
      description        = "Phase 3 study cohort indicator (1 = Study CS2514-2017-0004)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 studies, when STUDY_SULDUR_PHASE2 is also 0)",
      notes              = paste(
        "Selects the durlobactam Phase 3 proportional residual magnitude,",
        "carried over unchanged from the adult parent model. See the",
        "STUDY_SULDUR_PHASE2 notes."
      ),
      source_name        = "study phase"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 8000L,
    n_studies        = 0L,
    age_range        = "birth (28 weeks of gestation) to < 18 years",
    weight_range     = "not reported; CDC growth-chart body size by age and sex",
    disease_state    = paste(
      "Simulated pediatric patients with Acinetobacter",
      "baumannii-calcoaceticus complex infection. This is a SIMULATION",
      "population, not an estimation data set: no pediatric clinical trial",
      "data existed when the analysis was done. The parameter values were",
      "estimated in 373 adults contributing 5,188 plasma concentrations",
      "pooled from six Phase 1, one Phase 2, and one Phase 3 study",
      "(Cammarata 2024); the pediatric adaptation changed only the body-size",
      "terms and re-used every estimate unchanged. The simulated cohort was",
      "assigned an infection type of bacteremia and a non-East-Asian region."
    ),
    dose_range       = paste(
      "Simulations supporting the Phase 1b regimens of poster P-444 Table 2:",
      "25 mg/kg sulbactam with 25 mg/kg durlobactam q6h (12 to < 18 y, 6 to",
      "< 12 y, 1 to < 6 y, 3 mo to < 1 y, and term infants 2 to < 3 mo, with",
      "a cap of 1 g of each drug in the two oldest cohorts), 25 mg/kg q8h for",
      "term infants from birth to < 2 mo, 20 mg/kg q8h for preterm infants 2",
      "to < 3 mo, and 20 mg/kg q12h for preterm infants from birth to < 2 mo.",
      "The sulbactam-durlobactam ratio is 1:1 and every regimen is a 3-hour",
      "intravenous infusion. Birth is defined as 7 days post-natal."
    ),
    renal_function   = paste(
      "Normal renal function was assumed throughout. Poster P-444 Methods:",
      "the simulated dataset incorporated 'age- and sex-specific body size",
      "distributions based on the published Centers for Disease Control and",
      "Prevention (CDC) growth chart for height and weight and assuming",
      "normal renal function'. The Larson 2026 abstract adds that renal",
      "function was estimated using the Rhodin formula from fat-free mass;",
      "that formula's constants are not reported in either source."
    ),
    notes            = paste(
      "8,000 hypothetical pediatric patients, 1,000 per cohort across the",
      "eight age / maturity cohorts of poster P-444 Table 2. PK/PD target",
      "attainment was assessed against the poster's stated drivers: >= 50%",
      "fT>MIC for sulbactam and a free-drug AUC/MIC ratio of 10 for",
      "durlobactam. Dosing regimens were selected to give PTA >= 90% with a",
      "median exposure similar to adults, judged against the adult Day 1 and",
      "Day 3 exposures of poster P-444 Table 1."
    )
  )

  ini({
    # =====================================================================
    # DURLOBACTAM structural parameters. Every value is inherited verbatim
    # from the adult parent model (Cammarata 2024 Table 1, 'Final model /
    # Estimate' column, Durlobactam block); the poster re-used the adult
    # estimates and changed only the body-weight exponents below. Values
    # refer to the reference subject: CRCL = 100 mL/min/1.73 m^2, WT = 75 kg,
    # non-East-Asian region, no infection, CLcr >= 30 mL/min/1.73 m^2, no
    # hemodialysis running.
    # =====================================================================
    lcl <- log(9.33); label("Durlobactam total CL at CRCL = 100 mL/min/1.73 m^2, WT = 75 kg (L/h)")   # Cammarata 2024 Table 1 Durlobactam: CL = 9.33 L/h (%SEM 3.24)
    lvc <- log(12.5); label("Durlobactam central volume Vc at WT = 75 kg, no infection (L)")          # Cammarata 2024 Table 1 Durlobactam: Vc = 12.5 L (%SEM 2.93)
    lq  <- log(4.43); label("Durlobactam inter-compartmental clearance Q at WT = 75 kg (L/h)")        # Cammarata 2024 Table 1 Durlobactam: Q = 4.43 L/h (%SEM 3.77)
    lvp <- log(5.83); label("Durlobactam peripheral volume Vp at WT = 75 kg (L)")                     # Cammarata 2024 Table 1 Durlobactam: Vp = 5.83 L (%SEM 3.44)

    f_renal <- fixed(0.479); label("Durlobactam fraction of total CL excreted renally (unitless)")
    # Cammarata 2024 Table 1 Durlobactam: FE (%) = 0.479, reported in the
    # shaded 'fixed or not estimated' style with no %SEM.

    e_crcl_cl_renal <- 0.875; label("Durlobactam CLcr power exponent on the renal CL arm (unitless)") # Cammarata 2024 Table 1 Durlobactam: 'CL R, CLcr power' = 0.875 (%SEM 12)

    # ---------------------------------------------------------------------
    # ALLOMETRIC BODY-SIZE EXPONENTS - the ONLY change the pediatric
    # adaptation made to the adult model.
    #
    # Poster P-444 Methods: 'The first approach utilized the full covariate
    # model from adults with the only change being the use of allometric
    # scaling with fixed exponents of 0.75 for clearance and 1.0 for volume
    # of distribution.' Poster P-444 Figure 1 notes: 'WTKG effects added to
    # Q and Vp for both drugs to account for likely changes in these
    # parameters in children.'
    #
    # Q is a clearance and takes 0.75; Vp is a volume and takes 1.0. The
    # adult parent model carried ESTIMATED exponents on CL and Vc only
    # (0.646 and 0.521 for durlobactam) and no weight effect on Q or Vp.
    # ---------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Durlobactam allometric exponent on total CL (unitless)") # Poster P-444 Methods: fixed exponent 0.75 for clearance
    e_wt_vc <- fixed(1.0);  label("Durlobactam allometric exponent on Vc (unitless)")       # Poster P-444 Methods: fixed exponent 1.0 for volume of distribution
    e_wt_q  <- fixed(0.75); label("Durlobactam allometric exponent on Q (unitless)")        # Poster P-444 Methods (0.75 for clearance) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'
    e_wt_vp <- fixed(1.0);  label("Durlobactam allometric exponent on Vp (unitless)")       # Poster P-444 Methods (1.0 for volume) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'

    e_renalimp_sev_cl     <- -0.58;  label("Durlobactam proportional shift in total CL for CLcr < 30 mL/min/1.73 m^2 (fraction)")  # Cammarata 2024 Table 1 Durlobactam: BCLCRNLT30 = -0.58 (%SEM 8.15)
    e_region_eastasia_cl  <- -0.199; label("Durlobactam proportional shift in total CL for East Asian region (fraction)")          # Cammarata 2024 Table 1 Durlobactam: CLEASIAFL1 = -0.199 (%SEM 31.7)
    e_region_eastasia_vc  <- -0.263; label("Durlobactam proportional shift in Vc for East Asian region (fraction)")                # Cammarata 2024 Table 1 Durlobactam: V1EASIAFL1 = -0.263 (%SEM 49.2)

    e_habp_vabp_vc   <- 1.52;  label("Durlobactam proportional shift in Vc for HABP or VABP (fraction)")  # Cammarata 2024 Table 1 Durlobactam: V1INFTYPN1&2 = 1.52 (%SEM 26.5)
    e_cuti_vc        <- 0.343; label("Durlobactam proportional shift in Vc for cUTI (fraction)")          # Cammarata 2024 Table 1 Durlobactam: V1INFTYPN3 = 0.343 (%SEM 62.3)
    e_bacteremia_vc  <- 3.32;  label("Durlobactam proportional shift in Vc for bacteremia (fraction)")    # Cammarata 2024 Table 1 Durlobactam: V1INFTYPN4 = 3.32 (%SEM 53.7)

    # =====================================================================
    # SULBACTAM structural parameters (Cammarata 2024 Table 1, Sulbactam
    # block), same reference subject.
    # =====================================================================
    lcl_sbt <- log(13.5); label("Sulbactam total CL at CRCL = 100 mL/min/1.73 m^2, WT = 75 kg (L/h)")  # Cammarata 2024 Table 1 Sulbactam: CL = 13.5 L/h (%SEM 14.0)
    lvc_sbt <- log(12.0); label("Sulbactam central volume Vc at WT = 75 kg, no infection (L)")         # Cammarata 2024 Table 1 Sulbactam: Vc = 12 L (%SEM 8.90)
    lq_sbt  <- log(7.88); label("Sulbactam inter-compartmental clearance Q at WT = 75 kg (L/h)")       # Cammarata 2024 Table 1 Sulbactam: Q = 7.88 L/h (%SEM 18.9)
    lvp_sbt <- log(6.99); label("Sulbactam peripheral volume Vp at WT = 75 kg (L)")                    # Cammarata 2024 Table 1 Sulbactam: Vp = 6.99 L (%SEM 9.29)

    f_renal_sbt <- fixed(0.479); label("Sulbactam fraction of total CL excreted renally (unitless)")
    # Cammarata 2024 Table 1 Sulbactam: FE (%) = 0.479, shaded 'fixed or not estimated'.

    e_crcl_cl_renal_sbt <- 1.14; label("Sulbactam CLcr power exponent on the renal CL arm (unitless)") # Cammarata 2024 Table 1 Sulbactam: 'CL R, CLcr power' = 1.14 (%SEM 20.8)

    e_wt_cl_sbt <- fixed(0.75); label("Sulbactam allometric exponent on total CL (unitless)") # Poster P-444 Methods: fixed exponent 0.75 for clearance
    e_wt_vc_sbt <- fixed(1.0);  label("Sulbactam allometric exponent on Vc (unitless)")       # Poster P-444 Methods: fixed exponent 1.0 for volume of distribution
    e_wt_q_sbt  <- fixed(0.75); label("Sulbactam allometric exponent on Q (unitless)")        # Poster P-444 Methods (0.75 for clearance) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'
    e_wt_vp_sbt <- fixed(1.0);  label("Sulbactam allometric exponent on Vp (unitless)")       # Poster P-444 Methods (1.0 for volume) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'

    e_renalimp_sev_cl_sbt <- -0.635; label("Sulbactam proportional shift in total CL for CLcr < 30 mL/min/1.73 m^2 (fraction)") # Cammarata 2024 Table 1 Sulbactam: BCLCRNLT30 = -0.635 (%SEM 11)

    e_habp_cl_sbt       <- -0.424; label("Sulbactam proportional shift in total CL for HABP (fraction)")       # Cammarata 2024 Table 1 Sulbactam: CLINFTYPN1 = -0.424 (%SEM 19.4)
    e_vabp_cl_sbt       <- -0.298; label("Sulbactam proportional shift in total CL for VABP (fraction)")       # Cammarata 2024 Table 1 Sulbactam: CLINFTYPN2 = -0.298 (%SEM 36.5)
    e_cuti_cl_sbt       <- -0.157; label("Sulbactam proportional shift in total CL for cUTI (fraction)")       # Cammarata 2024 Table 1 Sulbactam: CLINFTYPN3 = -0.157 (%SEM 96.8)
    e_bacteremia_cl_sbt <- -0.444; label("Sulbactam proportional shift in total CL for bacteremia (fraction)") # Cammarata 2024 Table 1 Sulbactam: CLINFTYPN4 = -0.444 (%SEM 20.9)
    e_ap_cl_sbt         <- -0.382; label("Sulbactam proportional shift in total CL for acute pyelonephritis (fraction)") # Cammarata 2024 Table 1 Sulbactam: CLINFTYPN5 = -0.382 (%SEM 28.8)

    e_habp_vc_sbt       <-  0.836; label("Sulbactam proportional shift in Vc for HABP (fraction)")       # Cammarata 2024 Table 1 Sulbactam: V3INFTYPN1 = 0.836 (%SEM 25.3)
    e_vabp_vc_sbt       <-  1.43;  label("Sulbactam proportional shift in Vc for VABP (fraction)")       # Cammarata 2024 Table 1 Sulbactam: V3INFTYPN2 = 1.43 (%SEM 19.5)
    e_cuti_vc_sbt       <-  0.17;  label("Sulbactam proportional shift in Vc for cUTI (fraction)")       # Cammarata 2024 Table 1 Sulbactam: V3INFTYPN3 = 0.17 (%SEM 89.8)
    e_bacteremia_vc_sbt <-  1.85;  label("Sulbactam proportional shift in Vc for bacteremia (fraction)") # Cammarata 2024 Table 1 Sulbactam: V3INFTYPN4 = 1.85 (%SEM 44.7)
    e_ap_vc_sbt         <- -0.704; label("Sulbactam proportional shift in Vc for acute pyelonephritis (fraction)") # Cammarata 2024 Table 1 Sulbactam: V3INFTYPN5 = -0.704 (%SEM 13.2)

    # =====================================================================
    # HEMODIALYSIS SUB-MODEL (Cammarata 2024 Table S4), retained as part of
    # 'the full covariate model from adults'. Inert whenever
    # RRT_HEMODIAL_ACTIVE = 0, which is the case throughout the poster's
    # pediatric simulations.
    # =====================================================================
    e_hemodial_active_cl     <- log(6.24); label("Log fold-increase in durlobactam total CL while hemodialysis is running (unitless)") # Cammarata 2024 Table S4: CL-HDEFFECT (durlobactam) = 6.24 (%SEM 14.5)
    e_hemodial_active_cl_sbt <- log(8.19); label("Log fold-increase in sulbactam total CL while hemodialysis is running (unitless)")   # Cammarata 2024 Table S4: CL-HDEFFECT (sulbactam) = 8.19 (%SEM 23.2)

    # =====================================================================
    # EPITHELIAL LINING FLUID SUB-MODEL (Cammarata 2024 Table S6). Retained:
    # poster P-444 Figures 2[C] and 3[C] present 'total-drug ELF' exposures
    # for the pediatric cohort, so the ELF observables are part of the
    # published pediatric analysis. ELF concentration is an instantaneous
    # ratio of the plasma concentration.
    # =====================================================================
    lrelf     <- log(0.372); label("Durlobactam ELF / total-drug plasma concentration ratio (unitless)") # Cammarata 2024 Table S6: PLASMA-ELF ratio (durlobactam) = 0.372 (%SEM 3.6)
    lrelf_sbt <- log(0.533); label("Sulbactam ELF / total-drug plasma concentration ratio (unitless)")   # Cammarata 2024 Table S6: PLASMA-ELF ratio (sulbactam) = 0.533 (%SEM 5.41)

    # =====================================================================
    # Inter-individual variability, inherited verbatim from the adult parent
    # model. The paper reports omega^2 directly (the parenthetical %CV is
    # sqrt(omega^2) x 100), so the variances below are the published values.
    # =====================================================================
    etalcl + etalvc ~ c(0.0778,
                        0.0494, 0.0757)
    # Cammarata 2024 Table 1 Durlobactam: omega^2 CL = 0.0778 (27.9 %CV);
    # omega^2 Vc = 0.0757 (27.5 %CV); Covariance = 0.0494 (r^2 = 0.415)
    etalvp ~ 0.0773
    # Cammarata 2024 Table 1 Durlobactam: omega^2 Vp = 0.0773 (27.8 %CV)

    etalcl_sbt + etalvc_sbt ~ c(0.221,
                                0.0727, 0.0967)
    # Cammarata 2024 Table 1 Sulbactam: omega^2 CL = 0.221 (47.0 %CV);
    # omega^2 Vc = 0.0967 (31.1 %CV); Covariance = 0.0727 (r^2 = 0.2448)
    etalvp_sbt ~ 0.196
    # Cammarata 2024 Table 1 Sulbactam: omega^2 Vp = 0.196 (44.3 %CV)

    etae_hemodial_active_cl     ~ 0.124
    # Cammarata 2024 Table S4: omega^2 CL-HDEFFECT (durlobactam) = 0.124 (35.2 %CV)
    etae_hemodial_active_cl_sbt ~ 0.316
    # Cammarata 2024 Table S4: omega^2 CL-HDEFFECT (sulbactam) = 0.316 (56.2 %CV)

    # =====================================================================
    # Residual variability, inherited verbatim from the adult parent model.
    # Both tables report sigma^2; the SDs below are sqrt(sigma^2).
    # =====================================================================
    propSdPhase1 <- sqrt(0.019);   label("Durlobactam proportional residual SD, Phase 1 plasma (fraction)")   # Cammarata 2024 Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 1 = 0.019 (13.8 %CV)
    addSdPhase1  <- sqrt(0.00136); label("Durlobactam additive residual SD, Phase 1 plasma (mg/L)")           # Cammarata 2024 Table 1 Durlobactam: sigma^2 plasma, Additive Phase 1 = 0.00136 (0.0369 mg/L)
    propSdPhase2 <- sqrt(0.0794);  label("Durlobactam proportional residual SD, Phase 2 plasma (fraction)")   # Cammarata 2024 Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 2 = 0.0794 (28.2 %CV)
    propSdPhase3 <- sqrt(0.203);   label("Durlobactam proportional residual SD, Phase 3 plasma (fraction)")   # Cammarata 2024 Table 1 Durlobactam: sigma^2 plasma, Proportional Phase 3 = 0.203 (45.0 %CV)

    propSd_sbt <- sqrt(0.0433); label("Sulbactam proportional residual SD, plasma (fraction)") # Cammarata 2024 Table 1 Sulbactam: sigma^2 plasma, Proportional = 0.0433 (20.8 %CV)
    addSd_sbt  <- sqrt(0.0054); label("Sulbactam additive residual SD, plasma (mg/L)")         # Cammarata 2024 Table 1 Sulbactam: sigma^2 plasma, Additive = 0.0054 (0.0735 mg/L)

    propSd_Celf     <- sqrt(0.0322); label("Durlobactam proportional residual SD, ELF (fraction)") # Cammarata 2024 Table S6: sigma^2 ELF, Proportional (durlobactam) = 0.0322 (17.9 %CV)
    propSd_Celf_sbt <- sqrt(0.0628); label("Sulbactam proportional residual SD, ELF (fraction)")   # Cammarata 2024 Table S6: sigma^2 ELF, Proportional (sulbactam) = 0.0628 (25.1 %CV)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariate terms. The five infection-type indicators are
    #    mutually exclusive, so the proportional shifts are summed inside a
    #    single (1 + ...) bracket; the shared reference category is the
    #    uninfected Phase 1 subject.
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

    # Hemodialysis fold-change on total CL; collapses to exp(0) = 1 between
    # sessions and in every subject who is not on dialysis.
    hdeffect     <- exp((e_hemodial_active_cl     + etae_hemodial_active_cl)     * RRT_HEMODIAL_ACTIVE)
    hdeffect_sbt <- exp((e_hemodial_active_cl_sbt + etae_hemodial_active_cl_sbt) * RRT_HEMODIAL_ACTIVE)

    # ------------------------------------------------------------------
    # 2. Durlobactam individual PK parameters. Identical to the adult parent
    #    model except that the body-weight power terms now use the fixed
    #    allometric exponents and are additionally applied to Q and Vp.
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

    q  <- exp(lq)          * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp

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

    q_sbt  <- exp(lq_sbt)              * (WT / 75)^e_wt_q_sbt
    vp_sbt <- exp(lvp_sbt + etalvp_sbt) * (WT / 75)^e_wt_vp_sbt

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
