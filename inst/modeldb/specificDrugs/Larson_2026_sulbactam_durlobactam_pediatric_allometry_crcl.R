Larson_2026_sulbactam_durlobactam_pediatric_allometry_crcl <- function() {
  description <- paste(
    "Pediatric 'Allometry + CLCR' adaptation of the adult",
    "sulbactam-durlobactam population PK model, the second of the two",
    "approaches explored to select dosing regimens for the Phase 1b pediatric",
    "trial NCT06801223 (Larson 2026 / Cammarata IDWeek 2025 poster P-444).",
    "Starting from the sibling 'Adult with Allometry' model, every covariate",
    "relationship is removed except body weight on all four disposition",
    "parameters of both drugs (fixed allometric exponents 0.75 on CL and Q,",
    "1.0 on Vc and Vp) and renal function on the renal clearance arm only.",
    "The infection-type and East Asian region proportional shifts, the",
    "severe-renal-impairment shift on total clearance, and the hemodialysis",
    "gate are all dropped; the four-compartment structure, the renal /",
    "non-renal clearance split, the inter-individual variability, the",
    "residual-error model, and the epithelial lining fluid ratios are",
    "retained. This approach predicted lower exposures than the sibling model",
    "and would have supported higher doses, so the sibling was chosen as the",
    "more conservative basis for initial dose selection. Durlobactam uses the",
    "unsuffixed canonical compartment / parameter set; sulbactam carries the",
    "sibling-drug suffix _sbt throughout. See",
    "modellib('Larson_2026_sulbactam_durlobactam_pediatric_allometry') and",
    "modellib('Cammarata_2024_sulbactam_durlobactam')."
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
    "IDWeek 2025, Atlanta, GA. The defining sentence for this model is the",
    "poster Methods bullet: 'The second approach removed all covariate",
    "relationships from the Adult Allometry model except for body weight (on",
    "all parameters) and renal function (on renal clearance only)'.",
    "All retained structural, IIV and residual-error values are inherited",
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
  # adult parent model. The poster's second approach removed covariate
  # relationships on the PK PARAMETERS; it reports no re-estimation of the
  # residual-error model, so the adult phase-stratified residual structure is
  # carried over unchanged, exactly as in the sibling model.
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
        "One of the only two covariates this model retains. Poster P-444",
        "Methods: 'The second approach removed all covariate relationships",
        "from the Adult Allometry model except for body weight (on all",
        "parameters) and renal function (on renal clearance only).' Weight",
        "therefore enters ALL FOUR disposition parameters of BOTH drugs with",
        "the fixed allometric exponents of the Adult Allometry model - CL and",
        "Q at 0.75, Vc and Vp at 1.0 - per the poster Methods ('allometric",
        "scaling with fixed exponents of 0.75 for clearance and 1.0 for volume",
        "of distribution') and the Figure 1 note ('WTKG effects added to Q and",
        "Vp for both drugs to account for likely changes in these parameters",
        "in children'). Reference weight 75 kg is inherited from the adult",
        "parent model, whose reference subject is defined in Cammarata 2024",
        "Results as having 'a body weight of 75 kg'; neither the poster nor",
        "the abstract restates it. Baseline (time-fixed) per subject.",
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
        "The other retained covariate, and the one this model is named for.",
        "Poster P-444 Methods keeps 'renal function (on renal clearance",
        "only)', so the CLcr power function on the RENAL clearance arm is",
        "retained with the adult exponents while the severe-renal-impairment",
        "proportional shift on TOTAL clearance - which is a renal-function",
        "effect but NOT on renal clearance - is dropped. Poster P-444",
        "Figure 1 notes give the units: 'CLCR represents GFR normalized to BSA",
        "in units of mL/min/1.73m2'. Reference value 100 mL/min/1.73 m^2, from",
        "Cammarata 2024 Results ('normal renal function (CLcr of 100",
        "mL/min/1.73 m^2)'). The Larson 2026 abstract states that renal",
        "function for the simulated pediatric cohort was estimated using the",
        "Rhodin formula from fat-free mass; that formula's constants are not",
        "reported in either source on disk, so the covariate VALUES are a",
        "simulation input rather than part of this model - see the vignette",
        "Errata."
      ),
      source_name        = "CLcr"
    ),
    STUDY_SULDUR_PHASE2 = list(
      description        = "Phase 2 study cohort indicator (1 = Study CS2514-2017-0003)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 studies, when STUDY_SULDUR_PHASE3 is also 0)",
      notes              = paste(
        "Selects the durlobactam Phase 2 proportional residual magnitude,",
        "carried over unchanged from the adult parent model. The poster's",
        "second approach removed covariate relationships on the PK",
        "PARAMETERS; it reports no change to the residual-error model, so the",
        "adult phase-stratified residual structure is retained. Paired with",
        "STUDY_SULDUR_PHASE3; both 0 selects Phase 1."
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

  covariatesDataExcluded <- list(
    REGION_EASTASIA = list(
      description = "East Asian region of origin indicator (1 = China, Taiwan, or South Korea)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Present in the adult parent model and in the sibling Adult with",
        "Allometry model, but deliberately removed here: poster P-444 Methods",
        "state that this approach 'removed all covariate relationships from",
        "the Adult Allometry model except for body weight ... and renal",
        "function'. Documented so the provenance of the removal is visible."
      )
    ),
    RENALIMP_SEV = list(
      description = "Severe renal impairment indicator (1 = baseline CLcr < 30 mL/min/1.73 m^2)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Removed here. It is a renal-function covariate, but it acts as a",
        "proportional shift on TOTAL clearance, whereas poster P-444 retains",
        "renal function 'on renal clearance only'. The continuous CRCL power",
        "function on the renal arm is what survives."
      )
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemodialysis-session gate (1 while a session is running)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Removed here along with the other non-retained covariate",
        "relationships. The simulated pediatric cohort had normal renal",
        "function and no dialysis, so the term would have been inert in any",
        "case."
      )
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia cohort indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Infection-type shift removed per poster P-444 Methods, second approach."
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia cohort indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Infection-type shift removed per poster P-444 Methods, second approach."
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection cohort indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Infection-type shift removed per poster P-444 Methods, second approach."
    ),
    DIS_BACTEREMIA = list(
      description = "Bacteremia / bloodstream-infection cohort indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Infection-type shift removed per poster P-444 Methods, second",
        "approach. Note that the sibling Adult with Allometry model was",
        "simulated with infection type set to bacteremia, whose large",
        "durlobactam Vc shift (+3.32) and sulbactam CL shift (-0.444) are",
        "absent here - the main reason this model predicts lower exposures."
      )
    ),
    DIS_AP = list(
      description = "Acute pyelonephritis cohort indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Infection-type shift removed per poster P-444 Methods, second approach."
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
      "data existed when the analysis was done. The retained parameter values",
      "were estimated in 373 adults contributing 5,188 plasma concentrations",
      "pooled from six Phase 1, one Phase 2, and one Phase 3 study",
      "(Cammarata 2024). Because every infection-type covariate is removed,",
      "this model predicts the exposure of a body-size- and",
      "renal-function-matched subject with no infection-related shift."
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
      "eight age / maturity cohorts of poster P-444 Table 2. Poster P-444",
      "Results: 'Compared to the adult with allometry model, the adult with",
      "allometry + CLCR model predicted lower SUD exposures, which would have",
      "supported the selection of higher doses. Thus, the adult with allometry",
      "model represents a more conservative approach for initial dose",
      "selection for the pediatric clinical trial.' PK/PD target attainment",
      "was assessed against >= 50% fT>MIC for sulbactam and a free-drug",
      "AUC/MIC ratio of 10 for durlobactam."
    )
  )

  ini({
    # =====================================================================
    # DURLOBACTAM structural parameters, inherited verbatim from the adult
    # parent model (Cammarata 2024 Table 1, Durlobactam block). Values refer
    # to the reference subject: CRCL = 100 mL/min/1.73 m^2, WT = 75 kg.
    # With every infection-type, region and severe-renal-impairment shift
    # removed, those are the only two covariate conditions that remain.
    # =====================================================================
    lcl <- log(9.33); label("Durlobactam total CL at CRCL = 100 mL/min/1.73 m^2, WT = 75 kg (L/h)") # Cammarata 2024 Table 1 Durlobactam: CL = 9.33 L/h (%SEM 3.24)
    lvc <- log(12.5); label("Durlobactam central volume Vc at WT = 75 kg (L)")                      # Cammarata 2024 Table 1 Durlobactam: Vc = 12.5 L (%SEM 2.93)
    lq  <- log(4.43); label("Durlobactam inter-compartmental clearance Q at WT = 75 kg (L/h)")      # Cammarata 2024 Table 1 Durlobactam: Q = 4.43 L/h (%SEM 3.77)
    lvp <- log(5.83); label("Durlobactam peripheral volume Vp at WT = 75 kg (L)")                   # Cammarata 2024 Table 1 Durlobactam: Vp = 5.83 L (%SEM 3.44)

    f_renal <- fixed(0.479); label("Durlobactam fraction of total CL excreted renally (unitless)")
    # Cammarata 2024 Table 1 Durlobactam: FE (%) = 0.479, shaded 'fixed or not estimated'.

    e_crcl_cl_renal <- 0.875; label("Durlobactam CLcr power exponent on the renal CL arm (unitless)") # Cammarata 2024 Table 1 Durlobactam: 'CL R, CLcr power' = 0.875 (%SEM 12); RETAINED per poster P-444 'renal function (on renal clearance only)'

    e_wt_cl <- fixed(0.75); label("Durlobactam allometric exponent on total CL (unitless)") # Poster P-444 Methods: fixed exponent 0.75 for clearance
    e_wt_vc <- fixed(1.0);  label("Durlobactam allometric exponent on Vc (unitless)")       # Poster P-444 Methods: fixed exponent 1.0 for volume of distribution
    e_wt_q  <- fixed(0.75); label("Durlobactam allometric exponent on Q (unitless)")        # Poster P-444 Methods (0.75 for clearance) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'; poster Methods 'body weight (on all parameters)'
    e_wt_vp <- fixed(1.0);  label("Durlobactam allometric exponent on Vp (unitless)")       # Poster P-444 Methods (1.0 for volume) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'; poster Methods 'body weight (on all parameters)'

    # =====================================================================
    # SULBACTAM structural parameters (Cammarata 2024 Table 1, Sulbactam
    # block), same reference subject.
    # =====================================================================
    lcl_sbt <- log(13.5); label("Sulbactam total CL at CRCL = 100 mL/min/1.73 m^2, WT = 75 kg (L/h)") # Cammarata 2024 Table 1 Sulbactam: CL = 13.5 L/h (%SEM 14.0)
    lvc_sbt <- log(12.0); label("Sulbactam central volume Vc at WT = 75 kg (L)")                      # Cammarata 2024 Table 1 Sulbactam: Vc = 12 L (%SEM 8.90)
    lq_sbt  <- log(7.88); label("Sulbactam inter-compartmental clearance Q at WT = 75 kg (L/h)")      # Cammarata 2024 Table 1 Sulbactam: Q = 7.88 L/h (%SEM 18.9)
    lvp_sbt <- log(6.99); label("Sulbactam peripheral volume Vp at WT = 75 kg (L)")                   # Cammarata 2024 Table 1 Sulbactam: Vp = 6.99 L (%SEM 9.29)

    f_renal_sbt <- fixed(0.479); label("Sulbactam fraction of total CL excreted renally (unitless)")
    # Cammarata 2024 Table 1 Sulbactam: FE (%) = 0.479, shaded 'fixed or not estimated'.

    e_crcl_cl_renal_sbt <- 1.14; label("Sulbactam CLcr power exponent on the renal CL arm (unitless)") # Cammarata 2024 Table 1 Sulbactam: 'CL R, CLcr power' = 1.14 (%SEM 20.8); RETAINED per poster P-444 'renal function (on renal clearance only)'

    e_wt_cl_sbt <- fixed(0.75); label("Sulbactam allometric exponent on total CL (unitless)") # Poster P-444 Methods: fixed exponent 0.75 for clearance
    e_wt_vc_sbt <- fixed(1.0);  label("Sulbactam allometric exponent on Vc (unitless)")       # Poster P-444 Methods: fixed exponent 1.0 for volume of distribution
    e_wt_q_sbt  <- fixed(0.75); label("Sulbactam allometric exponent on Q (unitless)")        # Poster P-444 Methods (0.75 for clearance) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'
    e_wt_vp_sbt <- fixed(1.0);  label("Sulbactam allometric exponent on Vp (unitless)")       # Poster P-444 Methods (1.0 for volume) with the Figure 1 note 'WTKG effects added to Q and Vp for both drugs'

    # =====================================================================
    # EPITHELIAL LINING FLUID SUB-MODEL (Cammarata 2024 Table S6). Retained:
    # poster P-444 Figure 3[C] presents 'total-drug ELF' exposures for this
    # model as well. An ELF-to-plasma ratio is not a covariate relationship,
    # so it is unaffected by the covariate removal.
    # =====================================================================
    lrelf     <- log(0.372); label("Durlobactam ELF / total-drug plasma concentration ratio (unitless)") # Cammarata 2024 Table S6: PLASMA-ELF ratio (durlobactam) = 0.372 (%SEM 3.6)
    lrelf_sbt <- log(0.533); label("Sulbactam ELF / total-drug plasma concentration ratio (unitless)")   # Cammarata 2024 Table S6: PLASMA-ELF ratio (sulbactam) = 0.533 (%SEM 5.41)

    # =====================================================================
    # Inter-individual variability, inherited verbatim from the adult parent
    # model. The paper reports omega^2 directly (the parenthetical %CV is
    # sqrt(omega^2) x 100). The hemodialysis etas are absent because the
    # hemodialysis term is not part of this model.
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

    # =====================================================================
    # Residual variability, inherited verbatim from the adult parent model.
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
    # 1. Durlobactam individual PK parameters. Only body weight and the
    #    renal-arm CLcr power function survive the covariate removal:
    #
    #      CL_R  = CL * FE       * (CLcr / 100)^theta_CLcr
    #      CL_NR = CL * (1 - FE)
    #      CL_T  = (CL_R + CL_NR) * (WT / 75)^0.75 * exp(eta_CL)
    # ------------------------------------------------------------------
    cl_renal  <- exp(lcl) * f_renal * (CRCL / 100)^e_crcl_cl_renal
    cl_nonren <- exp(lcl) * (1 - f_renal)

    cl <- (cl_renal + cl_nonren) * (WT / 75)^e_wt_cl * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    q  <- exp(lq)           * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp

    # ------------------------------------------------------------------
    # 2. Sulbactam individual PK parameters (same construction).
    # ------------------------------------------------------------------
    cl_renal_sbt  <- exp(lcl_sbt) * f_renal_sbt * (CRCL / 100)^e_crcl_cl_renal_sbt
    cl_nonren_sbt <- exp(lcl_sbt) * (1 - f_renal_sbt)

    cl_sbt <- (cl_renal_sbt + cl_nonren_sbt) * (WT / 75)^e_wt_cl_sbt * exp(etalcl_sbt)
    vc_sbt <- exp(lvc_sbt + etalvc_sbt) * (WT / 75)^e_wt_vc_sbt
    q_sbt  <- exp(lq_sbt)               * (WT / 75)^e_wt_q_sbt
    vp_sbt <- exp(lvp_sbt + etalvp_sbt) * (WT / 75)^e_wt_vp_sbt

    # ------------------------------------------------------------------
    # 3. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_sbt <- cl_sbt / vc_sbt
    k12_sbt <- q_sbt  / vc_sbt
    k21_sbt <- q_sbt  / vp_sbt

    # ------------------------------------------------------------------
    # 4. Four-compartment system: two independent two-compartment IV
    #    dispositions. Doses are 3-hour intravenous infusions into the
    #    respective central compartments.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    d/dt(central_sbt)     <- -(kel_sbt + k12_sbt) * central_sbt + k21_sbt * peripheral1_sbt
    d/dt(peripheral1_sbt) <-   k12_sbt * central_sbt            - k21_sbt * peripheral1_sbt

    # ------------------------------------------------------------------
    # 5. Observations. Dose in mg, volumes in L -> concentrations in mg/L.
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
