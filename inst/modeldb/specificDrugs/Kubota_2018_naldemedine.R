Kubota_2018_naldemedine <- function() {
  description <- paste0(
    "Population pharmacokinetic model of naldemedine, a peripherally acting ",
    "mu-opioid receptor antagonist for opioid-induced constipation (OIC), ",
    "pooled across 18 phase 1 to phase 3 studies in healthy subjects, ",
    "patients with chronic non-cancer pain and OIC, and cancer patients with ",
    "OIC (Kubota 2018, n = 949, 8,146 plasma concentrations). Two-compartment ",
    "disposition with first-order absorption and an absorption lag time. ",
    "Apparent oral clearance carries median-centred power effects of age and ",
    "creatinine clearance plus multiplicative non-White-race and female-sex ",
    "factors; apparent central volume carries a linear body-weight ratio plus ",
    "multiplicative health-status (chronic non-cancer pain OIC, cancer OIC) ",
    "and fed-state factors; the absorption rate constant carries a ",
    "median-centred power effect of age. The reference subject is a ",
    "52-year-old, 76 kg, White, male patient with chronic non-cancer pain and ",
    "OIC, creatinine clearance 108 mL/min, dosed fasted: CL/F 9.10 L/h, ",
    "Vc/F 91.1 L, Ka 2.94 1/h, Q/F 4.77 L/h, Vp/F 41.8 L, ALAG 0.195 h. ",
    "Seven companion landmark exposure-response models in the ",
    "Kubota_2018_naldemedine_* family relate the steady-state AUC this model ",
    "produces to the probability of a spontaneous-bowel-movement response ",
    "and to the probability of gastrointestinal adverse events."
  )
  reference <- paste(
    "Kubota R, Fukumura K, Wajima T.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Naldemedine.",
    "Pharm Res. 2018;35(11):225.",
    "doi:10.1007/s11095-018-2501-7. PMCID: PMC6182381.",
    "Final-model parameter estimates and the covariate equations are Table III;",
    "the covariate functional forms are Eqs. (1) and (2);",
    "baseline covariate distributions are Table II.",
    sep = " "
  )
  vignette <- "Kubota_2018_naldemedine"
  units <- list(
    time = "h",
    dosing = "mg (naldemedine free base; the clinical dose is 0.2 mg once daily)",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amount units follow the `units$dosing` entry above;
  # the ng/mL observation scale is produced by the 1000x factor in model().
  compartmentData <- list(
    depot = list(
      analyte = "naldemedine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "naldemedine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "naldemedine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters BOTH CL/F and Ka as a median-centred power term, (AGE / 52)^theta,",
        "per Kubota 2018 Eq. (1) 'PKP = theta1 x (COV / median of COV)^theta2'.",
        "The centring constant 52 years is the cohort median from Table II",
        "(mean 51, SD 14, range 18-90). The two exponents have the same sign but",
        "very different magnitudes: -0.195 on CL/F (a 90-year-old has CL/F about",
        "10% below the 52-year-old reference) and -1.16 on Ka (the same subject",
        "absorbs roughly 1.9-fold more slowly). Age and creatinine clearance are",
        "correlated in this cohort (r = -0.45, Supplemental Table S6a), so the two",
        "clearance terms are not independently identified; see the vignette."
      ),
      source_name = "Age"
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline, calculated with the Cockcroft-Gault equation. RAW mL/min, NOT body-surface-area normalised.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as the median-centred power term (CRCL / 108)^0.0739 per",
        "Kubota 2018 Eq. (1). The centring constant 108 mL/min is the cohort median",
        "from Table II (mean 109.5, SD 42.4, range 5.8-311.8).",
        "IMPORTANT units caveat for reuse: this column is a raw Cockcroft-Gault",
        "creatinine clearance in mL/min, whereas the CRCL canonical's primary",
        "definition is BSA-normalised mL/min/1.73 m^2. Supplying a BSA-normalised",
        "value to this model biases CL/F by (ratio)^0.0739; the exponent is small",
        "enough that the error is modest, but it is still an error. The same raw",
        "Cockcroft-Gault convention is carried by Delattre_2010_amikacin.R and",
        "Jin_2026_colistinSulfate.R.",
        "The 95% CI of this exponent includes zero (-0.0133 to 0.161). Kubota 2018",
        "retained the term anyway (Discussion) because a renal-impairment study had",
        "shown a real AUC difference; that is a source-stated modelling decision,",
        "not an oversight, and it is reproduced here unchanged."
      ),
      source_name = "CLcr"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "486 of 949 subjects (51.2%) male, 463 (48.8%) female (Table II).",
        "Kubota 2018's source indicator is 'Gender = 1 for female, Gender = 0 for",
        "male' (Table II footnote), which is exactly the SEXF canonical polarity,",
        "so no value transformation is needed. Enters CL/F as the multiplicative",
        "power-of-indicator factor 0.902^SEXF per Eq. (2): female CL/F is 9.8%",
        "below male."
      ),
      source_name = "Gender"
    ),
    RACE_WHITE = list(
      description = "White race indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-White) for the canonical column; note that the SOURCE's reference level is White -- see notes.",
      notes = paste(
        "558 of 949 subjects (58.8%) White, 391 (41.2%) non-White (Table II).",
        "POLARITY INVERSION: Kubota 2018's indicator is named 'White' but is coded",
        "the other way round -- the Table II footnote states 'White = 1 for",
        "non-White, White = 0 for White'. The published factor 0.870 therefore",
        "applies to NON-White subjects, whose CL/F is 13% below the White",
        "reference. The RACE_WHITE canonical has the opposite polarity, so model()",
        "raises the published factor to the power (1 - RACE_WHITE); the estimate",
        "itself is carried unchanged and the inversion is visible at the call site.",
        "The Results text confirms the footnote independently: the reference",
        "population -- the subject at which every indicator is zero, so CL/F is",
        "THETA(1) = 9.10 L/h unmodified -- is described as '52-year-old, 76 kg",
        "male, white, non-cancer, OIC patients'. That subject is White and",
        "carries no race factor, so White = 0 for White.",
        "Reading the footnote the other way would move the effect onto the wrong",
        "group and invert the direction of a 13% clearance difference.",
        "Kubota 2018 tested White/non-White, Japanese/non-Japanese and",
        "Hispanic-or-Latino/non-Hispanic-or-non-Latino in SEPARATE models to avoid",
        "confounding; only White/non-White was retained."
      ),
      source_name = "White"
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Vc/F as the median-centred ratio (WT / 76). The centring constant",
        "76 kg is the cohort median from Table II (mean 79.5, SD 23.4, range",
        "34.4-188.1). Note that Table III writes this term as '* (Body weight/76)'",
        "with NO exponent, unlike the age and creatinine-clearance terms which are",
        "written as powers -- i.e. the exponent was held at 1 rather than",
        "estimated, so ini() carries it as fixed(1)."
      ),
      source_name = "Body weight"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant cohort indicator: 1 = healthy subject, 0 = patient (chronic non-cancer pain with OIC, or cancer with OIC).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient; the complement is the union of the two OIC patient cohorts)",
      notes = paste(
        "Kubota 2018 codes health status with TWO indicator variables over three",
        "levels (Table II footnote): 'non-Cancer = 1 and Cancer = 0 for patients",
        "with chronic non-cancer pain and OIC, non-Cancer = 0 and Cancer = 1 for",
        "cancer patients with OIC, non-Cancer = 0 and Cancer = 0 for healthy",
        "subjects'. Healthy subjects are therefore the level at which BOTH source",
        "indicators vanish.",
        "This model carries the three levels on the two ratified canonicals",
        "DIS_HEALTHY and DIS_CANCER rather than minting a new one, and recovers",
        "the source's chronic-non-cancer-pain indicator inside model() as",
        "dis_cncp_oic = (1 - DIS_HEALTHY) * (1 - DIS_CANCER). Encode a subject as",
        "healthy with DIS_HEALTHY = 1, DIS_CANCER = 0; as chronic-non-cancer-pain",
        "OIC with both 0; as cancer OIC with DIS_HEALTHY = 0, DIS_CANCER = 1.",
        "DIS_HEALTHY = 1 together with DIS_CANCER = 1 is not a state the source",
        "data can take and is not defined."
      ),
      source_name = "Health status (healthy subjects level)"
    ),
    DIS_CANCER = list(
      description = "Cancer-cohort indicator: 1 = cancer patient with opioid-induced constipation, 0 = healthy subject or patient with chronic non-cancer pain and OIC.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject or chronic-non-cancer-pain OIC patient)",
      notes = paste(
        "The oncology arm of the pooled analysis: studies 1108V9222 (phase 2b dose",
        "finding in cancer patients) and 1331V9236 (phase 3 in Japanese cancer",
        "patients). Enters Vc/F as the multiplicative power-of-indicator factor",
        "1.27^DIS_CANCER per Eq. (2): cancer patients have Vc/F 27% above the",
        "healthy reference.",
        "Note that Kubota 2018's cohort is 'cancer patients with OIC' generally,",
        "not specifically the advanced/metastatic solid-tumour population named in",
        "the DIS_CANCER register entry's primary description; this is within the",
        "entry's paper-defined-complement semantics, which is why its Scope is",
        "'specific'. See DIS_HEALTHY for the full three-level coding scheme."
      ),
      source_name = "Health status (cancer level)"
    ),
    FED = list(
      description = "Fed-versus-fasted state at the dose record: 1 = fed, 0 = fasted.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "817 of 949 subjects (86.1%) fasted, 132 (13.9%) fed (Table II).",
        "The Table II footnote states 'Food = 1 for fed condition, Food = 0 for",
        "fasted condition', which matches the FED canonical polarity directly.",
        "Enters Vc/F as the multiplicative factor 1.12^FED per Eq. (2).",
        "Two source conventions are worth carrying forward when assembling data",
        "for this model (Methods): food condition was ASSUMED fasted in the two",
        "phase 2b studies (1107V9221, 1108V9222), where food-intake times were not",
        "recorded; and in the phase 3 studies it was DEFINED as fasted when",
        "naldemedine was given more than 1 h after food."
      ),
      source_name = "Food condition"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  # Documented here for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index at baseline.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested on CL/F and Vc/F; not retained. Table II: mean 28.1, SD 7.2, median 26.9, range 14.4-58.8."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units = "g/dL",
      type = "continuous",
      notes = "Tested on CL/F only; not retained. Table II: mean 4.3, SD 0.5, median 4.3, range 2.1-5.4."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline.",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F only; not retained. Table II: mean 22, SD 13, median 19, range 6-223."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline.",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F only; not retained. Table II: mean 22, SD 15, median 18, range 2-212."
    ),
    TBILI = list(
      description = "Total bilirubin at baseline.",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested on CL/F only; not retained. Table II: mean 0.5, SD 0.3, median 0.4, range 0.04-2.4."
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator.",
      units = "(binary)",
      type = "binary",
      notes = "249 of 949 (26.2%) Japanese. Tested on CL/F and Vc/F in a SEPARATE model from the White/non-White split, to avoid confounding the race-effect estimate; not retained in the final model."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator.",
      units = "(binary)",
      type = "binary",
      notes = "96 of 949 (10.1%) Hispanic or Latino. Tested in a SEPARATE model from the other two race/ethnicity splits; not retained."
    ),
    AGE_GE65 = list(
      description = "Categorical age indicator: 1 = 65 years or older, 0 = under 65.",
      units = "(binary)",
      type = "binary",
      notes = "166 of 949 (17.5%) aged 65 or over. Tested as a categorical alternative to continuous AGE; the continuous power form was retained instead."
    ),
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F; not retained. Data from the concomitant-treatment periods of the phase 1 DDI studies (1202V9218 cyclosporine, 1403V921D rifampin, 1502V921E itraconazole/fluconazole) were EXCLUDED from the population PK analysis as worst-case scenarios not clinically relevant to OIC treatment, so this indicator is 1 for only 58 of 949 subjects (6.1%)."
    ),
    CONMED_CYP3A_INH = list(
      description = "Concomitant CYP3A inhibitor indicator (strong or moderate).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F; not retained. 14 subjects (1.5%) with a strong and 43 (4.5%) with a moderate inhibitor. See CONMED_PGP_INH for the DDI-period exclusion."
    ),
    CONMED_CYP3A_IND = list(
      description = "Concomitant CYP3A inducer indicator (strong or moderate).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F; not retained. 10 subjects (1.1%) with a strong and 6 (0.6%) with a moderate inducer. See CONMED_PGP_INH for the DDI-period exclusion."
    ),
    FORM_NALDEMEDINE = list(
      description = "Formulation category: solution or suspension, phase 1/2 tablet, or phase 3 tablet.",
      units = "(categorical)",
      type = "categorical",
      notes = "Tested on Ka only; not retained. Table II: solution or suspension 54 (5.7%), phase 1 or 2 tablet 130 (13.7%), phase 3 tablet 765 (80.6%)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 949L,
    n_studies = 18L,
    n_observations = paste(
      "8,146 naldemedine plasma concentrations from 949 subjects.",
      "9,077 samples were drawn from 1,026 subjects; 39 were excluded (14 not",
      "measured, 4 with an unidentified sampling or dosing time, 2 detectable",
      "before the first dose, 19 unexpected outliers) and 892 below-limit-of-",
      "quantification records (402 pre-dose, 490 post-dose) were treated as",
      "missing. LLOQ 0.01 ng/mL."
    ),
    age_range = "18-90 years (mean 51, SD 14, median 52); 783 (82.5%) under 65 years, 166 (17.5%) 65 or older",
    weight_range = "34.4-188.1 kg (mean 79.5, SD 23.4, median 76.0)",
    bmi_range = "14.4-58.8 kg/m^2 (mean 28.1, SD 7.2, median 26.9)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) 5.8-311.8 mL/min (mean 109.5, SD 42.4, median 108.0); includes a dedicated renal-impairment study (1401V921B, n = 38)",
    hepatic_function = "Total bilirubin 0.04-2.4 mg/dL, AST 6-223 U/L, ALT 2-212 U/L, albumin 2.1-5.4 g/dL; includes a dedicated hepatic-impairment study (1402V921C, n = 24)",
    sex_female_pct = 48.8,
    race_ethnicity = c(
      White = 58.8,
      Asian = 26.4,
      `Black or African American` = 13.5,
      `American Indian or Alaska Native` = 0.9,
      `Native Hawaiian or other Pacific Islander` = 0.2,
      Other = 0.1
    ),
    ethnicity_detail = "Japanese 249 (26.2%); Hispanic or Latino 96 (10.1%)",
    disease_state = paste(
      "Pooled across three populations: healthy subjects (phase 1), patients",
      "with chronic non-cancer pain and opioid-induced constipation (phase 2",
      "and phase 3), and cancer patients with opioid-induced constipation",
      "(phase 2b study 1108V9222 and phase 3 study 1331V9236). All patient",
      "cohorts were on chronic opioid therapy."
    ),
    dose_range = "0.01 to 3 mg oral naldemedine, single and multiple dose; the marketed clinical dose is 0.2 mg once daily",
    regions = "Global; phase 1 and phase 3 studies in Japan and internationally, with a large Japanese subgroup (26.2%)",
    notes = paste(
      "Analysis dataset pooled from 18 studies spanning phase 1 to phase 3",
      "(Table I). Data from the concomitant-treatment periods of the three",
      "phase 1 drug-drug-interaction studies (1202V9218 cyclosporine,",
      "1403V921D rifampin, 1502V921E itraconazole/fluconazole) were excluded",
      "because those regimens were designed as worst-case perpetrator exposures",
      "and are not clinically relevant to OIC treatment.",
      "Fitted in NONMEM 7.3 with FOCE-I; 200-replicate nonparametric bootstrap",
      "(98 runs, 49.0%, completed successfully) and a 1000-replicate",
      "prediction-corrected VPC were used for evaluation, with 7.4% of",
      "observations outside the 90% prediction interval.",
      "Kubota 2018 notes that observations above 100 h post-dose tended to",
      "exceed predictions; those concentrations are below 1/100 of Cmax and the",
      "authors judged the misfit not clinically meaningful."
    )
  )

  ini({
    # ======================================================================
    # Final-model estimates, Kubota 2018 Table III. The covariate functional
    # forms are the paper's Eq. (1) for continuous covariates
    #
    #   PKP = theta1 * (COV / median of COV)^theta2
    #
    # and Eq. (2) for categorical covariates
    #
    #   PKP = theta_(CAT=0) * (theta_CAT_i)^CAT_i
    #
    # so every categorical effect below is a MULTIPLICATIVE FACTOR raised to
    # a 0/1 indicator, not a log-scale shift. The centring constants (52
    # years, 108 mL/min, 76 kg) are the cohort medians of Table II, exactly
    # as Eq. (1) prescribes.
    #
    # Reference subject (all indicators at their source-zero level): a
    # 52-year-old, 76 kg, White, male patient with chronic non-cancer pain and
    # OIC, creatinine clearance 108 mL/min, dosed fasted. The paper reports
    # CL/F 9.10 L/h and Vc/F 91.1 L for that subject; note that Vc/F 91.1 is
    # NOT theta(6) but theta(6) * 1.20 = 75.9 * 1.20, because the reference
    # subject is a chronic-non-cancer-pain patient rather than a healthy
    # volunteer. The vignette re-derives both anchors mechanically.
    # ======================================================================

    # ----- Apparent oral clearance CL/F (L/h) -----
    lcl <- log(9.10); label("Apparent oral clearance CL/F (L/h) in the reference subject")  # Table III THETA(1) = 9.10 L/h (95% CI 8.73-9.47; bootstrap median 9.15)
    e_age_cl <- -0.195; label("Power exponent of median-centred age (AGE/52) on CL/F (unitless)")  # Table III THETA(2) = -0.195 (95% CI -0.291 to -0.0986)
    e_crcl_cl <- 0.0739; label("Power exponent of median-centred creatinine clearance (CRCL/108) on CL/F (unitless)")  # Table III THETA(3) = 0.0739 (95% CI -0.0133 to 0.161; CI spans zero, retained per Discussion)
    e_nonwhite_cl <- 0.870; label("Multiplicative factor on CL/F for NON-White race, applied as a power of (1 - RACE_WHITE) (unitless)")  # Table III THETA(4) = 0.870 (95% CI 0.820-0.920); Table II footnote codes the source indicator 'White' as 1 for non-White
    e_female_cl <- 0.902; label("Multiplicative factor on CL/F for female sex, applied as a power of SEXF (unitless)")  # Table III THETA(5) = 0.902 (95% CI 0.857-0.947); Table II footnote 'Gender = 1 for female'

    # ----- Apparent central volume Vc/F (L) -----
    lvc <- log(75.9); label("Apparent central volume Vc/F (L) at 76 kg in a HEALTHY subject dosed fasted")  # Table III THETA(6) = 75.9 L (95% CI 73.2-78.6)
    e_wt_vc <- fixed(1); label("Power exponent of median-centred body weight (WT/76) on Vc/F (unitless, held at 1)")  # Table III writes this term as '* (Body weight/76)' with NO exponent, unlike the age and CLcr power terms; the exponent was held at 1 rather than estimated
    e_cncp_vc <- 1.20; label("Multiplicative factor on Vc/F for chronic-non-cancer-pain OIC patients versus healthy (unitless)")  # Table III THETA(7) = 1.20 (95% CI 1.12-1.28), the source's 'non-Cancer' indicator
    e_cancer_vc <- 1.27; label("Multiplicative factor on Vc/F for cancer OIC patients versus healthy (unitless)")  # Table III THETA(8) = 1.27 (95% CI 1.05-1.49), the source's 'Cancer' indicator
    e_fed_vc <- 1.12; label("Multiplicative factor on Vc/F for the fed state, applied as a power of FED (unitless)")  # Table III THETA(9) = 1.12 (95% CI 1.05-1.19)

    # ----- Absorption -----
    lka <- log(2.94); label("First-order absorption rate constant Ka (1/h) at the reference age")  # Table III THETA(10) = 2.94 1/h (95% CI 2.32-3.56)
    e_age_ka <- -1.16; label("Power exponent of median-centred age (AGE/52) on Ka (unitless)")  # Table III THETA(11) = -1.16 (95% CI -1.26 to -1.06)
    ltlag <- log(0.195); label("Absorption lag time ALAG (h)")  # Table III ALAG = 0.195 h (95% CI 0.188-0.202); no covariates, no IIV

    # ----- Distribution -----
    lq <- log(4.77); label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table III Q/F = 4.77 L/h (95% CI 4.16-5.38); no covariates tested (model unstable)
    lvp <- log(41.8); label("Apparent peripheral volume Vp/F (L)")  # Table III Vp/F = 41.8 L (95% CI 38.4-45.2); no covariates tested (model unstable)

    # ----- Inter-individual variability (exponential / log-normal) -----
    # Table III reports IIV as a percent CV. The reported CV is the
    # APPROXIMATE one, CV% = 100 * sqrt(omega^2), NOT the exact log-normal
    # CV% = 100 * sqrt(exp(omega^2) - 1). The discriminator is the symmetry of
    # the printed 95% CIs: NONMEM's covariance step gives a symmetric interval
    # on the estimated variance, so back-transforming the CI endpoints with the
    # correct formula must restore that symmetry. On the two widest intervals
    # -- Ka (161.2%, 142.7-177.9) and Q/F (46.3%, 29.9-58.2), where the two
    # readings diverge most -- omega^2 = (CV/100)^2 gives upper/lower half-width
    # ratios of 1.007 and 0.995, while omega^2 = log(1 + CV^2) gives 0.860 and
    # 0.898. The narrow intervals (CL/F, Vc/F, Vp/F) cannot discriminate because
    # both transforms are near-linear over them. The vignette re-runs this test
    # mechanically.
    #
    # Kubota 2018 reports no correlation between the random effects
    # (Supplemental Fig. S3: "There was no clear correlation between the
    # inter-individual variabilities"), so the omega matrix is left diagonal.
    # IIV on ALAG was removed by the authors to allow the estimation and
    # covariance routines to converge, so no etaltlag is defined.
    etalcl ~ 0.143641  # Table III IIV CL/F 37.9% CV (95% CI 35.2-40.5), shrinkage 6.6% -> 0.379^2
    etalvc ~ 0.064009  # Table III IIV Vc/F 25.3% CV (95% CI 20.8-29.2), shrinkage 40.5% -> 0.253^2
    etalka ~ 2.598544  # Table III IIV Ka 161.2% CV (95% CI 142.7-177.9), shrinkage 32.6% -> 1.612^2
    etalq ~ 0.214369  # Table III IIV Q/F 46.3% CV (95% CI 29.9-58.2), shrinkage 60.6% -> 0.463^2
    etalvp ~ 0.131769  # Table III IIV Vp/F 36.3% CV (95% CI 30.2-41.6), shrinkage 57.1% -> 0.363^2

    # ----- Residual unexplained variability -----
    propSd <- 0.257; label("Proportional residual error (fraction)")  # Table III intra-individual variability, proportional 25.7% CV (95% CI 24.5-26.9), shrinkage 10.8%
  })

  model({
    # ------------------------------------------------------------------
    # Health status is three-level in the source (healthy / chronic
    # non-cancer pain with OIC / cancer with OIC) and is carried here on the
    # two ratified canonicals DIS_HEALTHY and DIS_CANCER. The source's
    # 'non-Cancer' indicator -- 1 only for chronic-non-cancer-pain OIC
    # patients -- is recovered as the both-negative cell.
    # ------------------------------------------------------------------
    dis_cncp_oic <- (1 - DIS_HEALTHY) * (1 - DIS_CANCER)

    # ------------------------------------------------------------------
    # Individual parameters. Continuous covariates enter as median-centred
    # power terms (Eq. 1); categorical covariates as multiplicative factors
    # raised to their 0/1 indicator (Eq. 2). RACE_WHITE is inverted at the
    # call site because the source's indicator is 1 for NON-White.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (AGE / 52)^e_age_cl *
      (CRCL / 108)^e_crcl_cl *
      e_nonwhite_cl^(1 - RACE_WHITE) *
      e_female_cl^SEXF

    vc <- exp(lvc + etalvc) *
      (WT / 76)^e_wt_vc *
      e_cncp_vc^dis_cncp_oic *
      e_cancer_vc^DIS_CANCER *
      e_fed_vc^FED

    ka <- exp(lka + etalka) * (AGE / 52)^e_age_ka
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    tlag <- exp(ltlag)

    # ------------------------------------------------------------------
    # Two-compartment disposition with first-order absorption and an
    # absorption lag time.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # Observation. Doses are in mg and vc is in L, so central / vc is in
    # mg/L = ug/mL; the 1000x factor converts to the ng/mL scale on which
    # naldemedine concentrations and the residual error were reported.
    # Sanity anchor: the reference subject's steady-state daily AUC is
    # Dose / (CL/F) = 0.2 mg / 9.10 L/h = 0.021978 mg*h/L = 21.98 ng*h/mL,
    # which is the value Kubota 2018 reports in Results.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
