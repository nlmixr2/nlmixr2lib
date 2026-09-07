Sato_2024_sglt2_hba1c_mbma <- function() {
  description <- paste0(
    "MBMA. Unified class-level dose-response model-based meta-analysis of ",
    "HbA1c reduction across six sodium-glucose co-transporter-2 (SGLT2) ",
    "inhibitors -- canagliflozin, dapagliflozin, empagliflozin, ",
    "ipragliflozin, luseogliflozin and tofogliflozin -- in type 2 diabetes. ",
    "The paper's central device is a UGE-NORMALIZED DOSE: each drug's mg ",
    "dose is divided by a drug-specific 'reference dose', defined as the ",
    "dose producing 51.4 g/day of urinary glucose excretion (UGE) in healthy ",
    "phase I volunteers, which is the geometric mean of the six drugs' UGE at ",
    "their clinical doses. After that normalization a SINGLE sigmoid Emax ",
    "curve describes all six drugs. Emax carries five covariates (baseline ",
    "HbA1c, body weight, eGFR, drug-naive status and diabetes duration) plus ",
    "a canagliflozin-specific 1.33-fold potentiation, which the authors ",
    "attribute to canagliflozin's comparatively weak SGLT2-over-SGLT1 ",
    "selectivity; the placebo term carries a concomitant-antihyperglycemic ",
    "-medication effect. Fitted to 295 study-arm means from 83 published ",
    "phase II/III randomized trials of at least 12 weeks. The model is purely ",
    "ALGEBRAIC and time-independent -- it predicts the end-of-treatment ",
    "study-arm mean HbA1c change from baseline, has no ODE states and ",
    "consumes no rxode2 dose events. Variability is BETWEEN-STUDY ",
    "(inter-study, ISV), encoded as correlated study-level etas, so the model ",
    "simulates study-arm mean outcomes and is NOT suitable for ",
    "individual-subject simulation. Companion SGLT2-inhibitor MBMA on the ",
    "same endpoint: modellib('Yao_2023_sglt2_endpoints_mbma'), which instead ",
    "drives HbA1c through exposure (AUC) and an FPG turnover model."
  )

  reference <- paste(
    "Sato H, Ishikawa A, Yoshioka H, Jin R, Sano Y, Hisaka A.",
    "Model-based meta-analysis of HbA1c reduction across SGLT2 inhibitors",
    "using dose adjusted by urinary glucose excretion.",
    "Sci Rep. 2024 Oct 21;14(1):24695.",
    "doi:10.1038/s41598-024-76256-6.",
    sep = " "
  )
  vignette <- "Sato_2024_sglt2_hba1c_mbma"
  units <- list(
    time          = paste0(
      "week (placeholder). The model is a time-independent end-of-treatment ",
      "dose-response; Sato 2024 restricted the pool to trials of at least 12 ",
      "weeks and found study duration non-significant as a covariate."
    ),
    dosing        = paste0(
      "mg/day. Per-arm daily doses enter as the six DOSE_<drug>_MGD covariate ",
      "columns, NOT as rxode2 dose events; the model has no compartments."
    ),
    concentration = paste0(
      "%/arm (study-arm mean change in HbA1c from baseline, NGSP percentage ",
      "points, SIGNED so that Cc = -0.7 means a 0.7-point HbA1c REDUCTION). ",
      "Output Cc is NOT a drug concentration; the slash in this unit string ",
      "is required by checkModelConventions parsing."
    )
  )

  covariateData <- list(
    DOSE_CANA_MGD = list(
      description        = "Per-arm total daily canagliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level, not individual-level. Exactly one of the six ",
        "DOSE_<drug>_MGD columns is non-zero on any active record; a placebo ",
        "arm sets all six to 0, the normalized dose collapses to 0 and the ",
        "prediction reduces to the placebo term. Setting two columns non-zero ",
        "simultaneously is outside the source's calibration and would make the ",
        "model ADD the two normalized doses. This column also selects the ",
        "canagliflozin Emax potentiation (e_cana_emax): the model derives the ",
        "indicator as DOSE_CANA_MGD > 0. Sato 2024 instead keyed the CANA ",
        "indicator off the study ID, so placebo arms inside canagliflozin ",
        "trials also carried CANA = 1; the two definitions give IDENTICAL ",
        "predictions because the Emax term is multiplied by a normalized dose ",
        "of 0 on every placebo arm. Doses in the pooled trials: 50-300 mg/day ",
        "(Sato 2024 Supplementary Table S2)."
      ),
      source_name        = "Dose (mg), rows with Drug = canagliflozin (Sato 2024 Table S2)"
    ),
    DOSE_DAPA_MGD = list(
      description        = "Per-arm total daily dapagliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Study-arm-level. See DOSE_CANA_MGD for the mutual-exclusivity rule. Doses in the pooled trials: 1-50 mg/day (Sato 2024 Supplementary Table S2).",
      source_name        = "Dose (mg), rows with Drug = dapagliflozin (Sato 2024 Table S2)"
    ),
    DOSE_EMPA_MGD = list(
      description        = "Per-arm total daily empagliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level. See DOSE_CANA_MGD for the mutual-exclusivity rule. ",
        "Doses in the pooled trials: 1-50 mg/day (Sato 2024 Supplementary ",
        "Table S2). Same canonical column already used by ",
        "Baron_2016_empagliflozin.R and Riggs_2014_empagliflozin.R, where it ",
        "is an individual-level daily dose; here it is a study-arm mean."
      ),
      source_name        = "Dose (mg), rows with Drug = empagliflozin (Sato 2024 Table S2)"
    ),
    DOSE_IPRA_MGD = list(
      description        = "Per-arm total daily ipragliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Study-arm-level. See DOSE_CANA_MGD for the mutual-exclusivity rule. Doses in the pooled trials: 12.5-300 mg/day (Sato 2024 Supplementary Table S2).",
      source_name        = "Dose (mg), rows with Drug = ipragliflozin (Sato 2024 Table S2)"
    ),
    DOSE_LUSEO_MGD = list(
      description        = "Per-arm total daily luseogliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Study-arm-level. See DOSE_CANA_MGD for the mutual-exclusivity rule. Doses in the pooled trials: 0.5-10 mg/day (Sato 2024 Supplementary Table S2).",
      source_name        = "Dose (mg), rows with Drug = luseogliflozin (Sato 2024 Table S2)"
    ),
    DOSE_TOFO_MGD = list(
      description        = "Per-arm total daily tofogliflozin dose (mg/day); 0 for placebo arms and for arms given a different SGLT2 inhibitor.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Study-arm-level. See DOSE_CANA_MGD for the mutual-exclusivity rule. Doses in the pooled trials: 2.5-40 mg/day (Sato 2024 Supplementary Table S2).",
      source_name        = "Dose (mg), rows with Drug = tofogliflozin (Sato 2024 Table S2)"
    ),
    HBA1C = list(
      description        = "Per-arm mean BASELINE HbA1c before randomization (NGSP %).",
      units              = "% (NGSP)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level. Enters the Emax covariate model in ",
        "exponential-linear-deviation form exp(0.438 * (HBA1C - 8)), centered ",
        "on 8%; a higher baseline HbA1c gives a LARGER reduction. Observed ",
        "range across the 295 arms 7.2-9.1%, over which f(COV) spans ",
        "0.704-1.62 (Sato 2024 Table 2B). NOTE: Sato 2024 Table 2B misprints ",
        "the unit of this covariate as 'mg dL-1'; the tabulated 7.2-9.1 range ",
        "and Supplementary Table S2 column header ('Baseline HbA1c (%)') both ",
        "confirm NGSP percent. Distinct from the model OUTPUT Cc, which is the ",
        "CHANGE in HbA1c from this baseline."
      ),
      source_name        = "Baseline HbA1c (%) (Sato 2024 Table S2; Table 2B)"
    ),
    WT = list(
      description        = "Per-arm mean body weight (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level. Enters the Emax covariate model as the power form ",
        "(WT / 81.6)^-0.661, normalized to 81.6 kg; a HEAVIER arm gets a ",
        "SMALLER reduction. Observed range 61.0-96.7 kg, over which f(COV) ",
        "spans 1.21-0.894 (Sato 2024 Table 2B)."
      ),
      source_name        = "Body weight (kg) (Sato 2024 Table S2; Table 2B)"
    ),
    CRCL = list(
      description        = "Per-arm mean baseline estimated glomerular filtration rate, body-surface-area normalized (mL/min/1.73 m^2).",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level. Sato 2024 reports eGFR, which the canonical CRCL ",
        "column covers (CRCL is the canonical for BSA-normalized renal ",
        "function from either a creatinine-based estimate or a tracer-measured ",
        "GFR). Enters the Emax covariate model as the power form ",
        "(CRCL / 85.9)^0.821, normalized to 85.9 mL/min/1.73 m^2; POORER renal ",
        "function gives a SMALLER reduction, consistent with the known ",
        "attenuation of SGLT2-inhibitor efficacy in renal impairment. Observed ",
        "range 38.5-154.5, over which f(COV) spans 0.517-1.62 (Sato 2024 ",
        "Table 2B)."
      ),
      source_name        = "eGFR (mL/min/1.73 m2) (Sato 2024 Table S2; Table 2B 'GFR')"
    ),
    T_DIAG_DIAB = list(
      description        = "Per-arm mean time since type 2 diabetes diagnosis (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Study-arm-level. Enters the Emax covariate model as the linear form ",
        "1 - 0.025 * (T_DIAG_DIAB - 6.6), centered on 6.6 years; LONGER ",
        "disease duration gives a SMALLER reduction. Observed range 0.25-18.2 ",
        "years, over which f(COV) spans 1.16-0.71 (Sato 2024 Table 2B). The ",
        "form is linear, not power or exponential, so it is only valid inside ",
        "roughly the observed range -- it crosses zero at 46.6 years and turns ",
        "the drug effect backwards beyond that."
      ),
      source_name        = "Diabetic duration (years) (Sato 2024 Table S2; Table 2B)"
    ),
    TRT_T2DM_NAIVE = list(
      description        = "1 = the study arm enrolled patients naive to antihyperglycemic therapy; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = paste0(
        "Study-arm-level. Sato 2024 Table 2B reports this covariate as ",
        "'Pre-treatment: Yes: 1; No: 1.230', i.e. arms with NO prior treatment ",
        "get a 1.230-fold LARGER Emax. The indicator in the source data is ",
        "drug-naive, despite the Supplementary Table S2 column being headed ",
        "'Pre-treatment': the authors' own analysis script (Supplementary File ",
        "MOESM1, line reading NAIVE = `Drug naive`) reads that column as the ",
        "naive flag, the arms coded 1 are drug-naive trials (e.g. Rosenstock ",
        "2016, 'Initial Combination Therapy ... for Drug-Naive Type 2 ",
        "Diabetes') while insulin add-on trials are coded 0, and the Results ",
        "text states the greater effect went with 'no prior treatment'. ",
        "Together with TRT_T2DM_ADDON this reproduces three of the four levels ",
        "of the canonical TRT_T2DM_* family: naive (this column 1, ADDON 0), ",
        "non-naive/washed-out (both 0) and add-on (ADDON 1); the 295 arms ",
        "split 51 / 53 / 191 and the naive and add-on flags never co-occur. ",
        "TRT_T2DM_MIXED is unused because every Sato 2024 arm resolves its ",
        "treatment history."
      ),
      source_name        = "Pre-treatment (Sato 2024 Table S2, read as drug-naive per Supplementary File MOESM1)"
    ),
    TRT_T2DM_ADDON = list(
      description        = "1 = the study arm received the SGLT2 inhibitor (or placebo) on top of ongoing background antihyperglycemic therapy; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other treatment-history stratum)",
      notes              = paste0(
        "Study-arm-level. Sato 2024 Table 2B reports this covariate on the ",
        "PLACEBO term as 'Concomitant medications: Yes: 0; No: 0.156', i.e. ",
        "arms WITHOUT background therapy get +0.156 added to Base, so their ",
        "placebo response is a slight HbA1c RISE (-0.124 + 0.156 = +0.032) ",
        "whereas add-on arms fall by 0.124 points. Derived in the authors' ",
        "script as ADD = 0 when the Supplementary Table S2 'Concomitant drug' ",
        "cell is empty and 1 otherwise (metformin, sulfonylurea, insulin, ",
        "DPP-4 inhibitors, thiazolidinediones, glinides, alpha-glucosidase ",
        "inhibitors and combinations). See TRT_T2DM_NAIVE for how the two ",
        "flags reconstruct the canonical four-level family."
      ),
      source_name        = "Concomitant drug (Sato 2024 Table S2; Table 2B 'Concomitant medications')"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Per-arm mean baseline age. Screened in the Sato 2024 stepwise covariate search but NOT retained in the final model.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists baseline age among the candidate covariates; it did not survive forward inclusion at p < 0.01 / backward elimination at p < 0.001. No point estimate is published, so no effect parameter appears in ini(). Tabulated per arm in Supplementary Table S2.",
      source_name        = "Age (years) (Sato 2024 Table S2)"
    ),
    SEXF = list(
      description        = "Per-arm proportion female. Screened but NOT retained.",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists sex among the candidate covariates; not retained. Supplementary Table S2 tabulates the complement, 'Male (%)'.",
      source_name        = "Male (%) (Sato 2024 Table S2)"
    ),
    BMI = list(
      description        = "Per-arm mean baseline body mass index. Screened but NOT retained.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists BMI among the candidate covariates; body weight was retained instead.",
      source_name        = "BMI (kg/m2) (Sato 2024 Table S2)"
    ),
    FPG = list(
      description        = "Per-arm mean baseline fasting plasma glucose. Screened but NOT retained.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists fasting plasma glucose among the candidate covariates; baseline HbA1c was retained instead.",
      source_name        = "FPG (mg/dL) (Sato 2024 Table S2)"
    ),
    SBP = list(
      description        = "Per-arm mean baseline systolic blood pressure. Screened but NOT retained.",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists systolic blood pressure among the candidate covariates; not retained.",
      source_name        = "SBP (mmHg) (Sato 2024 Table S2)"
    ),
    DBP = list(
      description        = "Per-arm mean baseline diastolic blood pressure. Screened but NOT retained.",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists diastolic blood pressure among the candidate covariates; not retained.",
      source_name        = "DBP (mmHg) (Sato 2024 Table S2)"
    ),
    TRT_DURATION = list(
      description        = "Trial treatment duration. Screened but NOT retained.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sato 2024 Methods lists study duration among the candidate covariates; not retained. This is why the final model is time-independent: the pool is restricted to trials of at least 12 weeks (range 12-104 weeks) and duration did not explain inter-study variability.",
      source_name        = "Study duration (weeks) (Sato 2024 Table S2)"
    ),
    REGION_JAPAN = list(
      description        = "1 = the trial was conducted in Japan. Screened but NOT retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (conducted outside Japan)",
      notes              = "Sato 2024 Methods lists 'whether the study was conducted in Japan' among the candidate covariates; not retained. Relevant because ipragliflozin, luseogliflozin and tofogliflozin are studied almost exclusively in Japan.",
      source_name        = "Study in Japan (Sato 2024 Table S2)"
    )
  )

  population <- list(
    species        = "human",
    n_studies      = 83L,
    n_data_points  = 295L,
    age_range      = "per-arm mean age approximately 51-66 years (Sato 2024 Supplementary Table S2)",
    weight_range   = "per-arm mean body weight 61.0-96.7 kg (Sato 2024 Table 2B covariate range)",
    disease_state  = paste0(
      "Type 2 diabetes mellitus. Per-arm mean baseline HbA1c 7.2-9.1% ",
      "(NGSP), per-arm mean eGFR 38.5-154.5 mL/min/1.73 m^2 (so the pool ",
      "spans normal renal function through moderate-to-severe impairment, ",
      "including dedicated renal-impairment trials), per-arm mean diabetes ",
      "duration 0.25-18.2 years. Treatment history across the 295 arms: 51 ",
      "drug-naive, 53 previously treated but on monotherapy in the trial, 191 ",
      "add-on to background antihyperglycemic therapy."
    ),
    dose_range     = paste0(
      "Placebo plus canagliflozin 50-300, dapagliflozin 1-50, empagliflozin ",
      "1-50, ipragliflozin 12.5-300, luseogliflozin 0.5-10 and tofogliflozin ",
      "2.5-40 mg/day; equivalently UGE-normalized doses of roughly 0.07-3.1. ",
      "Treatment durations 12-104 weeks."
    ),
    regions        = paste0(
      "International. Canagliflozin, dapagliflozin and empagliflozin trials ",
      "were mostly conducted outside Japan; ipragliflozin, luseogliflozin and ",
      "tofogliflozin trials were mostly Japanese (Sato 2024 Supplementary ",
      "Table S2)."
    ),
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "study-arm mean HbA1c change from baseline, not an individual ",
      "measurement. 83 trials contributed 295 arms after a PubMed search of ",
      "phase II/III trials of at least 12 weeks reporting HbA1c (137 studies ",
      "screened). Per-drug study counts: canagliflozin 17, dapagliflozin 29, ",
      "empagliflozin 22, ipragliflozin 7, luseogliflozin 5, tofogliflozin 3. ",
      "The dose-normalization layer is calibrated on a SEPARATE population -- ",
      "healthy Japanese volunteers in six single-dose phase I UGE studies ",
      "(n = 32-57 per study, mean age 23-27 years, mean body weight 61-65 kg; ",
      "Sato 2024 Supplementary Table S1) -- so the reference doses carry the ",
      "assumption, defended at length in the paper's Discussion, that the ",
      "dose-UGE relationship ranks the six drugs the same way in healthy ",
      "volunteers as in patients."
    )
  )

  ini({
    # ========================================================================
    # Sato 2024 Eq. 6 (sigmoid Emax on the UGE-normalized dose):
    #
    #   HbA1c_change = Base + Emax * nDose^n / (ED50^n + nDose^n)
    #                  + SD / sqrt(nsub) * eps
    #
    # with the parameter / covariate models of Eq. 7-10:
    #
    #   Base = theta1 + sum f(COV) + eta1        (ADDITIVE covariates + eta)
    #   Emax = theta2 * prod f(COV) * exp(eta2)  (MULTIPLICATIVE covariates,
    #                                             EXPONENTIAL eta)
    #   ED50 = theta3
    #   n    = theta4
    #
    # Point estimates are Sato 2024 Table 2A / 2B; Supplementary Table S3(A)
    # gives the same values with standard errors, and S3(B) the 488-sample
    # bootstrap summary. Sign convention: HbA1c change is SIGNED, negative =
    # reduction, so both Base and Emax are negative.
    # ========================================================================

    # ----- Placebo term (Sato 2024 Table 2A "Base") -------------------------
    pmax <- -0.124
    label("Placebo HbA1c change from baseline (% NGSP; signed, negative = reduction)")  # Sato 2024 Table 2A Base = -0.124 (bootstrap 95% CI -0.170 to -0.0718; SE 0.0261 per Table S3A)

    # ----- Maximum drug effect (Sato 2024 Table 2A "Emax") ------------------
    # Bare signed name rather than the canonical log form lemax: Emax is
    # NEGATIVE (it is a reduction on a signed change scale), so it cannot be
    # log-transformed. Same convention as the signed MBMA parameters in
    # modellib('Li_2015_taspoglutide_mbma') (dmax_hb) and
    # modellib('Yao_2023_sglt2_endpoints_mbma') (slopefd).
    emax <- -0.796
    label("Maximum HbA1c reduction at saturating normalized dose (% NGSP; signed, negative = reduction)")  # Sato 2024 Table 2A Emax = -0.796 (bootstrap 95% CI -1.17 to -0.689; SE 0.071 per Table S3A)

    # ----- Potency and sigmoidicity (Sato 2024 Table 2A) --------------------
    # Log-transformed for positivity; exp(led50) recovers 0.251.
    # The paper's Discussion flags ED50 as the one unstable parameter in the
    # bootstrap (mean 0.495 vs median 0.258, SE 3.22, 95% CI 0.191-1.09), but
    # notes that outlying ED50 draws paired with small Hill coefficients so
    # the resulting curves were similar; the point estimate is used here.
    led50 <- log(0.251)
    label("Normalized dose giving half-maximal HbA1c reduction, ED50 (dimensionless, dose / reference dose)")  # Sato 2024 Table 2A ED50 = 0.251 (bootstrap 95% CI 0.191 to 1.09; SE 0.0539 per Table S3A)

    lhill <- log(0.662)
    label("Hill coefficient n of the normalized-dose response (unitless)")  # Sato 2024 Table 2A n = 0.662 (bootstrap 95% CI 0.388 to 0.981; SE 0.149 per Table S3A)

    # ----- Emax covariates (Sato 2024 Table 2B, "f(COV)" column) ------------
    # Every one of these multiplies Emax (Eq. 8).
    e_hba1c_emax <- 0.438
    label("Coefficient on (HBA1C - 8) in the exponential baseline-HbA1c effect on Emax (per % NGSP)")  # Sato 2024 Table 2B f(COV) = exp(0.438 * (COV - 8)); S3A SE 0.0719

    e_wt_emax <- -0.661
    label("Power exponent on (WT / 81.6) for Emax (unitless)")  # Sato 2024 Table 2B f(COV) = (COV / 81.6)^-0.661; S3A SE 0.175

    e_cana_emax <- 0.328
    label("Fractional increase in Emax for canagliflozin arms (unitless)")  # Sato 2024 Table 2B f(COV) = "Cana: 1.33; Others: 1"; S3A theta = 0.328 (SE 0.045) entering as (1 + 0.328) = 1.328, and the authors' script MOESM1 confirms EMAX * (1 + CANA)

    e_crcl_emax <- 0.821
    label("Power exponent on (CRCL / 85.9) for Emax (unitless)")  # Sato 2024 Table 2B f(COV) = (COV / 85.9)^0.821; S3A SE 0.102

    e_trt_t2dm_naive_emax <- 0.230
    label("Fractional increase in Emax for drug-naive arms (unitless)")  # Sato 2024 Table 2B f(COV) = "Yes: 1; No: 1.230" for Pre-treatment, i.e. (1 + 0.230) when the arm is drug-naive; S3A theta = 0.230 (SE 0.081)

    e_t_diag_diab_emax <- -0.025
    label("Linear coefficient on (T_DIAG_DIAB - 6.6) for Emax (per year)")  # Sato 2024 Table 2B f(COV) = 1 - 0.025 * (COV - 6.6), encoded here as 1 + (-0.025) * (COV - 6.6); see the sign note below

    # NOTE on the diabetes-duration sign. Table S3 reports this theta with
    # OPPOSITE signs in its two halves: +0.025 in (A) "Original data" and
    # -0.0247 (median -0.0253, 95% CI -0.041 to -0.0054) in (B) the bootstrap
    # summary. The FUNCTION is unambiguous either way: Table 2B prints
    # f(COV) = 1 - 0.025 * (COV - 6.6), and its printed range 1.16-0.71 over
    # COV 0.25-18.2 reproduces exactly under that form (1 - 0.025 * (0.25 -
    # 6.6) = 1.159; 1 - 0.025 * (18.2 - 6.6) = 0.710). This model encodes the
    # bootstrap's signed parameterization 1 + theta * (COV - 6.6) with
    # theta = -0.025, which is numerically identical to the printed f(COV).

    # ----- Placebo covariate (Sato 2024 Table 2B) ---------------------------
    # The only covariate on Base, and it is ADDITIVE (Eq. 7).
    e_trt_t2dm_addon_pmax <- 0.156
    label("Additive shift in the placebo HbA1c change for arms with NO background antihyperglycemic therapy (% NGSP)")  # Sato 2024 Table 2B f(COV) = "Yes: 0; No: 0.156"; S3A SE 0.0348

    # ========================================================================
    # Drug-specific UGE reference doses (Sato 2024 Fig. 1, DIGITIZED).
    #
    # PROVENANCE -- NOT PAPER-PRINTED. These six values are the only
    # quantities in this file that are not printed in the paper or its
    # supplement. Sato 2024 defines the reference dose as the dose at which
    # the drug's fitted linear-logarithmic dose-UGE curve reaches the
    # reference UGE of 51.4 g/day (Eq. 1-4), but publishes the values only as
    # the red dashed vertical lines in Figure 1; the underlying data file
    # (data_uge.csv, read by the authors' Supplementary File MOESM1 script)
    # is not distributed. They were recovered from the 600-dpi Figure 1
    # panels by locating the red dashed line against the x-axis tick
    # calibration, and cross-checked two independent ways:
    #   (1) extracting each panel's fitted black curve and solving
    #       a + b*log10(dose) = 51.4 -- agrees to within 1.7% on every drug
    #       (96.7 / 7.02 / 11.90 / 80.33 / 7.38 / 14.59 mg), with the
    #       log-linear form of Eq. 1 confirmed at R^2 = 0.982-0.999;
    #   (2) locating the clinical-dose markers by their exact palette colours
    #       from the authors' script -- the recovered x-positions reproduce
    #       the known clinical doses (100, 10, 100, 5, 20 mg) to within 0.1-1%,
    #       validating the x-axis calibration itself.
    # Digitization error is therefore of order 1-2%, which propagates to the
    # normalized dose and is far smaller than the parameter uncertainty on
    # ED50 (bootstrap 95% CI 0.191-1.09, a better-than-fourfold span).
    # See the vignette Assumptions and deviations section.
    # ========================================================================
    dref_cana <- fixed(97.7)
    label("Canagliflozin UGE reference dose, i.e. dose giving 51.4 g/day urinary glucose excretion (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 canagliflozin panel, red dashed line

    dref_dapa <- fixed(7.05)
    label("Dapagliflozin UGE reference dose (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 dapagliflozin panel, red dashed line

    dref_empa <- fixed(11.8)
    label("Empagliflozin UGE reference dose (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 empagliflozin panel, red dashed line

    dref_ipra <- fixed(79.0)
    label("Ipragliflozin UGE reference dose (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 ipragliflozin panel, red dashed line

    dref_luseo <- fixed(7.38)
    label("Luseogliflozin UGE reference dose (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 luseogliflozin panel, red dashed line

    dref_tofo <- fixed(15.0)
    label("Tofogliflozin UGE reference dose (mg/day; digitized, not paper-printed)")  # Sato 2024 Fig. 1 tofogliflozin panel, red dashed line

    # ========================================================================
    # Inter-study variability (Sato 2024 Table 2A, "ISV (omega)" column).
    #
    # omega is a STANDARD DEVIATION, not a variance: the paper's Eq. 7-10
    # narrative states "eta: ISV that follows the normal distribution with
    # mean 0 and variance omega^2". nlmixr2 ini() takes the variance, so the
    # tabulated 0.191 and 0.214 are squared here. The correlation coefficient
    # 0.785 is reported on its own row and gives the off-diagonal covariance
    # 0.785 * 0.191 * 0.214 = 0.0320861.
    #
    # These are BETWEEN-STUDY effects, not popPK between-subject effects; the
    # eta_study_* prefix marks that, following
    # modellib('Yao_2023_sglt2_endpoints_mbma') and
    # modellib('Yang_2010_rosuvastatin_mbma').
    #
    # The two etas enter on DIFFERENT scales, per Eq. 7 and 8: eta_study_pmax
    # is ADDITIVE on Base (units of % HbA1c) while eta_study_emax is
    # EXPONENTIAL on Emax (a 21.4% CV, preserving Emax's negative sign).
    # ========================================================================
    # Block entries below, in nlmixr2 lower-triangle order:
    #   var(eta_study_pmax) = 0.191^2 = 0.036481   (Sato 2024 Table 2A ISV Base)
    #   cov                 = 0.785 * 0.191 * 0.214 = 0.0320861
    #   var(eta_study_emax) = 0.214^2 = 0.045796   (Sato 2024 Table 2A ISV Emax)
    # Comments must stay OUTSIDE the c() -- a comment inside an omega c()
    # breaks the rxode2 comment-to-label parser.
    eta_study_pmax + eta_study_emax ~ c(0.036481, 0.0320861, 0.045796)

    # ========================================================================
    # Residual error (Sato 2024 Eq. 6).
    #
    # The paper does NOT estimate a residual magnitude: "the standard
    # deviation of the residual error (eps) was fixed at 1 and scaled by the
    # standard error of the observed HbA1c change", so the residual term is
    # SD_arm / sqrt(n_arm) * eps with eps ~ N(0, 1). That per-record weight is
    # a property of the DATA (each published arm's own reported SD and sample
    # size), not of the model, so it cannot live in ini().
    #
    # This file therefore encodes the unit-weight residual sigma = 1 exactly
    # as the paper fixes it -- the same convention every other size-weighted
    # MBMA in this package uses (see the "unit study weight" note in
    # modellib('Yao_2023_sglt2_endpoints_mbma')).
    #
    # CONSEQUENCE FOR SIMULATION: a bare stochastic rxSolve() applies an
    # additive residual SD of 1 HbA1c point, which is roughly 10-40x too
    # large for a real trial arm. To simulate a specific arm, multiply by that
    # arm's standard error -- e.g. an arm of 250 patients with a within-arm SD
    # of 1.0% has SE = 1.0 / sqrt(250) = 0.063 -- or drop the residual
    # entirely with rxode2::zeroRe(mod, "sigma") when reproducing the paper's
    # population curves, which is what the vignette does.
    # ========================================================================
    addSd <- fixed(1)
    label("Residual SD on the standardized scale; the source scales it per record by each arm's standard error SD/sqrt(n)")  # Sato 2024 Eq. 6 text: "the standard deviation of the residual error was fixed at 1 and scaled by the standard error of the observed HbA1c change"
  })

  model({
    # ======================================================================
    # 1. UGE-normalized dose (Sato 2024 Eq. 5: nDose = Dose / Dose_reference)
    #
    # Exactly one of the six dose columns is non-zero on an active arm, so
    # this sum picks out that drug's normalized dose; a placebo arm sets all
    # six to 0 and the sum is 0. Same shape as the multi-drug sum in
    # modellib('Yao_2023_sglt2_endpoints_mbma').
    # ======================================================================
    ndose <-
      DOSE_CANA_MGD  / dref_cana +
      DOSE_DAPA_MGD  / dref_dapa +
      DOSE_EMPA_MGD  / dref_empa +
      DOSE_IPRA_MGD  / dref_ipra +
      DOSE_LUSEO_MGD / dref_luseo +
      DOSE_TOFO_MGD  / dref_tofo

    # Canagliflozin indicator for the Emax potentiation. Equivalent to Sato
    # 2024's study-level CANA flag on every record that can affect the
    # prediction; see the DOSE_CANA_MGD covariate note.
    canaFlag <- (DOSE_CANA_MGD > 0)

    ed50 <- exp(led50)
    hill <- exp(lhill)

    # ======================================================================
    # 2. Placebo term (Sato 2024 Eq. 7: Base = theta1 + sum f(COV) + eta1)
    #
    # ADDITIVE covariate and ADDITIVE eta. The covariate is keyed on the
    # ABSENCE of background therapy, matching Table 2B's "Yes: 0; No: 0.156".
    # ======================================================================
    base <- pmax +
      e_trt_t2dm_addon_pmax * (1 - TRT_T2DM_ADDON) +
      eta_study_pmax

    # ======================================================================
    # 3. Maximum effect (Sato 2024 Eq. 8:
    #    Emax = theta2 * prod f(COV) * exp(eta2))
    #
    # MULTIPLICATIVE covariates and an EXPONENTIAL eta, so the product
    # preserves the negative sign of theta2 for every covariate combination
    # inside the observed ranges.
    # ======================================================================
    emaxArm <- emax *
      exp(e_hba1c_emax * (HBA1C - 8)) *
      (WT / 81.6)^e_wt_emax *
      (1 + e_cana_emax * canaFlag) *
      (CRCL / 85.9)^e_crcl_emax *
      (1 + e_trt_t2dm_naive_emax * TRT_T2DM_NAIVE) *
      (1 + e_t_diag_diab_emax * (T_DIAG_DIAB - 6.6)) *
      exp(eta_study_emax)

    # ======================================================================
    # 4. Sigmoid Emax prediction (Sato 2024 Eq. 6)
    #
    # Cc is the study-arm mean HbA1c CHANGE FROM BASELINE in NGSP percentage
    # points, signed so that a negative value is a reduction. It is not a
    # concentration; the canonical output name is kept per package convention
    # and per modellib('Yang_2010_rosuvastatin_mbma'), the sibling algebraic
    # dose-response MBMA.
    # ======================================================================
    Cc <- base + emaxArm * ndose^hill / (ed50^hill + ndose^hill)

    Cc ~ add(addSd)
  })
}
