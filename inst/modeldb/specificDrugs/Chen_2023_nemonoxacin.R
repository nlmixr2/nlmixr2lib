Chen_2023_nemonoxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral nemonoxacin capsules in",
    "Chinese adults, pooled across phase I to III trials (Chen 2023;",
    "n = 161 subjects / 195 cases; 2007 plasma concentrations). First-order",
    "absorption with a lag time and estimated relative bioavailability.",
    "Clearance is an additive linear function of Cockcroft-Gault creatinine",
    "clearance (10.3 L/h intercept + 0.026 L/h per mL/min), and allometric",
    "body-weight scaling is fixed at 0.75 on CL / Q and 1 on Vc / Vp with a",
    "70 kg reference. Covariate effects: female sex lowers Vc by 11%;",
    "community-acquired pneumonia raises Vp by 23%; food slows absorption",
    "(ka x 0.44), lengthens the lag time (T_lag x 1.6), and lowers",
    "bioavailability (F x 0.88). Inter-occasion variability on CL",
    "distinguishes the single-dose occasion from steady state (72 h after",
    "the first of multiple doses). The paper's companion PK/PD analysis is",
    "an exposure-response logistic regression on AUC0-24/MIC, Cmax/MIC and",
    "%T > MIC; its regression coefficients are not reported, so only the",
    "population PK layer is encoded here (the published PK/PD targets are",
    "reproduced in the validation vignette)."
  )
  reference <- paste(
    "Chen Y, Wu X, Tsai C, Chang L, Yu J, Cao G, Guo B, Shi Y, Zhu D, Hu F,",
    "Yuan J, Liu Y, Zhao X, Zhang Y, Wu J, Zhang J (2023). Integrative",
    "population pharmacokinetic/pharmacodynamic analysis of nemonoxacin",
    "capsule in Chinese patients with community-acquired pneumonia.",
    "Frontiers in Pharmacology 14:912962. doi:10.3389/fphar.2023.912962.",
    sep = " "
  )
  vignette <- "Chen_2023_nemonoxacin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives allometric scaling of CL and Q (exponent 0.75) and of Vc",
        "and Vp (exponent 1), both fixed a priori rather than estimated,",
        "on a 70 kg reference weight. Chen 2023 Methods 'Fixed effect",
        "model' final paragraph: 'body weight (BW) was added on CL, V2, V3",
        "and Q using power model ... The power exponent was set as 0.75, 1,",
        "1 and 0.75, respectively, according to the literature (Germovsek",
        "et al., 2017; Holford & Anderson, 2017; Sinha et al., 2019). The",
        "power base was 70 kg (Germovsek et al., 2019).' Table 2 footnote:",
        "'Parameter values are for a 70 kg adult.' Cohort weight 61.0 +/-",
        "9.0 kg (median 60, range 42-90) per Table 1."
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CLcr. Chen 2023 Methods 'Data processing': 'The",
        "serum creatinine was measured using enzymatic method, and",
        "creatinine clearance (CLcr) was calculated by Cockcroft Gault",
        "formula'; Table 1 footnote a: 'CLcr was estimated by",
        "Cockcroft-Gault formula'. The paper enters CLcr in raw mL/min",
        "without BSA normalization, matching the raw-mL/min variant of the",
        "canonical CRCL column (see the CRCL entry in",
        "inst/references/covariate-columns.md and the",
        "Georges_2009_ceftazidime.R / Delattre_2010_amikacin.R",
        "precedents). The effect on CL is additive linear with no",
        "centering: Eq. 8 'CL = (10.3 + 0.026 x CLcr) x (BW/70)^0.75 x",
        "exp(eta1 + kappa)'. Cohort CLcr 114.4 +/- 29.2 mL/min (median",
        "116.2, range 50.7-200.7) per Table 1. CAP patients with moderate",
        "or severe renal impairment were NOT enrolled, so CLcr values far",
        "below ~50 mL/min are an extrapolation (Discussion 'Some",
        "limitations', third point)."
      ),
      source_name        = "CLcr"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; structural reference for the Table 2 typical Vc of 103 L)",
      notes              = paste(
        "Time-fixed per subject. Multiplicative power-form effect on Vc:",
        "0.89^SEXF per Eq. 8 'V2 = 103 x (BW/70) x 0.89^Sex x exp(eta2)'.",
        "Chen 2023 states the coding explicitly after Eq. 8: 'Value of sex",
        "was 0 for male and 1 for female', so the source column maps",
        "one-to-one onto the canonical SEXF orientation with no value",
        "transformation. Results 'PPK model': 'central volume of",
        "distribution (V2) decreased by 11% in female subjects compared to",
        "male subjects.' Cohort 111 male / 84 female (Table 1)."
      ),
      source_name        = "SEX"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (1 = healthy subject, 0 = CAP patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (community-acquired pneumonia patient); the model file shifts the",
        "structural lvp onto this reference so DIS_HEALTHY = 0 reproduces",
        "Eq. 8's CAP-patient Vp of 28 x 1.23 = 34.4 L"
      ),
      notes              = paste(
        "Time-fixed per subject. Chen 2023 uses the reverse-coded DisStat",
        "indicator (stated after Eq. 8: 'DisStat indicates disease status,",
        "0 for healthy subjects and 1 for CAP patients'), which is",
        "re-expressed here as DIS_HEALTHY = 1 - DisStat per the canonical",
        "orientation (same re-expression as Yoneyama_2017_emicizumab.R,",
        "Taubert_2018_finafloxacin.R and Bulitta_2010_ceftazidime.R).",
        "Multiplicative power-form effect on Vp: Eq. 8 'V3 = 28 x (BW/70) x",
        "1.23^DisStat x exp(eta3)', i.e. e_dis_healthy_vp = 1 / 1.23 =",
        "0.813 applied as 0.813^DIS_HEALTHY on the patient-referenced lvp.",
        "Results 'PPK model': 'peripheral volume of distribution (V3)",
        "increased by 23% in CAP patients compared to healthy subjects.'",
        "The complement reference category is the pooled Chinese CAP",
        "cohort from the phase II (TG-873870-C-3) and phase III",
        "(TG-873870-C-4) trials; the healthy stratum is the phase I",
        "(TG-873870-C-1) cohort. 125 CAP patients and 36 healthy subjects",
        "(Table 1)."
      ),
      source_name        = "DisStat"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (fasted; structural reference for the Table 2 typical ka of 2.2",
        "1/h, T_lag of 0.19 h, and F anchored at 1)"
      ),
      notes              = paste(
        "Per-dose-record indicator. Chen 2023 defines the food window",
        "explicitly after Eq. 8: 'Value of food was 1 if taking food within",
        "2 h before administration or within 30 min after administration,",
        "otherwise the value was 0' - a general fed-vs-fasted flag rather",
        "than a high-fat-meal challenge, so the canonical FED column",
        "applies (not FED_HIGHFAT). Three multiplicative power-form effects",
        "per Eq. 8: ka x 0.44^Food, T_lag x 1.6^Food, F1 x 0.88^Food.",
        "Results 'PPK model': 'absorption rate of nemonoxacin under fed",
        "condition was only 44% of that under fasting condition. Meanwhile,",
        "the lag in absorption time increased by 60%, and the absolute",
        "bioavailability reduced by 12%.' Cohort 91 fed / 104 fasted",
        "(Table 1)."
      ),
      source_name        = "Food"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for inter-occasion variability on CL",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Two occasions. Chen 2023 Methods 'Base model': 'Since PK",
        "parameters of nemonoxacin may change at steady state (72 h",
        "following multiple doses) compared to the single dose",
        "administration, parameter kappa was introduced to describe the",
        "variation of PK parameters over time (IOV) (Karlsson & Sheiner,",
        "1993).' OCC = 1 identifies the single-dose occasion (before 72 h",
        "of multiple dosing) and OCC = 2 the steady-state occasion (72 h",
        "onward); Table 2 reports two eta-shrinkage values ('45 or 20') for",
        "the single pi^2 CL term, consistent with two occasions. Decomposed",
        "inside model() into binary indicators oc1 / oc2 that multiplex the",
        "IOV etas on log-CL, following the Jonsson_2011_ethambutol.R",
        "precedent. For single-occasion simulation pass OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 161L,
    n_studies      = 3L,
    age_range      = "18 to 70 years",
    age_median     = "27.5 years (mean 34.0 +/- 14.4)",
    weight_range   = "42 to 90 kg",
    weight_median  = "60 kg (mean 61.0 +/- 9.0)",
    sex_female_pct = 100 * 84 / 195,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Pooled Chinese cohort of 125 adults with community-acquired",
      "pneumonia (CAP) and 36 healthy subjects. Among CAP patients 38",
      "(30.4%) had an underlying disease (hypertension 6, COPD 3, chronic",
      "bronchitis 3, emphysema 3 of 49 total underlying-disease",
      "occurrences), 12 (9.6%) had a relevant disease history, and 55",
      "(44.0%) took a concomitant medication (most often ambroxol 36,",
      "licorice 5, codeine oral liquid 4 of 94 total generic-name",
      "occurrences). Patients >= 75 years of age or with moderate to",
      "severe renal impairment were excluded from all three trials."
    ),
    renal_function = paste(
      "Cockcroft-Gault CLcr 114.4 +/- 29.2 mL/min (median 116.2, range",
      "50.7-200.7) overall; 103.8 +/- 27.9 (median 100.4, range",
      "50.7-187.5) in CAP patients and 130.0 +/- 17.3 (median 127.2,",
      "range 100.9-185.3) in healthy subjects. No subject had moderate or",
      "severe renal impairment."
    ),
    dose_range     = paste(
      "Oral nemonoxacin capsules (250 mg capsules, Huayu Wuxi",
      "Pharmaceutical Co., Ltd.). Doses represented in the pooled PK",
      "dataset span 250, 500 and 750 mg (Results 'PPK model': 'VPC results",
      "for 750 mg and 250 mg of nemonoxacin were similar to that of 500 mg",
      "dose'), including both single-dose and multiple-dose (q24 h)",
      "administration; the CAP dosing regimens evaluated by Monte Carlo",
      "simulation were 500 mg and 750 mg q24 h for 10 days."
    ),
    regions        = "China.",
    notes          = paste(
      "Data pooled from three trials of oral nemonoxacin capsule: phase I",
      "(TG-873870-C-1), phase II in CAP (TG-873870-C-3), and phase III vs",
      "levofloxacin in CAP (TG-873870-C-4). 161 subjects contributed 195",
      "cases and 2027 plasma concentrations, of which 2007 (99.0%) were",
      "retained after removing below-LLOQ values, Chauvenet-criterion",
      "outliers, diagnostic-plot outliers, and records from",
      "poorly-compliant subjects. Assay: LC-MS/MS, LLOQ 0.005 mg/L, linear",
      "range 0.005-1 mg/L. NONMEM 7.4 with FO then FOCEI; evaluated by 300",
      "bootstrap replicates, VPC (1000 simulations per subject) and a",
      "1000-shuffle covariate randomization test. Baseline demographics",
      "from Table 1; final parameter estimates from Table 2; final",
      "covariate model from Eq. 8. Albumin and hemoglobin in Table 1 are",
      "reported as the paper's normal-range-scaled ratio measured",
      "value / [0.5 x (upper + lower limit of normal)] and were screened",
      "but not retained in the final model."
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Identified as a significant covariate during FO screening but",
        "dropped under FOCEI and absent from the final model. Results 'PPK",
        "model': 'Age, bilirubin and concomitant medication were not",
        "included in the final model due to their variable effects on PK",
        "parameters.'"
      )
    ),
    ALB = list(
      description = "Serum albumin, expressed as the paper's normal-range-scaled ratio",
      units       = "(ratio to normal-range midpoint)",
      type        = "continuous",
      notes       = paste(
        "Screened as significant under FO ('seven covariates (CLcr,",
        "weight, sex, albumin, hemoglobin, food, and age) were preliminarily",
        "identified as significant covariates') but not retained after",
        "FOCEI screening, which confirmed only CLcr, body weight, DisStat,",
        "sex and food. Reported in Table 1 as measured value / [0.5 x",
        "(upper + lower limit of normal range)]: 1.04 +/- 0.12 overall."
      )
    ),
    HGB = list(
      description = "Hemoglobin, expressed as the paper's normal-range-scaled ratio",
      units       = "(ratio to normal-range midpoint)",
      type        = "continuous",
      notes       = paste(
        "Screened as significant under FO but not retained after FOCEI",
        "screening. Reported in Table 1 as the normal-range-scaled ratio:",
        "1.00 +/- 0.10 overall."
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "(ratio to normal-range midpoint)",
      type        = "continuous",
      notes       = paste(
        "Explicitly excluded from the final model: 'Age, bilirubin and",
        "concomitant medication were not included in the final model due to",
        "their variable effects on PK parameters.'"
      )
    ),
    NEUT = list(
      description = "Absolute neutrophil count",
      units       = "(ratio to normal-range midpoint)",
      type        = "continuous",
      notes       = paste(
        "Reduced OBJ when added on V3 (FO, dOBJ 8.1) and on CL and T_lag",
        "(FOCEI, dOBJ 15.3 and 8.3) but deliberately excluded because",
        "'neutrophil was often used as the PD index in clinical studies'",
        "(Discussion)."
      )
    ),
    RESPRATE = list(
      description = "Respiratory rate",
      units       = "breaths/min",
      type        = "continuous",
      notes       = paste(
        "Reduced OBJ by 19.3 on CL and 11.0 on V2 under FO but excluded",
        "because 'the relationship between respiratory rate and drug",
        "clearance could not be explained' (Discussion)."
      )
    ),
    GLUC = list(
      description = "Serum glucose",
      units       = "(ratio to normal-range midpoint)",
      type        = "continuous",
      notes       = paste(
        "Reduced OBJ significantly when added on V1 or Ka but excluded",
        "because 'to our knowledge, glucose was not regarded as a covariate",
        "on PK parameter' (Discussion)."
      )
    )
  )

  ini({
    # Structural parameters -- Chen 2023 Table 2 'Estimate based on original dataset'
    # column, with the covariate model as printed in Eq. 8. Table 2 footnote:
    # 'Parameter values are for a 70 kg adult.'
    #
    # Eq. 8 verbatim:
    #   CL   = (10.3 + 0.026 x CLcr) x (BW/70)^0.75 x exp(eta1 + kappa)   L/h
    #   V2   = 103 x (BW/70)         x 0.89^Sex     x exp(eta2)           L
    #   Q    = 2.0 x (BW/70)^0.75                                        L/h
    #   V3   = 28  x (BW/70)         x 1.23^DisStat x exp(eta3)           L
    #   KA   = 2.2 x 0.44^Food                     x exp(eta4)           1/h
    #   TLAG = 0.19 x 1.6^Food                     x exp(eta5)           h
    #   F1   = 1    x 0.88^Food                    x exp(eta6)

    # CL intercept: the renal-function-independent component of clearance, i.e. CL at
    # CLcr = 0 for a 70 kg adult. The additive linear CRCL term below carries the renal
    # component; this follows the Georges_2009_ceftazidime.R / Delattre_2010_amikacin.R
    # encoding of an intercept-plus-slope CL covariate model.
    lcl <- log(10.3);  label("CL intercept at CLcr = 0 and WT = 70 kg (L/h)")   # Table 2 row 'CL' = 10.3 L/h (bootstrap 10.4, RSE 1%); Eq. 8
    lvc <- log(103);   label("Central volume of distribution Vc at WT = 70 kg, male (L)") # Table 2 row 'V2' = 103 L (bootstrap 103, RSE 1%); Eq. 8
    lq  <- log(2.0);   label("Intercompartmental clearance Q at WT = 70 kg (L/h)")        # Table 2 row 'Q' = 2.0 L/h (bootstrap 2.0, RSE 2%); Eq. 8

    # Vp is referenced to the CAP-patient state (DIS_HEALTHY = 0) so the canonical
    # DIS_HEALTHY effect below is oriented as 'effect of being healthy'. Eq. 8 prints
    # the healthy typical V3 = 28 L with a 1.23-fold CAP-patient factor, so the
    # patient-referenced value is 28 x 1.23 = 34.44 L. Same re-expression as
    # Yoneyama_2017_emicizumab.R / Taubert_2018_finafloxacin.R / Bulitta_2010_ceftazidime.R.
    lvp <- log(28 * 1.23); label("Peripheral volume of distribution Vp at WT = 70 kg, CAP patient (L)") # Table 2 row 'V3' = 28 L (healthy) x Eq. 8 DisStat factor 1.23

    lka   <- log(2.2);  label("Absorption rate constant ka, fasted (1/h)")   # Table 2 row 'KA' = 2.2 1/h (bootstrap 2.3, RSE 11%); Eq. 8
    ltlag <- log(0.19); label("Absorption lag time T_lag, fasted (h)")       # Table 2 row 'TLAG' = 0.19 h (bootstrap 0.18, RSE 8%); Eq. 8

    # Bioavailability anchor: Eq. 8 prints 'F1 = 1 x 0.88^Food x exp(eta6)'. The
    # leading 1 is a structural anchor -- Table 2 reports the Food-on-F exponent and
    # the IIV on F but no estimated typical F1 -- so it is encoded as fixed.
    lfdepot <- fixed(log(1)); label("Relative bioavailability anchor F1, fasted (fraction)") # Eq. 8 'F1 = 1 x 0.88^Food'; no typical-value row in Table 2

    # Allometric exponents, fixed a priori (not estimated) per Methods 'Fixed effect
    # model': 'The power exponent was set as 0.75, 1, 1 and 0.75, respectively' for
    # CL, V2, V3 and Q, 'according to the literature'. 'The power base was 70 kg.'
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)") # Methods 'Fixed effect model', priori-information paragraph 2
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent of body weight on Vc and Vp (unitless)") # Methods 'Fixed effect model', priori-information paragraph 2

    # Covariate effects (Chen 2023 Table 2 'Estimate based on original dataset'; the
    # applied forms are as printed in Eq. 8).
    e_crcl_cl        <- 0.026;      label("Slope of the additive linear CLcr effect on CL (L/h per mL/min)") # Table 2 row 'CLcr (CL)' = 0.026 (bootstrap 0.026, RSE 3%); Eq. 8 additive term
    e_sexf_vc        <- 0.89;       label("Female-sex multiplicative factor on Vc (power-form base)")         # Table 2 row 'SEX (V2)' = 0.89 (bootstrap 0.88, RSE 1%); Eq. 8 '0.89^Sex'
    e_dis_healthy_vp <- 1 / 1.23;   label("Healthy-subject multiplicative factor on Vp (power-form base)")    # Eq. 8 '1.23^DisStat' re-expressed onto DIS_HEALTHY = 1 - DisStat (Table 2 row 'DisStat (V3)' rounds this to 1.2)
    e_fed_ka         <- 0.44;       label("Fed-state multiplicative factor on ka (power-form base)")          # Table 2 row 'Food (Ka)' = 0.44 (bootstrap 0.42, RSE 8%); Eq. 8 '0.44^Food'
    e_fed_tlag       <- 1.6;        label("Fed-state multiplicative factor on T_lag (power-form base)")       # Table 2 row 'Food (TLAG)' = 1.6 (bootstrap 1.6, RSE 9%); Eq. 8 '1.6^Food'
    e_fed_fdepot     <- 0.88;       label("Fed-state multiplicative factor on F1 (power-form base)")          # Table 2 row 'Food(F)' = 0.88 (bootstrap 0.88, RSE 2%); Eq. 8 '0.88^Food'

    # Inter-individual variability. Chen 2023 reports the omega^2 variances directly
    # in Table 2 (rows 'omega^2 CL' etc.) on the log-normal exponential-IIV scale of
    # Eq. 1 'Para_ind = Para_pop x exp(eta + kappa)'. Results 'PPK model' confirms the
    # SD reading: 'the IIV (SD) and IOV (SD) of CL were 3.6% and 13%, and the IIV of
    # V2 and V3 was 5.7% and 8.8%, respectively. The IIV of Ka and Tlag were 81% and
    # 83%. IIV of bioavailability (F) was 19%.' -- e.g. sqrt(0.0013) = 3.6%,
    # sqrt(0.64) = 80%, sqrt(0.035) = 18.7%. No IIV was estimated on Q (no Table 2 row).
    etalcl     ~ 0.0013  # Table 2 'omega^2 CL'   = 0.0013 (shrinkage 75%) -> SD 3.6%
    etalvc     ~ 0.0034  # Table 2 'omega^2 V2'   = 0.0034 (shrinkage 71%) -> SD 5.8%
    etalvp     ~ 0.0080  # Table 2 'omega^2 V3'   = 0.0080 (shrinkage 64%) -> SD 8.9%
    etalka     ~ 0.64    # Table 2 'omega^2 KA'   = 0.64   (shrinkage 29%) -> SD 80%
    etaltlag   ~ 0.71    # Table 2 'omega^2 TLAG' = 0.71   (shrinkage 39%) -> SD 84%
    etalfdepot ~ 0.035   # Table 2 'omega^2 F'    = 0.035  (shrinkage 27%) -> SD 19%

    # Inter-occasion variability on log-CL (kappa in Eq. 1 / Eq. 8), single shared
    # pi^2 across the two occasions. nlmixr2 has no NONMEM `$OMEGA BLOCK(1) SAME`
    # shortcut, so occasion 2 gets its own eta with the variance fixed equal to the
    # occasion-1 estimate (Jonsson_2011_ethambutol.R precedent).
    etaiov_cl_1 ~ 0.017        # Table 2 'pi^2 CL' = 0.017 (shrinkage 45 or 20) -> IOV SD 13% (estimated)
    etaiov_cl_2 ~ fix(0.017)   # SAME-equivalent: equal to the occasion-1 IOV variance

    # Residual error. Eq. 2 is the mixed model 'Y = Ypred x (1 + eps1) + eps2', so the
    # Table 2 sigma^2 entries are variances and the nlmixr2 SDs are their square
    # roots. Table 2 note: the additive row carries unit 10^-4, i.e. sigma^2_add =
    # 0.34e-4.
    propSd <- sqrt(0.032);    label("Proportional residual error (fraction)") # Table 2 'sigma^2 prop' = 0.032 (shrinkage 16%) -> propSd = 0.179
    addSd  <- sqrt(0.34e-4);  label("Additive residual error (mg/L)")         # Table 2 'sigma^2 add'  = 0.34 x 10^-4 (shrinkage 14%) -> addSd = 0.00583
  })

  model({
    # Decompose the integer occasion column into binary indicators multiplexing the
    # IOV etas on log-CL. OCC = 1 is the single-dose occasion, OCC = 2 the
    # steady-state occasion (from 72 h of multiple dosing onward).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Individual PK parameters, Eq. 8. Body-weight allometry uses the 70 kg reference
    # with exponents fixed at 0.75 (CL, Q) and 1 (Vc, Vp). CL carries an additive
    # linear CLcr term (raw Cockcroft-Gault mL/min, no centering), and the log-normal
    # BSV + IOV multiply the combined intercept + renal slope, matching Eq. 8's
    # placement of exp(eta1 + kappa) outside the parenthesized CL expression.
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * (WT / 70)^e_wt_cl_q  * exp(etalcl + iov_cl)
    vc <- exp(lvc + etalvc)             * (WT / 70)^e_wt_vc_vp * e_sexf_vc^SEXF
    q  <- exp(lq)                       * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp)             * (WT / 70)^e_wt_vc_vp * e_dis_healthy_vp^DIS_HEALTHY

    ka     <- exp(lka     + etalka)     * e_fed_ka^FED
    tlag   <- exp(ltlag   + etaltlag)   * e_fed_tlag^FED
    fdepot <- exp(lfdepot + etalfdepot) * e_fed_fdepot^FED

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with first-order absorption from a depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Plasma nemonoxacin concentration in mg/L (dose in mg, volumes in L). Combined
    # proportional + additive residual error per Eq. 2.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
