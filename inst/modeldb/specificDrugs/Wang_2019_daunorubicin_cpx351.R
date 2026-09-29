Wang_2019_daunorubicin_cpx351 <- function() {
  description <- paste(
    "Two-compartment linear population PK model for TOTAL daunorubicin",
    "(encapsulated plus released) following a 90-minute intravenous",
    "infusion of the dual-drug liposome CPX-351 in adults with hematologic",
    "malignancies (Wang 2019; 2023 plasma concentrations from 195 patients",
    "pooled across the phase 1 study 101, the phase 2 study 206 and the",
    "phase 3 study 301). Because more than 99% of circulating drug remains",
    "inside the liposome, these parameters describe the CPX-351 liposome",
    "rather than free daunorubicin: clearance is roughly 880-fold lower and",
    "the central volume roughly 320-fold smaller than for non-liposomal",
    "daunorubicin (Table S4). Body surface area scales all four disposition",
    "parameters on a 1.95 m^2 reference, with estimated exponents on CL and",
    "Vc and exponents fixed at 1 on Q and Vp. Total bilirubin has a very",
    "shallow positive effect on clearance (exponent 0.0829 on a 0.60 mg/dL",
    "reference), and the frozen formulation used in study 101 lowers all",
    "four disposition parameters relative to the lyophilized reference;",
    "the authors judged neither effect clinically meaningful. The",
    "log-additive residual error itself carries between-subject variability",
    "on its magnitude."
  )
  reference <- paste(
    "Wang Q, Banerjee K, Vasilinin G, Marier JF, Gibbons JA (2019).",
    "Population Pharmacokinetics and Exposure-Response Analyses for CPX-351",
    "in Patients With Hematologic Malignancies. The Journal of Clinical",
    "Pharmacology 59(5):748-762. doi:10.1002/jcph.1366.",
    sep = " "
  )
  vignette <- "Wang_2019_cpx351"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "daunorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "daunorubicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling on a 1.95 m^2 reference applied to CL, Vc, Q and Vp,",
        "per the Table 2 footnote b: the tabulated population estimates are",
        "'calculated for a typical patient with BSA of 1.95 m^2 who received",
        "the lyophilized formulation'. Note 1.95 m^2 is NOT the cohort median,",
        "which Table 1 gives as 1.94 m^2 (range 1.26-2.80); it is the",
        "standardisation constant the authors chose. The exponents on CL",
        "(0.829) and Vc (1.12) were estimated; those on Q and Vp were fixed at",
        "1.00 (Table S3, SE column reads 'fixed'). Results: the CL, Vc and Vp",
        "of daunorubicin 'were mainly dependent on BSA and consistent with the",
        "current dosing paradigm, where CPX-351 is dosed based on mg/m2'. The",
        "BSA formula used is unspecified in the source. Body weight, body mass",
        "index and BSA were all screened; BSA was the size descriptor",
        "retained."
      ),
      source_name = "BSA"
    ),
    TBILI = list(
      description = "Baseline total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source reports bilirubin in US-convention mg/dL; the canonical",
        "register unit is SI umol/L, so model() converts inline with",
        "tbili_mgdL <- TBILI / 17.1 before applying the power term. Power",
        "scaling on a 0.60 mg/dL (10.3 umol/L) reference with an estimated",
        "exponent of 0.0829 on CL only (Table 2 and Table S3 row BILI_CL,",
        "SE 0.0524, RSE 63.2%, 95% CI -0.0198 to 0.186). Note the confidence",
        "interval spans zero: the effect is retained because the paper's full",
        "model performed no reduction step, not because it is well",
        "identified. Cohort median 0.60 mg/dL, range 0.1-2.5 (Table 1).",
        "Results: 'The effect of bilirubin on CL was very shallow, with an",
        "exponent of 0.0829. These results suggest that patients with higher",
        "bilirubin values were associated with a slightly faster CL of the",
        "liposome and/or daunorubicin, as daunorubicin is primarily eliminated",
        "through hepatic pathways.' Table 3 shows the resulting AUCtau",
        "difference between the < 1.2 and 1.2-3 mg/dL strata is smaller than",
        "the within-stratum %CV. The cohort contains no patient above 3 mg/dL,",
        "so the term must not be extrapolated beyond that."
      ),
      source_name = "Bilirubin"
    ),
    FORM_CPX351_FROZEN = list(
      description = "Frozen (versus lyophilized) CPX-351 formulation indicator",
      units = "unitless",
      type = "categorical",
      reference_category = "0 (lyophilized formulation; the typical-value reference)",
      notes = paste(
        "1 = the patient received the frozen CPX-351 formulation, 0 = the",
        "lyophilized formulation. Study-fixed: all 38 patients of study 101",
        "received the frozen formulation and all 157 patients of studies 206",
        "and 301 the lyophilized (Table 1), so the covariate is perfectly",
        "confounded with study and cannot be separated from a study effect.",
        "Retained on all four disposition parameters in the final daunorubicin",
        "model and on none in the companion cytarabine model. Table S3 gives",
        "the coefficients on the log scale as 'x exp(Form_X) if Frozen':",
        "Form_CL -0.275, Form_Vc -0.221, Form_Q -0.551, Form_Vp -0.686. Table",
        "2 prints the same effects already exponentiated -- exp(-0.275) =",
        "0.760, exp(-0.221) = 0.802, exp(-0.551) = 0.576, exp(-0.686) = 0.504",
        "-- and the log-scale form is used here because it is the estimated",
        "parameterisation. Discussion: 'Formulation had only minor effects on",
        "the pharmacokinetics of total daunorubicin, but not cytarabine (Table",
        "2), and thus is unlikely to be of clinical importance.'"
      ),
      source_name = "Formulation"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final
  # daunorubicin model. Documentation only: checkModelConventions() does not
  # require these to be referenced in model(). Results: 'dose and other
  # intrinsic covariates did not exert a significant effect on the
  # pharmacokinetic parameters of daunorubicin.'
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a size descriptor; BSA was retained instead. Cohort median 79.8 kg, range 38.9-156.5 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed among the screened intrinsic factors (Methods, 'Population Pharmacokinetic Modeling'); not retained. No summary statistics are tabulated."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 67 years, range 24-81 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "unitless",
      type = "categorical",
      notes = "Screened; not retained. 119 of 195 patients (61.0%) were male, i.e. 39.0% female (Table 1)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "unitless",
      type = "categorical",
      notes = "Race was screened as a single categorical factor; not retained. 165 of 195 patients (84.6%) were White (Table 1). The paper does not state which level was the reference."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "unitless",
      type = "categorical",
      notes = "Screened as part of the race factor; not retained. 9 of 195 patients (4.6%) (Table 1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "unitless",
      type = "categorical",
      notes = "Screened as part of the race factor; not retained. 10 of 195 patients (5.1%) (Table 1)."
    ),
    WBC = list(
      description = "Baseline white blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = paste(
        "Screened; not retained. Cohort median 3.4 x 10^9/L, range 0.2-110.9",
        "(Table 1). Discussion: 'The inclusion of WBC counts in the models",
        "was based on the observation by Krogh-Madsen et al that baseline WBC",
        "count was a significant covariate for cytarabine and daunorubicin",
        "pharmacokinetics.'"
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened both as a continuous covariate and as the categorical renal",
        "function strata (normal >= 90, mild 60-89, moderate 30-59, severe",
        "15-29 mL/min); not retained. Cohort median 85.3 mL/min, range",
        "27.5-211.7 (Table 1). Table 3 shows mean AUCtau 18% higher in",
        "moderate impairment than in normal function, a difference smaller",
        "than the 33.5-45.5% within-stratum %CV. Only one patient had severe",
        "impairment and none had end-stage disease, so the paper explicitly",
        "declines to extrapolate there."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 0.90 mg/dL, range 0.34-2.02 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a marker of hepatic function; not retained. Cohort median 23.0 U/L, range 5-115 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a marker of hepatic function; not retained. Cohort median 24.0 U/L, range 3-153 (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 195L,
    n_studies = 3L,
    age_median = "67 years",
    age_range = "24-81 years",
    weight_median = "79.8 kg",
    weight_range = "38.9-156.5 kg",
    bsa_median = "1.94 m^2",
    bsa_range = "1.26-2.80 m^2",
    sex_female_pct = 100 * (195 - 119) / 195,
    race_ethnicity = c(
      White = 100 * 165 / 195,
      Black = 100 * 9 / 195,
      Asian = 100 * 10 / 195,
      `Native American` = 100 * 1 / 195,
      Other = 100 * 10 / 195
    ),
    disease_state = paste(
      "Adults with hematologic malignancies. Study 101 (NCT00389428, phase 1",
      "dose escalation, n = 38) enrolled relapsed or refractory acute",
      "myeloid leukemia, acute lymphocytic leukemia or myelodysplastic",
      "syndrome; study 206 (NCT02238925, phase 2 QTc study, n = 26) enrolled",
      "acute myeloid leukemia, relapsed or refractory acute lymphocytic",
      "leukemia or myelodysplastic syndrome in adults 18-80 years; study 301",
      "(phase 3 randomized, n = 131) enrolled patients 60-75 years with newly",
      "diagnosed high-risk or secondary acute myeloid leukemia."
    ),
    renal_function = paste(
      "Creatinine clearance median 85.3 mL/min, range 27.5-211.7 (Table 1).",
      "At baseline 83 patients (42.6%) had normal renal function (>= 90",
      "mL/min), 83 (42.6%) mild impairment (60-89), 28 (14.4%) moderate",
      "impairment (30-59) and 1 (0.5%) severe impairment (15-29). No patient",
      "had end-stage renal disease."
    ),
    hepatic_function = paste(
      "Total bilirubin median 0.60 mg/dL, range 0.1-2.5 (Table 1); 179",
      "patients (91.8%) below 1.2 mg/dL and 16 (8.2%) between 1.2 and 3",
      "mg/dL. No patient exceeded 3 mg/dL. ALT median 24.0 U/L (3-153), AST",
      "median 23.0 U/L (5-115), alkaline phosphatase median 72.0 U/L",
      "(21-319)."
    ),
    dose_range = paste(
      "CPX-351 given as 90-minute intravenous infusions on days 1, 3 and 5 of",
      "an induction cycle. One unit of CPX-351 = 1 mg cytarabine + 0.44 mg",
      "daunorubicin. Study 101 escalated from 3 units/m^2 (1.3 mg/m^2",
      "daunorubicin) to 134 units/m^2 (59 mg/m^2 daunorubicin); studies 206",
      "and 301 used the recommended 100 units/m^2 induction dose (44 mg/m^2",
      "daunorubicin) with 65 units/m^2 (29 mg/m^2 daunorubicin) for",
      "consolidation. 157 of 195 patients (80.5%) received 100-101 units/m^2",
      "in the first induction cycle."
    ),
    formulation = paste(
      "38 patients (19.5%), all from study 101, received the frozen",
      "formulation; 157 (80.5%), all from studies 206 and 301, received the",
      "lyophilized formulation (Table 1)."
    ),
    notes = paste(
      "2176 plasma samples were assayed; 2033 entered the initial",
      "daunorubicin analysis and 2023 the final one after removing 10",
      "outliers with absolute conditional weighted residual > 4. 100 samples",
      "(4.6%) were below the limit of quantification and set to missing. PK",
      "data were collected during the first induction cycle only. Total",
      "(encapsulated plus released) daunorubicin was measured after a",
      "liposomal rupture step, by HPLC with tandem mass spectrometry. NONMEM",
      "7.3 with the mu-referencing SAEM/IMP methods and ITS pre-estimation;",
      "PsN 4.4.8; R 3.2.5. All covariates screened into the full model were",
      "retained without a reduction step, per Harrell. Baseline demographics",
      "from Table 1; final parameter estimates from Table 2 and Supplemental",
      "Table S3. Note the paper states two different assay limits: the",
      "Bioanalytical Methods section gives a range of detection of 1000 to",
      "100 000 ng/mL, while Results reports samples below a limit of",
      "quantification of 5 ng/mL. The two are not reconciled in the source."
    )
  )

  ini({
    # Structural parameters -- Wang 2019 Table 2 and Supplemental Table S3
    # (final daunorubicin model). Typical values are for a patient with
    # BSA = 1.95 m^2, baseline total bilirubin = 0.60 mg/dL and the
    # lyophilized formulation (Table 2 footnote b). Table 2 writes the full
    # covariate model inline as
    #   CL = 0.147 x (BSA/1.95)^0.829 x (Bilirubin/0.60)^0.0829
    #        x 0.760 if frozen formulation
    # and Table S3 decomposes the same expression into the rows used below.
    # Methods give exponential IIV: theta_i = theta x exp(eta_i).
    lcl <- log(0.147); label("Clearance at BSA = 1.95 m^2, bilirubin = 0.60 mg/dL, lyophilized formulation (L/h)") # Table S3 row 'CL (L/h)' = 0.147 (SE 0.00683, RSE 4.6%, 95% CI 0.134-0.161); Table 2
    lvc <- log(4.29); label("Central volume of distribution at BSA = 1.95 m^2, lyophilized formulation (L)") # Table S3 row 'Vc (L)' = 4.29 (SE 0.105, RSE 2.4%, 95% CI 4.08-4.49); Table 2
    lq <- log(0.0294); label("Intercompartmental clearance at BSA = 1.95 m^2, lyophilized formulation (L/h)") # Table S3 row 'Q (L/h)' = 0.0294 (SE 0.00326, RSE 11.1%, 95% CI 0.0230-0.0358); Table 2
    lvp <- log(0.593); label("Peripheral volume of distribution at BSA = 1.95 m^2, lyophilized formulation (L)") # Table S3 row 'Vp (L)' = 0.593 (SE 0.0936, RSE 15.8%, 95% CI 0.410-0.777); Table 2

    # BSA power exponents on a 1.95 m^2 reference. Those on CL and Vc were
    # estimated; those on Q and Vp carry 'fixed' in the Table S3 SE column and
    # a dash for RSE and CI. Results describe the base model as carrying 'an
    # allometric model of BSA on CL, Vc, Q, and Vp (estimated exponents)'.
    e_bsa_cl <- 0.829; label("Body surface area exponent on CL (unitless)") # Table S3 row 'BSA_CL' = 0.829 (SE 0.302, RSE 36.4%, 95% CI 0.237-1.42); Table 2
    e_bsa_vc <- 1.12; label("Body surface area exponent on Vc (unitless)") # Table S3 row 'BSA_Vc' = 1.12 (SE 0.162, RSE 14.5%, 95% CI 0.802-1.44); Table 2
    e_bsa_q <- fixed(1.00); label("Body surface area exponent on Q (unitless)") # Table S3 row 'BSA_Q' = 1.00, SE column reads 'fixed'
    e_bsa_vp <- fixed(1.00); label("Body surface area exponent on Vp (unitless)") # Table S3 row 'BSA_Vp' = 1.00, SE column reads 'fixed'

    # Total-bilirubin power effect on CL, referenced to 0.60 mg/dL. Its 95% CI
    # spans zero; it is retained because the paper ran a full model with no
    # reduction step, per Harrell.
    e_tbili_cl <- 0.0829; label("Total bilirubin exponent on CL (unitless)") # Table S3 row 'BILI_CL' = 0.0829 (SE 0.0524, RSE 63.2%, 95% CI -0.0198 to 0.186); Table 2

    # Frozen-formulation effects, on the log scale as Table S3 estimates them
    # ('x exp(Form_X) if Frozen'). Exponentiating reproduces the multipliers
    # Table 2 prints: exp(-0.275) = 0.760 on CL, exp(-0.221) = 0.802 on Vc,
    # exp(-0.551) = 0.576 on Q, exp(-0.686) = 0.504 on Vp.
    e_form_cpx351_frozen_cl <- -0.275; label("Log-scale frozen-formulation effect on CL (unitless)") # Table S3 row 'Form_CL' = -0.275 (SE 0.107, RSE 39.0%, 95% CI -0.484 to -0.0649)
    e_form_cpx351_frozen_vc <- -0.221; label("Log-scale frozen-formulation effect on Vc (unitless)") # Table S3 row 'Form_Vc' = -0.221 (SE 0.057, RSE 25.7%, 95% CI -0.333 to -0.110)
    e_form_cpx351_frozen_q <- -0.551; label("Log-scale frozen-formulation effect on Q (unitless)") # Table S3 row 'Form_Q' = -0.551 (SE 0.285, RSE 51.9%, 95% CI -1.11 to 0.00902)
    e_form_cpx351_frozen_vp <- -0.686; label("Log-scale frozen-formulation effect on Vp (unitless)") # Table S3 row 'Form_Vp' = -0.686 (SE 0.173, RSE 25.2%, 95% CI -1.03 to -0.347)

    # Between-subject variability. Table S3 reports these as omega^2 with the
    # matching %CV in parentheses, and its footnote states the conversion
    # 'omega^2 % were calculated as sqrt(exp(omega^2)-1)'. The values below are
    # therefore variances, not CVs: sqrt(exp(0.218)-1) = 0.493 reproduces the
    # 49.3% printed in Table 2, sqrt(exp(0.0547)-1) = 0.237 reproduces 23.7%,
    # sqrt(exp(0.392)-1) = 0.693 reproduces 69.3%, sqrt(exp(0.563)-1) = 0.869
    # reproduces 86.9%, and sqrt(exp(0.720)-1) = 1.027 reproduces 102.7%.
    #
    # The paper fitted an OMEGA block on CL, Vc, Q and Vp but reports only the
    # diagonal; the off-diagonal covariances are not published, so these are
    # encoded as independent etas.
    etalcl ~ 0.218 # Table S3 row 'BSV_CL' = 0.218 (49.3%), RSE 11.3%, shrinkage 1.2%
    etalvc ~ 0.0547 # Table S3 row 'BSV_Vc' = 0.0547 (23.7%), RSE 13.0%, shrinkage 5.4%
    etalq ~ 0.392 # Table S3 row 'BSV_Q' = 0.392 (69.3%), RSE 21.1%, shrinkage 34.4%
    etalvp ~ 0.563 # Table S3 row 'BSV_Vp' = 0.563 (86.9%), RSE 21.0%, shrinkage 27.2%
    etaexpSd ~ 0.720 # Table S3 row 'BSV_LogErr' = 0.720 (102.7%), RSE 10.2%, shrinkage 4.4%

    # Residual error: additive on log-transformed concentrations, which
    # nlmixr2 writes as ~ lnorm(). Table S3 gives the equation verbatim:
    #   ln(Cobs) = ln(Cpred) + LogErr x exp(etaLogErr)
    # so the residual magnitude itself carries exponential between-subject
    # variability (etaexpSd above). Methods: 'Between-subject variability was
    # also included on the residual unknown error model to reflect the fact
    # that data from different individuals may have different information
    # content, considering the very wide range of dose and concentrations
    # observed in the current analysis.'
    expSd <- 0.143; label("Log-additive residual error SD, typical subject (log-scale SD)") # Table S3 row 'LogErr' = 0.143 (SE 0.0131, RSE 9.2%, 95% CI 0.117-0.169, shrinkage 11.7%); Table 2 'Error model' 14.3%
  })

  model({
    # The canonical register carries TBILI in SI umol/L; Wang 2019 reports
    # bilirubin in US-convention mg/dL and the 0.60 reference below is on that
    # scale, so convert before applying the power term.
    tbili_mgdL <- TBILI / 17.1 # SI umol/L -> US-convention mg/dL (1 mg/dL = 17.1 umol/L)

    # Table 2, daunorubicin column, written out term by term. The formulation
    # term is an indicator-scaled log-scale shift, so it is exactly 1 for the
    # lyophilized reference (FORM_CPX351_FROZEN = 0).
    cl <- exp(lcl + etalcl + e_form_cpx351_frozen_cl * FORM_CPX351_FROZEN) *
      (BSA / 1.95)^e_bsa_cl * (tbili_mgdL / 0.60)^e_tbili_cl
    vc <- exp(lvc + etalvc + e_form_cpx351_frozen_vc * FORM_CPX351_FROZEN) *
      (BSA / 1.95)^e_bsa_vc
    q <- exp(lq + etalq + e_form_cpx351_frozen_q * FORM_CPX351_FROZEN) *
      (BSA / 1.95)^e_bsa_q
    vp <- exp(lvp + etalvp + e_form_cpx351_frozen_vp * FORM_CPX351_FROZEN) *
      (BSA / 1.95)^e_bsa_vp

    # Per-subject residual error magnitude: the log-scale residual SD is
    # expSd x exp(etaexpSd) (Table S3 equation column).
    expSdi <- expSd * exp(etaexpSd)

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # Amounts in mg and volumes in L give mg/L, which is ug/mL -- the units of
    # the Table 3 AUCtau column (ug*h/mL).
    Cc <- central / vc

    Cc ~ lnorm(expSdi)
  })
}
