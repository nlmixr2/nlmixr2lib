GallegoHernandez_2026_linezolid <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for linezolid in elderly hospitalized",
    "patients aged 65-87 years (Gallego-Hernandez 2026), developed from routine",
    "therapeutic-drug-monitoring records at a single Spanish tertiary hospital. Clearance",
    "(4.25 L/h typical) carries three centred power covariates: absolute CKD-EPI eGFR",
    "(exponent 0.29, referenced to the population median 59.56 mL/min), treatment duration in",
    "days (exponent -0.18, referenced to 3.5 days) and age (exponent -1.16, referenced to",
    "78 years). The treatment-duration term is the paper's novel finding: apparent clearance",
    "declines progressively over a course of therapy, which the authors describe as a",
    "phenomenological description of the observed data rather than evidence of a specific",
    "biological mechanism. Central volume is a single typical value of 25.6 L with no",
    "covariates and no interindividual variability, the predominantly trough-oriented sampling",
    "design having been unable to support either. Interindividual variability is carried on",
    "clearance alone (33.2% CV) and residual error is proportional (26.5% CV). The model is",
    "intended as a Bayesian TDM tool for detecting linezolid overexposure in the elderly, not",
    "as a fully descriptive structural model.",
    sep = " "
  )
  reference <- paste(
    "Gallego-Hernandez G, Albarran-Gomez A, Sanchez-Hernandez JG, Garcia-Casanueva JC,",
    "Otero MJ. Population Pharmacokinetics of Linezolid in Elderly Hospitalized Patients:",
    "Implications for Therapeutic Drug Monitoring.",
    "Pharmaceutics. 2026;18(5):528. doi:10.3390/pharmaceutics18050528",
    sep = " "
  )
  vignette <- "GallegoHernandez_2026_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "linezolid", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate, CKD-EPI equation, in ABSOLUTE mL/min (de-normalized, NOT mL/min/1.73 m^2).",
      units              = "mL/min (absolute, NOT BSA-normalized)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as the centred power term (CRCL / 59.56)^0.29 (Results 3.2 final-model",
        "equation). Three cautions carry forward. First, the normalisation: Methods 2.1 states",
        "that 'the CKD-EPI estimated glomerular filtration rate (eGFR) expressed in absolute",
        "values (mL/min) was used', so this column is NOT on the mL/min/1.73 m^2 scale that is",
        "this canonical's default, and the paper argues in Discussion that absolute eGFR is the",
        "physiologically coherent pairing for a drug clearance. Supplying a BSA-normalized value",
        "silently rescales the renal term. Second, the values were CAPPED AT 130 mL/min before",
        "fitting, because reduced muscle mass in the elderly depresses serum creatinine and so",
        "inflates any creatinine-based estimate (Methods 2.1); a user supplying uncapped data is",
        "extrapolating past what was fitted. Third, the column is TIME-VARYING: Methods 2.1",
        "specifies renal function was 'handled as a time-varying covariate, using the value",
        "closest to each concentration measurement'. Development-cohort median 48.8 mL/min",
        "(IQR 25.9-75.4, Table 1), so the 59.56 centring constant sits ABOVE the cohort median",
        "and the typical patient's clearance is below the printed 4.25 L/h. The centring constant",
        "appears as 59.56 in the Results 3.2 prose and as the rounded 59.6 inside the printed",
        "equation; the two agree to the precision shown and the full-precision form is used here.",
        "The paper tested Cockcroft-Gault on several weight descriptors as well and found",
        "absolute eGFR the strongest predictor."
      ),
      source_name        = "eGFR"
    ),
    AGE = list(
      description        = "Age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as the centred power term (AGE / 78)^-1.16 (Results 3.2 final-model",
        "equation), so clearance FALLS with advancing age. The 78-year centring constant is the",
        "development-cohort median (Table 1). Two cautions. First, the exponent is steep",
        "(-1.16) but rests on a narrow fitted range of 65-87 years and is the least precisely",
        "estimated fixed effect in the model (RSE 43%, bootstrap 95% CI -1.90 to -0.46); across",
        "the full observed age range the age multiplier spans only about 1.24 to 0.86, so the",
        "steepness is an artefact of extrapolating a power function whose argument barely moves",
        "away from 1. It must NOT be extrapolated below 65 years, where it diverges rapidly.",
        "Second, age is itself an input to the CKD-EPI equation that produces CRCL, so the two",
        "retained covariates are not mechanistically independent; the authors addressed this",
        "directly and report only weak collinearity (Spearman rho = -0.203, variance inflation",
        "factor 1.05), concluding age carries explanatory information beyond eGFR in this",
        "dataset. Treated as time-fixed."
      ),
      source_name        = "AGE"
    ),
    T_FIRSTDOSE = list(
      description        = "Time elapsed since the first linezolid dose of the treatment course, i.e. treatment duration.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The paper's DAY covariate, 'treatment duration in days'. Carried in canonical HOURS and",
        "divided by 24 inside model() to recover the paper's days, exactly as the T_FIRSTDOSE",
        "register entry prescribes and as Eechoute_2012_imatinib.R and",
        "Chen_2021_lorlatinib_hypertriglyceridemia.R do. Enters clearance as the centred power",
        "term (DAY / 3.5)^-0.18 (Results 3.2 final-model equation), so apparent clearance",
        "declines as therapy continues. THIS TERM IS UNDEFINED AT T_FIRSTDOSE = 0: a negative",
        "power of zero is infinite, so a simulation must supply a strictly positive value. That",
        "is a property of the published model, not an encoding choice, and it is harmless in",
        "practice because the fitted data contain no record at day 0 -- the local TDM protocol",
        "drew the first sample before the third to fifth dose, giving a median 3 days",
        "(IQR 2-4) to the first measurement (Table 1), and treatment duration ran to a median",
        "7 days (range 3-26). The paper's own simulations (Figures 4-6) hold treatment duration",
        "FIXED per scenario at 3, 7 and 10 days rather than letting it run as a clock, and the",
        "vignette reproduces them that way. Users fitting real TDM records should supply the",
        "per-record treatment duration, which makes the column time-varying and monotonically",
        "increasing. The authors caution that in a retrospective real-world dataset this",
        "covariate may capture concurrent changes in clinical status, inflammation, renal",
        "function, TDM-driven dose adaptation or survivor bias as much as any time-dependent",
        "pharmacokinetic process (Discussion)."
      ),
      source_name        = "DAY"
    )
  )

  # Covariates collected and screened in the PsN stepwise covariate analysis but
  # NOT retained in the final model (Methods 2.1 lists the collected variables;
  # Methods 2.3 describes the forward-inclusion / backward-elimination screen;
  # Table 2 shows the three retained effects). Documented here for provenance;
  # they are deliberately absent from model() and from covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.", units = "kg", type = "continuous",
      notes = "Development-cohort median 70 kg (IQR 60-80), Table 1. Screened and not retained. Ideal body weight (median 58.5 kg, IQR 54-64.3), adjusted body weight (0.4 correction factor), body mass index (median 26.6 kg/m^2, IQR 23.8-30.0) and body surface area (Mosteller) were also collected and screened; note that weight descriptors entered the screen chiefly as inputs to alternative Cockcroft-Gault renal-function estimators, all of which lost to absolute CKD-EPI eGFR (Results 3.2)."
    ),
    SEXF = list(
      description = "Female sex indicator.", units = "1 = female, 0 = male", type = "categorical",
      notes = "Development cohort 65.0% male, i.e. 35.0% female (Table 1). Screened and not retained."
    ),
    ALB = list(
      description = "Serum albumin.", units = "g/dL", type = "continuous",
      notes = "Development-cohort median 3.0 g/dL (IQR 2.7-3.3), Table 1. Screened and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine.", units = "mg/dL", type = "continuous",
      notes = "Collected (Methods 2.1) as the input to both the Cockcroft-Gault and CKD-EPI renal-function estimators. Not tabulated as a standalone row in Table 1 and not retained in its raw form; renal function entered the model through CRCL instead."
    ),
    ALT = list(
      description = "Alanine aminotransferase (hepatic-function marker).", units = "U/L", type = "continuous",
      notes = "Development-cohort median 24 U/L (IQR 11.5-36.5), Table 1. Screened and not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase (hepatic-function marker).", units = "U/L", type = "continuous",
      notes = "Development-cohort median 44 U/L (IQR 30-77), Table 1. Screened and not retained."
    ),
    TBILI = list(
      description = "Total bilirubin (hepatic-function marker).", units = "mg/dL", type = "continuous",
      notes = "Development-cohort median 0.4 mg/dL (IQR 0.3-0.6), Table 1, where it is abbreviated BLT. Screened and not retained."
    ),
    LDH = list(
      description = "Lactate dehydrogenase.", units = "U/L", type = "continuous",
      notes = "Development-cohort median 212.0 U/L (IQR 174.5-285.5), Table 1. Screened and not retained."
    ),
    TPRO = list(
      description = "Total serum protein.", units = "g/dL", type = "continuous",
      notes = "Development-cohort median 5.5 g/dL (IQR 5.1-6.0), Table 1. Screened and not retained."
    ),
    CRP = list(
      description = "C-reactive protein (inflammation marker).", units = "mg/L", type = "continuous",
      notes = "Development-cohort median 9.4 mg/L (IQR 4.0-19.2), Table 1. Screened and not retained."
    ),
    PROCALCITONIN = list(
      description = "Serum procalcitonin (inflammation / sepsis marker).", units = "ng/mL", type = "continuous",
      notes = "Development-cohort median 0.5 ng/mL (IQR 0.2-1.3), Table 1. Screened and not retained."
    ),
    HGB = list(
      description = "Haemoglobin.", units = "g/dL", type = "continuous",
      notes = "Development-cohort median 10.1 g/dL (IQR 9.0-11.6), Table 1. Collected as a linezolid haematological-toxicity marker as well as a screened covariate. Not retained."
    ),
    PLT = list(
      description = "Platelet count.", units = "10^9/L", type = "continuous",
      notes = "Development-cohort median 241 x10^9/L (IQR 172.5-349.5), Table 1. Collected chiefly as the linezolid thrombocytopenia marker. Not retained."
    ),
    CONMED_RIF = list(
      description = "Concomitant rifampicin indicator (a known inducer).", units = "1 = yes, 0 = no", type = "categorical",
      notes = "3 of 103 development-cohort patients (2.9%), Table 1. Methods 2.1 states that concomitant medications with potential pharmacokinetic interaction were recorded 'with particular attention to known enzyme or transporter inducers and inhibitors, including rifampicin and macrolides'. Not retained; at n = 3 the cohort carries almost no information about this interaction."
    ),
    CONMED_MACROLIDE = list(
      description = "Concomitant macrolide indicator.", units = "1 = yes, 0 = no", type = "categorical",
      notes = "7 of 103 development-cohort patients (6.8%), Table 1. Screened per Methods 2.1 and not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor indicator.", units = "1 = yes, 0 = no", type = "categorical",
      notes = "88 of 103 development-cohort patients (85.4%), Table 1. Recorded as part of the polypharmacy profile; not retained. Note the near-universal prevalence leaves little contrast to estimate an effect from."
    ),
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal indicator (CYP3A4 / P-gp inhibitor).", units = "1 = yes, 0 = no", type = "categorical",
      notes = "13 of 103 development-cohort patients (12.6%), Table 1. Screened per Methods 2.1 and not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 103,
    n_studies      = 1,
    age_median     = "78 years (range 65-87)",
    age_range      = "65-87 years",
    weight_median  = "70 kg (IQR 60-80)",
    bmi_median     = "26.6 kg/m^2 (IQR 23.8-30.0)",
    sex_female_pct = 35.0,
    race_ethnicity = "Single-centre Spanish cohort; race/ethnicity not reported in the source.",
    disease_state  = "Elderly hospitalized adults receiving intravenous linezolid as targeted or empirical therapy for Gram-positive infection. Diagnoses in the development cohort were respiratory infection 21.4%, skin and soft tissue infection 22.3%, urinary tract infection 23.3% and other infections 33.0% (Table 1). Patients on oral linezolid, with active oncological or haematological disease, on renal replacement therapy, or critically ill requiring ICU admission were EXCLUDED by design, so the model carries no information about those groups.",
    renal_function = "Broad and skewed towards impairment: absolute CKD-EPI eGFR median 48.8 mL/min (IQR 25.9-75.4). Strata: >90 mL/min 11.7%, 60-89 mL/min 26.2%, 30-59 mL/min 32.0%, <30 mL/min 30.1% (Table 1). Values were capped at 130 mL/min for modelling. Renal replacement therapy was an exclusion criterion.",
    dose_range     = "All patients started intravenous linezolid 600 mg every 12 h as a 1-h infusion, after which dosing was individualized on TDM results. Development-cohort daily dose median 1200 mg/day (range 300-1800), maximum daily dose median 1200 mg/day (range 600-2400), daily dose per body weight median 14.7 mg/kg/day (IQR 10.0-17.7). Treatment duration median 7 days (range 3-26).",
    regions        = "Spain (University Hospital of Salamanca).",
    notes          = paste(
      "Retrospective, single-centre study of routine TDM records, January 2024 to September",
      "2025. 149 patients contributing 293 quantifiable serum concentrations were randomly",
      "split about 2:1 into a development cohort of 103 patients / 198 concentrations (the fit",
      "this model reproduces) and an independent validation cohort of 46 patients / 95",
      "concentrations. A further 15 measurements below the 0.8 mg/L assay LLOQ were excluded",
      "rather than handled by an M-method. Sampling is predominantly TROUGH-ORIENTED and sparse",
      "-- median 2 samples per patient (range 1-5), the first drawn before the third to fifth",
      "dose -- which is the stated reason the model carries no peripheral compartment and no",
      "IIV on volume, and why the authors call it a pragmatic TDM tool rather than a fully",
      "descriptive structural model. Concentrations were measured by ENZYME IMMUNOASSAY",
      "(ARK Linezolid Assay on an Abbott Architect ci4100), not LC-MS/MS, with limited reported",
      "cross-reactivity against the inactive linezolid metabolites; the authors flag this as a",
      "contributor to residual variability. Independent-validation performance was a mean",
      "prediction error of -10.54% and a mean absolute prediction error of 47.8%, which the",
      "authors contextualise against an external evaluation in which 25 published linezolid",
      "models returned median MAPEs of 58-73%. Estimation was FOCE-I in NONMEM 7.5 with PsN",
      "5.3.1; the final model was confirmed by a 1000-replicate bootstrap (990 successful).",
      "Development-cohort mean observed concentration 6.5 mg/L (SD 4.0)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural typical values. Table 2 'Final Model / Estimate' column.
    # Linezolid is given intravenously here, so CL and V are absolute, not
    # apparent, and there is no absorption compartment and no F term.
    # ------------------------------------------------------------------
    lcl <- log(4.25); label("Clearance at the reference covariate values (L/h)")  # Table 2 CLpop = 4.25 L/h (RSE 9%; bootstrap mean 4.22, 95% CI 3.76-4.82). Reference patient: eGFR 59.56 mL/min, 3.5 days of treatment, age 78 years.
    lvc <- log(25.60); label("Central volume of distribution, V (L)")             # Table 2 Vpop = 25.60 L (RSE 17%; bootstrap mean 25.40, 95% CI 21.8-30.35). No covariates; Discussion compares it with the 27.6 L of Matsumoto 2014 in elderly Japanese patients.

    # ------------------------------------------------------------------
    # Covariate effects on clearance. All three are CENTRED POWER terms
    # multiplying CLpop, per the final-model equation printed in
    # Results 3.2:
    #   CL_i (L/h) = 4.25 * (eGFR_i/59.6)^0.29 * (DAY_i/3.5)^-0.18
    #                     * (AGE_i/78)^-1.16 * exp(eta_CL,i)
    # The centring constants are carried as literals in model() rather
    # than as estimated parameters because the paper fixed them to
    # population medians; see the model() block.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.29; label("Power exponent for absolute CKD-EPI eGFR on CL (unitless)")        # Table 2 'eGFR-CL' = 0.29 (RSE 16%; bootstrap mean 0.28, 95% CI 0.18-0.38)
    e_tfirstdose_cl <- -0.18; label("Power exponent for treatment duration on CL (unitless)")    # Table 2 'DAY-CL' = -0.18 (RSE 14%; bootstrap mean -0.18, 95% CI -0.25 to -0.12)
    e_age_cl <- -1.16; label("Power exponent for age on CL (unitless)")                          # Table 2 'AGE-CL' = -1.16 (RSE 43%; bootstrap mean -1.17, 95% CI -1.90 to -0.46). The least precisely estimated fixed effect in the model.

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential / log-normal on CL, per the
    # exp(eta_CL,i) factor of the Results 3.2 equation and the statement
    # that 'IIV on CL followed a log-normal distribution'.
    #
    # VARIANCE CONVENTION. Table 2 heads the row 'IIVCL (CV, %)' = 33.20
    # and the residual row 'RUVprop (CV, %)' = 26.50 with the SAME '(CV, %)'
    # label. For a proportional residual error the only sensible reading of
    # that header is sqrt(sigma^2) * 100, so the shared header pins the
    # omega row to the same convention: 33.20 is sqrt(omega^2) * 100, and
    # omega^2 = 0.332^2 = 0.110224. The exact log-normal alternative,
    # omega^2 = log(1 + 0.332^2) = 0.10455, differs by 2.6% on the SD and
    # is recorded in the vignette Errata. No supplementary control stream
    # is published for this paper, so the header is the only evidence.
    #
    # IIV on V is NOT carried: 'IIV on V was initially evaluated but was
    # not supported by the limited information content of the predominantly
    # late sampling design and was consequently excluded from the final
    # model' (Results 3.2).
    # ------------------------------------------------------------------
    etalcl ~ 0.110224  # Table 2 'IIVCL (CV, %)' = 33.20 -> 0.332^2 (RSE 9%, bootstrap 95% CI 28.3-36.10; eta-shrinkage 4.3%). Base structural model was 42.9% before the three covariates entered.

    # ------------------------------------------------------------------
    # Residual unexplained variability. Proportional: 'RUV was best
    # described by a proportional error model. Additive and combined error
    # structures were evaluated but did not improve OFV or diagnostic
    # performance' (Results 3.2).
    # ------------------------------------------------------------------
    propSd <- 0.265; label("Proportional residual SD for Cc (fraction)")  # Table 2 'RUVprop (CV, %)' = 26.50 -> 0.265 (RSE 17%, bootstrap 95% CI 24.50-30.00)
  })

  model({
    # 1. Covariate centring constants, all fixed by the authors to
    #    population medians rather than estimated (Results 3.2).
    #      - 59.56 mL/min is named in the Results 3.2 prose as the
    #        'population median value'; the printed equation shows the
    #        rounded 59.6. The two agree to the precision shown and the
    #        difference moves CL by 0.02% at the 0.29 exponent.
    #      - 3.5 days is the centring of the continuous time-dependent
    #        power model.
    #      - 78 years is the development-cohort median age (Table 1).
    crcl_ref <- 59.56  # mL/min, absolute CKD-EPI eGFR
    day_ref <- 3.5     # days of treatment
    age_ref <- 78      # years

    # 2. T_FIRSTDOSE is carried in canonical hours; the paper's DAY
    #    covariate is in days, so divide by 24. Undefined at
    #    T_FIRSTDOSE = 0 because the exponent is negative -- see
    #    covariateData$T_FIRSTDOSE$notes.
    day_trt <- T_FIRSTDOSE / 24

    # 3. Individual parameters. Results 3.2 final-model equation.
    cl <- exp(lcl + etalcl) *
      (CRCL / crcl_ref)^e_crcl_cl *
      (day_trt / day_ref)^e_tfirstdose_cl *
      (AGE / age_ref)^e_age_cl
    vc <- exp(lvc)

    # 4. Micro-constant
    kel <- cl / vc

    # 5. One compartment with first-order elimination. Linezolid is
    #    administered intravenously (1-h infusion), so there is no
    #    absorption compartment and bioavailability is structurally 1.
    #    A two-compartment structure was evaluated and rejected: it gave
    #    only dOFV ~11.7 and the peripheral volume was imprecise
    #    (V2 RSE 69%) under the trough-oriented design (Results 3.2).
    d/dt(central) <- -kel * central

    # 6. Observation. Dose in mg, volume in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
