GarciaHervalejo_2026_teicoplanin <- function() {
  description <- paste(
    "One-compartment IV-infusion population PK model for teicoplanin in 100",
    "adult patients with haematological malignancies treated for febrile",
    "neutropenia at a single Spanish centre (Garcia-Hervalejo 2026).",
    "Clearance carries three covariates, all printed in the paper's Equation 1:",
    "CL (L/h) = 1.28 * (1 + 0.012 * (AGE - 62)) * (eGFR / 92.15)^0.35 *",
    "(IBW / 61)^3.2. IDEAL body weight -- not total, adjusted, BMI or BSA -- was",
    "the anthropometric descriptor that best explained clearance variability,",
    "and its exponent of 3.2 is far steeper than any allometric value: clearance",
    "spans 0.33 to 3.05 L/h across IBW 40 to 80 kg at otherwise typical",
    "covariates, a 9-fold range that dominates the model and drives every dosing",
    "conclusion the paper reaches. Between-subject variability is exponential on",
    "CL and Vc; residual error is purely additive at 2.61 mg/L, the paper having",
    "tested and rejected proportional and combined forms because the observed",
    "concentrations are almost all troughs over a narrow range. The dataset is",
    "trough-only therapeutic drug monitoring, so a two-compartment model was",
    "evaluated but not identifiable (dOFV about 0.2, distribution-parameter RSE",
    "over 100%) and the parsimonious one-compartment structure was retained.",
    "The paper's Monte Carlo simulations conclude that conventional regimens",
    "(6 mg/kg or 600 mg q12h x 3 then once daily) rarely reach a trough of",
    "15-20 mg/L, and propose an intensified five-dose 12 mg/kg q12h loading",
    "phase followed by 12 mg/kg once daily."
  )
  reference <- paste(
    "Garcia-Hervalejo M, Sanchez-Hernandez JG, Conde-Gonzalez I,",
    "Avendano Pita A, Otero MJ. Population Pharmacokinetics and",
    "Model-Informed Dose Optimization of Teicoplanin in Adults with",
    "Hematological Malignancies.",
    "Pharmaceutics. 2026;18(1):100.",
    "doi:10.3390/pharmaceutics18010100"
  )
  vignette <- "GarciaHervalejo_2026_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "teicoplanin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    IBW = list(
      description        = "Ideal body weight by the Devine formula; the anthropometric descriptor retained on clearance in preference to total body weight, adjusted body weight, BMI and BSA",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column IBW. Enters clearance as the POWER term printed in",
        "Garcia-Hervalejo 2026 Equation 1, normalised to the development-cohort",
        "median of 61 kg: (IBW / 61)^3.2. Table 1 gives the development-cohort",
        "median IBW as 61.4 kg (range 46.25-77); the equation's normalising",
        "constant is the rounded 61, and the Results text states the term as",
        "'(IBW/61)^3.2 ... normalized to the median population value'.",
        "THE EXPONENT 3.2 IS THE DEFINING FEATURE OF THIS MODEL AND IS NOT AN",
        "ALLOMETRIC EXPONENT. It is roughly four times the classical 0.75, so",
        "clearance is extraordinarily sensitive to IBW: relative to the 61 kg",
        "reference the multiplier is 0.26 at IBW 40 kg, 0.95 at 60 kg, 1.42 at",
        "70 kg and 2.38 at 80 kg. The estimate is nevertheless well identified --",
        "RSE 21%, bootstrap mean 3.25, 95% CI 2.13-4.78 (Table 2) -- so the CI",
        "excludes 1 comfortably and cannot be reconciled with any allometric",
        "reading. It is also the term that produces every dosing conclusion in",
        "the paper: the Figure 4 findings that IBW <= 40 kg patients are",
        "predicted to exceed the 60 mg/L safety threshold at steady state while",
        "IBW >= 80 kg patients cannot maintain 15 mg/L on once-daily dosing both",
        "follow directly from this exponent and are reproduced in the validation",
        "vignette.",
        "PRACTICAL CONSEQUENCE FOR REUSE: because the covariate is so steep,",
        "supplying a DIFFERENT ideal-body-weight formula than the source's, or",
        "accidentally supplying total body weight, biases clearance severely. At",
        "the cohort medians total body weight is 68 kg against an IBW of 61.4 kg,",
        "a ratio of only 1.11, yet substituting it inflates clearance by",
        "1.11^3.2 = 1.39, i.e. 39%. Always supply Devine ideal body weight.",
        "The paper prints the Devine formula with a typographical error, writing",
        "'IBW (kg) = 50.0 kg + 2.3 kg x ((height in cm - 152.4) - 2.54)' for men",
        "and the same form with 45.5 for women. The final operator must be a",
        "DIVISION, not a subtraction: (height_cm - 152.4) / 2.54 converts the",
        "excess height over 60 inches from centimetres to inches, which is the",
        "published Devine 1974 formula. The printed subtraction is arithmetically",
        "impossible -- it would give a 170 cm man an IBW of 84.6 kg, whereas the",
        "correct division gives 65.9 kg, and only the latter is compatible with a",
        "cohort whose median total body weight is 68 kg and median IBW 61.4 kg.",
        "Use men IBW = 50.0 + 2.3 * (height_cm - 152.4) / 2.54 and women",
        "IBW = 45.5 + 2.3 * (height_cm - 152.4) / 2.54.",
        "TIME-FIXED at baseline. The observed range in the development cohort is",
        "46.25-77 kg (Table 1); the paper's own simulations extrapolate to IBW 80",
        "kg and it flags this explicitly -- 'As patients with IBW >= 80 kg were",
        "not represented in the model development dataset, these findings rely on",
        "model extrapolation and should be interpreted with caution'. Given the",
        "exponent, extrapolation beyond the fitted range is unusually hazardous."
      ),
      source_name        = "IBW"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate by the CKD-EPI equation, BSA-normalised; power effect on clearance",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. CKD-EPI-estimated glomerular filtration rate",
        "(Methods, 'Study Design and Population': 'Renal function was estimated",
        "using the CKD-EPI (Chronic Kidney Disease Epidemiology Collaboration)",
        "equation to calculate the estimated glomerular filtration rate').",
        "Enters clearance as the POWER term printed in Equation 1, normalised to",
        "the development-cohort median: (eGFR / 92.15)^0.35.",
        "UNITS: BSA-normalised mL/min/1.73 m^2. The paper is internally",
        "inconsistent on this point -- Table 1 heads the row 'Estimated",
        "glomerular filtration rate *, mL/min/1.73m2' with a median of 92.15,",
        "while the Table 2 abbreviation footnote writes 'eGFR, estimated",
        "glomerular filtration rate (CKD-EPI formula) in mL/min.' Table 1",
        "governs: CKD-EPI is by construction a BSA-normalised equation and",
        "returns mL/min/1.73 m^2 unless explicitly de-normalised, which the",
        "Methods never mention. The same header-versus-Methods hazard is recorded",
        "against Taylor_2026_methotrexate.R in this register. In practice the",
        "distinction is partly self-cancelling here because the covariate is a",
        "ratio to the cohort median, but a user supplying a raw (de-normalised)",
        "mL/min value for a patient whose BSA is far from 1.73 m^2 will bias the",
        "renal term; supply the BSA-normalised value.",
        "The exponent 0.35 is weak, so the multiplier spans only 0.72 at an eGFR",
        "of 26 to 1.16 at 141 mL/min/1.73 m^2 -- the extremes of the development",
        "cohort. It is also the least precisely estimated fixed effect in the",
        "model (RSE 47%, bootstrap 95% CI 0.06-0.67), which is worth carrying",
        "forward: the CI's lower bound is close to zero, so the renal term is",
        "directionally supported but poorly pinned. The IBW term above dominates",
        "clearance by a wide margin.",
        "TIME-FIXED at baseline in this encoding. The Methods state that",
        "'Time-dependent factors, such as fluctuations in renal function during",
        "treatment, were also explored through additional analyses in NONMEM',",
        "but the paper reports no time-varying renal term in the final model and",
        "Equation 1 carries a single eGFR value per subject.",
        "Patients with end-stage renal disease or on renal replacement therapy",
        "were EXCLUDED by design (Methods), so the model carries no information",
        "about dialysis and the observed range is 26-141 mL/min/1.73 m^2."
      ),
      source_name        = "eGFR"
    ),
    AGE = list(
      description        = "Subject age; centred linear effect on clearance",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column AGE. Enters clearance as the median-CENTRED LINEAR",
        "multiplier printed in Equation 1: (1 + 0.012 * (AGE - 62)), with 62",
        "years the development-cohort median age (Table 1, range 18-85).",
        "Note the form: this is a linear multiplier on the raw deviation in",
        "years, NOT a power term and NOT normalised by the median. The",
        "coefficient therefore carries units of per year.",
        "The effect is positive -- clearance RISES with age, which is the",
        "opposite of the usual direction and which the Discussion addresses",
        "directly: 'Although age remained in the model as a covariate with a",
        "positive coefficient on CL, this does not necessarily reflect a true",
        "physiological increase in renal elimination. More plausibly, it",
        "represents an apparent increase in total CL driven by lower albumin",
        "concentrations and altered body composition in older or cachectic",
        "patients.' Teicoplanin is 90-95% protein bound, so hypoalbuminaemia",
        "raises the unbound fraction and hence apparent total clearance. Albumin",
        "itself did not reach significance in this cohort (see",
        "covariatesDataExcluded[[ALB]]), and the paper attributes that partly to",
        "collinearity with age -- so this AGE term is best read as a surrogate",
        "for the albumin/body-composition axis rather than a renal-ageing effect.",
        "Over the observed 18-85 year range the multiplier spans 0.472 to 1.276,",
        "so the term is not negligible, though it is far smaller than the IBW",
        "effect. Note that the multiplier would turn NEGATIVE below an age of",
        "62 - 1/0.012 = 8.7 years; the model is for adults (>= 18 years) only and",
        "must not be applied paediatrically.",
        "VALUE OF THE COEFFICIENT: 0.012 is taken from the printed Equation 1.",
        "Table 2 reports the same parameter as 0.01 and the Results prose",
        "restates it as 0.01, but Table 2 formats every estimate to two decimal",
        "places, so 0.012 and 0.010 are indistinguishable there. The bootstrap",
        "column discriminates between them: it reports a 95% CI of 0.01-0.02",
        "alongside an RSE of 16%. An estimate of 0.012 with RSE 16% implies a CI",
        "of about 0.0082-0.0158, which rounds to exactly the printed 0.01-0.02;",
        "an estimate of 0.010 with the same RSE implies 0.0069-0.0131, which",
        "would have printed as 0.01-0.01. The printed equation and the bootstrap",
        "interval therefore agree on 0.012, and standing policy prefers the",
        "printed equation over prose in any case. The practical difference is",
        "small (at age 85 the multiplier is 1.276 versus 1.230, a 3.7% shift in",
        "clearance).",
        "TIME-FIXED at baseline."
      ),
      source_name        = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 68 kg (range 41-127) in the development cohort (Table 1). Screened as one of five anthropometric descriptors on CL and not retained: 'Among the anthropometric descriptors evaluated (total body weight, adjusted body weight, IBW, BMI, and BSA), IBW was the descriptor that best explained interindividual variability in CL' (Results). The Discussion adds that total body weight 'may overestimate the metabolically active mass contributing to drug clearance' in patients with cachexia, fluid overload or malnutrition. Per-descriptor statistics are in Supplementary Table S1, which is not on disk."
    ),
    ABW = list(
      description = "Adjusted body weight, computed with a 0.4 correction factor applied to the excess of total over ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Collected as a demographic variable (Methods: 'adjusted body weight (using a 0.4 correction factor for excess body weight)') and screened as an anthropometric descriptor on CL, losing to IBW (Results). No cohort summary is tabulated. The name ABW is NOT registered in inst/references/covariate-columns.md because the model does not use it; it is documented here only to preserve the provenance of the paper's covariate screen."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Collected as weight/height^2 (Methods) and screened as an anthropometric descriptor on CL, losing to IBW (Results). No cohort summary is tabulated."
    ),
    BSA = list(
      description = "Body surface area by the Du Bois formula",
      units       = "m^2",
      type        = "continuous",
      notes       = "Collected by the Du Bois formula (Methods) and screened as an anthropometric descriptor on CL, losing to IBW (Results). No cohort summary is tabulated."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Median 3.3 g/dL (range 2.4-4.3) in the development cohort (Table 1) --",
        "uniformly low, the normal adult range being roughly 3.5-5.0 g/dL.",
        "Screened on CL and NOT retained, which the Discussion treats as a",
        "notable negative result because albumin is a retained covariate in",
        "several other teicoplanin models: 'In our cohort, this variable did not",
        "reach statistical significance, likely due to consistently low albumin",
        "levels, its partial collinearity with age, and the use of total rather",
        "than unbound STC, which may limit the ability to detect its true",
        "contribution.' The AGE term retained in the final model is best read as",
        "a partial surrogate for this axis (see covariateData[[AGE]]$notes)."
      )
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Median 5.4 g/dL (range 3.7-8.7) in the development cohort (Table 1). Screened as a biochemical covariate and not retained; per-covariate statistics are in Supplementary Table S1, which is not on disk."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Median 0.71 mg/dL (range 0.30-2.11) in the development cohort (Table 1). Collected as the input to the CKD-EPI equation; renal function entered the final model through the derived eGFR (canonical CRCL) rather than through raw creatinine."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Median 7.82 mg/dL in the development cohort (Table 1; the printed range '7.82-43.54' shares its lower bound with the median, which appears to be a transcription slip in the source table -- the validation cohort's range is given as 0.13-47.84). Screened as a biochemical covariate and not retained."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Median 9.2 g/dL (range 6.9-13.3) in the development cohort (Table 1), reflecting the anaemia typical of this population. Screened as a biochemical covariate and not retained."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "cells/uL",
      type        = "continuous",
      notes       = "Median 36,500/uL (range 3000-565,000) in the development cohort (Table 1), reflecting the profound thrombocytopenia of neutropenic haematology patients. Screened as a biochemical covariate and not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "40% of the development cohort were female (Table 1). Sex is not a covariate in the final model, but it is an input to the Devine ideal-body-weight formula that produces IBW (intercept 50.0 kg for men, 45.5 kg for women), so it acts on clearance indirectly through that derivation. Categorical covariates were screened by ANOVA on the eta values (Methods)."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Collected as a demographic variable (Methods). No cohort summary is tabulated. Not a covariate in the final model, but it is the other input to the Devine ideal-body-weight formula, so a user starting from raw demographics needs HT and SEXF to derive the IBW column this model consumes."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 100L,
    n_studies        = 1L,
    n_concentrations = 168L,
    age_median       = "62 years",
    age_range        = "18-85 years (development cohort, Table 1); the study enrolled adults aged >= 18 years",
    weight_median    = "68 kg total body weight; 61.4 kg ideal body weight",
    weight_range     = "41-127 kg total body weight; 46.25-77 kg ideal body weight (development cohort, Table 1)",
    sex_female_pct   = 40,
    race_ethnicity   = "Not reported by category; single-centre Spanish cohort (University Hospital of Salamanca)",
    disease_state    = paste(
      "Hospitalised adults (>= 18 years) with malignant haematological",
      "neoplasms receiving intravenous teicoplanin as part of treatment for",
      "febrile neutropenia. Diagnoses in the development cohort (Table 1):",
      "non-Hodgkin's lymphoma 20%, acute myeloblastic leukaemia 20%, multiple",
      "myeloma 20%, acute lymphoblastic leukaemia 9%, other haematological",
      "malignancies 31%. The population is characterised by hypoalbuminaemia",
      "(median albumin 3.3 g/dL), anaemia (median haemoglobin 9.2 g/dL) and",
      "profound thrombocytopenia (median platelets 36,500/uL). Patients with",
      "end-stage renal disease or on renal replacement therapy were EXCLUDED,",
      "as were those with incomplete data or no valid concentration",
      "measurement."
    ),
    dose_range       = paste(
      "Teicoplanin by 30-minute intravenous infusion. All patients started on",
      "600 mg every 12 h, with subsequent doses individualised by therapeutic",
      "drug monitoring (Methods). Realised exposure in the development cohort",
      "(Table 1): median 9.1 mg/kg overall (range 3-14.6), 8.8 mg/kg in the",
      "loading phase (6.3-14.6) and 9.2 mg/kg in the maintenance phase",
      "(3-14.5). NOTE that these mg/kg figures are per kilogram of TOTAL body",
      "weight, not ideal: 600 mg / 68 kg median total body weight = 8.8 mg/kg,",
      "matching the printed loading-phase median exactly. The published Monte",
      "Carlo simulations compare 6 mg/kg q12h x 3 then 6 mg/kg q24h (summary of",
      "product characteristics), 600 mg q12h x 3 then 600 mg q24h (local",
      "protocol), and the proposed intensified 12 mg/kg q12h x 5 then 12 mg/kg",
      "q24h."
    ),
    regions          = "Spain (University Hospital of Salamanca; patients treated February 2021 to December 2023)",
    renal_function   = paste(
      "Preserved to augmented. Development cohort CKD-EPI eGFR median 92.15",
      "mL/min/1.73 m^2 (range 26-141) and serum creatinine median 0.71 mg/dL",
      "(range 0.30-2.11). Patients with end-stage renal disease or on renal",
      "replacement therapy were excluded, so the model carries no information",
      "about dialysis. The Introduction highlights augmented renal clearance as",
      "'particularly common in younger patients with good functional status who",
      "receive cytotoxic chemotherapy or intensive supportive care', and the",
      "simulations identify preserved renal function as one of the two drivers",
      "of target-attainment failure (the other being high ideal body weight)."
    ),
    screened_covariates = paste(
      "All collected demographic, clinical and biochemical variables were",
      "screened (Methods): age, sex, total body weight, height, ideal body",
      "weight, adjusted body weight, BMI, BSA, haematological diagnosis, dosing",
      "regimen, treatment duration, serum albumin, total protein, serum",
      "creatinine, C-reactive protein, haemoglobin, platelet count and CKD-EPI",
      "eGFR. Screening combined physiological plausibility, visual inspection of",
      "eta-versus-covariate relationships, stepwise linear regression for",
      "continuous covariates and ANOVA for categorical ones, with retention",
      "requiring p < 0.05 and r^2 > 0.10 on at least one PK parameter; survivors",
      "went through NONMEM stepwise covariate modelling with forward inclusion",
      "at dOFV > 3.84 and backward elimination at dOFV > 6.63. Only IBW, eGFR",
      "and AGE survived, all on CL; no covariate was retained on Vd. Their",
      "combined inclusion reduced IIV on CL from 47.8% to 34.1% and on Vd from",
      "35.8% to 31.0% (Results; the Discussion restates the CL figure as 47.7%).",
      "The per-covariate univariate statistics are in Supplementary Table S1,",
      "which is not on disk -- see the vignette Errata."
    ),
    external_validation = paste(
      "An independent cohort of the 51 patients recruited after the development",
      "cohort closed, contributing 95 concentrations (median age 58 years,",
      "median total body weight 70 kg, median eGFR 98.77 mL/min/1.73 m^2,",
      "49.01% female). Predictive performance was assessed with goodness-of-fit",
      "plots and with the mean prediction error (bias) and mean absolute",
      "prediction error (precision); the paper reports 'acceptable accuracy and",
      "precision and no evidence of systematic bias' but does NOT print numeric",
      "MPE or MAPE values."
    ),
    notes            = paste(
      "Retrospective single-centre study; 151 patients and 263 serum",
      "concentrations in total, split 100 patients / 168 concentrations for",
      "model development and 51 patients / 95 concentrations for external",
      "validation. SAMPLING IS TROUGH-DOMINATED, which is the single most",
      "important limitation for reuse: the first sample was drawn 24 h after",
      "treatment initiation and subsequent samples every 72-96 h, 'always",
      "immediately before dose administration'. 60% of development-cohort",
      "samples came from the loading phase and 40% from maintenance; mean",
      "concentration 14.37 mg/L (SD 7.20). Because no sample was drawn during",
      "the distribution phase, a two-compartment model was not identifiable",
      "(dOFV about 0.2, distribution-parameter RSE > 100%) and the",
      "one-compartment structure was retained deliberately. The model should",
      "therefore not be used to predict peak or early post-infusion",
      "concentrations, only troughs and cumulative exposure.",
      "Assay: QMS Teicoplanin turbidimetric immunoassay (Thermo Fisher) on an",
      "Abbott Architect ci4100; LLOQ 3.0 mg/L, calibration range 3.0-50 mg/L,",
      "intra- and inter-day precision below 10%. All measured concentrations",
      "were above the LLOQ, so the dataset has no censored observations. Note",
      "that the assay measures TOTAL teicoplanin; the paper lists the absence of",
      "an unbound-fraction measurement as a limitation, teicoplanin being about",
      "90-95% protein bound.",
      "Estimation: NONMEM 7.5 with FOCE-I; diagnostics and simulation in R",
      "4.5.0. Qualification: goodness-of-fit plots, a 1000-replicate",
      "non-parametric bootstrap in PsN 4.9 (994 successful runs; 6 discarded for",
      "sitting close to a parameter-space boundary), a prediction-corrected VPC",
      "from 1000 simulations over 48 h, and the external validation above.",
      "Eta-shrinkage was 23.9% on CL and 24.7% on Vd."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects. Garcia-Hervalejo 2026 Table 2, 'Final Model /
    # Estimate' column, and the typical-value equation printed as Equation 1
    # in Results section 3.2:
    #
    #   CL (L/h) = 1.28 * (1 + 0.012 * (AGE - 62))
    #                   * (eGFR / 92.15)^0.35
    #                   * (IBW / 61)^3.2
    #   Vd (L)   = 92.10
    #
    # lcl is therefore the clearance of the REFERENCE subject: age 62 years,
    # eGFR 92.15 mL/min/1.73 m^2 and IBW 61 kg -- all three the
    # development-cohort medians (Table 1). No covariate was retained on Vd.
    #
    # Dose in mg with volume in L gives central / vc in mg/L, matching the
    # units$concentration declared above and the paper's mg/L targets.
    # ------------------------------------------------------------------
    lcl <- log(1.28);  label("Clearance in the reference subject, age 62 y / eGFR 92.15 mL/min/1.73 m^2 / IBW 61 kg (L/h)")  # Garcia-Hervalejo 2026 Table 2 row 'CLpop (L/h)' = 1.28 (RSE 7%; bootstrap mean 1.26, 95% CI 1.11-1.42), and the leading constant of Equation 1
    lvc <- log(92.10); label("Volume of distribution (L)")                                                                    # Garcia-Hervalejo 2026 Table 2 row 'Vpop (L)' = 92.10 (RSE 5%; bootstrap mean 92.36, 95% CI 85.25-100.50). The Results text rounds it to 92.1 L and notes it is 'consistent with those reported in the literature for teicoplanin in adult patients, reflecting its high plasma protein binding and extensive tissue distribution'

    # ------------------------------------------------------------------
    # Covariate effects on clearance. All three are printed in Equation 1 and
    # tabulated in Table 2; applied in model() exactly as printed.
    #
    # e_age_cl is a per-YEAR slope inside a centred LINEAR multiplier, not an
    # exponent -- note the different functional form from the two power terms.
    # Its value is 0.012 from Equation 1 rather than the 0.01 printed in
    # Table 2 and the Results prose; see covariateData[[AGE]]$notes for the
    # bootstrap-interval arithmetic that discriminates between the two, and the
    # vignette Errata.
    # ------------------------------------------------------------------
    e_ibw_cl  <- 3.20;  label("Power exponent of ideal body weight on CL (unitless)")            # Garcia-Hervalejo 2026 Table 2 row 'IBW-CL' = 3.20 (RSE 21%; bootstrap mean 3.25, 95% CI 2.13-4.78), and the (IBW/61) exponent of Equation 1
    e_crcl_cl <- 0.35;  label("Power exponent of CKD-EPI eGFR on CL (unitless)")                 # Garcia-Hervalejo 2026 Table 2 row 'eGFR-CL' = 0.35 (RSE 47%; bootstrap mean 0.33, 95% CI 0.06-0.67), and the (eGFR/92.15) exponent of Equation 1
    e_age_cl  <- 0.012; label("Linear slope of age on CL, centred at 62 years (1/year)")         # Garcia-Hervalejo 2026 Equation 1: (1 + 0.012 x (AGE - 62)). Table 2 row 'AGE-CL' prints 0.01 at its uniform two-decimal precision (RSE 16%; bootstrap mean 0.01, 95% CI 0.01-0.02); the bootstrap CI of 0.01-0.02 is reproduced by 0.012 at RSE 16% but not by 0.010, which would print 0.01-0.01

    # ------------------------------------------------------------------
    # Between-subject variability. Methods: 'Interindividual variability in
    # pharmacokinetic parameters was assumed to follow a log-normal
    # distribution' and 'The magnitude of both IIV and RUV was expressed as a
    # coefficient of variation (CV%)'.
    #
    # Table 2 reports IIV as CV%, so the omega variances below are the exact
    # log-normal back-transform omega^2 = log(1 + CV^2):
    #   CL: log(1 + 0.341^2) = 0.110003
    #   Vd: log(1 + 0.310^2) = 0.091758
    # Reading the percentages instead as approximate SDs (omega = CV) would
    # give 0.116281 and 0.096100 -- omega 0.341 and 0.310 rather than 0.332 and
    # 0.303, a difference of under 3% on the SD that is immaterial for any
    # simulation use. Reading them as VARIANCES is not tenable: an omega^2 of
    # 34.1 is an omega of 5.8 on the log scale, which would place the 95%
    # interval of individual clearance across five orders of magnitude.
    # No off-diagonal covariance between the two etas was published.
    # ------------------------------------------------------------------
    etalcl ~ 0.110003  # Garcia-Hervalejo 2026 Table 2 row 'IIVCL (CV, %)' = 34.10 (RSE 23%; bootstrap mean 31.60, 95% CI 23.2-39.80) -> log(1 + 0.341^2). Reduced from 47.8% in the base model by adding IBW, eGFR and AGE
    etalvc ~ 0.091758  # Garcia-Hervalejo 2026 Table 2 row 'IIVV (CV, %)' = 31.0 (RSE 26%; bootstrap mean 30.10, 95% CI 21.0-39.7) -> log(1 + 0.310^2). Reduced from 35.8% in the base model, although no covariate was retained on Vd

    # ------------------------------------------------------------------
    # Residual variability: purely ADDITIVE, in mg/L. Results section 3.2:
    # 'Proportional and combined residual error models were evaluated but did
    # not result in an improvement in OFV or GOF diagnostics compared with the
    # additive error model. Given the narrow range of observed concentrations
    # and the predominance of trough samples, proportional error components
    # were poorly identifiable and did not meaningfully contribute to model
    # performance.'
    #
    # At the cohort mean concentration of 14.37 mg/L this corresponds to an
    # effective 18% coefficient of variation, but it is flat in absolute terms,
    # so simulated concentrations near the LLOQ of 3.0 mg/L carry proportionally
    # much more noise and can go negative. That is a property of the published
    # model, not of this encoding.
    # ------------------------------------------------------------------
    addSd <- 2.61; label("Additive residual error (mg/L)")  # Garcia-Hervalejo 2026 Table 2 row 'RUVadi (mg/L)' = 2.61 (RSE 21%; bootstrap mean 2.59, 95% CI 1.93-3.03)
  })

  model({
    # Individual PK parameters. All three covariates act on clearance and are
    # applied exactly as printed in Garcia-Hervalejo 2026 Equation 1: two power
    # terms normalised to the cohort medians, and one centred linear multiplier
    # on age. Vd carries log-normal between-subject variability but no
    # covariate.
    #
    # IBW is Devine ideal body weight in kg (see covariateData[[IBW]]$notes for
    # the formula, and for the source's typographical error in printing it);
    # CRCL is CKD-EPI eGFR in mL/min/1.73 m^2; AGE is in years.
    cl <- exp(lcl + etalcl) *
      (1 + e_age_cl * (AGE - 62)) *
      (CRCL / 92.15)^e_crcl_cl *
      (IBW / 61)^e_ibw_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One compartment with first-order elimination. Teicoplanin was given as a
    # 30-minute intravenous infusion directly into the systemic circulation
    # (Methods), so there is no absorption compartment and no bioavailability
    # term; supply the infusion via the event table's rate or duration column.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
