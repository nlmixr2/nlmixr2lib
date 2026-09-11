Zhou_2025_teicoplanin <- function() {
  description <- paste(
    "Two-compartment IV-infusion population PK model for teicoplanin in 79",
    "adult renal transplant recipients at a single Chinese centre (Zhou 2025).",
    "Real-time Cockcroft-Gault creatinine clearance is the only retained",
    "covariate and enters clearance as a POWER term normalised to 17 mL/min:",
    "CL = 0.711 * (CrCL / 17)^0.198 L/h. Between-subject variability is",
    "exponential on CL and Vc only; the inter-compartmental clearance CLd and",
    "the peripheral volume Vp carry no random effect. Residual error is",
    "proportional. The cohort is severely renally impaired early after",
    "transplantation (median baseline CrCL 9 mL/min), so the reported CL is",
    "lower and the total volume of distribution (Vc + Vp = 46.6 L) smaller",
    "than in the critically ill teicoplanin cohorts of Wang 2023 and Wi 2017.",
    "The paper's Monte Carlo simulations use the model to recommend CrCL-",
    "stratified loading and maintenance regimens against a trough target of",
    "15 mg/L and an AUC24/MIC target of 610.4, with a toxicity threshold of",
    "40 mg/L."
  )
  reference <- paste(
    "Zhou Y, Peng J, Xu P, Wang F, Xi J, Zhang H, Hu S, Yan H, Tan L, Cai H,",
    "Zhang B, Lan G. Population pharmacokinetics and dosing optimization of",
    "teicoplanin in renal transplant patients.",
    "Antimicrob Agents Chemother. 2025;69(6):e01568-24.",
    "doi:10.1128/aac.01568-24"
  )
  vignette <- "Zhou_2025_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Real-time creatinine clearance calculated with the Cockcroft-Gault formula, the only covariate retained in the final model (power effect on clearance)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CrCL. RAW Cockcroft-Gault creatinine clearance in mL/min,",
        "NOT normalised to 1.73 m^2 body surface area (Zhou 2025 Methods, 'Study",
        "design and patients': 'Creatinine clearance (CrCL) was calculated by the",
        "Cockcroft-Gault formula'). Supply the value on that raw scale; a",
        "BSA-normalised eGFR silently rescales the renal term. Precedent for raw",
        "mL/min under the canonical CRCL column: Delattre_2010_amikacin.R and",
        "Chen_2023_nemonoxacin.R.",
        "TIME-VARYING. The paper screened BOTH a time-fixed baseline CrCL (source",
        "column BCRCL) and a REAL-TIME CrCL updated through the treatment course,",
        "and retained the real-time version: Table S2 gives CL-BCRCL dOFV -5.5",
        "(below the -3.84 forward-inclusion threshold) versus CL-CRCL dOFV -10.2,",
        "restated in the Results text as dOFV 10.203. That distinction is central",
        "to the paper: renal function recovers rapidly in the first 2-4 weeks after",
        "transplantation (Introduction), so a subject's CrCL rises substantially",
        "within a single teicoplanin course and a baseline-only covariate cannot",
        "track it. Only the time-varying column is used here, so the model needs a",
        "single CRCL column and no companion CRCL_BASE column.",
        "Enters the final model as the POWER term printed in the Zhou 2025 Results",
        "equation block: CL (L/h) = 0.711 * (CrCL / 17)^0.198 * e^eta1. The",
        "normalising constant 17 mL/min is printed in that equation; it is NOT the",
        "cohort median BASELINE CrCL, which is 9 mL/min in the development set and",
        "6.8 mL/min overall (Table 1), and it sits above both because the real-time",
        "values rise as the graft recovers.",
        "The exponent 0.198 is weak, so the multiplier spans only 0.55 at",
        "CrCL = 1 mL/min to 1.47 at CrCL = 120 mL/min. The Monte Carlo simulations",
        "cover CrCL <=10, 10-30, 30-60, 60-90 and 90-120 mL/min (Methods, 'Monte",
        "Carlo simulations and dosing optimization'), which is the range over which",
        "the model was exercised; the observed cohort itself is concentrated at the",
        "low end. None of the enrolled patients received renal replacement therapy",
        "or mechanical ventilation (Results, 'Patient characteristics'), so CrCL is",
        "the sole renal descriptor and no dialysis-clearance term is present.",
        "The Results prose describes the retained relationship as 'a proportional",
        "model', which conflicts with the power form printed in the equation block",
        "on the same page. The printed equation governs (standing policy: text vs",
        "printed-equation conflict, trust the equation), and it is corroborated",
        "three ways. (i) Table 2 labels the estimate 'CrCLCL: the influence",
        "coefficient of CrCL on CL', which is exponent-like language. (ii) The",
        "linear reading is arithmetically impossible in this cohort: a centred",
        "proportional form 0.711 * (1 + 0.198 * (CrCL - 17)) turns NEGATIVE below",
        "CrCL = 17 - 1/0.198 = 11.9 mL/min, and the development set's median",
        "baseline CrCL is 9 mL/min, so more than half the modelled subjects would",
        "have negative clearance (the ODE system does not solve at all under that",
        "reading). (iii) The power reading reproduces the paper's own printed",
        "Monte Carlo table (Table 3, 160 cells) to a median absolute difference",
        "under 2 percentage points; see vignettes/articles/Zhou_2025_teicoplanin.Rmd."
      ),
      source_name        = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "23 of 79 development-set subjects were female (56/79 male, Zhou 2025 Table 1). Screened on CL (dOFV -0.6) and on Vc (dOFV -0.7) and not retained (Table S2)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 42.5 years (IQR 36-52.8) in the development set (Zhou 2025 Table 1). Screened on CL (dOFV -3.2) and on Vc (dOFV -2.5); both fall short of the -3.84 forward-inclusion threshold (Table S2)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Median 64.8 kg (IQR 54-74.8) in the development set (Zhou 2025 Table 1).",
        "Screened on Vc (dOFV -3.7) and on CL (dOFV -2.5) and not retained",
        "(Table S2). The Discussion attributes the null result to the limited",
        "sample size and the narrow weight distribution, explicitly contrasting it",
        "with earlier teicoplanin models that did retain weight, and states that",
        "the Monte Carlo simulations were therefore run with FIXED (not",
        "weight-banded) doses."
      )
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Median 37.0 g/L (IQR 34.4-39.5) in the development set (Zhou 2025",
        "Table 1). Screened on CL (dOFV -0.6) and on Vc (dOFV -1.3) and not",
        "retained (Table S2, source column BALB). The Discussion notes that",
        "albumin acts on teicoplanin volume in critically ill and haematological",
        "cohorts but that this cohort's albumin sat within the normal range."
      )
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Median 6.3 umol/L (IQR 4.9-9) in the development set (Zhou 2025 Table 1). Screened on CL (dOFV -2.3) and on Vc (dOFV -1.3) and not retained (Table S2, source column BTBIL). Direct bilirubin was collected but dropped before covariate screening for collinearity with total bilirubin (r = 0.8 > 0.6, Table S2 footnote a)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Median 12.9 U/L (IQR 8.7-16.6) in the development set (Zhou 2025 Table 1). Screened on CL (dOFV -3.3) and on Vc (dOFV -0.5) and not retained (Table S2, source column BALT). Aspartate aminotransferase was collected but dropped before covariate screening for collinearity with ALT (r = 0.8 > 0.6, Table S2 footnote b)."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault), time-fixed",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Median 9 mL/min (IQR 6.6-13.8) in the development set (Zhou 2025",
        "Table 1). Screened on CL (dOFV -5.5) and on Vc (dOFV -0.5) under source",
        "column BCRCL. Although the CL effect cleared the -3.84 forward-inclusion",
        "threshold, the REAL-TIME CrCL gave a substantially larger drop on the same",
        "parameter (dOFV -10.2) and was carried forward instead (Table S2), so the",
        "final model has no baseline-CrCL term. Documented here only as a screened",
        "covariate."
      )
    ),
    T_POSTTRANSPLANT_GT1MO = list(
      description = "Post-transplantation interval indicator: 1 = more than 1 month since transplantation, 0 = within 1 month",
      units       = "(binary)",
      type        = "binary",
      notes       = "Source column DFLAG, described in the Table S2 abbreviations as 'a categorical covariate distinguishing postoperative duration intervals, dichotomized into two temporal strata: the short-term postoperative phase (<=1 month) and extended postoperative period (>1 month)'. 65 of 79 development-set subjects were within 1 month of transplantation (Zhou 2025 Table 1). Screened on CL (dOFV 0) and on Vc (dOFV 0) and not retained (Table S2). Listed here as documentation only; the name is not registered in inst/references/covariate-columns.md because the model does not use it."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 79L,
    n_studies        = 1L,
    n_concentrations = 306L,
    age_median       = "42.5 years (development set); 42.0 years overall",
    age_range        = "IQR 36-52.8 years in the development set (Zhou 2025 Table 1); the study enrolled adults aged >= 18 years",
    weight_median    = "64.8 kg (development set); 62.7 kg overall",
    weight_range     = "IQR 54-74.8 kg in the development set (Zhou 2025 Table 1)",
    height_median    = "168 cm (development set)",
    sex_female_pct   = 29.1,
    race_ethnicity   = "Not reported by category; single-centre Chinese cohort (Second Xiangya Hospital, Central South University, Changsha, Hunan)",
    disease_state    = paste(
      "Adult renal transplant recipients (age >= 18 years) receiving",
      "intravenous teicoplanin for at least 72 h. 80 of the 99 enrolled",
      "patients received it as perioperative antimicrobial prophylaxis and 19",
      "for confirmed or probable Gram-positive infection (Results, 'Patient",
      "characteristics'; Table S1 details the 19: pneumonia 14, urinary tract",
      "infection 5; pathogens MRSA 9, Enterococcus faecium 5; teicoplanin MICs",
      "0.5-1 mg/L for MRSA and 1 mg/L for E. faecium). Most were early",
      "post-transplant: 65 of 79 development-set subjects were within 1 month of",
      "surgery, 4 between 1 and 12 months and 10 beyond 12 months. Comorbidities",
      "were hypertension 81.0% and diabetes 25.3%. All 99 patients were on",
      "tacrolimus and mycophenolate mofetil. NONE were on mechanical ventilation",
      "or renal replacement therapy. Median teicoplanin treatment duration was 9",
      "days (IQR 8-12) in the development set."
    ),
    dose_range       = paste(
      "Teicoplanin (Zhejiang Medicine Co., Ltd) by 1 h intravenous infusion. All",
      "patients received the same empirical regimen: three loading doses of 400",
      "mg every 12 h followed by a maintenance dose of 400 mg once daily for at",
      "least 3 days (Methods, 'Dosing regimen and blood sampling'). For",
      "perioperative prophylaxis the first dose was given 2 h before surgery.",
      "The published Monte Carlo simulations (Table 3) explore loading doses of",
      "400, 600, 800 and 1000 mg q12h given either three or five times, followed",
      "by maintenance doses of 200-1000 mg q24h."
    ),
    regions          = "China (Second Xiangya Hospital, Changsha, Hunan; development cohort enrolled January 2022 to December 2023)",
    renal_function   = paste(
      "Severely impaired at baseline and recovering during treatment.",
      "Development set (Zhou 2025 Table 1): baseline Cockcroft-Gault creatinine",
      "clearance median 9 mL/min (IQR 6.6-13.8); baseline blood urea nitrogen",
      "median 18.8 mmol/L (IQR 12.8-26.1). No patient received renal replacement",
      "therapy. The Introduction notes that serum creatinine falls rapidly and",
      "renal function typically normalises within 2-4 weeks of transplantation,",
      "which is why the final model uses REAL-TIME rather than baseline CrCL; the",
      "17 mL/min normalising constant in the CL equation reflects the higher",
      "real-time values rather than the baseline median."
    ),
    screened_covariates = paste(
      "Tested-but-not-retained (Zhou 2025 Table S2, dOFV against the",
      "two-compartment zero-order-input base model at OFV 1246.3; forward",
      "inclusion threshold dOFV <= -3.84): on CL - sex -0.6, post-transplant",
      "interval (DFLAG) 0, age -3.2, body weight -2.5, baseline total bilirubin",
      "-2.3, baseline albumin -0.6, baseline ALT -3.3, baseline CrCL -5.5;",
      "on Vc - sex -0.7, DFLAG 0, age -2.5, body weight -3.7, baseline total",
      "bilirubin -1.3, baseline albumin -1.3, baseline ALT -0.5, baseline CrCL",
      "-0.5, real-time CrCL -0.1. Only real-time CrCL on CL (dOFV -10.2,",
      "final-model OFV 1236.1) was retained. Direct bilirubin and AST were",
      "dropped before screening for collinearity (r = 0.8) with total bilirubin",
      "and ALT respectively."
    ),
    external_validation = paste(
      "A prospectively recruited independent cohort of 20 patients contributing",
      "80 samples, enrolled November 2024 to March 2025 (Methods, 'Study design",
      "and patients'). Against the final model the root mean square error was",
      "8.49% and the mean prediction error 0.704% (Results, 'Model evaluation and",
      "validation'); goodness-of-fit plots and a prediction-corrected VPC on the",
      "validation set are Figures S2 and 3."
    ),
    notes            = paste(
      "Single-centre prospective two-phase study; 99 patients and 386 plasma",
      "samples in total, split 79 patients / 306 samples for model development",
      "and 20 patients / 80 samples for external validation. Sampling: two to",
      "four 2 mL samples per patient, drawn immediately before the seventh dose",
      "and at 0, 1, 3, 5, 7, 9, 11, 13, 14, 17, 18, 21 and 22 h after the end of",
      "an infusion, so each development-set patient contributed about four",
      "samples. Assay: automated two-dimensional LC-MS/MS quantifying the",
      "teicoplanin A2-2 isoform as the surrogate for teicoplanin, with",
      "teicoplanin itself as internal standard; linear 5.0-100 mg/L,",
      "LLOD 0.5 mg/L, intra- and inter-day CV <= 3.72%. Estimation: Phoenix NLME",
      "8.1 with FOCE-ELS. Qualification: goodness-of-fit plots, a 1000-run",
      "bootstrap (992 successful), a prediction-corrected VPC, and the external",
      "validation above. The authors note as a limitation that the short sampling",
      "window prevented the model from capturing teicoplanin's very long terminal",
      "phase (published terminal half-life 83-163 h), so this two-compartment",
      "model describes the alpha and beta phases only and should not be used to",
      "project exposure far beyond the simulated 192 h horizon."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects. Zhou 2025 Table 2, 'Final model / Estimate'
    # column, and the typical-value equation block printed in Results
    # ('Population pharmacokinetic modeling'):
    #   CL (L/h)  = 0.711 * (CrCL / 17)^0.198 * e^eta1
    #   Vc (L)    = 11.3 * e^eta2
    #   CLd (L/h) = 4.22,  Vp (L) = 35.3
    #
    # lcl is the clearance of the REFERENCE subject at CrCL = 17 mL/min, the
    # normalising constant printed in that equation. At the cohort median
    # BASELINE CrCL of 9 mL/min the model gives
    # 0.711 * (9 / 17)^0.198 = 0.624 L/h.
    # ------------------------------------------------------------------
    lcl <- log(0.711); label("Clearance at CrCL = 17 mL/min (L/h)")             # Zhou 2025 Table 2: CL = 0.711 L/h (SE 6%; bootstrap median 0.714, 90% CI 0.614-0.802)
    lvc <- log(11.3);  label("Central volume Vc (L)")                           # Zhou 2025 Table 2: Vc = 11.3 L (SE 7%; bootstrap median 11.3, 90% CI 9.93-12.8)
    lq  <- log(4.22);  label("Inter-compartmental clearance CLd (L/h)")         # Zhou 2025 Table 2: CLd = 4.22 L/h (SE 11%; bootstrap median 4.22, 90% CI 3.36-5.24)
    lvp <- log(35.3);  label("Peripheral volume Vp (L)")                        # Zhou 2025 Table 2 and the Results equation block: Vp = 35.3 L (SE 10%; bootstrap median 35, 90% CI 28.4-42.8). The Abstract and the sentence under Table 2 both round this to 35.2; 35.3 is used because it is the value in the parameter table AND in the printed equation, and because Vc + Vp = 11.3 + 35.3 = 46.6 L is the total volume the Discussion quotes.

    # ------------------------------------------------------------------
    # Covariate effect: POWER exponent of real-time Cockcroft-Gault CrCL on
    # CL, normalised to 17 mL/min. Applied in model() as (CRCL / 17)^e_crcl_cl,
    # matching the printed equation exactly.
    #
    # The Results prose calls this 'a proportional model'; the equation printed
    # three lines later is a power model. The equation governs -- see the
    # covariateData[[CRCL]]$notes for why the magnitude of 0.198 also rules out
    # a per-mL/min linear slope.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.198; label("Power exponent of real-time CrCL on CL (unitless)")  # Zhou 2025 Table 2 row 'CrCLCL', the influence coefficient of CrCL on CL = 0.198 (SE 30%; bootstrap median 0.196, 90% CI 0.073-0.319)

    # ------------------------------------------------------------------
    # Between-subject variability, exponential (Methods: 'The inter-individual
    # variability of the PK parameters was evaluated using an exponential
    # model'). Zhou 2025 Table 2 reports the random effects under a
    # 'Random effects (%)' header, as omega (the SD on the log scale) x 100.
    #
    # The column is SD-like, not variance-like: the residual row in the same
    # block reads 14.6% and the external validation returned an RMSE of 8.49%
    # (Results). Read as a variance, 14.6% would imply a residual SD of 38%,
    # which is irreconcilable with that RMSE and with a 3.72%-CV LC-MS/MS
    # assay. The variances below are therefore the squared percentages. Reading
    # the same percentages as approximate CVs instead -- omega =
    # sqrt(log(1 + CV^2)) -- would give 0.386 and 0.361 rather than 0.401 and
    # 0.373, a difference too small to matter for any simulation use.
    # ------------------------------------------------------------------
    etalcl ~ 0.160801  # Zhou 2025 Table 2: omega_CL = 40.1% (SE 12%; bootstrap median 39.6%, 90% CI 30.6-49.2) -> omega^2 = 0.401^2
    etalvc ~ 0.139129  # Zhou 2025 Table 2: omega_Vc = 37.3% (SE 21%; bootstrap median 36.2%, 90% CI 16.7-53.0) -> omega^2 = 0.373^2

    # Zhou 2025 Table 2 reports no random effect on CLd or on Vp, and the
    # printed equation block writes both as bare constants
    # ('CLd (L/h) = 4.22, Vp (L) = 35.3') with no e^eta factor, in contrast to
    # the CL and Vc lines which both carry one. No eta is declared on either.
    # No off-diagonal covariances were published.

    # ------------------------------------------------------------------
    # Residual variability: proportional only (Results: 'a two-compartment
    # model with first-order elimination with proportional residual variability
    # and interindividual variability').
    #
    # Phoenix parameterises its proportional error on the log scale as
    # ln Y_ij = ln F_ij + eps_1,ij (the second of the three candidate forms
    # printed in Methods), which is formally exponential rather than
    # proportional. At an SD of 0.146 the two forms differ by well under 1% in
    # the simulated concentration distribution, and the canonical nlmixr2lib
    # encoding of a residual the source labels 'proportional' is prop(propSd).
    # ------------------------------------------------------------------
    propSd <- 0.146; label("Proportional residual error (fraction)")  # Zhou 2025 Table 2: proportional residual = 14.6% (SE 14%; bootstrap median 14.3%, 90% CI 10.4-18.2), encoded as a fraction
  })

  model({
    # Individual PK parameters. Real-time Cockcroft-Gault CrCL (canonical
    # column CRCL, raw mL/min) enters CL only, as the power term normalised to
    # 17 mL/min printed in the Zhou 2025 Results equation block. Vc carries
    # log-normal IIV with no covariate; CLd and Vp carry neither.
    cl <- exp(lcl + etalcl) * (CRCL / 17)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Intravenous infusion into the central compartment (1 h infusions in the
    # study; Zhou 2025 Methods, 'Dosing regimen and blood sampling', and the
    # Table S2 base-model description 'two-compartment model with a zero-order
    # input rate'), first-order distribution to peripheral1 and first-order
    # elimination from central. Dose in mg and volumes in L give central / vc
    # in mg/L.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
