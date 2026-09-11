Wang_2024_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 120 hospitalized",
    "Chinese patients aged 60 years or older (Wang 2024), including",
    "critically ill patients. Imipenem-cilastatin was given empirically as",
    "an IV infusion of 250-1000 mg every 6-12 h and monitored by",
    "therapeutic drug monitoring. Clearance scales as a power of raw",
    "Cockcroft-Gault creatinine clearance (reference 71 mL/min); central",
    "volume, peripheral volume and intercompartmental clearance carry no",
    "covariate effects. Inter-individual variability was retained only on",
    "clearance; the omegas on Vc, Q and Vp were constrained to zero for",
    "lack of precision. Residual variability is additive. Fitted in NONMEM",
    "7.3.0."
  )
  reference <- paste(
    "Wang J, Fang Q, Luo X, Jin L, Zhu H.",
    "Population pharmacokinetics and dosing optimization of imipenem in",
    "Chinese elderly patients.",
    "Front Pharmacol. 2025;15:1524272.",
    "doi:10.3389/fphar.2024.1524272.",
    sep = " "
  )
  vignette <- "Wang_2024_imipenem"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated by the Cockcroft-Gault equation,",
        "reported raw in mL/min and NOT BSA-normalised to mL/min/1.73 m^2"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Paper Methods 'Data collection and sampling schedule': 'The",
        "creatinine clearance rate (CLCR) was calculated by Cockcroft-Gault",
        "Equation.' Reported raw in mL/min; the cohort median is 58.9",
        "mL/min (IQR 35.54-97.33, paper Table 1). Table 1 lists a separate",
        "BSA-normalised eGFR row (median 86.9 mL/min/1.73 m^2) which is NOT",
        "the covariate used in the model -- supply raw Cockcroft-Gault",
        "mL/min, not the eGFR value. The reference constant in the power",
        "model is 71 mL/min (paper Equation 1), which is NOT the cohort",
        "median 58.9 mL/min; the paper does not state how 71 was chosen, so",
        "it must be taken from Equation 1 verbatim rather than recomputed",
        "from Table 1. Used in the multiplicative power model",
        "CL = 13.1 * (CLCR/71)^0.263 * exp(eta_CL) (paper Equation 1).",
        "The Monte Carlo dosing simulations stratify virtual patients into",
        "CLCR bands 0-30, 30-60, 60-90 and 90-120 mL/min (paper Results",
        "'Monte Carlo simulations of dosage regimens'). Stored under the",
        "canonical CRCL column per inst/references/covariate-columns.md,",
        "which accepts raw measured / Cockcroft-Gault CrCL in mL/min when",
        "the source paper does not BSA-normalise, provided the per-model",
        "description records the assay form -- precedent: Delattre 2010",
        "amikacin, Couffignal 2014 imipenem, Lamoth 2009 imipenem."
      ),
      source_name        = "CLCR"
    )
  )

  # Covariates screened during model building but NOT retained in the final
  # model. Documented here (rather than in covariateData) because they are
  # never referenced in model(); see Supplementary Table S1, which reports
  # the full forward-inclusion / backward-elimination trace. Each entered CL
  # multiplicatively during forward inclusion and was dropped again in
  # recursive backward elimination because removal raised the OFV by less
  # than the 6.63 (p < 0.01) retention threshold. The paper prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = paste(
        "Cohort median 63.8 mg/L (IQR 26.92-115.14, paper Table 1). Added",
        "to CL in forward inclusion (model 3, dOFV -6.02, p < 0.05) and",
        "removed in backward elimination (model 8, dOFV +6.02 < 6.63;",
        "Supplementary Table S1). Not in the final model; no coefficient",
        "is reported."
      )
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = paste(
        "Cohort median 8.6 x10^9/L (IQR 5.67-13.07, paper Table 1). Added",
        "to CL in forward inclusion (model 4, dOFV -4.47, p < 0.05) and",
        "removed first in backward elimination (model 6, dOFV +6.42 < 6.63;",
        "Supplementary Table S1). Not in the final model; no coefficient",
        "is reported."
      )
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy during imipenem treatment",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "24 of 120 modeling-group patients (20%) received CRRT during",
        "imipenem therapy (paper Table 1). Added to CL in forward inclusion",
        "(model 5, dOFV -6.64, p < 0.01) and removed in backward elimination",
        "(model 7, dOFV +4.69 < 6.63; Supplementary Table S1). The paper's",
        "Discussion gives a third, inconsistent figure for the same",
        "screening step -- 39 patients on CRRT and an OFV decrease of 5.81",
        "-- against Supplementary Table S1's 6.64; both agree the covariate",
        "was dropped. The Discussion argues the CRRT effect is already",
        "carried implicitly by CLCR, since CLCR is derived from serum",
        "creatinine which reflects CRRT solute removal. Not in the final",
        "model; no coefficient is reported."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 120L,
    n_studies        = 1L,
    age_range        = ">= 60 years (inclusion criterion)",
    age_median       = "72 years (IQR 68-81)",
    weight_range     = "35-93.5 kg",
    weight_median    = "65 kg (IQR 59.00-65.33)",
    sex_female_pct   = 35.0,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Hospitalized Chinese patients aged 60 years or above treated",
      "empirically with imipenem-cilastatin sodium for injection and",
      "monitored by therapeutic drug monitoring, including patients in",
      "critical condition. Median C-reactive protein 63.8 mg/L, white",
      "blood cell count 8.6 x10^9/L, procalcitonin 0.593 ng/mL and serum",
      "albumin 33 g/L indicate an acutely infected, hypoalbuminaemic",
      "cohort (paper Table 1). 24 patients (20%) received continuous",
      "renal replacement therapy during imipenem therapy. Patients on",
      "extracorporeal membrane oxygenation (ECMO) were excluded."
    ),
    dose_range       = paste(
      "IV infusion of 250-1000 mg imipenem, dosing interval every 6 h to",
      "every 12 h, chosen empirically by the treating clinicians (paper",
      "Methods 'Quantification of imipenem concentrations'). The infusion",
      "rate was set from the actual infusion duration recorded in the",
      "Electronic Health Record; the paper does not report the duration",
      "distribution. The Monte Carlo dosing simulations evaluate 0.25 g",
      "q6h, 0.5 g q6h, 0.5 g q8h, 1 g q6h, 1 g q8h and 1 g q12h."
    ),
    regions          = "China (single centre: Nanjing Drum Tower Hospital, Nanjing, Jiangsu)",
    renal_function   = paste(
      "Cockcroft-Gault creatinine clearance median 58.9 mL/min (IQR",
      "35.54-97.33, raw mL/min, not BSA-normalised); BSA-normalised eGFR",
      "median 86.9 mL/min/1.73 m^2 (IQR 50.12-132.20); serum creatinine",
      "median 232 umol/L (IQR 150.25-354.00) (paper Table 1). Table 1",
      "also carries a separate 'CREA' row with median 75.00 umol/L, which",
      "is inconsistent with the 'SCR' row; the paper does not reconcile",
      "them and neither value enters the model directly."
    ),
    n_concentrations = 370L,
    notes            = paste(
      "Retrospective single-centre observational study, October 2021 to",
      "April 2024 (paper Methods 'Study design and ethics'; ethics",
      "approval No. 2023-380-02). A total of 142 patients contributing 370",
      "plasma concentration records met the inclusion criteria; 120 were",
      "used to develop the model and a separate 22 formed an external",
      "validation cohort (paper Results 'Patient demographics'). The",
      "abstract instead states that all 370 observations from 142 patients",
      "were incorporated in the PPK model; Table 1 and the Results text",
      "give the 120 / 22 split, and the per-group observation counts are",
      "not reported. Sampling was sparse and predominantly trough",
      "concentrations, which the paper lists as a limitation. Imipenem in",
      "plasma was measured by HPLC-UV at 300 nm after ultrafiltration,",
      "with 3-morpholine propyl sulfonic acid as stabiliser; calibration",
      "linear over 0.5-50 ug/mL with a limit of quantitation of 0.5",
      "ug/mL. Observed concentrations had median 1.8 ug/mL (IQR",
      "0.3-2.775). 1 ug/mL = 1 mg/L, the units declared above. Model",
      "fitted in NONMEM 7.3.0. Final-model 95% confidence intervals from",
      "1000 nonparametric bootstrap resamples, of which 967 (96.7%)",
      "converged. External validation gave MPE 7.6%, MAPE 39.4%, F20 34%",
      "and F30 52.8%."
    )
  )

  ini({
    # ===== Structural PK (Wang 2024 Table 2 'Final model' column, and
    # Equations 1-4) =====
    # Reference subject for CL: CLCR = 71 mL/min (paper Equation 1).
    lcl <- log(13.1); label("Typical clearance at CLCR=71 mL/min (L/h)")  # Wang 2024 Table 2 final model: CL = 13.1 L/h (RSE 4.80%; bootstrap median 13.3, 95% CI 9.11-14.4); Equation 1
    lvc <- log(11.7); label("Typical central volume Vc (L)")  # Wang 2024 Table 2 final model: Vc = 11.7 L (RSE 5.20%; bootstrap median 11.7, 95% CI 3.41-12.7); Equation 2
    lq  <- log(11.9); label("Typical intercompartmental clearance Q (L/h)")  # Wang 2024 Table 2 final model: Q = 11.9 L/h (RSE 24.5%; bootstrap median 12.1, 95% CI 5.67-18.0); Equation 3 prints the unit as "(L)", which is a typesetting error -- Table 2 and the abstract both give L/h, and L/h is the only unit dimensionally consistent with a two-compartment model
    lvp <- log(29.3); label("Typical peripheral volume Vp (L)")  # Wang 2024 Table 2 final model: Vp = 29.3 L (RSE 12.3%; bootstrap median 30.3, 95% CI 17.5-41.4); Equation 4

    # ===== Covariate effect =====
    # Multiplicative power model on CL, paper Equation 1:
    #   CL (L/h) = 13.1 * (CLCR/71)^0.263 * EXP(eta_CL)
    # CLCR was the only covariate retained after forward inclusion and
    # backward elimination (paper Results 'Population pharmacokinetic
    # modeling'; Supplementary Table S1). Estimated, not fixed.
    e_crcl_cl <- 0.263; label("Power exponent on (CRCL/71) for CL (unitless)")  # Wang 2024 Table 2 final model: theta5 (CLCR on CL) = 0.263 (RSE 14.0%; bootstrap median 0.264, 95% CI 0.196-0.341); Equation 1 exponent

    # ===== Inter-individual variability (Wang 2024 Table 2) =====
    # Exponential IIV: CL_i = CL_typ * exp(eta_CL,i), eta ~ N(0, omega^2)
    # (paper Methods 'Base model': "Exponential error models were employed
    # to define inter-individual variability (IIV, eta) ... with a mean of
    # zero and a variance of omega^2"). Table 2 labels the row "omega^2 CL",
    # so 0.0832 is the VARIANCE and is used directly -- nlmixr2 takes
    # variances on the eta line. omega^2 = 0.0832 corresponds to
    # sqrt(exp(0.0832) - 1) = 29.5% CV.
    #
    # No eta is declared on Vc, Q or Vp. Paper Results 'Population
    # pharmacokinetic modeling': "In the final model, the inter-individual
    # variation values for Vc, Q, and Vp were constrained due to inadequate
    # precision in the omega estimates." Table 2 reports no omega for them.
    # They are deliberately omitted rather than written as `~ fixed(0)`,
    # which would make OMEGA singular and break simulation.
    etalcl ~ 0.0832  # Wang 2024 Table 2 final model: omega^2 CL = 0.0832 (RSE 16.7%, shrinkage 10.6%; bootstrap median 0.0802, 95% CI 0.0545-0.111)

    # ===== Residual error (Wang 2024 Table 2) =====
    # Additive-only error model, selected over exponential, proportional
    # and combined (paper Results 'Population pharmacokinetic modeling':
    # "the additive residual error model assuming a normal distribution was
    # determined to best fit the data").
    #
    # Table 2 tabulates 0.575 under the heading "Residual variability
    # (RSV)", row "Additive error". That number is a VARIANCE, not an SD:
    # paper Methods 'Base model' defines every residual error model as
    # "assuming a symmetric distribution around a mean of zero with a
    # variance represented by sigma^2", and the additive form is written
    # Cobs = Cpred + epsilon -- so the estimated quantity for that row is
    # sigma^2, which is also what NONMEM's $SIGMA block reports by default.
    # The companion IIV row in the same table is likewise a variance
    # ("omega^2 CL"). nlmixr2's add() takes a STANDARD DEVIATION, hence
    # sqrt(0.575) = 0.7583 mg/L. See the vignette's "Assumptions and
    # deviations" section, which records this reading and the checks run
    # against Supplementary Figure S1.
    addSd <- sqrt(0.575); label("Additive residual error (mg/L)")  # Wang 2024 Table 2 final model: additive sigma^2 = 0.575 (RSE 13.8%, shrinkage 8.40%; bootstrap median 0.572, 95% CI 0.422-0.721); sqrt() converts the reported variance to the SD that add() expects
  })

  model({
    # ----- Individual PK parameters -----
    # Paper Equations 1-4. Only CL carries a covariate and an eta; Vc, Q
    # and Vp are typical values with no covariate and no IIV.
    cl <- exp(lcl + etalcl) * (CRCL / 71)^e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # Imipenem-cilastatin is given as an IV infusion straight into the
    # central compartment; the infusion duration comes from the event
    # table's rate / dur column (paper Methods 'Quantification of imipenem
    # concentrations': "The infusion rate was established based on the
    # actual infusion duration documented in the Electronic Health
    # Record"). No absorption compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    # Dose in mg, vc in L -> mg/L, which equals the ug/mL reporting units
    # of the HPLC-UV assay (paper Methods) and matches the "Serum imipenem
    # concentration (mg/L)" axis of Supplementary Figures S1 and S2.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
