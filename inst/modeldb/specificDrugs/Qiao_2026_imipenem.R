Qiao_2026_imipenem <- function() {
  description <- "Two-compartment IV population PK model for imipenem in 46 Chinese critically ill adults on continuous renal replacement therapy (CVVH or CVVHDF) (Qiao 2026). Total clearance is the sum of an estimated endogenous clearance of 6.08 L/h, scaled by Cockcroft-Gault creatinine clearance through a power term centred on the 40.07 mL/min cohort median, and the individually calculated CRRT clearance supplied as the data column QEFF. Inter-individual variability on the endogenous clearance and the central volume; proportional residual error."
  reference <- "Qiao Z, Li D, Ai S, Zhang Z, Gan Y, Hu H, Luo C, Mo R, Huang J, He Y, Tao Y, Zheng X. Population Pharmacokinetic and Pharmacodynamic Analysis of Imipenem in Critically Ill Patients Receiving Continuous Renal Replacement Therapy. Drug Des Devel Ther. 2026;20:606942. doi:10.2147/DDDT.S606942"
  vignette <- "Qiao_2026_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation on total body weight; raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Qiao 2026 Table 1 footnote prints the estimating equation:",
        "CrCl = [(140 - age) x weight (kg)] / [0.818 x SCr], multiplied by 0.85 for women,",
        "with serum creatinine in umol/L (the footnote writes 'mmol/L', but the tabulated",
        "SCr median of 126 and range 34-673 are umol/L values; the 0.818 divisor is the",
        "umol/L form of the familiar 72). Cohort median 40.07 mL/min, range 6.70-157.98",
        "(Table 1); 40.07 is also the centring constant printed in Equation 10.",
        "Measured during CRRT, when extracorporeal removal lowers serum creatinine, so the",
        "value is an index of residual native renal function rather than a true GFR.",
        "The MDRD-type Chinese eGFR was screened on clearance as an alternative and not",
        "retained (dOFV -2.841, Discussion)."
      ),
      source_name = "CrCl"
    ),
    QEFF = list(
      description = "Individually calculated CRRT (CVVH or CVVHDF) extracorporeal clearance of imipenem",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Qiao 2026 Methods 'Quantification of CLCRRT', Equations 1-5: for CVVH",
        "CL_CRRT = SC x Q_E x [Q_b / (Q_b + Q_R)] / 1000 and for CVVHDF",
        "CL_CRRT = SA x Q_E x [Q_b / (Q_b + Q_R)] / 1000, where SC (sieving coefficient)",
        "and SA (saturation coefficient) are the measured effluent- or ultrafiltrate-to-plasma",
        "imipenem concentration ratios, Q_E the effluent flow (mL/h, = Q_d + Q_uf),",
        "Q_b the blood flow and Q_R the replacement-fluid flow; the bracket is the",
        "predilution correction. Supplied to the model as a per-subject data column",
        "(CL_CRRT,i of Equation 10) and ADDED to the endogenous clearance outside the",
        "random effect. Table 1: CL_CRRT median 2.48 L/h (range 1.37-3.72); Q_b median",
        "150 mL/min (115-200); Q_uf median 2362 mL/h (2029-3840). The paper reports no",
        "per-patient or median SC / SA values, so to simulate a prescribed CRRT dose",
        "(mL/h/kg) a user must supply a sieving coefficient and the flow settings to",
        "Equation 4 or 5 themselves and divide by 1000 to get L/h."
      ),
      source_name = "CL_CRRT"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Qiao 2026 Supplementary Table S1: entered on V1 (dOFV -16.17) and on V2 (dOFV -32.57) in forward inclusion, removed in backward elimination (p < 0.001 threshold). Table 1 median 60 kg (range 35-85). The 60 kg median is load-bearing OUTSIDE the model: the paper's Monte Carlo simulations convert the weight-normalised CRRT dose (20 / 37 / 74 mL/h/kg) into an absolute effluent flow at 60 kg."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Qiao 2026 Supplementary Table S1: entered on V1 in forward inclusion (dOFV -7.86), removed in backward elimination. Table 1 median 68 years (range 18-96)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Qiao 2026 Supplementary Table S1: entered on V1 in forward inclusion (dOFV -7.35), removed in backward elimination. Table 1 median 49.5 U/L (range 2.00-4281)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Qiao 2026 Supplementary Table S1 ('UREA'): entered on Q (dOFV -22.58) and on V2 (dOFV -13.30) in forward inclusion, removed in backward elimination. Table 1 median 9.55 mmol/L (range 3.50-50.30)."
    ),
    RBC = list(
      description = "Red blood cell count",
      units = "10^12 cells/L",
      type = "continuous",
      notes = "Qiao 2026 Supplementary Table S1: entered on V2 in forward inclusion (dOFV -24.84), removed in backward elimination. Table 1 median 2.84 x 10^12/L (range 1.57-4.95). Neutrophil percentage on Q (dOFV -6.54) and septic shock on Q (dOFV -8.78) and on V2 (dOFV -11.54) were likewise entered and then removed."
    )
  )

  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    age_range = "18-96 years",
    age_median = "68 years",
    weight_range = "35-85 kg",
    weight_median = "60 kg",
    sex_female_pct = 47.83,
    race_ethnicity = "Chinese (single-centre cohort, Chongqing)",
    disease_state = "Critically ill adults in the intensive care unit receiving imipenem-cilastatin together with continuous renal replacement therapy for at least 24 h. Sepsis or septic shock 40 (86.96%); acute kidney injury 17 (36.96%).",
    renal_function = "Cockcroft-Gault creatinine clearance median 40.07 mL/min (range 6.70-157.98), measured during CRRT. CRRT modality CVVH 38 (82.61%), CVVHDF 8 (17.39%); Baxter Prismaflex or NIKKISO Aquarius, citrate or heparin anticoagulation, filter run time 24-48 h. Blood flow median 150 mL/min (115-200); ultrafiltration flow median 2362 mL/h (2029-3840); calculated CRRT clearance median 2.48 L/h (1.37-3.72).",
    dose_range = "Imipenem-cilastatin 0.5 g q6h (22 patients, 47.82%) or 1 g q8h (24 patients, 52.17%) as a 2-h intravenous infusion (100 or 200 mL at 50 or 100 mL/h).",
    regions = "China (Chongqing University Central Hospital). Enrolled January-December 2025.",
    n_concentrations = 560L,
    notes = "Retrospective single-centre study; Table 1 gives the demographics. 560 plasma and 280 effluent samples; prefilter, postfilter and effluent drawn pre-dose and 0.5, 1, 1.5, 2, 3, 4, 5, 6 and 8 h after a dose at or near steady state. UPLC-MS/MS, LLOQ 0.02 mg/L. NONMEM 7.5.0 FOCE-I with PsN and Pirana; BQL handled by M1. Evaluated by GOF, NPDE, VPC (1000 simulations) and a 1000-sample bootstrap (998 converged)."
  )

  ini({
    # Structural parameters: Qiao 2026 Table 2, 'Final Model Estimate (RSE %)'
    # column. Equation 10 is the final clearance model:
    #
    #   CL_TOTAL = 6.08 * (CrCl_i / 40.07)^0.383 * exp(eta_CL,i) + CL_CRRT,i
    #
    # so 6.08 L/h is the endogenous clearance at the 40.07 mL/min cohort-median
    # creatinine clearance, and the per-patient CRRT clearance (data column
    # QEFF) is added after the random effect.
    lcl <- log(6.08); label("Endogenous clearance at CRCL = 40.07 mL/min (L/h)") # Table 2 CL 6.080 L/h (RSE 7%)
    lvc <- log(24.4); label("Central volume of distribution V1 (L)") # Table 2 V1 24.400 L (RSE 9%)
    lq <- log(3.19); label("Intercompartmental clearance Q (L/h)") # Table 2 Q 3.190 L/h (RSE 30%)
    lvp <- log(12.0); label("Peripheral volume of distribution V2 (L)") # Table 2 V2 12.000 L (RSE 36%)

    # Power exponent of CrCl on the endogenous clearance. Table 2 prints the
    # row as 'CLCR-CL, L/h', but the value is the dimensionless exponent of
    # Equation 10 (Results: 'the power coefficient for CrCl was 0.383').
    e_crcl_cl <- 0.383; label("Power exponent of CRCL on endogenous clearance (unitless)") # Table 2 CLCR-CL 0.383 (RSE 26%); Eq. 10

    # IIV: exponential model (Methods Eq. 6) with eta ~ N(0, omega^2). Table 2's
    # 'IIV-CL' and 'IIV-V1' rows are the NONMEM omega^2 VARIANCES. Reading
    # them as standard deviations instead flattens the paper's own Figure 6
    # target-attainment curves: across 24 digitised cells the variance
    # reading lands within an RMSE of about 2 percentage points and the SD
    # reading about 13 (e.g. 0.5 g q12h, MIC 4 mg/L, CrCl <= 30 mL/min:
    # published 75%, variance 82%, SD 97%). The vignette repeats the
    # comparison.
    # Table 2 IIV-CL 0.192 (RSE 29%, shrinkage 5%)
    etalcl ~ 0.192
    # Table 2 IIV-V1 0.161 (RSE 36%, shrinkage 24%)
    etalvc ~ 0.161

    # Proportional residual error (Methods Eq. 7, Y = IPRED * (1 + eps), eps ~
    # N(0, sigma^2)). Table 2 'PROB-erro' 0.0849 is sigma^2 on the same NONMEM
    # variance scale as the IIV rows, so the SD is sqrt(0.0849) = 0.291.
    propSd <- sqrt(0.0849); label("Proportional residual error (fraction)") # Table 2 PROB-erro 0.0849 (RSE 14%, shrinkage 10%)
  })
  model({
    # Endogenous clearance (Equation 10): power function of creatinine
    # clearance normalised to the 40.07 mL/min cohort median; the eta
    # multiplies this arm only.
    cl_body <- exp(lcl + etalcl) * (CRCL / 40.07)^e_crcl_cl

    # Total clearance: the per-patient CRRT clearance QEFF (L/h) is added
    # outside the random effect, as in Equation 10.
    cl <- cl_body + QEFF

    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give mg/L, the unit of the UPLC-MS/MS plasma
    # assay. Total (not unbound) plasma concentration: the paper applies no
    # protein-binding correction.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
