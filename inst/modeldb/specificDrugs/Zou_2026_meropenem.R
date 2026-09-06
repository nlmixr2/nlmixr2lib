Zou_2026_meropenem <- function() {
  description <- "One-compartment IV population PK model for meropenem in 44 Chinese adults with severe postoperative infections receiving therapeutic drug monitoring in a general ICU (Zou 2026). Postoperative Cockcroft-Gault creatinine clearance acts as a power covariate on CL, centred at the cohort mean of 47.7 mL/min. Inter-individual variability on CL and Vc is strongly correlated (rho = 0.94). Residual error is additive; see the vignette Errata for the reconciliation of Table 2's residual-error row, which reports a variance where the bootstrap column reports a standard deviation."
  reference <- paste(
    "Zou Y, Ren J, Lei H, Chen G, Li C, He X, Hu Y, Liu X.",
    "Population pharmacokinetic modeling and Monte Carlo simulation to",
    "optimize meropenem dosing in patients with severe postoperative",
    "infections.",
    "Front Pharmacol. 2026;17:1778552.",
    "doi:10.3389/fphar.2026.1778552.",
    sep = " "
  )
  vignette <- "Zou_2026_meropenem"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Postoperative Cockcroft-Gault creatinine clearance (raw, not BSA-normalised)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Postoperative creatinine clearance, calculated with the",
        "Cockcroft-Gault equation (Zou 2026 Table 1 footnote a; Table 3",
        "footnote 3), was the only covariate retained by stepwise covariate",
        "modelling (forward alpha = 0.05, backward alpha = 0.001).",
        "It enters CL as a power function centred at the cohort mean",
        "47.7 mL/min (Zou 2026 Section 3.2 equation):",
        "CL = 6.472 * (CrCL/47.7)^0.3834 * exp(eta_CL).",
        "Cohort mean 47.7 mL/min, range 11.9-136 mL/min (Table 1).",
        "Note that Section 3.1 prose instead reports a mean CRCL of",
        "54.3 +/- 29.6 mL/min over the same 11.9-136 range; the 47.7 value",
        "is the one used as the model centring constant and is the value",
        "printed inside the covariate equation, so it is the value encoded",
        "here. Monte Carlo simulations were stratified into CRCL bands of",
        "<10, 10-25, 26-50, 51-90 and 91-140 mL/min (Section 2.6).",
        "Stored under canonical CRCL with raw mL/min (the source-aliases",
        "entry in inst/references/covariate-columns.md permits the raw",
        "Cockcroft-Gault form when the source paper does not apply",
        "BSA-normalisation), following AbdulAziz_2016_doripenem.R."
      ),
      source_name        = "CRCL"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 44L,
    n_studies        = 1L,
    n_concentrations = 135L,
    age_range        = "20.0-91.0 years (Section 3.1 prose); Table 1 reports mean 66.0, range 20.0-83.0 years",
    age_mean         = "68.5 +/- 14.8 years (Section 3.1 prose); 66.0 years (Table 1)",
    weight_range     = "40.0-81.0 kg (Section 3.1 prose); Table 1 reports mean 58.0, range 45.0-81.0 kg",
    weight_mean      = "56.7 +/- 9.04 kg (Section 3.1 prose); 58.0 kg (Table 1)",
    sex_female_pct   = 29.55,
    race_ethnicity   = "Chinese (single-centre ICU, The Central Hospital of Xiangtan, Hunan, China)",
    disease_state    = paste(
      "Adults (>18 years) with severe postoperative infections requiring",
      "meropenem for more than 3 days, enrolled in a general ICU. Intra-",
      "abdominal infection combined with hospital-acquired pneumonia was",
      "the most common presentation (47.73%), followed by hospital-acquired",
      "pneumonia alone (22.73%), intra-abdominal infection alone (20.45%)",
      "and other (9.09%). 93.18% were classified as critically ill and",
      "15.91% were receiving haemodialysis. Patients with preoperative",
      "meropenem-treated infection, pregnancy or lactation, carbapenem",
      "allergy, or concomitant sodium valproate were excluded."
    ),
    dose_range       = "1.0 g meropenem IV every 8 h as a prolonged infusion over 2 h or 2.5 h; dose adjustment permitted at the treating physician's discretion based on TDM",
    regions          = "China (single tertiary-hospital general ICU in Xiangtan, Hunan)",
    crcl_range       = "11.9-136 mL/min (Cockcroft-Gault; cohort mean 47.7 mL/min used as the model centring constant)",
    albumin_range    = "17.7-44.4 g/L (mean 31.9)",
    notes            = paste(
      "Prospective therapeutic drug monitoring cohort enrolled between",
      "March 2023 and January 2025 (Zou 2026 Section 2.1). Samples were",
      "drawn at steady state after at least four consecutive doses, at",
      "30 min before dosing and at 2 h (or 2.5 h), 4 h and 6 h after the",
      "start of infusion; approximately three sampling times per patient.",
      "Observed meropenem concentrations: mean 16.183 mg/L, range",
      "0.477-71.507 mg/L (Table 1). Assay: validated HPLC-UV at 299 nm,",
      "calibrated 0.5-80.0 mg/L, intra-batch precision 4.95-10.07% and",
      "inter-batch precision 4.94-9.07% (Section 2.4). Model fitted in",
      "Phoenix NLME 8.1 using FOCE-ELS; OFV 685.07 (Table 2). Internal",
      "validation by 1,000-replicate pcVPC and 500-replicate bootstrap.",
      "Protein binding of meropenem is approximately 2%, so the authors",
      "treated total plasma concentrations as free concentrations in the",
      "Monte Carlo simulations (Section 2.6); the model is therefore",
      "equally a total- and a free-concentration model."
    )
  )

  ini({
    # ===== Structural PK (Zou 2026 Table 2 'Final model' column) =====
    # Reference subject: postoperative Cockcroft-Gault CRCL = 47.7 mL/min.
    lcl <- log(6.472); label("Typical clearance at CRCL = 47.7 mL/min (L/h)")  # Zou 2026 Table 2 final model: CL = 6.472 L/h, RSE 8.53% (bootstrap median 6.46, 95% CI 5.52-7.81)
    lvc <- log(26.69); label("Typical central volume of distribution (L)")     # Zou 2026 Table 2 final model: Vc = 26.69 L, RSE 13.55% (bootstrap median 26.6, 95% CI 20.3-38.1)

    # ===== Covariate effect =====
    # Power effect of postoperative Cockcroft-Gault creatinine clearance on
    # CL, centred at the cohort mean 47.7 mL/min (Zou 2026 Section 3.2):
    #   CL = 6.472 * (CrCL/47.7)^0.3834 * exp(eta_CL)
    e_crcl_cl   <- 0.3834;      label("Power exponent of postoperative CRCL on CL (unitless)")  # Zou 2026 Table 2 final model, "CRCL on CL" = 0.3834, RSE 16.54% (bootstrap median 0.396, 95% CI 0.254-0.545)
    crcl_ref_cl <- fixed(47.7); label("Reference postoperative Cockcroft-Gault CRCL for the power CL covariate (mL/min, cohort mean)")  # Zou 2026 Section 3.2 covariate equation denominator; equals the Table 1 cohort mean creatinine clearance

    # ===== Inter-individual variability (Zou 2026 Table 2 final model) =====
    # Table 2 reports omega as a percentage and its footnote defines it as
    # the "square root of inter-individual variance", i.e. omega itself (the
    # SD of the eta) rather than a %CV requiring the log(1 + CV^2) transform.
    # Taken at face value: omega_CL = 0.494, omega_Vc = 0.8199, so
    #   var(etalcl) = 0.494^2   = 0.244036
    #   var(etalvc) = 0.8199^2  = 0.672236
    # The reported correlation rho = 0.94 gives the off-diagonal covariance
    #   cov = 0.94 * 0.494 * 0.8199 = 0.380727
    # The resulting 2x2 matrix is positive definite (eigenvalues 0.8949 and
    # 0.0213), but is strongly correlated; see the vignette for the check.
    etalcl + etalvc ~ c(0.244036,
                        0.380727, 0.672236)

    # ===== Residual error (Zou 2026 Table 2 final model) =====
    # Additive residual error, SD = 4.17 mg/L.
    #
    # Table 2's residual row is labelled "Additional error", unit mg/L, with
    # a final-model value of 17.4 and a bootstrap median of 4.64 (95% CI
    # 2.99-6.34). Those two columns cannot both be standard deviations: 17.4
    # lies far outside the bootstrap CI, contradicting Section 3.3's claim
    # that "all final parameter values fell within the 95% confidence
    # intervals". They reconcile exactly if the final-model column reports
    # the residual VARIANCE and the bootstrap column reports the SD, since
    # sqrt(17.4) = 4.17 mg/L, which is inside 2.99-6.34.
    #
    # Section 2.5/3.2 prose calls the retained model "proportional", but
    # digitising the 102 resolvable markers of Figure 1b (observed vs
    # individual-predicted) shows the residual SD is flat in concentration
    # (3.1, 3.5, 6.4, 5.2, 3.2, 2.4 and 0.7 mg/L in IPRED bins 0-5, 5-10,
    # 10-15, 15-20, 20-30, 30-45 and 45-75 mg/L) while the RELATIVE SD falls
    # monotonically from 84% to 1.3%. That is an additive error; a
    # proportional error of any magnitude is refuted at both ends of the
    # range. The overall digitised residual SD is 4.44 mg/L, bracketing
    # sqrt(17.4) = 4.17 and the bootstrap median 4.64. The table label
    # ("Additional" = additive) and the unit (mg/L) agree. Full reasoning
    # and the digitisation are in the vignette Errata.
    addSd <- 4.17; label("Additive residual error (mg/L)")  # Zou 2026 Table 2 final model: residual variance 17.4 (mg/L)^2, RSE 12.0%; SD = sqrt(17.4) = 4.17 mg/L (bootstrap median SD 4.64, 95% CI 2.99-6.34)
  })

  model({
    # ----- Individual PK parameters -----
    # Power effect of postoperative Cockcroft-Gault CRCL on CL, centred at
    # the cohort mean 47.7 mL/min (Zou 2026 Section 3.2 equation).
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref_cl)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Meropenem is given as an IV infusion directly into the central
    # compartment; one compartment with linear elimination (Section 3.2).
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Plasma meropenem concentration: dose in mg, vc in L -> mg/L. Protein
    # binding is approximately 2%, so total and free concentrations were
    # treated as interchangeable by the authors (Section 2.6).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
