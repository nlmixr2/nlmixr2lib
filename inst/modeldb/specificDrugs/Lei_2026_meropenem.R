Lei_2026_meropenem <- function() {
  description <- "Two-compartment IV population PK model for meropenem in 28 Chinese adults with sepsis or septic shock (Lei 2026). Clearance carries a power effect of Cockcroft-Gault creatinine clearance and a piecewise-linear (hinge) decrease with baseline plasma 4-pyridoxic acid, an endogenous OAT1/OAT3 biomarker, above an estimated knot of 26.5 ng/mL. Central volume carries a power effect of serum albumin and a multiplicative renal-status factor (acute kidney injury, continuous renal replacement therapy). Combined proportional and additive residual error."
  reference <- paste(
    "Lei Z, Liang H, Cheng Q, Kong W, Du Y, Yao X, Feng F, Jin Y, Wang W,",
    "Li H, Lu M, Liu D, Shen N.",
    "Pyridoxic Acid as a Novel Endogenous Biomarker for OAT-Mediated",
    "Clearance of Meropenem: A Population Pharmacokinetic Modeling and",
    "PK/PD-Response Analysis in Patients with Sepsis.",
    "Pharmaceutics. 2026;18(9):1113.",
    "doi:10.3390/pharmaceutics18091113.",
    sep = " "
  )
  vignette <- "Lei_2026_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cockcroft-Gault creatinine clearance at baseline (Lei 2026 Table 1",
        "footnote 2; Section 3.2.2). The paper's covariate equation (Eq. 3)",
        "uses CLcr in L/h centred at the PopPK-cohort median 2.62 L/h",
        "(Section 3.2.3: 'the constants 2.62 L/h and 31.4 g/L correspond to",
        "the median values of CLcr and ALB observed in the PopPK study",
        "population'); the model converts the canonical mL/min column with",
        "CRCL * 60 / 1000 before centring. 2.62 L/h = 43.7 mL/min.",
        "Full-cohort (n = 49) median 47.1 mL/min, IQR 32.8-80.8 (Table 1).",
        "Stored under canonical CRCL with raw mL/min, as in",
        "Zou_2026_meropenem.R and Delattre_2010_amikacin.R."
      ),
      source_name = "CLcr"
    ),
    PYRIDOXIC_ACID_BL = list(
      description = "Baseline plasma 4-pyridoxic acid (PDA) concentration, an endogenous OAT1/OAT3 substrate",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured by LC-MS/MS in the pre-dose (baseline) sample only",
        "(Lei 2026 Section 2.2), so time-fixed per subject. Enters CL as a",
        "hinge: no effect at or below the estimated knot of 26.5 ng/mL, a",
        "linear decrease of 0.042 L/h per ng/mL above it (Eq. 3).",
        "Full-cohort median 6.8 ng/mL, IQR 4.3-17.3 (Table 1); median of",
        "the above-knot subpopulation 49.6 ng/mL and maximum 121.3 ng/mL",
        "(Section 4.1). Because the PDA term is subtracted rather than",
        "multiplied, Eq. 3 gives a non-positive clearance when CLcr is very",
        "low and PDA very high (e.g. below about 13 mL/min at the observed",
        "maximum PDA of 121.3 ng/mL); the model is encoded as published and",
        "is not meant to be simulated in that corner."
      ),
      source_name = "PDA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline serum albumin. Power covariate on Vc centred at the PopPK-cohort median 31.4 g/L (Eq. 4; Section 3.2.3). Full-cohort median 31.0 g/L, IQR 29.1-34.8 (Table 1).",
      source_name = "ALB"
    ),
    DIS_ARF = list(
      description = "Acute kidney injury at baseline (KDIGO criteria), with or without continuous renal replacement therapy",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "One of two indicators decomposing the paper's three-level GROUP",
        "covariate (sepsis with normal renal function / acute kidney injury",
        "/ continuous renal replacement therapy, classified by KDIGO",
        "criteria at baseline; Section 2.1). DIS_ARF = 1 for every subject",
        "with acute kidney injury, including those on CRRT; DIS_ARF = 0",
        "selects the normal-renal-function reference. The levels map as",
        "normal (DIS_ARF = 0, RRT_CRRT_STATUS = 0), AKI without CRRT",
        "(1, 0) and CRRT (1, 1). The combination (0, 1) is not a level of",
        "the source covariate."
      ),
      source_name = "GROUP"
    ),
    RRT_CRRT_STATUS = list(
      description = "Receiving continuous renal replacement therapy at baseline",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Second indicator of the decomposed GROUP covariate (see DIS_ARF); a CRRT subject also carries DIS_ARF = 1. The CRRT level multiplies Vc by 4.08 in Eq. 4, which equals the AKI factor squared (2.02^2 = 4.0804); the model encodes it that way with the single estimated 'Vc-GROUP' parameter of Table 2, as e_dis_arf_vc^(DIS_ARF + RRT_CRRT_STATUS).",
      source_name = "GROUP"
    )
  )

  covariatesDataExcluded <- list(
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Added on Vc in forward inclusion (Supplementary Table S2 step 2) and removed in backward elimination (step 8)."
    ),
    URINE_VOL_24H = list(
      description = "24-hour urine output",
      units = "mL",
      type = "continuous",
      notes = "Added on CL in forward inclusion (Supplementary Table S2 step 5) and removed in backward elimination (step 9)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28L,
    n_studies = 1L,
    n_concentrations = 138L,
    age_range = "full cohort (n = 49): median 74.0 years, IQR 64.0-85.0 (Table 1); adults >= 18 years",
    sex_female_pct = 63.3,
    race_ethnicity = "Han Chinese (100% of the full cohort, Table 1)",
    bmi_range = "full cohort median 21.5 kg/m^2, IQR 18.6-24.7 (Table 1)",
    disease_state = paste(
      "Adults with sepsis or septic shock (Sepsis-3) receiving IV meropenem",
      "in a single-centre ICU (Peking University Third Hospital, Beijing),",
      "January 2022 to December 2024. Full cohort renal status: normal",
      "renal function 21, AKI without CRRT 21, CRRT 7 (KDIGO). Median",
      "APACHE II 17, SOFA 7; 83.7% pulmonary infection (Table 1)."
    ),
    dose_range = "0.5 g q12h (3 patients), 1.0 g q8h (19) or 1.0 g q12h (27), all as 1-h IV infusions (Section 3.1; full cohort of 49)",
    regions = "China",
    crcl_range = "full cohort median 47.1 mL/min, IQR 32.8-80.8 (Cockcroft-Gault); PopPK-subset median 2.62 L/h (43.7 mL/min)",
    albumin_range = "full cohort median 31.0 g/L, IQR 29.1-34.8; PopPK-subset median 31.4 g/L",
    pyridoxic_acid_range = "full cohort median 6.8 ng/mL, IQR 4.3-17.3; maximum 121.3 ng/mL",
    notes = paste(
      "The PopPK model was fitted to the intensively sampled subset of 28",
      "patients (138 samples: pre-dose, end of infusion, 0.5 h post-",
      "infusion, 3.5 h (q8h) or 5.5 h (q12h) post-infusion, and pre-next-",
      "dose) out of a 49-patient exposure-response cohort; demographics",
      "are tabulated only for the 49 (Table 1). NONMEM 7.2, FOCE-I.",
      "Meropenem assay LC-MS/MS, 50.0-5000 ng/mL calibration range;",
      "observed concentrations 0.28-223.12 mg/L (Section 3.1)."
    )
  )

  ini({
    # ===== Structural PK (Lei 2026 Table 2, 'Final Model' column) =====
    # Reference subject: CLcr = 2.62 L/h, baseline PDA <= 26.5 ng/mL,
    # ALB = 31.4 g/L, sepsis with normal renal function.
    lcl <- log(7.21); label("Clearance coefficient at CLcr = 2.62 L/h, before the PDA term (L/h)")  # Table 2 final model CL = 7.21 L/h (95% CI -7.25, 21.68); Eq. 3 leading coefficient
    lvc <- log(7.41); label("Central volume at ALB = 31.4 g/L, normal renal function (L)")         # Table 2 final model Vc = 7.41 L (95% CI 3.29, 11.53)
    lq <- log(6.45); label("Intercompartmental clearance (L/h)")                                   # Table 2 final model Q = 6.45 L/h (95% CI 0.002, 12.90)
    lvp <- log(9.33); label("Peripheral volume (L)")                                               # Table 2 final model Vp = 9.33 L (95% CI 3.43, 15.23)

    # ===== Covariate effects (Lei 2026 Table 2; Eqs. 3 and 4) =====
    e_crcl_cl <- 0.29; label("Power exponent of CLcr on CL (unitless)")                                     # Table 2 'CL-CLcr' = 0.29 (95% CI -0.167, 0.747); Eq. 3 exponent
    crcl_ref_cl <- fixed(2.62); label("Reference CLcr for the CL power term (L/h, PopPK-cohort median)")     # Eq. 3 denominator; Section 3.2.3 'the constants 2.62 L/h ... correspond to the median values of CLcr'
    e_pyridoxic_acid_bl_cl <- -0.042; label("Linear slope of baseline PDA on CL above the knot (L/h per ng/mL)")  # Table 2 'CL-PDA' = -0.042 (95% CI -0.309, 0.225); Eq. 3 '- 0.042 x PDA'
    lpyridoxic_acid_bl_hinge <- log(26.5); label("Baseline PDA knot below which CL does not depend on PDA (ng/mL)")  # Section 3.2.3 'The PDA cutoff value was estimated as a structural parameter ... identified at 26.5 ng/mL'; Eq. 3 branch limit
    e_alb_vc <- -3.67; label("Power exponent of ALB on Vc (unitless)")                                      # Table 2 'Vc-ALB' = -3.67 (95% CI -6.12, -1.22); Eq. 4 exponent
    alb_ref_vc <- fixed(31.4); label("Reference ALB for the Vc power term (g/L, PopPK-cohort median)")       # Eq. 4 denominator; Section 3.2.3 '31.4 g/L correspond to the median values of ... ALB'
    e_dis_arf_vc <- 2.02; label("Multiplicative factor on Vc per renal-status step: AKI = factor, CRRT = factor^2 (unitless)")  # Table 2 'Vc-GROUP' = 2.02 (95% CI 0.922, 3.118); Eq. 4 AKI row x 2.02, CRRT row x 4.08 = 2.02^2

    # ===== Inter-individual variability (Lei 2026 Table 2, final model) =====
    # Table 2 prints IIV as a percentage. Section 3.2.3 quotes the CL and Vc
    # variances as 0.052 and 0.069, which match sqrt(omega^2) = 22.7% and
    # 26.4% (the text attaches them to the opposite parameters; the table
    # and the abstract both give CL = 26.4%). The percentages are therefore
    # sqrt(omega^2) x 100, and omega^2 = (IIV%/100)^2.
    etalcl ~ 0.069696  # Table 2 'IIV on CL' = 26.4%, shrinkage 2.6%: 0.264^2
    etalvc ~ 0.051529  # Table 2 'IIV on Vc' = 22.7%, shrinkage 45.5%: 0.227^2
    etalq ~ 1.179396   # Table 2 'IIV on Q' = 108.6%, shrinkage 29.2%: 1.086^2

    # ===== Residual error (Lei 2026 Table 2, final model) =====
    # Table 2 prints the raw NONMEM SIGMA variances (0.0578 proportional,
    # 0.105 additive), not SDs; see the vignette Assumptions section for
    # why (Supplementary Figure S1A scatter is consistent with a 24%
    # proportional SD and not with 5.8%).
    propSd <- 0.2404163; label("Proportional residual error (fraction)")  # Table 2 'Proportional' = 0.0578 (variance): sqrt(0.0578)
    addSd <- 0.3240370; label("Additive residual error (mg/L)")           # Table 2 'Additive' = 0.105 (variance, mg/L^2): sqrt(0.105)
  })

  model({
    # ----- Covariate terms -----
    # Eq. 3 is written with CLcr in L/h; the canonical CRCL column is mL/min.
    crcl_lh <- CRCL * 60 / 1000
    pyridoxic_acid_bl_hinge <- exp(lpyridoxic_acid_bl_hinge)

    # ----- Individual PK parameters -----
    # Eq. 3, both branches in one line:
    #   PDA <= 26.5: CL = 7.21 * (CLcr/2.62)^0.29 - 0.042 * 26.5
    #   PDA >  26.5: CL = 7.21 * (CLcr/2.62)^0.29 - 0.042 * PDA
    # The flat arm subtracts the slope times the knot, so the two arms meet
    # at the knot. Exponential IIV (Section 2.3) multiplies the whole
    # typical value.
    cl <- (exp(lcl) * (crcl_lh / crcl_ref_cl)^e_crcl_cl +
      e_pyridoxic_acid_bl_cl * max(PYRIDOXIC_ACID_BL, pyridoxic_acid_bl_hinge)) * exp(etalcl)
    # Eq. 4: renal-status factor 1 (normal), 2.02 (AKI), 4.08 = 2.02^2 (CRRT,
    # where DIS_ARF = 1 and RRT_CRRT_STATUS = 1).
    vc <- exp(lvc + etalvc) * (ALB / alb_ref_vc)^e_alb_vc *
      e_dis_arf_vc^(DIS_ARF + RRT_CRRT_STATUS)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system (Eqs. 1-2; IV infusion into central) -----
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
