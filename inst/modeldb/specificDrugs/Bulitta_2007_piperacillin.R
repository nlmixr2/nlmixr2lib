Bulitta_2007_piperacillin <- function() {
  description <- "Two-compartment first-order IV population PK model for piperacillin in 8 adult cystic-fibrosis patients and 26 adult healthy volunteers receiving 4 g piperacillin as a 5-min intravenous infusion (Bulitta 2007). Lean body mass (LBM) is the size descriptor with allometric scaling (exponents 0.75 on CL and Q, 1.0 on V1 and V2; LBM_STD = 53 kg). A cystic-fibrosis disease-state indicator multiplicatively scales V1 and V2 via fcyf_vss^DIS_CF (fcyf_vss = 0.926), with fcyf_cl^DIS_CF retained on CL at its boundary estimate of 1.00 for model-form traceability."
  reference <- paste(
    "Bulitta JB, Duffull SB, Kinzig-Schippers M, Holzgrabe U, Stephan U,",
    "Drusano GL, Sorgel F. Systematic comparison of the population",
    "pharmacokinetics and pharmacodynamics of piperacillin in cystic fibrosis",
    "patients and healthy volunteers. Antimicrob Agents Chemother.",
    "2007;51(7):2497-2507. doi:10.1128/AAC.01477-06."
  )
  vignette <- "Bulitta_2007_piperacillin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Computed by the Cheymol / James formulas",
        "(Bulitta 2007 Materials and Methods 'Subjects', citing refs 9 and 33).",
        "Allometric scaling on CL and Q (exponent 0.75 fixed) and on V1 and V2",
        "(exponent 1.0 fixed) with reference LBM_STD = 53 kg (Bulitta 2007",
        "Materials and Methods 'Size models')."
      ),
      source_name        = "LBM"
    ),
    DIS_CF = list(
      description        = "Cystic-fibrosis disease-state indicator (1 = CF patient, 0 = healthy participant)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Time-fixed per subject. Multiplicative power-form effect on V1 and V2",
        "via fcyf_vss^DIS_CF (Bulitta 2007 Table 4 LBM-allometric row:",
        "FCYF_VSS = 0.926, 90% CI 0.82-1.02) and on CL via fcyf_cl^DIS_CF",
        "(Bulitta 2007 Table 4 LBM-allometric row: FCYF_CL = 1.00, 90% CI",
        "0.92-1.09; estimated at the model boundary, retained for model-form",
        "traceability with no numerical effect at the published point estimate)."
      ),
      source_name        = "CF"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "CF patients 21 +/- 4 years (mean +/- SD); healthy participants 25 +/- 4 years",
    weight_range   = "CF patients 43.1 +/- 7.8 kg total body weight; healthy participants 71.1 +/- 11.8 kg",
    weight_median  = NULL,
    sex_female_pct = 47.1,
    race_ethnicity = "Caucasian (all 34 subjects)",
    disease_state  = "8 adult cystic-fibrosis patients (CFTR-confirmed by sweat test and clinical history; studied during an infection-free period) and 26 adult healthy volunteers (control group)",
    dose_range     = "Single 4 g piperacillin as a 5-min IV infusion (one CF patient received 3 g)",
    regions        = "Germany (single centre, Nurnberg-Heroldsberg)",
    lbm_range      = "CF patients 37.2 +/- 6.9 kg LBM; healthy participants 56.4 +/- 7.2 kg LBM",
    bsa_range      = "CF patients 1.38 +/- 0.18 m^2; healthy participants 1.85 +/- 0.18 m^2",
    bmi_range      = "CF patients 16.7 +/- 1.1 kg/m^2; healthy participants 23.6 +/- 3.7 kg/m^2",
    height_range   = "CF patients 160 +/- 11.9 cm; healthy participants 174 +/- 8.4 cm",
    notes          = paste(
      "Open-label single-dose parallel-group study. Plasma samples collected",
      "before infusion, at end of infusion (5 min), and at frequent intervals",
      "to 24 h post-infusion (21 nominal time points). Piperacillin assayed by",
      "HPLC with mezlocillin internal standard (LOQ 0.200 mg/L; linearity",
      "0.200-150 mg/L; interday precision 3.5-9.2 %). Protein binding fixed at",
      "30 % for both groups in the Monte Carlo PD simulations. Demographics",
      "from Bulitta 2007 Table 1; structural parameters from Table 3 (final",
      "two-compartment LBM-allometric model); disease-state scale factors from",
      "Table 4 (LBM-allometric row)."
    )
  )

  ini({
    # ===== Structural PK (Bulitta 2007 Table 3, LBM-allometric final model) =====
    # Typical-value reference subject: healthy adult with LBM = LBM_STD = 53 kg.
    # CL and V1 / V2 columns for "CF patients" in Table 3 are derived from
    # these healthy-typical estimates by the disease-state factors below.
    lcl <- log(11.3); label("Typical CL at LBM = 53 kg, healthy adult (L/h)")           # Bulitta 2007 Table 3: CL = 11.3 L/h (90% CI 10.5-12.3 CF, 10.9-11.7 healthy bootstrap)
    lvc <- log(7.01); label("Typical V1 (central) at LBM = 53 kg, healthy adult (L)")    # Bulitta 2007 Table 3: V1 healthy = 7.01 L (90% CI 6.36-8.02)
    lvp <- log(3.37); label("Typical V2 (peripheral) at LBM = 53 kg, healthy adult (L)") # Bulitta 2007 Table 3: V2 healthy = 3.37 L (90% CI 2.74-3.86)
    lq  <- log(12.8); label("Typical inter-compartmental clearance at LBM = 53 kg (L/h)") # Bulitta 2007 Table 3: CLic = 12.8 L/h (90% CI 7.31-15.4)

    # Disease-state scale factors (Bulitta 2007 Table 4, LBM-allometric row).
    # fcyf_cl was estimated at the boundary 1.00 (no CF effect at the
    # published point estimate); retained for model-form traceability so the
    # encoded equations match the paper's Eq. CL_i / V_i exactly.
    fcyf_cl  <- 1.00;  label("CF / healthy ratio on CL (unitless)")                       # Bulitta 2007 Table 4 LBM-allometric: FCYF_CL = 1.00 (90% CI 0.92-1.09)
    fcyf_vss <- 0.926; label("CF / healthy ratio on V1 and V2 (unitless)")                # Bulitta 2007 Table 4 LBM-allometric: FCYF_VSS = 0.926 (90% CI 0.82-1.02)

    # Allometric exponents (Bulitta 2007 Materials and Methods 'Size models':
    # "The allometric exponent was fixed at 1.0 for all volume terms and at
    # 0.75 for all clearance terms.").
    e_lbm_cl <- fixed(0.75); label("Allometric exponent on CL (unitless, fixed)")        # Bulitta 2007 Methods: 0.75 fixed
    e_lbm_q  <- fixed(0.75); label("Allometric exponent on Q (unitless, fixed)")         # Bulitta 2007 Methods: 0.75 fixed
    e_lbm_vc <- fixed(1.00); label("Allometric exponent on V1 (unitless, fixed)")        # Bulitta 2007 Methods: 1.0 fixed
    e_lbm_vp <- fixed(1.00); label("Allometric exponent on V2 (unitless, fixed)")        # Bulitta 2007 Methods: 1.0 fixed

    # Reference LBM (Bulitta 2007 Materials and Methods 'Size models':
    # "a standard lean body mass (LBMSTD) of 53 kg").
    lbm_std <- fixed(53); label("Reference lean body mass (kg, fixed)")                  # Bulitta 2007 Methods: LBM_STD = 53 kg

    # ===== Between-subject variability (Bulitta 2007 Table 3, log-normal IIV) =====
    # The paper reports the apparent CV on a log scale; the variance is
    # related to CV via omega^2 = log(1 + CV^2):
    #   10.4% CV -> log(1 + 0.104^2) = 0.010758       (CL)
    #   26.0% CV -> log(1 + 0.260^2) = 0.065396       (V1)
    #   34.2% CV -> log(1 + 0.342^2) = 0.110667       (V2)
    # Off-diagonal: cov(V1, V2) = r * sd(V1) * sd(V2)
    #   r = -0.80 (Table 3 footnote e); cov = -0.80 * sqrt(0.065396 * 0.110667)
    #   = -0.068061. Block specified as c(var_V1, cov_V1V2, var_V2).
    etalcl ~ 0.010758                                                                     # Bulitta 2007 Table 3: BSV CL = 10.4% CV
    etalvc + etalvp ~ c(0.065396,
                       -0.068061, 0.110667)                                               # Bulitta 2007 Table 3: BSV V1 = 26.0% CV, BSV V2 = 34.2% CV, r(V1,V2) = -0.80

    # ===== Residual error (Bulitta 2007 Table 3: combined proportional + additive) =====
    propSd <- 0.132; label("Proportional residual error (fraction)")                      # Bulitta 2007 Table 3: CVC = 13.2%
    addSd  <- 1.88;  label("Additive residual error (mg/L)")                              # Bulitta 2007 Table 3: SDC = 1.88 mg/L
  })

  model({
    # ----- Allometric and disease-state factors -----
    size_cl <- (LBM / lbm_std)^e_lbm_cl
    size_q  <- (LBM / lbm_std)^e_lbm_q
    size_vc <- (LBM / lbm_std)^e_lbm_vc
    size_vp <- (LBM / lbm_std)^e_lbm_vp

    # ----- Individual PK parameters (Bulitta 2007 Eq. CL_i / V_i with FCYF) -----
    cl <- exp(lcl + etalcl) * size_cl * fcyf_cl^DIS_CF
    vc <- exp(lvc + etalvc) * size_vc * fcyf_vss^DIS_CF
    vp <- exp(lvp + etalvp) * size_vp * fcyf_vss^DIS_CF
    q  <- exp(lq) * size_q

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- Two-compartment IV ODE system -----
    # Piperacillin dosed directly into the central compartment (IV infusion;
    # the 5-min zero-order input from the source study is supplied via the
    # event table, not as a structural parameter).
    d/dt(central)      <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-                   k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Plasma piperacillin concentration. Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
