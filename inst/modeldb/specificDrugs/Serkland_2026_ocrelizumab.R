Serkland_2026_ocrelizumab <- function() {
  description <- "Two-compartment population PK of intravenous ocrelizumab with linear clearance, coupled to an indirect-response CD19+ B-lymphocyte model in which ocrelizumab removes B cells by a second-order (concentration x count) kill term, in adults with newly diagnosed relapsing-remitting multiple sclerosis given repeated 600 mg infusions"
  reference <- paste(
    "Serkland TT, Oma AKF, Hallin EI, Overas MH, Nygaard GO, Myhr KM,",
    "Torkildsen O, Roblitz S, Skrede S (2026).",
    "Pharmacokinetic-Pharmacodynamic Modelling of Repeated Ocrelizumab Dosing",
    "in Relapsing-Remitting Multiple Sclerosis.",
    "Clin Pharmacokinet 65(10):1541-1552. doi:10.1007/s40262-026-01692-8.",
    sep = " "
  )
  vignette <- "Serkland_2026_ocrelizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")
  dosing <- "central"

  # Final model ("model 1 without TMDD"; Serkland 2026 Table 2 and Supplementary
  # Table S1 middle column). PK from Eq. 5-6 with the target-mediated term
  # removed (Results: "simplified ... by removing the target-mediated
  # elimination from the PK model"); B cells from Eq. 12 (Eq. 7 with
  # kin = B0 * kout):
  #   dC/dt = -Q/V1 * C + Q/V1 * P - CL/V1 * C
  #   dP/dt =  Q/V2 * C - Q/V2 * P
  #   dB/dt =  B0 * kout - kdeg * C * B - kout * B,   B(0) = B0
  # C and P are concentrations (mg/L); the amount form used below is the same
  # system multiplied through by V1 and V2.

  compartmentData <- list(
    central = list(analyte = "ocrelizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ocrelizumab", units = "mg", specimen = "tissue", verified = TRUE),
    bcell = list(analyte = "CD19+ B-lymphocytes", units = "cells/mm^3", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    BLBCELL = list(
      description = "Observed baseline (pre-infusion) CD19+ B-lymphocyte count",
      units = "cells/mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not a covariate effect: the individual observed baseline count IS the",
        "B-cell baseline of the final model (Serkland 2026 Eq. 10, B0,i = B0,m,",
        "with no random effect). It sets the initial condition bcell(0) and the",
        "zero-order production rate kin = BLBCELL * kout (Eq. 12). Measured by",
        "flow cytometry (BD Multitest 6-color TBNK) and reported per mm^3",
        "(= per uL). Cohort mean 297, range 127-726 cells/mm^3 (Table 1). For",
        "the two patients with a missing baseline the authors used the",
        "population median.",
        sep = " "
      ),
      source_name = "B0,m"
    )
  )

  # Covariates tested by Monolix stepwise forward addition / backward deletion
  # on all parameters; none was retained ("not supported as it introduced
  # large parameter uncertainties", Results).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on all parameters (Methods, Sect. 2.4); not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on all parameters (Methods, Sect. 2.4); not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on all parameters (Methods, Sect. 2.4); not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened on all parameters (Methods, Sect. 2.4); not retained. Cohort 8 F / 3 M."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11,
    n_studies = 1,
    age_range = "22-50 years",
    age_mean = "36 years",
    weight_range = "61-96 kg",
    weight_mean = "73 kg",
    height_range = "156-187 cm",
    height_mean = "169 cm",
    sex_female_pct = 72.7,
    disease_state = "newly diagnosed, treatment-naive relapsing-remitting multiple sclerosis",
    dose_range = "600 mg IV infusion (about 4 h) every 6 months; first and second infusions analysed",
    regions = "Norway (OVERLORD-MS trial, ROS-MS sub-study)",
    baseline_cd19_count = "mean 297 cells/mm^3 (range 127-726)",
    n_observations = "114 serum ocrelizumab concentrations plus CD19+ lymphocyte counts over 0-24 weeks after each of two infusions",
    notes = paste(
      "Baseline demographics from Serkland 2026 Table 1. Ocrelizumab was",
      "measured by LC-MS/MS (LLOQ 4 ug/mL, ULOQ 256 ug/mL); CD19+ counts below",
      "the 7 per mm^3 quantification threshold that were reported as zero were",
      "treated as censored (M3) in the fit. Race was not reported.",
      sep = " "
    )
  )

  ini({
    # --- PK (Table 2) -------------------------------------------------------
    lvc <- log(2.58); label("Central volume of distribution V1 (L)") # Table 2 'V1 (L)' 2.58 (RSE 4.49%)
    lvp <- log(1.95); label("Peripheral volume of distribution V2 (L)") # Table 2 'V2 (L)' 1.95 (RSE 5.54%)
    lcl <- log(0.153); label("Target-independent clearance CL (L/day)") # Table 2 'CL (L/day)' 0.153 (RSE 7.21%)
    lq <- log(0.220); label("Intercompartmental clearance Q (L/day)") # Table 2 'Q (L/day)' 0.220 (RSE 15.7%)

    # --- CD19+ B-lymphocyte turnover (Table 2) ------------------------------
    lkout <- fixed(log(0.02)); label("First-order CD19+ lymphocyte loss rate constant kout (1/day)") # Table 2 'kout (1/day)' 0.02 (fixed), taken from Pan et al. (Results)
    lkirrev <- log(1.70); label("Second-order ocrelizumab-dependent CD19+ lymphocyte elimination rate constant kdeg (L/(day*mg))") # Table 2 'kdeg (L/(day mg))' 1.70 (RSE 32.6%)

    # --- IIV: Table 2 reports the SD of the log-normal random effects; the
    # variances below are omega^2.
    etalvc ~ 0.019044 # Table 2 'omega V1' 0.138 (SD); 0.138^2
    etalvp ~ 0.011881 # Table 2 'omega V2' 0.109 (SD); 0.109^2
    etalcl ~ 0.054756 # Table 2 'omega CL' 0.234 (SD); 0.234^2

    # --- Residual error (Table 2; Supplementary Material 2 Eq. 1:
    # y = C + (sigma_add + sigma_prop * C) * e, i.e. Monolix combined1)
    addSd <- 0.545; label("Additive residual error on ocrelizumab concentration (mg/L)") # Table 2 'sigma add,OCR' 0.545 (RSE 17.4%)
    propSd <- 0.0662; label("Proportional residual error on ocrelizumab concentration (fraction)") # Table 2 'sigma prop,OCR' 0.0662 (RSE 13.5%)
    addSd_Bcell <- 2.28; label("Additive residual error on CD19+ lymphocyte count (cells/mm^3)") # Table 2 'sigma add,CD19+' 2.28 (RSE 14.6%)
  })

  model({
    # Individual parameters
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    kout <- exp(lkout)
    kirrev <- exp(lkirrev)

    # Steady-state production from the observed baseline (Eq. 12)
    kin <- BLBCELL * kout

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # PK: Eq. 5-6 without the target-mediated term, in amounts (mg)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc # mg/L

    # CD19+ B lymphocytes (Eq. 12); kirrev (L/(day*mg)) x Cc (mg/L) is 1/day
    d/dt(bcell) <- kin - kirrev * Cc * bcell - kout * bcell
    bcell(0) <- BLBCELL

    Bcell <- bcell # CD19+ lymphocyte count, cells/mm^3

    Cc ~ add(addSd) + prop(propSd) + combined1()
    Bcell ~ add(addSd_Bcell)
  })
}
