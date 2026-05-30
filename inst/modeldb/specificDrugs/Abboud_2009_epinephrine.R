Abboud_2009_epinephrine <- function() {
  description <- "One-compartment population PK model for intravenous epinephrine (adrenaline) infusion in adults with septic shock, with a constant endogenous epinephrine production rate (R0) feeding the central compartment and body weight and SAPS II severity score as power covariates on clearance (Abboud 2009)."
  reference <- paste(
    "Abboud I, Lerolle N, Urien S, Tadie JM, Leviel F, Fagon JY, Faisy C.",
    "Pharmacokinetics of epinephrine in patients with septic shock:",
    "modelization and interaction with endogenous neurohormonal status.",
    "Crit Care 2009;13(4):R120. doi:10.1186/cc7972.",
    sep = " "
  )
  vignette <- "Abboud_2009_epinephrine"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )
  # Conversion: epinephrine molecular weight 183.2 g/mol, so 1 mg = 5460 nmol
  # (used to convert the paper's mg/h and ug/kg/min rates to nmol/h for use
  # with this model; the conversion is recorded in the validation vignette).

  covariateData <- list(
    WT = list(
      description        = "Body weight at study inclusion",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL: CL_i = CL * (WT/70)^e_wt_cl with reference 70 kg per Abboud 2009 Table 3 footnote; the body-weight exponent is estimated at 0.60 (not fixed at the canonical allometric 0.75).",
      source_name        = "BW"
    ),
    SAPS_II = list(
      description        = "New Simplified Acute Physiology Score II at ICU admission",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL: CL_i = CL * (SAPS_II/50)^e_saps_ii_cl with reference 50 units per Abboud 2009 Table 3 footnote; the exponent is -0.67 so higher disease severity is associated with lower epinephrine clearance.",
      source_name        = "SAPS II"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 1L,
    age_range      = "Adults (mean 64 +/- 15 years)",
    age_median     = "mean 64 years",
    weight_range   = "Adults (mean 68 +/- 19 kg)",
    weight_median  = "mean 68 kg",
    sex_female_pct = 34.3,
    saps_ii_mean   = "64 +/- 23",
    disease_state  = "Adults with septic shock (mean arterial pressure below 65 mmHg refractory to fluid resuscitation; SAPS II 64 +/- 23) requiring vasopressor therapy in a tertiary medical intensive care unit; epinephrine was the exclusive first-line catecholamine.",
    dose_range     = "Intravenous epinephrine infusion started at 0.15 ug/kg/min and titrated to mean arterial pressure 65-75 mmHg. Observed infusion rate at the per-protocol C1 sampling time (cumulative dose 0.15 mg/kg) had a median of 2 mg/h (range 0.1-7) = 0.52 ug/kg/min (range 0.026-1.7) per Table 2; expressed in molar units that is ~10,090 nmol/h (range 545-38,208).",
    regions        = "France (Paris)",
    notes          = "Single-center prospective study at Hopital Europeen Georges Pompidou, Assistance Publique - Hopitaux de Paris (January-June 2006). 73 plasma epinephrine concentrations across 38 patients: a baseline sample C0 within 15 minutes before infusion onset plus one steady-state sample C1 at cumulative epinephrine dose 0.15 mg/kg (median delay from start 415 minutes, range 90-1260); three C1 samples were lost. Concentrations measured by HPLC with coulometric detection, limit of quantification 0.10 nmol/L. Exclusions: pregnancy, renal replacement therapy during the study, catecholamines in the 24 hours preceding enrolment. Hydrocortisone and recombinant human activated protein C were not used. ICU mortality 65.7%. Causes of septic shock per Table 1: community-acquired pneumonia 10, nosocomial pneumonia 12, mediastinitis 4, intra-abdominal infection 6, other 4, not documented 2. Baseline demographics in Table 1; baseline hemodynamics and plasma hormones in Table 2."
  )

  ini({
    # Structural parameters from Abboud 2009 Table 3 final population PK
    # estimates. Reference subject: 70 kg body weight, SAPS II = 50 units.
    # Table 3 reports CL in L/h, V in L, and R0 (baseline endogenous
    # epinephrine input rate) in nmol/h.
    lcl <- log(127);  label("Clearance for reference subject (CL, L/h)")              # Table 3 (127 L/h at 70 kg BW, SAPS II 50)
    lvc <- log(7.9);  label("Volume of distribution (V, L)")                          # Table 3 (final-model V = 7.9 L)
    lrbase <- log(43.5); label("Endogenous epinephrine production rate (R0, nmol/h)")    # Table 3

    # Covariate effects from Abboud 2009 Table 3 (final relationship for CL
    # in the Results section): CL_i = 127 * (BW/70)^0.60 * (SAPS II/50)^-0.67.
    e_wt_cl      <-  0.60; label("Power exponent of (WT/70) on CL (theta_BW, unitless)")          # Table 3
    e_saps_ii_cl <- -0.67; label("Power exponent of (SAPS_II/50) on CL (theta_SAPSII, unitless)") # Table 3

    # Inter-individual variability (lognormal, exponential per Methods:
    # "BSVs were assumed to be exponential"). Table 3 reports the square
    # root of omega^2 for each term; variance = SD^2. No covariance was
    # retained in the final model (Methods: "if the correlation between
    # terms was low, it was fixed at 0"), so the two etas are independent.
    etalcl ~ 0.1089  # Table 3 BSV(CL) = 0.33; 0.33^2 = 0.1089
    etalrbase ~ 1.5129  # Table 3 BSV(R0) = 1.23; 1.23^2 = 1.5129

    # Residual error: both components held fixed at the HPLC assay-
    # quantification values per Methods ("the residual variability
    # parameters were fixed as follows: 10% and 0.1 nmol/L for the
    # proportional and additive components, according to the assay
    # quantification"). Table 3 footnote "b" marks these as fixed values.
    propSd <- fixed(0.10); label("Proportional residual SD (fraction)")  # Table 3 (fixed)
    addSd  <- fixed(0.10); label("Additive residual SD (nmol/L)")        # Table 3 (fixed)
  })
  model({
    # Individual PK parameters. Reference subject: WT = 70 kg, SAPS_II = 50.
    # Final CL relationship per Abboud 2009 Results:
    #   CL_i = 127 * (WT/70)^0.60 * (SAPS_II/50)^-0.67  (L/h)
    # IIV is on CL and on the endogenous production rate R0; V has no IIV
    # (the Methods retained only omega^2_CL initially, and only BSV for CL and
    # R0 could be accurately estimated once the residual error was fixed).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (SAPS_II / 50)^e_saps_ii_cl
    vc <- exp(lvc)
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc

    # One-compartment open model with first-order elimination. Exogenous
    # epinephrine is administered as an IV infusion (cmt = central, rate
    # column on the dose row); rxode2 adds the event-driven infusion rate
    # to d/dt(central) automatically. The constant endogenous epinephrine
    # production rbase (nmol/h) is added as a baseline input so the
    # steady-state plateau equals (rate + R0) / CL, matching the explicit
    # formula in Abboud 2009 Results:
    #   C_plateau (nmol/L) = (rate + R0)
    #                        / ( 127 * (BW/70)^0.60 * (SAPS_II/50)^-0.67 )
    d/dt(central) <- rbase - kel * central

    # Initial condition: pre-infusion endogenous steady state, central(0) =
    # R0 / kel = R0 * V / CL, so the observed baseline Cc(0) = R0 / CL
    # reproduces the C0 sample drawn within 15 minutes before infusion
    # onset. For a typical subject this yields 43.5 / 127 = 0.343 nmol/L,
    # which matches the observed median baseline epinephrine of
    # 0.34 nmol/L in Table 2.
    central(0) <- rbase / kel

    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
