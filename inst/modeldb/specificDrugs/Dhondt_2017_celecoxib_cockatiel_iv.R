Dhondt_2017_celecoxib_cockatiel_iv <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order elimination for",
    "celecoxib after a single 10 mg/kg intravenous bolus of an analytical",
    "standard solution (STD) to cockatiels (Nymphicus hollandicus).",
    "Fitted in Phoenix NLME (FOCE-ELS) as the intravenous arm of a",
    "three-drug comparative study; Table 1, 'IV STD' block. Every volume",
    "and clearance term in the source is normalised to body weight, so the",
    "model is coded per kilogram: the dosed amount is ug/kg and the volume",
    "is L/kg, which makes central/vc land directly in ng/mL, the assay",
    "units of Figure 2a. Body weight and sex were screened as covariates",
    "and neither was retained. See Dhondt_2017_celecoxib_cockatiel_oral_std",
    "and Dhondt_2017_celecoxib_cockatiel_oral_cf for the separately fitted",
    "oral arms.",
    sep = " "
  )
  reference <- paste(
    "Dhondt L, Devreese M, Croubels S, De Baere S, Haesendonck R,",
    "Goessens T, Gehring R, De Backer P, Antonissen G.",
    "Comparative population pharmacokinetics and absolute oral",
    "bioavailability of COX-2 selective inhibitors celecoxib, mavacoxib and",
    "meloxicam in cockatiels (Nymphicus hollandicus).",
    "Sci Rep. 2017;7(1):12043.",
    "doi:10.1038/s41598-017-12159-z.",
    sep = " "
  )
  vignette <- "Dhondt_2017_cox2inhibitors_cockatiel"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  # Per-kg normalisation: Table 1 reports Vd in L/kg and the dose in
  # mg/kg, so the state amount is ug/kg and central/vc is ug/L == ng/mL.
  compartmentData <- list(
    central = list(analyte = "celecoxib", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Body weight and sex were formally screened by stepwise forward-backward
  # selection (Methods, 'Pharmacokinetic analysis') and neither reached the
  # p < 0.01 inclusion threshold for any of the three drugs, so neither
  # appears in the final model (Results, 'Pharmacokinetic analysis').
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a continuous covariate by stepwise forward-backward",
        "selection against the -2LL criterion; not significant for any of",
        "the drugs and therefore not retained. Methods 'Pharmacokinetic",
        "analysis'; Results 'Pharmacokinetic analysis'."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "The source screened 'gender' as a categorical covariate and did",
        "not retain it; the source does not state which sex was the",
        "reference category. Cohorts were balanced 17 male / 17 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "91 g",
    weight_range = "91 +/- 10 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 10 mg/kg body weight intravenous bolus into the vena cutanea",
      "ulnaris (wing vein) of a celecoxib analytical standard solution",
      "(5 mg/mL in polyethylene glycol 400:physiological saline, 75:25 v/v)."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The 34 birds (17 male / 17 female) received celecoxib STD both",
      "intravenously and orally in a two-way crossover with a one-month",
      "washout; the two arms were fitted as separate models and the oral",
      "arm is Dhondt_2017_celecoxib_cockatiel_oral_std. A sparse sampling",
      "protocol was used because of the limited blood volume of cockatiels:",
      "sampling times were randomly allocated across birds with a maximum",
      "of two samples per bird, drawn before dosing and at 5, 15, 30 and",
      "45 min and 1, 2, 4, 6, 8 and 12 h. Celecoxib was quantified in",
      "plasma by LC-MS/MS over 5-5000 ng/mL (LOQ 5 ng/mL, LOD 0.22 ng/mL);",
      "values below the LOQ were excluded before fitting. Plasma protein",
      "binding measured separately was 98.98 +/- 0.07% and is not a model",
      "parameter. Immediately after dosing each bird received a 2 mL",
      "intra-crop feed bolus. See Dhondt 2017 Methods 'Animals and",
      "experimental procedure' and 'Celecoxib PK study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 1, 'IV STD' block.
    # Methods equation (2): C(t) = C0 * exp(-Cl/Vd * t), a one-compartment
    # model with first-order elimination.
    #
    # Both values are confirmed by the table's own computed secondary
    # parameters: C0 = D/Vd = 10000/4.67 = 2141 ng/mL (Table 1 reports
    # 2141.54) and AUC(0-inf) = D/Cl = 10000/2.42 = 4132 ng.h/mL (Table 1
    # reports 4140.36).
    # ------------------------------------------------------------------
    lvc <- log(4.67); label("Volume of distribution Vd (log L/kg)") # Table 1, IV STD: Vd = 4.67 L/kg (RSE 8.06%)
    lcl <- log(2.42); label("Total body clearance Cl (log L/h/kg)") # Table 1, IV STD: Cl = 2.42 L/h.kg (RSE 7.93%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 1, 'omega' column. Methods equation (4)
    # declares exponential IIV, P_i = theta_P * exp(eta_Pi), with eta of
    # mean zero and variance omega^2, and states that 'Interindividual
    # variability is reported as omega'.
    #
    # SCALE OF THE omega COLUMN. The table captions call omega 'variance of
    # the interindividual variability', which contradicts Methods equation
    # (4) plus its 'reported as omega' sentence. The equations win: the
    # same captions also call sigma_mult and sigma_add 'variances', yet
    # Methods equation (7) divides one by the other and multiplies the
    # result by a concentration, which is only dimensionally coherent if
    # both are standard deviations. The tabulated second moments are
    # therefore read as SDs throughout and squared here, because nlmixr2
    # omega entries are variances. See the vignette 'Assumptions and
    # deviations' section.
    #   Vd: omega < 0.001 -> encoded at the printed upper bound,
    #       0.001^2 = 1e-06, which keeps OMEGA positive definite.
    #   Cl: omega = 0.0418 -> 0.0418^2 = 0.00174724
    # ------------------------------------------------------------------
    etalvc ~ 1e-06 # Table 1, IV STD omega for Vd, reported as '< 0.001'
    etalcl ~ 0.00174724 # Table 1, IV STD omega for Cl = 0.0418

    # ------------------------------------------------------------------
    # Residual error -- Methods: 'For the IV and PO STD dosing, a
    # multiplicative (equation (5)) and additive (equation (6)) error model
    # was used ... respectively'. Equation (5) is
    # C_obs = C_pred * (1 + epsilon), which is nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.42; label("Proportional residual error (fraction)") # Table 1, IV STD: Res. Error = 0.42 (RSE 12.32%)
  })

  model({
    # 1. Individual parameters (exponential IIV, Methods equation (4)).
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system -- Methods equation (2) written as a differential
    # equation. The intravenous dose enters central directly.
    d/dt(central) <- -kel * central

    # 4. Observation. Amount in ug/kg over an apparent volume in L/kg
    # gives ug/L, which is ng/mL -- the assay units of Figure 2a.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
