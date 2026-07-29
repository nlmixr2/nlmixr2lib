Hall_2017_dapsone <- function() {
  description <- "One-compartment population PK model with first-order oral absorption for dapsone in healthy US adults across a wide weight range; covariate effects on Ka, CL, and Vc are encoded via the published MARS piecewise-linear basis functions of weight, age, and blood urea nitrogen (Hall 2017)."
  reference <- "Hall RG II, Pasipanodya JG, Swancutt MA, Meek C, Leff R, Gumbo T. Supervised machine-learning reveals that old and obese people achieve low dapsone concentrations. CPT Pharmacometrics Syst Pharmacol. 2017;6(8):552-559. doi:10.1002/psp4.12208"
  vignette <- "Hall_2017_dapsone"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives all three MARS hinge models (Ka, CL, Vc). In-cohort range 58-138 kg (Table 1); the MARS hinges have no physical constraint outside that range and the linear-scale typical value is clamped at a small positive floor in model() to guard against extrapolation.",
      source_name        = "Weight"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the MARS hinge model on CL via the mirror basis functions BF3_CL = max(0, AGE - 27) * BF2_CL and BF4_CL = max(0, 27 - AGE) * BF2_CL. In-cohort range 21-77 years (Table 1).",
      source_name        = "Age"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the MARS hinge model on Ka via BF1_Ka = max(0, BUN - 7), which interacts with a weight hinge to drive Ka. In-cohort range 7-28 mg/dL (Table 1).",
      source_name        = "BUN"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35,
    n_studies      = 1,
    age_range      = "21-77 years",
    age_median     = "36 years",
    weight_range   = "58-138 kg",
    weight_median  = "85.1 kg",
    bmi_range      = "23.6-43.9 kg/m^2",
    bmi_median     = "30.2 kg/m^2",
    sex_female_pct = 51,
    race_ethnicity = c(White = 83, NonWhite = 17),
    disease_state  = "Healthy adult volunteers (NCT01165840)",
    dose_range     = "100 mg single oral dose under directly observed therapy",
    regions        = "USA (Texas)",
    notes          = "Stratified recruitment by BMI: 12 normal-weight (BMI < 25), 12 overweight / class I-II obese (BMI 25-40), 12 class III obese (BMI > 40); equal numbers of men and women within each stratum. Eight blood draws per subject over 72 h, sampling times chosen by D-optimal design in ADAPT 5. Baseline demographics in Table 1 of Hall 2017."
  )

  ini({
    # Structural population intercepts from the MARS final model (Hall 2017
    # Table 2, Eqs. 2-4). Each intercept is the MARS typical value when all
    # hinge basis functions for that parameter evaluate to zero (e.g., the Ka
    # intercept applies when WT > 74.8 and (WT > 63.7 or BUN <= 7); the Vc
    # intercept applies when WT <= 69.8). Log-transformed for log-normal IIV.
    lka <- log(3.967);  label("Absorption-rate MARS intercept (Ka, 1/h)")   # Table 2 Eq. (2)
    lcl <- log(2.048);  label("Clearance MARS intercept (CL, L/h)")         # Table 2 Eq. (3)
    lvc <- log(36.648); label("Central volume MARS intercept (Vc, L)")      # Table 2 Eq. (4)

    # MARS hinge-function coefficients. Held FIXED at the Hall 2017 Table 2
    # values; these come from a MARS regression on individual posthoc PK
    # estimates, not from the popPK fit itself.
    e_bf3ka_ka <- fixed(-0.062); label("MARS coef on BF3_Ka -> Ka (1/h per kg*mg/dL)")  # Table 2 Eq. (2)
    e_bf5ka_ka <- fixed( 0.162); label("MARS coef on BF5_Ka -> Ka (1/h per kg)")        # Table 2 Eq. (2)
    e_bf3cl_cl <- fixed(-0.002); label("MARS coef on BF3_CL -> CL (L/h per kg*yr)")     # Table 2 Eq. (3)
    e_bf4cl_cl <- fixed(-0.005); label("MARS coef on BF4_CL -> CL (L/h per kg*yr)")     # Table 2 Eq. (3)
    e_bf8v_vc  <- fixed( 8.894); label("MARS coef on BF8_V -> Vc (L per kg)")           # Table 2 Eq. (4)
    e_bf10v_vc <- fixed(-8.818); label("MARS coef on BF10_V -> Vc (L per kg)")          # Table 2 Eq. (4)

    # IIV: omega^2 = log(CV^2 + 1). Hall 2017 reports base 1-compartment %CVs
    # (94.7% Ka, 37.6% CL, 50.5% Vc) on p. 553 / 554 but does NOT report
    # covariate-adjusted IIV after the MARS hinges are applied. Using the
    # base-model %CV here overstates residual IIV (see vignette Errata).
    etalka ~ 0.6402  # CV 94.7% base model                                  # Results p. 553-554
    etalcl ~ 0.1322  # CV 37.6% base model
    etalvc ~ 0.2271  # CV 50.5% base model

    # Residual error: Hall 2017 does not explicitly report a residual-error
    # model. propSd is approximated from Figure 2b (one-compartment obs vs
    # predicted r^2 = 0.91, observed mean 1.63 mg/L SD 1.03 mg/L), giving
    # sqrt(1 - r^2) * SD / mean ~= 19% CV. Held fixed.
    propSd <- fixed(0.20); label("Proportional residual error (fraction)")  # Fig. 2b derivation
  })
  model({
    # MARS hinge basis functions (Hall 2017 Table 2). Knots at WT 63.7, 69.8,
    # 74.8, 77.2 kg; AGE 27 yr; BUN 7 mg/dL.
    BF1_Ka <- max(0, BUN - 7)
    BF3_Ka <- max(0, 63.7 - WT) * BF1_Ka
    BF5_Ka <- max(0, 74.8 - WT)
    BF2_CL <- max(0, 77.2 - WT)
    BF3_CL <- max(0, AGE - 27) * BF2_CL
    BF4_CL <- max(0, 27 - AGE) * BF2_CL
    BF8_V  <- max(0, WT - 69.8)
    BF10_V <- max(0, WT - 74.8)

    # MARS typical values (linear scale): intercept + additive hinge shifts.
    ka_typ <- exp(lka) + e_bf3ka_ka * BF3_Ka + e_bf5ka_ka * BF5_Ka
    cl_typ <- exp(lcl) + e_bf3cl_cl * BF3_CL + e_bf4cl_cl * BF4_CL
    vc_typ <- exp(lvc) + e_bf8v_vc  * BF8_V  + e_bf10v_vc * BF10_V

    # MARS hinge extrapolation outside the cohort (WT 58-138 kg, AGE 21-77 yr,
    # BUN 7-28 mg/dL) can drive the linear-scale typical value to <= 0, which
    # is non-physical. Clamp at a small positive floor before applying IIV so
    # the ODE remains well-posed under simulation.
    ka <- max(0.01, ka_typ) * exp(etalka)
    cl <- max(0.01, cl_typ) * exp(etalcl)
    vc <- max(0.1,  vc_typ) * exp(etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
