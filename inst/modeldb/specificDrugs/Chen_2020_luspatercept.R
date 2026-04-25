Chen_2020_luspatercept <- function() {
  description <- "One-compartment population PK model for luspatercept (activin receptor type IIB / IgG1 Fc-fusion) in adults with anemia due to myelodysplastic syndromes (Chen 2020), with first-order subcutaneous absorption, first-order linear elimination parameterised in CL/F and V1/F, body weight + age + baseline albumin power covariates on CL/F, and body weight + baseline albumin power covariates on V1/F."
  reference <- "Chen N, Kassir N, Laadem A, Giuseppi AC, Shetty JK, Maxwell SE, Sriraman P, Ritland S, Linde PG, Budda B, Reynolds J, Ramji P, Palmisano M, Zhou S. Population Pharmacokinetics and Exposure-Response of Luspatercept, an Erythroid Maturation Agent, in Anemic Patients With Myelodysplastic Syndromes. CPT Pharmacometrics Syst Pharmacol. 2020 Oct;9(10):395-404. doi:10.1002/psp4.12515"
  vignette <- "Chen_2020_luspatercept"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg per Chen 2020 final-model covariate equations for CL/F and V1/F. Power exponents 0.769 (CL/F) and 0.877 (V1/F) per Table 2.",
      source_name        = "Weight"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 72 years (the dataset median) per Chen 2020 final-model CL/F equation. Power exponent -0.534 per Table 2.",
      source_name        = "Age"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 44 g/L (the dataset median) per Chen 2020 final-model equations for CL/F and V1/F. Power exponents -1.17 (CL/F) and -0.610 (V1/F) per Table 2.",
      source_name        = "Albumin"
    )
  )

  population <- list(
    n_subjects     = 260L,
    n_observations = 2403L,
    n_studies      = 2L,
    age_range      = "27-95 years",
    age_median     = "72 years",
    weight_range   = "46-124 kg",
    weight_median  = "76.3 kg",
    sex_female_pct = 38.8,
    race_ethnicity = c(White = 82.3),
    disease_state  = "Anemia due to lower-risk myelodysplastic syndromes (adults). 72.7% Very-low/Low IPSS-R risk, 23.1% Intermediate, 4.2% High/Very-high. 83.1% with positive ring sideroblasts.",
    dose_range     = "Subcutaneous luspatercept 0.125-1.75 mg/kg every 3 weeks (q3w). 91.0% of patients started at 1.0 mg/kg with stepwise titration to 1.33 or 1.75 mg/kg as needed; the remaining 9.0% received a constant 0.125-0.75 mg/kg.",
    regions        = "Multinational pooled data from a phase II dose-finding/expansion study (A536-03 'PACE-MDS', NCT01749514, n=107) and a pivotal phase III study (ACE-536-MDS-001 'MEDALIST', NCT02631070, n=153).",
    notes          = "Baseline demographics from Chen 2020 Table 1. Renal: 26.9% no impairment, 51.5% mild (eGFR 60-89), 21.5% moderate (eGFR 30-59); none severe. Hepatic: 59.2% none, 31.5% mild, 8.8% moderate, 0.4% severe. Albumin median 44 g/L (31.0-52.6). eGFR median 73.1 mL/min/1.73 m^2 (29.6-150). 38.5% on concurrent iron chelation therapy. 2,403 quantifiable luspatercept concentrations collected 4-784 days after the first dose; ELISA range 50-600 ng/mL; 0.6% of postdose samples were below the limit of quantitation and excluded."
  )

  ini({
    # Structural PK parameters - Chen 2020 Table 2 final-model NONMEM
    # estimates. Reference subject: 70 kg, 72 years, 44 g/L baseline albumin.
    lcl <- log(0.469); label("Apparent clearance CL/F (L/day) at reference covariates")        # Chen 2020 Table 2: CL/F = 0.469 L/day
    lvc <- log(9.22);  label("Apparent central volume V1/F (L) at reference covariates")       # Chen 2020 Table 2: V1/F = 9.22 L
    lka <- log(0.456); label("Absorption rate Ka (1/day)")                                     # Chen 2020 Table 2: Ka = 0.456 1/day

    # Covariate exponents - Chen 2020 Table 2 and the final-model covariate
    # equations in Results: CL/F = 0.469 * (WT/70)^0.769 * (AGE/72)^-0.534 *
    # (ALB/44)^-1.17 and V1/F = 9.22 * (WT/70)^0.877 * (ALB/44)^-0.610.
    e_wt_cl  <-  0.769;  label("Power exponent of (WT/70 kg) on CL/F (unitless)")              # Chen 2020 Table 2: Weight on CL/F = 0.769
    e_age_cl <- -0.534;  label("Power exponent of (AGE/72 yr) on CL/F (unitless)")             # Chen 2020 Table 2: Age on CL/F = -0.534
    e_alb_cl <- -1.17;   label("Power exponent of (ALB/44 g/L) on CL/F (unitless)")            # Chen 2020 Table 2: Albumin on CL/F = -1.17
    e_wt_vc  <-  0.877;  label("Power exponent of (WT/70 kg) on V1/F (unitless)")              # Chen 2020 Table 2: Weight on V1/F = 0.877
    e_alb_vc <- -0.610;  label("Power exponent of (ALB/44 g/L) on V1/F (unitless)")            # Chen 2020 Table 2: Albumin on V1/F = -0.610

    # Inter-individual variability - Chen 2020 Table 2 reports IIV as the
    # square root of omega^2 expressed as a percentage (NONMEM exponential
    # IIV). The cross-check is the AUC IIV: with sqrt(omega^2_CL) = 0.364,
    # the predicted IIV(AUC_ss) is 100*sqrt(exp(0.364^2)-1) = 37.7%, which
    # matches the 38.0% reported in Results for AUC_ss IIV. No IIV on Ka
    # (Methods: 'Inclusion of IIV for absorption rate constant led to large
    # shrinkage'). Diagonal Omega; no correlation reported in Table 2.
    etalcl ~ 0.1325; label("IIV variance on log-CL/F (Chen 2020 Table 2: 36.4% as sqrt(omega^2))")
    etalvc ~ 0.0506; label("IIV variance on log-V1/F (Chen 2020 Table 2: 22.5% as sqrt(omega^2))")

    # Residual error - Chen 2020 Table 2 reports 22.4% for "Residual
    # variability". Methods state luspatercept concentrations were natural-
    # log-transformed prior to analysis with an additive error model, which
    # is equivalent to a proportional error in linear space; the reported
    # 22.4% is sigma on the log scale (= proportional SD in nlmixr2).
    propSd <- 0.224; label("Proportional residual error (fraction)")                            # Chen 2020 Table 2: residual = 22.4% (log-additive)
  })

  model({
    # Individual PK parameters with Chen 2020 final-model covariate equations
    # (Results: 'The final covariate model for CL/F at the population level
    # is described as follows: CL/F (L/day) = 0.469 * (Weight/70)^0.769 *
    # (Age/72)^-0.534 * (Albumin/44)^-1.17' and 'V1/F (L) = 9.22 *
    # (Weight/70)^0.877 * (Albumin/44)^-0.610'). Reference subject: 70 kg,
    # 72 years, 44 g/L baseline albumin.
    cl <- exp(lcl + etalcl) *
          (WT  / 70)^e_wt_cl *
          (AGE / 72)^e_age_cl *
          (ALB / 44)^e_alb_cl
    vc <- exp(lvc + etalvc) *
          (WT  / 70)^e_wt_vc *
          (ALB / 44)^e_alb_vc
    ka <- exp(lka)

    # One-compartment SC model with first-order absorption and first-order
    # elimination (Results: 'A one-compartment model with first-order
    # absorption and elimination best described the concentration-time
    # profiles of luspatercept after subcutaneous injection.'). Apparent
    # CL/F and V1/F absorb the bioavailability into the typical values; the
    # full SC dose enters the depot and no separate F is estimated.
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration in ug/mL (= mg/L); ELISA assay range 50-600 ng/mL =
    # 0.05-0.6 ug/mL. Steady-state mean AUC_ss values predicted for the
    # 1.0-1.75 mg/kg dose range are 151-264 day*ug/mL per Results.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
