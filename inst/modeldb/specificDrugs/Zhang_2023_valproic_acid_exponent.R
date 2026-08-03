Zhang_2023_valproic_acid_exponent <- function() {
  description <- "One-compartment population PK model with first-order absorption for total plasma valproic acid in Chinese children with epilepsy (Zhang 2023 Model V, the simple exponent model). Apparent clearance follows a power function of the patient's own daily dose per kilogram, CL/F = CLp/F * (DD/25)^0.658, with DD in mg/kg/day and a reference daily dose of 25 mg/kg/day; this is the empirical dose-dependence strategy the authors contrast against mechanistic protein-binding models. Formulation-specific absorption rate constants FIXED from the literature (syrup 2.64 1/h reference, conventional tablet 1.57 1/h, sustained-release tablet 0.46 1/h). Best prediction-based performance of the five strategies (MDPE 1.50%, MAPE 17.68%) but the authors conclude it does not describe the underlying non-linearity; see modellib('Zhang_2023_valproic_acid_base') for the reference model."
  reference <- "Zhang L, Liu M, Qin W, Shi D, Mao J, Li Z. Modeling the protein binding non-linearity in population pharmacokinetic model of valproic acid in children with epilepsy: a systematic evaluation study. Front Pharmacol. 2023;14:1228641. doi:10.3389/fphar.2023.1228641. PMID 37860114. Model V structure from Eq. 7; parameter estimates from Supplementary Table S3."
  vignette <- "Zhang_2023_valproic_acid_protein_binding"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    DOSE_VPA_MGKGD = list(
      description        = "Patient's own total daily valproic acid dose per kilogram body weight",
      units              = "mg/kg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters as a self-dose-rate power regressor on CL/F: cl = CLp/F * (DOSE_VPA_MGKGD / 25)^0.658 (Zhang 2023 Eq. 7). The reference value 25 mg/kg/day is the rounded cohort mean daily dose (24.50 mg/kg/day, median 23.44; Zhang 2023 Table 2). Per-dose-record covariate, constant within an inter-dose interval and updated when the prescriber alters the daily dose. Observed range 8.70-57.69 mg/kg/day. Must be strictly positive - a zero or NA value makes the power term undefined.",
      source_name        = "DD"
    ),
    FORM_TABLET = list(
      description        = "Conventional (immediate-release) oral tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral syrup, the reference formulation in this cohort)",
      notes              = "Selects the FIXED conventional-tablet absorption rate constant Ka = 1.57 1/h. Oral syrup (Ka = 2.64 1/h) is the reference formulation: FORM_TABLET = 0 and FORM_VPA_SR = 0. The Zhang 2023 evaluation cohort itself received only syrup (194 records) or sustained-release tablet (61 records); the conventional-tablet level is retained so the full literature Ka set the authors quote is reachable. Zhang 2023 Supplementary Table S3 footnote.",
      source_name        = "conventional tablet"
    ),
    FORM_VPA_SR = list(
      description        = "Sustained-release valproic acid tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral syrup, the reference formulation in this cohort)",
      notes              = "Selects the FIXED sustained-release-tablet absorption rate constant Ka = 0.46 1/h. 61 of 255 observation records in the Zhang 2023 cohort were on the sustained-release tablet (Table 2). Zhang 2023 Supplementary Table S3 footnote.",
      source_name        = "sustained release tablet"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Collected in the evaluation cohort (median 19.00 kg, range 4.00-70.00; Zhang 2023 Table 2) but NOT retained on CL/F or V/F in any of the authors' own models. Weight enters this model only indirectly, through the per-kilogram normalisation of the daily dose covariate DOSE_VPA_MGKGD."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected in the evaluation cohort (mean 42.21 g/L, median 42.10, range 29.70-70.50; Zhang 2023 Table 2). Albumin enters only the paper's Model I one-binding-site strategy (Eq. 3), not this simple exponent model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 202,
    n_studies      = 1,
    age_range      = "0.17-15.00 years",
    age_median     = "4.92 years",
    weight_range   = "4.00-70.00 kg",
    weight_median  = "19.00 kg",
    sex_female_pct = 31.2,
    race_ethnicity = "Chinese (single-centre Han-predominant cohort; sub-ethnicity not reported)",
    disease_state  = "Childhood epilepsy on maintenance valproic acid with routine therapeutic drug monitoring",
    dose_range     = "60-1250 mg/day (median 480); 8.70-57.69 mg/kg/day (median 23.44); oral syrup or sustained-release tablet given once, twice or three times daily",
    regions        = "China (Wuhan Children's Hospital, single centre, January 2016 - November 2018)",
    notes          = "255 total plasma valproic acid trough concentrations in 202 children (139 male / 63 female), measured by gas chromatography (LOQ 1 mg/L, calibration range 12.5-150 mg/L, CV < 10%). Observed concentrations 22.60-118.50 mg/L (median 50.40). All samples were troughs collected under steady-state conditions, so absorption and distribution parameters are only weakly identified. The authors caution in the Discussion that because the daily dose is itself the quantity a therapeutic-drug-monitoring model is meant to predict, using it as a CL/F covariate is circular - this model is included for fidelity to the published comparison, not as a recommended dosing tool. Baseline demographics: Zhang 2023 Table 2."
  )

  ini({
    # ----------------------------------------------------------------
    # Absorption - FIXED to the literature Ka set the authors adopted.
    # Zhang 2023 Supplementary Table S3 footnote: "Ka was fixed to
    # 2.64, 1.57, 0.46 for syrup, conventional tablet and SR tablet,
    # respectively". Oral syrup is the reference formulation.
    # ----------------------------------------------------------------
    lka <- fixed(log(2.64)); label("Absorption rate constant, oral syrup reference (1/h)")            # Zhang 2023 Supplementary Table S3 footnote (syrup Ka = 2.64 1/h)
    e_form_tablet_ka <- fixed(log(1.57 / 2.64)); label("Log-ratio shift on Ka for conventional tablet vs syrup") # Zhang 2023 Supplementary Table S3 footnote (conventional tablet Ka = 1.57 1/h)
    e_form_vpa_sr_ka <- fixed(log(0.46 / 2.64)); label("Log-ratio shift on Ka for sustained-release tablet vs syrup") # Zhang 2023 Supplementary Table S3 footnote (SR tablet Ka = 0.46 1/h)

    # ----------------------------------------------------------------
    # Structural parameters (Zhang 2023 Supplementary Table S3, "Model
    # V: The simple exponent model" column; estimate (% RSE)). Both are
    # apparent (oral) parameters on the TOTAL plasma concentration.
    # ----------------------------------------------------------------
    lcl <- log(0.331); label("Apparent plasma clearance CLp/F at 25 mg/kg/day (L/h)")  # Zhang 2023 Supplementary Table S3 Model V (CLp/F = 0.331 L/h, 2.1% RSE)
    lvc <- log(17.2);  label("Apparent volume of distribution V/F (L)")                # Zhang 2023 Supplementary Table S3 Model V (V/F = 17.2 L, 20.0% RSE)

    # ----------------------------------------------------------------
    # Daily-dose power exponent on CL/F (Zhang 2023 Eq. 7:
    # CL/F = CLp * (DD/25)^k, DD in mg/kg/day). Reported in
    # Supplementary Table S3 as the row "CL_DD". Estimated on the
    # authors' own dataset (the only parameter in Eq. 7 not FIXED).
    # ----------------------------------------------------------------
    e_dose_vpa_mgkgd_cl <- 0.658; label("Power exponent of daily dose (mg/kg/day) on CL/F (unitless)")  # Zhang 2023 Supplementary Table S3 Model V (CL_DD = 0.658, 7.7% RSE); form from Eq. 7

    # ----------------------------------------------------------------
    # IIV - diagonal (no correlation reported). Zhang 2023 reports
    # between-subject variability as CV%, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #   BSV_CL/F = 24.3% -> omega^2 = log(0.243^2 + 1) = 0.05737
    #   BSV_V/F  = 45.1% -> omega^2 = log(0.451^2 + 1) = 0.18515
    # ----------------------------------------------------------------
    etalcl ~ 0.05737  # Zhang 2023 Supplementary Table S3 Model V (BSV_CL/F = 24.3%, 11.5% RSE)
    etalvc ~ 0.18515  # Zhang 2023 Supplementary Table S3 Model V (BSV_V/F = 45.1%, 15.9% RSE)

    # ----------------------------------------------------------------
    # Residual error - combined proportional + additive
    # ----------------------------------------------------------------
    propSd <- 0.141; label("Proportional residual SD (fraction; 14.1%)")  # Zhang 2023 Supplementary Table S3 Model V (residual error 14.1%, 11.9% RSE)
    addSd  <- 4.2;   label("Additive residual SD (mg/L)")                 # Zhang 2023 Supplementary Table S3 Model V (residual error 4.2 mg/L, 35.7% RSE)
  })

  model({
    # 1. Formulation-specific absorption rate constant. Syrup is the
    #    reference (both indicators 0).
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Individual PK parameters. CL/F is a power function of the
    #    patient's own daily dose per kg, referenced to 25 mg/kg/day
    #    (Zhang 2023 Eq. 7).
    cl <- exp(lcl + etalcl) * (DOSE_VPA_MGKGD / 25)^e_dose_vpa_mgkgd_cl
    vc <- exp(lvc + etalvc)

    # 3. Micro-constant
    kel <- cl / vc

    # 4. One-compartment ODE system with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation (total plasma valproic acid) and residual error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
