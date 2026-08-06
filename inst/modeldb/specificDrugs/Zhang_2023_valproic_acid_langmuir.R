Zhang_2023_valproic_acid_langmuir <- function() {
  description <- "One-compartment population PK model with first-order absorption for valproic acid in Chinese children with epilepsy (Zhang 2023 Model II, the Langmuir protein-binding model). The disposition is carried on the UNBOUND concentration Cu = central/(V/F) with linear unbound clearance; the observed TOTAL plasma concentration is reconstructed as Cc = Cu + Cb, where the bound concentration follows the single-site Langmuir isotherm Cb = Bm * Cu / (Kd + Cu) with Bm and Kd FIXED to the literature values of Ueshima 2008. Formulation-specific absorption rate constants FIXED from the literature (syrup 2.64 1/h reference, conventional tablet 1.57 1/h, sustained-release tablet 0.46 1/h). One of five protein-binding non-linearity strategies the authors compare; see modellib('Zhang_2023_valproic_acid_nonsaturable') for the linear non-saturable extension the authors ultimately preferred."
  reference <- "Zhang L, Liu M, Qin W, Shi D, Mao J, Li Z. Modeling the protein binding non-linearity in population pharmacokinetic model of valproic acid in children with epilepsy: a systematic evaluation study. Front Pharmacol. 2023;14:1228641. doi:10.3389/fphar.2023.1228641. PMID 37860114. Model II binding isotherm from Eq. 4 with Kd and Bm quoted from Ueshima S, et al. (2008); parameter estimates from Supplementary Table S3."
  vignette <- "Zhang_2023_valproic_acid_protein_binding"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = TRUE,
                   notes = "Holds administered valproic acid amount. Dividing by V/F yields the UNBOUND plasma concentration Cu, not the total: V/F is the unbound-referenced apparent volume (341 L, about 12-fold the base model's 27.8 L, matching the 9.8% unbound fraction the Langmuir isotherm produces). The observed total concentration is the algebraic sum Cu + Cb.")
  )

  covariateData <- list(
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
      notes       = "Collected in the evaluation cohort (median 19.00 kg, range 4.00-70.00; Zhang 2023 Table 2) but NOT retained on CL/F or V/F in any of the authors' own six models. Zhang 2023 Supplementary Table S3 reports no weight term."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected in the evaluation cohort (mean 42.21 g/L, median 42.10, range 29.70-70.50; Zhang 2023 Table 2). The Langmuir isotherm of Eq. 4 folds the binding capacity into the single constant Bm and does not scale it by albumin; albumin enters explicitly only in the Model I one-binding-site strategy (Eq. 3)."
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
    notes          = "255 total plasma valproic acid trough concentrations in 202 children (139 male / 63 female), measured by gas chromatography (LOQ 1 mg/L, calibration range 12.5-150 mg/L, CV < 10%). Observed concentrations 22.60-118.50 mg/L (median 50.40). All samples were troughs collected under steady-state conditions, so absorption and distribution parameters are only weakly identified. Prediction-based performance MDPE 0.82%, MAPE 20.49% (Zhang 2023 Table 3). Baseline demographics: Zhang 2023 Table 2."
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
    # II: Langmuir equation" column; estimate (% RSE)). Both are
    # UNBOUND-referenced apparent (oral) parameters: elimination acts on
    # Cu, and Cu = central / (V/F). Their ratio to the base model's
    # total-concentration parameters (0.311 L/h and 27.8 L) recovers the
    # unbound fraction - 9.15% on CL and 8.15% on V - which agrees with
    # the 9.8% the Langmuir isotherm itself gives at the steady-state
    # unbound concentration for the cohort mean daily dose. That triple
    # agreement is what identifies the unbound reading, because the
    # paper prints Eq. 4 without saying how it couples to the PK.
    # ----------------------------------------------------------------
    lcl <- log(3.4);   label("Apparent unbound plasma clearance CLp/F (L/h)")    # Zhang 2023 Supplementary Table S3 Model II (CLp/F = 3.4 L/h, 4.6% RSE)
    lvc <- log(341.0); label("Apparent unbound volume of distribution V/F (L)")  # Zhang 2023 Supplementary Table S3 Model II (V/F = 341.0 L, 5.5% RSE)

    # ----------------------------------------------------------------
    # Langmuir isotherm constants (Zhang 2023 Eq. 4), FIXED to the
    # values the authors quote from Ueshima 2008 and never estimated.
    # ----------------------------------------------------------------
    lbmax_pb <- fixed(log(130.0)); label("Maximum binding-site concentration Bm (mg/L)")  # Zhang 2023 Eq. 4 text (Bm = 130 mg/L, fixed from Ueshima 2008)
    lkd_pb   <- fixed(log(7.8));   label("Dissociation constant Kd for the saturable site (mg/L)")  # Zhang 2023 Eq. 4 text (Kd = 7.8, fixed from Ueshima 2008)

    # ----------------------------------------------------------------
    # IIV - diagonal (no correlation reported). Zhang 2023 reports
    # between-subject variability as CV%, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #   BSV_CL/F = 51.9% -> omega^2 = log(0.519^2 + 1) = 0.23851
    #   BSV_V/F  = 72.2% -> omega^2 = log(0.722^2 + 1) = 0.41955
    # ----------------------------------------------------------------
    etalcl ~ 0.23851  # Zhang 2023 Supplementary Table S3 Model II (BSV_CL/F = 51.9%, 8.1% RSE)
    etalvc ~ 0.41955  # Zhang 2023 Supplementary Table S3 Model II (BSV_V/F = 72.2%, 5.5% RSE)

    # ----------------------------------------------------------------
    # Residual error - combined proportional + additive, on the TOTAL
    # plasma concentration (the measured quantity).
    # ----------------------------------------------------------------
    propSd <- 0.125; label("Proportional residual SD (fraction; 12.5%)")  # Zhang 2023 Supplementary Table S3 Model II (residual error 12.5%, 11.2% RSE)
    addSd  <- 4.3;   label("Additive residual SD (mg/L)")                 # Zhang 2023 Supplementary Table S3 Model II (residual error 4.3 mg/L, 37.3% RSE)
  })

  model({
    # 1. Formulation-specific absorption rate constant. Syrup is the
    #    reference (both indicators 0).
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Individual PK parameters (unbound-referenced; no covariates)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 3. Fixed binding-isotherm constants
    bmax_pb <- exp(lbmax_pb)
    kd_pb   <- exp(lkd_pb)

    # 4. Micro-constant. Elimination acts on the unbound concentration,
    #    CL * Cu = (CL / V) * central, so the amount kinetics stay
    #    first-order; the non-linearity sits entirely in the observation
    #    equation below.
    kel <- cl / vc

    # 5. One-compartment ODE system with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation. Cu is the unbound plasma concentration; Cb is the
    #    protein-bound concentration from the Langmuir isotherm
    #    (Zhang 2023 Eq. 4); the measured quantity is their sum.
    Cu <- central / vc
    Cb <- bmax_pb * Cu / (kd_pb + Cu)
    Cc <- Cu + Cb
    Cc ~ add(addSd) + prop(propSd)
  })
}
