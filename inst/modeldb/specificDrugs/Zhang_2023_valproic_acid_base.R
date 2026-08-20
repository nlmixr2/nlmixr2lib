Zhang_2023_valproic_acid_base <- function() {
  description <- "One-compartment population PK model with first-order absorption for total plasma valproic acid in Chinese children with epilepsy (Zhang 2023 base model). Linear clearance with no covariates on CL/F or V/F; formulation-specific absorption rate constants FIXED from the literature (syrup 2.64 1/h reference, conventional tablet 1.57 1/h, sustained-release tablet 0.46 1/h). This is the reference model against which the paper's five protein-binding non-linearity strategies are compared; see modellib('Zhang_2023_valproic_acid_exponent') for the daily-dose power model."
  reference <- "Zhang L, Liu M, Qin W, Shi D, Mao J, Li Z. Modeling the protein binding non-linearity in population pharmacokinetic model of valproic acid in children with epilepsy: a systematic evaluation study. Front Pharmacol. 2023;14:1228641. doi:10.3389/fphar.2023.1228641. PMID 37860114. Base-model parameter estimates from Supplementary Table S3."
  vignette <- "Zhang_2023_valproic_acid_protein_binding"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = FALSE)
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
      notes       = "Collected in the evaluation cohort (median 19.00 kg, range 4.00-70.00; Zhang 2023 Table 2) but NOT retained on CL/F or V/F in the base model or in any of the five protein-binding models. Zhang 2023 Supplementary Table S3 reports no weight term. The Discussion notes body weight was the single most common covariate across the ten reviewed literature models and questions whether the 3/4 allometric exponent is reliable in children, but the authors' own models are not weight-scaled."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected in the evaluation cohort (mean 42.21 g/L, median 42.10, range 29.70-70.50; Zhang 2023 Table 2) and screened out of the base model. Albumin enters only the paper's Model I one-binding-site strategy (Eq. 3), not this base model."
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
    notes          = "255 total plasma valproic acid trough concentrations in 202 children (139 male / 63 female), measured by gas chromatography (LOQ 1 mg/L, calibration range 12.5-150 mg/L, CV < 10%). Observed concentrations 22.60-118.50 mg/L (median 50.40). All samples were troughs collected under steady-state conditions, so absorption and distribution parameters are only weakly identified - the authors state this explicitly as a limitation, which is why V/F is large relative to the paediatric literature. Patients with hepatic or renal impairment, abnormal albumin, or poor adherence were excluded. 77 of 202 received a concomitant antiseizure medication, most commonly levetiracetam (27), oxcarbazepine (24) and topiramate (19); no co-medication effect was retained. Baseline demographics: Zhang 2023 Table 2."
  )

  ini({
    # ----------------------------------------------------------------
    # Absorption - FIXED to the literature Ka set the authors adopted.
    # Zhang 2023 Supplementary Table S3 footnote: "Ka was fixed to
    # 2.64, 1.57, 0.46 for syrup, conventional tablet and SR tablet,
    # respectively". Oral syrup is the reference formulation, so lka is
    # the syrup value and the two indicators are log-ratio shifts
    # (Aruldhas 2021 methadone / Kleideiter 2017 cebranopadol pattern
    # for a multi-level formulation-specific Ka).
    # ----------------------------------------------------------------
    lka <- fixed(log(2.64)); label("Absorption rate constant, oral syrup reference (1/h)")            # Zhang 2023 Supplementary Table S3 footnote (syrup Ka = 2.64 1/h)
    e_form_tablet_ka <- fixed(log(1.57 / 2.64)); label("Log-ratio shift on Ka for conventional tablet vs syrup") # Zhang 2023 Supplementary Table S3 footnote (conventional tablet Ka = 1.57 1/h)
    e_form_vpa_sr_ka <- fixed(log(0.46 / 2.64)); label("Log-ratio shift on Ka for sustained-release tablet vs syrup") # Zhang 2023 Supplementary Table S3 footnote (SR tablet Ka = 0.46 1/h)

    # ----------------------------------------------------------------
    # Structural parameters (Zhang 2023 Supplementary Table S3, "Base
    # model" column; estimate (% relative standard error)). Both are
    # apparent (oral) parameters on the TOTAL plasma concentration.
    # No covariates were retained on CL/F or V/F.
    # ----------------------------------------------------------------
    lcl <- log(0.311); label("Apparent plasma clearance CLp/F (L/h)")            # Zhang 2023 Supplementary Table S3 base model (CLp/F = 0.311 L/h, 4.6% RSE)
    lvc <- log(27.8);  label("Apparent volume of distribution V/F (L)")          # Zhang 2023 Supplementary Table S3 base model (V/F = 27.8 L, 8.8% RSE)

    # ----------------------------------------------------------------
    # IIV - diagonal (no correlation reported). Zhang 2023 reports
    # between-subject variability as CV%, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #   BSV_CL/F = 45.7% -> omega^2 = log(0.457^2 + 1) = 0.18966
    #   BSV_V/F  = 68.6% -> omega^2 = log(0.686^2 + 1) = 0.38567
    # ----------------------------------------------------------------
    etalcl ~ 0.18966  # Zhang 2023 Supplementary Table S3 base model (BSV_CL/F = 45.7%, 8.1% RSE)
    etalvc ~ 0.38567  # Zhang 2023 Supplementary Table S3 base model (BSV_V/F = 68.6%, 8.9% RSE)

    # ----------------------------------------------------------------
    # Residual error - combined proportional + additive
    # ----------------------------------------------------------------
    propSd <- 0.089; label("Proportional residual SD (fraction; 8.9%)")          # Zhang 2023 Supplementary Table S3 base model (residual error 8.9%, 29.1% RSE)
    addSd  <- 6.3;   label("Additive residual SD (mg/L)")                        # Zhang 2023 Supplementary Table S3 base model (residual error 6.3 mg/L, 35.6% RSE)
  })

  model({
    # 1. Formulation-specific absorption rate constant. Syrup is the
    #    reference (both indicators 0); exactly one of the three
    #    formulation levels applies to any dose record.
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Individual PK parameters (no covariate effects in the base model)
    cl <- exp(lcl + etalcl)
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
