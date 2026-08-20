Zhang_2023_valproic_acid_onebindingsite <- function() {
  description <- "One-compartment population PK model with first-order absorption for valproic acid in Chinese children with epilepsy (Zhang 2023 Model I, the one-binding-site protein-binding model). The disposition is carried on the UNBOUND concentration Cu = central/(V/F) with linear unbound clearance; the observed TOTAL plasma concentration is reconstructed as Cc = Cu + Cb, where the albumin-bound concentration follows the one-binding-site isotherm Cb = N * K * Cu * ALB / (1 + K * Cu) with N and K FIXED to the adult literature values of Dutta 2007. Formulation-specific absorption rate constants FIXED from the literature (syrup 2.64 1/h reference, conventional tablet 1.57 1/h, sustained-release tablet 0.46 1/h). One of five protein-binding non-linearity strategies the authors compare; see modellib('Zhang_2023_valproic_acid_base') for the linear reference model."
  reference <- "Zhang L, Liu M, Qin W, Shi D, Mao J, Li Z. Modeling the protein binding non-linearity in population pharmacokinetic model of valproic acid in children with epilepsy: a systematic evaluation study. Front Pharmacol. 2023;14:1228641. doi:10.3389/fphar.2023.1228641. PMID 37860114. Model I binding isotherm from Eq. 3 with N and K quoted from Dutta S, Reed RC. Distinct absorption characteristics of oral formulations of valproic acid/divalproex available in the United States. Epilepsy Res. 2007;73(3):275-83; parameter estimates from Supplementary Table S3."
  vignette <- "Zhang_2023_valproic_acid_protein_binding"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = TRUE,
                   notes = "Holds administered valproic acid amount. Dividing by V/F yields the UNBOUND plasma concentration Cu, not the total: V/F is the unbound-referenced apparent volume (13000 L, about 470-fold the base model's 27.8 L, matching the 0.27% unbound fraction the binding isotherm produces). The observed total concentration is the algebraic sum Cu + Cb.")
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Scales the saturable binding capacity in the one-binding-site isotherm (Zhang 2023 Eq. 3): Cb = N * K * Cu * ALB / (1 + K * Cu). Cohort mean 42.21 g/L, median 42.10, range 29.70-70.50 (Zhang 2023 Table 2); patients with abnormal albumin were excluded from the cohort. Albumin is the only covariate in any of the six models the authors fit. IMPORTANT UNITS NOTE: the paper labels K as 15.5 per mM, but the published parameter estimates only reconcile when Eq. 3 is evaluated with Cu in mg/L, ALB in g/L and K taken at its printed numeric value; see the model file comment on lkassoc_pb and the vignette Errata.",
      source_name        = "ALB"
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
      notes       = "Collected in the evaluation cohort (median 19.00 kg, range 4.00-70.00; Zhang 2023 Table 2) but NOT retained on CL/F or V/F in any of the authors' own six models. Zhang 2023 Supplementary Table S3 reports no weight term."
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
    notes          = "255 total plasma valproic acid trough concentrations in 202 children (139 male / 63 female), measured by gas chromatography (LOQ 1 mg/L, calibration range 12.5-150 mg/L, CV < 10%). Observed concentrations 22.60-118.50 mg/L (median 50.40). All samples were troughs collected under steady-state conditions, so absorption and distribution parameters are only weakly identified. The authors caution that the one-binding-site constants come from an ADULT population and that age is positively correlated with the valproate free fraction, so this strategy may be poorly suited to children - it was nonetheless the third-best of the five on prediction-based metrics (MDPE 3.48%, MAPE 19.38%). Baseline demographics: Zhang 2023 Table 2."
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
    # I: One-binding site model" column; estimate (% RSE)). Both are
    # UNBOUND-referenced apparent (oral) parameters: elimination acts on
    # Cu, and Cu = central / (V/F). Their ratio to the base model's
    # total-concentration parameters (0.311 L/h and 27.8 L) recovers the
    # unbound fraction - 0.236% on CL and 0.214% on V - which agrees
    # with the 0.27% the isotherm itself gives at the steady-state
    # unbound concentration for the cohort mean daily dose. That triple
    # agreement is what identifies the unbound reading, because the
    # paper prints Eq. 3 without saying how it couples to the PK.
    # ----------------------------------------------------------------
    lcl <- log(132.0);   label("Apparent unbound plasma clearance CLp/F (L/h)")        # Zhang 2023 Supplementary Table S3 Model I (CLp/F = 132.0 L/h, 7.3% RSE)
    lvc <- log(13000);   label("Apparent unbound volume of distribution V/F (L)")      # Zhang 2023 Supplementary Table S3 Model I (V/F = 13000 L, 13.8% RSE)

    # ----------------------------------------------------------------
    # One-binding-site isotherm constants (Zhang 2023 Eq. 3), FIXED to
    # the values the authors quote from Dutta 2007 and never estimated.
    #   N = 1.98 binding sites per unit of single-site binding material
    #   K = 15.5 binding-site association constant
    # The paper labels K "mM-1", but evaluating Eq. 3 in molar units
    # gives a total concentration of 3.3 mg/L at the cohort mean daily
    # dose - 17-fold below the observed median of 50.4 mg/L - whereas
    # evaluating it with Cu in mg/L and ALB in g/L gives 59.9 mg/L and
    # reproduces the CL and V ratios above. The as-implemented reading
    # is therefore mg/L + g/L and the printed "mM-1" label is
    # dimensionally inconsistent with the authors' own estimates. See
    # vignette Errata.
    # ----------------------------------------------------------------
    lnsite_pb  <- fixed(log(1.98)); label("Binding sites per unit of single-site binding material (unitless)")  # Zhang 2023 Eq. 3 text (N = 1.98, fixed from Dutta 2007)
    lkassoc_pb <- fixed(log(15.5)); label("Binding-site association constant K (L/mg as implemented)")           # Zhang 2023 Eq. 3 text (K = 15.5, fixed from Dutta 2007; printed as mM-1)

    # ----------------------------------------------------------------
    # IIV - diagonal (no correlation reported). Zhang 2023 reports
    # between-subject variability as CV%, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #   BSV_CL/F = 47.2% -> omega^2 = log(0.472^2 + 1) = 0.20113
    #   BSV_V/F  = 68.9% -> omega^2 = log(0.689^2 + 1) = 0.38847
    # ----------------------------------------------------------------
    etalcl ~ 0.20113  # Zhang 2023 Supplementary Table S3 Model I (BSV_CL/F = 47.2%, 13.0% RSE)
    etalvc ~ 0.38847  # Zhang 2023 Supplementary Table S3 Model I (BSV_V/F = 68.9%, 6.5% RSE)

    # ----------------------------------------------------------------
    # Residual error - combined proportional + additive, on the TOTAL
    # plasma concentration (the measured quantity).
    # ----------------------------------------------------------------
    propSd <- 0.261; label("Proportional residual SD (fraction; 26.1%)")  # Zhang 2023 Supplementary Table S3 Model I (residual error 26.1%, 5.1% RSE)
    addSd  <- 2.6;   label("Additive residual SD (mg/L)")                 # Zhang 2023 Supplementary Table S3 Model I (residual error 2.6 mg/L, 107.3% RSE)
  })

  model({
    # 1. Formulation-specific absorption rate constant. Syrup is the
    #    reference (both indicators 0).
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Individual PK parameters (unbound-referenced; no covariates)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 3. Fixed binding-isotherm constants
    nsite_pb  <- exp(lnsite_pb)
    kassoc_pb <- exp(lkassoc_pb)

    # 4. Micro-constant. Elimination acts on the unbound concentration,
    #    CL * Cu = (CL / V) * central, so the amount kinetics stay
    #    first-order; the non-linearity sits entirely in the observation
    #    equation below.
    kel <- cl / vc

    # 5. One-compartment ODE system with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation. Cu is the unbound plasma concentration; Cb is the
    #    albumin-bound concentration from the one-binding-site isotherm
    #    (Zhang 2023 Eq. 3); the measured quantity is their sum.
    Cu <- central / vc
    Cb <- nsite_pb * kassoc_pb * Cu * ALB / (1 + kassoc_pb * Cu)
    Cc <- Cu + Cb
    Cc ~ add(addSd) + prop(propSd)
  })
}
