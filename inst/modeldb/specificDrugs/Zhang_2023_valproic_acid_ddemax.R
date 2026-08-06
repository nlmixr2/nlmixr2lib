Zhang_2023_valproic_acid_ddemax <- function() {
  description <- "One-compartment population PK model with first-order absorption for total plasma valproic acid in Chinese children with epilepsy (Zhang 2023 Model III, the dose-dependent maximum effect model). Apparent clearance rises with the patient's own total daily dose through a sigmoid-Emax term, CL/F = CLp/F * (1 + Emax * DD^gamma / (DD50^gamma + DD^gamma)), with Emax, gamma and DD50 all FIXED to the values of Ding 2015. Formulation-specific absorption rate constants FIXED from the literature (syrup 2.64 1/h reference, conventional tablet 1.57 1/h, sustained-release tablet 0.46 1/h). The weakest of the five non-linearity strategies the authors compare - its objective function value is worse than the linear base model's - and it is packaged for fidelity to the published comparison; see modellib('Zhang_2023_valproic_acid_nonsaturable') for the strategy the authors recommend. DD50 is not reported by Zhang 2023 and is carried over from Ding 2015; see the vignette Errata."
  reference <- "Zhang L, Liu M, Qin W, Shi D, Mao J, Li Z. Modeling the protein binding non-linearity in population pharmacokinetic model of valproic acid in children with epilepsy: a systematic evaluation study. Front Pharmacol. 2023;14:1228641. doi:10.3389/fphar.2023.1228641. PMID 37860114. Model III structure from Eq. 5; Emax, gamma and DD50 quoted from Ding J, Wang Y, Lin W, et al. A population pharmacokinetic model of valproic acid in pediatric patients with epilepsy: a non-linear pharmacokinetic model based on protein-binding saturation. Clin Pharmacokinet. 2015;54(3):305-17. doi:10.1007/s40262-014-0212-8 (tabulated in Zhang 2023 Table 1); CLp/F, V/F, BSV and residual error from Supplementary Table S3."
  vignette <- "Zhang_2023_valproic_acid_protein_binding"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = TRUE,
                   notes = "Holds total valproic acid amount; central/(V/F) is the TOTAL plasma concentration. Unlike the paper's Models I, II and IV, this strategy carries no explicit binding isotherm - the protein-binding non-linearity is absorbed empirically into the dose dependence of CL/F.")
  )

  covariateData <- list(
    DOSE_VPA_MGD = list(
      description        = "Patient's own total daily valproic acid dose",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the sigmoid-Emax term on apparent clearance (Zhang 2023 Eq. 5): cl = CLp/F * (1 + Emax * DOSE_VPA_MGD^gamma / (DD50^gamma + DOSE_VPA_MGD^gamma)). Per-dose-record covariate, constant within an inter-dose interval and updated when the prescriber alters the daily dose. Observed range 60-1250 mg/day (median 480, mean 513.04; Zhang 2023 Table 2). Must be strictly positive. Eq. 5's text defines DD in mg/kg/day, but DD50 is only reconcilable with the paper's reported clearance and predictive metrics on the mg/day scale - operator-ratified reading, see the model file comment on dd_50 and the vignette Errata.",
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
      notes       = "Collected in the evaluation cohort (median 19.00 kg, range 4.00-70.00; Zhang 2023 Table 2) but NOT retained on CL/F or V/F in any of the authors' own six models. Zhang 2023 Supplementary Table S3 reports no weight term, and on the ratified mg/day reading of Eq. 5 weight does not enter this model even indirectly."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected in the evaluation cohort (mean 42.21 g/L, median 42.10, range 29.70-70.50; Zhang 2023 Table 2). Albumin enters explicitly only in the Model I one-binding-site strategy (Eq. 3); this model represents binding saturation empirically through the dose dependence of CL/F."
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
    notes          = "255 total plasma valproic acid trough concentrations in 202 children (139 male / 63 female), measured by gas chromatography (LOQ 1 mg/L, calibration range 12.5-150 mg/L, CV < 10%). Observed concentrations 22.60-118.50 mg/L (median 50.40). All samples were troughs collected under steady-state conditions, so absorption and distribution parameters are only weakly identified. This model was the poorest performer of the five strategies: its objective function value (1773.0) is 20.9 points WORSE than the linear base model's (1752.1), and it had the highest MAPE (27.67%) and lowest F30 (52.94%) of the five (Zhang 2023 Table 3, Supplementary Table S3). The authors nonetheless note it is the only empirical strategy whose tendency plot reproduces a genuinely non-linear dose-clearance relationship. Baseline demographics: Zhang 2023 Table 2."
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
    # III: Dose-dependent maximum effect model" column; estimate
    # (% RSE)). CLp/F is the zero-dose intercept of Eq. 5, NOT the
    # typical clearance: the sigmoid-Emax multiplier is close to its
    # 1 + Emax = 3.8 ceiling across the observed dose range, so typical
    # CL/F is about 0.307 L/h - essentially the base model's 0.311.
    # ----------------------------------------------------------------
    lcl <- log(0.0815); label("Apparent plasma clearance CLp/F at zero daily dose (L/h)")  # Zhang 2023 Supplementary Table S3 Model III (CLp/F = 0.0815 L/h, 33.7% RSE)
    lvc <- log(52.3);   label("Apparent volume of distribution V/F (L)")                   # Zhang 2023 Supplementary Table S3 Model III (V/F = 52.3 L, 12.6% RSE)

    # ----------------------------------------------------------------
    # Sigmoid-Emax dose effect on CL/F (Zhang 2023 Eq. 5). All three
    # constants are FIXED, none is estimated on the Zhang cohort.
    #
    # Emax = 2.8 and gamma = 1.68 are stated in the Eq. 5 text, quoted
    # from Ding 2015.
    #
    # DD50 is NOT reported anywhere in Zhang 2023 - it is absent from
    # Supplementary Table S3 (whose CL_DD row is "/" for this model) and
    # from the body text. It is carried over from the same Ding 2015
    # source as Emax and gamma, where Zhang's own Table 1 tabulates
    # Ding's CL/F term as (1 + 2.8 * DDW^1.68 / (37.4^1.68 + DDW^1.68)).
    # OPERATOR-RATIFIED READING (sidecar oare_PMC10587682 q2, option C):
    # DD50 = 37.4 evaluated on the mg/DAY scale, even though Ding
    # tabulates it against DDW in mg/kg/day. Rationale: only the mg/day
    # reading reproduces Zhang's own numbers - typical CL/F comes out at
    # 0.307 L/h against the base model's 0.311 L/h and a typical
    # steady-state concentration of 69.7 mg/L, in family with the
    # 60-70 mg/L that all five other models give. On the mg/kg/day
    # reading typical CL/F is 0.153 L/h and the typical concentration
    # about 137 mg/L, incompatible with the -1.10% MDPE the paper
    # reports (Table 3). Numerical corroboration over printed units,
    # per the standing text-versus-equation policy. Supplementary
    # Figure S5A is the one piece of evidence that favours the
    # mg/kg/day reading; the conflict is documented in the vignette
    # Errata rather than resolved.
    # ----------------------------------------------------------------
    dd_emax  <- fixed(2.8);  label("Maximum fractional increase in CL/F at saturating daily dose (unitless)")  # Zhang 2023 Eq. 5 text (Emax = 2.8, fixed from Ding 2015)
    dd_50    <- fixed(37.4); label("Daily dose giving half-maximal increase in CL/F, DD50 (mg/day)")           # Ding 2015 via Zhang 2023 Table 1 (DD50 = 37.4); NOT reported in Zhang 2023 - see comment above
    dd_gamma <- fixed(1.68); label("Hill coefficient of the daily-dose effect on CL/F (unitless)")             # Zhang 2023 Eq. 5 text (gamma = 1.68, fixed from Ding 2015)

    # ----------------------------------------------------------------
    # IIV - diagonal (no correlation reported). Zhang 2023 reports
    # between-subject variability as CV%, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #   BSV_CL/F = 68.6% -> omega^2 = log(0.686^2 + 1) = 0.38567
    #   BSV_V/F  = 45.9% -> omega^2 = log(0.459^2 + 1) = 0.19118
    # ----------------------------------------------------------------
    etalcl ~ 0.38567  # Zhang 2023 Supplementary Table S3 Model III (BSV_CL/F = 68.6%, 20.5% RSE)
    etalvc ~ 0.19118  # Zhang 2023 Supplementary Table S3 Model III (BSV_V/F = 45.9%, 12.4% RSE)

    # ----------------------------------------------------------------
    # Residual error - combined proportional + additive
    # ----------------------------------------------------------------
    propSd <- 0.117; label("Proportional residual SD (fraction; 11.7%)")  # Zhang 2023 Supplementary Table S3 Model III (residual error 11.7%, 17.9% RSE)
    addSd  <- 6.8;   label("Additive residual SD (mg/L)")                 # Zhang 2023 Supplementary Table S3 Model III (residual error 6.8 mg/L, 31.9% RSE)
  })

  model({
    # 1. Formulation-specific absorption rate constant. Syrup is the
    #    reference (both indicators 0).
    ka <- exp(lka + e_form_tablet_ka * FORM_TABLET + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Individual PK parameters. CL/F carries the sigmoid-Emax daily
    #    dose effect of Zhang 2023 Eq. 5.
    ddeff <- 1 + dd_emax * DOSE_VPA_MGD^dd_gamma / (dd_50^dd_gamma + DOSE_VPA_MGD^dd_gamma)
    cl <- exp(lcl + etalcl) * ddeff
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
