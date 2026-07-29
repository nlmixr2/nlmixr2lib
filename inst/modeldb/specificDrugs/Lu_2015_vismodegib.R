Lu_2015_vismodegib <- function() {
  description <- "Semi-mechanism-based one-compartment population pharmacokinetic model for vismodegib (GDC-0449, oral Hedgehog pathway inhibitor) in adults with advanced solid tumors and healthy volunteers. First-order absorption, first-order elimination of unbound drug, and saturable fast-equilibrium binding to alpha-1-acid glycoprotein (AAG) jointly describe total and unbound plasma vismodegib concentrations. AAG is supplied as a time-varying covariate (uM); covariates retained on disposition are age (power on CLunbound, reference 60 years) and body weight (power on Vc, reference 75 kg); formulation (Phase I dry-blend capsule vs Phase II wet-granulation commercial capsule) and population (healthy volunteer vs patient) shift Ka and relative bioavailability F (Lu 2015)."
  reference   <- "Lu T, Wang B, Gao Y, Dresser M, Graham RA, Jin JY. Semi-Mechanism-Based Population Pharmacokinetic Modeling of the Hedgehog Pathway Inhibitor Vismodegib. CPT Pharmacometrics Syst Pharmacol. 2015;4(11):680-689. doi:10.1002/psp4.12039"
  vignette    <- "Lu_2015_vismodegib"
  units       <- list(time = "day", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the apparent clearance of unbound drug: CLunbound = exp(lcl + etalcl) * (AGE / 60)^e_age_cl. Reference 60 years is the median patient age in the Lu 2015 cohort (Table 2 footnote / Eq. 5).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the central volume of distribution: Vc = exp(lvc + etalvc) * (WT / 75)^e_wt_vc. Reference 75 kg is the patient-cohort weight reference in Lu 2015 Eq. 5 / Table 2 footnote.",
      source_name        = "WT"
    ),
    AAG = list(
      description        = "Plasma alpha-1-acid glycoprotein (AAG; orosomucoid) concentration. Time-varying covariate that drives the saturable fast-equilibrium binding of vismodegib in plasma.",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Lu 2015 reports AAG in micromolar (uM) throughout (typical baseline 31.1 uM in patients, 20.2 uM in healthy volunteers; Methods + Results paragraph 1 of Study Population). The canonical register lists AAG in g/L; users converting from g/L should multiply by 1000 / AAG_MW where AAG MW is ~41 kDa (so 1 g/L ~ 24.4 uM). The covariate enters as a STRUCTURAL term inside the saturable-binding mass balance (Methods, Structural model + Discussion paragraph on Widmer 2006 parameterisation), NOT as a multiplicative effect on a PK parameter. Missing AAG values were imputed by Lu 2015 using Last Observation Carried Forward (LOCF); when supplying simulated covariate columns, hold AAG constant within a dosing interval unless time-varying data are available.",
        sep = " "),
      source_name        = "AAG"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator (1 = healthy volunteer, 0 = patient with advanced solid tumor / locally-advanced or metastatic basal cell carcinoma).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patient pooled across SHH3925g, SHH4610g, SHH4476g)",
      notes              = "Two effects in Lu 2015 Eq. 4: (a) multiplicative shift on ka (estimable only in combination with formulation); typical ka in the Phase II / patient reference 9.025 day^-1 is multiplied by exp(0.671) = 1.96 for HV. (b) multiplicative shift on F applied only to the Phase I formulation (FORM_VISMO_PHASEI = 1); F = exp(log(0.346) + 0.881 * DIS_HEALTHY) for the Phase I formulation, so Phase I in patients has F = 0.346 and Phase I in HV has F = 0.836 (Lu 2015 Table 1). HV cohort enrolled in SHH4433g + SHH4683g.",
      source_name        = "HV"
    ),
    FORM_VISMO_PHASEI = list(
      description        = "Vismodegib formulation indicator (1 = Phase I development formulation, dry-blend capsules; 0 = Phase II / commercial formulation, wet-granulation capsules).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase II / commercial wet-granulation capsule; F fixed at 1 and ka = exp(lka) day^-1 reference per Lu 2015 Eq. 4)",
      notes              = "Per-subject formulation indicator (each subject received a single formulation across their PK profile). Two effects on absorption (Lu 2015 Eq. 4): (a) multiplicative shift on ka via exp(-0.602) = 0.55x; (b) multiplicative shift on F: F = 1 (Phase II reference) vs F = exp(log(0.346) + 0.881 * DIS_HEALTHY) (Phase I, value depends on healthy-volunteer status). The Phase I formulation was used in SHH3925g and the Phase I cohort of SHH4610g (24.4% of subjects, Lu 2015 Results paragraph 1); the Phase II formulation was used in SHH4476g, SHH4433g, and SHH4683g (75.6% of subjects).",
      source_name        = NULL
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 225L,
    n_studies      = 5L,
    age_range      = "26-89 years (patient cohort median 59.7 years); 47-65 years (HV cohort)",
    age_median     = "60 years (patient reference; Lu 2015 Eq. 5 typical value)",
    weight_range   = "Not reported in detail; reference 75 kg (Lu 2015 Eq. 5)",
    weight_median  = "75 kg (reference covariate value, Lu 2015 Eq. 5)",
    sex_female_pct = 46.2,
    race_ethnicity = "Patients: predominantly Caucasian (97.1%); HV cohort comprised of women of non-childbearing potential (WONCBP), all Caucasian.",
    disease_state  = "Adults with advanced solid tumors including locally-advanced or metastatic basal cell carcinoma (n = 204 patients across SHH3925g, SHH4610g, SHH4476g) pooled with healthy women of non-childbearing potential (n = 21 across SHH4433g, SHH4683g).",
    dose_range     = "Oral vismodegib once daily 150-540 mg, with single-dose and multiple-dose data; 150 mg QD is the approved dose for the Phase II / commercial formulation.",
    regions        = "United States (FDA-approved indication is metastatic or locally-advanced basal cell carcinoma in adults).",
    n_observations = 4942L,
    cohort_split   = "Patients 204 (90.7%) / healthy volunteers 21 (9.3%); Phase I formulation 55 (24.4%) / Phase II formulation 170 (75.6%).",
    studies        = c("SHH3925g (NCT00607724)", "SHH4610g (NCT00968981)", "SHH4476g (NCT00833417)", "SHH4433g", "SHH4683g (NCT00991718)"),
    notes          = "Population pooled across five Phase I / Phase II clinical studies (Lu 2015 Methods / Results). External validation used SHH4871g (NCT01173536), a dedicated QT study with 21 HVs receiving 150 mg QD Phase II formulation (470 total concentrations). 47.9% of AAG records were missing across all data (11.4% after excluding insensitive sub-day-sampling and HV data) and were imputed by LOCF. Mean baseline AAG: 31.1 uM in patients, 20.2 uM in HVs."
  )

  ini({
    # Structural PK parameters. Reference covariate set: AGE = 60 y, WT = 75 kg,
    # DIS_HEALTHY = 0 (cancer patient), FORM_VISMO_PHASEI = 0 (Phase II / commercial
    # formulation; F fixed at 1). AAG enters as a structural input, not as a
    # multiplicative covariate on a PK parameter.
    lka      <- log(9.025);  label("Absorption rate constant (1/day) for Phase II formulation in patients (reference)")     # Lu 2015 Table 2: exp(h4) = 9.025 day^-1
    lcl      <- log(1326);   label("Apparent clearance of unbound vismodegib (L/day) at AGE = 60 y (reference)")             # Lu 2015 Table 2: exp(h1) = 1326 L/day
    lvc      <- log(58.0);   label("Apparent central volume of distribution (L) at WT = 75 kg (reference)")                  # Lu 2015 Table 2: exp(h2) = 58.0 L
    lkdAAG   <- log(0.056);  label("Dissociation constant for fast-equilibrium binding to AAG (uM)")                          # Lu 2015 Table 2: exp(h3) = 0.056 uM

    # Covariate effects on disposition (Lu 2015 Eq. 5)
    e_age_cl <- -0.527;      label("Power exponent of AGE on CLunbound (unitless; reference 60 y)")                           # Lu 2015 Table 2: h9 = -0.527
    e_wt_vc  <-  0.660;      label("Power exponent of WT on Vc (unitless; reference 75 kg)")                                  # Lu 2015 Table 2: h10 = 0.660

    # Covariate effects on absorption (Lu 2015 Eq. 4):
    # ka = exp(lka + e_healthy_ka * DIS_HEALTHY + e_form_phaseI_ka * FORM_VISMO_PHASEI)
    e_healthy_ka      <-  0.671;  label("Linear-on-log-scale effect of HV on ka (additive on log(ka); unitless)")             # Lu 2015 Table 2: h5 =  0.671
    e_form_phaseI_ka  <- -0.602;  label("Linear-on-log-scale effect of Phase I formulation on ka (additive on log(ka))")      # Lu 2015 Table 2: h8 = -0.602

    # Covariate effects on relative bioavailability F (Lu 2015 Eq. 4):
    # F = 1 for FORM_VISMO_PHASEI = 0 (Phase II / commercial reference); F = exp(e_form_phaseI_f + e_healthy_f * DIS_HEALTHY)
    # for FORM_VISMO_PHASEI = 1.
    e_form_phaseI_f   <- -1.061;  label("Linear-on-log-scale effect of Phase I formulation on F in patients (= log(0.346))")  # Lu 2015 Table 2: h6 = log(0.346) = -1.061
    e_healthy_f       <-  0.881;  label("Additional linear-on-log-scale effect of HV on F when Phase I formulation")          # Lu 2015 Table 2: h7 =  0.881

    # Inter-individual variability. Lu 2015 Table 2 reports IIV as %CV
    # (intersubject variability); the stored log-normal variance follows
    # omega^2 = log(CV^2 + 1). CLunbound 48.7% CV -> 0.21278;
    # Vc 45.5% CV -> 0.18804; correlation(CLunbound, Vc) = 0.55 (Results
    # paragraph after Eq. 6: "...pairwise correlation of 0.55 for CLunbound
    # vs. Vc...") so cov = 0.55 * sqrt(0.21278 * 0.18804) = 0.11003.
    # KDAAG 19.7% CV -> 0.03810. No IIV on ka or F per Lu 2015
    # Results: "Intersubject variability ... was not estimated for the
    # absorption parameters (F and ka) due to the limited data in the
    # absorption phase."
    etalcl + etalvc ~ c(0.21278, 0.11003, 0.18804)                # Lu 2015 Table 2: IIV CLunbound 48.7% CV, Vc 45.5% CV, correlation 0.55
    etalkdAAG       ~ 0.03810                                     # Lu 2015 Table 2: IIV KDAAG 19.7% CV

    # Residual error. Lu 2015 Methods: "Natural log-transformed data were
    # used for modeling. ... An additive error model on the log-transformed
    # data was applied." Additive on log scale ~ proportional in linear
    # space (nlmixr2 prop()). Lu 2015 Table 2 reports the linear-scale CV%
    # for each output: 26.7% for total vismodegib plasma concentration and
    # 42.4% for unbound vismodegib plasma concentration.
    propSd    <- 0.267;   label("Proportional residual error on total vismodegib Cc (fraction)")   # Lu 2015 Table 2: total 26.7% CV
    propSd_Cu <- 0.424;   label("Proportional residual error on unbound vismodegib Cu (fraction)") # Lu 2015 Table 2: unbound 42.4% CV
  })

  model({
    # Vismodegib molecular weight (g/mol). Used to convert input dose
    # records (mg) into the depot compartment in micromoles so that
    # concentrations come out in uM (matching the units of KDAAG and AAG
    # reported by Lu 2015). 421.30 g/mol from PubChem CID 24776
    # (GDC-0449; C19H14Cl2N2O3S).
    mw_vismo <- 421.30

    # Individual PK parameters with covariate effects (Lu 2015 Eqs 4, 5):
    #   ka          = exp(h4 + h5 * HV + h8 * Phase_I)
    #   CLunbound_i = exp(h1 + h9 * log(AGE/60)) * exp(etalcl)
    #   Vc_i        = exp(h2 + h10 * log(WT/75)) * exp(etalvc)
    #   KDAAG_i     = exp(h3) * exp(etalkdAAG)
    ka         <- exp(lka + e_healthy_ka * DIS_HEALTHY + e_form_phaseI_ka * FORM_VISMO_PHASEI)
    cl_unbound <- exp(lcl + etalcl) * (AGE / 60)^e_age_cl
    vc         <- exp(lvc + etalvc) * (WT  / 75)^e_wt_vc
    kdAAG      <- exp(lkdAAG + etalkdAAG)

    # Relative bioavailability F (Lu 2015 Eq. 4):
    #   F = 1                                          if FORM_VISMO_PHASEI = 0 (Phase II reference)
    #   F = exp(h6 + h7 * DIS_HEALTHY)                  if FORM_VISMO_PHASEI = 1 (Phase I formulation)
    # The single closed-form expression below yields F = 1 when
    # FORM_VISMO_PHASEI = 0 (because exp(0) = 1) and the Phase-I-specific
    # value otherwise.
    F_rel <- exp(FORM_VISMO_PHASEI * (e_form_phaseI_f + e_healthy_f * DIS_HEALTHY))

    # Total drug concentration in central compartment (uM = umol/L; the
    # depot/central amounts carry umol because f(depot) converts the input
    # mg dose to umol via the molecular-weight scaling below).
    Cc <- central / vc

    # Saturable fast-equilibrium binding to AAG (Lu 2015 Methods, Structural
    # model; Discussion paragraph citing Widmer 2006 imatinib popPK):
    #   Cc = Cu + (AAG * Cu) / (KDAAG + Cu)
    # Solve the resulting quadratic Cu^2 + (KDAAG + AAG - Cc) * Cu - KDAAG * Cc = 0
    # for the unbound concentration Cu. Both AAG and KDAAG are in uM, so Cu
    # comes out in uM. The +sqrt(...) branch is the physically meaningful
    # (non-negative) root.
    bind_b <- kdAAG + AAG - Cc
    Cu     <- 0.5 * (sqrt(bind_b * bind_b + 4 * kdAAG * Cc) - bind_b)

    # ODE system: depot -> central, with elimination acting on UNBOUND drug
    # (cl_unbound is the apparent clearance of unbound vismodegib; Cu is
    # the unbound concentration in uM). The amounts in depot and central
    # carry umol because of the f(depot) unit conversion below.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl_unbound * Cu

    # Bioavailability + unit conversion. Input AMT is in mg; f(depot) scales
    # by F_rel and by (1000 / MW_vismo) to convert mg -> umol so that depot
    # carries umol throughout. 1000 / 421.30 ~ 2.373 umol/mg.
    f(depot) <- F_rel * 1000 / mw_vismo

    # Residual error. The paper applied an additive error on log-transformed
    # data, which is equivalent to a proportional error in linear space for
    # small CV; for the unbound CV of 42.4% the linear-scale and log-scale
    # forms diverge slightly (see vignette Assumptions and deviations).
    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
