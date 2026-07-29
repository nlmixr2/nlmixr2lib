Sassen_2017_crisantaspase <- function() {
  description <- "Two-compartment population PK model for intravenous Erwinia asparaginase (crisantaspase; Erwinase) in pediatric acute lymphoblastic leukemia patients, with allometric scaling on clearance and volumes and a higher first-month clearance (Sassen 2017)."
  reference   <- "Sassen SDT, Mathot RAA, Pieters R, Kloos RQH, de Haas V, Kaspers GJL, et al. Population pharmacokinetics of intravenous Erwinia asparaginase in pediatric acute lymphoblastic leukemia patients. Haematologica. 2017;102(3):552-561. doi:10.3324/haematol.2016.149195"
  vignette    <- "Sassen_2017_crisantaspase"
  units       <- list(time = "hour", dosing = "IU", concentration = "IU/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Treated as a time-fixed baseline in Sassen 2017 (pediatric ALL cohort, median 24.5 kg, range 11.7-99.0 kg; Table 1). Drives the allometric scaling on CL, Q (exponent 0.75) and Vc, Vp (exponent 1) with reference 70 kg.",
      source_name        = "WT"
    ),
    MONTH1 = list(
      description        = "Indicator for the first 30 days of Erwinia asparaginase treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (month 2 onwards)",
      notes              = "Time-varying within subject; 1 if the observation falls within the first 30 days of treatment, 0 otherwise. Multiplicative effect on CL via (1 + e_month1_cl * MONTH1) gives a 14% higher CL during the first month of treatment (Sassen 2017 Results, Pharmacokinetic model section; abstract and Discussion). Data assemblers derive `MONTH1 = as.integer(time_post_treatment_start_days < 30)`.",
      source_name        = "MONTH1"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 51L,
    n_studies           = 1L,
    n_observations      = 714L,
    age_range           = "1.9-17.7 years",
    age_median          = "6 years",
    weight_range        = "11.7-99.0 kg",
    weight_median       = "24.5 kg",
    sex_female_pct      = 37.3,
    race_ethnicity      = "Not reported (multicenter Dutch cohort).",
    disease_state       = "Pediatric acute lymphoblastic leukemia treated per DCOG ALL-10 (Nov 2004-Apr 2012) or ALL-11 (since Apr 2012) protocols, after the development of an allergy to E. coli-derived asparaginase or silent inactivation of E. coli-derived asparaginase.",
    dose_range          = "Starting dose 20,000 IU/m^2 IV over 1 h thrice weekly (Mon/Wed/Fri), adjusted by therapeutic drug monitoring (TDM) from week 3 onwards based on 72 h trough concentrations (target >= 100 IU/L).",
    regions             = "Netherlands (7 pediatric oncology centers, Dutch Childhood Oncology Group).",
    bioanalytic_methods = "Erwinia asparaginase activity in serum, LLOQ < 5 IU/L.",
    notes               = "Baseline demographics from Sassen 2017 Table 1; sample summary from Table 2. NONMEM, 2-compartment linear elimination with allometric scaling (fixed exponents 0.75 / 1) on a 70 kg reference. IIV on CL only (Vc, Q, Vp IIV not estimable per paper). Inter-occasion variability (IOV) on CL of 13% CV on monthly intervals was reported in Table 3 but is NOT encoded in this ini() (nlmixr2lib popPK convention is to carry IPV only); see vignette Assumptions and deviations. ClinicalTrials.gov / NTR 3379 Dutch Trial Register."
  )

  ini({
    # Reference 70 kg subject; clearances in L/h, volumes in L. Sassen 2017
    # Table 3 final-model estimates: CL 0.44 L/h/70kg, Vc 3.2 L/70kg,
    # Q 0.15 L/h/70kg, Vp 2.9 L/70kg.
    lcl <- log(0.44); label("Clearance for a 70 kg subject (CL, L/h)")                        # Sassen 2017 Table 3 (CL 0.44 L/h/70kg)
    lvc <- log(3.2);  label("Central volume of distribution for a 70 kg subject (Vc, L)")     # Sassen 2017 Table 3 (Vc 3.2 L/70kg)
    lq  <- log(0.15); label("Intercompartmental clearance for a 70 kg subject (Q, L/h)")      # Sassen 2017 Table 3 (Q 0.15 L/h/70kg)
    lvp <- log(2.9);  label("Peripheral volume of distribution for a 70 kg subject (Vp, L)")  # Sassen 2017 Table 3 (Vp 2.9 L/70kg)

    # Allometric exponents held fixed at canonical values (Sassen 2017 Methods,
    # Pharmacokinetic analysis: "standard fixed exponent values of 0.75 for
    # the flow dependent physiological process parameters and 1 for volume-
    # related parameters").
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")   # Sassen 2017 Methods (Pharmacokinetic analysis)
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")  # Sassen 2017 Methods (Pharmacokinetic analysis)

    # First-month covariate effect on CL. Sassen 2017 Table 3 reports a
    # multiplier of 1.1 (rounded; bootstrap median 1.12; 95%CI 1.06-1.22;
    # RSE 3%) and the abstract / Discussion explicitly describe this as
    # "14% higher" CL in the first month. Encoded here as a fractional
    # change e_month1_cl = 0.14 applied via (1 + e_month1_cl * MONTH1) so
    # that CL_month1 = 1.14 * TVCL.
    e_month1_cl <- 0.14; label("Fractional increase in CL during the first month of treatment (unitless)")  # Sassen 2017 Table 3 ("CL month 1 diff" 1.14 multiplier) and abstract / Discussion ("14% higher")

    # IIV on CL only (Sassen 2017 Table 3: IIV CL = 33% CV; IIV on Vc, Q, Vp
    # was not estimable, per Results, Pharmacokinetic model section).
    # omega^2 = log(1 + CV^2) = log(1 + 0.33^2) = 0.10337.
    etalcl ~ 0.10337  # Sassen 2017 Table 3 (IIV CL 33% CV converted as log(1 + 0.33^2))

    # Residual error. Sassen 2017 Methods (Pharmacokinetic analysis) state
    # the analysis used log-transformed concentrations with an additive
    # error model. In NONMEM this is Y = LOG(IPRED) + EPS(1) with
    # EPS(1) ~ N(0, sigma^2); the equivalent nlmixr2 residual is a log-
    # normal error Cc ~ lnorm(expSd). Table 3 prints the residual SD as
    # 0.57 with a misleading "IU/L" unit header; on log scale it is
    # unitless.
    expSd <- 0.57; label("Log-scale residual SD (unitless)")  # Sassen 2017 Table 3 (residual error 0.57 on log scale)
  })

  model({
    # Multiplicative first-month CL factor: in the first 30 days of
    # treatment CL is increased by 14% (MONTH1 = 1); otherwise the factor
    # collapses to 1.
    cl_month1_factor <- 1 + e_month1_cl * MONTH1

    # Allometric weight scaling on a 70 kg reference (exponent 0.75 on CL
    # and Q, exponent 1 on Vc and Vp; Sassen 2017 Methods).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * cl_month1_factor
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment model with IV input (1 h infusion into the central
    # compartment; no extravascular depot; Sassen 2017 Methods).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in IU, vc in L -> central/vc in IU/L (matching
    # the assay output of asparaginase activity in IU/L).
    Cc <- central / vc

    # Log-normal residual error (Sassen 2017 fitted log-transformed
    # concentrations with an additive residual error on the log scale).
    Cc ~ lnorm(expSd)
  })
}
