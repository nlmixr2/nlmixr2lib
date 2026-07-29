Budha_2015_rg7652 <- function() {
  description <- "Population PK/PD model for RG7652 (an anti-PCSK9 monoclonal antibody) in healthy hypercholesterolemic subjects (Budha 2015): one-compartment PK with first-order SC absorption and combined linear plus Michaelis-Menten elimination from the central compartment, linked to a Type 3 indirect-response model for serum low-density lipoprotein cholesterol (LDL-C) in which RG7652 stimulates LDL-C degradation through an Emax function."
  reference   <- "Budha NR, Leabman M, Jin JY, Wada DR, Baruch A, Peng K, Tingley WG, Davis JD. Modeling and Simulation to Support Phase 2 Dose Selection for RG7652, a Fully Human Monoclonal Antibody Against Proprotein Convertase Subtilisin/Kexin Type 9. AAPS J. 2015;17(4):881-890. doi:10.1208/s12248-015-9750-8"
  vignette    <- "Budha_2015_rg7652"
  paper_specific_compartments <- c("LDL")

  units       <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Power-form effect on Ka with reference 47 years (Budha 2015 Table I footnote: 'Typical value of Ka = 0.348 * (AGE/47)^-0.886'). Older subjects have a lower absorption rate.",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power-form (allometric) scaling on linear CL and central volume of distribution V with reference 80 kg (Budha 2015 Table I footnotes: 'Typical value of CL = 0.426 * (BW/80)^0.813' and 'Typical value of V1 = 8.38 * (BW/80)^0.288').",
      source_name        = "BW"
    ),
    LDLC = list(
      description        = "Subject's natural (pre-treatment, screening) baseline serum LDL-C concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Two roles in this model: (1) initialises the indirect-response LDL state (LDL(0) <- LDLC) and anchors the zero-order production Ksyn = Kdeg * LDLC so the system is at LDLC at steady state in the absence of drug and statin; (2) drives a power-form effect on Emax with reference 145 mg/dL (Budha 2015 Results: 'Baseline (pre-dose) LDLc values ranged from 58 to 242 mg/dL (median, 145 mg/dL)' and Table II 'Baseline ~ Emax' covariate). Subjects with a higher natural baseline LDL-C have a lower Emax. For statin-pretreated subjects in the Phase 1 study the LDLC entered into the model is the screening LDL-C measured 3 to 4 weeks prior to RG7652 and immediately prior to atorvastatin administration (Budha 2015 Methods: 'screening LDLc values collected prior to atorvastatin dosing were used as the initial condition for the indirect response model'); the time-varying CONMED_STATIN_MONO covariate then dynamically pulls the LDL state to the lower equilibrium during the atorvastatin run-in and accounts for the LDL rebound after statin cessation.",
      source_name        = "LDLC"
    ),
    CONMED_STATIN_MONO = list(
      description        = "Concomitant atorvastatin monotherapy indicator, 1 = subject is receiving atorvastatin 40 mg daily and no other lipid-lowering comedication, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant statin)",
      notes              = "Time-varying within subject in the Phase 1 study: atorvastatin 40 mg daily started >=3 weeks before the first RG7652 dose and continued 5 weeks after, then stopped on day 35 (Budha 2015 Methods, Phase I Study). Multiplicative exponential effect on the indirect-response zero-order production rate Ksyn: Ksyn(t) = Kdeg * LDLC * exp(theta * CONMED_STATIN_MONO(t)) with theta = -0.648 (Budha 2015 Table II 'Statin on screening LDLc'). The same parameter accounts for the ~45% lower pre-dose baseline LDL-C in the statin cohorts and the LDL rebound after statin cessation (Budha 2015 Discussion). Reference patient is not on statin therapy.",
      source_name        = "STATIN"
    )
  )

  population <- list(
    species          = "human",
    n_subjects_pk    = 60L,
    n_subjects_pd    = 80L,
    n_studies        = 1L,
    phase            = "Phase 1 (NCT-not-reported in main text; Phase 2 NCT01609140 was the simulation target)",
    age_range        = "19-64 years (median 46.5 years)",
    weight_range     = "52.6-114.5 kg (median 82.6 kg)",
    sex_female_pct   = 47.5,
    race_ethnicity   = "Not reported in the main paper.",
    disease_state    = "Healthy adults with elevated LDL-C (130-220 mg/dL) at screening; not patients with coronary heart disease (the Phase 2 simulation target population).",
    dose_range       = "Single-dose cohorts: 10, 40, 150, 300, 600, or 800 mg SC. Multiple-dose cohorts: 40 or 150 mg SC once weekly for 4 weeks (QWx4) with or without concomitant atorvastatin 40 mg daily.",
    regions          = "Not reported in the main paper.",
    baseline_ldlc_screening_range = "125-220 mg/dL (median 163 mg/dL) at screening (pre-statin).",
    baseline_ldlc_predose_range   = "58-242 mg/dL (median 145 mg/dL) at pre-dose (post-statin for statin-treated cohorts).",
    concomitant      = "CONMED_STATIN_MONO = 1 only for cohorts receiving atorvastatin 40 mg daily as the sole lipid-lowering therapy.",
    samples_pk       = "687 evaluable RG7652 serum concentrations measured by validated target-binding ELISA (biotin-rhuPCSK9 capture, anti-CDR detection).",
    samples_pd       = "1070 serum LDL-C concentrations measured by the Roche LDL Cholesterol Plus 2nd Generation direct method.",
    notes            = "Sequential modelling: individual post-hoc PK parameters from the PK fit were used to drive the PD model. PK fit n = 60 (RG7652-treated subjects only); PD fit n = 80 (all subjects including placebo). The Table II PD parameter set tabulated here is the model in which the statin covariate acts on the screening-LDLc / Ksyn anchor (the main-text published estimates). A sensitivity 'statin model' that adds a second statin effect directly on Emax is described in the paper's Discussion and supplemental material but is not coded here."
  )

  ini({
    # ----- PK structural parameters (Budha 2015 Table I final-model typical values) -----
    # Reference covariates: AGE = 47 years (median rounded from 46.5), WT = 80 kg
    # (rounded from median 82.6 kg per Table I footnote).
    # Units: time = day, dose = mg, concentration = ug/mL (= mg/L).
    lka   <- log(0.348);  label("First-order SC absorption rate Ka at 47-year reference (1/day)")               # Budha 2015 Table I final-model Ka = 0.348
    lcl   <- log(0.426);  label("Apparent linear clearance CL/F at 80-kg reference (L/day)")                    # Budha 2015 Table I final-model CL/F = 0.426
    lvc   <- log(8.38);   label("Apparent central volume of distribution V/F at 80-kg reference (L)")           # Budha 2015 Table I final-model V/F = 8.38
    lvmax <- log(1.91);   label("Maximum saturable elimination rate Vmax (mg/day)")                              # Budha 2015 Table I final-model Vmax = 1.91
    lkm   <- log(4.29);   label("Michaelis-Menten constant Km (ug/mL)")                                          # Budha 2015 Table I final-model Km = 4.29

    # ----- PK covariate effects (Budha 2015 Table I covariate block, power form) -----
    e_age_ka <- -0.886;   label("Power exponent of (AGE / 47 years) on Ka (unitless)")                          # Budha 2015 Table I Age ~ Ka = -0.886
    e_wt_cl  <-  0.813;   label("Power exponent of (WT / 80 kg) on linear CL (unitless)")                        # Budha 2015 Table I BW ~ CL/F = 0.813
    e_wt_vc  <-  0.288;   label("Power exponent of (WT / 80 kg) on central V (unitless)")                        # Budha 2015 Table I BW ~ V1 = 0.288

    # ----- PK inter-individual variability (Budha 2015 Table I BSV column, variances) -----
    # No IIV reported on Vmax or Km ('NE = not estimated' in Table I).
    etalka ~ 0.257                                                                                              # Budha 2015 Table I BSV on Ka = 0.257
    etalcl ~ 0.143                                                                                              # Budha 2015 Table I BSV on CL/F = 0.143
    etalvc ~ 0.0311                                                                                             # Budha 2015 Table I BSV on V/F = 0.0311

    # ----- PK residual error (Budha 2015 Table I) -----
    propSd <- 0.178; label("Proportional PK residual error on Cc (fraction)")                                   # Budha 2015 Table I Proportional error = 17.8%

    # ----- PD structural parameters (Budha 2015 Table II final-model typical values) -----
    # Type 3 indirect-response model: drug stimulates LDL-C degradation rate Kdeg.
    # LDL-C concentration units: mg/dL.
    lkdeg <- log(0.0476); label("LDL-C first-order degradation rate Kdeg (1/day)")                              # Budha 2015 Table II Kdeg = 0.0476
    lemax <- log(1.95);   label("Maximum drug-induced stimulation factor on Kdeg at 145 mg/dL baseline-LDLC reference (unitless)")  # Budha 2015 Table II Emax = 1.95
    lec50 <- log(13.8);   label("RG7652 serum concentration producing 50% of Emax (ug/mL)")                     # Budha 2015 Table II EC50 = 13.8
    lldlc <- fixed(log(145)); label("Reference (median pre-dose baseline) LDL-C used as the Emax-covariate reference (mg/dL)")  # Budha 2015 Results: 'Baseline (pre-dose) LDLc values ranged from 58 to 242 mg/dL (median, 145 mg/dL)'

    # ----- PD covariate effects (Budha 2015 Table II covariate block) -----
    # CONMED_STATIN_MONO is time-varying (atorvastatin pretreatment + 5 weeks; off thereafter).
    # The log-effect captures both the ~45% lower pre-dose baseline LDL-C during the
    # statin run-in (paper: 'about 45% lower'; exp(-0.648) ~= 0.523) and the LDL rebound
    # after statin cessation (paper: 'rebound phenomenon ... approaching the screening levels').
    e_conmed_statin_mono_ldlc <- -0.648; label("Log-effect of CONMED_STATIN_MONO on the LDL-C Ksyn anchor (unitless)")  # Budha 2015 Table II Statin on screening LDLc = -0.648
    e_ldlc_emax               <- -1.49;  label("Power exponent of (LDLC / 145 mg/dL) on Emax (unitless)")              # Budha 2015 Table II Baseline ~ Emax = -1.49

    # ----- PD inter-individual variability (Budha 2015 Table II OMEGA block, variances) -----
    etalkdeg ~ 0.162                                                                                            # Budha 2015 Table II omega^2_Kdeg = 0.162
    etalemax ~ 0.159                                                                                            # Budha 2015 Table II omega^2_Emax = 0.159
    etalldlc ~ 0.0209                                                                                           # Budha 2015 Table II omega^2_Baseline = 0.0209

    # ----- PD residual error (Budha 2015 Table II) -----
    propSd_LDL <- 0.131; label("Proportional PD residual error on LDL-C (fraction)")                            # Budha 2015 Table II Proportional residual error = 13.1%
  })
  model({
    # ----- Individual PK parameters (Budha 2015 Table I reference AGE = 47 y, WT = 80 kg) -----
    ka   <- exp(lka + etalka) * (AGE / 47)^e_age_ka
    cl   <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # ----- PK ODE system (Budha 2015 Equations in Methods 'Population PK and PD Models') -----
    # One-compartment SC-input disposition with combined linear (CL) and saturable
    # (Vmax / Km) elimination from the central compartment. Cc is the algebraic central
    # concentration (ug/mL = mg/L). Vmax has units mg/day, Km has units ug/mL, so the
    # saturable term Vmax * Cc / (Km + Cc) carries units of mg/day (a mass elimination rate).
    Cc <- central / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl * Cc - vmax * Cc / (km + Cc)

    # ----- Individual PD parameters (Budha 2015 Table II reference LDLC = 145 mg/dL) -----
    # LDLC is the subject's natural (pre-treatment, screening) baseline LDL-C value.
    # IIV on baseline is applied multiplicatively to the input covariate.
    kdeg     <- exp(lkdeg + etalkdeg)
    emax_ind <- exp(lemax + etalemax) * (LDLC / exp(lldlc))^e_ldlc_emax
    ec50     <- exp(lec50)
    ldlc_bl  <- LDLC * exp(etalldlc)

    # Zero-order LDL production rate Ksyn anchored to natural baseline; the time-varying
    # CONMED_STATIN_MONO multiplier reproduces the lower pre-dose baseline LDL during the
    # atorvastatin run-in and the rebound after cessation. RG7652 stimulates Kdeg via an
    # Emax function (Budha 2015 indirect-response equation in Methods).
    ksyn <- kdeg * ldlc_bl * exp(e_conmed_statin_mono_ldlc * CONMED_STATIN_MONO)

    LDL(0) <- ldlc_bl
    d/dt(LDL) <- ksyn - kdeg * (1 + emax_ind * Cc / (ec50 + Cc)) * LDL

    # ----- Observation and error models -----
    # Cc reported in ug/mL (= mg/L); LDL reported in mg/dL.
    Cc  ~ prop(propSd)
    LDL ~ prop(propSd_LDL)
  })
}
