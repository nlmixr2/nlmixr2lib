Gulati_2014_brownSnakeVenom <- function() {
  description <- paste(
    "QSP (reduced). Five-state 'extracted' coagulation model of Gulati 2014,",
    "obtained by proper lumping of the 62-state Wajima 2009 / Gulati 2013",
    "humoral coagulation network, fitted as a population PK-PD model to",
    "fibrinogen concentrations after Australian brown snake (Pseudonaja spp.)",
    "envenoming with complete venom-induced consumption coagulopathy (VICC).",
    "Venom is absorbed first order from a bite-site state (depot) into a",
    "plasma venom state (central) and eliminated first order; antivenom",
    "multiplies the venom elimination rate 100-fold from the time it is given",
    "(time-varying CONMED_ANTIVENOM). Plasma venom activates a lumped",
    "procoagulant pool (lumped, assumed dominated by prothrombin) to thrombin",
    "by a saturable venom-driven term, and thrombin consumes fibrinogen by a",
    "saturable thrombin-driven term; fibrinogen and the lumped pool have",
    "zero-order production and first-order loss. The venom 'dose' is a fixed",
    "0.0075 nmol/L initial amount in the bite-site state for every patient,",
    "with between-subject variability on its bioavailable fraction",
    "(typical value fixed to 1). IIV on thrombin loss, the venom-driven Vmax",
    "and Km, fibrinogen loss and venom bioavailability; combined additive +",
    "proportional residual error on fibrinogen in g/L."
  )
  reference <- paste(
    "Gulati A, Isbister GK, Duffull SB. Scale reduction of a systems",
    "coagulation model with an application to modeling",
    "pharmacokinetic-pharmacodynamic data. CPT Pharmacometrics Syst Pharmacol.",
    "2014;3(1):e90. doi:10.1038/psp.2013.67.",
    "Structure: Eq 1-5 and Figure 1; estimates: Table 2; original 62-state",
    "model: Supplementary Data original_model.m and Supplementary Table S1."
  )
  vignette <- "Gulati_2014_brownSnakeVenom"
  units <- list(time = "h", dosing = "nmol/L", concentration = "g/L")

  compartmentData <- list(
    depot = list(
      analyte = "brown snake venom (procoagulant activity)",
      units = "nmol/L",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "brown snake venom (procoagulant activity)",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    fibrinogen = list(analyte = "fibrinogen", units = "nmol/L", specimen = "plasma", verified = TRUE),
    thrombin = list(analyte = "thrombin (factor IIa)", units = "nmol/L", specimen = "plasma", verified = TRUE),
    lumped = list(
      analyte = "lumped procoagulant pool (58 original states, assumed dominated by prothrombin)",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONMED_ANTIVENOM = list(
      description = "Time-varying indicator that antivenom has been administered: 0 before the antivenom dose, 1 from the time of the antivenom dose onward.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no antivenom yet)",
      notes = "Gulati 2014 Eq 2: the venom plasma elimination term is multiplied by 'antivenom' = 1 if antivenom is absent and 100 if present; the deposited original_model.m switches it at 't > time_of_antivenom'. Every complete-VICC patient in the analysis received antivenom (Methods, 'The data'); the lumping simulations placed it at 4 h after the bite, the average in the brown snake data (Methods, 'Simplification of the blood coagulation systems pharmacology model'). Supply the column as a per-record value that switches from 0 to 1 at the subject's antivenom time.",
      source_name = "antivenom (Eq 2); time_of_antivenom (original_model.m)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61,
    n_studies = 1,
    age_range = "41 +/- 18 years (mean +/- SD, Table 3)",
    sex_female_pct = 26,
    disease_state = "Australian brown snake (Pseudonaja spp.) envenoming with complete venom-induced consumption coagulopathy (VICC); all patients received antivenom.",
    dose_range = "Single snakebite; venom dose unknown and set to a nominal 0.0075 (Eq 1 initial condition, nmol/L) for every patient, with between-subject variability on the bioavailable fraction.",
    regions = "Australia (Australian Snakebite Project, over 100 hospitals, January 2004 to May 2008).",
    notes = "Table 3: 73 patients bitten by brown snakes, 61 with complete VICC (the analysed fibrinogen data), 74% male. Fibrinogen limit of quantification 0.2 g/L; below-limit data handled with Beal's M3 in NONMEM 7.2 (FOCE-I). The wider Australian Snakebite Project cohort held 138 snakebite cases."
  )

  ini({
    # Gulati 2014 Table 2 (final model). Eq 1-5 give the structure.
    lka <- log(0.854); label("Venom absorption rate constant k_abs^brown (1/h)") # Table 2 k_a^Brown = 0.854 1/h (SE 24.2%)
    lkel <- log(0.744); label("Venom plasma elimination rate constant d^brown without antivenom (1/h)") # Table 2 d^Brown = 0.744 1/h (SE 0.87%)
    kel_antivenom_fold <- fixed(100); label("Fold increase of venom elimination once antivenom is given (unitless)") # Eq 2: antivenom = 100 if present, 1 if absent
    lfdepot <- fixed(log(1)); label("Bioavailable fraction of the nominal venom dose F_venom (unitless)") # Methods: fixed-effects value fixed to 1, BSV estimated

    lkdeg_thrombin <- log(58.3); label("Thrombin (factor IIa) first-order loss rate d^IIa (1/h)") # Table 2 d^IIa = 58.3 1/h (SE 4.1%)
    lvmax_lumped <- log(1940); label("Vmax of venom-driven lumped-state activation to thrombin (1/h as used in Eq 4-5)") # Table 2 vmax Lumped_state->IIa (venom) = 1940 (SE 10.5%)
    lkm_lumped <- log(1.38); label("Km (venom concentration) of lumped-state activation to thrombin (nmol/L)") # Table 2 Km Lumped_state->IIa (venom) = 1.38 nmol/L (SE 4.8%)
    kin_lumped <- fixed(147); label("Zero-order production rate of the lumped state p^Lumped_state (nmol/L/h)") # Results 'Identifiability': fixed at 147 nmol/l/h (= pII of the original model)
    kdeg_lumped <- fixed(0.01); label("First-order loss rate of the lumped state d^Lumped_state (1/h)") # Results 'Identifiability': fixed at 0.01/h (= dII of the original model)

    lkin_fibrinogen <- log(216); label("Zero-order fibrinogen production rate p^Fg (nmol/L/h)") # Table 2 p^Fg = 216 nmol/l/h (SE 6.8%)
    lkdeg_fibrinogen <- log(0.0174); label("First-order fibrinogen loss rate d^Fg (1/h)") # Table 2 d^Fg = 0.0174 1/h (SE 9.4%)
    lvmax_fibrinogen <- log(8.88); label("Vmax of thrombin-driven fibrinogen conversion to fibrin (1/h as used in Eq 3)") # Table 2 vmax Fg->F (IIa) = 8.88 (SE 4.5%)
    lkm_fibrinogen <- log(22.9); label("Km (thrombin concentration) of fibrinogen conversion to fibrin (nmol/L)") # Table 2 Km Fg->F (IIa) = 22.9 nmol/l (SE 5.7%)

    bl_fibrinogen <- fixed(8945.5); label("Fibrinogen initial condition C3L(0) (nmol/L)") # Eq 3: C3L(0) = 8945.5
    bl_lumped <- fixed(1394.4); label("Lumped-state initial condition C5L(0) (nmol/L)") # Eq 5: C5L(0) = 1394.4
    mw_fibrinogen <- fixed(340); label("Fibrinogen molecular weight converting nmol/L to g/L (kDa)") # not printed in the source; standard fibrinogen MW, see vignette Assumptions

    # Table 2 BSV rows are labelled '(%)' but carry the NONMEM OMEGA variances
    # of exponential random effects (see vignette Assumptions and deviations).
    etalkdeg_thrombin ~ 0.005 # Table 2 BSV_dIIa = 0.00500 (SE 10.0%, shrinkage 75%)
    etalvmax_lumped ~ 2.10 # Table 2 BSV_vmax Lumped_state->IIa = 2.10 (SE 23.0%, shrinkage 43%)
    etalkm_lumped ~ 0.990 # Table 2 BSV_Km Lumped_state->IIa = 0.990 (SE 24.5%, shrinkage 61%)
    etalkdeg_fibrinogen ~ 0.530 # Table 2 BSV_dFg = 0.530 (SE 17.8%, shrinkage 46%)
    etalfdepot ~ 1.40 # Table 2 BSV_Fvenom = 1.40 (SE 4.5%, shrinkage 54%)

    # Table 2 prints sigma^2_prop = 0.390. Read as an SD (39% CV) because only
    # that reading reproduces the Figure 3 VPC prediction interval; read as a
    # variance (SD 0.624) the simulated 10th percentile at 40 h is 0.4 g/L
    # against about 0.95 g/L in Figure 3. See vignette Assumptions.
    propSd <- 0.390; label("Proportional residual error on fibrinogen (fraction)") # Table 2 sigma^2_prop row = 0.390, used as the SD

    addSd <- sqrt(0.00100); label("Additive residual error on fibrinogen (g/L)") # Table 2 sigma^2_add = 0.00100 -> SD = sqrt(0.001)
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    fdepot <- exp(lfdepot + etalfdepot)
    kdeg_thrombin <- exp(lkdeg_thrombin + etalkdeg_thrombin)
    vmax_lumped <- exp(lvmax_lumped + etalvmax_lumped)
    km_lumped <- exp(lkm_lumped + etalkm_lumped)
    kin_fibrinogen <- exp(lkin_fibrinogen)
    kdeg_fibrinogen <- exp(lkdeg_fibrinogen + etalkdeg_fibrinogen)
    vmax_fibrinogen <- exp(lvmax_fibrinogen)
    km_fibrinogen <- exp(lkm_fibrinogen)

    # Eq 2: antivenom multiplies venom elimination by 100 once given.
    antivenom <- 1 + (kel_antivenom_fold - 1) * CONMED_ANTIVENOM

    # Venom-driven activation of the lumped pool to thrombin (Eq 4 and 5).
    r_thrombin <- vmax_lumped * lumped * central / (km_lumped + central)
    # Thrombin-driven conversion of fibrinogen to fibrin (Eq 3).
    r_fibrin <- vmax_fibrinogen * fibrinogen * thrombin / (km_fibrinogen + thrombin)

    d/dt(depot) <- -ka * depot # Eq 1
    d/dt(central) <- ka * depot - kel * central * antivenom # Eq 2
    d/dt(fibrinogen) <- kin_fibrinogen - kdeg_fibrinogen * fibrinogen - r_fibrin # Eq 3
    d/dt(thrombin) <- -kdeg_thrombin * thrombin + r_thrombin # Eq 4
    d/dt(lumped) <- kin_lumped - kdeg_lumped * lumped - r_thrombin # Eq 5

    # Eq 1: C1L(0) = 0.0075 is supplied as a dose of 0.0075 into depot.
    f(depot) <- fdepot
    fibrinogen(0) <- bl_fibrinogen
    lumped(0) <- bl_lumped

    # Fibrinogen was measured in g/L (Figures 2-4, LOQ 0.2 g/L).
    Cfibrinogen <- fibrinogen * mw_fibrinogen / 1e6
    Cfibrinogen ~ add(addSd) + prop(propSd)
  })
}
