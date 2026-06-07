Chang_2011_PF04455242_rat <- function() {
  description <- "Preclinical (Sprague-Dawley rat). Competitive antagonism PK-PD model of PF-04455242 (kappa opioid receptor antagonist) on spiradoline-induced plasma prolactin elevation. One-compartment first-order absorption PK for both spiradoline (KOR agonist challenge) and PF-04455242, with a dose-dependent absorption rate constant for PF-04455242 (1.64 /h at 3.2 mg/kg SC, 0.385 /h at 10 mg/kg SC). Direct-response sigmoid Emax PD: prolactin = baseline + Emax * Csp^gamma / ((EC50 * (1 + Cpf/Ki))^gamma + Csp^gamma) with competitive antagonism of the spiradoline-induced rise by PF-04455242. Spiradoline plasma compartments are declared via paper_specific_compartments rather than registering a new sibling-drug suffix; see Errata in the vignette for the rationale."
  reference <- paste(
    "Chang C, Byon W, Lu Y, Jacobsen LK, Badura LL, Sawant-Basak A, Miller E,",
    "Liu J, Grimwood S, Wang EQ, Maurer TS.",
    "Quantitative PK-PD Model-Based Translational Pharmacology of a Novel",
    "Kappa Opioid Receptor Antagonist Between Rats and Humans.",
    "AAPS J. 2011;13(4):565-575. doi:10.1208/s12248-011-9296-3.",
    sep = " "
  )
  vignette <- "Chang_2011_PF04455242_prolactin"
  paper_specific_compartments <- c("depot_spiradoline", "central_spiradoline")

  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Per-subject assigned PF-04455242 SC dose level (mg/kg). Used as a threshold-based switch for the dose-dependent PF-04455242 absorption rate constant: DOSE > 5 mg/kg selects the high-dose Ka (0.385 /h), DOSE <= 5 selects the low-dose Ka (1.64 /h). Spiradoline dose level does NOT enter this switch (spiradoline has its own depot/central compartments and a single dose-independent Ka).",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Set to the rat's PF-04455242 dose level (0.32, 3.2, or 10 mg/kg SC in Chang 2011 study 3). For simulations of the PD study, the typical levels are 3.2 mg/kg (low) and 10 mg/kg (high) per Fig. 4b. The eta on Ka was reported only at 3.2 mg/kg (NE at 10 mg/kg per Table I footnote b) and is applied uniformly here.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 9L,
    n_studies      = 3L,
    age_range      = "adult",
    weight_range   = "300-400 g",
    sex_female_pct = 0,
    disease_state  = "healthy male rats",
    dose_range     = paste(
      "Study 1 (spiradoline dose-response for prolactin): spiradoline 0, 0.1, 0.32, or 1.0 mg/kg SC (n=6-9 per group).",
      "Study 2 (spiradoline PK satellite): spiradoline 0.1, 0.32, or 1.0 mg/kg SC, single timepoint at 30 min (n=5).",
      "Study 3 (PF-04455242 + spiradoline challenge): PF-04455242 0 (vehicle), 3.2, or 10 mg/kg SC at t=-30 min followed by spiradoline 0.32 mg/kg SC at t=0 (n=6-9 per group). Samples at -30, 10, 20, 30, 45, 60, 90, 120 min post-spiradoline.",
      sep = " "
    ),
    regions        = "Pfizer Inc., USA (Institutional Animal Care and Use Committee protocols)",
    notes          = "Three independent rat studies pooled for the sequential PK-PD fit: spiradoline PK estimated from study 2, PF-04455242 PK estimated from study 3 PK samples, then PD parameters estimated from studies 1 and 3 prolactin observations with the PK parameters fixed at their Stage 1/2 estimates (Methods 'Preclinical PK-PD Modeling')."
  )

  ini({
    # ------------------------------------------------------------------
    # PF-04455242 PK -- one-compartment first-order absorption.
    # Apparent parameters per kg body weight (CL/F, Vc/F). Dose-dependent
    # absorption rate constant: Ka = 1.64 /h at the 3.2 mg/kg PF-04455242
    # dose and 0.385 /h at the 10 mg/kg dose (Chang 2011 Table I, footnotes
    # a and b). Both Ka values are estimated; the IIV omega-Ka was reported
    # only at 3.2 mg/kg (NE at 10 mg/kg) and is applied uniformly below.
    # ------------------------------------------------------------------
    lka_pf_low  <- log(1.64);  label("PF-04455242 absorption rate constant at 3.2 mg/kg SC (1/h)")  # Table I, footnote a
    lka_pf_high <- log(0.385); label("PF-04455242 absorption rate constant at 10 mg/kg SC (1/h)")   # Table I, footnote b
    lcl         <- log(2.5);   label("PF-04455242 apparent oral clearance CL/F per kg (L/h/kg)")    # Table I
    lvc         <- log(1.31);  label("PF-04455242 apparent central volume Vc/F per kg (L/kg)")      # Table I

    # IIV on PF-04455242 PK. omega translated from reported CV% via
    # omega^2 = log(CV^2 + 1) (lognormal exponential random-effect model).
    etalka_pf_low ~ 0.0969  # Table I, omega-Ka at 3.2 mg/kg = 31.9%CV -> log(0.319^2+1) = 0.0969
    etalcl        ~ 0.1193  # Table I, omega-CL  = 35.6%CV -> log(0.356^2+1) = 0.1193
    etalvc        ~ 0.0659  # Table I, omega-Vc  = 26.1%CV -> log(0.261^2+1) = 0.0659

    propSd <- 0.141; label("PF-04455242 proportional residual error (fraction)")  # Table I, 14.1%

    # ------------------------------------------------------------------
    # Spiradoline PK -- one-compartment first-order absorption. KOR
    # agonist used as a challenge drug. Apparent parameters per kg body
    # weight. Compartment names depot_spiradoline / central_spiradoline
    # are declared via paper_specific_compartments above; corresponding
    # parameter names use the _spiradoline suffix.
    # ------------------------------------------------------------------
    lka_spiradoline <- log(4.75);  label("Spiradoline absorption rate constant (1/h)")                 # Table I
    lcl_spiradoline <- log(6.06);  label("Spiradoline apparent oral clearance CL/F per kg (L/h/kg)")   # Table I
    lvc_spiradoline <- log(7.15);  label("Spiradoline apparent central volume Vc/F per kg (L/kg)")     # Table I

    etalka_spiradoline ~ 0.2992  # Table I, omega-Ka = 59.1%CV -> log(0.591^2+1) = 0.2992
    etalcl_spiradoline ~ 0.0319  # Table I, omega-CL = 18.0%CV -> log(0.180^2+1) = 0.0319
    etalvc_spiradoline ~ 0.1957  # Table I, omega-Vc = 46.5%CV -> log(0.465^2+1) = 0.1957

    propSd_Cspir <- 0.4472; label("Spiradoline plasma proportional residual error (fraction)")  # Table I, 44.72%

    # ------------------------------------------------------------------
    # Prolactin PD -- direct-response sigmoid Emax with competitive
    # antagonism by PF-04455242 (Chang 2011 Eq. 1, Fig. 2). Baseline,
    # Emax, EC50 (spiradoline), Ki (PF-04455242), and the Hill exponent
    # gamma. PD was fit with the PK parameters fixed at the spiradoline
    # and PF-04455242 PK estimates above.
    # ------------------------------------------------------------------
    lrbase <- log(1.34); label("Baseline plasma prolactin (ng/mL)")                                       # Table II
    lemax  <- log(39.7); label("Maximal spiradoline-induced prolactin elevation above baseline (ng/mL)")  # Table II
    lec50  <- log(34.3); label("Spiradoline concentration giving half-maximal prolactin elevation (ng/mL)")  # Table II
    lki    <- log(414);  label("PF-04455242 in vivo Ki at rat KOR (ng/mL)")                               # Table II
    lhill  <- log(4.15); label("Hill coefficient (paper gamma) on the spiradoline stimulation (unitless)") # Table II

    etalrbase ~ 0.1422  # Table II, omega-BL = 39.1%CV -> log(0.391^2+1) = 0.1422

    propSd_PRL <- 1.31; label("Rat prolactin proportional residual error (fraction)")  # Table II, 131%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (typical-value * exp(eta) form).
    # PF-04455242 Ka switches on DOSE > 5 mg/kg via the difference of
    # the two log-Ka typical values; the same eta_low is applied to
    # both dose levels (paper reported NE at 10 mg/kg, so this is a
    # documented deviation - see vignette Assumptions and deviations).
    # ------------------------------------------------------------------
    ka_pf <- exp(lka_pf_low + (lka_pf_high - lka_pf_low) * (DOSE > 5) + etalka_pf_low)
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    kel   <- cl / vc

    ka_sp  <- exp(lka_spiradoline + etalka_spiradoline)
    cl_sp  <- exp(lcl_spiradoline + etalcl_spiradoline)
    vc_sp  <- exp(lvc_spiradoline + etalvc_spiradoline)
    kel_sp <- cl_sp / vc_sp

    # ------------------------------------------------------------------
    # 2. PD typical-value parameters (typical-value * exp(eta)).
    # ------------------------------------------------------------------
    rbase <- exp(lrbase + etalrbase)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)
    ki    <- exp(lki)
    hill  <- exp(lhill)

    # ------------------------------------------------------------------
    # 3. PK ODEs. PF-04455242 uses canonical depot / central. Spiradoline
    # uses paper-specific depot_spiradoline / central_spiradoline.
    # Subjects receive a PF-04455242 SC dose at t=-0.5 h (cmt=depot,
    # amt in mg/kg) and a spiradoline SC dose at t=0 (cmt=depot_spiradoline,
    # amt=0.32 mg/kg in Chang 2011 study 3).
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka_pf * depot
    d/dt(central) <-  ka_pf * depot - kel * central

    d/dt(depot_spiradoline)   <- -ka_sp * depot_spiradoline
    d/dt(central_spiradoline) <-  ka_sp * depot_spiradoline - kel_sp * central_spiradoline

    # ------------------------------------------------------------------
    # 4. Plasma concentrations and PD response. Internal compartment
    # state is in mg/kg (dose unit) and Vc/F is in L/kg, so central/vc
    # has units mg/L = ug/mL. Multiply by 1000 to display in ng/mL,
    # matching the paper's reporting (Tables I/II, Figures 3-4) and
    # keeping Cc / Cspir on the same scale as the PD parameters
    # (rbase, emax, ec50, ki are all expressed in ng/mL).
    # Cc      -- PF-04455242 plasma concentration (parent observation, ng/mL).
    # Cspir   -- spiradoline plasma concentration (paper-named output, ng/mL).
    # PRL     -- prolactin response (paper-named PD output, ng/mL): Chang
    #            2011 Eq. 1. PF-04455242 acts as a competitive antagonist
    #            by shifting the apparent spiradoline EC50 upward via the
    #            (1 + Cpf/Ki) factor.
    # ------------------------------------------------------------------
    Cc    <- central / vc * 1000
    Cspir <- central_spiradoline / vc_sp * 1000

    ec50_apparent <- ec50 * (1 + Cc / ki)
    stim          <- emax * Cspir^hill / (ec50_apparent^hill + Cspir^hill)
    PRL           <- rbase + stim

    Cc    ~ prop(propSd)
    Cspir ~ prop(propSd_Cspir)
    PRL   ~ prop(propSd_PRL)
  })
}
