Suksai_2026_levodopa <- function() {
  description <- "Deterministic simulation PK-PD framework for oral levodopa in Parkinson's disease (Suksai 2026). A gastrointestinal depot feeds a central plasma compartment that exchanges by first-order micro-constants with a peripheral tissue compartment (k12 / k21) and, in parallel, with a brain compartment (k13 / k31); elimination is first-order from plasma. An effect-site concentration equilibrates with the plasma concentration at rate ke0 and drives a sigmoidal Emax motor response on a TRS-like scale. The paper fits nothing: every value is a literature-informed nominal baseline used for in silico regimen comparison, so all parameters are fixed and the model carries no IIV and no residual error. NOTE: the executed model reproduced here is the rate-constant parameterisation evidenced by the baseline configuration and Figures 3-4, NOT the clearance-parameterised system printed as Equations 1-4, which references an unreported CLint and does not conserve mass; see the vignette Errata for the full deviation list and its verification."
  reference <- paste(
    "Suksai S, Suantai S, Wandee P.",
    "A computational pharmacokinetic-pharmacodynamic framework for simulating and comparing",
    "individualized levodopa dosing in Parkinson's disease.",
    "Front Pharmacol. 2026 May 29;17:1817435.",
    "doi:10.3389/fphar.2026.1817435",
    sep = " "
  )
  vignette <- "Suksai_2026_levodopa"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. `effect` holds a CONCENTRATION rather than an amount --
  # Suksai 2026 Equation 4 is written directly in Ce, with no effect-site
  # amount and no volume dividing it -- so its units are mg/L, not mg.
  compartmentData <- list(
    depot       = list(analyte = "levodopa", units = "mg",   specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levodopa", units = "mg",   specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "levodopa", units = "mg",   specimen = "tissue",              verified = TRUE),
    brain       = list(analyte = "levodopa", units = "mg",   specimen = "tissue",              verified = TRUE),
    effect      = list(analyte = "levodopa", units = "mg/L", specimen = "not applicable",      verified = TRUE)
  )

  # Suksai 2026 models NO covariate effects. Age, body weight and disease stage
  # are named as sources of cohort heterogeneity (Methods, "Simulation design")
  # and Table 3 varies ka and CL by illustrative case, but no covariate ever
  # enters an equation and no covariate coefficient is reported anywhere in the
  # paper. They are therefore documented here rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Suksai 2026 Methods ('Simulation design') sampled the virtual cohort over 50-80 kg and",
        "the narrative states that 'variations in body weight can modify distribution-related",
        "parameters', but no weight-on-parameter relationship is written down and no allometric",
        "exponent or centring weight is reported. Table 3's three illustrative cases differ in",
        "weight (60 / 75 / 70 kg) yet their differing ka and CL are presented as scenario",
        "selections, not as a weight function -- Case 3 is the LIGHTEST of the three cases with",
        "moderate weight but carries the HIGHEST ka and CL, so the table cannot be read as a",
        "weight relationship even qualitatively. Screened but not retained; not implemented.",
        sep = " "
      ),
      source_name        = "Weight (Suksai 2026 Table 3)"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort range 50-80 years (Suksai 2026 Methods, 'Simulation design'); Table 3 cases are",
        "60 / 72 / 50 years. No age effect on any parameter is reported or implemented.",
        sep = " "
      ),
      source_name        = "Age (Suksai 2026 Table 3)"
    ),
    HY_STAGE = list(
      description        = "Hoehn-Yahr Parkinson's disease severity stage",
      units              = "(stage)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Suksai 2026 sampled Hoehn-Yahr stages II-IV across the 100 virtual profiles (Methods,",
        "'Simulation design') and Table 3 labels the three illustrative cases by disease stage,",
        "but disease stage never enters an equation and no stage-specific parameter estimate is",
        "reported. Screened but not retained; not implemented. NOTE: HY_STAGE is documented here",
        "only as a screened-and-excluded covariate and is NOT registered as a canonical covariate",
        "column by this extraction -- registering it would require a model that actually uses it.",
        sep = " "
      ),
      source_name        = "Disease stage (Suksai 2026 Table 3)"
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 100,
    n_studies     = 0,
    age_range     = "50-80 years",
    weight_range  = "50-80 kg",
    disease_state = "Parkinson's disease, Hoehn-Yahr stages II-IV, including both stable responders and profiles with motor fluctuations",
    dose_range    = "100 mg orally every 4 h for six doses over 24 h (nominal baseline regimen); 100-300 mg per dose at 4-8 h intervals across the illustrative scenarios of Table 3",
    regions       = "not applicable (fully in silico)",
    notes         = paste(
      "NO PATIENT DATA. Suksai 2026 Methods ('Literature-based parameter selection'): 'No",
      "individual patient-level data were used; all simulations were conducted entirely in",
      "silico.' The 100 profiles are a virtual cohort generated by independent sampling from",
      "literature-informed parameter ranges; a full covariance structure was not modelled and",
      "the sampling ranges themselves are never reported, so this model file carries the",
      "nominal baseline typical values only (no IIV -- see the ini() note). The ten source",
      "studies that informed the parameter ranges are listed in Suksai 2026 Table 1 (Triggs",
      "1996, Troconiz 1998, Contin 2001, Chan 2005, Westin 2011, Mao 2013, Othman & Dutta 2014,",
      "Othman 2015, Senek 2018, Senek 2020) and the pharmacodynamic values they contributed are",
      "in Table 2; none of those studies is the source of the specific baseline numbers below,",
      "which Suksai 2026 reports as its own nominal configuration.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # ALL PARAMETERS ARE fixed(). Suksai 2026 estimates nothing: there is no
    # observed data, no objective-function fit and no reported uncertainty on
    # any value. Methods, "Model development": "a baseline parameter
    # configuration was defined from literature-informed values and used as the
    # nominal reference setting for computational experiments." Every value
    # below is that nominal configuration, quoted verbatim from the sentence
    # beginning "The baseline PK-PD parameters were ka = 1.2 h-1, ...".
    #
    # No IIV and no residual error are encoded, because the paper reports
    # neither. The 100-profile virtual cohort was drawn from "literature-informed
    # parameter ranges" that the paper never states (Methods, "Population
    # variability and virtual cohort generation"), so any variance here would be
    # invented. See the vignette Errata.
    # ------------------------------------------------------------------

    # -- Absorption ----------------------------------------------------
    lka   <- fixed(log(1.2))   ; label("First-order absorption rate constant from the gastrointestinal depot (ka, 1/h)")             # Methods, Model development: ka = 1.2 h-1

    # -- Distribution micro-constants ----------------------------------
    # Suksai 2026 parameterises distribution directly in first-order rate
    # constants, not in an inter-compartmental clearance. k12 / k21 exchange
    # with the peripheral tissue compartment; k13 / k31 exchange with the brain
    # compartment. Both pairs hang off plasma IN PARALLEL (see model()).
    lk12  <- fixed(log(0.6))   ; label("Plasma-to-peripheral-tissue distribution rate constant (k12, 1/h)")                          # Methods, Model development: k12 = 0.6 h-1
    lk21  <- fixed(log(0.5))   ; label("Peripheral-tissue-to-plasma distribution rate constant (k21, 1/h)")                          # Methods, Model development: k21 = 0.5 h-1
    lk13  <- fixed(log(0.25))  ; label("Plasma-to-brain distribution rate constant (k13, 1/h)")                                      # Methods, Model development: k13 = 0.25 h-1
    lk31  <- fixed(log(0.20))  ; label("Brain-to-plasma distribution rate constant (k31, 1/h)")                                      # Methods, Model development: k31 = 0.20 h-1

    # -- Elimination ---------------------------------------------------
    # The paper's symbol is `ke`. Equation 4 DEFINES ke as "the effect-site
    # equilibration rate constant", but the executed model uses it as the
    # first-order elimination rate constant from plasma and uses the separately
    # listed ke0 for effect-site equilibration. Three independent lines of
    # evidence, all reproduced in the vignette:
    #   (i) the baseline list contains BOTH ke and ke0, so they cannot be the
    #       same quantity;
    #   (ii) Figure 4 gives ke a large NEGATIVE sensitivity (-0.71) on peak Ce,
    #       which an equilibration constant cannot have but an elimination
    #       constant must;
    #   (iii) simulating with this assignment reproduces Figure 3 and every bar
    #       of Figure 4 (see the vignette's sensitivity-reproduction table).
    # Encoded under the canonical name for the role it actually plays.
    lkel  <- fixed(log(0.35))  ; label("First-order elimination rate constant from plasma (paper symbol ke, 1/h)")                   # Methods, Model development: ke = 0.35 h-1

    # -- Effect-site equilibration -------------------------------------
    lke0  <- fixed(log(0.6))   ; label("Effect-site equilibration rate constant (ke0, 1/h)")                                         # Methods, Model development: ke0 = 0.6 h-1

    # -- Volumes -------------------------------------------------------
    # Only vc is load-bearing: it converts the plasma amount into the
    # concentration that drives the effect site, so peak Ce is exactly
    # proportional to 1 / vc. Figure 4's Vp bar overshoots its -1 gridline,
    # and the +/-20% central-difference index of an exact inverse
    # proportionality is the analytic constant (1/1.2 - 1/0.8)/0.4 =
    # -1.0416667 -- which the implementation reproduces to six decimals.
    # Under Equation 1's clearance parameterisation Vplasma would ALSO govern
    # the distribution flux and that identity would not hold, so this is a
    # falsifier of the printed form and not merely a corroboration.
    # vp and v_brain are used ONLY to render the peripheral and brain
    # concentration curves of Figure 3's top panel; consistently, Figure 4
    # shows a sensitivity of exactly zero for both Vper and Vb.
    lvc      <- fixed(log(20)) ; label("Central (plasma) volume of distribution (paper symbol Vp, L)")                               # Methods, Model development: Vp = 20 L
    lvp      <- fixed(log(25)) ; label("Peripheral tissue volume of distribution (paper symbol Vper, L)")                            # Methods, Model development: Vper = 25 L
    lv_brain <- fixed(log(10)) ; label("Brain compartment volume of distribution (paper symbol Vb, L)")                              # Methods, Model development: Vb = 10 L

    # -- Pharmacodynamics (Equation 5) ---------------------------------
    # E = E0 + Emax * Ce^gamma / (C50^gamma + Ce^gamma). Equation 5 writes the
    # potency term as C50 while the baseline configuration and the Figure 4
    # axis both call it EC50; they are the same quantity and the canonical
    # name is used. E0 and Emax are on the paper's dimensionless "TRS-like"
    # motor-response scale (Figure 3, bottom panel y-axis: "Effect (a.u. /
    # TRS-like)"), so Emax = 10 is an ADDITIVE increment above E0 = 30, not a
    # fraction of baseline. This is settled by the figure: the bottom panel
    # spans 30 at baseline to ~36.3 at peak, i.e. a maximum span of 10 units
    # and not of 30 * 10.
    le0   <- fixed(log(30))    ; label("Baseline motor status (E0, TRS-like units)")                                                 # Methods, Model development: E0 = 30
    lemax <- fixed(log(10))    ; label("Maximum achievable drug effect above baseline (Emax, TRS-like units)")                       # Methods, Model development: Emax = 10
    lec50 <- fixed(log(2.5))   ; label("Effect-site concentration producing 50 percent of Emax (paper symbols EC50 and C50, mg/L)")  # Methods, Model development: EC50 = 2.5 mg/L
    lhill <- fixed(log(2.0))   ; label("Hill coefficient of the sigmoidal Emax response (paper symbol gamma, unitless)")             # Methods, Model development: gamma = 2.0
  })

  model({
    # 1. Individual parameters. No IIV terms: Suksai 2026 reports none
    #    (see the ini() note), so every subject solves the typical value.
    ka      <- exp(lka)
    k12     <- exp(lk12)
    k21     <- exp(lk21)
    k13     <- exp(lk13)
    k31     <- exp(lk31)
    kel     <- exp(lkel)
    ke0     <- exp(lke0)
    vc      <- exp(lvc)
    vp      <- exp(lvp)
    v_brain <- exp(lv_brain)

    e0      <- exp(le0)
    emax    <- exp(lemax)
    ec50    <- exp(lec50)
    hill    <- exp(lhill)

    # 2. ODE system.
    #
    #    depot -(ka)-> central; central <-(k12/k21)-> peripheral1 and
    #    central <-(k13/k31)-> brain IN PARALLEL; elimination kel from central.
    #
    #    The peripheral and brain compartments are two independent distribution
    #    arms off plasma. They are NOT in series: Suksai 2026's printed
    #    Equation 2 puts the effect site downstream of the peripheral
    #    compartment, but the executed model (Figures 3-4) does not, and the
    #    printed form additionally double-counts the peripheral efflux flux
    #    (subtracting it from peripheral in Equation 1 while adding it to the
    #    effect site in Equation 2, delivering the same drug to two places at
    #    once). The parallel topology is what the k13 / k31 pair -- which
    #    appears nowhere in Equations 1-4 -- requires, and it is what reproduces
    #    the published figures. Full justification in the vignette Errata.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (k12 + k13 + kel) * central +
                          k21 * peripheral1 +
                          k31 * brain
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(brain)       <-  k13 * central - k31 * brain

    # 3. Concentrations.
    Cc     <- central     / vc
    Cper   <- peripheral1 / vp
    Cbrain <- brain       / v_brain

    # 4. Effect-site equilibration.
    #
    #    The effect state holds a CONCENTRATION (mg/L), matching Suksai 2026
    #    Equation 4, which is written directly in Ce with no effect-site amount
    #    and no dividing volume.
    #
    #    Ce equilibrates with the PLASMA concentration Cc. Equation 4 as printed
    #    drives it from Aeffect / Veffect instead, but the executed model does
    #    not: Figure 4 gives Vb a sensitivity of exactly zero (had Ce chased the
    #    brain concentration, a +/-20% change in Vb would have moved peak Ce by
    #    ~17%), and Figure 3 shows Ce peaking at ~3.3 mg/L while the brain curve
    #    reaches ~7 mg/L, so Ce demonstrably does not track the brain curve.
    d/dt(effect) <- ke0 * (Cc - effect)
    Ce           <- effect

    # 5. Pharmacodynamic response (Equation 5), on the paper's TRS-like scale.
    E <- e0 + emax * Ce^hill / (ec50^hill + Ce^hill)

    # 6. No residual error model. Suksai 2026 has no observed data and reports
    #    no residual-error structure of any kind, so none is invented here; the
    #    model is a deterministic simulation framework. Cc, Cper, Cbrain, Ce and
    #    E are all returned by rxode2 as simulated outputs.
  })
}
