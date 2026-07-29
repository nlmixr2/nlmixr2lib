Hong_2013_glucose_insulin_MTT <- function() {
  description <- paste(
    "Integrated glucose-insulin population pharmacodynamic model for the",
    "meal tolerance test (MTT) in adults with type 2 diabetes mellitus",
    "(Hong 2013). Companion model to Hong_2013_glucose_insulin_HGC: shares",
    "the glucose-disposition and effect-compartment structure but replaces",
    "the intravenous glucose-clamp input with an oral-meal absorption arm",
    "and replaces the biphasic insulin-secretion equation with a",
    "power-function plus incretin Emax form. Glucose absorption from the",
    "meal is modelled as a continuous Savic 2007 analytical transit",
    "compartment chain (rxode2 transit(n, mtt, BIO)) with non-integer n =",
    "0.781, mean transit time MTT = 62.5 min, and a 'bioavailability' BIO",
    "(0.252) multiplied by a dummy meal dose of 100 g of glucose. The",
    "incretin effect (insulin secretion triggered by oral glucose intake",
    "via GLP-1) is captured by an Emax function of the instantaneous",
    "absorption rate ABSG: Emax = 2.02, ABSG50 = 14.8 mg/min FIXED to the",
    "Jauslin 2007 literature value because Emax / ABSG50 were correlated.",
    "Insulin secretion has no first-phase pulse (in contrast to the HGC",
    "model); the second-phase response is a power function of",
    "(glucose_concentration / glucose_baseline) with exponent IPRG = 3.06,",
    "modulated by the incretin factor. The disposition parameters CLG, VG,",
    "CLI, VI, kIE are FIXED at the HGC point estimates per the paper's",
    "Results section, on the assumption that after accounting for the",
    "incretin effect the disposition parameters should not markedly differ",
    "between HGC and MTT. Palosuran 125 mg b.i.d. was the investigational",
    "drug and was found to have no clinically meaningful effect, so the",
    "published final estimates set the palosuran treatment effect to zero",
    "-- the structural model below is the drug-free glucose-insulin",
    "homeostasis model in DIS_DIAB under MTT conditions."
  )
  reference <- paste(
    "Hong Y, Dingemanse J, Sidharta P, Mager DE (2013).",
    "Population Pharmacodynamic Modeling of Hyperglycemic Clamp and Meal",
    "Tolerance Tests in Patients with Type 2 Diabetes Mellitus.",
    "The AAPS Journal 15(4):1051-1063.",
    "doi:10.1208/s12248-013-9512-4.",
    "PMID 23913136; PMCID PMC3787234.",
    sep = " "
  )
  vignette <- "Hong_2013_glucose_insulin"

  units <- list(
    time          = "min",
    dosing        = "mg of glucose (dummy meal dose of 100000 mg = 100 g into the glucose compartment; absorbed via the Savic transit chain)",
    concentration = "mg/L for glucose (G / VG) and mU/L for insulin (I / VI); convert glucose to mg/dL via /10 and to mmol/L via /18.02"
  )

  covariateData <- list(
    FPG = list(
      description        = "Baseline fasting plasma glucose concentration (GCss in the paper notation). Used to derive the constant endogenous glucose production rate GP at steady state and to set the reference glucose for the power-function second-phase insulin secretion (G/(VG*GCss))^IPRG. Time-fixed per subject.",
      units              = "mg/L (paper reports baseline glucose in mg/dL with range 110-180 mg/dL across the DIS_DIAB cohort; multiply by 10 to convert to the mg/L scale used internally by the ODEs)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 1300 mg/L (= 130 mg/dL, mid-range for the Hong 2013 DIS_DIAB cohort with fasting plasma glucose 110-180 mg/dL).",
      source_name        = "GCss"
    ),
    INS_BL = list(
      description        = "Baseline fasting plasma insulin concentration (ICss in the paper notation). Used to derive the constant endogenous glucose production rate GP and the baseline insulin secretion ICss*CLI at steady state. Time-fixed per subject.",
      units              = "mU/L (paper internal units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 13 mU/L (representative DIS_DIAB fasting insulin). Companion canonical to FPG.",
      source_name        = "ICss"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "40-65 years (mean 53.7, SD 7.3)",
    weight_range   = NA_character_,
    sex_female_pct = 20,
    disease_state  = "Type 2 diabetes mellitus (DIS_DIAB) managed by diet only. Fasting blood glucose 110-180 mg/dL; HbA1c 5.4-8.3% (mean 6.4%, SD 0.8%). Patients excluded if treated with an antidiabetic drug in the 2 months prior to screening or with severe diabetes complications.",
    dose_range     = "Meal tolerance test: standardised breakfast approximately 618 kcal (65% carbohydrates, 17% proteins, 18% lipids) given after an overnight fast; modelled as a dummy oral dose of 100 g glucose into the gut absorbed via the Savic transit chain (n / MTT / BIO). Palosuran 125 mg b.i.d. (or placebo) was administered orally twice daily for 4 weeks and the MTT was performed 1 h after drug administration on day 28 of each treatment period; palosuran had no clinically meaningful effect and is set to zero in the published final model.",
    regions        = "Germany (Ethikkommission der Aerztekammer Nordrhein, Germany).",
    notes          = "Subject demographics from Hong 2013 Methods 'Patients and Study Design'. Two-way crossover design (palosuran vs placebo) with 4-week treatment periods and 4-week washout; MTT performed on day 28 and HGC on day 29 of each period. Published final estimates re-fit the MTT after fixing palosuran effect to zero; CLG / VG / CLI / VI / kIE are FIXED to the HGC point estimates from Table I. The MTT IOV table reflects intra-individual variability across the four MTT occasions (two periods x palosuran/placebo). ABSG50 was fixed to 14.8 mg/min (Jauslin 2007 lit value) because Emax / ABSG50 were correlated; sensitivity analysis fixing ABSG50 to 0.148 mg/min did not materially change estimates except for n (1.42 vs 0.781) and MTT (44.2 vs 62.5)."
  )

  ini({
    # ---------------------------------------------------------------------
    # MTT-specific structural parameters (Hong 2013 Table II Original-data
    # Estimate column). ABSG50 is fixed to the Jauslin 2007 literature
    # value because Emax / ABSG50 were correlated in the within-study fit.
    # ---------------------------------------------------------------------
    ln          <- log(0.781)            ; label("Number of absorption transit compartments n (continuous; dimensionless)")           # Hong 2013 Table II: N = 0.781
    lmtt        <- log(62.5)             ; label("Mean transit time MTT (min)")                                                       # Hong 2013 Table II: Mean transit time = 62.5
    lbio        <- log(0.252)            ; label("Bioavailability factor BIO on the dummy 100 g meal-glucose dose (fraction)")        # Hong 2013 Table II: BIO = 0.252
    lclgi_mtt   <- log(0.00425)          ; label("Insulin-dependent glucose clearance for MTT CLGI_MTT (L/min/(mU/L))")               # Hong 2013 Table II: CL_GI_MTT = 0.00425
    lemax       <- log(2.02)             ; label("Maximal incretin-driven insulin-secretion stimulation Emax (unitless)")             # Hong 2013 Table II: E_max = 2.02
    labsg50     <- fixed(log(14.8))      ; label("Glucose absorption rate giving 50 percent of maximal incretin effect ABSG50 (mg/min; FIXED to Jauslin 2007)")  # Hong 2013 Table II: ABSG_50 = 14.8 (FIXED, footnote a)
    liprg       <- log(3.06)             ; label("Power coefficient on (G/(VG*GCss)) for second-phase insulin secretion IPRG (unitless)")  # Hong 2013 Table II: IPRG = 3.06

    # ---------------------------------------------------------------------
    # Disposition parameters fixed at HGC point estimates (paper Results:
    # "The population parameters including CLG, VG, CLI, VI, and kIE were
    # fixed to those estimated from the HGC model with the assumption
    # that these disposition parameters should not markedly differ
    # between HGC and MTT after the incretin effect had been taken into
    # account").
    # ---------------------------------------------------------------------
    lclg        <- fixed(log(0.164))     ; label("Insulin-independent glucose clearance CLG (L/min; FIXED from HGC Table I)")        # Hong 2013 Table I (HGC) -> MTT FIXED
    lvg         <- fixed(log(23.7))      ; label("Glucose volume of distribution VG (L; FIXED from HGC Table I)")                    # Hong 2013 Table I (HGC) -> MTT FIXED
    lcli        <- fixed(log(1.54))      ; label("Insulin clearance CLI (L/min; FIXED from HGC Table I)")                            # Hong 2013 Table I (HGC) -> MTT FIXED
    lvi         <- fixed(log(6.09))      ; label("Insulin volume of distribution VI (L; FIXED from HGC Table I)")                    # Hong 2013 Table I (HGC) -> MTT FIXED
    lkie        <- fixed(log(0.00291))   ; label("Insulin effect-compartment equilibration rate kIE (1/min; FIXED from HGC Table I)") # Hong 2013 Table I (HGC) -> MTT FIXED

    # ---------------------------------------------------------------------
    # Inter-individual variability (Table II "Random effects model IIV"
    # column). IIV fixed to zero on n and ABSG50 per the paper.
    # ---------------------------------------------------------------------
    etalmtt      ~ log(0.306^2 + 1)             # MTT     IIV 30.6% CV -> var = log(1.09364) = 0.0895
    etalbio      ~ log(0.415^2 + 1)             # BIO     IIV 41.5% CV -> var = log(1.17223) = 0.159
    etalclgi_mtt ~ log(0.476^2 + 1)             # CLGI_MTT IIV 47.6% CV -> var = log(1.22658) = 0.204
    etalemax     ~ log(0.568^2 + 1)             # Emax    IIV 56.8% CV -> var = log(1.32262) = 0.280
    etaliprg     ~ log(0.389^2 + 1)             # IPRG    IIV 38.9% CV -> var = log(1.15132) = 0.141

    # ---------------------------------------------------------------------
    # Inter-occasion variability NOT structurally encoded. Hong 2013
    # Table II reports IOV on CLG (65.7% CV) and VG (19.5% CV) across
    # the four MTT occasions (two crossover periods x
    # palosuran/placebo); these are NOT encoded here for the same
    # reason as the HGC companion model -- rxode2's mu-reference parser
    # cannot combine IIV and IOV etas on a single mu-referenced line,
    # and the model-library use case has no operational occasion
    # column. See vignette Assumptions and deviations. Note that the
    # MTT-table IOV values are larger than the HGC-table values
    # because they span four occasions and reflect intra-individual
    # variability across the placebo / palosuran arms.
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Residual error (Table II "Residual proportional error" column).
    # As with the HGC model, the paper's "log-additive" form is
    # implemented as proportional in linear space per the paper's own
    # interpretation.
    # ---------------------------------------------------------------------
    propSd_Gc <- 0.0702                ; label("Proportional residual SD on MTT glucose concentration (fraction)")  # Hong 2013 Table II: sigma_G_MTT = 7.02%
    propSd_Ic <- 0.302                 ; label("Proportional residual SD on MTT insulin concentration (fraction)")  # Hong 2013 Table II: sigma_I_MTT = 30.2%
  })

  model({
    # Individual structural parameters.
    nn        <- exp(ln)
    mtt       <- exp(lmtt       + etalmtt)
    bio       <- exp(lbio       + etalbio)
    clgi_mtt  <- exp(lclgi_mtt  + etalclgi_mtt)
    emax      <- exp(lemax      + etalemax)
    absg50    <- exp(labsg50)
    iprg      <- exp(liprg      + etaliprg)

    # Disposition fixed-from-HGC. IOV on CLG / VG reported in Table II
    # but not encoded structurally; see ini() comment above.
    clg       <- exp(lclg)
    vg        <- exp(lvg)
    cli       <- exp(lcli)
    vi        <- exp(lvi)
    kie       <- exp(lkie)

    # Per-subject baseline glucose (mg/L) and insulin (mU/L) from the
    # covariate columns FPG (= GCss) and INS_BL (= ICss).
    gcss <- FPG
    icss <- INS_BL

    # Endogenous glucose production (constant; see HGC companion for the
    # steady-state derivation). The MTT uses CLGI_MTT for the
    # insulin-dependent arm of the baseline balance because all
    # baseline-rate quantities in the MTT model use the MTT-specific
    # disposition assumption (CLGI_MTT replaces CLGI_HGC).
    gp <- gcss * (clg + clgi_mtt * icss)

    # Glucose absorption rate from the Savic 2007 analytical transit
    # compartment chain (rxode2 builtin `transit(n, mtt, bio)` computes
    # the per-time absorption rate from the prior-dose amount and time
    # since dose). The dummy 100 g (= 100000 mg) meal dose is supplied
    # via the dataset as a dose into the glucose compartment with f
    # (glucose) = 0 below so the bolus is suppressed and the transit()
    # chain drives the absorption rate. The paper Eq 5 form is exactly
    # the Savic analytical Erlang form with continuous n.
    absg <- transit(nn, mtt, bio)

    # Insulin-secretion rate (paper Eq 6-7). Power function in
    # (glucose_concentration / GCss) for the second-phase response,
    # modulated by the Emax incretin factor in ABSG. Small epsilons
    # prevent zero divisions when glucose_concentration = 0 or
    # ABSG = 0 at t = 0.
    glucose_conc <- glucose / vg
    glu_ratio    <- glucose_conc / (gcss + 1e-12)
    incretin     <- 1 + emax * absg / (absg50 + absg + 1e-12)
    ir_base      <- icss * cli
    ir           <- ir_base * (glu_ratio^iprg) * incretin

    # ODE system (paper Eqs 4 and 6, plus shared insulin-effect Eq 3).
    # Glucose dynamics: constant production, transit-absorbed input,
    # insulin-independent and insulin-dependent elimination. Insulin
    # dynamics: secretion (power x incretin) minus first-order
    # elimination. Effect compartment as in HGC.
    d/dt(glucose) <- gp + absg - (clg / vg) * glucose - (clgi_mtt / vg) * glucose * effect
    d/dt(insulin) <- ir - (cli / vi) * insulin
    d/dt(effect)  <- kie * (insulin / vi - effect)

    # Initial conditions at steady state.
    glucose(0) <- gcss * vg
    insulin(0) <- icss * vi
    effect(0)  <- icss

    # Suppress the dose bolus on the glucose compartment so the
    # transit() chain drives the absorption rate from the prior-dose
    # amount and time. Downstream users supply the meal dose via the
    # dataset (e.g., et(amt = 100000, cmt = "glucose")).
    f(glucose) <- 0

    # Observation variables.
    Gc <- glucose_conc
    Ic <- insulin / vi

    Gc ~ prop(propSd_Gc)
    Ic ~ prop(propSd_Ic)
  })
}
