Hong_2013_glucose_insulin_HGC <- function() {
  description <- paste(
    "Integrated glucose-insulin population pharmacodynamic model for the",
    "hyperglycemic glucose clamp (HGC) test in adults with type 2 diabetes",
    "mellitus (Hong 2013). Built on the Silber 2007 glucose-insulin framework",
    "with two modifications: (1) endogenous glucose production GP is constant",
    "and derived at steady state from baseline glucose (GCss) and insulin",
    "(ICss) so that feedback control of GP is suppressed in DIS_DIAB (Silber's",
    "GP-feedback estimate was close to zero in DIS_DIAB); (2) the first-phase",
    "insulin response is captured by an empirical Gaussian secretion pulse",
    "(amplitude Amplitude, peak time Tsec fixed at 3.54 min, width Tdur fixed",
    "at 1.76 min) rather than a NONMEM bolus. The second-phase insulin",
    "response rises linearly with time through a constant gamma multiplied by",
    "the elevation of glucose above its steady state (rectified at zero).",
    "Glucose is eliminated via insulin-independent clearance (CLG) and",
    "insulin-dependent clearance (CLGI_HGC) gated by an effect-compartment",
    "insulin concentration ICE (rate constant kIE). Insulin follows",
    "first-order elimination (CLI / VI). Palosuran 125 mg b.i.d. was the",
    "investigational drug; the paper concludes it has no clinically",
    "meaningful effect on insulin secretion or sensitivity, so the published",
    "final estimates set the palosuran treatment effect to zero -- the",
    "structural model below is the drug-free glucose-insulin homeostasis",
    "model in DIS_DIAB. VI is fixed at 6.09 L (Silber 2007 DIS_DIAB literature",
    "value) because the within-study estimate (0.52 L) was not",
    "physiologically meaningful and biased CLI."
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
    dosing        = "mg of glucose (intravenous loading dose plus Biostator-controlled GIR into the glucose compartment)",
    concentration = "mg/L for glucose (G / VG) and mU/L for insulin (I / VI); convert to mg/dL via /10 and to mmol/L via /18.02 for glucose"
  )

  covariateData <- list(
    FPG = list(
      description        = "Baseline fasting plasma glucose concentration (GCss in the paper notation). Used to derive the constant endogenous glucose production rate GP at steady state via GP = GCss * (CLG + CLGI_HGC * ICss). Time-fixed per subject.",
      units              = "mg/L (paper reports baseline glucose in mg/dL with range 110-180 mg/dL across the DIS_DIAB cohort; multiply by 10 to convert to the mg/L scale used internally by the ODEs)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 1300 mg/L (= 130 mg/dL, mid-range for the Hong 2013 cohort of DIS_DIAB patients with fasting plasma glucose 110-180 mg/dL). Supply per-subject baseline glucose via this column to anchor each subject's drug-free steady state.",
      source_name        = "GCss"
    ),
    INS_BL = list(
      description        = "Baseline fasting plasma insulin concentration (ICss in the paper notation). Used to derive the constant endogenous glucose production rate GP and the baseline insulin secretion ICss*CLI at steady state. Time-fixed per subject.",
      units              = "mU/L (paper internal units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default reference value 13 mU/L (representative DIS_DIAB fasting insulin; paper does not report a single typical value because subject-level baseline insulin was used as a per-subject covariate). Companion canonical to FPG. NA_NA_paracetamol's INS_BL is in pmol/L with a 1/6.945 rescale; Hong 2013 uses mU/L directly so no rescaling is applied here.",
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
    dose_range     = "Hyperglycemic clamp procedure: intravenous glucose loading dose 150 mg/kg followed by a Biostator-controlled 20% glucose infusion regulated to clamp blood glucose at 240 mg/dL for 120 min. Palosuran 125 mg b.i.d. (or placebo) was administered orally twice daily for 4 weeks and the HGC was performed 1 h after drug administration on day 29 of each treatment period; palosuran had no clinically meaningful effect and is set to zero in the published final model.",
    regions        = "Germany (Ethikkommission der Aerztekammer Nordrhein, Germany).",
    notes          = "Subject demographics from Hong 2013 Methods 'Patients and Study Design'. Two-way crossover design (palosuran vs placebo) with 4-week treatment periods and 4-week washout; MTT performed on day 28 and HGC on day 29 of each period. Subject-level palosuran treatment effects (theta_pal in Eqs 8-9) were estimated as -0.122 with 95% CI including zero and -7.44% on glucose elimination (not clinically relevant), so the published final estimates fix the palosuran effect to zero. Glucose volume of distribution VI was fixed to 6.09 L from the Silber 2007 DIS_DIAB literature value because the within-study estimate (0.52 L) was not physiologically meaningful."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters (Hong 2013 Table I, Original-data Estimate
    # column). VI is fixed to the Silber 2007 DIS_DIAB literature value
    # (6.09 L) because the within-study estimate was not physiologically
    # meaningful; Tsec and Tdur are fixed to literature values from
    # reference 16 (Lima 2010 / Mari 2003 first-phase secretion-pulse
    # parameters) because the within-study estimates (3.27 and 0.98 min)
    # were too narrow relative to the 120-min HGC time-span and
    # destabilised the variance-covariance matrix.
    # ---------------------------------------------------------------------
    lclg       <- log(0.164)            ; label("Insulin-independent glucose clearance CLG (L/min)")                          # Hong 2013 Table I: CL_G = 0.164
    lclgi_hgc  <- log(0.0111)           ; label("Insulin-dependent glucose clearance for HGC CLGI_HGC (L/min/(mU/L))")        # Hong 2013 Table I: CL_GI_HGC = 0.0111
    lvg        <- log(23.7)             ; label("Glucose volume of distribution VG (L)")                                       # Hong 2013 Table I: V_G = 23.7
    lgamma     <- log(0.000431)         ; label("Second-phase insulin-secretion slope gamma (mU/(min^2*(mg/L)))")             # Hong 2013 Table I: gamma = 0.000431
    lamplitude <- log(32.2)             ; label("Amplitude of the Gaussian first-phase insulin pulse (mU)")                    # Hong 2013 Table I: Amplitude = 32.2
    lcli       <- log(1.54)             ; label("Insulin clearance CLI (L/min)")                                              # Hong 2013 Table I: CL_I = 1.54
    lvi        <- fixed(log(6.09))      ; label("Insulin volume of distribution VI (L; FIXED to Silber 2007 DIS_DIAB)")           # Hong 2013 Table I: V_I = 6.09 (FIXED, footnote a)
    lkie       <- log(0.00291)          ; label("Insulin effect-compartment equilibration rate kIE (1/min)")                  # Hong 2013 Table I: k_IE = 0.00291
    ltsec      <- fixed(log(3.54))      ; label("Time of peak first-phase insulin secretion Tsec (min; FIXED to lit value)")  # Hong 2013 Results: Tsec FIXED at 3.54 (ref 16)
    ltdur      <- fixed(log(1.76))      ; label("Width of the Gaussian first-phase pulse Tdur (min; FIXED to lit value)")     # Hong 2013 Results: Tdur FIXED at 1.76 (ref 16)

    # ---------------------------------------------------------------------
    # Inter-individual variability (Table I "Random effects model IIV"
    # column). Reported as %CV; converted to log-scale variance with
    # omega^2 = log(CV^2 + 1). Hong 2013 fixed IIV to zero on Tsec, Tdur,
    # gamma, and CLGI_HGC.
    # ---------------------------------------------------------------------
    etalclg       ~ log(0.0957^2 + 1)             # CLG IIV 9.57% CV  -> var = log(1.00916) = 0.00912
    etalvg        ~ log(0.158^2  + 1)             # VG  IIV 15.8% CV  -> var = log(1.02496) = 0.0247
    etalamplitude ~ log(1.97^2   + 1)             # Amplitude IIV 197% CV -> var = log(4.8809) = 1.586
    etalcli       ~ log(0.857^2  + 1)             # CLI IIV 85.7% CV  -> var = log(1.7344)  = 0.551
    etalvi        ~ log(0.776^2  + 1)             # VI  IIV 77.6% CV  -> var = log(1.6022)  = 0.471
    etalkie       ~ log(0.902^2  + 1)             # kIE IIV 90.2% CV  -> var = log(1.8136)  = 0.595

    # ---------------------------------------------------------------------
    # Inter-occasion variability NOT structurally encoded. Hong 2013
    # Table I reports IOV on CLG (31.8% CV) and VG (9.22% CV) across the
    # two HGC occasions (one per treatment period); these are NOT
    # encoded here because nlmixr2lib follows the Andrews 2017 /
    # Brooks 2021 tacrolimus precedent: the rxode2 mu-reference parser
    # does not accept `theta + eta_iiv + eta_iov` on a single line, and
    # the model-library use case has no operational occasion column.
    # Downstream users who want to simulate IOV can add an OCC indicator
    # and a per-occasion eta in their own rxode2 model. See vignette
    # Assumptions and deviations.
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Residual error (Table I "Residual proportional error" column).
    # Hong 2013 fits an additive error on log-transformed observations,
    # which the paper notes "approximately corresponds to a proportional
    # error model on nontransformed data". We implement it as proportional
    # in linear space, matching the paper's interpretation. Separate
    # parameters per output.
    # ---------------------------------------------------------------------
    propSd_Gc <- 0.103                  ; label("Proportional residual SD on HGC glucose concentration (fraction)")  # Hong 2013 Table I: sigma_G_HGC = 10.3%
    propSd_Ic <- 0.257                  ; label("Proportional residual SD on HGC insulin concentration (fraction)")  # Hong 2013 Table I: sigma_I_HGC = 25.7%
  })

  model({
    # Individual structural parameters. IIV-bearing parameters use the
    # standard mu-reference form `exp(lX + etalX)`; IOV (reported by
    # paper but not encoded; see ini() comment) would add an additional
    # per-occasion eta term that the parser does not support on the
    # same line.
    clg       <- exp(lclg       + etalclg)
    clgi_hgc  <- exp(lclgi_hgc)
    vg        <- exp(lvg        + etalvg)
    gamma     <- exp(lgamma)
    amplitude <- exp(lamplitude + etalamplitude)
    cli       <- exp(lcli       + etalcli)
    vi        <- exp(lvi        + etalvi)
    kie       <- exp(lkie       + etalkie)
    tsec      <- exp(ltsec)
    tdur      <- exp(ltdur)

    # Per-subject baseline glucose (mg/L) and insulin (mU/L) from the
    # covariate columns FPG (= GCss) and INS_BL (= ICss).
    gcss <- FPG
    icss <- INS_BL

    # Endogenous glucose production (constant, derived at steady state
    # from baseline glucose and insulin per paper Methods "The
    # endogenous glucose production was expressed as a function of
    # glucose concentrations at steady-state (GCss) and the elimination
    # rate of glucose at baseline"). At baseline, dG/dt = 0 and ICE =
    # ICss, so GP balances total elimination: GP = GCss*CLG +
    # CLGI_HGC*GCss*ICss (mg/min). The factor VG cancels because
    # CLG/VG*G_ss = CLG*GCss when G_ss = GCss*VG.
    gp <- gcss * (clg + clgi_hgc * icss)

    # Insulin-secretion components (paper Methods, "Model for
    # Hyperglycemia Glucose Clamp Test"). Baseline insulin secretion is
    # the steady-state product ICss*CLI. The first-phase Gaussian pulse
    # uses literature-fixed Tsec / Tdur; small epsilon avoids division
    # by zero when Tdur is at its lower bound.
    ir_base  <- icss * cli
    ir_first <- amplitude / (tdur * sqrt(2 * 3.14159265358979) + 1e-12) *
                exp(-(t - tsec)^2 / (2 * tdur^2 + 1e-12))

    # Second-phase rate: gamma * t * elevation_above_baseline,
    # rectified at zero so the term vanishes whenever G/VG <= GCss
    # (matches the paper's "set to zero when the glucose concentration
    # was equal to GCss"). Glucose concentration is G/VG (mg/L).
    glucose_conc <- glucose / vg
    delta_g      <- glucose_conc - gcss
    ir_second    <- gamma * t * delta_g * (delta_g > 0)

    # ODE system (paper Eqs 1-3 with constant GP). Glucose is
    # eliminated insulin-independent (CLG / VG) and insulin-dependent
    # (CLGI_HGC / VG * ICE); ICE is the effect-compartment insulin
    # concentration in mU/L. The Biostator-controlled glucose infusion
    # rate (GIR) enters the glucose compartment as a dose / infusion via
    # the dataset (amt + rate columns on cmt = glucose); no infusion
    # term is hard-coded in the ODE so downstream users can supply any
    # GIR profile.
    d/dt(glucose) <- gp - (clg / vg) * glucose - (clgi_hgc / vg) * glucose * effect
    d/dt(insulin) <- ir_base + ir_first + ir_second - (cli / vi) * insulin
    d/dt(effect)  <- kie * (insulin / vi - effect)

    # Initial conditions at steady state. Glucose amount G_ss = GCss *
    # VG; insulin amount I_ss = ICss * VI; effect-compartment insulin
    # concentration ICE_ss = ICss.
    glucose(0) <- gcss * vg
    insulin(0) <- icss * vi
    effect(0)  <- icss

    # Observation variables: glucose concentration (mg/L) and insulin
    # concentration (mU/L) in plasma. Convert glucose to mg/dL by
    # dividing by 10 in downstream analyses if desired.
    Gc <- glucose_conc
    Ic <- insulin / vi

    Gc ~ prop(propSd_Gc)
    Ic ~ prop(propSd_Ic)
  })
}
