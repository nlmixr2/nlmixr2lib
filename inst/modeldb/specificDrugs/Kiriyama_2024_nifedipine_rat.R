Kiriyama_2024_nifedipine_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment IV-infusion PK of nifedipine in spontaneously hypertensive rats (SHR/Izm)",
    "linked to four effect-compartment pharmacodynamic endpoints:",
    "mean blood pressure (BP) via a hypotensive Emax arm plus a homeostatic",
    "endogenous-hypertensive-substance (EHS) counter-regulatory arm,",
    "and heart rate (HR), QT interval (QT) and RR-corrected QT (QTc) via ordinary sigmoid Emax models.",
    "Each endpoint carries its own effect-compartment equilibration rate ke0.",
    "The BP model's hypotensive and hypertensive arms are parameterised with Emax as the",
    "ASYMPTOTIC BP LEVEL rather than as a change from baseline (see the ini() notes).",
    "Fit by naive pooling to mean profiles in WinNonlin; the paper reports no IIV and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Kiriyama A, Kimura S, Yamashita S.",
    "Exploring the multiple effects of nifedipine and captopril administration in",
    "spontaneously hypertensive rats through pharmacokinetic-pharmacodynamic analyses.",
    "Pharmacol Res Perspect. 2024;12(4):e1249.",
    "doi:10.1002/prp2.1249.",
    "Sister model file from the same paper:",
    "modellib('Kiriyama_2024_captopril_rat').",
    "The BP homeostasis (EHS feedback) structure was first reported by the same group in",
    "Kiriyama A et al. (reference 5 of the 2024 paper), a nifedipine + propranolol SHR study;",
    "the 2024 paper restates the full equation set, so no upstream source is required.",
    sep = " "
  )
  vignette <- "Kiriyama_2024_nifedipine_captopril_rat"
  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  # No covariates: the paper fits mean profiles with absolute (not weight-normalised)
  # V1 and micro-constants, so body weight enters only when converting the reported
  # mg/kg dose into the absolute infused amount. See population$notes.
  covariateData <- list()

  compartmentData <- list(
    central     = list(analyte = "nifedipine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nifedipine", units = "ng", specimen = "plasma", verified = TRUE),
    effect_bp   = list(analyte = "nifedipine", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    ehs         = list(analyte = "endogenous hypertensive substance", units = "arbitrary unit", specimen = "not applicable", verified = TRUE),
    effect_hr   = list(analyte = "nifedipine", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_qt   = list(analyte = "nifedipine", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_qtc  = list(analyte = "nifedipine", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  paper_specific_compartments <- c(
    # Four parallel (NOT chained) effect compartments, one per PD endpoint, because
    # the paper estimates a separate ke0 for each endpoint (Table 3). The blessed
    # `effect<n>` numbered form would falsely imply a chain, so each carries its
    # endpoint suffix instead. Precedent for an `effect_<suffix>` compartment:
    # Kleijn_2011_sugammadex_rocuronium.R (`effect_roc`).
    "effect_bp", "effect_hr", "effect_qt", "effect_qtc",
    # Endogenous hypertensive substance pool driving the BP homeostatic rebound.
    # Precedent for paper-mechanistic cardiovascular states: Fu_2022_atenolol_qsp.R
    # (`hr`, `edv`, `tpr`, `ctr`).
    "ehs"
  )

  population <- list(
    species        = "rat (spontaneously hypertensive rat, SHR/Izm; male)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "not reported",
    weight_range   = "215-315 g",
    sex_female_pct = 0,
    disease_state  = "spontaneous hypertension (SHR/Izm strain) under urethane anaesthesia (1.0 g/kg intraperitoneal)",
    dose_range     = "nifedipine 1.0 mg/kg (monotherapy) or 0.5 mg/kg (coadministered with captopril), 30-minute IV infusion into the right femoral vein in 0.5 mL polyethylene glycol 400",
    regions        = "Japan (Doshisha Women's College of Liberal Arts, Kyoto)",
    biomarkers     = paste(
      "Mean arterial BP via a left-carotid-artery catheter; HR and QT via ECG;",
      "both recorded continuously from -15 to 180 min after infusion start.",
      "QTc = QT / sqrt(RR) (Bazett). Plasma nifedipine sampled from the right jugular vein",
      "5 to 300 min after infusion start and assayed by LC-MS/MS (m/z 348 -> 316).",
      sep = " "
    ),
    notes = paste(
      "Per-analysis group sizes differ: PK parameters (Tables 1-2) come from 3 to 6 rats,",
      "PD parameters (Table 3) from 8 to 10 rats; each experiment used 3-10 rats.",
      "n_subjects records the PD maximum (10).",
      "DOSE UNITS: the paper reports mg/kg doses but fits absolute V1 (mL) and absolute",
      "micro-constants, so a body weight is needed to build an event table. Back-solving the",
      "Figure 2 fitted peak concentrations against the Table 2 monotherapy parameters gives",
      "247 g from the nifedipine panel and 249 g from the captopril panel -- two independent",
      "estimates agreeing to within 1% and both inside the stated 215-315 g range. The",
      "validation vignette therefore doses a 250 g rat (1.0 mg/kg = 0.25 mg = 2.5e5 ng).",
      "The reported body weight is NOT used as a model covariate.",
      "Urethane anaesthesia inhibits CYP3A4 and the paper flags it as a potential confounder",
      "on nifedipine PK (Discussion); BP baselines under anaesthesia are lower than in awake SHRs.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETICS -- Kiriyama 2024 Table 2, nifedipine "Monotherapy"
    # column (1.0 mg/kg). The 2-compartment structure and the ODEs are
    # Methods section 2.4 / Figure 1. All four values are WinNonlin
    # least-squares estimates on the MEAN concentration-time profile
    # (AIC 185.3); the paper reports no standard errors, so nothing is
    # wrapped in fixed().
    #
    # The paper also reports a separate coadministration PK fit
    # (V1 = 26.7 mL, k10 = 0.0768, k12 = 0.132, k21 = 0.0764 /min;
    # AIC 145.0) used for its drug-interaction simulations. Those values
    # are NOT encoded as a covariate effect because the paper explicitly
    # concludes that none of the PK differences reached significance
    # ("no significant PK interactions were observed"); the vignette
    # applies them through rxSolve(params = ) instead.
    # ==================================================================
    lvc  <- log(54.9);    label("Central volume V1 (mL)")                                   # Table 2 nifedipine monotherapy V1 = 54.9 mL
    lkel <- log(0.0218);  label("Elimination rate constant k10 (1/min)")                    # Table 2 nifedipine monotherapy k10 = 0.0218 1/min
    lk12 <- log(0.0121);  label("Central-to-peripheral rate constant k12 (1/min)")          # Table 2 nifedipine monotherapy k12 = 0.0121 1/min
    lk21 <- log(0.0165);  label("Peripheral-to-central rate constant k21 (1/min)")          # Table 2 nifedipine monotherapy k21 = 0.0165 1/min

    # ==================================================================
    # BLOOD PRESSURE -- Kiriyama 2024 Table 3, nifedipine "BP" column,
    # model (A) of Figure 1 (AIC 151.4).
    #
    # Structure (Methods 2.5). An effect compartment C3 drives a
    # hypotensive Emax arm; the resulting BP deficit stimulates
    # production of an endogenous hypertensive substance (EHS, C4) whose
    # own Emax arm pushes BP back up:
    #
    #   dE/dt = -(E0 - Emax1)*EC50_1/(EC50_1 + C3)^2 * dC3/dt
    #           +(Emax2 - E0)*EC50_2/(EC50_2 + C4)^2 * dC4/dt
    #
    # with C3 = 0, C4 = 0, E = E0 at t = 0. That dE/dt is an exact
    # differential; integrating it under those initial conditions gives
    # the closed form used in model():
    #
    #   E = (E0*EC50_1 + Emax1*C3)/(EC50_1 + C3)
    #     + (E0*EC50_2 + Emax2*C4)/(EC50_2 + C4)
    #     - E0
    #
    # NOTE ON Emax PARAMETERISATION: because the printed dE/dt carries
    # the factors (E0 - Emax1) and (Emax2 - E0), Emax1 and Emax2 are the
    # ASYMPTOTIC BP LEVELS each arm drives toward, not changes from
    # baseline. Emax1 = 67.1 mmHg is therefore the floor of the
    # hypotensive arm (consistent with the measured nadir of
    # 70.9 +/- 8.6 mmHg reported in Results); reading it as a 67.1 mmHg
    # DECREASE would put the nadir near 37 mmHg, which Figure 3A refutes.
    # ==================================================================
    lke0_bp      <- log(0.135);   label("Effect-compartment equilibration rate ke0 for BP (1/min)")                        # Table 3 nifedipine BP ke0 = 0.135 1/min
    le0_bp       <- log(96.8);    label("Baseline BP before dosing, E0 (mmHg)")                                            # Table 3 nifedipine BP E0 = 96.8 mmHg
    lemax_bp     <- log(67.1);    label("Asymptotic BP level of the hypotensive arm, Emax,1 (mmHg)")                       # Table 3 nifedipine BP Emax,1 = 67.1 mmHg
    lec50_bp     <- log(86.3);    label("Effect-compartment concentration at half-maximal hypotensive effect, EC50,1 (ng/mL)") # Table 3 nifedipine BP EC50,1 = 86.3 ng/mL

    # EHS turnover pool. The paper's K_EHS,in / k_EHS,out are the
    # zero-order production and first-order degradation rate constants of
    # the EHS pool, i.e. the canonical indirect-response kin / kout roles
    # with the pool named by the `_ehs` suffix.
    lkin_ehs     <- log(0.0217);  label("Zero-order EHS production rate constant K_EHS,in (unit/min)")                     # Table 3 nifedipine BP K_EHS,in = 0.0217 unit/min
    lkout_ehs    <- log(0.0119);  label("First-order EHS degradation rate constant k_EHS,out (1/min)")                     # Table 3 nifedipine BP k_EHS,out = 0.0119 1/min

    # Feedback gain. The paper's symbol is alpha, described as "a slope of
    # a linear equation" for the homeostatic feedback system. Its role --
    # a linear amplification factor applied to a normalised deviation of a
    # state from its baseline -- is the registered `lgamma` role (the same
    # shape as the Tetschke 2018 EPO feedback factor
    # gamma * (THB_MASS - thb) / THB_MASS), so it takes the `lgamma`
    # canonical rather than `lhill`. Table 3 prints the unit as mmHg/min,
    # but the equation multiplies the dimensionless ratio (E0 - E)/E0, so
    # alpha is dimensionless; the printed unit is a typographical slip.
    lgamma_bp    <- log(1.17);    label("Linear feedback gain alpha from the fractional BP deficit onto EHS production (unitless)") # Table 3 nifedipine BP alpha = 1.17

    lemax_bp_ehs <- log(49338);   label("Asymptotic BP level of the EHS-driven hypertensive arm, Emax,2 (mmHg)")           # Table 3 nifedipine BP Emax,2 = 49 338 mmHg
    lec50_bp_ehs <- log(1450);    label("EHS level at half-maximal hypertensive effect, EC50,2 (arbitrary unit)")          # Table 3 nifedipine BP EC50,2 = 1450 unit

    # ==================================================================
    # HEART RATE -- Table 3, nifedipine "HR" column. Ordinary inhibitory
    # sigmoid Emax with an effect compartment (Figure 1B, AIC 268.9):
    #   E = E0 - Emax * C3^gamma / (EC50^gamma + C3^gamma)
    # Here Emax IS a change from baseline, so the floor is
    # 306 - 63.3 = 242.7 beat/min, matching the ~250 beat/min plateau
    # reported in Results and drawn in Figure 3B.
    # ==================================================================
    lke0_hr   <- log(0.006589); label("Effect-compartment equilibration rate ke0 for HR (1/min)")          # Table 3 nifedipine HR ke0 = 0.006589 1/min
    le0_hr    <- log(306);      label("Baseline heart rate before dosing, E0 (beat/min)")                  # Table 3 nifedipine HR E0 = 306 beat/min
    lemax_hr  <- log(63.3);     label("Maximum heart-rate reduction, Emax (beat/min)")                     # Table 3 nifedipine HR Emax = 63.3 beat/min
    lec50_hr  <- log(115);      label("Effect-compartment concentration at half-maximal HR effect, EC50 (ng/mL)") # Table 3 nifedipine HR EC50 = 115 ng/mL
    lhill_hr  <- log(0.931);    label("Hill coefficient gamma for the HR sigmoid (unitless)")              # Table 3 nifedipine HR gamma = 0.931

    # ==================================================================
    # QT INTERVAL -- Table 3, nifedipine "QT" column. Stimulatory sigmoid
    # Emax with an effect compartment (Figure 1B, AIC 191.1):
    #   E = E0 + Emax * C3^gamma / (EC50^gamma + C3^gamma)
    # ==================================================================
    lke0_qt   <- log(0.0117);   label("Effect-compartment equilibration rate ke0 for QT (1/min)")          # Table 3 nifedipine QT ke0 = 0.0117 1/min
    le0_qt    <- log(75.8);     label("Baseline QT interval before dosing, E0 (msec)")                     # Table 3 nifedipine QT E0 = 75.8 msec
    lemax_qt  <- log(76.3);     label("Maximum QT prolongation, Emax (msec)")                              # Table 3 nifedipine QT Emax = 76.3 msec
    lec50_qt  <- log(3227);     label("Effect-compartment concentration at half-maximal QT effect, EC50 (ng/mL)") # Table 3 nifedipine QT EC50 = 3227 ng/mL
    lhill_qt  <- log(0.857);    label("Hill coefficient gamma for the QT sigmoid (unitless)")              # Table 3 nifedipine QT gamma = 0.857

    # ==================================================================
    # RR-CORRECTED QT (QTc, Bazett) -- these five values are reported in
    # the Results prose (section 3, paragraph on Figure 3D), NOT in
    # Table 3: "The calculated PD parameters of nifedipine were as
    # follows: ke0 = 0.0130 min-1, EC50 = 522 ng/mL, E0 = 169 msec,
    # Emax = 29.5 msec, and gamma = 1.92." Same stimulatory sigmoid Emax
    # form as QT. The companion captopril QTc model could not be fitted
    # (transient HR increase), so QTc appears only in this file.
    # ==================================================================
    lke0_qtc  <- log(0.0130);   label("Effect-compartment equilibration rate ke0 for QTc (1/min)")         # Results text: QTc ke0 = 0.0130 1/min
    le0_qtc   <- log(169);      label("Baseline Bazett-corrected QT before dosing, E0 (msec)")             # Results text: QTc E0 = 169 msec
    lemax_qtc <- log(29.5);     label("Maximum QTc prolongation, Emax (msec)")                             # Results text: QTc Emax = 29.5 msec
    lec50_qtc <- log(522);      label("Effect-compartment concentration at half-maximal QTc effect, EC50 (ng/mL)") # Results text: QTc EC50 = 522 ng/mL
    lhill_qtc <- log(1.92);     label("Hill coefficient gamma for the QTc sigmoid (unitless)")             # Results text: QTc gamma = 1.92

    # ==================================================================
    # RESIDUAL ERROR. The paper fits mean profiles by naive pooling in
    # WinNonlin and reports only AIC -- no residual-error magnitude and no
    # inter-individual variability for any endpoint. Every residual SD is
    # therefore fixed(0) rather than invented, and the model has no eta
    # terms: it reproduces the published typical-value curves exactly.
    # ==================================================================
    addSd     <- fixed(0); label("Additive residual SD on plasma nifedipine (ng/mL; not reported by the paper)")
    addSd_bp  <- fixed(0); label("Additive residual SD on BP (mmHg; not reported by the paper)")
    addSd_hr  <- fixed(0); label("Additive residual SD on HR (beat/min; not reported by the paper)")
    addSd_qt  <- fixed(0); label("Additive residual SD on QT (msec; not reported by the paper)")
    addSd_qtc <- fixed(0); label("Additive residual SD on QTc (msec; not reported by the paper)")
  })

  model({
    # ---- 1. Individual parameters ----
    vc  <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    ke0_bp      <- exp(lke0_bp)
    e0_bp       <- exp(le0_bp)
    emax_bp     <- exp(lemax_bp)
    ec50_bp     <- exp(lec50_bp)
    kin_ehs     <- exp(lkin_ehs)
    kout_ehs    <- exp(lkout_ehs)
    gamma_bp    <- exp(lgamma_bp)
    emax_bp_ehs <- exp(lemax_bp_ehs)
    ec50_bp_ehs <- exp(lec50_bp_ehs)

    ke0_hr   <- exp(lke0_hr)
    e0_hr    <- exp(le0_hr)
    emax_hr  <- exp(lemax_hr)
    ec50_hr  <- exp(lec50_hr)
    hill_hr  <- exp(lhill_hr)

    ke0_qt   <- exp(lke0_qt)
    e0_qt    <- exp(le0_qt)
    emax_qt  <- exp(lemax_qt)
    ec50_qt  <- exp(lec50_qt)
    hill_qt  <- exp(lhill_qt)

    ke0_qtc  <- exp(lke0_qtc)
    e0_qtc   <- exp(le0_qtc)
    emax_qtc <- exp(lemax_qtc)
    ec50_qtc <- exp(lec50_qtc)
    hill_qtc <- exp(lhill_qtc)

    # ---- 2. Two-compartment PK (Methods 2.4 / Figure 1) ----
    # The paper writes the PK ODEs in CONCENTRATION units:
    #   dC1/dt = k0/V1 - (k12 + k10)*C1 + k21*C2
    #   dC2/dt = k12*C1 - k21*C2                    (C2 scaled to V1)
    # Multiplying through by V1 gives the amount form below, which lets
    # rxode2's own zero-order infusion machinery supply k0 = D / Tinf
    # exactly (dose the `central` compartment with dur = 30 min).
    d/dt(central)     <- -(k12 + kel) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    Cc <- central / vc

    # ---- 3. Blood pressure: hypotensive Emax + EHS homeostatic rebound ----
    d/dt(effect_bp) <- ke0_bp * (Cc - effect_bp)

    # Closed-form integral of the paper's dE/dt (derivation in the ini()
    # notes above). Computed before d/dt(ehs) because the EHS production
    # term feeds back on the current BP.
    bp <- (e0_bp * ec50_bp     + emax_bp     * effect_bp) / (ec50_bp     + effect_bp) +
          (e0_bp * ec50_bp_ehs + emax_bp_ehs * ehs)       / (ec50_bp_ehs + ehs) -
          e0_bp

    # EHS pool, expressed as the DEVIATION from its own drug-free steady
    # state. The paper prints
    #   dC4/dt = K_EHS,in * (1 + (E0 - E)/E0 * alpha) - k_EHS,out * C4
    # which is the ODE for TOTAL EHS, whose drug-free steady state is
    # K_EHS,in / k_EHS,out = 1.824 units. The paper also states the initial
    # condition "at t = 0, C3 = 0, C4 = 0, and E = E0" -- and E = E0 at
    # C4 = 0 only holds if C4 is measured from that steady state.
    # Substituting delta = C4_total - K_EHS,in/k_EHS,out into the printed
    # ODE cancels the constant term exactly:
    #   d(delta)/dt = K_EHS,in * (E0 - E)/E0 * alpha - k_EHS,out * delta,
    #   delta(0) = 0
    # so the line below IS the printed equation, written for the state the
    # printed initial condition refers to. Retaining the "1 +" while also
    # starting at 0 would double-count the baseline production: BP would
    # then drift to 131 mmHg with NO drug present, reach a nadir of only
    # 78.5 mmHg (the paper measured 70.9 +/- 8.6), and end at 119 mmHg at
    # 180 min -- 22 mmHg ABOVE the untreated baseline that Figure 3A shows
    # recovering only to about 85-88 mmHg. The deviation form reproduces
    # all three (drug-free BP held at 96.8 mmHg, nadir 71.7 mmHg at
    # 15 min, 84.3 mmHg at 180 min). The vignette asserts the drug-free
    # case as a runtime gate.
    d/dt(ehs) <- kin_ehs * gamma_bp * (e0_bp - bp) / e0_bp - kout_ehs * ehs

    # ---- 4. Secondary endpoints: ordinary sigmoid Emax (Figure 1B) ----
    d/dt(effect_hr) <- ke0_hr * (Cc - effect_hr)
    hr <- e0_hr - emax_hr * effect_hr^hill_hr / (ec50_hr^hill_hr + effect_hr^hill_hr)

    d/dt(effect_qt) <- ke0_qt * (Cc - effect_qt)
    qt <- e0_qt + emax_qt * effect_qt^hill_qt / (ec50_qt^hill_qt + effect_qt^hill_qt)

    d/dt(effect_qtc) <- ke0_qtc * (Cc - effect_qtc)
    qtc <- e0_qtc + emax_qtc * effect_qtc^hill_qtc / (ec50_qtc^hill_qtc + effect_qtc^hill_qtc)

    # ---- 5. Observations ----
    Cc  ~ add(addSd)
    bp  ~ add(addSd_bp)
    hr  ~ add(addSd_hr)
    qt  ~ add(addSd_qt)
    qtc ~ add(addSd_qtc)
  })
}
