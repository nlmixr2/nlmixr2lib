Kiriyama_2024_captopril_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment IV-infusion PK of captopril in spontaneously hypertensive rats (SHR/Izm)",
    "linked to three effect-compartment pharmacodynamic endpoints -- mean blood pressure (BP),",
    "heart rate (HR) and QT interval (QT) -- each an ordinary sigmoid Emax model with its own",
    "effect-compartment equilibration rate ke0.",
    "Unlike its nifedipine sister model, captopril's BP response needed no homeostatic",
    "counter-regulation arm; the HR and QT ke0 values are so small (1e-5 /min) that those two",
    "effect compartments behave as cumulative-exposure integrators.",
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
    "modellib('Kiriyama_2024_nifedipine_rat').",
    sep = " "
  )
  vignette <- "Kiriyama_2024_nifedipine_captopril_rat"
  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  # No covariates: the paper fits mean profiles with absolute (not weight-normalised)
  # V1 and micro-constants, so body weight enters only when converting the reported
  # mg/kg dose into the absolute infused amount. See population$notes.
  covariateData <- list()

  compartmentData <- list(
    central     = list(analyte = "captopril", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "captopril", units = "ng", specimen = "plasma", verified = TRUE),
    effect_bp   = list(analyte = "captopril", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_hr   = list(analyte = "captopril", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_qt   = list(analyte = "captopril", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  paper_specific_compartments <- c(
    # Three parallel (NOT chained) effect compartments, one per PD endpoint,
    # because the paper estimates a separate ke0 for each endpoint (Table 3).
    # The blessed `effect<n>` numbered form would falsely imply a chain, so each
    # carries its endpoint suffix instead. Precedent for an `effect_<suffix>`
    # compartment: Kleijn_2011_sugammadex_rocuronium.R (`effect_roc`).
    "effect_bp", "effect_hr", "effect_qt"
  )

  population <- list(
    species        = "rat (spontaneously hypertensive rat, SHR/Izm; male)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "not reported",
    weight_range   = "215-315 g",
    sex_female_pct = 0,
    disease_state  = "spontaneous hypertension (SHR/Izm strain) under urethane anaesthesia (1.0 g/kg intraperitoneal)",
    dose_range     = "captopril 15.0 mg/kg (monotherapy) or 5.0 mg/kg (coadministered with nifedipine), 30-minute IV infusion into the right femoral vein in 0.5 mL polyethylene glycol 400",
    regions        = "Japan (Doshisha Women's College of Liberal Arts, Kyoto)",
    biomarkers     = paste(
      "Mean arterial BP via a left-carotid-artery catheter; HR and QT via ECG;",
      "both recorded continuously from -15 to 180 min after infusion start.",
      "Plasma captopril sampled from the right jugular vein 5 to 300 min after infusion start",
      "and assayed by LC-MS/MS after derivatisation with 2,4'-dibromoacetophenone to stabilise",
      "the thiol (captopril-BPB, m/z 416 -> 255).",
      sep = " "
    ),
    notes = paste(
      "Per-analysis group sizes differ: PK parameters (Tables 1-2) come from 3 to 6 rats,",
      "PD parameters (Table 3) from 8 to 10 rats; each experiment used 3-10 rats.",
      "n_subjects records the PD maximum (10).",
      "DOSE UNITS: the paper reports mg/kg doses but fits absolute V1 (mL) and absolute",
      "micro-constants, so a body weight is needed to build an event table. Back-solving the",
      "Figure 2 fitted peak concentrations against the Table 2 monotherapy parameters gives",
      "249 g from the captopril panel and 247 g from the nifedipine panel -- two independent",
      "estimates agreeing to within 1% and both inside the stated 215-315 g range. The",
      "validation vignette therefore doses a 250 g rat (15 mg/kg = 3.75 mg = 3.75e6 ng).",
      "The reported body weight is NOT used as a model covariate.",
      "The transient HR RISE seen immediately after captopril infusion (measured",
      "Emax = 325 +/- 35 beat/min against a 312.8 beat/min baseline) is NOT part of the model:",
      "the sigmoid Emax fit captures only the subsequent HR reduction, and the paper states",
      "the QTc relationship could not be quantified at all because of that transient.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETICS -- Kiriyama 2024 Table 2, captopril "Monotherapy"
    # column (15.0 mg/kg). The 2-compartment structure and the ODEs are
    # Methods section 2.4 / Figure 1. All four values are WinNonlin
    # least-squares estimates on the MEAN concentration-time profile
    # (AIC 218.2); the paper reports no standard errors, so nothing is
    # wrapped in fixed().
    #
    # The paper also reports a separate coadministration PK fit
    # (V1 = 41.7 mL, k10 = 0.0842, k12 = 0.0631, k21 = 0.0529 /min;
    # AIC 188.6) used for its drug-interaction simulations. Those values
    # are NOT encoded as a covariate effect because the paper explicitly
    # concludes that none of the PK differences reached significance
    # ("no significant PK interactions were observed"); the vignette
    # applies them through rxSolve(params = ) instead.
    # ==================================================================
    lvc  <- log(54.2);    label("Central volume V1 (mL)")                                   # Table 2 captopril monotherapy V1 = 54.2 mL
    lkel <- log(0.0521);  label("Elimination rate constant k10 (1/min)")                    # Table 2 captopril monotherapy k10 = 0.0521 1/min
    lk12 <- log(0.153);   label("Central-to-peripheral rate constant k12 (1/min)")          # Table 2 captopril monotherapy k12 = 0.153 1/min
    lk21 <- log(0.0694);  label("Peripheral-to-central rate constant k21 (1/min)")          # Table 2 captopril monotherapy k21 = 0.0694 1/min

    # ==================================================================
    # PHARMACODYNAMICS -- Kiriyama 2024 Table 3, "Captopril" rows. All
    # three endpoints use the ordinary sigmoid Emax model with an effect
    # compartment of Figure 1B (Methods 2.5):
    #   dC3/dt = ke0 * (C1 - C3)
    #   BP, HR:  E = E0 - Emax * C3^gamma / (EC50^gamma + C3^gamma)
    #   QT:      E = E0 + Emax * C3^gamma / (EC50^gamma + C3^gamma)
    # Emax here is a CHANGE from baseline (contrast the nifedipine BP
    # model, where the printed dE/dt makes Emax an asymptotic level).
    # ==================================================================

    # Blood pressure (AIC 213.2). Floor = 105 - 41.8 = 63.2 mmHg; the
    # simulated nadir is 81.6 mmHg at 41 min, matching Figure 3A.
    lke0_bp  <- log(0.0400);  label("Effect-compartment equilibration rate ke0 for BP (1/min)")                        # Table 3 captopril BP ke0 = 0.0400 1/min
    le0_bp   <- log(105);     label("Baseline BP before dosing, E0 (mmHg)")                                            # Table 3 captopril BP E0 = 105 mmHg
    lemax_bp <- log(41.8);    label("Maximum BP reduction, Emax (mmHg)")                                               # Table 3 captopril BP Emax = 41.8 mmHg
    lec50_bp <- log(7885);    label("Effect-compartment concentration at half-maximal BP effect, EC50 (ng/mL)")        # Table 3 captopril BP EC50 = 7885 ng/mL
    lhill_bp <- log(0.662);   label("Hill coefficient gamma for the BP sigmoid (unitless)")                            # Table 3 captopril BP gamma = 0.662

    # Heart rate (AIC 299). ke0 = 1.21e-5 /min corresponds to an
    # effect-compartment equilibration half-life of 5.7e4 min, i.e. about
    # 40 days -- over the 180 min experiment the effect compartment
    # integrates plasma exposure rather than equilibrating with it. Paired
    # with EC50 = 7.6 ng/mL the sigmoid is near-saturated by 180 min,
    # giving 255 beat/min against the ~250 beat/min plateau in Figure 3B.
    # This is a faithful transcription of Table 3, not a unit slip: the
    # combination reproduces the published curve.
    lke0_hr  <- log(1.21e-5); label("Effect-compartment equilibration rate ke0 for HR (1/min)")                        # Table 3 captopril HR ke0 = 1.21e-5 1/min
    le0_hr   <- log(312.8);   label("Baseline heart rate before dosing, E0 (beat/min)")                                # Table 3 captopril HR E0 = 312.8 beat/min
    lemax_hr <- log(80.4);    label("Maximum heart-rate reduction, Emax (beat/min)")                                   # Table 3 captopril HR Emax = 80.4 beat/min
    lec50_hr <- log(7.6);     label("Effect-compartment concentration at half-maximal HR effect, EC50 (ng/mL)")        # Table 3 captopril HR EC50 = 7.6 ng/mL
    lhill_hr <- log(1.4);     label("Hill coefficient gamma for the HR sigmoid (unitless)")                            # Table 3 captopril HR gamma = 1.4

    # QT interval (AIC 220.0). Same integrator behaviour as HR
    # (ke0 = 4.05e-5 /min); the very shallow gamma = 0.1347 is what makes
    # captopril's QT prolongation rise faster than nifedipine's, the
    # qualitative contrast the paper highlights in its Results.
    lke0_qt  <- log(4.05e-5); label("Effect-compartment equilibration rate ke0 for QT (1/min)")                        # Table 3 captopril QT ke0 = 4.05e-5 1/min
    le0_qt   <- log(81.2);    label("Baseline QT interval before dosing, E0 (msec)")                                   # Table 3 captopril QT E0 = 81.2 msec
    lemax_qt <- log(41.9);    label("Maximum QT prolongation, Emax (msec)")                                            # Table 3 captopril QT Emax = 41.9 msec
    lec50_qt <- log(4873);    label("Effect-compartment concentration at half-maximal QT effect, EC50 (ng/mL)")        # Table 3 captopril QT EC50 = 4873 ng/mL
    lhill_qt <- log(0.1347);  label("Hill coefficient gamma for the QT sigmoid (unitless)")                            # Table 3 captopril QT gamma = 0.1347

    # ==================================================================
    # RESIDUAL ERROR. The paper fits mean profiles by naive pooling in
    # WinNonlin and reports only AIC -- no residual-error magnitude and no
    # inter-individual variability for any endpoint. Every residual SD is
    # therefore fixed(0) rather than invented, and the model has no eta
    # terms: it reproduces the published typical-value curves exactly.
    # ==================================================================
    addSd    <- fixed(0); label("Additive residual SD on plasma captopril (ng/mL; not reported by the paper)")
    addSd_bp <- fixed(0); label("Additive residual SD on BP (mmHg; not reported by the paper)")
    addSd_hr <- fixed(0); label("Additive residual SD on HR (beat/min; not reported by the paper)")
    addSd_qt <- fixed(0); label("Additive residual SD on QT (msec; not reported by the paper)")
  })

  model({
    # ---- 1. Individual parameters ----
    vc  <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    ke0_bp  <- exp(lke0_bp)
    e0_bp   <- exp(le0_bp)
    emax_bp <- exp(lemax_bp)
    ec50_bp <- exp(lec50_bp)
    hill_bp <- exp(lhill_bp)

    ke0_hr  <- exp(lke0_hr)
    e0_hr   <- exp(le0_hr)
    emax_hr <- exp(lemax_hr)
    ec50_hr <- exp(lec50_hr)
    hill_hr <- exp(lhill_hr)

    ke0_qt  <- exp(lke0_qt)
    e0_qt   <- exp(le0_qt)
    emax_qt <- exp(lemax_qt)
    ec50_qt <- exp(lec50_qt)
    hill_qt <- exp(lhill_qt)

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

    # ---- 3. Effect compartments and sigmoid Emax responses ----
    d/dt(effect_bp) <- ke0_bp * (Cc - effect_bp)
    bp <- e0_bp - emax_bp * effect_bp^hill_bp / (ec50_bp^hill_bp + effect_bp^hill_bp)

    d/dt(effect_hr) <- ke0_hr * (Cc - effect_hr)
    hr <- e0_hr - emax_hr * effect_hr^hill_hr / (ec50_hr^hill_hr + effect_hr^hill_hr)

    d/dt(effect_qt) <- ke0_qt * (Cc - effect_qt)
    qt <- e0_qt + emax_qt * effect_qt^hill_qt / (ec50_qt^hill_qt + effect_qt^hill_qt)

    # ---- 4. Observations ----
    Cc ~ add(addSd)
    bp ~ add(addSd_bp)
    hr ~ add(addSd_hr)
    qt ~ add(addSd_qt)
  })
}
