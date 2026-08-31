Darwish_2012_armodafinil <- function() {
  description <- "Population pharmacokinetic and MSLT pharmacodynamic model for oral armodafinil in patients with excessive sleepiness associated with shift work disorder (SWD) (Darwish 2012). One-compartment first-order-absorption disposition with a study-specific absorption rate constant and an apparent volume expressed per kg of body weight (Vd/F 0.47 L/kg, ke 0.054 1/h), plus a linear centred body-weight effect on Vd/F and a Study-1 volume offset. The pharmacodynamic layer treats the Multiple Sleep Latency Test (MSLT) as a time-to-event readout: sleep latency is exponentially distributed with a hazard whose logarithm is a cubic polynomial in clock time over the 0000-0800 h testing window, and the predicted plasma concentration reduces that hazard through an inhibitory Emax model with EC50 4.6 ug/mL shared with modafinil. Derived outputs `hazard` (risk of falling asleep, 1/min) and `mslt` (expected sleep latency truncated at the 20-minute session length, min). Companion model from the same paper: modellib('Darwish_2012_modafinil')."
  reference <- paste(
    "Darwish M, Bond M, Ezzet F.",
    "Armodafinil and modafinil in patients with excessive sleepiness",
    "associated with shift work disorder: a pharmacokinetic/pharmacodynamic",
    "model for predicting and comparing their concentration-effect",
    "relationships.",
    "J Clin Pharmacol. 2012;52(9):1328-1342.",
    "doi:10.1177/0091270011417825.",
    "Sister model file from the same paper:",
    "modellib('Darwish_2012_modafinil').",
    sep = " "
  )
  vignette <- "Darwish_2012_armodafinil_modafinil_shift_work_disorder"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "armodafinil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "armodafinil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the model twice. (1) Pharmacokinetics: Darwish 2012 Eq. 3 expresses the apparent volume of distribution per kg of body weight, Vd_i = [Vd * BW] * exp(eta_Vd), and Darwish 2012 Table II adds a linear weight effect on the per-kg volume itself (-0.003 L/kg per kg). The Methods state that 'Observed BW and the other covariates (ie, sex, age, and race) were centered using their median values', but the median is not printed; 86 kg is used as the centring weight (see the model file's Implementation notes for the derivation). (2) Pharmacodynamics: Darwish 2012 Table III carries a body-weight coefficient of 0.004 on the log baseline hazard, footnoted 'Coefficient of body weight >86 kg'. Body-weight range in the armodafinil efficacy trial was 45-153 kg (Darwish 2012 Results, 'Predicted Steady-State Plasma Concentration-Versus-Time Profiles').",
      source_name        = "BW"
    ),
    TCLOCK = list(
      description        = "Wall-clock time of day of the MSLT session, in decimal hours on a 0-24 scale.",
      units              = "hour of day (0-24)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per record. This is the argument `t` of the cubic baseline-hazard polynomial of Darwish 2012 Eq. 5, where t is clock time within the MSLT testing interval, 'usually from 0000 to 0800 h'. In the armodafinil efficacy trial the MSLT was run as five 20-minute sessions every 2 hours from 0000 to 0800 h, so the polynomial is supported on TCLOCK = 0, 2, 4, 6, 8. Darwish 2012 Discussion states that 'Extrapolation to times later than 0800 h was not considered appropriate'; do not evaluate the hazard outside 0-8 h. Distinct from the model's integration-time axis, which carries time after dose (dosing was around 2200 h, i.e. 2 hours before TCLOCK = 0).",
      source_name        = "t (clock time in Darwish 2012 Eq. 5)"
    ),
    STUDY_MODAF = list(
      description        = "Integer study identifier of the Darwish 2012 pooled analysis; the levels armodafinil was administered in are 1, 2, 6 (pharmacokinetic studies) and 7 (the MSLT efficacy trial).",
      units              = "(integer 1-8)",
      type               = "categorical",
      reference_category = "6 for the pharmacokinetic sub-model (SWD-patient study); 7 for the MSLT hazard sub-model (the armodafinil efficacy trial, which is the baseline-hazard reference)",
      notes              = "Time-fixed per subject. Selects the study-specific absorption rate constant of Darwish 2012 Table II ('the population pharmacokinetic model detected no difference in ka between armodafinil and modafinil, [but] ka differed across studies'). Armodafinil was dosed in Study 1 (ka 2.04 1/h), Study 2 (ka 0.41 1/h) and Study 6 (ka 0.79 1/h); Study 1 additionally carries a +0.13 L/kg apparent-volume offset (Table II footnote a, 'Difference from population average'). Study 7 is the armodafinil MSLT efficacy trial, which contributed no plasma concentrations -- its patients' concentrations were predicted with the SWD-patient (Study 6) absorption rate, so levels 6, 7 and 8 all select ka = 0.79 1/h. Use STUDY_MODAF = 7 when simulating the pharmacodynamic layer, and set it to 7 for the sister modafinil model too when reproducing the paper's like-for-like drug comparison (Darwish 2012 Figure 8) under a common baseline hazard.",
      source_name        = "Study no. (Darwish 2012 Table I; Studies 7 and 8 per the paper's Methods text)"
    )
  )

  # Screened by the source paper but NOT retained in either the population
  # pharmacokinetic or the pharmacokinetic/pharmacodynamic model, and
  # therefore not referenced in model(). Darwish 2012 Methods: 'Observed BW
  # and the other covariates (ie, sex, age, and race) were centered using
  # their median values'; Results: 'Previously established and validated
  # population models (data on file, Cephalon) indicated that only body
  # weight (BW) was a significant covariate influencing the drugs'
  # pharmacokinetics'; Discussion: 'Testing of subject characteristics
  # revealed no significant effect, except for a modest effect of BW.'
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a covariate in both the population pharmacokinetic and the MSLT hazard model; not retained. No point estimate is reported."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a covariate in both sub-models; not retained. No point estimate is reported."
    ),
    RACE_BLACK = list(
      description = "Race indicator used to represent the paper's unspecified race covariate screen.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Darwish 2012 names 'race' among the screened covariates but neither defines its categories nor reports an estimate; no race effect is retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 284L,
    n_studies      = 6L,
    age_range      = "not reported in the source paper",
    weight_range   = "45-153 kg in the armodafinil efficacy trial (Study 7); 52-138 kg in the modafinil efficacy trial (Study 8)",
    weight_median  = "not reported; 86 kg is used as the centring weight (see Implementation notes)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Pharmacokinetic model: 190 healthy volunteers (Studies 1-5) pooled with 94 patients with excessive sleepiness associated with shift work disorder (Study 6). Pharmacodynamic (MSLT) model: 449 patients with excessive sleepiness associated with shift work disorder from two 3-month randomised, double-blind, placebo-controlled efficacy trials (Study 7, armodafinil 150 mg/d, 245 patients; Study 8, modafinil 200 mg/d, 204 patients), contributing close to 7000 MSLT observations.",
    dose_range     = "Armodafinil single doses of 50, 100, 200, 300 or 400 mg (Study 1) and 100, 150, 200 or 300 mg/d (Study 2); multiple doses of 150, 200 or 250 mg/d (Study 6); 150 mg/d in the efficacy trial (Study 7). All efficacy-trial doses were given around 2200 h, 30 to 60 minutes before each night shift.",
    regions        = "not reported in the source paper",
    biomarkers     = "Multiple Sleep Latency Test (MSLT) sleep latency in minutes, recorded as five 20-minute sessions every 2 hours from 0000 to 0800 h at baseline and at weeks 4, 8 and 12. Values are right-censored at the 20-minute session length.",
    notes          = "Pharmacokinetic dataset: 3773 plasma concentration observations from 284 subjects across the six studies of Darwish 2012 Table I. Population pharmacokinetic model fitted by maximum-likelihood NLME in S-PLUS 8.0 (TIBCO); pharmacokinetic/pharmacodynamic model fitted in NONMEM VI level 1 with the first-order conditional (Laplace) method, accounting for right censoring at 20 minutes. Dropouts: 38 of 254 randomised (15%) in Study 7 and 16 of 209 (7.7%) in Study 8; last observation carried forward."
  )

  # Implementation notes (see the vignette's 'Assumptions and deviations'
  # section for the full justification of each item):
  #
  # * Centring weight 86 kg. Darwish 2012 Methods state that body weight was
  #   centred at its median but do not print the median. 86 kg is the only
  #   body-weight landmark the paper prints (Table III footnote b,
  #   'Coefficient of body weight >86 kg'), and it is corroborated
  #   independently by the paper's own reported exposure: with Vd/F 0.47 L/kg
  #   and the Study-6 ka of 0.79 1/h, a steady-state armodafinil 200 mg daily
  #   profile stays above the EC50 of 4.6 ug/mL for 8.5 h at 86 kg, matching
  #   the Results statement that 'the plasma concentration of armodafinil
  #   200 mg was maintained ... above the EC50 for up to 9 hours after
  #   dosing' (12.3 h at 70 kg and 7.6 h at 90 kg both miss it).
  #
  # * Body-weight effect on the baseline hazard. Darwish 2012 Table III
  #   footnote b reads 'Coefficient of body weight >86 kg'. This is encoded as
  #   a threshold (hockey-stick) term, e_wt_hazard * max(0, WT - 86), so the
  #   effect is zero at or below 86 kg. The alternative reading -- plain
  #   centring at 86 kg, e_wt_hazard * (WT - 86), which would also *lower* the
  #   hazard below 86 kg -- is not excluded by the paper text; the two agree
  #   above 86 kg and diverge below it. The threshold reading is used because
  #   the paper's median-centring statement is made only about the population
  #   pharmacokinetic model, and the '>' in the footnote is the sole guidance
  #   given for the hazard model. See the vignette Errata.
  #
  # * Apparent-parameter framework. Darwish 2012 Methods: 'Since the Vd of the
  #   central compartment in a 2-compartment model (V1) and the bioavailability
  #   fraction (F) cannot be estimated simultaneously, estimates of CL and V1
  #   were expressed as CL/F and V1/F.' Every volume and clearance here is
  #   therefore apparent (per unit bioavailability).
  #
  # * Common size / IIV scaling. Darwish 2012 Eq. 2 puts the random effect on
  #   the parameter itself (theta_i = theta * exp(eta_i)) and Eq. 3 applies it
  #   to Vd; Table II reports random effects only on ka and Vd/F. The
  #   elimination rate constant ke is a subject-independent fixed effect, so
  #   the individual apparent clearance scales with the individual apparent
  #   volume. model() therefore forms one scaling factor `fvd` from body
  #   weight, the centred weight effect and the volume random effect, and
  #   multiplies both the per-kg volume and the per-kg clearance by it. This
  #   reproduces ke = cl/vc = 0.054 1/h exactly for every subject.
  #
  # * Clearance value. Table II reports ke and Vd/F rather than CL/F; the
  #   product is entered here so both published numbers stay visible in the
  #   source trace. It reproduces the Table II footnote 'Armodafinil
  #   CL = 0.025 L/h/kg' (0.054 * 0.47 = 0.0254).
  #
  # * Residual error. Darwish 2012 Eq. 1 models concentration with an ADDITIVE
  #   error, C_ij(t) = C_pred,ij(t) + eps_ij(t), and Table II reports
  #   SD(eps) = 0.7 ug/mL. There is no proportional component.
  #
  # * Correlated random effects. The Methods state that the random effects of
  #   CL and Vd 'were assumed to be correlated, but independent of eta of the
  #   absorption rate constant (ka)'. Table II of the final model reports only
  #   SD(eta_ka) and SD(eta_V/F) and no covariance, and carries no random
  #   effect on clearance at all, so the random effects are encoded as
  #   uncorrelated. No covariance is invented.
  #
  # * MSLT likelihood. Darwish 2012 Eq. 4 gives sleep latency an exponential
  #   density with hazard h, fitted with right censoring at the 20-minute
  #   session length. That event likelihood has no nlmixr2 residual-error
  #   analogue, so this file is written for forward simulation: `hazard` and
  #   `mslt` are exposed as derived outputs, where `mslt` is the mean of the
  #   exponential distribution truncated at 20 minutes,
  #   E[min(X, 20)] = (1 - exp(-20 h)) / h, and the only error model attaches
  #   to the plasma-concentration output Cc.
  ini({
    # ---- Absorption: study-specific first-order rate constant (Table II) ----
    # Stratum-suffixed because ka is estimated separately in every study; the
    # paper found no armodafinil-versus-modafinil difference in ka. Only the
    # studies in which armodafinil was administered are carried here.
    lka_s1 <- log(2.04); label("First-order absorption rate constant ka in Study 1 (1/h)")                       # Darwish 2012 Table II 'ka Study 1' = 2.04 (SE 0.33, 16% CV)
    lka_s2 <- log(0.41); label("First-order absorption rate constant ka in Study 2 (1/h)")                       # Darwish 2012 Table II 'ka Study 2' = 0.41 (SE 0.033, 8% CV)
    lka_s6 <- log(0.79); label("First-order absorption rate constant ka in Study 6, SWD patients (1/h)")         # Darwish 2012 Table II 'ka Study 6' = 0.79 (SE 0.071, 9% CV); also applied to the efficacy trials (Studies 7 and 8), which had no PK sampling

    # ---- Disposition: apparent volume and clearance, per kg of body weight ----
    lvc <- log(0.47);        label("Apparent volume of distribution Vd/F per kg of body weight (L/kg)")          # Darwish 2012 Table II 'Vd/F Armodafinil' = 0.47 L/kg (SE 0.011, 2% CV)
    lcl <- log(0.054 * 0.47); label("Apparent clearance CL/F per kg of body weight (L/h/kg)")                    # Darwish 2012 Table II 'ke Armodafinil' = 0.054 1/h x 'Vd/F Armodafinil' = 0.47 L/kg; reproduces the Table II footnote 'Armodafinil CL = 0.025 L/h/kg'

    # ---- Covariate effects on the apparent volume (Table II) ----
    e_wt_vc     <- -0.003; label("Slope of body weight on Vd/F per kg, centred at 86 kg (L/kg per kg)")          # Darwish 2012 Table II 'Weight effect on Vd/F' = -0.003 (SE 0.0003, 10% CV); centring weight per Implementation notes
    e_study1_vc <-  0.13;  label("Study-1 offset on Vd/F per kg (L/kg)")                                         # Darwish 2012 Table II 'Vd/F Armodafinil (Study 1)' = 0.13 (SE 0.025, 19% CV), footnote a 'Difference from population average'

    # ---- Inter-individual variability (exponential, Darwish 2012 Eq. 2) ----
    # Table II reports the random effects as standard deviations of eta;
    # ini() takes variances, so each entry is SD^2.
    etalka ~ 0.4356  # Darwish 2012 Table II 'SD(eta_ka)'  = 0.66; 0.66^2 = 0.4356
    etalvc ~ 0.0324  # Darwish 2012 Table II 'SD(eta_V/F)' = 0.18; 0.18^2 = 0.0324

    # ---- Residual error (additive, Darwish 2012 Eq. 1) ----
    addSd <- 0.7; label("Additive residual error on plasma armodafinil concentration (ug/mL)")                   # Darwish 2012 Table II 'SD(eps)' = 0.7 ug/mL

    # ---- MSLT baseline (placebo) hazard: cubic polynomial in clock time ----
    # Darwish 2012 Eq. 5: h_p,i(t) = exp(b0 + b1 t + b2 t^2 + b3 t^3 + eta_i),
    # with t the clock time of the MSLT session within the 0000-0800 h window.
    b0 <- -1.67;   label("Intercept of the log baseline hazard at 0000 h (log 1/min)")                           # Darwish 2012 Table III 'b0 log, min-1' = -1.67 (SE 0.06, 3% CV)
    b1 <-  0.15;   label("Linear clock-time coefficient of the log baseline hazard (log 1/min per h)")           # Darwish 2012 Table III 'b1 log, min-1.h-1' = 0.15 (SE 0.04, 30% CV)
    b2 <-  0.048;  label("Quadratic clock-time coefficient of the log baseline hazard (log 1/min per h^2)")      # Darwish 2012 Table III 'b2 log, min-1.h-2' = 0.048 (SE 0.01, 28% CV)
    b3 <- -0.0065; label("Cubic clock-time coefficient of the log baseline hazard (log 1/min per h^3)")          # Darwish 2012 Table III 'b3 log, min-1.h-3' = -0.0065 (SE 0.001, 17% CV)

    e_wt_hazard <- 0.004; label("Slope of body weight above 86 kg on the log baseline hazard (log 1/min per kg)") # Darwish 2012 Table III 'Body weight log, min.kg-1' = 0.004 (SE 0.002, 50% CV), footnote b 'Coefficient of body weight >86 kg'

    # ---- MSLT drug effect (Darwish 2012 Eq. 6) ----
    lec50 <- 1.53; label("Log of the concentration producing 50% of the maximum hazard reduction (log ug/mL)")   # Darwish 2012 Table III 'Log EC50' = 1.53 (SE 0.15, 10% CV); Table III also prints the back-transformed 'EC50 = 4.6 ug/mL' and the Results state the EC50 is assumed equal for armodafinil and modafinil

    # ---- Inter-individual variability in the baseline hazard ----
    etab0 ~ 0.77  # Darwish 2012 Table III 'sigma_eta^2 log, min-2' = 0.77 (SE 0.083, 11% CV), reported as a VARIANCE
  })

  model({
    # Centring / threshold body weight shared by the pharmacokinetic and the
    # pharmacodynamic weight effects (see Implementation notes).
    wt_ref <- 86

    # ---- Study-specific absorption rate constant (Table II) ----
    # Studies 6, 7 and 8 all take the Study-6 (SWD-patient) value: the two
    # efficacy trials contributed MSLT data only, and their patients'
    # concentrations were predicted with the patient-study absorption rate.
    lka_study <- lka_s1 * (STUDY_MODAF == 1) +
      lka_s2 * (STUDY_MODAF == 2) +
      lka_s6 * (STUDY_MODAF >= 6)
    ka <- exp(lka_study + etalka)

    # ---- Apparent volume per kg: population value, linear centred weight
    #      effect, and the Study-1 offset (Table II) ----
    vdkg <- exp(lvc) + e_wt_vc * (WT - wt_ref) + e_study1_vc * (STUDY_MODAF == 1)

    # One scaling factor shared by the apparent volume and the apparent
    # clearance, so that the subject-independent elimination rate constant
    # ke = cl / vc is preserved exactly (see Implementation notes).
    fvd <- (vdkg / exp(lvc)) * WT * exp(etalvc)
    vc  <- exp(lvc) * fvd
    cl  <- exp(lcl) * fvd
    kel <- cl / vc

    # ---- One-compartment disposition with first-order absorption ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    # ---- MSLT hazard model ----
    # Baseline (placebo) hazard, Darwish 2012 Eq. 5, in units of 1/min. The
    # body-weight term is the threshold form of Table III footnote b.
    hp <- exp(b0 + b1 * TCLOCK + b2 * TCLOCK^2 + b3 * TCLOCK^3 +
                e_wt_hazard * max(0, WT - wt_ref) + etab0)

    # Inhibitory drug effect, Darwish 2012 Eq. 6: the predicted plasma
    # concentration reduces the hazard, to zero for large Cc and to hp for
    # Cc = 0.
    ec50   <- exp(lec50)
    hazard <- hp * (1 - Cc / (ec50 + Cc))

    # Expected MSLT: the mean of the Eq. 4 exponential sleep-latency
    # distribution truncated at the 20-minute session length.
    mslt <- (1 - exp(-20 * hazard)) / hazard

    Cc ~ add(addSd)
  })
}
