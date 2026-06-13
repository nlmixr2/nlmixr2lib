LiemMoolenaar_2011_scopolamine <- function() {
  description <- paste(
    "Two-compartment population PK model for IV scopolamine in healthy adults",
    "(Liem-Moolenaar 2011, Table 2) with ten parallel effect-compartment",
    "linear-concentration-effect PK/PD models for central-nervous-system",
    "endpoints (Table 3): heart rate, saccadic peak velocity, adaptive tracking,",
    "VAS external perception, body sway, VAS alertness, VAS internal perception,",
    "smooth pursuit, VAS feeling high, and finger tapping (with an additive",
    "time-trend on finger tapping). PK was fit to 90 healthy male volunteers",
    "given a single 0.5 mg scopolamine i.v. infusion over 15 minutes; the ten",
    "PD endpoints were fit independently as effect-compartment linear-slope",
    "models on the empirical-Bayes individual PK profiles. PD parameter sets",
    "are grouped in Table 3 by equilibration half-life (heart rate <0.5 h;",
    "saccadic peak velocity and adaptive tracking 1-1.5 h; VAS external,",
    "body sway, VAS alertness, VAS internal, smooth pursuit 2.5-3.5 h;",
    "VAS feeling high and finger tapping >8 h)."
  )
  reference <- paste(
    "Liem-Moolenaar M, de Boer P, Timmers M, Schoemaker RC, van Hasselt JGC,",
    "Schmidt S, van Gerven JMA. Pharmacokinetic-pharmacodynamic relationships",
    "of central nervous system effects of scopolamine in healthy subjects.",
    "Br J Clin Pharmacol. 2011;71(6):886-898.",
    "doi:10.1111/j.1365-2125.2011.03936.x"
  )
  vignette <- "LiemMoolenaar_2011_scopolamine"
  units <- list(time = "minute", dosing = "mg", concentration = "pg/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 85L,
    n_studies      = 2L,
    age_range      = "18-55 years",
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = NULL,
    disease_state  = "Healthy male volunteers (90 enrolled across two studies; 85 evaluable for PK/PD).",
    dose_range     = "0.5 mg scopolamine IV infusion over 15 min, single dose (placebo period within a four-way crossover; placebo arm not modelled here).",
    regions        = "Netherlands (Centre for Human Drug Research, Leiden)",
    notes          = paste(
      "BMI inclusion range 18-28.5 kg/m^2. Subjects were medically screened",
      "within 3 weeks prior to enrolment and were not allowed to use",
      "CNS-active medication, smoke, or consume caffeine/alcohol during the",
      "study. Demographics per Methods, 'Subjects' section. The PK fit used",
      "the same population that supplied the per-endpoint PD records to the",
      "Table 3 PK/PD analyses."
    )
  )

  ini({
    # =================================================================
    # PK structural parameters (Liem-Moolenaar 2011 Table 2)
    # Time unit: minute. Concentration unit: pg/mL.
    # =================================================================
    lcl  <- log(2.53);   label("Clearance (L/min)")                                # Table 2: CL = 2.53 L/min
    lvc  <- log(66.3);   label("Central volume of distribution (L)")               # Table 2: Vc = 66.3 L
    lq   <- log(4.78);   label("Inter-compartmental clearance (L/min)")            # Table 2: Q = 4.78 L/min
    lvp  <- log(183.7);  label("Peripheral volume of distribution (L)")            # Table 2: Vss = 250 L; derived Vp = Vss - Vc = 250 - 66.3 = 183.7 L

    # PK IIV (CV%; omega^2 = log(1 + CV^2))
    etalcl ~ 0.023568    # CV 15.4% on CL                                          # Table 2: IIVar CL = 15.4%
    etalvc ~ 0.125855    # CV 36.6% on Vc                                          # Table 2: IIVar Vc = 36.6%
    etalq  ~ 0.007885    # CV  8.9% on Q                                           # Table 2: IIVar Q = 8.9%
    # The paper reports a 7.2% IIV on Vss (=Vc+Vp) rather than on Vp directly;
    # see vignette Assumptions and deviations for the Vp IIV treatment.

    propSd <- 0.102; label("Proportional residual error on scopolamine plasma concentration (fraction)") # Table 2: residual error SD/mean = 0.102

    # =================================================================
    # PK/PD: ten effect-compartment linear-slope models (Table 3)
    # Effect compartments are numbered effect1..effect10 in the order of
    # Table 3 (ascending equilibration half-life):
    #   effect1  = heart rate (HR)
    #   effect2  = saccadic peak velocity (SPV)
    #   effect3  = adaptive tracking
    #   effect4  = VAS external perception (log mm)
    #   effect5  = body sway (log mm)
    #   effect6  = VAS alertness (mm)
    #   effect7  = VAS internal perception (log mm)
    #   effect8  = smooth pursuit (%)
    #   effect9  = VAS feeling high (log mm)
    #   effect10 = finger tapping (taps per 10 s)
    # =================================================================

    # ---- effect1: heart rate (HR), t1/2,keo 16.8 min -----------------
    lintercept_hr <- log(55.2);          label("Heart rate baseline intercept (beats/min)")            # Table 3 HR Intercept
    slope_hr      <- -0.00675;           label("HR slope vs effect-compartment Ce (beats/min per pg/mL)") # Table 3 HR Slope (footnote on sign parameterisation)
    lke0_hr       <- log(log(2) / 16.8); label("HR ke0 (1/min); t1/2,keo = 16.8 min")                  # Table 3 HR t1/2,keo
    etalintercept_hr ~ 0.0002956   # CV 1.72%                                                          # Table 3 HR intercept IICV
    etaslope_hr      ~ 0.0207100   # CV 14.46%                                                         # Table 3 HR slope IICV
    etalke0_hr       ~ 0.0260887   # CV 16.25%                                                         # Table 3 HR t1/2,keo IICV
    propSd_hr <- 0.103; label("Proportional residual error on HR (fraction; LTBS-style interpretation; see vignette)") # Table 3 HR residual error

    # ---- effect2: saccadic peak velocity (SPV), t1/2,keo 65.1 min ----
    lintercept_spv <- log(485);           label("SPV baseline intercept (deg/s)")                       # Table 3 SPV Intercept
    slope_spv      <- -0.0737;            label("SPV slope vs Ce (deg/s per pg/mL)")                    # Table 3 SPV Slope
    lke0_spv       <- log(log(2) / 65.1); label("SPV ke0 (1/min); t1/2,keo = 65.1 min")                  # Table 3 SPV t1/2,keo
    etalintercept_spv ~ 0.265716   # CV 55.05%                                                         # Table 3 SPV intercept IICV
    etaslope_spv      ~ 0.186753   # CV 45.3%                                                          # Table 3 SPV slope IICV
    # No IIV on SPV ke0 (Table 3 reports IICV 0.0%)
    addSd_spv <- 28.4; label("Additive residual error on SPV (deg/s)")                                  # Table 3 SPV residual error

    # ---- effect3: adaptive tracking (placebo-subtracted %), t1/2,keo 86.6 min ----
    lintercept_tracking <- log(0.0479);          label("Adaptive tracking intercept (% change from baseline, placebo-subtracted)") # Table 3 adaptive tracking Intercept
    slope_tracking      <- -0.0217;              label("Adaptive tracking slope vs Ce (% per pg/mL)")                     # Table 3 adaptive tracking Slope
    lke0_tracking       <- log(log(2) / 86.6);   label("Adaptive tracking ke0 (1/min); t1/2,keo = 86.6 min")               # Table 3 adaptive tracking t1/2,keo
    etalintercept_tracking ~ 0.001877   # CV 4.336%                                                                       # Table 3 adaptive tracking intercept IICV
    etaslope_tracking      ~ 0.049917   # CV 22.6%                                                                        # Table 3 adaptive tracking slope IICV
    etalke0_tracking       ~ 0.077385   # CV 28.4%                                                                        # Table 3 adaptive tracking t1/2,keo IICV
    addSd_tracking <- 2.98; label("Additive residual error on adaptive tracking (%)")                                     # Table 3 adaptive tracking residual error

    # ---- effect4: VAS external perception (log mm), t1/2,keo 161 min ----
    lintercept_vasext <- log(0.34);           label("VAS external perception intercept (log mm)")                      # Table 3 VAS external perception Intercept
    slope_vasext      <- 0.000633;            label("VAS external perception slope vs Ce (log mm per pg/mL)")          # Table 3 VAS external perception Slope
    lke0_vasext       <- log(log(2) / 161);   label("VAS external perception ke0 (1/min); t1/2,keo = 161 min")          # Table 3 VAS external perception t1/2,keo
    etalintercept_vasext ~ 0.078594   # CV 28.6%                                                                        # Table 3 VAS external perception intercept IICV
    etaslope_vasext      ~ 0.844218   # CV 107%                                                                         # Table 3 VAS external perception slope IICV
    etalke0_vasext       ~ 0.540045   # CV 86.5%                                                                        # Table 3 VAS external perception t1/2,keo IICV
    addSd_vasext <- 0.0935; label("Additive residual error on VAS external perception (log mm)")                       # Table 3 VAS external perception residual error

    # ---- effect5: body sway (log mm), t1/2,keo 181 min ----
    lintercept_bodysway <- log(2.4);           label("Body sway intercept (log mm)")                                     # Table 3 body sway Intercept
    slope_bodysway      <- 0.00147;            label("Body sway slope vs Ce (log mm per pg/mL)")                        # Table 3 body sway Slope
    lke0_bodysway       <- log(log(2) / 181);  label("Body sway ke0 (1/min); t1/2,keo = 181 min")                        # Table 3 body sway t1/2,keo
    etalintercept_bodysway ~ 0.006693   # CV 8.2%                                                                       # Table 3 body sway intercept IICV
    etaslope_bodysway      ~ 0.099369   # CV 32.4%                                                                      # Table 3 body sway slope IICV
    etalke0_bodysway       ~ 0.089164   # CV 30.6%                                                                      # Table 3 body sway t1/2,keo IICV
    addSd_bodysway <- 0.102; label("Additive residual error on body sway (log mm)")                                     # Table 3 body sway residual error

    # ---- effect6: VAS alertness (mm), t1/2,keo 199 min ----
    lintercept_vasalert <- log(53.6);          label("VAS alertness intercept (mm)")                                    # Table 3 VAS alertness Intercept
    slope_vasalert      <- -0.0622;            label("VAS alertness slope vs Ce (mm per pg/mL)")                        # Table 3 VAS alertness Slope
    lke0_vasalert       <- log(log(2) / 199);  label("VAS alertness ke0 (1/min); t1/2,keo = 199 min")                    # Table 3 VAS alertness t1/2,keo
    etalintercept_vasalert ~ 0.021727   # CV 14.8%                                                                      # Table 3 VAS alertness intercept IICV
    etaslope_vasalert      ~ 0.329597   # CV 62.9%                                                                      # Table 3 VAS alertness slope IICV
    etalke0_vasalert       ~ 0.253522   # CV 53.8%                                                                      # Table 3 VAS alertness t1/2,keo IICV
    addSd_vasalert <- 4.05; label("Additive residual error on VAS alertness (mm)")                                     # Table 3 VAS alertness residual error

    # ---- effect7: VAS internal perception (log mm), t1/2,keo 200 min ----
    lintercept_vasint <- log(0.336);          label("VAS internal perception intercept (log mm)")                       # Table 3 VAS internal perception Intercept
    slope_vasint      <- 0.000331;            label("VAS internal perception slope vs Ce (log mm per pg/mL)")           # Table 3 VAS internal perception Slope
    lke0_vasint       <- log(log(2) / 200);   label("VAS internal perception ke0 (1/min); t1/2,keo = 200 min")           # Table 3 VAS internal perception t1/2,keo
    etalintercept_vasint ~ 0.075853   # CV 28.1%                                                                         # Table 3 VAS internal perception intercept IICV
    etaslope_vasint      ~ 1.047319   # CV 130%                                                                          # Table 3 VAS internal perception slope IICV
    etalke0_vasint       ~ 0.000256   # CV 1.6%                                                                          # Table 3 VAS internal perception t1/2,keo IICV
    addSd_vasint <- 0.083; label("Additive residual error on VAS internal perception (log mm)")                          # Table 3 VAS internal perception residual error

    # ---- effect8: smooth pursuit (%), t1/2,keo 221 min ----
    lintercept_smoothpursuit <- log(3.1);           label("Smooth pursuit intercept (%)")                               # Table 3 smooth pursuit Intercept
    slope_smoothpursuit      <- -0.0264;            label("Smooth pursuit slope vs Ce (% per pg/mL)")                   # Table 3 smooth pursuit Slope
    lke0_smoothpursuit       <- log(log(2) / 221);  label("Smooth pursuit ke0 (1/min); t1/2,keo = 221 min")              # Table 3 smooth pursuit t1/2,keo
    etalintercept_smoothpursuit ~ 0.014908   # CV 12.25%                                                                # Table 3 smooth pursuit intercept IICV
    etaslope_smoothpursuit      ~ 0.179251   # CV 44.4%                                                                 # Table 3 smooth pursuit slope IICV
    etalke0_smoothpursuit       ~ 0.655339   # CV 97.8%                                                                 # Table 3 smooth pursuit t1/2,keo IICV
    addSd_smoothpursuit <- 6.28; label("Additive residual error on smooth pursuit (%)")                                  # Table 3 smooth pursuit residual error

    # ---- effect9: VAS feeling high (log mm), t1/2,keo 483 min ----
    lintercept_vashigh <- log(0.328);           label("VAS feeling high intercept (log mm)")                            # Table 3 VAS feeling high Intercept
    slope_vashigh      <- 0.00313;              label("VAS feeling high slope vs Ce (log mm per pg/mL)")                # Table 3 VAS feeling high Slope
    lke0_vashigh       <- log(log(2) / 483);    label("VAS feeling high ke0 (1/min); t1/2,keo = 483 min")                # Table 3 VAS feeling high t1/2,keo
    etalintercept_vashigh ~ 0.042989   # CV 21.0%                                                                       # Table 3 VAS feeling high intercept IICV
    etaslope_vashigh      ~ 0.631968   # CV 93.9%                                                                       # Table 3 VAS feeling high slope IICV
    etalke0_vashigh       ~ 1.652984   # CV 207%                                                                        # Table 3 VAS feeling high t1/2,keo IICV
    addSd_vashigh <- 0.216; label("Additive residual error on VAS feeling high (log mm)")                                # Table 3 VAS feeling high residual error

    # ---- effect10: finger tapping (taps per 10 s), t1/2,keo 649 min ----
    # Finger tapping uniquely carries an additive time-trend (paper:
    # "Finger tapping has a time slope of 19.3 number per day"). The slope
    # captures learning/practice effects during the visit; it is reported as
    # taps per 10 s per day.
    lintercept_tapping <- log(63.4);          label("Finger tapping intercept (taps per 10 s)")                          # Table 3 finger tapping Intercept
    slope_tapping      <- -0.0797;            label("Finger tapping slope vs Ce (taps per 10 s per pg/mL)")             # Table 3 finger tapping Slope
    lke0_tapping       <- log(log(2) / 649);  label("Finger tapping ke0 (1/min); t1/2,keo = 649 min")                    # Table 3 finger tapping t1/2,keo
    time_slope_tapping <- 19.3;               label("Finger tapping time slope (taps per 10 s per day)")                 # Table 3 finger tapping footnote ** (time slope 19.3 number per day)
    etalintercept_tapping ~ 0.015729   # CV 12.6%                                                                        # Table 3 finger tapping intercept IICV
    etaslope_tapping      ~ 0.061810   # CV 25.3%                                                                        # Table 3 finger tapping slope IICV
    etalke0_tapping       ~ 0.329597   # CV 62.9%                                                                        # Table 3 finger tapping t1/2,keo IICV
    addSd_tapping <- 3.09; label("Additive residual error on finger tapping (taps per 10 s)")                            # Table 3 finger tapping residual error
  })

  model({
    # =================================================================
    # PK: two-compartment IV (time in minutes; concentration in pg/mL)
    # The reported peripheral V is derived as Vss - Vc = 250 - 66.3 L.
    # The published IIV on Vss is not directly reproducible on Vp; only
    # the Vc IIV is applied here (the Vp typical value is fixed and the
    # PK uses no Vp eta) -- see vignette Assumptions and deviations.
    # =================================================================
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Convert central (mg) / vc (L) = mg/L; 1 mg/L = 1e6 pg/mL
    Cc <- central / vc * 1e6
    Cc ~ prop(propSd)

    # =================================================================
    # PK/PD: ten parallel effect-compartment linear-slope models.
    # Each compartment effectN integrates Cc - effectN at rate ke0_<endpoint>.
    # Slope IIV is multiplicative (slope_i = slope_pop * exp(eta)) which
    # preserves the sign of the population slope.
    # =================================================================

    # ---- effect1: heart rate (HR) ----
    intercept_hr_i <- exp(lintercept_hr + etalintercept_hr)
    slope_hr_i     <- slope_hr * exp(etaslope_hr)
    ke0_hr         <- exp(lke0_hr + etalke0_hr)
    d/dt(effect1)  <- ke0_hr * (Cc - effect1)
    hr             <- intercept_hr_i + slope_hr_i * effect1
    hr ~ prop(propSd_hr)

    # ---- effect2: saccadic peak velocity (SPV) ----
    intercept_spv_i <- exp(lintercept_spv + etalintercept_spv)
    slope_spv_i     <- slope_spv * exp(etaslope_spv)
    ke0_spv         <- exp(lke0_spv)
    d/dt(effect2)   <- ke0_spv * (Cc - effect2)
    spv             <- intercept_spv_i + slope_spv_i * effect2
    spv ~ add(addSd_spv)

    # ---- effect3: adaptive tracking ----
    intercept_tracking_i <- exp(lintercept_tracking + etalintercept_tracking)
    slope_tracking_i     <- slope_tracking * exp(etaslope_tracking)
    ke0_tracking         <- exp(lke0_tracking + etalke0_tracking)
    d/dt(effect3)        <- ke0_tracking * (Cc - effect3)
    tracking             <- intercept_tracking_i + slope_tracking_i * effect3
    tracking ~ add(addSd_tracking)

    # ---- effect4: VAS external perception ----
    intercept_vasext_i <- exp(lintercept_vasext + etalintercept_vasext)
    slope_vasext_i     <- slope_vasext * exp(etaslope_vasext)
    ke0_vasext         <- exp(lke0_vasext + etalke0_vasext)
    d/dt(effect4)      <- ke0_vasext * (Cc - effect4)
    vasext             <- intercept_vasext_i + slope_vasext_i * effect4
    vasext ~ add(addSd_vasext)

    # ---- effect5: body sway ----
    intercept_bodysway_i <- exp(lintercept_bodysway + etalintercept_bodysway)
    slope_bodysway_i     <- slope_bodysway * exp(etaslope_bodysway)
    ke0_bodysway         <- exp(lke0_bodysway + etalke0_bodysway)
    d/dt(effect5)        <- ke0_bodysway * (Cc - effect5)
    bodysway             <- intercept_bodysway_i + slope_bodysway_i * effect5
    bodysway ~ add(addSd_bodysway)

    # ---- effect6: VAS alertness ----
    intercept_vasalert_i <- exp(lintercept_vasalert + etalintercept_vasalert)
    slope_vasalert_i     <- slope_vasalert * exp(etaslope_vasalert)
    ke0_vasalert         <- exp(lke0_vasalert + etalke0_vasalert)
    d/dt(effect6)        <- ke0_vasalert * (Cc - effect6)
    vasalert             <- intercept_vasalert_i + slope_vasalert_i * effect6
    vasalert ~ add(addSd_vasalert)

    # ---- effect7: VAS internal perception ----
    intercept_vasint_i <- exp(lintercept_vasint + etalintercept_vasint)
    slope_vasint_i     <- slope_vasint * exp(etaslope_vasint)
    ke0_vasint         <- exp(lke0_vasint + etalke0_vasint)
    d/dt(effect7)      <- ke0_vasint * (Cc - effect7)
    vasint             <- intercept_vasint_i + slope_vasint_i * effect7
    vasint ~ add(addSd_vasint)

    # ---- effect8: smooth pursuit ----
    intercept_smoothpursuit_i <- exp(lintercept_smoothpursuit + etalintercept_smoothpursuit)
    slope_smoothpursuit_i     <- slope_smoothpursuit * exp(etaslope_smoothpursuit)
    ke0_smoothpursuit         <- exp(lke0_smoothpursuit + etalke0_smoothpursuit)
    d/dt(effect8)             <- ke0_smoothpursuit * (Cc - effect8)
    smoothpursuit             <- intercept_smoothpursuit_i + slope_smoothpursuit_i * effect8
    smoothpursuit ~ add(addSd_smoothpursuit)

    # ---- effect9: VAS feeling high ----
    intercept_vashigh_i <- exp(lintercept_vashigh + etalintercept_vashigh)
    slope_vashigh_i     <- slope_vashigh * exp(etaslope_vashigh)
    ke0_vashigh         <- exp(lke0_vashigh + etalke0_vashigh)
    d/dt(effect9)       <- ke0_vashigh * (Cc - effect9)
    vashigh             <- intercept_vashigh_i + slope_vashigh_i * effect9
    vashigh ~ add(addSd_vashigh)

    # ---- effect10: finger tapping (with additive time trend) ----
    # Time t is in minutes; convert the per-day time slope to per-minute (1440 min/day).
    intercept_tapping_i <- exp(lintercept_tapping + etalintercept_tapping)
    slope_tapping_i     <- slope_tapping * exp(etaslope_tapping)
    ke0_tapping         <- exp(lke0_tapping + etalke0_tapping)
    d/dt(effect10)      <- ke0_tapping * (Cc - effect10)
    tapping             <- intercept_tapping_i + slope_tapping_i * effect10 + time_slope_tapping * (t / 1440)
    tapping ~ add(addSd_tapping)
  })
}
