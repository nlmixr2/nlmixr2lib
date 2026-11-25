Cirincione_2017_mealQT_24h <- function() {
  description <- paste(
    "Semimechanistic model of meal-driven plasma glucose, heart rate and the",
    "uncorrected QT interval over a full 24-hour baseline day in healthy",
    "adults, with a two-harmonic circadian rhythm on heart-rate production",
    "(Cirincione 2017, 24-hour model). Meal carbohydrate (mg) is absorbed",
    "first-order from a gut depot into a one-compartment glucose pool with",
    "zero-order endogenous production at steady state. The change in glucose",
    "from its steady state linearly stimulates the zero-order production of",
    "an indirect-response heart-rate state, which is further multiplied by",
    "(1 + CIRC), CIRC being the sum of a 24-h and a 12-h cosine. QT follows",
    "a direct power-law QT-RR relationship on the predicted RR interval plus",
    "a delayed glucose effect through one signal-transduction compartment",
    "that shortens QT. Model time 0 is 8 pm, the start of the baseline day."
  )
  reference <- paste(
    "Cirincione B, Sager PT, Mager DE. Influence of Meals and Glycemic",
    "Changes on QT Interval Dynamics. J Clin Pharmacol. 2017;57(8):966-976.",
    "doi:10.1002/jcph.933"
  )
  vignette <- "Cirincione_2017_mealQT"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/dL"
  )

  compartmentData <- list(
    depot = list(
      analyte = "meal carbohydrate (apparent glucose dose)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    glucose = list(
      analyte = "glucose",
      units = "mg/dL",
      specimen = "plasma",
      verified = TRUE
    ),
    hr = list(
      analyte = "heart rate",
      units = "beats/min",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect = list(
      analyte = "glucose-driven QT signal-transduction state",
      units = "mg/dL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "mean 40.0 (SD 12.4) years",
    weight_range = "mean 91.2 (SD 12.1) kg",
    sex_female_pct = 5,
    race_ethnicity = c(White = 85, Other = 15),
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Standardized 2500 kcal/day diet. Breakfast carbohydrate 81, 86 and 99 g",
      "(cohorts A, B, C) at 7 am (A) or 6 am (B, C); lunch carbohydrate 190,",
      "144 and 104 g at 2 pm. Carbohydrate amounts were scaled by the",
      "recorded fraction of the meal eaten (0/25/50/75/100%)."
    ),
    regions = "United States (Covance Clinical Research Unit, Evansville, Indiana)",
    notes = paste(
      "Open-label phase 1 pilot study run to design a thorough QT study of an",
      "antidiabetic agent; only the drug-free 24-hour baseline day (8 pm to",
      "8 pm) was analysed. 21 subjects enrolled in three cohorts (A n = 6,",
      "B n = 6, C n = 9); one was removed for an abnormal postprandial",
      "glucose profile. Demographics from Cirincione 2017 Results 'Baseline",
      "Subject PD and Demographics'. Holter ECGs hourly (every 30 min around",
      "breakfast), 26 HR and QT measurements; plasma glucose at 14 times."
    )
  )

  ini({
    # Glucose sub-model (Eq. 1). Table 1, '24-Hour Window' column.
    lvg <- log(653)
    label("Glucose volume of distribution Vg (dL)") # Table 1 Vg = 653 dL (RSE 22.5%)
    lka_glucose <- fixed(log(1.72))
    label("First-order absorption rate of meal glucose ka_glu (1/h)") # Table 1 ka_glu = 1.72 1/h, RSE 'NE' (fixed; Methods cites ref 20)
    lkout_glucose <- log(1.94)
    label("First-order glucose elimination rate constant kout_glu (1/h)") # Table 1 kout_glu = 1.94 1/h (RSE 22.3%)
    lrbase_glucose <- log(92.8)
    label("Steady-state plasma glucose Gluss (mg/dL)") # Table 1 Gluss = 92.8 mg/dL (RSE 1.35%)

    # Heart-rate indirect response (Eqs. 3 and 9)
    lkout_hr <- log(14.3)
    label("First-order loss rate constant of heart rate kout_HR (1/day)") # Table 1 kout_HR = 14.3 1/days (RSE 26.2%); the Discussion's 14.8 is a typo
    lrbase_hr <- log(64.1)
    label("Steady-state (baseline) heart rate HRss (beats/min)") # Table 1 HRss = 64.1 bpm (RSE 2.2%)
    lsstim_glucose_hr <- log(0.00103)
    label("Linear stimulation of heart-rate production by the change in glucose SHR (dL/mg)") # Table 1 SHR = 0.00103 ([mg/dL]-1) (RSE 80%)

    # Circadian rhythm on heart-rate production (Eq. 10, two harmonics)
    kin_circ_famp1 <- 0.062
    label("Fractional amplitude of the 24-h cosine on heart-rate production AMP1 (unitless)") # Table 1 AMP1 = 0.062 (RSE 13.6%)
    tacro1 <- 17.3
    label("Acrophase of the 24-h cosine ph1 (h after 8 pm)") # Table 1 ph1 = 17.3 h (RSE 3.64%)
    kin_circ_famp2 <- 0.039
    label("Fractional amplitude of the 12-h cosine on heart-rate production AMP2 (unitless)") # Table 1 AMP2 = 0.039 (RSE 21.3%)
    tacro2 <- 10.7
    label("Acrophase of the 12-h cosine ph2 (h after 8 pm)") # Table 1 ph2 = 10.7 h (RSE 5.51%)

    # QT-RR relationship and delayed glucose effect on QT (Eqs. 5, 7, 8b)
    lqtb <- log(409)
    label("Intercept of the QT-RR relationship QTint (ms at RR = 1 s)") # Table 1 QTint = 409 ms (RSE 0.7%)
    lbrr <- log(0.425)
    label("QT-RR correction exponent alpha (unitless)") # Table 1 alpha = 0.425 (RSE 7.76%)
    slope_effect_qt <- fixed(-1.27)
    label("Slope of QT on the glucose signal-transduction state STR (ms per mg/dL)") # Table 1 STR = -1.27, RSE 'NE'; fixed to the breakfast-window estimate (Results 'PD Modeling: 24-Hour Model')
    lktr <- fixed(log(0.12))
    label("First-order rate constant of the glucose-to-QT signal delay ktr (1/h)") # Table 1 ktr = 0.12 1/h, RSE 'NE'; fixed to the breakfast-window estimate

    # IIV: log-normal. Table 1 'IIV (%)' is 100 * sqrt(omega^2), so
    # omega^2 = (IIV / 100)^2. The log(1 + CV^2) reading gives a mean
    # postprandial HR rise of 5.5 bpm, outside the published 6.0-7.6 bpm
    # 95% CI of the 500-subject simulation (Results); this reading gives 6.9.
    etalkout_glucose ~ 0.0784 # Table 1 kout_glu IIV 28% -> 0.28^2
    etalkout_hr ~ 0.573049 # Table 1 kout_HR IIV 75.7% -> 0.757^2
    etalrbase_hr ~ 0.008464 # Table 1 HRss IIV 9.2% -> 0.092^2
    etalsstim_glucose_hr ~ 1.971216 # Table 1 SHR IIV 140.4% -> 1.404^2
    etalqtb ~ 0.00076729 # Table 1 QTint IIV 2.77% -> 0.0277^2

    # Residual error: proportional on each endpoint (Methods 'Interindividual
    # and Residual Error Models').
    propSd_glucose <- 0.113
    label("Proportional residual error on glucose (fraction)") # Table 1 second 'RV glucose' row = 11.3 %CV (final-model residual block); first-stage row = 11%
    propSd_hr <- 0.0796
    label("Proportional residual error on heart rate (fraction)") # Table 1 RV HR = 7.96 %CV
    propSd_qt <- 0.0301
    label("Proportional residual error on QT (fraction)") # Table 1 RV QT = 3.01 %CV
  })

  model({
    vg <- exp(lvg)
    ka_glucose <- exp(lka_glucose)
    kout_glucose <- exp(lkout_glucose + etalkout_glucose)
    rbase_glucose <- exp(lrbase_glucose)
    # Steady-state production: Kin_glu = kout_glu * Gluss
    kin_glucose <- kout_glucose * rbase_glucose

    # kout_HR is reported per day; the model runs in hours
    kout_hr <- exp(lkout_hr + etalkout_hr) / 24
    rbase_hr <- exp(lrbase_hr + etalrbase_hr)
    kin_hr <- kout_hr * rbase_hr
    sstim_glucose_hr <- exp(lsstim_glucose_hr + etalsstim_glucose_hr)

    qtb <- exp(lqtb + etalqtb)
    brr <- exp(lbrr)
    ktr <- exp(lktr)

    # Eq. 10 with two harmonics; t is hours after 8 pm (the start of the
    # baseline day), so ph1 = 17.3 h is 1:18 pm and ph2 = 10.7 h is 6:42 am
    # and 6:42 pm.
    circ <- kin_circ_famp1 * cos(2 * pi * (t - tacro1) / 24) +
      kin_circ_famp2 * cos(2 * pi * (t - tacro2) / 12)
    circ0 <- kin_circ_famp1 * cos(2 * pi * (0 - tacro1) / 24) +
      kin_circ_famp2 * cos(2 * pi * (0 - tacro2) / 12)

    glucose(0) <- rbase_glucose
    hr(0) <- rbase_hr * (1 + circ0)

    d/dt(depot) <- -ka_glucose * depot
    d/dt(glucose) <- ka_glucose * depot / vg + kin_glucose - kout_glucose * glucose
    d/dt(hr) <- kin_hr * (1 + sstim_glucose_hr * (glucose - rbase_glucose)) * (1 + circ) - kout_hr * hr
    d/dt(effect) <- ktr * (glucose - rbase_glucose - effect)

    # RR interval in seconds (Eq. 5)
    rr <- 60 / hr
    # Eq. 8b: uncorrected QT
    qt <- qtb * rr^brr + slope_effect_qt * effect
    # Eq. 6: heart-rate-corrected QT (mQTc)
    qtc_hr <- qt / rr^brr
    # Heart-rate- and glucose-corrected QT (mQTcG); Eq. 8a is printed with
    # '+ STR * TR', but only the '-' form removes the glucose term of Eq. 8b
    # and gives the flat mQTcG profile of Figure 5D.
    qtc_hr_glucose <- (qt - slope_effect_qt * effect) / rr^brr

    glucose ~ prop(propSd_glucose)
    hr ~ prop(propSd_hr)
    qt ~ prop(propSd_qt)
  })
}
