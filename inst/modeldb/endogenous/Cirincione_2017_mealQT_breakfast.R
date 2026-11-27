Cirincione_2017_mealQT_breakfast <- function() {
  description <- paste(
    "Semimechanistic model of meal-driven plasma glucose, heart rate and the",
    "uncorrected QT interval around breakfast in healthy adults (Cirincione",
    "2017, breakfast-window model, 6 am to 2 pm). Meal carbohydrate (mg) is",
    "absorbed first-order from a gut depot into a one-compartment glucose",
    "pool with zero-order endogenous production at steady state. The change",
    "in glucose from its steady state linearly stimulates the zero-order",
    "production of an indirect-response heart-rate state. QT follows a",
    "direct power-law QT-RR relationship on the predicted RR interval plus a",
    "delayed glucose effect through one signal-transduction compartment that",
    "shortens QT. No circadian term; model time 0 is 6 am."
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
      "Standardized breakfast with 81, 86 and 99 g carbohydrate (cohorts A,",
      "B, C) at 7 am (A) or 6 am (B, C), scaled by the recorded fraction of",
      "the meal eaten (0/25/50/75/100%)."
    ),
    regions = "United States (Covance Clinical Research Unit, Evansville, Indiana)",
    notes = paste(
      "Open-label phase 1 pilot study run to design a thorough QT study of an",
      "antidiabetic agent; this model covers the richly sampled breakfast",
      "window (6 am to 2 pm) of the drug-free baseline day. 21 subjects",
      "enrolled in three cohorts (A n = 6, B n = 6, C n = 9); one was removed",
      "for an abnormal postprandial glucose profile. Demographics from",
      "Cirincione 2017 Results 'Baseline Subject PD and Demographics'. Holter",
      "ECGs every 30 min around breakfast."
    )
  )

  ini({
    # Glucose sub-model (Eq. 1). Table 1, 'Breakfast Window' column.
    lvg <- log(781)
    label("Glucose volume of distribution Vg (dL)") # Table 1 Vg = 781 dL (RSE 20.6%)
    lka_glucose <- fixed(log(1.72))
    label("First-order absorption rate of meal glucose ka_glu (1/h)") # Table 1 ka_glu = 1.72 1/h, RSE 'NE'; Methods: fixed to 1.72 h-1 (ref 20)
    lkout_glucose <- log(1.52)
    label("First-order glucose elimination rate constant kout_glu (1/h)") # Table 1 kout_glu = 1.52 1/h (RSE 24.5%)
    lrbase_glucose <- log(90)
    label("Steady-state plasma glucose Gluss (mg/dL)") # Table 1 Gluss = 90 mg/dL (RSE 1.83%)

    # Heart-rate indirect response (Eqs. 2 and 3)
    lkout_hr <- log(7.51)
    label("First-order loss rate constant of heart rate kout_HR (1/day)") # Table 1 kout_HR = 7.51 1/days (RSE 33.5%)
    lrbase_hr <- log(64.3)
    label("Steady-state (baseline) heart rate HRss (beats/min)") # Table 1 HRss = 64.3 bpm (RSE 2.67%)
    lsstim_glucose_hr <- log(0.00407)
    label("Linear stimulation of heart-rate production by the change in glucose SHR (dL/mg)") # Table 1 SHR = 0.00407 ([mg/dL]-1) (RSE 25.5%)

    # QT-RR relationship and delayed glucose effect on QT (Eqs. 4, 5, 7, 8b)
    lqtb <- log(406)
    label("Intercept of the QT-RR relationship QTint (ms at RR = 1 s)") # Table 1 QTint = 406 ms (RSE 0.68%)
    lbrr <- log(0.362)
    label("QT-RR correction exponent alpha (unitless)") # Table 1 alpha = 0.362 (RSE 10.3%)
    slope_effect_qt <- -1.27
    label("Slope of QT on the glucose signal-transduction state STR (ms per mg/dL)") # Table 1 STR = -1.27 (RSE 42.8%)
    lktr <- log(0.12)
    label("First-order rate constant of the glucose-to-QT signal delay ktr (1/h)") # Table 1 ktr = 0.12 1/h (RSE 62.5%); Results: 0.120 h-1

    # IIV: log-normal. Table 1 'IIV (%)' is 100 * sqrt(omega^2), so
    # omega^2 = (IIV / 100)^2 (the reading that reproduces the published
    # 24-hour simulation; see Cirincione_2017_mealQT_24h).
    etalkout_glucose ~ 0.164025 # Table 1 kout_glu IIV 40.5% -> 0.405^2
    etalkout_hr ~ 1.2769 # Table 1 kout_HR IIV 113% -> 1.13^2
    etalrbase_hr ~ 0.00743044 # Table 1 HRss IIV 8.62% -> 0.0862^2
    etalqtb ~ 0.00081225 # Table 1 QTint IIV 2.85% -> 0.0285^2

    # Residual error: proportional on each endpoint (Methods 'Interindividual
    # and Residual Error Models').
    propSd_glucose <- 0.0878
    label("Proportional residual error on glucose (fraction)") # Table 1 second 'RV glucose' row = 8.78 %CV (final-model residual block); first-stage row = 9.18%
    propSd_hr <- 0.07
    label("Proportional residual error on heart rate (fraction)") # Table 1 RV HR = 7 %CV
    propSd_qt <- 0.0248
    label("Proportional residual error on QT (fraction)") # Table 1 RV QT = 2.48 %CV
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
    sstim_glucose_hr <- exp(lsstim_glucose_hr)

    qtb <- exp(lqtb + etalqtb)
    brr <- exp(lbrr)
    ktr <- exp(lktr)

    glucose(0) <- rbase_glucose
    hr(0) <- rbase_hr

    d/dt(depot) <- -ka_glucose * depot
    d/dt(glucose) <- ka_glucose * depot / vg + kin_glucose - kout_glucose * glucose
    d/dt(hr) <- kin_hr * (1 + sstim_glucose_hr * (glucose - rbase_glucose)) - kout_hr * hr
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
