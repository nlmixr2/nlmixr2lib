Cirincione_2017_exenatide_qt <- function() {
  description <- "Semimechanistic PK/PD model of exenatide effects on plasma glucose, heart rate and the QT interval in healthy adults on a stepped continuous IV infusion (Cirincione 2017 J Clin Pharmacol): one-compartment exenatide PK with parallel linear and Michaelis-Menten elimination; an indirect-response glucose model fed by first-order meal-carbohydrate absorption whose disposition rate is stimulated by exenatide (Emax), with a study-day-specific meal absorption rate; an indirect-response heart-rate model whose production carries a two-harmonic circadian rhythm and is stimulated by the glucose excursion (linear) and by exenatide (Emax); and an uncorrected QT interval given by a power function of the RR interval plus a delayed glucose-transduction term. The glucose-HR-QT parameters are fixed from the same group's baseline-day analysis (Supplemental Table S1)."
  reference <- paste(
    "Cirincione B, LaCreta F, Sager P, Mager DE.",
    "Model-Based Evaluation of Exenatide Effects on the QT Interval in Healthy Subjects",
    "Following Continuous IV Infusion.",
    "J Clin Pharmacol. 2017;57(8):956-965. doi:10.1002/jcph.882.",
    "Fixed glucose-heart rate-QT parameters (Supplemental Table S1) from:",
    "Cirincione B, Sager PT, Mager DE. Influence of Meals and Glycemic Changes on",
    "QT Interval Dynamics. J Clin Pharmacol. 2017;57(8):966-976. doi:10.1002/jcph.933."
  )
  vignette <- "Cirincione_2017_exenatide_qt"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "pg/mL",
    meal = "mg carbohydrate (dosed into gut)",
    glucose = "mg/dL",
    hr = "beats/min",
    qt = "ms"
  )

  paper_specific_compartments <- c("transduction", "qt")

  covariateData <- list(
    OCC = list(
      description = "Study day (occasion) selecting the day-specific first-order meal-glucose absorption rate constant",
      units = "(count)",
      type = "categorical",
      reference_category = "1 (baseline / placebo day)",
      notes = "1 = 24-h baseline day without exenatide; 2, 3, 4 = the first, second and third days of the stepped IV infusion (target 200, 300 and 500 pg/mL). Table 1 estimates ka_glu separately for each study day (days 1-4); values outside 1-4 give a zero absorption rate. Time-varying; a study day runs 8 pm to 8 pm.",
      source_name = "study day"
    )
  )

  compartmentData <- list(
    central = list(analyte = "exenatide", units = "ug", specimen = "plasma", verified = TRUE),
    gut = list(
      analyte = "meal carbohydrate (glucose equivalent)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    glucose = list(analyte = "glucose", units = "mg/dL", specimen = "plasma", verified = TRUE),
    hr = list(analyte = "heart rate", units = "beats/min", specimen = "not applicable", verified = TRUE),
    transduction = list(
      analyte = "delayed glucose-deviation signal",
      units = "mg/dL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "Healthy adults (mean 40.0 +/- 12.4 years in the 20 baseline-day subjects of the companion analysis)",
    weight_range = "Mean 91.2 +/- 12.1 kg in the 20 baseline-day subjects of the companion analysis",
    sex_female_pct = 5,
    race_ethnicity = "Primarily white (85%) in the baseline-day cohort of the companion analysis",
    disease_state = "Healthy volunteers",
    dose_range = "Stepped continuous IV exenatide infusion, 1.0 to 6.3 ug/h, targeting steady-state 200, 300 and 500 pg/mL on infusion days 1-3 (the 700 pg/mL day 4 was not analysed). Cohort A: 1.0 (0-6 h), 1.8 (6-30 h), 2.7 (30-48 h), 3.6 (48-54 h), 4.5 ug/h (54-72 h); cohorts B and C: 1.0, 1.8, 2.5, 3.2, 3.9 ug/h over the same intervals (Figure 1).",
    regions = "United States (single inpatient clinic)",
    notes = "21 healthy subjects were enrolled in an open-label phase 1 pilot thorough-QT design study (cohort A n = 6, B n = 6, C n = 9); 18 had exenatide exposure and glycemic measurements beyond the 24-h baseline period. Standardised 2500 kcal/day diet: breakfast 7 am (cohort A) or 6 am (cohorts B and C), lunch 2 pm, dinner 8 pm; carbohydrate content and percent consumed in Figure 1. Plasma glucose, insulin and exenatide at 14 times per day; HR and QT from continuous Holter ECG (26 baseline extractions, 14 per treatment day). Demographics are not tabulated in this paper; the companion baseline-day analysis (Cirincione, Sager, Mager 2017, 20 evaluable subjects, 1 woman) reports them."
  )

  ini({
    # Exenatide PK (Table 1; Equation 1). One compartment, IV infusion,
    # linear (kel) plus Michaelis-Menten (Vmax, Km in amount units) elimination.
    lkel <- log(0.424)
    label("First-order linear elimination rate constant kel (1/h)") # Table 1: Kel = 0.424 1/h (RSE 7.57%)
    lkm <- log(5.3)
    label("Amount at half-maximal nonlinear elimination Km (ug)") # Table 1: Km = 5.3 ug (RSE 29.2%); Results: Km/Vc = 453 pg/mL
    lvmax <- log(1.79)
    label("Maximal nonlinear elimination rate Vmax (ug/h)") # Table 1: Vmax = 1.79 ug/h (RSE 9.66%)
    lvc <- log(11.7)
    label("Volume of distribution V (L)") # Table 1: V = 11.7 L (RSE 7.17%)

    # Glucose indirect-response model. System parameters fixed from the
    # baseline-day analysis (Supplemental Table S1, 24-hour window).
    lvc_glucose <- fixed(log(653))
    label("Volume of distribution of glucose Vglu (dL)") # Supplemental Table S1: Vg = 653 dL
    lkout_glucose <- fixed(log(1.94))
    label("First-order glucose disposition rate constant kout_glu (1/h)") # Supplemental Table S1: kout_glu = 1.94 1/h
    lrbase_glucose <- fixed(log(92.8))
    label("Steady-state (fasting) plasma glucose Gluss (mg/dL)") # Supplemental Table S1: Gluss = 92.8 mg/dL
    lka_glucose <- fixed(log(1.72))
    label("First-order meal-glucose absorption rate constant ka_glu, study day 1 / baseline (1/h)") # Table 1: ka_glu day 1 = 1.72 1/h (NE); Supplemental Table S1: 1.72 (NE)
    lka_glucose_occ2 <- log(0.162)
    label("First-order meal-glucose absorption rate constant ka_glu, study day 2 / 200 pg/mL target (1/h)") # Table 1: ka_glu day 2 = 0.162 1/h (RSE 29.8%)
    lka_glucose_occ3 <- log(0.121)
    label("First-order meal-glucose absorption rate constant ka_glu, study day 3 / 300 pg/mL target (1/h)") # Table 1: ka_glu day 3 = 0.121 1/h (RSE 41.4%)
    lka_glucose_occ4 <- log(0.194)
    label("First-order meal-glucose absorption rate constant ka_glu, study day 4 / 500 pg/mL target (1/h)") # Table 1: ka_glu day 4 = 0.194 1/h (RSE 43.8%)
    lemax_glucose <- log(0.77)
    label("Maximal fold stimulation of glucose disposition by exenatide Smax_glu/exen (unitless)") # Table 1: Smaxglu/exen = 0.77 (RSE 22.7%)
    lec50_glucose <- fixed(log(90))
    label("Exenatide concentration at half-maximal stimulation of glucose disposition SC50_glu/exen (pg/mL)") # Table 1: SC50glu/exen = 90 (NE); PK/PD Modeling: fixed to 90 pg/mL

    # Heart-rate indirect-response model. System parameters fixed from
    # Supplemental Table S1; exenatide effect from Table 1.
    lkout_hr <- fixed(log(14.3 / 24))
    label("First-order loss rate constant of heart rate kout_HR (1/h)") # Supplemental Table S1: kout_HR = 14.3 1/day, converted to 1/h
    lrbase_hr <- fixed(log(64.1))
    label("Steady-state (baseline) heart rate HRss (beats/min)") # Supplemental Table S1: HRss = 64.1 bpm
    lslope_glucose_hr <- fixed(log(0.00103))
    label("Linear stimulation of heart-rate production per unit glucose excursion above Gluss, S_HR (dL/mg)") # Supplemental Table S1: S_HR = 0.00103 (mg/dL)^-1
    kin_circ_famp1 <- fixed(0.062)
    label("Fractional amplitude of the 24-h cosine on heart-rate production AMP1 (unitless)") # Supplemental Table S1: AMP1 = 0.062
    tacro1 <- fixed(17.3)
    label("Acrophase of the 24-h cosine on heart-rate production ph1 (h after 8 pm)") # Supplemental Table S1: ph1 = 17.3 h
    kin_circ_famp2 <- fixed(0.039)
    label("Fractional amplitude of the 12-h cosine on heart-rate production AMP2 (unitless)") # Supplemental Table S1: AMP2 = 0.039
    tacro2 <- fixed(10.7)
    label("Acrophase of the 12-h cosine on heart-rate production ph2 (h after 8 pm)") # Supplemental Table S1: ph2 = 10.7 h
    lemax_hr <- log(0.527)
    label("Maximal fold stimulation of heart-rate production by exenatide Smax_HR/exen (unitless)") # Table 1: SmaxHR/exen = 0.527 (RSE 23.7%)
    lec50_hr <- log(584)
    label("Exenatide concentration at half-maximal stimulation of heart-rate production SC50_HR/exen (pg/mL)") # Table 1: SC50HR/exen = 584 pg/mL (RSE 44.3%)

    # QT model (Figure 2 equation; Supplemental Table S1).
    lrbase_qt <- fixed(log(409))
    label("QT interval at an RR interval of 1 s, intercept of the QT-RR relationship QTint (ms)") # Supplemental Table S1: QTint = 409 msec
    e_rr_qt <- fixed(0.425)
    label("Exponent of the RR interval in the QT-RR relationship alpha (unitless)") # Supplemental Table S1: alpha = 0.425
    lktr <- fixed(log(0.12))
    label("First-order rate constant of the glucose signal-transduction delay ktr (1/h)") # Supplemental Table S1: ktr = 0.12 1/h (NE)
    e_transduction_qt <- fixed(-1.27)
    label("Change in QT per unit of the delayed glucose-transduction signal S_TR (ms per mg/dL)") # Supplemental Table S1: S_TR = -1.27 msec/(mg/dL) (NE)

    # IIV. 'IIV (%)' is sqrt(omega^2): only that reading makes the printed
    # cov(IIV V, IIV Km) = 0.23 a valid covariance (correlation 0.90); with
    # omega^2 = log(1 + CV^2) the correlation would be 1.11.
    etalkel ~ 0.0320 # Table 1: Kel IIV 17.9% -> 0.179^2
    etalvc + etalkm ~ c(0.0605, 0.23, 1.0816) # Table 1: V IIV 24.6% -> 0.246^2; cov(IIV V, IIV Km) = 0.23; Km IIV 104% -> 1.04^2
    etalemax_glucose ~ 0.8780 # Table 1: Smaxglu/exen IIV 93.7% -> 0.937^2
    etalec50_hr ~ 0.6956 # Table 1: SC50HR/exen IIV 83.4% -> 0.834^2
    etalkout_glucose ~ fixed(0.0784) # Supplemental Table S1: kout_glu IIV 28% -> 0.28^2
    etalkout_hr ~ fixed(0.5730) # Supplemental Table S1: kout_HR IIV 75.7% -> 0.757^2
    etalrbase_hr ~ fixed(0.008464) # Supplemental Table S1: HRss IIV 9.2% -> 0.092^2
    etalslope_glucose_hr ~ fixed(1.9712) # Supplemental Table S1: S_HR IIV 140.4% -> 1.404^2
    etalrbase_qt ~ fixed(0.000767) # Supplemental Table S1: QTint IIV 2.77% -> 0.0277^2

    # Residual error: proportional for every endpoint (Methods; Table 1).
    propSd <- 0.0903
    label("Proportional residual SD of exenatide concentration (fraction)") # Table 1: RVexen = 9.03 %CV (RSE 12.3%)
    propSd_glucose <- 0.104
    label("Proportional residual SD of plasma glucose (fraction)") # Table 1: RVGlu = 10.4 %CV (RSE 15.7%)
    propSd_hr <- 0.0782
    label("Proportional residual SD of heart rate (fraction)") # Table 1: RVHR = 7.82 %CV (RSE 9.55%)
    propSd_qt <- 0.029
    label("Proportional residual SD of the uncorrected QT interval (fraction)") # Table 1: RVQT = 2.9 %CV (RSE 9.64%)
  })

  model({
    # Individual parameters
    kel <- exp(lkel + etalkel)
    km <- exp(lkm + etalkm)
    vmax <- exp(lvmax)
    vc <- exp(lvc + etalvc)

    vc_glucose <- exp(lvc_glucose)
    kout_glucose <- exp(lkout_glucose + etalkout_glucose)
    rbase_glucose <- exp(lrbase_glucose)
    emax_glucose <- exp(lemax_glucose + etalemax_glucose)
    ec50_glucose <- exp(lec50_glucose)
    # Study-day-specific meal absorption rate (Table 1, ka_glu day 1-4)
    ka_glucose <- exp(lka_glucose) * (OCC == 1) +
      exp(lka_glucose_occ2) * (OCC == 2) +
      exp(lka_glucose_occ3) * (OCC == 3) +
      exp(lka_glucose_occ4) * (OCC == 4)

    kout_hr <- exp(lkout_hr + etalkout_hr)
    rbase_hr <- exp(lrbase_hr + etalrbase_hr)
    slope_glucose_hr <- exp(lslope_glucose_hr + etalslope_glucose_hr)
    emax_hr <- exp(lemax_hr)
    ec50_hr <- exp(lec50_hr + etalec50_hr)

    rbase_qt <- exp(lrbase_qt + etalrbase_qt)
    ktr <- exp(lktr)

    # Exenatide plasma concentration: ug / L * 1000 = pg/mL
    Cc <- 1000 * central / vc

    # Drug effects (Equations 2 and 3)
    stim_glucose <- 1 + emax_glucose * Cc / (ec50_glucose + Cc)
    stim_hr <- 1 + emax_hr * Cc / (ec50_hr + Cc)

    # Circadian rhythm on heart-rate production (Equation 4); t = 0 is 8 pm,
    # the start of the baseline day, as in the source analysis.
    circ <- kin_circ_famp1 * cos(2 * pi * (t - tacro1) / 24) +
      kin_circ_famp2 * cos(2 * pi * (t - tacro2) / 12)

    kin_hr <- kout_hr * rbase_hr

    # Equation 1
    d/dt(central) <- -kel * central - vmax * central / (km + central)
    # Meal carbohydrate (mg) absorbed first-order into the glucose space
    d/dt(gut) <- -ka_glucose * gut
    # Equation 2 with Kin_glu = kout_glu * stimulation * Gluss, so exenatide
    # speeds the return of a meal excursion to Gluss and leaves the fasting
    # level at Gluss (the paper's own VPC excludes a fasting fall; see vignette).
    d/dt(glucose) <- ka_glucose * gut / vc_glucose -
      kout_glucose * stim_glucose * (glucose - rbase_glucose)
    # Equation 3 with the glucose stimulation S_HR shown in Figure 2 and in the
    # baseline-day model
    d/dt(hr) <- kin_hr * (1 + slope_glucose_hr * (glucose - rbase_glucose)) *
      (1 + circ) * stim_hr - kout_hr * hr
    # Glucose signal-transduction delay on QT (single transit compartment)
    d/dt(transduction) <- ktr * (glucose - rbase_glucose - transduction)

    glucose(0) <- rbase_glucose
    hr(0) <- rbase_hr

    # QT model (Figure 2): QT = QTint * (60/HR)^alpha + S_TR * TR
    rr_s <- 60 / hr
    qt <- rbase_qt * rr_s^e_rr_qt + e_transduction_qt * transduction
    # Model-predicted heart-rate-corrected QT, mQTc = QT / RR^alpha (Figure 6)
    QTc <- qt / rr_s^e_rr_qt

    Cc ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    hr ~ prop(propSd_hr)
    qt ~ prop(propSd_qt)
  })
}
