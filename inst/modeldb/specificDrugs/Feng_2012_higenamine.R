Feng_2012_higenamine <- function() {
  description <- "Population PK/PD model for intravenous higenamine in 10 healthy Chinese subjects (Feng 2012). Two-compartment disposition with Michaelis-Menten (saturable) elimination from the central compartment, plus a direct-effect Emax sub-model for the cardiovascular-stress heart-rate response (E = E0 + Emax * Cc / (EC50 + Cc)). No demographic covariates were retained in the final model (sex, height, weight, BMI, and age were graphically screened but did not influence PK or PD)."
  reference   <- "Feng S, Jiang J, Hu P, Zhang JY, Liu T, Zhao Q, Li BL (2012). A phase I study on pharmacokinetics and pharmacodynamics of higenamine in healthy Chinese subjects. Acta Pharmacologica Sinica 33(11):1353-1358. doi:10.1038/aps.2012.114"
  vignette    <- "Feng_2012_higenamine"
  units       <- list(time = "min", dosing = "ug", concentration = "ug/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator (1 = female)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Feng 2012 Methods (page 1355) screened sex graphically against the basic-model individual parameters; not retained in the final model. Cohort: 4 male, 6 female (Table 1).",
      source_name        = "sex"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Feng 2012 Methods (page 1355) screened body weight graphically; not retained in the final model. Cohort mean (SD) 60.4 (4.2) kg, range 52.5-66 kg (Table 1).",
      source_name        = "weight"
    ),
    HT = list(
      description        = "Body height",
      units              = "m",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Feng 2012 Methods (page 1355) screened height graphically; not retained in the final model. Cohort mean (SD) 1.60 (0.07) m, range 1.5-1.73 m (Table 1).",
      source_name        = "height"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Feng 2012 Methods (page 1355) screened BMI graphically; not retained in the final model. Cohort mean (SD) 23.3 (0.81) kg/m^2, range 22.1-24.4 (Table 1).",
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Feng 2012 Methods (page 1355) screened age graphically; not retained in the final model. Cohort mean (SD) 30.2 (6.8) years, range 22-41 years (Table 1).",
      source_name        = "age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "22-41 years",
    age_median     = "mean 30.2 (SD 6.8) years",
    weight_range   = "52.5-66 kg",
    weight_median  = "mean 60.4 (SD 4.2) kg",
    sex_female_pct = 60.0,
    race_ethnicity = "Chinese (Han, single-centre cohort at Peking Union Medical College Hospital)",
    disease_state  = "Healthy adult volunteers (clinical-laboratory, ECG, and physical-examination screen normal; resting heart rate 45-90 bpm, systolic BP <= 140 mmHg, diastolic BP <= 90 mmHg).",
    dose_range     = "Continuous IV infusion of higenamine hydrochloride at escalating rates 0.5, 1.0, 2.0, 4.0 ug/kg/min, each rate maintained for 3 min (cumulative dose 22.5 ug/kg over 12 min).",
    regions        = "China (Peking Union Medical College Hospital, Beijing, single centre).",
    notes          = "Per Feng 2012 Table 1. Phase I open-label, non-controlled, single-dose, single-centre study. 4 male and 6 female (40 percent vs 60 percent). Higenamine is an active ingredient of Aconite root being developed as a pharmacologic stress-test agent; the cohort received a single IV stress-test dose with rich plasma (16 timepoints over 102 min postdose) and urinary sampling, plus 10 heart-rate measurements over 42 min postdose. Two non-related adverse events (one mild bilirubin rise, one moderate dizziness and nausea); both transient."
  )

  ini({
    # ---- Structural PK parameters (Feng 2012 Table 4) -------------------
    # Two-compartment disposition with Michaelis-Menten (saturable)
    # elimination from the central compartment. Time unit = min throughout.
    lvc   <- log(18.7);  label("Central volume of distribution (Vc, L)")                     # Table 4: Vc = 18.7 L (estimate %CV 6.6)
    lvp   <- log(43.0);  label("First peripheral volume of distribution (Vp, L)")            # Table 4: Vp = 43.0 L (estimate %CV 3.3)
    lq    <- log(3.8);   label("Inter-compartmental clearance (CLd, L/min)")                 # Table 4: CLd = 3.8 L/min (estimate %CV 5.2)
    lvmax <- log(48.3);  label("Maximum elimination rate (Vmax, ug/min)")                    # Table 4: Vmax = 48.3 (estimate %CV 1.4) -- see Errata in vignette: paper labels Vmax (L/min) but dimensional analysis of the standard Michaelis-Menten elimination rate = Vmax * Cc / (Km + Cc) with Vmax in ug/min reproduces the published Cmax and NCA CL (250 L/h).
    lkm   <- log(3.1);   label("Michaelis-Menten constant (Km, ug/L)")                       # Table 4: Km = 3.1 ug/L (estimate %CV 3.7)

    # ---- Structural PD parameters (Feng 2012 Table 4) -------------------
    # Direct-effect Emax model with baseline for heart rate:
    #   HR = E0 + Emax * Cc / (EC50 + Cc)
    # No effect compartment / no hysteresis -- per Methods page 1355: "a
    # simple direct effect model with baseline".
    le0    <- log(68);   label("Baseline heart rate (E0, bpm)")                              # Table 4: E0 = 68 bpm (estimate %CV 2.0)
    lemax  <- log(73);   label("Maximum drug-induced increase in heart rate (Emax, bpm)")    # Table 4: Emax = 73 bpm (estimate %CV 3.8)
    lec50  <- log(8.1);  label("Higenamine concentration producing 50% of Emax (EC50, ug/L)")# Table 4: EC50 = 8.1 ug/L (estimate %CV 9.1)

    # ---- IIV (Feng 2012 Table 4 'Inter-individual variability (%CV)') --
    # The column is reported in %CV; converted to log-normal omega^2 via
    # omega^2 = log(1 + (CV/100)^2). Parameters with a blank IIV cell in
    # Table 4 (Vp, EC50) are encoded without an eta in the model() block.
    # All IIVs are small (<= 8.7 percent CV) reflecting the rich-sampling
    # 10-subject design with a highly restricted demographic window.
    etalvc   ~ 0.007541    # Table 4: Vc IIV 8.7%   ; log(1 + 0.087^2) = 0.007541
    etalkm   ~ 0.001680    # Table 4: Km IIV 4.1%   ; log(1 + 0.041^2) = 0.001680
    etalvmax ~ 0.000289    # Table 4: Vmax IIV 1.7% ; log(1 + 0.017^2) = 0.000289
    etalq    ~ 0.000100    # Table 4: CLd IIV 1.0%  ; log(1 + 0.010^2) = 0.000100
    etale0   ~ 0.000049    # Table 4: E0 IIV 0.7%   ; log(1 + 0.007^2) = 0.000049
    etalemax ~ 0.000001    # Table 4: Emax IIV 0.1% ; log(1 + 0.001^2) = 0.000001

    # ---- Residual error (Feng 2012 Table 4) ----------------------------
    # Methods page 1355: "A residual error model with a multiplicative
    # component was applied in the final model" => proportional on the
    # linear-space observation.
    propSd    <- 0.260;  label("Proportional residual SD on higenamine plasma Cc (fraction)") # Table 4: Plasma residual error = 0.260 (estimate %CV 6.1)
    propSd_HR <- 0.089;  label("Proportional residual SD on heart rate (fraction)")            # Table 4: Heart rate residual error = 0.089 (estimate %CV 8.1)
  })

  model({
    # ---- Individual PK parameters --------------------------------------
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq + etalq)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm + etalkm)

    # ---- Individual PD parameters --------------------------------------
    e0    <- exp(le0 + etale0)
    emax  <- exp(lemax + etalemax)
    ec50  <- exp(lec50)

    # ---- PK ODE system -------------------------------------------------
    # Two-compartment disposition with Michaelis-Menten elimination from
    # the central compartment. Dosing is IV infusion directly into central
    # (no depot; see Methods 'Methods and drug administration').
    Cc <- central / vc

    d/dt(central)     <- -vmax * Cc / (km + Cc) - q * (Cc - peripheral1 / vp)
    d/dt(peripheral1) <-                          q * (Cc - peripheral1 / vp)

    # ---- PD: direct-effect Emax on heart rate --------------------------
    # E = E0 + Emax * Cc / (EC50 + Cc) per Feng 2012 Methods page 1355.
    HR <- e0 + emax * Cc / (ec50 + Cc)

    # ---- Observations and residual error -------------------------------
    Cc ~ prop(propSd)
    HR ~ prop(propSd_HR)
  })
}
