Bokor_2025_procalcitonin <- function() {
  description <- paste(
    "One-compartment turnover model for serum procalcitonin (PCT) in",
    "Japanese preterm and term neonates and infants after birth (Bokor",
    "2025). Birth is a nominal bolus of 1 unit into a first-order",
    "birth-stimulus compartment, released after a lag time; the",
    "stimulus increases PCT production linearly (SYNB = 1 + SLOPEB * CB).",
    "PCT clearance is allometrically scaled to 70 kg (exponent 3/4) and",
    "matures with postmenstrual age (sigmoid Hill function); volume",
    "scales linearly with weight. Typical values are the Scenario 3",
    "(clinical-practice clock-time imputation) estimates the authors",
    "selected as best; the random-effect variances are those of the",
    "Scenario 0 control stream (the only variances the paper reports).",
    sep = " "
  )
  reference <- paste(
    "Bokor AJ, Holford N, Hannam JA.",
    "Imputation of missing clock times - application to procalcitonin",
    "concentration time course after birth.",
    "J Pharmacokinet Pharmacodyn. 2025;52(2):20.",
    "doi:10.1007/s10928-025-09965-8.",
    sep = " "
  )
  vignette <- "Bokor_2025_procalcitonin"
  units <- list(
    time = "h",
    dosing = "unit (nominal birth-event bolus of 1 into the birth-stimulus compartment at the time of birth)",
    concentration = "ug/L"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Allometric scaling to 70 kg: exponent 3/4 on",
        "clearance and 1 on volume (fixed, theory based). The source",
        "control stream uses birth weight (WTB) for the first 5 days of",
        "life and the current (or Sumpter-Holford 2011 estimated) weight",
        "thereafter; supply WT on that basis to reproduce the fit. The",
        "production rate RateInPCT is NOT weight scaled."
      ),
      source_name = "WTB / CWT / EWT (grams in the source data; WT in kg here)"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying; postmenstrual age in WEEKS (the source equation and",
        "TM50PCT are in weeks). Drives the sigmoid maturation of PCT",
        "clearance Fmat = 1 / (1 + (PMA / TM50PCT)^-HillPCT) (Eq. 4)."
      ),
      source_name = "PMAW"
    )
  )

  compartmentData <- list(
    birth = list(
      analyte = "birth-event stimulus (latent, nominal 1-unit bolus at birth)",
      units = "unit",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(
      analyte = "procalcitonin",
      units = "ug",
      specimen = "serum",
      verified = TRUE
    )
  )

  # The birth-stimulus compartment (NONMEM COMP (BIRTH)) is a latent state
  # specific to this model: a unit amount released at birth after a lag that
  # decays first-order and drives PCT production.
  paper_specific_compartments <- c("birth")

  population <- list(
    species = "human",
    n_subjects = 282L,
    n_studies = 1L,
    n_observations = 1275L,
    age_range = "Birth to 19 weeks postnatal age",
    gestational_age_range = "22 to 41 weeks",
    weight_range = "Not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "Japanese",
    disease_state = paste(
      "Preterm and term neonates and infants with serum PCT measured as",
      "part of routine care (dataset of Fukuzumi 2016, Sci Rep 6:1-6,",
      "paper reference 16); the control stream is titled 'healthy",
      "neonates'."
    ),
    dose_range = "n/a (endogenous biomarker; birth event as nominal unit bolus)",
    regions = "Japan",
    notes = paste(
      "Methods, first paragraph. Sampling dates were recorded but clock",
      "times were not; the paper compared clock-time imputation",
      "scenarios (Scenario 0 none; 1A/1B sequence-preserving with 8 or",
      "12 h minimum spacing; 2 rise-then-fall reordering; 3 clinical",
      "practice, first sample at 09:00). Scenario 3 was best (average",
      "dOFV -62.6 vs Scenario 0, Table 1)."
    )
  )

  ini({
    # Typical values: Scenario 3 averages over 100 imputed datasets
    # (Table 2 column '3'; Supplementary Table 5 'Average'). Time in hours.
    lkin <- log(0.00466)
    label("PCT production rate RateInPCT (ug/h; not weight scaled)") # Table 2 Scenario 3 RateInPCT = 0.00466 ug/h (RSE 3.86%)
    lcl <- log(0.493)
    label("PCT clearance CLPCT at 70 kg, fully mature (L/h)") # Table 2 Scenario 3 CLPCT = 0.493 L/h/70kg (RSE 3.66%)
    lvc <- fixed(log(15))
    label("PCT volume of distribution VPCT at 70 kg (L)") # Table 2 VPCT = 15.0 L/70kg FIX (calcitonin volume, paper reference 17)
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of weight on CLPCT (unitless)") # Methods 'Modelling methods': fixed theory-based exponent 3/4; control stream PWRCL=3/4
    e_wt_vc <- fixed(1)
    label("Allometric exponent of weight on VPCT (unitless)") # Methods 'Modelling methods': fixed exponent 1; control stream VC=THETA(3)*FSZ
    ltm50_cl <- log(24.7)
    label("Postmenstrual age at 50% maturation of CLPCT, TM50PCT (weeks)") # Table 2 Scenario 3 TM50PCT = 24.7 weeks (RSE 3.04%)
    lhill_cl <- log(7.82)
    label("Hill exponent of the CLPCT maturation function, HillPCT (unitless)") # Table 2 Scenario 3 HillPCT = 7.82 (RSE 4.88%)
    lthalf_birth <- log(5.67)
    label("Elimination half-life of the birth-effect stimulus, TelB (h)") # Table 2 Scenario 3 TelB = 5.67 h (RSE 4.63%)
    ltlag <- log(12.0)
    label("Lag time between birth and the onset of the birth stimulus, TLAG (h)") # Table 2 Scenario 3 TLAG = 12.0 h (RSE 9.60%)
    lslope_birth <- log(75.4)
    label("Slope of PCT production on birth-stimulus amount, SLOPEB (1/unit)") # Table 2 Scenario 3 SLOPEB = 75.4 (RSE 6.68%)

    # Between-subject variability: $OMEGA BLOCK(5) of the Supplementary Text 1
    # control stream (Scenario 0 final estimates), order RATEIN, CL, SLOPEB,
    # TLAGB, TM50. The paper does not report variances for Scenario 3; the
    # random-effects model was the same in every scenario (Methods
    # 'Parameter estimation'). PPV on TELB was 0 FIX and is omitted.
    etalkin + etalcl + etalslope_birth + etaltlag + etaltm50_cl ~ c(
      0.276,
      0.125, 0.0699,
      0.144, 0.0416, 2.37,
      0.0102, 0.00467, 0.0609, 0.00439,
      0.0156, 0.019, -0.315, -0.00915, 0.112
    ) # Supplementary Text 1 $OMEGA BLOCK(5) PPV_RATEIN, PPV_CL, PPV_SLOPEB, PPV_TLAGB, PPV_TM50

    # Residual error: proportional only (RUV_SD 0 FIX in the control stream).
    propSd <- 0.294
    label("Proportional residual error on PCT (fraction)") # Table 2 Scenario 3 Proportional RUV = 0.294 (RSE 4.37%)
  })

  model({
    # Size and maturation (Eq. 4; control stream $PK FSZ, FMAT)
    fsz <- WT / 70
    tm50_cl <- exp(ltm50_cl + etaltm50_cl)
    hill_cl <- exp(lhill_cl)
    fmat <- 1 / (1 + (PAGE / tm50_cl)^(-hill_cl))

    # Individual parameters
    kin <- exp(lkin + etalkin)
    cl <- exp(lcl + etalcl) * fsz^e_wt_cl * fmat
    vc <- exp(lvc) * fsz^e_wt_vc
    thalf_birth <- exp(lthalf_birth)
    kel_birth <- log(2) / thalf_birth
    slope_birth <- exp(lslope_birth + etalslope_birth)
    tlag <- exp(ltlag + etaltlag)

    # Baseline PCT before the birth stimulus starts (Eq. 3)
    base <- kin / cl

    # Birth stimulus amount in a compartment of unit volume (VB = 1) drives
    # production linearly (Eq. 2); it is zero until the lagged birth bolus.
    synb <- 1 + slope_birth * birth

    d/dt(birth) <- -kel_birth * birth
    alag(birth) <- tlag
    # Eq. 1 in amounts, as coded in the control stream $DES:
    # DADT(2) = RATEIN_PCT * SYNB - DCP * CL
    d/dt(central) <- kin * synb - cl * central / vc
    central(0) <- base * vc

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
