Krause_2023_daridorexant <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral daridorexant (a dual orexin ",
    "receptor antagonist) in 1898 healthy subjects and patients with insomnia ",
    "disorder pooled from 13 phase 1, two phase 2 and two phase 3 studies ",
    "(Krause 2023). First-order absorption after a lag time; Michaelis-Menten ",
    "elimination from the central compartment. Relative bioavailability decreases ",
    "with dose (power, referenced to 50 mg). Food status (fasted or high-fat, ",
    "high-calorie meal versus a light meal / uncontrolled food intake) acts on the ",
    "lag time and ka, and morning versus bedtime dosing on ka. Lean body weight ",
    "scales the central volume, fat mass the peripheral volume and the ",
    "intercompartmental clearance, and age, alkaline phosphatase and lean body ",
    "weight scale Km. Log-normal IIV on all eight structural parameters, ",
    "inter-occasion variability on F, the lag time and ka, and a combined ",
    "additive + proportional residual error."
  )
  reference <- paste(
    "Krause A, Lott D, Brussee JM, Muehlan C, Dingemanse J. (2023).",
    "Population pharmacokinetic modeling of daridorexant, a novel dual orexin",
    "receptor antagonist. CPT: Pharmacometrics & Systems Pharmacology",
    "12(1):74-86. doi:10.1002/psp4.12877.",
    sep = " "
  )
  vignette <- "Krause_2023_daridorexant"

  # Doses are in mg and volumes in L, so central / vc is mg/L (= ug/mL, the unit
  # of Km); the 1000 factor in model() reports Cc in ng/mL, the unit of the
  # paper's concentrations and of the additive residual error.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description = "Administered daridorexant dose on the current dose record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on relative bioavailability, F = 0.41 * (DOSE / 50)^-0.47",
        "(Table 2 row 'Dose', reference dose 50 mg per the Table 2 footnote).",
        "A 200 mg dose has 4^-0.47 = 52% of the bioavailability of a 50 mg dose",
        "(Results; Table 3). Doses studied 5-200 mg. Must be carried on every",
        "record (dose and observation rows) and be > 0. Place the column after",
        "amt in the event table: rxode2 reads a DOSE column that precedes amt",
        "as the dose amount instead of as a covariate."
      ),
      source_name = "Dose"
    ),
    FFM = list(
      description = "Lean body weight (fat-free mass, Janmahasatian formula) at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper calls it 'lean body weight' and derives it with the",
        "Janmahasatian et al. (2005) formula: 9270 * WT / (6680 + 216 * BMI) for",
        "males and 9270 * WT / (8780 + 244 * BMI) for females (Methods). Power",
        "effects on Vc (exponent 0.49) and on Km (exponent -0.20), referenced",
        "to 55 kg (Table 2 footnote). Cohort median 47.2 kg, range 29.1-83.9 kg",
        "(Table 1)."
      ),
      source_name = "LBW"
    ),
    FM = list(
      description = "Fat mass at baseline, body weight minus lean body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "FM = WT - FFM with FFM by the Janmahasatian formula (Methods). Power",
        "effects on Vp (exponent 0.47) and Q (exponent 1.47), referenced to",
        "20 kg (Table 2 footnote). Cohort median 23.4 kg, range 7.0-58.6 kg",
        "(Table 1)."
      ),
      source_name = "FM"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on Km (exponent 0.30), referenced to 30 years per the",
        "Table 2 footnote. The 'typical subject' of the Figure 3 simulations is",
        "60 years old, which is NOT the normalisation reference: only the",
        "30-year reference reproduces the Figure 3 reference-subject AUC0-24,",
        "Cmax, C8h and tmax (see the vignette). Cohort median 54, range 18-88",
        "years (Table 1)."
      ),
      source_name = "Age"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase at baseline",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on Km (exponent 0.34), referenced to 60 U/L (Table 2",
        "footnote). Cohort median 61, range 16-244 U/L (Table 1)."
      ),
      source_name = "ALP"
    ),
    FED = list(
      description = "Not-fasted indicator at dosing (per dose record): 1 = any food (light meal, uncontrolled or high-fat meal), 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (light meal / uncontrolled food intake, together with FED_HIGHFAT = 0)",
      notes = paste(
        "Krause 2023 food status has three levels with light meal / uncontrolled",
        "as the reference: fasted, light meal / uncontrolled, and high-fat,",
        "high-calorie food. Encoded as FED = 0 for fasted, FED = 1 and",
        "FED_HIGHFAT = 0 for the light-meal / uncontrolled reference, and",
        "FED = 1 and FED_HIGHFAT = 1 for the high-fat meal. The fasted effects",
        "therefore enter as exp(beta * (1 - FED)): tlag x exp(-0.51) (15 min",
        "versus 25 min) and ka x exp(0.09) (+9%) (Table 2; Results)."
      ),
      source_name = "Food status"
    ),
    FED_HIGHFAT = list(
      description = "High-fat, high-calorie meal at dosing indicator (per dose record)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (light meal / uncontrolled food intake when FED = 1)",
      notes = paste(
        "1 = drug taken simultaneously with a high-fat, high-calorie meal",
        "(Results; Conclusions). Requires FED = 1. Multiplies tlag by",
        "exp(0.62) (25 to 46 min) and ka by exp(-1.19) (-70%) relative to the",
        "light-meal / uncontrolled reference (Table 2)."
      ),
      source_name = "Food status"
    ),
    DOSETIME_EVENING = list(
      description = "Bedtime dose indicator (per dose record): 1 = bedtime, 0 = morning",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (bedtime administration is the model reference)",
      notes = paste(
        "Krause 2023 contrasts morning administration (mostly the phase 1",
        "studies) with bedtime administration (the phase 2 and 3 insomnia",
        "trials and some phase 1 periods); bedtime is the reference category",
        "(Table 2 footnote, 'evening for time of administration'). The paper",
        "gives no clock-time window. Morning dosing multiplies ka by exp(1.05)",
        "= 2.86 (+186%, Results), entered as exp(1.05 * (1 - DOSETIME_EVENING))."
      ),
      source_name = "Time of administration"
    ),
    OCC = list(
      description = "Occasion index (1-4) for the inter-occasion variability on F, tlag and ka",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Table 2 reports one IOV standard deviation each for F, tlag and ka",
        "(Monolix occasion structure); the paper does not state how many",
        "occasions a subject contributes or how occasions were defined. The",
        "model carries four occasion slots with equal variances (occasions 2-4",
        "fixed to the occasion-1 value, the '$OMEGA BLOCK(1) SAME' idiom).",
        "Records with OCC outside 1-4 get no IOV. A typical-value or",
        "single-occasion simulation may use OCC = 1 throughout."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Primary body-size descriptor tested; replaced by lean body weight and fat mass, which described the data better (Results)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; lean body weight and fat mass were retained instead (Results). Enters the model only through the Janmahasatian FFM formula."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened; sex provides no information beyond body composition (Results, Figure 2)."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault, mL/min)",
      units = "mL/min",
      type = "continuous",
      notes = "Not statistically significant at the 1% level on F, absorption or elimination (Results)."
    ),
    RENALIMP_SEV = list(
      description = "Severe renal impairment indicator (dedicated study)",
      units = "(binary)",
      type = "binary",
      notes = "Tested as renal impairment yes/no; not significant (Results). 7 subjects (Table 1)."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (dedicated study)",
      units = "(binary)",
      type = "binary",
      notes = "Tested pooled with moderate impairment as hepatic impairment yes/no; not significant (Results). 8 subjects (Table 1)."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (dedicated study)",
      units = "(binary)",
      type = "binary",
      notes = "Tested pooled with mild impairment as hepatic impairment yes/no; not significant (Results). 8 subjects (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Not significant at the 1% level (Results). 175 subjects (Table 1)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "daridorexant",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "daridorexant",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "daridorexant",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1898,
    n_studies = 17,
    n_observations = 12627,
    age_range = "18-88 years",
    age_median = "54 years",
    weight_range = "42.0-119.5 kg",
    weight_median = "73.6 kg",
    bmi_range = "17.6-39.6 kg/m^2 (median 25.4)",
    lean_body_weight_range = "29.1-83.9 kg (median 47.2)",
    fat_mass_range = "7.0-58.6 kg (median 23.4)",
    sex_female_pct = 60.9,
    race_ethnicity = c(
      White = 86.7,
      `Black or African American` = 9.2,
      Asian = 1.9,
      Japanese = 0.9,
      `American Indian or Alaska Native` = 0.2,
      `Native Hawaiian or Other Pacific Islander` = 0.2,
      `Other / not reported` = 0.8
    ),
    disease_state = paste(
      "Healthy subjects and subjects with hepatic or renal impairment (phase 1)",
      "and patients with insomnia disorder (phases 2 and 3)"
    ),
    dose_range = "5-200 mg orally (clinical dose 25-50 mg at bedtime)",
    hepatic_function = "Healthy 1882, mild impairment 8, moderate impairment 8; ALP median 61 U/L (range 16-244)",
    renal_function = "Healthy 1891, severe impairment 7; creatinine clearance median 98 mL/min (range 16-241)",
    regions = "Not reported",
    notes = paste(
      "Demographics from Krause 2023 Table 1 (n = 1898). Phase 1: 9420",
      "concentrations from 412 subjects with dense sampling; phases 2 and 3:",
      "3207 next-morning concentrations from 1486 insomnia patients (Results,",
      "Data characteristics). Race and sex percentages are computed from the",
      "Table 1 counts."
    )
  )

  ini({
    # Final model, Krause 2023 Table 2 'All data, final model' column. Parameters
    # without an RSE in that column were held fixed at the 'All Phase I data'
    # estimates (Methods step 4; Table 2 note), hence fixed() below.
    # Reference subject (Table 2 footnote): age 30 years, fat mass 20 kg, lean
    # body weight 55 kg, ALP 60 U/L, dose 50 mg, light meal / uncontrolled food,
    # bedtime (evening) dosing.

    # --- Bioavailability and absorption ---------------------------------------
    lfdepot <- fixed(log(0.41)); label("Relative bioavailability F at the 50 mg reference dose (fraction)") # Table 2 'F' = 0.41 (no RSE in final column)
    e_dose_fdepot <- fixed(-0.47); label("Power exponent of dose on F, referenced to 50 mg (unitless)") # Table 2 'Dose' = -0.47 (no RSE in final column)
    ltlag <- fixed(log(0.41)); label("Absorption lag time, light meal / uncontrolled (h)") # Table 2 't lag (h)' = 0.41 (no RSE in final column)
    e_fasted_tlag <- fixed(-0.51); label("Fasted effect on log tlag (unitless)") # Table 2 'Food: fasted on t lag' = -0.51
    e_fed_highfat_tlag <- fixed(0.62); label("High-fat, high-calorie meal effect on log tlag (unitless)") # Table 2 'Food: fed on t lag' = 0.62
    lka <- fixed(log(1.05)); label("Absorption rate constant ka, light meal / uncontrolled, bedtime (1/h)") # Table 2 'k a (1/h)' = 1.05 (no RSE in final column)
    e_fasted_ka <- fixed(0.09); label("Fasted effect on log ka (unitless)") # Table 2 'Food: fasted on k a' = 0.09
    e_fed_highfat_ka <- fixed(-1.19); label("High-fat, high-calorie meal effect on log ka (unitless)") # Table 2 'Food: fed on k a' = -1.19
    e_morning_ka <- fixed(1.05); label("Morning administration effect on log ka (unitless)") # Table 2 'Morning administration on k a' = 1.05

    # --- Distribution -----------------------------------------------------------
    lvc <- fixed(log(14.6)); label("Central volume Vc (L)") # Table 2 'V c (L)' = 14.60 (no RSE in final column)
    e_ffm_vc <- 0.49; label("Power exponent of lean body weight on Vc, referenced to 55 kg (unitless)") # Table 2 'Lean body weight on V c' = 0.49 (RSE 13.86%)
    lvp <- fixed(log(13.7)); label("Peripheral volume Vp (L)") # Table 2 'V p (L)' = 13.70 (no RSE in final column)
    e_fm_vp <- 0.47; label("Power exponent of fat mass on Vp, referenced to 20 kg (unitless)") # Table 2 'Fat mass on V p' = 0.47 (RSE 13.83%)
    lq <- fixed(log(3.58)); label("Intercompartmental clearance Q (L/h)") # Table 2 'Q(L/h)' = 3.58 (no RSE in final column)
    e_fm_q <- 1.47; label("Power exponent of fat mass on Q, referenced to 20 kg (unitless)") # Table 2 'Fat mass on Q' = 1.47 (RSE 6.48%)

    # --- Michaelis-Menten elimination ----------------------------------------
    lvmax <- fixed(log(6.94)); label("Maximum elimination rate Vm (mg/h)") # Table 2 'V m (mg/h)' = 6.94 (no RSE in final column)
    lkm <- fixed(log(2.36)); label("Michaelis-Menten constant Km (ug/mL)") # Table 2 'K m (ug/ml)' = 2.36 (no RSE in final column)
    e_age_km <- 0.30; label("Power exponent of age on Km, referenced to 30 years (unitless)") # Table 2 'Age on K m' = 0.30 (RSE 6.87%)
    e_alp_km <- 0.34; label("Power exponent of ALP on Km, referenced to 60 U/L (unitless)") # Table 2 'ALP on K m' = 0.34 (RSE 10.54%)
    e_ffm_km <- -0.20; label("Power exponent of lean body weight on Km, referenced to 55 kg (unitless)") # Table 2 'Lean body weight on K m' = -0.20 (RSE 25.95%)

    # --- Inter-individual variability -------------------------------------------
    # Table 2 reports Monolix omegas as SD(.) of the log-normal random effect, so
    # each variance is SD^2. All held fixed at the phase 1 estimates (Methods
    # step 4; no RSE in the final column).
    etalfdepot ~ fixed(0.1764) # Table 2 'SD(F)' = 0.42 -> 0.42^2
    etaltlag ~ fixed(0.0676) # Table 2 'SD(t lag)' = 0.26 -> 0.26^2
    etalka ~ fixed(0.3481) # Table 2 'SD(k a)' = 0.59 -> 0.59^2
    etalvc ~ fixed(0.0289) # Table 2 'SD(V c)' = 0.17 -> 0.17^2
    etalvp ~ fixed(0.0529) # Table 2 'SD(V p)' = 0.23 -> 0.23^2
    etalq ~ fixed(0.25) # Table 2 'SD(Q)' = 0.50 -> 0.50^2
    etalvmax ~ fixed(0.0081) # Table 2 'SD(V m)' = 0.09 -> 0.09^2
    etalkm ~ fixed(0.1296) # Table 2 'SD(K m)' = 0.36 -> 0.36^2

    # --- Inter-occasion variability (four occasion slots, equal variances) ------
    etaiov_fdepot_1 ~ fixed(0.04) # Table 2 'IOV(F)' = 0.20 -> 0.20^2
    etaiov_fdepot_2 ~ fixed(0.04) # occasion 2, equal to occasion 1
    etaiov_fdepot_3 ~ fixed(0.04) # occasion 3, equal to occasion 1
    etaiov_fdepot_4 ~ fixed(0.04) # occasion 4, equal to occasion 1
    etaiov_tlag_1 ~ fixed(0.1764) # Table 2 'IOV(t lag)' = 0.42 -> 0.42^2
    etaiov_tlag_2 ~ fixed(0.1764) # occasion 2, equal to occasion 1
    etaiov_tlag_3 ~ fixed(0.1764) # occasion 3, equal to occasion 1
    etaiov_tlag_4 ~ fixed(0.1764) # occasion 4, equal to occasion 1
    etaiov_ka_1 ~ fixed(0.4356) # Table 2 'IOV(k a)' = 0.66 -> 0.66^2
    etaiov_ka_2 ~ fixed(0.4356) # occasion 2, equal to occasion 1
    etaiov_ka_3 ~ fixed(0.4356) # occasion 3, equal to occasion 1
    etaiov_ka_4 ~ fixed(0.4356) # occasion 4, equal to occasion 1

    # --- Residual error ---------------------------------------------------------
    addSd <- 17.59; label("Additive residual SD (ng/mL)") # Table 2 'Additive error' = 17.59 (RSE 1.97%)
    propSd <- 0.16; label("Proportional residual SD (fraction)") # Table 2 'Multiplicative error' = 0.16 (RSE 1.34%)
  })

  model({
    # 1. Occasion indicators for the inter-occasion variability ----------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 + oc3 * etaiov_tlag_3 + oc4 * etaiov_tlag_4
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4

    # 2. Individual parameters (Methods: continuous covariates as
    #    theta * (cov / ref)^beta, categorical as theta * exp(beta)) ------------
    # Food status: fasted is FED = 0; the light-meal / uncontrolled reference is
    # FED = 1, FED_HIGHFAT = 0; a high-fat meal is FED = 1, FED_HIGHFAT = 1.
    fasted <- 1 - FED
    morning <- 1 - DOSETIME_EVENING

    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) * (DOSE / 50)^e_dose_fdepot
    tlag <- exp(ltlag + etaltlag + iov_tlag + e_fasted_tlag * fasted + e_fed_highfat_tlag * FED_HIGHFAT)
    ka <- exp(lka + etalka + iov_ka + e_fasted_ka * fasted + e_fed_highfat_ka * FED_HIGHFAT + e_morning_ka * morning)
    vc <- exp(lvc + etalvc) * (FFM / 55)^e_ffm_vc
    vp <- exp(lvp + etalvp) * (FM / 20)^e_fm_vp
    q <- exp(lq + etalq) * (FM / 20)^e_fm_q
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm) * (AGE / 30)^e_age_km * (ALP / 60)^e_alp_km * (FFM / 55)^e_ffm_km

    # 3. ODE system --------------------------------------------------------------
    # Michaelis-Menten elimination on the central concentration in mg/L (= ug/mL,
    # the unit of Km); vmax is in mg/h.
    cmgl <- central / vc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - vmax * cmgl / (km + cmgl) - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # 4. Observation and error ---------------------------------------------------
    Cc <- 1000 * cmgl
    # The paper does not name the Monolix error form; the Monolix default
    # combined form 'combined1', y = f + (a + b * f) * e, is used.
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
