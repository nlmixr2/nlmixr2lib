Shoji_2011_pregabalin <- function() {
  description <- "One-compartment population PK model for pregabalin in adults (Shoji 2011 BJCP; pooled healthy volunteers, subjects with impaired renal function, and patients with post-herpetic neuralgia or diabetic peripheral neuropathy from 14 clinical trials). CL/F is proportional to Cockcroft-Gault creatinine clearance (capped at an estimated break point) with an additional ideal-body-weight power effect. V/F depends on ideal body weight, body mass index, age, and sex. Absorption rate and lag-time are reduced by a high-fat meal at the time of dosing. Combined proportional + additive residual error is stratified by healthy-vs-patient status."
  reference <- "Shoji S, Suzuki M, Tomono Y, Bockbrader HN, Matsui S. Population pharmacokinetics of pregabalin in healthy subjects and patients with post-herpetic neuralgia or diabetic peripheral neuropathy. Br J Clin Pharmacol. 2011;72(1):63-76. doi:10.1111/j.1365-2125.2011.03932.x"
  vignette <- "Shoji_2011_pregabalin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  paper_specific_residual_sds <- c(
    "propSdHealthy", "addSdHealthy",
    "propSdPatient", "addSdPatient"
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance by the Cockcroft-Gault equation. NOT BSA-normalized; raw Cockcroft-Gault mL/min. Capped at the estimated break point th_bp = 107 mL/min when entering CL/F (the paper's saturation knot above which CL/F no longer scales linearly with CLcr).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per visit; covariate values from the earliest visit are carried forward in the source dataset (Methods, Inclusion of covariates paragraph). Cohort range 10.0-230 mL/min, mean 86 mL/min (Table 1, All Total row). CL/F is modelled as a strict linear function of min(CRCL, th_bp); there is no intercept term (Methods, Base model development; the intercept 95% CI -0.00628 to 0.240 included zero so was dropped). Same canonical CRCL form used in Delattre 2010 amikacin (raw Cockcroft-Gault, NOT BSA-normalized).",
      source_name        = "CLcr"
    ),
    IBW = list(
      description        = "Ideal body weight derived from total body weight, height, and sex per the formula in Methods (the Cockcroft-Gault context formula). Time-fixed.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 62 kg (the population-mean IBW called out explicitly in the Discussion: 'even when IBW became extremely low (i.e. half of the mean 62 kg) CL/F deceased by 22%'). Enters CL/F and V/F as a power scalar (IBW / 62)^e_ibw_<param>. The paper's Methods equation for IBW is from reference [12] (Devine-family variant); the per-paper formula should be applied when assembling a virtual cohort if only TBW + HT + SEXF are available.",
      source_name        = "IBW"
    ),
    BMI = list(
      description        = "Body mass index (kg/m^2) from total body weight and height, time-fixed per subject. Enters V/F as a power scalar (BMI / 25)^e_bmi_vc.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 25 kg/m^2, approximately equal to the population-mean BMI (mean TBW 71 kg, mean HT ~167 cm -> BMI ~25.4 kg/m^2). The reference value 25 is the comparator used in the paper's own Discussion sensitivity test ('The V/F decreased to 84% when BMI changed from 25 to 18 kg/m^2'). The two-subject HT imputation to 167 cm (Methods, Demographic data) keeps the cohort-mean BMI ~25.",
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Subject age at study entry (years). Time-fixed. Enters V/F as a power scalar (AGE / 59)^e_age_vc.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 59 years (the population-mean age called out explicitly in the Discussion: 'it decreased only 6.5% from the mean age (59 years) to the maximal observed age (101 years)'). Cohort range 19-101 years (Table 1, All Total row).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed. Shoji 2011 codes SEX as 'male = 0, female = 1' in Figure 2 caption, matching the canonical SEXF orientation directly. Enters V/F as the categorical multiplicative factor e_sexf_vc^SEXF; the female-vs-male V/F ratio is 0.906 (Table 3, q_Gender on V/F final model). Cohort 37.2% female, 62.8% male (Table 2).",
      source_name        = "SEX"
    ),
    FED = list(
      description        = "Per-dose fed-vs-fasted indicator (1 = dose given with food within 2 h of the meal, or food but exact meal time unknown; 0 = fasted).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-dose-event indicator, time-varying within subject across multiple dosing events. Self-reported food status (Methods, Inclusion of covariates paragraph). Enters ka and tlag as multiplicative shifts: ka_fed = ka * (1 + e_food_ka) -> ~93% reduction (NONMEM (1 + theta * FED) parameterisation); tlag_fed = tlag * (1 + e_food_tlag) -> ~81% increase. Consistent with the published high-fat-meal slowing of pregabalin absorption.",
      source_name        = "FOOD"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator (1 = healthy subject from studies HV01-HV09 or impaired-renal-function / elderly substudies HV05 / HV07; 0 = patient with post-herpetic neuralgia or diabetic peripheral neuropathy from studies PT01-PT05).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient cohort)",
      notes              = "Time-fixed per subject. Used in model() to switch the proportional + additive residual error magnitudes between the two paper-defined groups; the healthy-volunteer arms had richly-sampled PK and lower residual variability (CV 22%, SD 0.0239 ug/mL), while patient arms had sparse outpatient sampling and higher residual variability (CV 28.5%, SD 0.236 ug/mL) (Table 3, residual variability rows final model). The patient-type variable was tested as a covariate on CL/F but did not remain in the final model (ratio patient/healthy 92.5% with 95% CI 85.2-99.9% fell within the 80-125% clinically-less-significant range; Discussion paragraph 5).",
      source_name        = "(derived from study identifier; HV studies = 1, PT studies = 0)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 616L,
    n_studies      = 14L,
    age_range      = "19-101 years",
    age_median     = "59.1 years (mean; Table 1, All Total row)",
    weight_range   = "31-142 kg",
    weight_median  = "71.0 kg (mean total body weight; Table 1, All Total row)",
    sex_female_pct = 37.2,
    race_ethnicity = c(White = 52.9, Black = 0.97, Asian = 41.6, Other = 4.55),
    disease_state  = "Pooled healthy-volunteer (n=195, 31.7%), post-herpetic neuralgia (n=267, 43.3%), and painful diabetic peripheral neuropathy (n=154, 25.0%) cohorts drawn from 14 clinical trials. Five phase-1 studies (HV01-HV09) used dense PK sampling in healthy adults including substudies in subjects with impaired renal function (HV05) and elderly subjects (HV07). Four post-herpetic neuralgia trials (PT01-PT04) and one diabetic peripheral neuropathy trial (PT05) used sparse outpatient PK sampling (1-2 samples per patient).",
    dose_range     = "Single doses 1-300 mg oral; multiple-dose 75-300 mg BID or 25-200 mg TID (up to 600 mg/day in efficacy studies, up to 900 mg/day in phase-1). PT05 included a CrCl-based dose reduction to 150 mg BID for patients allocated to 300 mg BID who had 30 <= CrCl < 60 mL/min.",
    regions        = "United States, Japan, European countries, Australia, Canada (per Table 1 study identifiers).",
    bioanalysis    = "Plasma pregabalin quantified by validated HPLC-UV (LLQ 0.005-0.05 ug/mL across studies) or LC-MS/MS (LLQ 0.025 ug/mL; HV08, HV09, PT04, PT05). Precision and accuracy 1.7-5.5% and -4.4 to +4.0% across assays (Methods, Clinical studies and assay methods).",
    notes          = "Total 5275 plasma pregabalin concentrations: 4650 (88.2%) from healthy subjects and 625 (11.8%) from patients. Sparse-sample distribution: 36 observations in the absorption phase (<= 0.75 h), 358 around Cmax (0.75-3.00 h), and 231 in the elimination phase (> 3.0 h). Two subjects with missing height were imputed to 167 cm (the cohort-mean height). Estimation used NONMEM V Level 1.1 with key results confirmed in NONMEM 7 Level 1.2; first-order conditional estimation with eta-epsilon interaction (FOCE-I)."
  )

  ini({
    # ---------------------------------------------------------------
    # Clearance: CL/F is parameterised in the paper as a strict slope on
    # min(CRCL, th_bp) with no intercept (Methods, Base model development:
    # the intercept 95% CI -0.00628 to 0.240 included zero so was dropped).
    # The Table 3 slope (q_CL/F = 0.0462 L/h per mL/min) implies a mean
    # CL/F = 0.0462 * 86 = 3.97 L/h at the population-mean CRCL = 86 mL/min
    # (Results, Base model development: "0.0460 multiplied by 86 = 3.96 L/h").
    # We re-express the same proportional relationship in canonical form
    # by anchoring CL/F at the population-mean CRCL: CL/F = exp(lcl) *
    # (min(CRCL, th_bp) / 86) * (IBW / 62)^e_ibw_cl. The slope on CLcr is
    # implicit in the linear ratio (CRCL / 86). The Table 3 slope is
    # recovered as exp(lcl) / 86 = 3.97 / 86 = 0.0462 (L/h per mL/min).
    # Above the break point th_bp the CL/F effect saturates; the null
    # hypothesis th_bp = 150 was rejected (Discussion paragraph 3).
    lcl      <- log(0.0462 * 86); label("Typical CL/F at reference CRCL=86 mL/min, IBW=62 kg (L/h)")  # Table 3, q_CL/F final model x Table 1 mean CRCL
    th_bp    <- 107;              label("CRCL break point above which CL/F saturates (mL/min)")       # Table 3, q_BP final model
    e_ibw_cl <- 0.354;            label("IBW power exponent on CL/F (unitless)")                       # Table 3, q_IBW on CL/F final model

    # ---------------------------------------------------------------
    # Apparent volume of distribution at the reference covariates:
    #   IBW = 62 kg, BMI = 25 kg/m^2, AGE = 59 years, SEXF = 0 (male).
    # All covariate effects enter V/F as power-of-ratio or categorical-power
    # scalars; the typical V/F at the reference covariates is exp(lvc) = 35.6 L.
    lvc       <- log(35.6); label("Apparent volume of distribution V/F at reference covariates (L)")  # Table 3, q_V/F final model
    e_ibw_vc   <- 0.819;     label("IBW power exponent on V/F (unitless)")                              # Table 3, q_IBW on V/F final model
    e_bmi_vc   <- 0.525;     label("BMI power exponent on V/F (unitless)")                              # Table 3, q_BMI on V/F final model
    e_age_vc   <- -0.125;    label("AGE power exponent on V/F (unitless)")                              # Table 3, q_Age on V/F final model
    e_sexf_vc  <- 0.906;     label("Female-vs-male V/F ratio applied as e_sexf_vc ^ SEXF (unitless)")    # Table 3, q_Gender on V/F final model

    # ---------------------------------------------------------------
    # Absorption: first-order rate constant and absorption lag-time. Both
    # carry a food covariate effect applied multiplicatively as
    # parameter_fed = parameter * (1 + e_food_<param> * FED) -- the
    # NONMEM (1 + theta * FED) parameterisation. At FED = 1 the food
    # multiplier becomes (1 - 0.930) = 0.070 on ka (~93% slower) and
    # (1 + 0.811) = 1.811 on tlag (~81% longer).
    lka         <- log(7.99);  label("Absorption rate constant ka at FED = 0 (1/h)")  # Table 3, q_ka final model
    e_food_ka   <- -0.930;     label("Food effect on ka (multiplicative: ka_fed = ka * (1 + e_food_ka))")  # Table 3, q_Food on ka final model
    ltlag       <- log(0.243); label("Absorption lag-time at FED = 0 (h)")  # Table 3, q_t_lag final model
    e_food_tlag <- 0.811;      label("Food effect on tlag (multiplicative: tlag_fed = tlag * (1 + e_food_tlag))")  # Table 3, q_Food on tlag final model

    # ---------------------------------------------------------------
    # Inter-individual variability: diagonal log-normal IIV on CL/F, V/F,
    # and ka. Methods paragraph 3 specifies "a variance-covariance matrix W
    # with independent correlation structure" -- diagonal omega matrix.
    # Final-model CV% values: 14.8% (CL/F), 8.54% (V/F), 93.2% (ka).
    # omega^2 = log(CV^2 + 1) for log-normal variance.
    etalcl ~ 0.02168   # Table 3, CV% (CL/F) final model 14.8%; log(1 + 0.148^2) = 0.02168
    etalvc ~ 0.00727   # Table 3, CV% (V/F)  final model  8.54%; log(1 + 0.0854^2) = 0.00727
    etalka ~ 0.62526   # Table 3, CV% (ka)   final model 93.2%; log(1 + 0.932^2)  = 0.62526

    # ---------------------------------------------------------------
    # Residual variability: combined proportional + additive error model,
    # Y_ij = C_ij * (1 + eps1) + eps2 (Methods paragraph 3). Magnitudes
    # are stratified between healthy and patient cohorts. Switched in
    # model() via DIS_HEALTHY. All four names are declared in
    # `paper_specific_residual_sds` because the canonical residual-error
    # matcher recognises only bare propSd / addSd / expSd -- not
    # per-cohort suffixes.
    propSdHealthy <- 0.220;  label("Proportional residual SD, healthy subjects (fraction)")  # Table 3, CV% (healthy) final model 22.0%
    addSdHealthy  <- 0.0239; label("Additive residual SD, healthy subjects (ug/mL)")          # Table 3, SD (healthy) final model
    propSdPatient <- 0.285;  label("Proportional residual SD, patients (fraction)")           # Table 3, CV% (patient) final model 28.5%
    addSdPatient  <- 0.236;  label("Additive residual SD, patients (ug/mL)")                  # Table 3, SD (patient) final model
  })

  model({
    # ---------------------------------------------------------------
    # Individual PK parameters from the final-model equations.
    # CL/F is linear in min(CRCL, th_bp) with an IBW power adjustment.
    # The linear-ratio form (CRCL / 86) is mathematically identical to the
    # paper's slope form q_CL/F * CRCL (with q_CL/F = exp(lcl) / 86 = 0.0462);
    # we use the ratio form so etalcl pairs with the canonical lcl.
    crcl_eff <- min(CRCL, th_bp)
    cl       <- exp(lcl + etalcl) * (crcl_eff / 86) * (IBW / 62)^e_ibw_cl

    # V/F: power-of-ratio scaling on IBW, BMI, AGE; categorical female
    # multiplier e_sexf_vc applied as e_sexf_vc^SEXF (1 for males, 0.906 for
    # females). Reference covariate values: IBW = 62 kg, BMI = 25 kg/m^2,
    # AGE = 59 years, SEXF = 0 (male).
    vc <- exp(lvc + etalvc) *
            (IBW / 62)^e_ibw_vc *
            (BMI / 25)^e_bmi_vc *
            (AGE / 59)^e_age_vc *
            e_sexf_vc^SEXF

    # Absorption rate and lag-time with food-effect multipliers.
    ka   <- exp(lka + etalka) * (1 + e_food_ka   * FED)
    tlag <- exp(ltlag)        * (1 + e_food_tlag * FED)

    # Elimination first-order rate.
    kel <- cl / vc

    # ODE system: depot -> central -> elimination. Standard 1-cmt FO oral.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag-time on the depot compartment.
    alag(depot) <- tlag

    # Observed plasma concentration: central (mg) / vc (L) = mg/L = ug/mL,
    # matching the bioanalytical assay units.
    Cc <- central / vc

    # ---------------------------------------------------------------
    # Residual error: combined proportional + additive, stratified by
    # healthy-vs-patient cohort via DIS_HEALTHY. Exactly one of the two
    # indicator products is 1 for any given subject; the corresponding
    # propSd / addSd pair is selected at runtime.
    propSd <- propSdHealthy * DIS_HEALTHY + propSdPatient * (1 - DIS_HEALTHY)
    addSd  <- addSdHealthy  * DIS_HEALTHY + addSdPatient  * (1 - DIS_HEALTHY)

    Cc ~ prop(propSd) + add(addSd)
  })
}
