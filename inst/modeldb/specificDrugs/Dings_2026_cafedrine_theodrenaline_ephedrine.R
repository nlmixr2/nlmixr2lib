Dings_2026_cafedrine_theodrenaline_ephedrine <- function() {
  description <- paste(
    "Population kinetic/pharmacodynamic (K/PD) model of maternal",
    "haemodynamics after intravenous cafedrine/theodrenaline (C/T, 20:1",
    "fixed combination) or ephedrine for spinal-anaesthesia-induced",
    "hypotension during elective caesarean section (Dings 2026, HYPOTENS",
    "study). No pharmacokinetic samples were taken, so exposure is",
    "represented by a virtual K/PD compartment with first-order",
    "elimination and the apparent volume fixed to 1 L for identifiability.",
    "The K/PD concentration feeds TWO PARALLEL first-order lag cascades of",
    "different speed - three slow stages (mean transit time 28.0 min,",
    "carrying persistence) and four fast stages (2.27 min, carrying onset)",
    "- whose terminal members sum to one driver concentration. That driver",
    "acts through three parallel Emax relationships on heart rate, mean",
    "arterial pressure and systolic blood pressure, each anchored on the",
    "subject's own at-diagnosis baseline. A separate decaying state carries",
    "the progressive blood-pressure fall from the spinal anaesthetic",
    "itself, and step effects capture surgical incision and uterotomy.",
    "Ephedrine raises the attainable heart-rate ceiling by 15% and lowers",
    "all three EC50 values by 71.1%, which is what makes 50 mg C/T and",
    "15 mg ephedrine comparable for blood pressure while ephedrine causes",
    "substantially more tachycardia. Companion file",
    "`Dings_2026_neonatal_acidosis.R` consumes this model's per-subject",
    "MAXHR and kel as covariates on the newborn acid-base endpoints."
  )
  reference <- paste(
    "Dings C, Lehr T, Kranke P, Vojnar B, Gaik C, Koch T, Eberhart L,",
    "Huljic-Lankinen S, Murst M, Kreuer S. Pharmacometric Analysis of",
    "Cafedrine/Theodrenaline Versus Ephedrine on Maternal Hemodynamics and",
    "Neonatal Acidosis During Cesarean Section. Pharmaceutics.",
    "2026;18(3):296. doi:10.3390/pharmaceutics18030296"
  )
  vignette <- "Dings_2026_cafedrine_theodrenaline_ephedrine"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # `effect_anaesthesia` is the paper's A9 (Dings 2026 Eqs A9 / A11): a
  # virtual state holding a hypothetical unit amount dosed at the time of
  # spinal anaesthesia and decaying at kdel. It is not a drug compartment
  # and is scoped to this paper's surgical-anaesthesia mechanism.
  paper_specific_compartments <- c("effect_anaesthesia")

  compartmentData <- list(
    depot_kpd = list(
      analyte = "cafedrine/theodrenaline or ephedrine", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    # The lag-cascade stages are driven by depot_kpd/vc (a concentration),
    # not by a mass transfer, so they carry concentration units (mg/L)
    # rather than amounts. Dings 2026 Eqs A2-A8.
    effect_slow1 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_slow2 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_slow3 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_fast1 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_fast2 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_fast3 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_fast4 = list(analyte = "cafedrine/theodrenaline or ephedrine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    effect_anaesthesia = list(analyte = "spinal anaesthesia effect (virtual unit amount)", units = "(unitless)", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Maternal body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on K/PD clearance, centred at the population median of 30 kg/m^2 (Dings 2026 Section 3.4 typical parturient; Section 2.3 defines the typical parturient as the median covariate values). Cohort mean (SD) 31.6 (6.58) kg/m^2. A higher BMI lowers clearance (exponent -0.948) and therefore lengthens the effect half-life.",
      source_name        = "BMI"
    ),
    DBP_BL = list(
      description        = "Maternal diastolic blood pressure at hypotension diagnosis, immediately before the first antihypotensive bolus",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the K/PD apparent volume, centred at 50 mmHg. Dings 2026 does not tabulate baseline DBP in Table 2; the 50 mmHg centring value is the typical-parturient value from Section 3.4, and is internally consistent with the tabulated baseline MAP and SBP via MAP = DBP + (SBP - DBP)/3 (50 + (92 - 50)/3 = 64 mmHg, matching the reported baseline MAP).",
      source_name        = "DBP (baseline row, Table 2 / Section 3.4)"
    ),
    HR_PRESURG = list(
      description        = "Maternal heart rate before spinal anaesthesia (undisturbed pre-surgical set-point)",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the attainable heart-rate ceiling, centred at the population median of 92 beats/min. Cohort means 92.6 (C/T) and 91.7 (ephedrine) beats/min, Dings 2026 Table 2 pre-surgery block. Distinct from HR_BL, which is measured after anaesthesia at hypotension diagnosis.",
      source_name        = "HR (pre-surgery row, Table 2)"
    ),
    HR_BL = list(
      description        = "Maternal heart rate at hypotension diagnosis, immediately before the first antihypotensive bolus",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Additive baseline anchor of the heart-rate Emax term (Dings 2026 Eq. A12: BL_HR). No covariate effect is estimated on it. Typical value 84 beats/min (Section 3.4); cohort means 85.4 (C/T) and 87.5 (ephedrine) beats/min. Because the estimated ceiling for C/T (77.8 beats/min) lies BELOW this baseline, the recovered Emax is negative for C/T and positive for ephedrine - this is the paper's central finding that C/T is essentially heart-rate-neutral.",
      source_name        = "HR (baseline row, Table 2)"
    ),
    MAP_PRESURG = list(
      description        = "Maternal mean arterial pressure before spinal anaesthesia",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the attainable MAP ceiling, centred at the population median of 100 mmHg. Cohort mean 101 mmHg in both arms (Dings 2026 Table 2 pre-surgery block). Enters with the OPPOSITE sign to MAP_BL on the same parameter, which is why the two anchors must be separate columns.",
      source_name        = "MAP (pre-surgery row, Table 2)"
    ),
    MAP_BL = list(
      description        = "Maternal mean arterial pressure at hypotension diagnosis, immediately before the first antihypotensive bolus",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Serves two roles: the additive baseline anchor of the MAP Emax term (Dings 2026 Eq. A13: BL_MAP) and a power covariate on the MAP ceiling with exponent -0.53, centred at 64 mmHg. Cohort means 64.2 (C/T) and 64.5 (ephedrine) mmHg.",
      source_name        = "MAP (baseline row, Table 2)"
    ),
    SBP_PRESURG = list(
      description        = "Maternal systolic blood pressure before spinal anaesthesia",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the attainable SBP ceiling, centred at the population median of 139 mmHg. Cohort means 139 (C/T) and 138 (ephedrine) mmHg.",
      source_name        = "SBP (pre-surgery row, Table 2)"
    ),
    SBP_BL = list(
      description        = "Maternal systolic blood pressure at hypotension diagnosis, immediately before the first antihypotensive bolus",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Additive baseline anchor of the SBP Emax term (Dings 2026 Eq. A14: BL_SBP) and a power covariate on the SBP ceiling with exponent -0.494, centred at 92 mmHg. Cohort means 92.3 (C/T) and 93 (ephedrine) mmHg.",
      source_name        = "SBP (baseline row, Table 2)"
    ),
    SPINAL_BLOCK = list(
      description        = "Highest anaesthetised dermatome of the spinal block, as a thoracic segment number (5 = T5)",
      units              = "thoracic segment number",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the additive anaesthesia blood-pressure effect, centred at the population median of 5 (T5). THE ORIENTATION IS INVERTED relative to the phrase 'block height': a LOWER segment number is an anatomically HIGHER, more extensive block. With exponent -2.05, a lower number therefore predicts a LARGER fall in blood pressure, matching Dings 2026 Section 3.2 ('larger decreases observed at lower thoracic block levels') and the positive per-segment coefficients in the companion neonatal model. Cohort means 5.79 (C/T) and 5.02 (ephedrine).",
      source_name        = "Spinal block height / SPINAL"
    ),
    TRT_EPHEDRINE = list(
      description        = "Indicator that the parturient was treated with ephedrine rather than cafedrine/theodrenaline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cafedrine/theodrenaline arm; every subject in this analysis received one drug or the other, so C/T is the reference)",
      notes              = "Enters via the paper's categorical form `1 + COV * theta` (Dings 2026 Eq. 2) on two parameters: the attainable heart-rate ceiling (+15%, Table A1 F_MAXHR = 0.15) and ALL THREE EC50 values (-71.1%, Table A1 F_EC50 = -0.711 - a single estimated parameter shared across the HR, MAP and SBP EC50s, confirmed by Section 3.2 quoting all three ephedrine EC50 values as 71.1% lower). n = 108 of 243 (44.4%) received ephedrine.",
      source_name        = "Treatment with E"
    ),
    T_ANAESTHESIA = list(
      description        = "Time from induction of spinal anaesthesia to the first antihypotensive bolus (model time zero)",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sets the initial condition of effect_anaesthesia (Dings 2026 Eqs A9 / A11). Typical value 7 min (Section 3.4). The construction makes the anaesthesia-driven blood-pressure fall already realised before treatment vanish at t = 0 - it is absorbed into the observed MAP_BL / SBP_BL anchors - so only the REMAINING fall is predicted forward.",
      source_name        = "TIMEAN"
    ),
    T_INCISION = list(
      description        = "Time from the first antihypotensive bolus to surgical incision",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; the 0/1 indicator I of Dings 2026 Eqs A12-A14 is derived inside model() as `t >= T_INCISION`, so no time-varying covariate column is needed. Typical value 5 min after treatment (Section 3.4).",
      source_name        = "I (derived indicator)"
    ),
    T_UTEROTOMY = list(
      description        = "Time from the first antihypotensive bolus to uterotomy",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; the 0/1 indicator U of Dings 2026 Eqs A12-A14 is derived inside model() as `t >= T_UTEROTOMY`. Typical value 9 min after treatment (Section 3.4). Each of the two events adds one increment of the stimulation effect, so `iu` runs 0 -> 1 -> 2.",
      source_name        = "U (derived indicator)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Maternal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the Dings 2026 covariate analysis (Section 2.2, maternal demographics) but NOT retained in the final haemodynamic model; BMI was retained instead. Cohort mean (SD) 87.6 (19.7) kg. Maternal weight IS retained in the companion `Dings_2026_neonatal_acidosis.R` base-excess regression."
    ),
    HT = list(
      description = "Maternal body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened (Section 2.2) but not retained; BMI carries the size effect. Cohort mean (SD) 166 (6.91) cm."
    ),
    DOSE_BUPIVACAINE_MG = list(
      description = "Intrathecal bupivacaine dose given for the spinal block",
      units       = "mg",
      type        = "continuous",
      notes       = "Screened as an anaesthesia-related covariate (Section 2.2) but not retained in the maternal haemodynamic model. Cohort means 10.5 (C/T) and 10.2 (ephedrine) mg. IS retained in the companion neonatal model."
    ),
    EGA = list(
      description = "Pregnancy duration at delivery",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened as 'pregnancy duration' (Section 2.2) but not retained in the maternal haemodynamic model. Cohort mean (SD) 262 (10.2) days = 37.4 weeks. IS retained in the companion neonatal pH regression."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 243L,
    n_studies      = 1L,
    age_range      = NA_character_,
    age_median     = "33.4 years (mean; SD 5.22)",
    weight_range   = NA_character_,
    weight_median  = "87.6 kg (mean; SD 19.7)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported in Dings 2026; multicentre German cohort.",
    disease_state  = "Spinal-anaesthesia-induced hypotension (SBP < 100 mmHg or < 90% of the pre-operative baseline SBP) during elective caesarean section in parturients at term. Prophylactically treated patients were excluded.",
    dose_range     = "Intravenous boluses at the attending anaesthesiologist's discretion. Cafedrine/theodrenaline 10-200 mg per dose (most common single dose 50 mg; median cumulative dose over 30 min 100 mg, range 20-400 mg), expressed as cafedrine equivalents. Ephedrine 5-40 mg per dose (most common 15 mg; median cumulative 30 mg, range 10-120 mg).",
    regions        = "Germany (multicentre)",
    notes          = "HYPOTENS study (Eberhart 2018, Curr Med Res Opin 34:953-961; primary results Kranke 2021, Eur J Anaesthesiol 38:1067-1076), a prospective, open-label, two-armed, non-interventional study conducted July 2016 - February 2018. Of 283 per-protocol patients, 40 (14.1%) were excluded for blood loss > 1000 mL, colloid / blood-product administration, conversion to general anaesthesia, or a local anaesthetic other than bupivacaine, leaving 243 analysed: 135 (55.6%) cafedrine/theodrenaline and 108 (44.4%) ephedrine. Haemodynamics were recorded immediately before the first bolus and at 1-10, 12, 15, 20, 25 and 30 min after it. Demographics from Dings 2026 Tables 1 and 2. NOTE: no pharmacokinetic samples were obtained, so every kinetic parameter here is a functional descriptor of the response time course, NOT a systemic pharmacokinetic quantity."
  )

  ini({
    # ---- K/PD kinetics (Dings 2026 Table A1, Eq. A1) --------------------
    # No PK data exist for cafedrine, theodrenaline or ephedrine, so
    # exposure is a virtual K/PD compartment. V is FIXED to 1 L purely for
    # identifiability and is a scaling constant, not a physiological
    # volume (Section 2.2). kel = cl/vc = 0.0971/min gives the reported
    # 7.14 min effect half-life (Section 3.2).
    lcl <- log(0.0971); label("K/PD clearance (L/min)")                        # Table A1 "Clearance", RSE 11.2%
    lvc <- fixed(log(1)); label("K/PD apparent volume of distribution (L)")    # Table A1 "Volume of distribution", 1 (fixed)

    # ---- Dual-rate effect-delay cascades (Eqs A2-A8) --------------------
    # Slow chain = 3 stages, mean transit time 3/0.107 = 28.0 min.
    # Fast chain = 4 stages, mean transit time 4/1.76  =  2.27 min.
    # Both transit times are reported in Section 3.2, which is what pins
    # the stage counts to the chains (the paper's Ktr1 is the SLOWER
    # constant, so do not map Ktr1/Ktr2 onto slow/fast by their digit).
    lktr_slow <- log(0.107); label("Transit rate of the slow (long-term) effect cascade (1/min)")  # Table A1 "Transit rate long-term effect" Ktr1, RSE 7.3%
    lktr_fast <- log(1.76); label("Transit rate of the fast (short-term) effect cascade (1/min)")  # Table A1 "Transit rate short-term effect" Ktr2, RSE 5.9%

    # ---- Emax ceilings and potencies (Eqs A12-A14) ---------------------
    # The paper estimates the maximum attainable LEVEL of each endpoint
    # (MAXi), and Emax is recovered as the difference MAXi - BLi. Hence
    # lrmax_* (plateau level) rather than lemax (increment): for C/T the
    # heart-rate ceiling 77.8 lies below the 84 beats/min baseline, so the
    # recovered Emax is NEGATIVE and heart rate falls.
    lrmax_hr <- log(77.8); label("Maximum attainable heart rate (beats/min)")              # Table A1 "Maximum HR" MAXHR, RSE 2.1%
    lrmax_map <- log(119); label("Maximum attainable mean arterial pressure (mmHg)")       # Table A1 "Maximum MAP" MAXMAP, RSE 4.8%
    lrmax_sbp <- log(169); label("Maximum attainable systolic blood pressure (mmHg)")      # Table A1 "Maximum SBP" MAXSBP, RSE 4.7%
    lec50_hr <- log(5.63); label("Driver concentration for half-maximal heart-rate effect (mg/L)")        # Table A1 "Half-effective concentration HR" EC50HR, RSE 46.9%
    lec50_map <- log(61.7); label("Driver concentration for half-maximal MAP effect (mg/L)")              # Table A1 "Half-effective concentration MAP" EC50MAP, RSE 23.3%
    lec50_sbp <- log(64.4); label("Driver concentration for half-maximal SBP effect (mg/L)")              # Table A1 "Half-effective concentration SBP" EC50SBP, RSE 23.4%

    # ---- Spinal-anaesthesia and intraoperative-event effects -----------
    # eff_anaes_bp is negative (a fall) and is applied to MAP and SBP
    # alike; eff_event_pp is a PULSE-PRESSURE effect, so SBP receives it
    # in full and MAP receives one third of it (Eqs A13-A14). Section 3.2
    # confirms the arithmetic: 2.39/3 = 0.797 ~ the reported 0.8 mmHg MAP
    # rise, against 2.39 mmHg for SBP.
    lkdel <- log(0.171); label("Rate constant of the anaesthesia blood-pressure effect (1/min)")  # Table A1 "Rate of BP effect of anesthesia" Kdel, RSE 47.3%
    eff_anaes_bp <- -9.78; label("Anaesthesia effect on MAP and SBP at a T5 block (mmHg)")        # Table A1 "Anesthesia effect on BP" EffBP, RSE 29.7%
    eff_event_hr <- 2.75; label("Incision / uterotomy effect on heart rate, per event (beats/min)")  # Table A1 "Incision/uterotomy effect on HR" EffHR, RSE 18.9%
    eff_event_pp <- 2.39; label("Incision / uterotomy effect on pulse pressure, per event (mmHg)")   # Table A1 "Incision/uterotomy effect on PP" EffPP, RSE 21.7%

    # ---- Covariate effects (Table A1 "Covariate Effects" block) --------
    # Continuous covariates use the paper's power form centred at the
    # population median (Eq. 1); categorical covariates use the fractional
    # form 1 + COV*theta (Eq. 2). Centring values are the typical-
    # parturient values of Section 3.4, which Section 2.3 defines as the
    # population medians.
    e_bmi_cl <- -0.948; label("Power exponent of BMI on K/PD clearance, centred at 30 kg/m^2")                          # Table A1 "BMI on clearance" Cl(BMI), p = 9.4e-5
    e_dbp_bl_vc <- -1.49; label("Power exponent of baseline DBP on K/PD volume, centred at 50 mmHg")                    # Table A1 "Baseline DBP on V" V(DBPBL), p = 1.6e-6
    e_hr_presurg_rmax_hr <- 0.399; label("Power exponent of pre-surgery HR on the HR ceiling, centred at 92 beats/min") # Table A1 "Pre-surgery HR on MAXHR", p = 5.9e-9
    e_map_presurg_rmax_map <- 0.499; label("Power exponent of pre-surgery MAP on the MAP ceiling, centred at 100 mmHg") # Table A1 "Pre-surgery MAP on MAXMAP", p = 1.9e-11
    e_map_bl_rmax_map <- -0.53; label("Power exponent of baseline MAP on the MAP ceiling, centred at 64 mmHg")          # Table A1 "Baseline MAP on MAXMAP", p = 5.7e-17
    e_sbp_presurg_rmax_sbp <- 0.577; label("Power exponent of pre-surgery SBP on the SBP ceiling, centred at 139 mmHg") # Table A1 "Pre-surgery SBP on MAXSBP", p = 9.6e-16
    e_sbp_bl_rmax_sbp <- -0.494; label("Power exponent of baseline SBP on the SBP ceiling, centred at 92 mmHg")         # Table A1 "Baseline SBP on MAXSBP", p = 1.8e-14
    e_spinal_block_effbp <- -2.05; label("Power exponent of spinal block segment number on the anaesthesia BP effect, centred at T5")  # Table A1 "Spinal block height effect on EffBP" EffBP(SPINAL), p = 3.8e-6
    e_trt_ephedrine_rmax_hr <- 0.15; label("Fractional change in the HR ceiling for ephedrine vs cafedrine/theodrenaline")            # Table A1 "Fractional MAXHR change for E" F_MAXHR, RSE 16.6%
    e_trt_ephedrine_ec50 <- -0.711; label("Fractional change in ALL THREE EC50 values for ephedrine vs cafedrine/theodrenaline")      # Table A1 "Fractional EC50 change for Ephed." F_EC50, RSE 4%

    # ---- IIV (exponential; Table A1 "IIV [%CV]" column) ---------------
    # Reported as %CV for an exponential (log-normal) random effect, so
    # the nlmixr2 variance is log(1 + CV^2). The conversion is written out
    # rather than pre-evaluated so it stays auditable against Table A1.
    etalcl ~ log(1 + 0.669^2)             # Table A1 Clearance IIV 66.9 %CV, RSE 9.15%
    etalvc ~ log(1 + 0.775^2)             # Table A1 Volume IIV 77.5 %CV, RSE 13%
    etalktr_fast ~ log(1 + 0.691^2)       # Table A1 Ktr2 IIV 69.1 %CV, RSE 5.77%
    etalrmax_hr ~ log(1 + 0.158^2)        # Table A1 MAXHR IIV 15.8 %CV, RSE 7.35%
    etalrmax_map ~ log(1 + 0.165^2)       # Table A1 MAXMAP IIV 16.5 %CV, RSE 15.4%
    etalrmax_sbp ~ log(1 + 0.132^2)       # Table A1 MAXSBP IIV 13.2 %CV, RSE 18.2%

    # ---- Residual variability (Table A1 "Residual variability") -------
    # Proportional error reported as %CV; entered as a linear-scale SD.
    propSd_hr <- 0.13; label("Proportional residual SD for heart rate")               # Table A1 PEHR 13 %CV, RSE 5.4%
    propSd_map <- 0.0937; label("Proportional residual SD for MAP")                   # Table A1 PEMAP 9.37 %CV, RSE 6.1%
    propSd_sbp <- 0.0855; label("Proportional residual SD for SBP")                   # Table A1 PESBP 8.55 %CV, RSE 6.8%
  })

  model({
    # ---- Individual K/PD parameters ------------------------------------
    # Power covariate model, Dings 2026 Eq. (1): Effect = (COV/median)^theta.
    cl <- exp(lcl + e_bmi_cl * log(BMI / 30) + etalcl)
    vc <- exp(lvc + e_dbp_bl_vc * log(DBP_BL / 50) + etalvc)
    kel <- cl / vc
    ktr_slow <- exp(lktr_slow)
    ktr_fast <- exp(lktr_fast + etalktr_fast)
    kdel <- exp(lkdel)

    # ---- Ceilings and potencies ----------------------------------------
    # Categorical covariate model, Dings 2026 Eq. (2): Effect = 1 + COV*theta.
    # One shared F_EC50 scales all three EC50 values.
    ec50_fac <- 1 + e_trt_ephedrine_ec50 * TRT_EPHEDRINE
    ec50_hr <- exp(lec50_hr) * ec50_fac
    ec50_map <- exp(lec50_map) * ec50_fac
    ec50_sbp <- exp(lec50_sbp) * ec50_fac

    rmax_hr <- exp(lrmax_hr + e_hr_presurg_rmax_hr * log(HR_PRESURG / 92) + etalrmax_hr) *
      (1 + e_trt_ephedrine_rmax_hr * TRT_EPHEDRINE)
    rmax_map <- exp(lrmax_map + e_map_presurg_rmax_map * log(MAP_PRESURG / 100) +
                      e_map_bl_rmax_map * log(MAP_BL / 64) + etalrmax_map)
    rmax_sbp <- exp(lrmax_sbp + e_sbp_presurg_rmax_sbp * log(SBP_PRESURG / 139) +
                      e_sbp_bl_rmax_sbp * log(SBP_BL / 92) + etalrmax_sbp)

    # Power effect on a NEGATIVE additive parameter: a lower segment
    # number (an anatomically higher block) magnifies the fall.
    eff_bp <- eff_anaes_bp * (SPINAL_BLOCK / 5)^e_spinal_block_effbp

    # ---- K/PD compartment (Eq. A1) -------------------------------------
    d/dt(depot_kpd) <- -kel * depot_kpd
    kpd_conc <- depot_kpd / vc

    # ---- Slow cascade, 3 stages (Eqs A2-A4) ----------------------------
    d/dt(effect_slow1) <- ktr_slow * (kpd_conc - effect_slow1)
    d/dt(effect_slow2) <- ktr_slow * (effect_slow1 - effect_slow2)
    d/dt(effect_slow3) <- ktr_slow * (effect_slow2 - effect_slow3)

    # ---- Fast cascade, 4 stages (Eqs A5-A8) ----------------------------
    d/dt(effect_fast1) <- ktr_fast * (kpd_conc - effect_fast1)
    d/dt(effect_fast2) <- ktr_fast * (effect_fast1 - effect_fast2)
    d/dt(effect_fast3) <- ktr_fast * (effect_fast2 - effect_fast3)
    d/dt(effect_fast4) <- ktr_fast * (effect_fast3 - effect_fast4)

    # ---- Spinal-anaesthesia effect (Eqs A9, A11) -----------------------
    # A unit amount is dosed at the time of anaesthesia, T_ANAESTHESIA
    # minutes before model time zero, so at t = 0 the state already sits
    # at exp(-kdel*T_ANAESTHESIA). Setting that as the initial condition
    # is equivalent to the paper's pre-time-zero dose and avoids negative
    # event times. `anaes_eff` then rises from 0 at t = 0 toward that
    # value, so only the anaesthesia-driven fall STILL TO COME is added -
    # the part already realised is absorbed in the observed BL anchors.
    d/dt(effect_anaesthesia) <- -kdel * effect_anaesthesia
    effect_anaesthesia(0) <- exp(-kdel * T_ANAESTHESIA)
    anaes_eff <- exp(-kdel * T_ANAESTHESIA) - effect_anaesthesia

    # ---- Effect driver (Eq. A10) ---------------------------------------
    # Printed as "CONC = A4 + A9", which is a typo for A4 + A8: A9 is the
    # dimensionless anaesthesia state and already enters via anaes_eff,
    # whereas A8 is the terminal fast-cascade stage and would otherwise be
    # dangling along with A5-A7. See the vignette Errata for the numeric
    # evidence (the A4+A8 reading reproduces the published 30.4%
    # ephedrine tachycardia incidence; A4+A9 gives 23.4%).
    conc <- effect_slow3 + effect_fast4

    # Intraoperative stimulation: 0 before incision, 1 after incision,
    # 2 after uterotomy (Eqs A12-A14 covariates I and U).
    iu <- (t >= T_INCISION) + (t >= T_UTEROTOMY)

    # ---- Haemodynamic endpoints (Eqs A12-A14) --------------------------
    hr <- HR_BL + (rmax_hr - HR_BL) * conc / (conc + ec50_hr) + eff_event_hr * iu
    map <- MAP_BL + (rmax_map - MAP_BL) * conc / (conc + ec50_map) +
      eff_bp * anaes_eff + eff_event_pp / 3 * iu
    sbp <- SBP_BL + (rmax_sbp - SBP_BL) * conc / (conc + ec50_sbp) +
      eff_bp * anaes_eff + eff_event_pp * iu

    hr ~ prop(propSd_hr)
    map ~ prop(propSd_map)
    sbp ~ prop(propSd_sbp)
  })
}
