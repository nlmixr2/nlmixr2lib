Dings_2026_neonatal_acidosis <- function() {
  description <- paste(
    "Neonatal acid-base outcome models from the HYPOTENS caesarean-section",
    "study (Dings 2026, Appendix A.2): three stepwise multiple linear",
    "regressions predicting umbilical-arterial pH, base excess and lactate",
    "in newborns whose mothers were treated for spinal-anaesthesia-induced",
    "hypotension with cafedrine/theodrenaline (C/T) or ephedrine. These are",
    "static algebraic predictors, not ODE models: there are no dosing",
    "events and no time course. Two of the base-excess predictors are",
    "per-subject parameter estimates from the companion maternal",
    "haemodynamic K/PD model `Dings_2026_cafedrine_theodrenaline_ephedrine.R`",
    "(the attainable heart-rate ceiling MAXHR and the K/PD elimination rate",
    "kel), and the strongest predictor of both pH and base excess is the",
    "cumulative minutes of maternal MAP below her pre-surgery MAP, which is",
    "derived from that model's simulated MAP-time profile. The three",
    "endpoints explain 6.06%, 22.1% and 36.7% of observed variance",
    "respectively, so they are hypothesis-generating rather than",
    "predictive at the individual level."
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

  covariateData <- list(
    DUR_MAP_BELOW_PRESURG = list(
      description        = "Cumulative minutes for which maternal MAP stayed below her own pre-surgery MAP, from antihypotensive treatment until umbilical-cord clamping",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper's TIMEMAP. Derived from the simulated (or observed) maternal MAP-time profile of the companion model, NOT measured directly. Dings 2026 Section 2.3 states the convention when MAP never recovers within the window: TIMEMAP is then set to the median interval between drug administration and cord clamping. Cohort means 10.1 min (C/T) and 10.2 min (ephedrine), Table 2. Retained in both the pH (Table A2) and base-excess (Table A3) regressions and is the strongest predictor of each; NOT retained for lactate (Table A4).",
      source_name        = "Duration MAP<pre-surgery MAP / TIMEMAP"
    ),
    EGA = list(
      description        = "Pregnancy duration at delivery",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dings 2026 reports this in DAYS (Table A2 coefficient is per day; cohort mean (SD) 262 (10.2) days). The canonical column is in weeks, so model() converts with EGA * 7 before applying the printed per-day coefficient - the coefficient value itself is unchanged from Table A2. Retained only in the pH regression.",
      source_name        = "Pregnancy duration"
    ),
    SPINAL_BLOCK = list(
      description        = "Highest anaesthetised dermatome of the maternal spinal block, as a thoracic segment number (5 = T5)",
      units              = "thoracic segment number",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear per-segment terms, POSITIVE in both Table A2 (+0.00534 pH per segment) and Table A3 (+0.3 mmol/L base excess per segment). Because a LOWER segment number is an anatomically HIGHER, more extensive block, a positive coefficient means a higher block predicts WORSE neonatal acid-base status - consistent with the negative power exponent on the anaesthesia blood-pressure effect in the companion maternal model. NOTE: the Section 3.3 prose describes both of these as decreases, which contradicts the signs printed in Tables A2 / A3; the tables are authoritative (see vignette Errata). Cohort means 5.79 (C/T) and 5.02 (ephedrine).",
      source_name        = "Spinal block height"
    ),
    TRT_CAFEDRINE_THEODRENALINE = list(
      description        = "Indicator that the mother was treated with cafedrine/theodrenaline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (untreated - cord clamping occurred before the mother received any antihypotensive)",
      notes              = "Reference here is UNTREATED, not the other drug: Dings 2026 Table A3 carries separate indicators for C/T (+3.26 mmol/L) and ephedrine (+2.64 mmol/L) against an untreated reference, and Figure A4 confirms the three-level coding (0 = none, 1 = ephedrine, 2 = C/T). Newborns counted as untreated: 10 of 228 pH (4.4%), 8 of 205 base excess (3.9%), 1 of 76 lactate (1.3%). Retained only in the base-excess regression.",
      source_name        = "Treatment with C/T"
    ),
    TRT_EPHEDRINE = list(
      description        = "Indicator that the mother was treated with ephedrine",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 - but the reference DIFFERS BY ENDPOINT: untreated for base excess (Table A3), versus C/T-or-untreated for lactate (Table A4)",
      notes              = "In the base-excess regression this is one of two treatment indicators against an untreated reference (+2.64 mmol/L). In the lactate regression it is the ONLY treatment indicator and Section 3.3 states the contrast explicitly as ephedrine 'compared with treatment with C/T or none' (+5.46 mg/dL), so its reference pools C/T with untreated. Encoding both endpoints from a single column therefore requires the two different contrasts to be respected, which model() does by using TRT_CAFEDRINE_THEODRENALINE only in the base-excess equation.",
      source_name        = "Treatment with E / Treatment with ephedrine"
    ),
    DOSE_BUPIVACAINE_MG = list(
      description        = "Intrathecal bupivacaine dose given to the mother for the spinal block",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the base-excess (-0.186 mmol/L per mg, Table A3) and lactate (+1.05 mg/dL per mg, Table A4) regressions; not retained for pH. Cohort means 10.5 (C/T) and 10.2 (ephedrine) mg.",
      source_name        = "Bupivacaine amount"
    ),
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MATERNAL weight, not the newborn's (Table A3 row label 'Weight (mother)'). +0.0176 mmol/L base excess per kg. Cohort mean (SD) 87.6 (19.7) kg. Retained only in the base-excess regression. The newborn's own weight was screened (Section 2.2) but not retained in any of the three models.",
      source_name        = "Weight (mother)"
    ),
    MAXHR_KPD = list(
      description        = "Mother's individual maximum attainable heart rate from the companion maternal haemodynamic K/PD fit",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A model-derived input, not a measurement: obtain it as the rmax_hr parameter of `Dings_2026_cafedrine_theodrenaline_ephedrine.R` for the same subject. Population typical values 77.8 (C/T) and 89.5 (ephedrine) beats/min. The -0.0256 mmol/L per beat/min coefficient (Table A3) is the quantitative basis for the paper's conclusion that a lower maternal peak heart rate benefits the newborn, and is the mechanism by which the choice of vasopressor reaches the neonatal endpoint. Retained only in the base-excess regression.",
      source_name        = "Model parameter MAXHR"
    ),
    KEL_KPD = list(
      description        = "Mother's individual K/PD elimination rate constant from the companion maternal haemodynamic K/PD fit",
      units              = "1/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A model-derived input: the kel = cl/vc of `Dings_2026_cafedrine_theodrenaline_ephedrine.R` for the same subject, population typical value 0.0971/min. Within a K/PD framework this describes how fast the drug EFFECT dissipates, not systemic clearance - no drug concentrations were measured - and the paper's Discussion cautions explicitly against a causal reading. A slower effect decay (smaller kel) predicts better base excess (-1.74 mmol/L per 1/min, Table A3). Retained only in the base-excess regression.",
      source_name        = "Model parameter kel"
    ),
    T_CORDCLAMP_SAMPLING = list(
      description        = "Delay between umbilical-cord clamping and drawing the umbilical-arterial blood-gas sample",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A pre-analytical covariate: continued anaerobic metabolism in the sampled blood raises measured lactate with elapsed time (+0.708 mg/dL per min, Table A4 - the largest single term in the lactate model). It is also a confounder for the between-treatment lactate comparison, because the delay itself differed significantly by arm (means 7.81 min for C/T vs 4.65 min for ephedrine, p < 0.001, Table 2). Retained only in the lactate regression.",
      source_name        = "Cord clamping to blood sampling"
    )
  )

  covariatesDataExcluded <- list(
    WT_BIRTH = list(
      description = "Newborn birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a candidate predictor (Dings 2026 Section 2.2, 'the neonate's weight and sex') but not retained in any of the three regressions. Cohort median (SD) 3370 (547) g; means 3440 (C/T) and 3290 (ephedrine) g, the only demographic difference between arms (p = 0.03, Table 1). Note that MATERNAL weight (WT) WAS retained, for base excess."
    ),
    SEXF = list(
      description = "Newborn female-sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Section 2.2) but not retained in any regression. Newborn sex distribution 42.2% female (C/T) and 47.2% female (ephedrine), p = 0.52 (Table 1)."
    ),
    BMI = list(
      description = "Maternal body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened via 'all variables considered in the maternal hemodynamic modeling' (Section 2.2) but not retained in any neonatal regression, although it IS retained on clearance in the companion maternal model. Cohort mean (SD) 31.6 (6.58) kg/m^2."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 243L,
    n_studies      = 1L,
    age_range      = "newborn (at delivery)",
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = "3370 g (newborn; SD 547)",
    sex_female_pct = 44.4,
    race_ethnicity = "Not reported in Dings 2026; multicentre German cohort.",
    disease_state  = "Newborns delivered by elective caesarean section under spinal anaesthesia, whose mothers developed and were treated for spinal-anaesthesia-induced hypotension. Acidosis was uncommon: 12 newborns (5.3%) had pH < 7.2 (4 of 135 C/T, 8 of 108 ephedrine; p = 0.133) and none exceeded the critical thresholds of -12 mmol/L base excess or 90 mg/dL lactate.",
    dose_range     = "Not applicable - the mother is dosed, not the newborn. Maternal doses are described in the companion model.",
    regions        = "Germany (multicentre)",
    notes          = "Endpoint availability differs by biomarker and sets each regression's n: umbilical-arterial pH 228 newborns (93.8%), base excess 205 (84.4%), lactate only 76 (31.3%), so lactate conclusions are the least generalisable. Model fits: pH F(3;224) = 5.88, p = 0.000703, adjusted R^2 = 0.0606; base excess F(8;196) = 8.23, p = 1.36e-9, adjusted R^2 = 0.221; lactate F(3;72) = 15.5, p = 6.99e-8, adjusted R^2 = 0.367. Newborn sex distribution 57.8% male (C/T) / 52.8% male (ephedrine); sex_female_pct above is the pooled female share. Estimates from Dings 2026 Tables A2, A3 and A4."
  )

  ini({
    # ---- Umbilical-arterial pH (Dings 2026 Table A2) --------------------
    # Stepwise forward selection (p = 0.1) / backward elimination
    # (p = 0.05) via the mixlm R package. Coefficients are on the natural
    # scale of each predictor and are NOT log-transformed: several are
    # negative, and these are linear-regression slopes rather than
    # strictly-positive structural parameters.
    ph_ua_int <- 7.52; label("Umbilical-arterial pH intercept")                                          # Table A2 Intercept, RSE 1.31%, t = 76.6
    e_dur_map_below_presurg_ph <- -0.00177; label("pH change per minute of maternal MAP below pre-surgery MAP")  # Table A2, RSE 36.2%, p = 0.00613
    e_ega_ph <- -0.000822; label("pH change per additional day of pregnancy duration")                   # Table A2 "Pregnancy duration", RSE 45.7%, p = 0.0298
    e_spinal_block_ph <- 0.00534; label("pH change per additional spinal segment number")                # Table A2, RSE 46.1%, p = 0.0309

    # ---- Umbilical-arterial base excess (Table A3) ----------------------
    # Reference for both treatment indicators is UNTREATED.
    be_ua_int <- -2.09; label("Umbilical-arterial base excess intercept (mmol/L)")                       # Table A3 Intercept, RSE 75.6%, p = 0.188 (not significant)
    e_trt_cafedrine_theodrenaline_be <- 3.26; label("Base-excess change for maternal C/T vs untreated (mmol/L)")   # Table A3 "Treatment with C/T", RSE 23.8%, p < 0.001
    e_trt_ephedrine_be <- 2.64; label("Base-excess change for maternal ephedrine vs untreated (mmol/L)")           # Table A3 "Treatment with E", RSE 29.7%, p < 0.001
    e_dur_map_below_presurg_be <- -0.105; label("Base-excess change per minute of maternal MAP below pre-surgery MAP (mmol/L/min)")  # Table A3, RSE 24.6%, p < 0.001
    e_spinal_block_be <- 0.3; label("Base-excess change per additional spinal segment number (mmol/L)")            # Table A3, RSE 30.7%, p = 0.00134
    e_dose_bupivacaine_mg_be <- -0.186; label("Base-excess change per mg maternal bupivacaine (mmol/L/mg)")        # Table A3 "Bupivacaine amount", RSE 35.5%, p = 0.00536
    e_wt_be <- 0.0176; label("Base-excess change per kg maternal body weight (mmol/L/kg)")                         # Table A3 "Weight (mother)", RSE 41.6%, p = 0.0171
    e_maxhr_kpd_be <- -0.0256; label("Base-excess change per beat/min of the maternal K/PD MAXHR (mmol/L/bpm)")    # Table A3 "Model parameter MAXHR", RSE 42.6%, p = 0.0202
    e_kel_kpd_be <- -1.74; label("Base-excess change per 1/min of the maternal K/PD kel (mmol/L per 1/min)")       # Table A3 "Model parameter kel", RSE 47.6%, p = 0.0368

    # ---- Umbilical-arterial lactate (Table A4) --------------------------
    # Here the ephedrine reference pools C/T with untreated (Section 3.3),
    # so no C/T indicator appears.
    lactate_ua_int <- 3.43; label("Umbilical-arterial lactate intercept (mg/dL)")                        # Table A4 Intercept, RSE 115%, p = 0.389 (not significant)
    e_t_cordclamp_sampling_lactate <- 0.708; label("Lactate change per minute from cord clamping to sampling (mg/dL/min)")  # Table A4, RSE 18.9%, p < 0.001
    e_trt_ephedrine_lactate <- 5.46; label("Lactate change for maternal ephedrine vs C/T or untreated (mg/dL)")             # Table A4 "Treatment with ephedrine", RSE 37.2%, p = 0.00883
    e_dose_bupivacaine_mg_lactate <- 1.05; label("Lactate change per mg maternal bupivacaine (mg/dL/mg)")                   # Table A4 "Bupivacaine amount", RSE 37.5%, p = 0.00967

    # ---- Residual variability -------------------------------------------
    # Dings 2026 reports adjusted R^2 and per-coefficient RSEs but no
    # residual standard error for any of the three regressions, so the
    # residual SDs cannot be sourced and are fixed to 0: these models
    # reproduce the published point predictions only. See the vignette
    # Errata. Do NOT substitute an invented value; a user needing
    # predictive intervals must refit or obtain the residual SE from the
    # authors.
    addSd_ph_ua <- fixed(0); label("Additive residual SD for umbilical-arterial pH (not reported in the source)")
    addSd_be_ua <- fixed(0); label("Additive residual SD for umbilical-arterial base excess (not reported in the source)")
    addSd_lactate_ua <- fixed(0); label("Additive residual SD for umbilical-arterial lactate (not reported in the source)")
  })

  model({
    # Pregnancy duration: the canonical EGA column is in weeks while the
    # Table A2 coefficient is per DAY, so convert the covariate rather
    # than rescaling the published coefficient.
    preg_days <- EGA * 7

    # Umbilical-arterial pH (Table A2).
    ph_ua <- ph_ua_int +
      e_dur_map_below_presurg_ph * DUR_MAP_BELOW_PRESURG +
      e_ega_ph * preg_days +
      e_spinal_block_ph * SPINAL_BLOCK

    # Umbilical-arterial base excess (Table A3). Both treatment
    # indicators are present because the reference is untreated.
    be_ua <- be_ua_int +
      e_trt_cafedrine_theodrenaline_be * TRT_CAFEDRINE_THEODRENALINE +
      e_trt_ephedrine_be * TRT_EPHEDRINE +
      e_dur_map_below_presurg_be * DUR_MAP_BELOW_PRESURG +
      e_spinal_block_be * SPINAL_BLOCK +
      e_dose_bupivacaine_mg_be * DOSE_BUPIVACAINE_MG +
      e_wt_be * WT +
      e_maxhr_kpd_be * MAXHR_KPD +
      e_kel_kpd_be * KEL_KPD

    # Umbilical-arterial lactate (Table A4). Only the ephedrine
    # indicator appears; its reference pools C/T with untreated.
    lactate_ua <- lactate_ua_int +
      e_t_cordclamp_sampling_lactate * T_CORDCLAMP_SAMPLING +
      e_trt_ephedrine_lactate * TRT_EPHEDRINE +
      e_dose_bupivacaine_mg_lactate * DOSE_BUPIVACAINE_MG

    ph_ua ~ add(addSd_ph_ua)
    be_ua ~ add(addSd_be_ua)
    lactate_ua ~ add(addSd_lactate_ua)
  })
}
