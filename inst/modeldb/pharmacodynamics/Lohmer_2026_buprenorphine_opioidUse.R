Lohmer_2026_buprenorphine_opioidUse <- function() {
  description <- paste(
    "Mixed-effects logistic exposure-response model of weekly opioid",
    "abstinence (negative urine drug screen AND negative self-report for",
    "past-week opioid use) under monthly subcutaneous extended-release",
    "buprenorphine (BUP-XR, SUBLOCADE) in adults with moderate-to-severe",
    "opioid use disorder, integrated across a pre-fentanyl phase 3 study",
    "(Study 1, NCT02357901, 2015-2016) and a high-risk fentanyl-era study",
    "(Study 2, NCT04995029, 2021-2024).",
    "logit P(abstinent) = alpha + Emax * C / (EC50 + C) + eta, with",
    "study-specific intercept alpha and Emax, a Study 1 300/100-mg arm",
    "factor and a baseline-employment factor on alpha, a Black/African",
    "American factor on Emax, and injection-route and daily-fentanyl-plus-",
    "amphetamine/methamphetamine factors on EC50. Study 1 uses the plasma",
    "buprenorphine concentration directly; Study 2 drives the effect",
    "through an effect compartment (ke0 half-life about 34 days) that",
    "captures its slower response. The model has no PK ODE: the",
    "buprenorphine plasma concentration CP_BPN_NGML is a required",
    "time-varying input (the authors used individual predictions from the",
    "Jones 2021 BUP-XR popPK model, packaged as",
    "modellib('Jones_2021_buprenorphine'))."
  )
  reference <- paste(
    "Lohmer LRL, Greenwald MK, Kakar R, Laffont CM. Optimal",
    "buprenorphine exposure to treat patients with opioid use disorder in",
    "the era of fentanyl and polysubstance use. Front Pharmacol.",
    "2026;17:1839067. doi:10.3389/fphar.2026.1839067.",
    "Parameter estimates from Table 2; model equations from the",
    "Supplementary Material ('Opioid Use Model Equations').",
    "Sister models from the same paper:",
    "modellib('Lohmer_2026_buprenorphine_craving'),",
    "modellib('Lohmer_2026_buprenorphine_dropout')."
  )
  vignette <- "Lohmer_2026_buprenorphine"

  depends <- c("CP_BPN_NGML")

  units <- list(
    time = "h",
    dosing = "N/A (PD-only; buprenorphine plasma concentration is a required input covariate)",
    concentration = "ng/mL (CP_BPN_NGML input covariate and the Study 2 effect compartment); the observation is a binary weekly opioid-use status"
  )

  compartmentData <- list(
    effect = list(
      analyte = "buprenorphine",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_BPN_NGML = list(
      description = "Buprenorphine plasma concentration driving the Emax effect on the logit of abstinence. Time-varying; supplied per record.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Methods, 'Modeling efficacy outcomes': individual predicted concentrations from the Jones 2021 BUP-XR / transmucosal popPK models (EBEs refitted to both studies). Set to 0 for placebo. In Study 2 the concentration feeds the `effect` compartment rather than the Emax term directly.",
      source_name = "C (Supplementary Material, fC = Emax model of BUP exposure)"
    ),
    STUDY_NCT04995029 = list(
      description = "Study 2 indicator: 1 = the 2021-2024 high-risk fentanyl-era study NCT04995029, 0 = the 2015-2016 phase 3 study NCT02357901 (Study 1).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Study 1, NCT02357901)",
      notes = "Selects the study-specific intercept and Emax (Table 2) and switches the Emax term from the plasma concentration (Study 1) to the effect-compartment concentration (Study 2, Results 'Opioid use model').",
      source_name = "Study 1 / Study 2"
    ),
    REGIMEN_MAINT_100MG = list(
      description = "Randomised to the 300/100-mg BUP-XR regimen (two 300-mg injections followed by 100-mg monthly maintenance doses); 0 = 300/300-mg regimen or placebo.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (300/300-mg regimen)",
      notes = "The intercept factor 0.764 applies only in Study 1 ('100-mg arm effect in Study 1', carried over from Laffont 2022); the model multiplies the indicator by (1 - STUDY_NCT04995029), so the Study 2 100-mg arm takes no intercept shift.",
      source_name = "100 mg (Table 2, theta alpha 100 mg Study 1)"
    ),
    EMPLOYED = list(
      description = "Employed at baseline (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not employed)",
      notes = "Multiplies the intercept alpha by 0.796 (Table 2); alpha is negative, so employment raises the drug-free probability of abstinence (Results).",
      source_name = "Employment"
    ),
    RACE_BLACK = list(
      description = "Black / African American race (1 = yes, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or other race)",
      notes = "Multiplies Emax by 0.843 ('Emax was 16% lower in Black/African American participants', Results).",
      source_name = "Black/AA"
    ),
    OPIOID_ROUTE_INJECTION = list(
      description = "Injection reported as the main route of opioid use at baseline (1 = yes, 0 = smoking, snorting, oral or other).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-injection main route)",
      notes = "Multiplies EC50 by 1.74 (Table 2 footnote: 'Injection use' means injection was reported as the main route of opioid use at baseline).",
      source_name = "Injection use"
    ),
    FENTANYL_DAILY = list(
      description = "Daily fentanyl use reported at baseline (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no daily fentanyl use)",
      notes = "Enters only through its interaction with AMPHETAMINE_USE: EC50 is multiplied by 8.25 when both are 1 (Table 2 'daily fentanyl and amphetamine/methamphetamine co-use at baseline'). Imputed as 0 for all Study 1 participants because fentanyl use was not assessed (Results, 'Patient characteristics').",
      source_name = "Daily fentanyl"
    ),
    AMPHETAMINE_USE = list(
      description = "Amphetamine / methamphetamine use reported at baseline (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no amphetamine/methamphetamine use)",
      notes = "Enters only through its interaction with FENTANYL_DAILY (see FENTANYL_DAILY).",
      source_name = "Amphetamines/Methamphetamines"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 989L,
    n_studies = 2L,
    n_observations = 19811L,
    age_range = "18-72 years (Study 1 median 38, range 19-64; Study 2 median 40, range 18-72)",
    age_mean = "Study 1 39.7 years (SD 11.1); Study 2 41.2 years (SD 10.6)",
    sex_female_pct = 37.9,
    race_ethnicity = c(White = 72.4, `Black/African American` = 23.7, Other = 3.9),
    disease_state = "Treatment-seeking adults with moderate-to-severe opioid use disorder (DSM-5); Study 2 additionally required high-risk opioid use (injection, >= 500 mg IV heroin equivalent, or any highly potent synthetic opioid).",
    dose_range = "Study 1: BUP-XR 300 mg x 6 (300/300), 300 mg x 2 then 100 mg x 4 (300/100), or placebo, every 4 weeks. Study 2: 300 mg at weeks 1 and 2, then 100 or 300 mg every 4 weeks x 8.",
    regions = "USA (Study 1); USA and Canada (Study 2)",
    co_medication = "Study 2 baseline: 49.4% daily fentanyl, 44.0% fentanyl >= 2 times/day, 60.0% amphetamine/methamphetamine, 29.2% daily fentanyl plus amphetamine/methamphetamine; fentanyl imputed absent in Study 1.",
    notes = "Table 1 (n = 489 Study 1, n = 500 Study 2). 19,811 weekly opioid-use observations (10,943 Study 1, 8,868 Study 2), Results 'Opioid use model'. Race and sex percentages pooled from the per-study Table 1 counts."
  )

  ini({
    # Intercepts on the logit of P(no opioid use), one per study (stratum suffix)
    logite0_study1 <- -4.28; label("Intercept alpha on the logit of abstinence, Study 1 (logit units)") # Table 2 'alpha Study 1' = -4.28 (RSE 9.7%)
    logite0_study2 <- -6.83; label("Intercept alpha on the logit of abstinence, Study 2 (logit units)") # Table 2 'alpha Study 2' = -6.83 (RSE 7.6%)
    e_regimen_maint_100mg_logite0 <- 0.764; label("Multiplicative factor on alpha for the 300/100-mg arm in Study 1 (unitless)") # Table 2 'theta alpha, 100 mg, Study 1' = 0.764 (RSE 8.5%)
    e_employed_logite0 <- 0.796; label("Multiplicative factor on alpha for employment at baseline (unitless)") # Table 2 'theta alpha, employment' = 0.796 (RSE 6.8%)

    # Emax relationship on the logit scale
    lemax_study1 <- log(6.68); label("Emax on the logit of abstinence, Study 1 (logit units)") # Table 2 'Emax, Study 1' = 6.68 (RSE 6.7%)
    lemax_study2 <- log(10.0); label("Emax on the logit of abstinence, Study 2 (logit units)") # Table 2 'Emax, Study 2' = 10.0 (RSE 3.8%)
    e_race_black_emax <- 0.843; label("Multiplicative factor on Emax for Black/African American participants (unitless)") # Table 2 'theta Emax, race' = 0.843 (RSE 4.8%)
    lec50 <- log(0.999); label("Buprenorphine concentration giving half of Emax (ng/mL)") # Table 2 'EC50' = 0.999 ng/mL (RSE 22.1%)
    e_opioid_route_injection_ec50 <- 1.74; label("Multiplicative factor on EC50 for injection as the main route of opioid use (unitless)") # Table 2 'theta EC50, injection' = 1.74 (RSE 35.7%)
    e_fentanyl_daily_amphetamine_use_ec50 <- 8.25; label("Multiplicative factor on EC50 for daily fentanyl plus amphetamine/methamphetamine co-use (unitless)") # Table 2 'theta EC50, fentanyl/amph' = 8.25 (RSE 18.5%)

    # Study 2 effect compartment
    lke0 <- log(0.000856); label("Effect-compartment equilibration rate constant, Study 2 (1/h)") # Table 2 'ke0, Study 2' = 0.000856 1/h (RSE 13.3%)

    # IIV: additive (normal) on the intercept, log-normal on EC50 and ke0.
    # Log-normal CV% converted as omega^2 = log(1 + CV^2), the convention of
    # the same group's Jones 2021 BUP-XR popPK model (see vignette).
    etalogite0 ~ 5.1529 # Table 2 'IIV alpha' SD = 2.27 (RSE 13.2%); variance 2.27^2
    etalec50 ~ 1.3672 # Table 2 'IIV EC50' = 171 CV% (RSE 14.1%); log(1 + 1.71^2)
    etalke0 ~ 1.1135 # Table 2 'IIV ke0' = 143 CV% (RSE 8.8%); log(1 + 1.43^2)
  })

  model({
    study2 <- STUDY_NCT04995029
    arm100_study1 <- REGIMEN_MAINT_100MG * (1 - study2)

    # Intercept: study-specific alpha scaled by the arm and employment factors
    # (Table 2 'coefficient on alpha'), plus the additive normal eta
    # (Supplementary Material: logit = alpha + fC + eta).
    alpha <- (logite0_study1 * (1 - study2) + logite0_study2 * study2) *
      e_regimen_maint_100mg_logite0^arm100_study1 *
      e_employed_logite0^EMPLOYED

    emax <- (exp(lemax_study1) * (1 - study2) + exp(lemax_study2) * study2) *
      e_race_black_emax^RACE_BLACK
    ec50 <- exp(lec50 + etalec50) *
      e_opioid_route_injection_ec50^OPIOID_ROUTE_INJECTION *
      e_fentanyl_daily_amphetamine_use_ec50^(FENTANYL_DAILY * AMPHETAMINE_USE)
    ke0 <- exp(lke0 + etalke0)

    # Effect compartment (Study 2 only; integrated in all subjects but read
    # only when STUDY_NCT04995029 = 1). Starts at 0 at the first injection.
    d/dt(effect) <- ke0 * (CP_BPN_NGML - effect)
    cdrive <- CP_BPN_NGML * (1 - study2) + effect * study2

    fc <- emax * cdrive / (ec50 + cdrive)
    logit_opioid_abstinence <- alpha + fc + etalogite0
    prob_opioid_abstinence <- expit(logit_opioid_abstinence)

    # Binary endpoint: rxode2 category 1 is no opioid use (paper Y = 0) and
    # category 2 is opioid use (Y = 1), so the category index is Y + 1.
    opioid_use ~ c(prob_opioid_abstinence)
  })
}
