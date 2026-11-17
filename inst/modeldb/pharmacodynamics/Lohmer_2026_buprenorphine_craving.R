Lohmer_2026_buprenorphine_craving <- function() {
  description <- paste(
    "Proportional-odds (cumulative-logit) mixed-effects exposure-response",
    "model of current opioid craving on the 0-100 mm Opioid Craving Visual",
    "Analog Scale (OC-VAS), categorised as 0, 1-5, 6-20 and >20 mm, under",
    "monthly subcutaneous extended-release buprenorphine (BUP-XR,",
    "SUBLOCADE) in adults with moderate-to-severe opioid use disorder,",
    "integrated across a pre-fentanyl phase 3 study (Study 1, NCT02357901)",
    "and a high-risk fentanyl-era study (Study 2, NCT04995029).",
    "logit P(category <= m) = alpha_m + Emax * C / (EC50 + C) +",
    "theta_placebo * (1 - exp(-ke0 * t)) + eta, with study-specific",
    "intercept alpha_1, threshold increments and Emax, a Study 1",
    "300/100-mg arm factor on alpha_1, and Black/African American,",
    "fentanyl >= 2 times/day and fentanyl-plus-injection factors on Emax.",
    "Study 1 uses the plasma buprenorphine concentration directly; Study 2",
    "drives the drug effect through an effect compartment sharing the",
    "placebo-onset rate ke0 (half-life about 13 days). The model has no PK",
    "ODE: the buprenorphine plasma concentration CP_BPN_NGML is a required",
    "time-varying input (the authors used individual predictions from the",
    "Jones 2021 BUP-XR popPK model, packaged as",
    "modellib('Jones_2021_buprenorphine'))."
  )
  reference <- paste(
    "Lohmer LRL, Greenwald MK, Kakar R, Laffont CM. Optimal",
    "buprenorphine exposure to treat patients with opioid use disorder in",
    "the era of fentanyl and polysubstance use. Front Pharmacol.",
    "2026;17:1839067. doi:10.3389/fphar.2026.1839067.",
    "Parameter estimates from Table 3; model equations from the",
    "Supplementary Material ('Opioid Craving Model Equations').",
    "Sister models from the same paper:",
    "modellib('Lohmer_2026_buprenorphine_opioidUse'),",
    "modellib('Lohmer_2026_buprenorphine_dropout')."
  )
  vignette <- "Lohmer_2026_buprenorphine"

  depends <- c("CP_BPN_NGML")

  units <- list(
    time = "h",
    dosing = "N/A (PD-only; buprenorphine plasma concentration is a required input covariate)",
    concentration = "ng/mL (CP_BPN_NGML input covariate and the Study 2 effect compartment); the observation is an ordinal OC-VAS category"
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
      description = "Buprenorphine plasma concentration driving the Emax effect on the cumulative logits of craving. Time-varying; supplied per record.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Methods, 'Modeling efficacy outcomes': individual predicted concentrations from the Jones 2021 BUP-XR / transmucosal popPK models. Set to 0 for placebo. In Study 2 the concentration feeds the `effect` compartment rather than the Emax term directly.",
      source_name = "C (Supplementary Material, fC = Emax model of BUP exposure)"
    ),
    STUDY_NCT04995029 = list(
      description = "Study 2 indicator: 1 = the 2021-2024 high-risk fentanyl-era study NCT04995029, 0 = the 2015-2016 phase 3 study NCT02357901 (Study 1).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Study 1, NCT02357901)",
      notes = "Selects the study-specific intercept, threshold increments and Emax (Table 3) and switches the Emax term from the plasma concentration (Study 1) to the effect-compartment concentration (Study 2).",
      source_name = "Study 1 / Study 2"
    ),
    REGIMEN_MAINT_100MG = list(
      description = "Randomised to the 300/100-mg BUP-XR regimen (two 300-mg injections followed by 100-mg monthly maintenance doses); 0 = 300/300-mg regimen or placebo.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (300/300-mg regimen)",
      notes = "The factor 0.971 on alpha_1 applies only in Study 1 (Table 3); the model multiplies the indicator by (1 - STUDY_NCT04995029).",
      source_name = "100 mg (Table 3, theta alpha 100 mg Study 1)"
    ),
    RACE_BLACK = list(
      description = "Black / African American race (1 = yes, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or other race)",
      notes = "Multiplies Emax by 1.43 ('43% higher Emax in Black/African Americans', Results).",
      source_name = "Black/AA"
    ),
    FENTANYL_GE2_DAILY = list(
      description = "Fentanyl use at least twice per day reported at baseline (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fentanyl less than twice daily or none)",
      notes = "Multiplies Emax by 0.736 without, and by 0.198 with, injection as the main route of opioid use (OPIOID_ROUTE_INJECTION; Table 3). Imputed as 0 for all Study 1 participants because fentanyl use was not assessed.",
      source_name = "Fentanyl >= 2x/day"
    ),
    OPIOID_ROUTE_INJECTION = list(
      description = "Injection reported as the main route of opioid use at baseline (1 = yes, 0 = smoking, snorting, oral or other).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-injection main route)",
      notes = "Enters only together with FENTANYL_GE2_DAILY (Table 3 'fentanyl use >= 2 times/day and injection use at baseline'); no main effect of injection on craving was retained.",
      source_name = "Injection use"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 987L,
    n_studies = 2L,
    n_observations = 18093L,
    age_range = "18-72 years (Study 1 median 38, range 19-64; Study 2 median 40, range 18-72)",
    age_mean = "Study 1 39.7 years (SD 11.1); Study 2 41.2 years (SD 10.6)",
    sex_female_pct = 37.9,
    race_ethnicity = c(White = 72.4, `Black/African American` = 23.7, Other = 3.9),
    disease_state = "Treatment-seeking adults with moderate-to-severe opioid use disorder (DSM-5); Study 2 additionally required high-risk opioid use.",
    dose_range = "Study 1: BUP-XR 300 mg x 6 (300/300), 300 mg x 2 then 100 mg x 4 (300/100), or placebo, every 4 weeks. Study 2: 300 mg at weeks 1 and 2, then 100 or 300 mg every 4 weeks x 8.",
    regions = "USA (Study 1); USA and Canada (Study 2)",
    biomarkers = "Current craving on the 0-100 mm Opioid Craving VAS (0 = no craving at all, 100 = strongest craving ever), categorised 0, 1-5, 6-20, >20 mm; >20 mm is clinically significant.",
    notes = "Table 1 demographics (n = 489 Study 1, n = 500 Study 2; the craving analysis used 489 + 498). 18,093 craving observations (11,347 Study 1, 6,746 Study 2), Results 'Opioid craving model'."
  )

  ini({
    # Cumulative logit of the lowest category (0 mm), one per study
    logite0_study1 <- -5.39; label("Intercept alpha_1, cumulative logit of OC-VAS 0 mm, Study 1 (logit units)") # Table 3 'alpha1, Study 1' = -5.39 (RSE 7.2%)
    logite0_study2 <- -8.61; label("Intercept alpha_1, cumulative logit of OC-VAS 0 mm, Study 2 (logit units)") # Table 3 'alpha1, Study 2' = -8.61 (RSE 16%)
    e_regimen_maint_100mg_logite0 <- 0.971; label("Multiplicative factor on alpha_1 for the 300/100-mg arm in Study 1 (unitless)") # Table 3 'theta alpha 100 mg, Study 1' = 0.971 (RSE 4.2%)
    # Threshold increments: alpha_2 = alpha_1 + delta_1, alpha_3 = alpha_2 + delta_2
    logitstep1_study1 <- 2.42; label("Increment delta_1 from alpha_1 to alpha_2, Study 1 (logit units)") # Table 3 'delta1, Study 1' = 2.42 (RSE 4.3%)
    logitstep1_study2 <- 3.30; label("Increment delta_1 from alpha_1 to alpha_2, Study 2 (logit units)") # Table 3 'delta1, Study 2' = 3.30 (RSE 6.5%)
    logitstep2_study1 <- 1.92; label("Increment delta_2 from alpha_2 to alpha_3, Study 1 (logit units)") # Table 3 'delta2, Study 1' = 1.92 (RSE 4.9%)
    logitstep2_study2 <- 1.87; label("Increment delta_2 from alpha_2 to alpha_3, Study 2 (logit units)") # Table 3 'delta2, Study 2' = 1.87 (RSE 6.4%)

    # Emax relationship on the cumulative-logit scale
    lemax_study1 <- log(2.65); label("Emax on the cumulative logits of craving, Study 1 (logit units)") # Table 3 'Emax, Study 1' = 2.65 (RSE 13%)
    lemax_study2 <- log(3.55); label("Emax on the cumulative logits of craving, Study 2 (logit units)") # Table 3 'Emax, Study 2' = 3.55 (RSE 30%)
    lec50 <- log(0.365); label("Buprenorphine concentration giving half of Emax (ng/mL)") # Table 3 'EC50' = 0.365 ng/mL (RSE 66%)
    e_race_black_emax <- 1.43; label("Multiplicative factor on Emax for Black/African American participants (unitless)") # Table 3 'theta Emax, race' = 1.43 (RSE 7.2%)
    e_fentanyl_ge2_daily_emax <- 0.736; label("Multiplicative factor on Emax for fentanyl >= 2 times/day without injection use (unitless)") # Table 3 'theta Emax, fentanyl, no injection' = 0.736 (RSE 12%)
    e_fentanyl_ge2_daily_injection_emax <- 0.198; label("Multiplicative factor on Emax for fentanyl >= 2 times/day with injection use (unitless)") # Table 3 'theta Emax, fentanyl+injection' = 0.198 (RSE 115%)

    # Placebo effect and the shared onset / equilibration rate
    emax_placebo <- 3.73; label("Maximal placebo effect on the cumulative logits (logit units)") # Table 3 'theta Placebo' = 3.73 (RSE 7.4%)
    lke0 <- log(0.00217); label("Placebo-onset and Study 2 effect-compartment rate constant ke0 (1/h)") # Table 3 'ke0' = 0.00217 1/h (RSE 7.6%)

    # IIV: additive (normal) on the intercept, log-normal on Emax, EC50, ke0.
    # Log-normal CV% converted as omega^2 = log(1 + CV^2) (see vignette).
    etalogite0 ~ 2.9241 # Table 3 'IIV alpha' SD = 1.71 (RSE 8.5%); variance 1.71^2
    etalemax ~ 0.2465 # Table 3 'IIV Emax' = 53.0 CV% (RSE 29%); log(1 + 0.53^2)
    etalec50 ~ 0.9604 # Table 3 'IIV EC50' = 127 CV% (RSE 67%); log(1 + 1.27^2)
    etalke0 ~ 1.5282 # Table 3 'IIV ke0' = 190 CV% (RSE 4.4%); log(1 + 1.90^2)
  })

  model({
    study2 <- STUDY_NCT04995029
    arm100_study1 <- REGIMEN_MAINT_100MG * (1 - study2)

    alpha1 <- (logite0_study1 * (1 - study2) + logite0_study2 * study2) *
      e_regimen_maint_100mg_logite0^arm100_study1
    delta1 <- logitstep1_study1 * (1 - study2) + logitstep1_study2 * study2
    delta2 <- logitstep2_study1 * (1 - study2) + logitstep2_study2 * study2

    fent_noinj <- FENTANYL_GE2_DAILY * (1 - OPIOID_ROUTE_INJECTION)
    fent_inj <- FENTANYL_GE2_DAILY * OPIOID_ROUTE_INJECTION
    emax <- (exp(lemax_study1) * (1 - study2) + exp(lemax_study2) * study2) *
      exp(etalemax) *
      e_race_black_emax^RACE_BLACK *
      e_fentanyl_ge2_daily_emax^fent_noinj *
      e_fentanyl_ge2_daily_injection_emax^fent_inj
    ec50 <- exp(lec50 + etalec50)
    ke0 <- exp(lke0 + etalke0)

    # Effect compartment (Study 2 only; integrated in all subjects but read
    # only when STUDY_NCT04995029 = 1). Starts at 0 at the first injection.
    d/dt(effect) <- ke0 * (CP_BPN_NGML - effect)
    cdrive <- CP_BPN_NGML * (1 - study2) + effect * study2
    fc <- emax * cdrive / (ec50 + cdrive)

    # Placebo effect, printed for the final model without a study or arm
    # restriction (Supplementary Material); t is time since the first
    # injection.
    fplacebo <- emax_placebo * (1 - exp(-ke0 * t))

    lp <- fc + fplacebo + etalogite0
    logit_ocvas_le0 <- alpha1 + lp
    logit_ocvas_le5 <- alpha1 + delta1 + lp
    logit_ocvas_le20 <- alpha1 + delta1 + delta2 + lp
    prob_ocvas_le0 <- expit(logit_ocvas_le0)
    prob_ocvas_le5 <- expit(logit_ocvas_le5)
    prob_ocvas_le20 <- expit(logit_ocvas_le20)

    prob_ocvas0 <- prob_ocvas_le0
    prob_ocvas1_5 <- prob_ocvas_le5 - prob_ocvas_le0
    prob_ocvas6_20 <- prob_ocvas_le20 - prob_ocvas_le5
    prob_ocvas_gt20 <- 1 - prob_ocvas_le20

    # Ordinal endpoint: category 1 = 0 mm, 2 = 1-5 mm, 3 = 6-20 mm and the
    # unlisted category 4 = >20 mm (the paper's Y = 1..4).
    ocvas_cat ~ c(prob_ocvas0, prob_ocvas1_5, prob_ocvas6_20)
  })
}
