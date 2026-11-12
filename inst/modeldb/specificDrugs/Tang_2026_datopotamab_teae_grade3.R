Tang_2026_datopotamab_teae_grade3 <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for any grade 3 or worse ",
    "treatment-emergent adverse event (TEAE) in adults with breast or ",
    "non-small cell lung cancer treated with datopotamab deruxtecan ",
    "(Dato-DXd) 0.27-10 mg/kg every 3 weeks (Tang 2026, n = 1081 pooled ",
    "from TROPION-PanTumor01, -Lung01, -Lung05 and -Breast01). The ",
    "probability is expit(-1.171 + 1.305 * CAV + region shifts + 0.454 ",
    "* ever-smoker), where CAV is the individual average plasma ",
    "concentration of the released PAYLOAD DXd (ng/mL) from the first ",
    "dose to the end of the cycle of the event (or of the last ",
    "treatment for patients without an event). United States ",
    "enrollment and never-smoking are the reference levels. There is ",
    "no PK layer and no ODE: exposure is supplied as data and was ",
    "derived in the source from the companion population PK model ",
    "packaged as Hong_2025_datopotamab. No between-subject random effect ",
    "and no residual error are estimated. One of ten exposure-response ",
    "models in the Tang_2026_datopotamab_* family."
  )
  reference <- paste(
    "Tang Z, Lim K, Jiang Y, Denduluri N, Rokutanda N, Pan Y, Hong Y, Ren S,",
    "Vajjah P, Zhou D. Exposure-Response Analyses of Datopotamab Deruxtecan",
    "(Dato-DXd) in Patients with Hormone Receptor (HR) Positive, Human",
    "Epidermal Growth Factor Receptor 2 (HER2)-Negative Breast Cancer.",
    "J Clin Pharmacol. 2026;66(4):245-258. doi:10.1002/jcph.70188.",
    "Coefficients are transcribed from Supporting Information Table S4.",
    sep = " "
  )
  vignette <- "Tang_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_teae_grade3 (probability of any grade 3 or worse TEAE, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average plasma concentration of the released PAYLOAD DXd, cumulative AUC divided by the time from the first dose to the end of the cycle of the event (or, without an event, to the cycle of the last treatment or the data cut-off, whichever came first)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "PAYLOAD (DXd), not the conjugate: the unit is ng/mL, and the",
        "companion Tang_2026_datopotamab_stomatitis / _ose / _ose_grade2",
        "models carry an intact Dato-DXd Cavg in ug/mL under the same",
        "column name. Enters LINEARLY and UNCENTRED. Derived in the source",
        "from empirical Bayes estimates of the Hong 2025 population PK",
        "model, accounting for dose interruptions and reductions; it is",
        "endpoint-specific because the averaging window ends at the event.",
        "Reproduce with modellib('Hong_2025_datopotamab'), whose DXd output",
        "Cc_dxd is in ug/mL (multiply by 1000). Tang 2026 Figure 2A shows",
        "the analysis range spanning roughly 0.3-1.4 ng/mL with medians",
        "near 0.5, 0.7 and 1.0 ng/mL at 4, 6 and 8 mg/kg."
      ),
      source_name = "DXd Cavg"
    ),
    REGION_JAPAN = list(
      description = "Japan region indicator (1 = enrolled in Japan, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is United States enrollment (all REGION_* indicators 0)",
      notes = "One of four non-reference levels of the geographical-region factor (Japan, Europe, mainland China, rest of world vs USA). Tang 2026 Table S2: 192 of 1081 (17.8%).",
      source_name = "Region - Japan"
    ),
    REGION_EUROPE = list(
      description = "Europe region indicator (1 = enrolled in Europe, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is United States enrollment",
      notes = "Tang 2026 Table S2: 301 of 1081 (27.8%).",
      source_name = "Region - Europe"
    ),
    REGION_CHINA = list(
      description = "Mainland China region indicator (1 = enrolled in mainland China, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is United States enrollment",
      notes = "Tang 2026 Table S2 'Mainland China': 49 of 1081 (4.5%). Taiwan and Hong Kong are not in this level; they fall in rest of world.",
      source_name = "Region - China"
    ),
    REGION_ROW = list(
      description = "Rest-of-world region indicator (1 = enrolled outside the USA, Japan, Europe and mainland China, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is United States enrollment",
      notes = "Tang 2026 Table S2 'Rest of the World': 218 of 1081 (20.2%); includes Korea, Taiwan, Hong Kong, Singapore, Canada, Australia, Latin America, India, South Africa and the Russian Federation. Differs from Hong 2025, whose rest of world also includes mainland China.",
      source_name = "Region - Rest of world"
    ),
    SMOKE_NEVER = list(
      description = "Never-smoker indicator (1 = never smoked, 0 = former or current smoker)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (never smoker) is the reference level of this model; the coefficient applies to the ever-smoker indicator 1 - SMOKE_NEVER",
      notes = "Tang 2026 dichotomises smoking as 'former or current smokers' versus 'patients who never smoked' (odds ratio 1.57), so the ever-smoker indicator is formed in model() as 1 - SMOKE_NEVER, as in the unpaired SMOKE_NEVER usage of the covariate register. Table S2: never 544 (50.3%), current 437, former 100.",
      source_name = "Former/Current Smoker"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1081L,
    n_studies = 4L,
    age_range = "26-86 years (median 60)",
    weight_range = "35.6-156 kg (median 64.2)",
    sex_female_pct = 67.2,
    race_ethnicity = "47.2% White, 40% Asian, 1.8% Black, 7.3% multiple or other, 2.9% unknown or not reported (Table S2)",
    disease_state = "advanced or metastatic non-small cell lung cancer (59.5%), HR-positive/HER2-negative breast cancer (36.4%) or triple-negative breast cancer (4.1%)",
    dose_range = "datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; 919 of 1081 (85%) at 6 mg/kg, 50 at 4 mg/kg and 82 at 8 mg/kg",
    regions = "USA 29.7%, Europe 27.8%, Japan 17.8%, mainland China 4.5%, rest of world 20.2%",
    notes = "Safety analysis set: patients with at least one Dato-DXd dose and at least one post-baseline PK record in TROPION-PanTumor01, TROPION-Lung01, TROPION-Lung05 and TROPION-Breast01 (data cut-off 29 April 2024). Observed grade 3 or worse TEAE incidence 30%, 41% and 60% at 4, 6 and 8 mg/kg (Table 2)."
  )

  ini({
    # ==================================================================
    # Tang 2026 Supporting Information Table S4, 'TEAE Gr3+' rows. Every
    # value is a printed Estimate on the log-odds scale, so nothing is
    # back-solved. The printed odds ratios in the Results are an internal
    # check: exp(-0.7647521) = 0.465, exp(-0.3995153) = 0.671,
    # exp(-0.3466507) = 0.707, exp(-0.6327399) = 0.531 and
    # exp(0.4538759) = 1.574 against the printed 0.465, 0.671, 0.707,
    # 0.531 and 1.57.
    # ==================================================================
    logit_ref <- -1.1706039; label("Logit of the event probability at zero DXd exposure for a never-smoker enrolled in the USA (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ (Intercept): -1.1706039 [-1.59; -0.755]
    e_cav_logit <- 1.3054368; label("Log-odds of a grade 3 or worse TEAE per 1 ng/mL increase in DXd Cavg (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ DXd Cavg: 1.3054368 [0.917; 1.69]
    e_region_japan_logit <- -0.6327399; label("Log-odds shift for Japan vs USA enrollment (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ Region - Japan: -0.6327399 [-1.02; -0.250]
    e_region_europe_logit <- -0.3466507; label("Log-odds shift for Europe vs USA enrollment (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ Region - Europe: -0.3466507 [-0.680; -0.0130]
    e_region_china_logit <- -0.3995153; label("Log-odds shift for mainland China vs USA enrollment (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ Region - China: -0.3995153 [-1.05; 0.250]
    e_region_row_logit <- -0.7647521; label("Log-odds shift for rest-of-world vs USA enrollment (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ Region - Rest of world: -0.7647521 [-1.14; -0.392]
    e_smoke_ever_logit <- 0.4538759; label("Log-odds shift for former or current smokers vs never smokers (unitless logit)") # Tang 2026 Table S4, TEAE Gr3+ Former/Current Smoker: 0.4538759 [0.194; 0.714]

    # Bernoulli-likelihood logistic regression (R 4.3): no variance
    # components. The placeholder additive term exists only because
    # rxode2 requires an observation declaration.
    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor. CAV enters uncentred, so logit_ref is the logit
    # at zero payload exposure for the reference patient (USA, never
    # smoker). The four REGION_* indicators are mutually exclusive; all
    # zero selects the USA.
    logit_teae_grade3 <- logit_ref +
      e_cav_logit * CAV +
      e_region_japan_logit * REGION_JAPAN +
      e_region_europe_logit * REGION_EUROPE +
      e_region_china_logit * REGION_CHINA +
      e_region_row_logit * REGION_ROW +
      e_smoke_ever_logit * (1 - SMOKE_NEVER)

    prob_teae_grade3 <- expit(logit_teae_grade3)

    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
