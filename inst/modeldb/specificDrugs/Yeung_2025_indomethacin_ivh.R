Yeung_2025_indomethacin_ivh <- function() {
  description <- paste(
    "Logistic-regression exposure-response model for severe",
    "intraventricular hemorrhage (IVH Grade III or IV vs none/I/II, by",
    "cranial ultrasound in the first week of life) in 75 preterm neonates",
    "given indomethacin IVH prophylaxis at two US NICUs. Linear predictor",
    "in the dose-normalised indomethacin AUC0-inf (from the CYP2C9-reduced",
    "popPK models), study site, gestational age, postnatal age at the",
    "first dose, sex and race, as tabulated in Yeung 2025 Table 2. Only",
    "gestational age was statistically significant; exposure had no",
    "detectable effect. Static per-subject regression: no ODE, no dosing",
    "events, no time dimension.",
    sep = " "
  )
  reference <- paste(
    "Yeung CHT, Sekulich DC, Scott A, Nolte WM, Gibson K, Su R,",
    "Alrifai MW, Lopata SM, Lewis T, Reese J, Maharaj A (2025).",
    "The Relationship of Indomethacin Exposure With Efficacy and Renal",
    "Toxicity Outcomes for Preterm Infants in the Neonatal Intensive Care",
    "Unit. Clin Transl Sci 18(5):e70251. doi:10.1111/cts.70251.",
    sep = " "
  )
  vignette <- "Yeung_2025_indomethacin"
  units <- list(
    time = "n/a (static per-subject landmark logistic regression; outcome is the worst IVH grade by ultrasound in the first week of life)",
    dosing = "n/a (no dose events; exposure enters through the AUC_INDO covariate supplied by Yeung_2025_indomethacin_cmh or Yeung_2025_indomethacin_vumc)",
    concentration = "prob_ivh_severe (probability of severe IVH, 0-1)"
  )

  covariateData <- list(
    AUC_INDO = list(
      description = "Individual model-predicted indomethacin AUC0-inf of the IVH-prophylaxis course, normalised by the total course dose in mg/kg.",
      units = "ug*h/mL per mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Yeung 2025 Methods 'Evaluation of IVH Prevention Efficacy': AUC0-inf",
        "from the CYP2C9-reduced popPK models (empirical Bayes estimates),",
        "divided by 0.3 mg/kg (CMH, 0.1 mg/kg x 3) or 0.2 mg/kg (VUMC). Total",
        "plasma concentration. Figure 2 medians are roughly 230 (CMH) and",
        "360 (VUMC) ug*h/mL per mg/kg."
      ),
      source_name = "AUC0-inf (dose-normalised)"
    ),
    STUDY_VUMC = list(
      description = "Study-site indicator; 1 = Vanderbilt University Medical Center cohort, 0 = Children's Mercy Hospital cohort.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CMH)",
      notes = "Table 2 'Site (0 = CMH, 1 = VUMC)'. The two sites also differ in regimen (VUMC single 0.2 mg/kg at 12 h of life; CMH 0.1 mg/kg q24h x3).",
      source_name = "Site"
    ),
    GA = list(
      description = "Gestational age at birth.",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters linearly and uncentred (Table 2).",
      source_name = "Gestational age"
    ),
    PNA = list(
      description = "Postnatal age at the first indomethacin dose.",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 coefficient is per HOUR of postnatal age at the first dose; model() recovers hours with PNA * 730.5 (1 month = 30.4375 days). Supply the value at the first dose, not a time-varying column.",
      source_name = "Postnatal age at first dose (h)"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 2 'Sex (0 = male, 1 = female)'; same orientation as the canonical.",
      source_name = "Sex"
    ),
    RACE_BLACK = list(
      description = "Black race indicator; 1 = Black, 0 = non-Black.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (Black) is the paper's reference; the effect is on non-Black = 1 - RACE_BLACK",
      notes = "Table 2 'Race (0 = Black, 1 = non-Black)'. White, other and not-reported were collated as non-Black (Table 1 footnote). model() uses (1 - RACE_BLACK).",
      source_name = "Race"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 75,
    n_studies = 2,
    age_range = "postnatal age at first dose 0.23-53.7 h (CMH) and 11.4-20.3 h (VUMC)",
    weight_range = "360-1290 g (CMH), 370-1040 g (VUMC) at first dose",
    ga_range = "22.1-29.7 weeks gestational age at birth",
    sex_female_pct = 46.7,
    race_ethnicity = c(Black = 42.7, White = 53.3, Other = 2.7, `Not reported` = 1.3),
    disease_state = "Preterm neonates receiving indomethacin IVH prophylaxis; 9 of 75 developed severe (Grade III/IV) IVH",
    dose_range = "CMH 0.1 mg/kg IV q24h x3; VUMC 0.2 mg/kg IV single dose",
    regions = "USA (Children's Mercy Hospital; Vanderbilt University Medical Center)",
    notes = "Table 1 CMH IVH (N = 45) and VUMC IVH (N = 30) columns pooled; percentages recomputed from the counts."
  )

  ini({
    logit_ref <- 14.5; label("Logit of severe IVH at zero covariate values (unitless)") # Table 2, IVH block 'Intercept' = 14.5 (95% CI 0.34, 31.8)
    e_auc_ivh <- -0.0005; label("Change in logit per ug*h/mL per mg/kg of dose-normalised AUC0-inf") # Table 2, IVH block 'AUC0-inf' = -0.0005 (95% CI -0.0037, 0.0003), p = 0.51
    e_study_vumc_ivh <- -1.79; label("Change in logit for the VUMC site vs CMH (unitless)") # Table 2, IVH block 'Site' = -1.79 (95% CI -4.17, 0.28)
    e_ga_ivh <- -0.66; label("Change in logit per week of gestational age (1/week)") # Table 2, IVH block 'Gestational age (weeks)' = -0.66 (95% CI -1.36, -0.11), p = 0.034
    e_pna_ivh <- 0.087; label("Change in logit per hour of postnatal age at first dose (1/h)") # Table 2, IVH block 'Postnatal age at first dose (h)' = 0.087 (95% CI -0.090, 0.23)
    e_sexf_ivh <- -0.78; label("Change in logit for female vs male (unitless)") # Table 2, IVH block 'Sex' = -0.78 (95% CI -2.57, 0.84)
    e_race_nonblack_ivh <- 0.91; label("Change in logit for non-Black vs Black (unitless)") # Table 2, IVH block 'Race' printed -0.91 with 95% CI -0.67, 2.74; sign taken from the CI (see vignette)

    addSd_prob_ivh_severe <- fixed(0.001); label("Placeholder additive residual SD on the probability; the source likelihood is Bernoulli") # not from source
  })

  model({
    pna_h <- PNA * 730.5
    logit_ivh <- logit_ref +
      e_auc_ivh * AUC_INDO +
      e_study_vumc_ivh * STUDY_VUMC +
      e_ga_ivh * GA +
      e_pna_ivh * pna_h +
      e_sexf_ivh * SEXF +
      e_race_nonblack_ivh * (1 - RACE_BLACK)
    prob_ivh_severe <- expit(logit_ivh)
    prob_ivh_severe ~ add(addSd_prob_ivh_severe)
  })
}
