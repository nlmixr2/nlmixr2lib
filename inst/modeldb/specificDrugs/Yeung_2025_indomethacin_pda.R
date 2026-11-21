Yeung_2025_indomethacin_pda <- function() {
  description <- paste(
    "Mixed-effects logistic-regression exposure-response model for patent",
    "ductus arteriosus (PDA) treatment failure after an indomethacin",
    "course (further indomethacin course or surgery needed vs no further",
    "treatment) in 18 preterm neonates at Children's Mercy Hospital.",
    "Linear predictor in the per-course indomethacin AUC (course start to",
    "12 h after the last dose, from the IND popPK model), gestational",
    "age, postnatal age, prior IVH prophylaxis, sex and race, plus a",
    "subject random intercept, as tabulated in Yeung 2025 Table 2. No",
    "predictor was statistically significant. Static per-course",
    "regression: no ODE, no dosing events.",
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
    time = "n/a (static per-course logistic regression; outcome assessed by echocardiogram after each course)",
    dosing = "n/a (no dose events; exposure enters through the AUC_INDO_COURSE covariate)",
    concentration = "prob_pda_failure (probability of PDA treatment failure after a course, 0-1)"
  )

  covariateData <- list(
    AUC_INDO_COURSE = list(
      description = "Individual model-predicted indomethacin AUC over one PDA treatment course, from the first dose of the course to 12 h after its last dose.",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Yeung 2025 Methods 'PopPK Modeling' and 'Evaluation of PDA",
        "Treatment Efficacy': simulated from empirical Bayes estimates of",
        "the published IND popPK model of Krzyzanski 2023 (with the",
        "CYP2C9 covariate), not dose-normalised. A course is three doses",
        "12 h apart, so the window is about 36 h. Figure 3 median roughly",
        "35-40 ug*h/mL."
      ),
      source_name = "AUCcourse"
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
      description = "Postnatal age at the course.",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 coefficient is per HOUR of postnatal age; model() recovers hours with PNA * 730.5. The paper does not state the reference time within a course; supply the value at the first dose of the course.",
      source_name = "Postnatal age (hours)"
    ),
    PRIOR_INDO = list(
      description = "Prior indomethacin IVH-prophylaxis indicator; 1 = the infant received indomethacin IVH prophylaxis after birth before PDA treatment, 0 = PDA treatment only.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PDA only)",
      notes = "Table 2 'Condition (0 = PDA only, 1 = IVH and PDA)'. Ten of the 18 PDA infants had prior IVH prophylaxis, median 8.16 (0.23-20.8) days before the first PDA dose (Results 'Patient Characteristics').",
      source_name = "Condition"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 2 'Sex (0 = male, 1 = female)'.",
      source_name = "Sex"
    ),
    RACE_BLACK = list(
      description = "Black race indicator; 1 = Black, 0 = non-Black.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (Black) is the paper's reference; the effect is on non-Black = 1 - RACE_BLACK",
      notes = "Table 2 'Race (0 = Black, 1 = non-Black)'; model() uses (1 - RACE_BLACK).",
      source_name = "Race"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18,
    n_studies = 1,
    age_range = "postnatal age at first PDA dose 3.1-27.4 days (median 11.2)",
    weight_range = "560-2630 g at first dose (median 835 g)",
    ga_range = "23-34 weeks gestational age at birth (median 25.9)",
    sex_female_pct = 33.3,
    race_ethnicity = c(Black = 38.9, White = 55.6, Other = 5.6),
    disease_state = "Preterm neonates with hemodynamically significant PDA treated with IV indomethacin (up to three courses)",
    dose_range = "0.1 mg/kg IV q12h x3 per course; dose could be increased after a failed course",
    regions = "USA (Children's Mercy Hospital, Kansas City)",
    notes = "Table 1 CMH PDA column. Final outcomes: closed 6, open without further treatment 6, open requiring surgery 6."
  )

  ini({
    logit_ref <- 2.79; label("Logit of PDA treatment failure at zero covariate values (unitless)") # Table 2, PDA block 'Intercept' = 2.79 (95% CI -8.77, 14.35)
    e_auc_course_pda <- -0.0065; label("Change in logit per ug*h/mL of course AUC") # Table 2, PDA block 'AUCcourse' = -0.0065 (95% CI -0.029, 0.016), p = 0.58
    e_ga_pda <- -0.23; label("Change in logit per week of gestational age (1/week)") # Table 2, PDA block 'Gestational age (weeks)' = -0.23 (95% CI -0.67, 0.22)
    e_pna_pda <- 0.009; label("Change in logit per hour of postnatal age (1/h)") # Table 2, PDA block 'Postnatal age (hours)' = 0.009 (95% CI -0.0015, 0.02)
    e_prior_indo_pda <- 0.24; label("Change in logit for prior IVH prophylaxis vs PDA only (unitless)") # Table 2, PDA block 'Condition' = 0.24 (95% CI -1.94, 2.42)
    e_sexf_pda <- 2.35; label("Change in logit for female vs male (unitless)") # Table 2, PDA block 'Sex' = 2.35 (95% CI -0.46, 5.16)
    e_race_nonblack_pda <- 0.99; label("Change in logit for non-Black vs Black (unitless)") # Table 2, PDA block 'Race' = 0.99 (95% CI -1.34, 3.34)

    etalogit_ref ~ fixed(0) # subject random intercept (Methods); variance not reported

    addSd_prob_pda_failure <- fixed(0.001); label("Placeholder additive residual SD on the probability; the source likelihood is Bernoulli") # not from source
  })

  model({
    pna_h <- PNA * 730.5
    logit_pda <- logit_ref + etalogit_ref +
      e_auc_course_pda * AUC_INDO_COURSE +
      e_ga_pda * GA +
      e_pna_pda * pna_h +
      e_prior_indo_pda * PRIOR_INDO +
      e_sexf_pda * SEXF +
      e_race_nonblack_pda * (1 - RACE_BLACK)
    prob_pda_failure <- expit(logit_pda)
    prob_pda_failure ~ add(addSd_prob_pda_failure)
  })
}
