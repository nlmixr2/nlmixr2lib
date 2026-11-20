Tanigawara_2012_garenoxacin_cure_cart <- function() {
  description <- paste(
    "Static landmark logistic exposure-response model for clinical cure with",
    "oral garenoxacin 400 mg once daily in Japanese patients with mild or",
    "moderate secondary infection of chronic respiratory disease (Tanigawara",
    "2012 phase III, n = 66), with the unbound AUC0-24/MIC ratio",
    "dichotomised at the CART-derived split of 28 h. The probability of cure",
    "is expit(logit_ref + e_faucmic28_success x I(fu x AUC_GARENOXACIN / MIC",
    "> 28)): 40.0% at or below the split and 91.8% above it. The",
    "continuous-index fit of the same data is",
    "Tanigawara_2012_garenoxacin_cure."
  )
  reference <- paste(
    "Tanigawara Y, Nozawa K, Tsuda H. Optimal dose finding of garenoxacin",
    "based on population pharmacokinetics/pharmacodynamics and Monte Carlo",
    "simulation. Eur J Clin Pharmacol. 2012;68(1):39-53.",
    "doi:10.1007/s00228-011-1095-3"
  )
  vignette <- "Tanigawara_2012_garenoxacin"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_GARENOXACIN covariate column)",
    concentration = "AUC_GARENOXACIN in ug*h/mL and MIC in ug/mL; prob_clinical_success is a probability (0-1)"
  )

  covariateData <- list(
    AUC_GARENOXACIN = list(
      description = "Steady-state total plasma garenoxacin AUC over the 24 h dosing interval (AUC0-24)",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical Bayes daily dose / CL/F from Tanigawara_2012_garenoxacin at",
        "400 mg once daily (Methods, 'PK/PD analysis of efficacy'). TOTAL",
        "concentration; the model multiplies by fu to form the unbound AUC."
      ),
      source_name = "AUC0-24"
    ),
    MIC = list(
      description = "Garenoxacin MIC against the causative pathogen isolated from the patient",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "CLSI method; when two or more pathogens were isolated the higher MIC was used (Methods).",
      source_name = "MIC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_studies = 1L,
    n_observations = "66 binary outcome records, one per patient",
    sex_female_pct = 34.6,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Mild or moderate secondary infection of chronic respiratory disease with a causative organism isolated (Table 2).",
    dose_range = "400 mg orally once daily for 10 days",
    regions = "Japan",
    notes = paste(
      "Tanigawara 2012 Tables 5-6: 61 patients above the split (56 cured) and",
      "5 at or below it (2 cured). sex_female_pct describes the 133-patient",
      "phase III PK set."
    )
  )

  ini({
    # Tanigawara 2012 Table 6, 'Logistic regression analysis by categorical
    # value'; the indicator is 1 for fAUC0-24/MIC > 28 h and 0 otherwise
    # (table footnote). The two coefficients reproduce Table 5 exactly:
    # expit(-0.4055) = 0.400 (2/5) and expit(-0.4055 + 2.8214) = 0.918
    # (56/61).
    logit_ref <- -0.4055; label("Logit of clinical cure at fAUC0-24/MIC <= 28 h (unitless logit)") # Table 6 categorical model: intercept -0.4055 (SE 0.9129)
    e_faucmic28_success <- 2.8214; label("Log-odds increment of clinical cure for fAUC0-24/MIC > 28 h (unitless logit)") # Table 6 categorical model: coefficient 2.8214 (SE 1.0253), odds ratio 16.80, p = 0.0059

    # CART split, a fixed breakpoint of the categorical regression rather
    # than one of its estimated coefficients.
    faucmic_split <- fixed(28); label("CART-derived fAUC0-24/MIC split (h)") # Table 6 footnote b and Results: 'an fAUC0-24/MIC ratio > 28 h was found to be a statistically significant (p = 0.0059) predictor'

    fu <- fixed(0.25); label("Unbound fraction of garenoxacin in serum (fraction)") # Discussion: 'The serum protein binding of garenoxacin is 75%'

    addSd_prob_clinical_success <- fixed(0.001); label("Placeholder additive residual SD on the clinical-cure probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })
  model({
    faucmic <- fu * AUC_GARENOXACIN / MIC
    above_split <- 0
    if (faucmic > faucmic_split) above_split <- 1

    logit_clinical_success <- logit_ref + e_faucmic28_success * above_split
    prob_clinical_success <- expit(logit_clinical_success)

    prob_clinical_success ~ add(addSd_prob_clinical_success)
  })
}
