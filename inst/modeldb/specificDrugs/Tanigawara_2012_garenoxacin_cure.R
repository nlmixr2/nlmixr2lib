Tanigawara_2012_garenoxacin_cure <- function() {
  description <- paste(
    "Static landmark logistic exposure-response model for clinical cure with",
    "oral garenoxacin 400 mg once daily in Japanese patients with mild or",
    "moderate secondary infection of chronic respiratory disease (Tanigawara",
    "2012 phase III, n = 66 patients with an isolated pathogen and MIC). The",
    "probability of cure is expit(logit_ref + e_log10faucmic_success x",
    "log10(fu x AUC_GARENOXACIN / MIC)), i.e. logistic in the common",
    "logarithm of the unbound AUC0-24/MIC ratio. There is no PK layer: the",
    "total steady-state AUC0-24 is a covariate column that can be generated",
    "from Tanigawara_2012_garenoxacin as daily dose / individual CL/F, and",
    "the unbound fraction 0.25 (75% serum protein binding) is fixed."
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
        "Tanigawara 2012 Methods, 'PK/PD analysis of efficacy': 'The AUC0-24",
        "was calculated for each patient using the final population PK",
        "model', i.e. the empirical Bayes estimate of daily dose / CL/F from",
        "Tanigawara_2012_garenoxacin at 400 mg once daily. TOTAL",
        "concentration; the model multiplies by fu to form the unbound AUC.",
        "Population mean at 400 mg once daily 104 ug*h/mL (Figure 2)."
      ),
      source_name = "AUC0-24"
    ),
    MIC = list(
      description = "Garenoxacin MIC against the causative pathogen isolated from the patient",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Phase III MICs were determined by the CLSI method (Methods,",
        "'Determination of MICs'). Polymicrobial rule: 'When two or more",
        "infecting bacteria were detected in a patient, the higher MIC was",
        "used.' Table 2 lists the phase III isolates, MIC 0.05 to 32 ug/mL",
        "(most at or below 0.06 ug/mL)."
      ),
      source_name = "MIC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_studies = 1L,
    n_observations = "66 binary outcome records, one per patient",
    age_range = "19-94 years (the 133-patient phase III PK set)",
    weight_range = "33-120 kg (the 133-patient phase III PK set)",
    sex_female_pct = 34.6,
    race_ethnicity = c(Japanese = 100),
    disease_state = paste(
      "Mild or moderate secondary infection of chronic respiratory disease",
      "(acute exacerbation of chronic bronchitis, diffuse panbronchiolitis,",
      "bronchiectasis, emphysema and pulmonary fibrosis, bronchial asthma, old",
      "pulmonary tuberculosis) with a causative organism isolated: MSSA 8,",
      "PSSP 6, PISP/PRSP 8, H. influenzae 29, M. catarrhalis 7, K. pneumoniae",
      "4, P. aeruginosa 5, and single isolates of S. pyogenes and S.",
      "maltophilia plus two H. parainfluenzae (Table 2)."
    ),
    dose_range = "400 mg orally once daily for 10 days",
    regions = "Japan",
    notes = paste(
      "Tanigawara 2012 phase III. Clinical efficacy evaluated 7-14 days after",
      "the start of treatment; cure = 1, failure = 0. 58 of 66 patients were",
      "cured (87.9%). The age, weight and sex figures describe the 133",
      "patients in the phase III PK analysis (87 men), of whom these 66 are a",
      "subset."
    )
  )

  ini({
    # Tanigawara 2012 Table 6, 'Logistic regression analysis by continuous
    # value': logit(P(cure)) = -0.4565 + 1.0865 x log10(fAUC0-24/MIC), with
    # footnote a: 'This value was transformed (common logarithmic
    # transformation)'. Figure 5 confirms the base-10 log: the curve passes
    # 39% at fAUC/MIC = 1 h, 75% at the 28 h CART split and 94% at 1000 h.
    logit_ref <- -0.4565; label("Logit of clinical cure at fAUC0-24/MIC = 1 h, i.e. log10 ratio of zero (unitless logit)") # Table 6 continuous model: intercept -0.4565 (SE 1.0317)
    e_log10faucmic_success <- 1.0865; label("Log-odds of clinical cure per unit increase in log10(fAUC0-24/MIC) (unitless logit)") # Table 6 continuous model: coefficient 1.0865 (SE 0.4526), odds ratio 2.964, p = 0.0164

    # Garenoxacin serum protein binding 75% (Discussion, citing Van Wart et
    # al.), so fu = 0.25; the paper uses it to turn the AUC0-24/MIC target of
    # 125 h into the fAUC0-24/MIC target of 30 h.
    fu <- fixed(0.25); label("Unbound fraction of garenoxacin in serum (fraction)") # Discussion: 'The serum protein binding of garenoxacin is 75%'

    # The source likelihood is Bernoulli; the placeholder additive residual
    # only satisfies rxode2's requirement for an observation declaration.
    addSd_prob_clinical_success <- fixed(0.001); label("Placeholder additive residual SD on the clinical-cure probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })
  model({
    # PK/PD index in hours: unbound AUC0-24 (ug*h/mL) over MIC (ug/mL).
    faucmic <- fu * AUC_GARENOXACIN / MIC

    logit_clinical_success <- logit_ref + e_log10faucmic_success * log10(faucmic)
    prob_clinical_success <- expit(logit_clinical_success)

    # Deterministic probability; sample outcomes downstream with
    # rbinom(n, 1, prob_clinical_success).
    prob_clinical_success ~ add(addSd_prob_clinical_success)
  })
}
