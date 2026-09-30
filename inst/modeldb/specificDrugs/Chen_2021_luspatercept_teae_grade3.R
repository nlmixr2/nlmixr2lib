Chen_2021_luspatercept_teae_grade3 <- function() {
  description <- paste0(
    "Landmark binomial logistic-regression exposure-SAFETY model relating the steady-state luspatercept AUC of the starting dose ",
    "to the probability of any treatment-emergent adverse event (TEAE) of grade 3 or higher during the first 2 treatment cycles in adults with beta-thalassemia ",
    "(Chen 2021, Figure 4A; n = 285 luspatercept-treated patients of the population PK analysis). ",
    "prob = expit(logit_ref + e_auc_logit * AUCss). Neither coefficient is printed: both are digitized from the fitted line in ",
    "Figure 4A, which is drawn as a vector path in the PDF. The paper describes the relationship as flat. There is no PK layer ",
    "and no ODE: exposure is supplied as the AUC_LUSP data column (starting dose / individual CL/F from Chen_2021_luspatercept). ",
    "No random effects and no residual error are estimated (Bernoulli likelihood). One of four exposure-response models in the ",
    "Chen_2021_luspatercept_* family."
  )
  reference <- paste(
    "Chen N, Kassir N, Laadem A, Giuseppi AC, Shetty J, Maxwell SE, Sriraman P, Ritland S, Linde PG, Budda B, Reynolds JG, Zhou S, Palmisano M.",
    "Population Pharmacokinetics and Exposure-Response Relationship of Luspatercept, an Erythroid Maturation Agent, in Anemic Patients With beta-Thalassemia.",
    "J Clin Pharmacol. 2021;61(1):52-63. doi:10.1002/jcph.1696. PMCID: PMC7754485.",
    "Logistic regression fit drawn in Figure 4A ('TEAES >= grade 3: cycles 1-2'); Results, 'Exposure-Response for TEAEs'.",
    sep = " "
  )
  vignette <- "Chen_2021_luspatercept"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_LUSP covariate column)",
    concentration = "prob_teae_grade3 (probability of any treatment-emergent adverse event (TEAE) of grade 3 or higher during treatment cycles 1-2, 0-1; also logit_teae_grade3)"
  )

  covariateData <- list(
    AUC_LUSP = list(
      description = "Individual luspatercept steady-state area under the serum concentration-time curve over the 21-day dosing interval for the starting dose (AUCss). Supplied as data: this model has no PK layer.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Chen 2021 Methods: AUC = dose / individual CL/F from the final population PK model (Chen_2021_luspatercept). During cycles 1-2 all patients were on their starting dose, so the dose-interval AUC at the event (AUC_TEAE) equals the starting-dose AUCss plotted on the Figure 4A axis. Untransformed and uncentred; the intercept is an extrapolation to zero exposure because placebo patients were shown for comparison but excluded from the fit (Methods). Observed range in Figure 4A about 13-301 ug*day/mL.",
      source_name = "AUCss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 285L,
    n_studies = 2L,
    n_observations = "285 binary event records, one per patient",
    age_range = "18-66 years",
    weight_range = "34.1-97.0 kg",
    sex_female_pct = 56.8,
    disease_state = "Adults with beta-thalassemia requiring regular red-blood-cell transfusions (the safety population equals the PK population).",
    dose_range = "Starting dose 0.2-1.25 mg/kg subcutaneously q3w (79.8% started at 1 mg/kg).",
    regions = "Phase 2 study A536-04 (NCT01749540) and phase 3 study ACE-536-B-THAL-001 BELIEVE (NCT02604433).",
    notes = "Observed events by AUCss quartile (Figure 4A): 7/71, 4/71, 5/71 and 7/72 (23/285 overall); placebo 1/109 (not part of the fit). The paper final exposure-safety analysis for this endpoint over all cycles is a Cox proportional-hazards model (Table 4) whose baseline hazard is not reported, so it is not packaged."
  )

  ini({
    # Chen 2021 Figure 4A: univariate logistic regression
    #
    #   logit(p) = logit_ref + e_auc_logit * AUCss
    #
    # Neither coefficient is printed (no odds ratio in the panel, text or
    # supplement). Both are recovered from the fitted line, which is a
    # vector path in the PDF: its two end nodes (AUCss 13.29 and 300.8)
    # read against the axis tick marks give p = 0.0861 and 0.0819, and the
    # logistic through those two points has the coefficients below. See the
    # vignette source trace.
    logit_ref <- -2.360; label("Logit of the probability of a grade >= 3 TEAE in cycles 1-2 extrapolated to zero luspatercept exposure (unitless logit)")  # digitized from Figure 4A fitted line (not printed); see vignette
    e_auc_logit <- -0.0001916; label("Log-odds of a grade >= 3 TEAE in cycles 1-2 per 1 ug*day/mL of luspatercept AUCss (unitless logit)")  # digitized from Figure 4A fitted line (not printed); see vignette

    # No between-subject variability and no residual error in the source
    # (Bernoulli likelihood). The tiny fixed additive residual exists only so
    # rxode2 has an error model to attach to the typical-value probability.
    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_teae_grade3 <- logit_ref + e_auc_logit * AUC_LUSP
    prob_teae_grade3 <- expit(logit_teae_grade3)
    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
