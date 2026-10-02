Brendel_2021_irinotecan_liposomal_neutropenia_sn38 <- function() {
  description <- paste(
    "Logistic-regression exposure-safety model for CTCAE grade 3 or worse neutropenia adverse event (neutropenia, decreased neutrophil count or febrile neutropenia) at first event",
    "in adults with metastatic pancreatic ductal adenocarcinoma treated with",
    "liposomal irinotecan-based regimens (Brendel 2021; N = 316 from",
    "NAPOLI-1 and the phase I/II first-line study NCT02551991). The",
    "probability of an event is expit(a + b * log10(CAV)), where CAV is the",
    "individual average steady-state SN-38 plasma concentration in",
    "ng/mL. There is no PK layer and no ODE: exposure is supplied as data",
    "and was obtained in the source from post hoc estimates of the",
    "companion population PK model Brendel_2021_irinotecan_liposomal. The",
    "exposure effect is NOT statistically significant; reported by the source and kept here for completeness. The slope is the printed odds ratio; the",
    "intercept is not printed and was recovered by digitising the fitted",
    "curve in Figure 5b (see the vignette).",
    sep = " "
  )
  reference <- paste(
    "Brendel K, Bekaii-Saab T, Boland PM, Dayyani F, Dean A, Macarulla T,",
    "Maxwell F, Mody K, Pedret-Dunn A, Wainberg ZA, Zhang B. Population",
    "pharmacokinetics of liposomal irinotecan in patients with cancer and",
    "exposure-safety analyses in patients with metastatic pancreatic cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10:1550-1563.",
    "doi:10.1002/psp4.12725. PMCID: PMC8674005.",
    "Odds ratio from the Results; intercept digitised from Figure 5b.",
    sep = " "
  )
  vignette <- "Brendel_2021_irinotecan_liposomal"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_neutropenia_grade3 (probability of CTCAE grade 3 or worse neutropenia adverse event (neutropenia, decreased neutrophil count or febrile neutropenia) at first event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state SN-38 plasma concentration, Cavg,ss = AUCss,tau / tau. Supplied as data: this model has no PK layer.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Brendel 2021 Methods (Drug exposure derivation): noncompartmental",
        "analysis of steady-state concentrations simulated from each",
        "patient's post hoc parameters of the final population PK model",
        "(reproduce with modellib('Brendel_2021_irinotecan_liposomal')),",
        "divided by the patient's dosing interval (2 or 3 weeks). Enters on",
        "the LOG10 scale and uncentred. Median about 10^-0.025 = 0.94 ng/mL (Figure 5b quartile boundaries at log10 -0.12, -0.025 and 0.095)."
      ),
      source_name = "Cavg,ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 316L,
    n_studies = 2L,
    n_observations = "316 binary records (one per patient; earliest event only)",
    disease_state = "Metastatic pancreatic ductal adenocarcinoma: second line after gemcitabine (NAPOLI-1, N = 260) or first line (NCT02551991, N = 56)",
    dose_range = "Liposomal irinotecan 70 mg/m^2 free base Q2W with 5-FU/LV or 100 mg/m^2 Q3W alone (NAPOLI-1); 50, 55 or 70 mg/m^2 Q2W with 5-FU/LV and oxaliplatin (NCT02551991)",
    notes = paste(
      "Event rate: grade 3 or worse diarrhea at first event 17.7% overall",
      "(16.9% NAPOLI-1, 21.4% NCT02551991); grade 3 or worse neutropenia",
      "adverse events at first event 24.7% overall (20.8% and 42.9%).",
      "Neutropenia adverse events pooled neutropenia, decreased neutrophil",
      "count and febrile neutropenia (umbrella approach). Univariable",
      "regressions: no other covariate enters this model."
    )
  )

  ini({
    # Provenance. The slope is the printed odds ratio per log10 unit of
    # Cavg,ss (Results: OR 3.14, p = 0.115 (Figure 5b)). No intercept is printed; it
    # was DIGITISED from the fitted curve in Figure 5b with the slope held
    # at log(3.14). Self-check of the digitisation: fitting BOTH
    # coefficients to the traced curve returns -1.167 and 1.155 (odds ratio 3.17), reproducing the
    # printed odds ratio. Observed events per exposure quartile in Figure 5b:
    # 13/79, 21/78, 20/80 and 24/79 (overall 24.7%).
    logit_ref <- -1.165
    label("Logit of the event probability at CAV = 1 ng/mL (unitless logit)") # digitised from Brendel 2021 Figure 5b; no printed value exists
    e_log10cav <- log(3.14)
    label("Log-odds of the event per unit increase in log10 CAV (unitless logit)") # Brendel 2021 Results, odds ratio 3.14 per log10 unit of SN-38 Cavg,ss

    # Binomial logistic regression (R glm): no random effect and no residual
    # error are estimated. The placeholder additive residual exists only so
    # that rxode2 accepts an observation declaration.
    addSd_prob_neutropenia_grade3 <- fixed(0.001)
    label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # CAV is uncentred, so logit_ref is the logit at CAV = 1 ng/mL
    # (log10 CAV = 0); it is not a reference-patient probability.
    logit_neutropenia_grade3 <- logit_ref + e_log10cav * log10(CAV)
    prob_neutropenia_grade3 <- expit(logit_neutropenia_grade3)

    prob_neutropenia_grade3 ~ add(addSd_prob_neutropenia_grade3)
  })
}
