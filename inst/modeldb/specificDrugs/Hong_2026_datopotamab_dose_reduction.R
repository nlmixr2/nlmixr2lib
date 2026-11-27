Hong_2026_datopotamab_dose_reduction <- function() {
  description <- paste(
    "Logistic-regression exposure-safety model for treatment-emergent",
    "adverse events (TEAEs) leading to dose reduction in adults with",
    "advanced non-small cell lung cancer or breast cancer receiving",
    "datopotamab deruxtecan (Dato-DXd) 0.27-10 mg/kg IV every 3 weeks (Hong",
    "2026, n = 729 pooled from TROPION-PanTumor01, -Lung05 and -Lung01).",
    "logit(P) = -3.263 + 0.1081 * AUC1, where AUC1 is the individual",
    "Cycle-1 intact Dato-DXd AUC in mg*h/mL (supplied as AUC_DATODXD in",
    "ug*day/mL and converted inside the model). This is the paper's",
    "DRUG-EFFECT model (exposure only), whose curve is plotted in Figure 4:",
    "the intercept and slope are back-solved from the digitised curve",
    "because no coefficient is printed. The final model also carries",
    "enrolment in Europe (odds ratio 2.20) and any renal impairment (1.90),",
    "but its intercept and slope are not published, so those effects are",
    "not encoded (see the vignette). There is no PK layer: exposure is",
    "supplied as data and was derived in the source from the companion",
    "population PK model packaged as Hong_2025_datopotamab. One of six",
    "exposure-response models in the Hong_2026_datopotamab_* family."
  )
  reference <- paste(
    "Hong Y, Mc Laughlin AM, Pan Y, Hietala SF, Peigne S, Chargualaf M,",
    "Tang Z, Zhou D, Tajima N, Abutarif M, Garimella T. Exposure-Response",
    "Analyses of Datopotamab Deruxtecan (Dato-DXd), a TROP2-Directed",
    "Antibody-Drug Conjugate, in Advanced Non-small Cell Lung Cancer.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70477. PMID 42733211.",
    "Intercept and slope back-solved from Figure 4 (bottom panel); model",
    "form from Supplementary Equation S2.",
    sep = " "
  )
  vignette <- "Hong_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_DATODXD covariate column)",
    concentration = "prob_dose_reduction (probability of a TEAE leading to dose reduction, 0-1)"
  )

  covariateData <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 (day 0 to day 21) area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Hong 2026 reports this exposure (AUC1) in mg/mL*h; the model",
        "converts with AUC_DATODXD * 24 / 1000. Enters LINEARLY and",
        "UNCENTRED. Figure 4 shows dose-group medians near 10, 15 and 21",
        "mg*h/mL (about 420, 625 and 875 ug*day/mL) at 4, 6 and 8 mg/kg.",
        "Derived in the source from empirical Bayes estimates of the Hong",
        "2025 population PK model, not observed."
      ),
      source_name = "Dato-DXd AUC1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 729L,
    n_studies = 3L,
    disease_state = "Advanced or metastatic non-small cell lung cancer (n = 644) or breast cancer (n = 85, TROPION-PanTumor01 only)",
    dose_range = "Datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; most patients at 6 mg/kg",
    regions = "USA, Europe, Japan, China and rest of world",
    notes = "Exposure-safety dataset: the 644-patient exposure-efficacy dataset (baseline characteristics in Table S3) plus 85 breast-cancer patients from TROPION-PanTumor01. Observed events by Dato-DXd AUC1 quartile 18/183, 28/182, 37/182 and 57/182 (Figure 4)."
  )

  ini({
    # Supplementary Equation S2 with a linear f(Exposure) for this endpoint
    # (Results). No coefficient is printed; nonlinear least squares of
    # expit(a + b * AUC1) on the 335 pixel columns of the Figure 4 bottom
    # model-predicted median gives the values below, with a largest
    # residual of 0.004 on the probability scale.
    logite0 <- -3.263; label("Logit of the dose-reduction TEAE probability at zero Dato-DXd exposure (unitless logit)") # back-solved from the Figure 4 bottom drug-effect curve
    e_auc_datodxd_logit <- 0.1081; label("Change in the logit per unit of Cycle-1 Dato-DXd AUC (1/(mg*h/mL))") # back-solved from the Figure 4 bottom drug-effect curve

    # The source likelihood is Bernoulli with no residual error; the
    # placeholder below only lets the nlmixr2 observation machinery accept
    # the model. Not from the source.
    addSd_prob_dose_reduction <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Cycle-1 AUC in the paper's unit: ug*day/mL * 24 h/day / 1000 ug/mg.
    auc1_mghml <- AUC_DATODXD * 24 / 1000
    logit_dose_reduction <- logite0 + e_auc_datodxd_logit * auc1_mghml
    prob_dose_reduction <- expit(logit_dose_reduction)
    prob_dose_reduction ~ add(addSd_prob_dose_reduction)
  })
}
