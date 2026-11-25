Hong_2026_datopotamab_teae_grade3 <- function() {
  description <- paste(
    "Logistic-regression exposure-safety model for treatment-related grade",
    "3 or worse treatment-emergent adverse events (TEAEs) in adults with",
    "advanced non-small cell lung cancer or breast cancer receiving",
    "datopotamab deruxtecan (Dato-DXd) 0.27-10 mg/kg IV every 3 weeks (Hong",
    "2026, n = 729 pooled from TROPION-PanTumor01, -Lung05 and -Lung01).",
    "logit(P) = -1.092 + 1.441 * CAV, where CAV is the individual average",
    "plasma concentration of the released PAYLOAD DXd (ng/mL) to the end of",
    "the cycle that included the event. This is the paper's DRUG-EFFECT",
    "model (exposure only), whose curve is plotted in Figure 4: the",
    "intercept and slope are back-solved from the digitised curve because",
    "no coefficient is printed. The final model also carries Asian race",
    "(odds ratio 0.566) and ECOG performance status > 0 (odds ratio 1.50),",
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
    "Intercept and slope back-solved from Figure 4 (top-left panel); model",
    "form from Supplementary Equation S2.",
    sep = " "
  )
  vignette <- "Hong_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_teae_grade3 (probability of a treatment-related grade 3 or worse TEAE, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average plasma concentration of the released PAYLOAD DXd from the first dose to the end of the cycle that included the event",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "PAYLOAD (DXd), not the conjugate: the companion",
        "Hong_2026_datopotamab_orr model carries an intact Dato-DXd average",
        "concentration in ug/mL under the same column name. Enters LINEARLY",
        "and UNCENTRED. Figure 4 shows the analysis range spanning roughly",
        "0.2-1.4 ng/mL, with dose-group medians near 0.4, 0.6 and 0.85 ng/mL",
        "at 4, 6 and 8 mg/kg. Derived in the source from empirical Bayes",
        "estimates of the Hong 2025 population PK model; reproduce with",
        "modellib('Hong_2025_datopotamab'), whose DXd output Cc_dxd is in",
        "ug/mL (multiply by 1000)."
      ),
      source_name = "DXd Cav"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 729L,
    n_studies = 3L,
    disease_state = "Advanced or metastatic non-small cell lung cancer (n = 644) or breast cancer (n = 85, TROPION-PanTumor01 only)",
    dose_range = "Datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; most patients at 6 mg/kg",
    regions = "USA, Europe, Japan, China and rest of world",
    notes = "Exposure-safety dataset: the 644-patient exposure-efficacy dataset (baseline characteristics in Table S3; median age 63 years, 46% female, 42% Asian) plus 85 breast-cancer patients from TROPION-PanTumor01; Hong 2026 prints no separate demographic table for the 729. Observed events by DXd Cav quartile 64/188, 75/180, 82/184 and 120/177 (Figure 4)."
  )

  ini({
    # Supplementary Equation S2 with a linear f(Exposure) for this endpoint
    # (Results). No coefficient is printed; nonlinear least squares of
    # expit(a + b * CAV) on the 327 pixel columns of the Figure 4 top-left
    # model-predicted median gives the values below, with a largest
    # residual of 0.004 on the probability scale.
    logite0 <- -1.092; label("Logit of the grade 3 or worse TEAE probability at zero DXd exposure (unitless logit)") # back-solved from the Figure 4 top-left drug-effect curve
    e_cav_logit <- 1.441; label("Change in the logit per unit of DXd Cav (1/(ng/mL))") # back-solved from the Figure 4 top-left drug-effect curve

    # The source likelihood is Bernoulli with no residual error; the
    # placeholder below only lets the nlmixr2 observation machinery accept
    # the model. Not from the source.
    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_teae_grade3 <- logite0 + e_cav_logit * CAV
    prob_teae_grade3 <- expit(logit_teae_grade3)
    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
