Hong_2026_datopotamab_orr <- function() {
  description <- paste(
    "Logistic-regression exposure-efficacy model for objective response",
    "(best objective response of complete or partial response by blinded",
    "independent central review) in adults with previously treated advanced",
    "non-small cell lung cancer receiving datopotamab deruxtecan (Dato-DXd)",
    "0.27-10 mg/kg IV every 3 weeks (Hong 2026, n = 592). logit(P) = -7.845",
    "+ 8.00 * CAV / (4.04 + CAV), where CAV is the individual average",
    "plasma concentration of INTACT Dato-DXd (ug/mL) to the end of the cycle",
    "in which the response occurred. This is the paper's DRUG-EFFECT model",
    "(exposure only), whose curve is plotted in Figure 3: the EC50 is the",
    "printed 0.00404 mg/mL and the intercept and Emax are back-solved from",
    "the digitised Figure 3 curve. The final model also carries squamous",
    "histology (odds ratio 0.281) and more than two prior lines of therapy",
    "(odds ratio 1.82), but its intercept is not published, so those effects",
    "are not encoded (see the vignette). There is no PK layer: exposure is",
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
    "EC50 from the Results; intercept and Emax back-solved from Figure 3;",
    "model form from Supplementary Equation S2.",
    sep = " "
  )
  vignette <- "Hong_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr_central (probability of objective response by blinded independent central review, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average plasma concentration of INTACT datopotamab deruxtecan from the first dose to the end of the cycle in which the response occurred",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "INTACT conjugate, not the DXd payload: the companion",
        "Hong_2026_datopotamab_teae_grade3 and _sae models carry a DXd",
        "average concentration in ng/mL under the same column name. Hong",
        "2026 reports it in mg/mL (multiply by 1000 for ug/mL). Figure S3",
        "percentiles in the ORR dataset: 13.0 / 22.5 / 36.0 / 52.6 ug/mL",
        "(5th / 25th / 75th / 95th); medians 20.7, 28.9 and 40.7 ug/mL at",
        "4, 6 and 8 mg/kg. Derived in the source from empirical Bayes",
        "estimates of the Hong 2025 population PK model; reproduce with",
        "modellib('Hong_2025_datopotamab'), whose Dato-DXd output Cc is in",
        "ug/mL."
      ),
      source_name = "Dato-DXd Cav"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 592L,
    n_studies = 3L,
    age_range = "26-84 years (median 63; efficacy dataset of 644, Table S3)",
    weight_range = "37.0-156 kg (median 65.0; efficacy dataset of 644, Table S3)",
    sex_female_pct = 45.8,
    race_ethnicity = "44% White, 42% Asian, 1.9% Black, 11% other (efficacy dataset of 644, Table S3)",
    disease_state = "Previously treated advanced or metastatic non-small cell lung cancer with a response-evaluable tumour assessment",
    dose_range = "Datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; most patients at 6 mg/kg",
    regions = "USA, Europe, Japan, China and rest of world",
    notes = "Response-evaluable subset (n = 592) of the 644-patient exposure-efficacy dataset (TROPION-PanTumor01 NSCLC, TROPION-Lung05, TROPION-Lung01). Observed responders by Cav quartile 24/148, 46/148, 49/148 and 55/148 (Figure 3)."
  )

  ini({
    # Supplementary Equation S2: logit(P) = logit(Pbase) + f(Exposure),
    # with an Emax f for this endpoint (Results). The EC50 is printed; the
    # intercept and Emax are not. Least squares on the logit of the
    # Figure 3 model-predicted median read at Cav = 0.01, 0.02, ..., 0.08
    # mg/mL (P = 0.105, 0.232, 0.310, 0.360, 0.393, 0.414, 0.431, 0.443),
    # with EC50 held at the printed value, gives the values below; the
    # largest residual is 0.0075 logit units.
    logite0 <- -7.845; label("Logit of the objective-response probability at zero Dato-DXd exposure (unitless logit)") # back-solved from the Figure 3 drug-effect curve, EC50 held at the Results value
    emax <- 8.001; label("Maximal logit shift of the objective-response probability produced by Dato-DXd Cav (unitless logit)") # back-solved from the Figure 3 drug-effect curve, EC50 held at the Results value
    lec50 <- log(4.04); label("Dato-DXd Cav producing half of the maximal logit shift (ug/mL)") # Results 'The estimated EC50 was 0.00404 mg/mL' = 4.04 ug/mL

    # The source likelihood is Bernoulli with no residual error; the
    # placeholder below only lets the nlmixr2 observation machinery accept
    # the model. Not from the source.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50)
    logit_orr <- logite0 + emax * CAV / (ec50 + CAV)
    prob_orr_central <- expit(logit_orr)
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
