Suri_2019_brentuximab_pn <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for grade 2 or worse ",
    "peripheral neuropathy in adults with previously untreated stage III ",
    "or IV classical Hodgkin lymphoma treated with brentuximab vedotin ",
    "1.2 mg/kg every 2 weeks plus AVD in ECHELON-1 (Suri 2019). The ",
    "probability of an event is expit(b0 + b1 * AUC_BV_ADC), where ",
    "AUC_BV_ADC is the individual time-averaged antibody-drug conjugate ",
    "(ADC) AUC up to the worst grade 2 or worse peripheral neuropathy ",
    "event (or over the treatment duration if none occurred), in ",
    "ug*day/mL. No covariate was retained. There is no PK layer: exposure ",
    "is supplied as data and was derived in the source from post hoc ",
    "estimates of the companion popPK model Suri_2019_brentuximab. NOTE: ",
    "Suri 2019 does not tabulate the regression coefficients; they were ",
    "recovered by digitising the fitted curve of Figure 4b (see the ",
    "vignette)."
  )
  reference <- paste(
    "Suri A, Mould DR, Song G, Collins GP, Endres CJ, Gomez-Navarro J,",
    "Venkatakrishnan K. Population Pharmacokinetic Modeling and",
    "Exposure-Response Assessment for the Antibody-Drug Conjugate",
    "Brentuximab Vedotin in Hodgkin's Lymphoma in the Phase III ECHELON-1",
    "Study. Clin Pharmacol Ther. 2019;106(6):1268-1279. doi:10.1002/cpt.1530.",
    "PMCID PMC6896233. The regression coefficients are not tabulated by the",
    "source; they were digitised from the fitted line in Figure 4b.",
    sep = " "
  )
  vignette <- "Suri_2019_brentuximab"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_BV_ADC covariate column)",
    concentration = "prob_peripheral_neuropathy_grade2 (probability of a grade 2 or worse peripheral neuropathy event, 0-1)"
  )

  covariateData <- list(
    AUC_BV_ADC = list(
      description = "Individual time-averaged brentuximab vedotin antibody-drug conjugate AUC up to the event (or over the treatment duration if no event occurred). Supplied as data: this model has no PK layer.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Suri 2019 Methods: individual AUCs from the popPK post hoc estimates were averaged over the time to the first occurrence of the worst-grade event, giving an exposure intensity that accounts for dose reductions and delays. Values span about 20-110 ug*day/mL (Figure 4b axis). Enters LINEARLY and uncentred. Reproduce with modellib('Suri_2019_brentuximab').",
      source_name = "ADC AUC/time"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 661L,
    n_studies = 1L,
    n_observations = "one binary event record per patient (landmark analysis, no repeated measures)",
    disease_state = "Previously untreated stage III or IV classical Hodgkin lymphoma (ECHELON-1, NCT01712490, A+AVD arm).",
    dose_range = "Brentuximab vedotin 1.2 mg/kg every 2 weeks (capped at 120 mg above 100 kg) with doxorubicin, vinblastine and dacarbazine, up to 6 cycles.",
    notes = "Grade 2 or worse peripheral neuropathy included peripheral sensory and motor neuropathy, paresthesia, hypoesthesia, polyneuropathy, muscular weakness and demyelinating polyneuropathy (CTCAE). ADC AUC/time was a significant predictor (P = 0.004); MMAE AUC/time was not (Figure S1c), so no MMAE-driven PN model exists. Events occurring from the first dose until 30 days after the last dose of frontline therapy."
  )

  ini({
    # ==================================================================
    # Provenance: DIGITISED FROM FIGURE 4b, not read from text or a table;
    # Suri 2019 tabulates no exposure-safety regression coefficients.
    # Method: the black fitted line was traced at 400 dpi over 963 pixel
    # columns (ADC AUC 21-111 ug*day/mL). The linear-exposure logit that
    # the supplement specifies fits with a root-mean-square residual of
    # 0.26 percentage points (0.02351 per ug*day/mL); a log-exposure logit
    # fits worse (1.78 points).
    # ==================================================================
    b0_pn <- -3.502;  label("Logit of grade 2 or worse peripheral neuropathy probability at zero ADC AUC/time (logit)") # digitised from Suri 2019 Figure 4b; no printed value exists
    e_auc_pn <- 0.02351; label("Log-odds of grade 2 or worse peripheral neuropathy per ug*day/mL of ADC AUC/time (logit per ug*day/mL)") # digitised from Suri 2019 Figure 4b; no printed value exists

    # Binomial logistic regression (Bernoulli likelihood); the placeholder
    # additive residual exists only so rxode2 accepts an endpoint.
    addSd_prob_peripheral_neuropathy_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_peripheral_neuropathy_grade2 <- b0_pn + e_auc_pn * AUC_BV_ADC
    prob_peripheral_neuropathy_grade2  <- expit(logit_peripheral_neuropathy_grade2)

    prob_peripheral_neuropathy_grade2 ~ add(addSd_prob_peripheral_neuropathy_grade2)
  })
}
