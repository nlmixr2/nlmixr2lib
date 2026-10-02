Suri_2019_brentuximab_fn_mmae <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for febrile neutropenia ",
    "driven by monomethyl auristatin E (MMAE) payload exposure in adults ",
    "with previously untreated stage III or IV classical Hodgkin lymphoma ",
    "treated with brentuximab vedotin 1.2 mg/kg every 2 weeks plus AVD in ",
    "ECHELON-1 (Suri 2019). The probability of an event is ",
    "expit(b0 + b1 * AUC_BV_MMAE + bG * CONMED_GCSF), where AUC_BV_MMAE is ",
    "the individual time-averaged MMAE AUC up to the first febrile ",
    "neutropenia event (ng*day/mL) and CONMED_GCSF flags primary ",
    "prophylactic G-CSF. No PK layer: exposure is supplied as data from ",
    "the companion popPK model Suri_2019_brentuximab. The paper fits a ",
    "separate ADC-driven model for the same endpoint ",
    "(Suri_2019_brentuximab_fn_adc). NOTE: the coefficients are not ",
    "tabulated; they were recovered by digitising both fitted curves of ",
    "Figure 4d (see the vignette)."
  )
  reference <- paste(
    "Suri A, Mould DR, Song G, Collins GP, Endres CJ, Gomez-Navarro J,",
    "Venkatakrishnan K. Population Pharmacokinetic Modeling and",
    "Exposure-Response Assessment for the Antibody-Drug Conjugate",
    "Brentuximab Vedotin in Hodgkin's Lymphoma in the Phase III ECHELON-1",
    "Study. Clin Pharmacol Ther. 2019;106(6):1268-1279. doi:10.1002/cpt.1530.",
    "PMCID PMC6896233. The regression coefficients are not tabulated by the",
    "source; they were digitised from the two fitted lines in Figure 4d.",
    sep = " "
  )
  vignette <- "Suri_2019_brentuximab"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_BV_MMAE covariate column)",
    concentration = "prob_febrile_neutropenia (probability of a febrile neutropenia event, 0-1)"
  )

  covariateData <- list(
    AUC_BV_MMAE = list(
      description = "Individual time-averaged AUC of monomethyl auristatin E released from brentuximab vedotin, up to the event (or over the treatment duration if no event occurred). Supplied as data: this model has no PK layer.",
      units = "ng*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Suri 2019 Methods: AUC averaged over the time to first occurrence of the event. Figure 4d axis spans about 2-41 ng*day/mL. Enters LINEARLY and uncentred. Reproduce with modellib('Suri_2019_brentuximab').",
      source_name = "MMAE AUC/time"
    ),
    CONMED_GCSF = list(
      description = "Primary prophylactic granulocyte colony-stimulating factor (1 = received primary prophylaxis, 0 = did not)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Fixed (time-invariant) covariate in Suri 2019. Additive shift on the logit: odds ratio exp(-0.8917) = 0.410.",
      source_name = "Prophylactic G-CSF"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 661L,
    n_studies = 1L,
    n_observations = "one binary event record per patient (landmark analysis, no repeated measures)",
    disease_state = "Previously untreated stage III or IV classical Hodgkin lymphoma (ECHELON-1, NCT01712490, A+AVD arm).",
    dose_range = "Brentuximab vedotin 1.2 mg/kg every 2 weeks (capped at 120 mg above 100 kg) with doxorubicin, vinblastine and dacarbazine, up to 6 cycles.",
    notes = "MMAE AUC/time was a significant predictor of febrile neutropenia (P < 0.001). Observed febrile neutropenia incidence in the A+AVD arm was 19% overall, 21% without and 11% with primary G-CSF prophylaxis."
  )

  ini({
    # ==================================================================
    # Provenance: DIGITISED FROM FIGURE 4d (red = no G-CSF, blue = G-CSF),
    # traced at 400 dpi over 918 pixel columns per curve (MMAE AUC 2-41
    # ng*day/mL). Separate fits give slopes 0.1093 and 0.1097 per
    # ng*day/mL (shared-slope model confirmed); joint-fit root-mean-square
    # residual 0.21 percentage points. Implied G-CSF odds reduction 59%,
    # inside the paper's stated '55-81%'.
    # ==================================================================
    b0_fn        <- -2.784; label("Logit of febrile neutropenia probability at zero MMAE AUC/time without G-CSF (logit)")  # digitised from Suri 2019 Figure 4d; no printed value exists
    e_auc_fn     <- 0.1095; label("Log-odds of febrile neutropenia per ng*day/mL of MMAE AUC/time (logit per ng*day/mL)") # digitised from Suri 2019 Figure 4d; no printed value exists
    e_conmed_gcsf_fn <- -0.8917; label("Log-odds shift of febrile neutropenia with primary G-CSF prophylaxis (logit)")   # digitised from Suri 2019 Figure 4d; no printed value exists

    addSd_prob_febrile_neutropenia <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_febrile_neutropenia <- b0_fn + e_auc_fn * AUC_BV_MMAE + e_conmed_gcsf_fn * CONMED_GCSF
    prob_febrile_neutropenia  <- expit(logit_febrile_neutropenia)

    prob_febrile_neutropenia ~ add(addSd_prob_febrile_neutropenia)
  })
}
