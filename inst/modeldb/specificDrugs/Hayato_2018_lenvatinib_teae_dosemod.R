Hayato_2018_lenvatinib_teae_dosemod <- function() {
  description <- paste0(
    "Weekly cumulative-logit (proportional-odds) exposure-safety model for ",
    "dose-altering adverse events of lenvatinib in radioiodine-refractory ",
    "differentiated thyroid cancer (Hayato 2018; phase 3 SELECT study, 260 ",
    "lenvatinib and 131 placebo patients, 24 weeks of weekly reports). ",
    "Each week is graded 0 (no dose-altering adverse event), 1 (a ",
    "treatment-related adverse event leading to dose interruption) or 2 (a ",
    "treatment-related adverse event leading to dose reduction or drug ",
    "withdrawal). logit P(grade >= i) = B_i + Emax * AUC / (EC50 + AUC) + ",
    "eta, with AUC the lenvatinib steady-state daily AUC at the highest ",
    "dose taken during the week (AUC_LEN; 0 on placebo) and eta a ",
    "subject-level sensitivity shared by both thresholds. Outputs the ",
    "ordinal grade teae_dosemod and the two cumulative weekly ",
    "probabilities prob_teae_dosemod (grade >= 1) and prob_teae_dosered ",
    "(grade 2). Time is in weeks. Companion model: ",
    "modellib('Hayato_2018_lenvatinib_tumor')."
  )
  reference <- paste(
    "Hayato S, Shumaker R, Ferry J, Binder T, Dutcus CE, Hussein Z.",
    "Exposure-response analysis and simulation of lenvatinib safety and",
    "efficacy in patients with radioiodine-refractory differentiated thyroid",
    "cancer. Cancer Chemother Pharmacol. 2018;82(6):971-978.",
    "doi:10.1007/s00280-018-3687-4.",
    sep = " "
  )
  vignette <- "Hayato_2018_lenvatinib"
  units <- list(
    time = "week (one record per patient-week; the model is algebraic and has no time dependence of its own)",
    dosing = "n/a (no dose events; lenvatinib exposure enters through the AUC_LEN covariate)",
    concentration = "(ordinal weekly dose-altering adverse-event grade 0-2 in teae_dosemod; prob_teae_dosemod and prob_teae_dosered are weekly probabilities, 0-1)"
  )

  covariateData <- list(
    AUC_LEN = list(
      description = "Lenvatinib steady-state daily AUC at the HIGHEST dose taken during the week of the adverse-event record.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying, one value per patient-week. Hayato 2018 Methods:",
        "'relating this to lenvatinib AUC based on the highest dose during",
        "the week'; the individual AUC comes from the Gupta 2016 population",
        "PK model (F * dose / (CL/F)). Not frozen at the starting dose:",
        "unlike the landmark Tamai 2017 model, this is a weekly",
        "longitudinal model, and in the source simulations the dose (and",
        "so AUC_LEN) steps down after each grade-2 week. The paper reports",
        "AUC and EC50 in ug*h/mL; this model carries the canonical ng*h/mL",
        "scale and converts the EC50 by x 1000. Set to 0 for placebo",
        "patients and for weeks with no lenvatinib, which leaves the",
        "baseline logits B1 and B2."
      ),
      source_name = "AUC (based on the highest dose during the week)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 391L,
    n_studies = 1L,
    n_observations = "7914 weekly adverse-event records over the first 24 weeks from 391 patients (260 lenvatinib, 131 placebo).",
    disease_state = "Radioiodine-refractory differentiated thyroid cancer (RR-DTC).",
    dose_range = "Lenvatinib 24 mg orally once daily starting dose with protocol-driven reductions (20, 14, 10 mg), or placebo.",
    regions = "Multinational phase 3 SELECT study (E7080-G000-303, NCT01321554).",
    notes = paste(
      "Of 392 SELECT patients, 1 lenvatinib patient had no PK data (Hayato",
      "2018 Methods). Demographics are not reported in Hayato 2018; they",
      "are those of SELECT (Schlumberger 2015). A dose-altering adverse",
      "event is any treatment-related adverse event leading to dose",
      "interruption, dose reduction or drug withdrawal; all adverse-event",
      "types are pooled into one event with one exposure-response",
      "relationship. No covariate analysis was performed. Estimation by",
      "the Laplacian method in NONMEM 7.2.0 with PDxPop 5.0."
    )
  )

  ini({
    # Hayato 2018 Methods equations and Table 2:
    #   EFF = Emax * AUC / (EC50 + AUC) + eta
    #   A1 = B1 + EFF ; A2 = B2 + EFF ; P(AE >= i) = expit(A_i)
    # The source parameterises the logits of the EXCEEDANCE probabilities
    # P(grade >= i), so the threshold increment B2 - B1 is negative (the
    # opposite sign to the P(Y <= m) convention of the logitstep family).
    logite0 <- -5.64; label("Logit of the weekly probability of a dose-altering adverse event (grade >= 1) at zero exposure, B1 (unitless)") # Hayato 2018 Table 2: B1 = -5.64 (%RSE 5.14)
    logitstep1 <- -0.830; label("Increment from logit P(grade >= 1) to logit P(grade >= 2), B2 - B1 (unitless)") # Hayato 2018 Table 2: B2-B1 = -0.830 (%RSE 6.07)
    emax <- 6.55; label("Maximum lenvatinib effect on both cumulative logits, Emax (unitless logit)") # Hayato 2018 Table 2: Emax = 6.55 (%RSE 11.4)
    lec50 <- log(2950); label("Lenvatinib daily AUC giving half of Emax, EC50 (ng*h/mL)") # Hayato 2018 Table 2: EC50 = 2.95 ug*h/mL (%RSE 26.3), x 1000 to ng*h/mL

    # The subject-specific sensitivity eta enters EFF additively and so
    # shifts both logits equally, i.e. it is a shared eta on the linear
    # predictor. Table 2 reports its SD.
    etalogite0 ~ 1.06^2 # Hayato 2018 Table 2: IIV SD = 1.06 -> 1.06^2
  })

  model({
    ec50 <- exp(lec50)
    eff <- emax * AUC_LEN / (ec50 + AUC_LEN)

    logit_ge1 <- logite0 + etalogite0 + eff
    logit_ge2 <- logit_ge1 + logitstep1

    prob_teae_dosemod <- expit(logit_ge1)
    prob_teae_dosered <- expit(logit_ge2)

    # Per-grade probabilities (Hayato 2018 Methods):
    # P(AE = 0) = 1 - P(AE >= 1); P(AE = 1) = P(AE >= 1) - P(AE >= 2);
    # P(AE = 2) = P(AE >= 2).
    p0 <- 1 - prob_teae_dosemod
    p1 <- prob_teae_dosemod - prob_teae_dosered
    p2 <- prob_teae_dosered

    # rxode2 ordinal endpoint: simulated category k is grade k - 1
    # (category 3, probability 1 - p0 - p1 = p2, is grade 2).
    teae_dosemod ~ c(p0, p1)
  })
}
