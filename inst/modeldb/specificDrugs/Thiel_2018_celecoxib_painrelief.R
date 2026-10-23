Thiel_2018_celecoxib_painrelief <- function() {
  description <- paste(
    "Exposure-response model linking the predicted decrease in prostaglandin",
    "E2 (PGE2) to placebo-adjusted mean pain relief after a single oral",
    "400 mg dose of celecoxib in adults with acute pain after oral (dental)",
    "surgery (Thiel 2018). Direct sigmoid Hill model without baseline:",
    "score = emax * x^hill / (ec50^hill + x^hill), where x = -PGE2_DELTA is",
    "the PGE2 decrease (umol/L) relative to the drug-free scenario and score",
    "is the mean pain relief score minus the placebo mean at the same time.",
    "Fit by least squares to the arm-mean 0-6 h time course (11 points);",
    "residual SD is the fit RMSE. The driver is supplied as a time-varying",
    "covariate: in the source it is the output of the Yang 2007",
    "arachidonic-acid network model driven by a PK-Sim whole-body PBPK",
    "prediction of celecoxib in venous blood cells, and neither of those",
    "layers is part of this model. Arm-mean simulator; no inter-individual",
    "variability was estimated.",
    sep = " "
  )
  reference <- paste(
    "Thiel C, Smit I, Baier V, Cordes H, Fabry B, Blank LM, Kuepfer L. (2018).",
    "Using quantitative systems pharmacology to evaluate the drug efficacy of",
    "COX-2 and 5-LOX inhibitors in therapeutic situations.",
    "npj Systems Biology and Applications 4:28.",
    "doi:10.1038/s41540-018-0062-3.",
    "Observed pain relief: Cheung R, Krishnaswami S, Kowalski K. Analgesic",
    "efficacy of celecoxib in postoperative oral surgery pain.",
    "Clin Ther 2007;29:2498-2510.",
    sep = " "
  )
  vignette <- "Thiel_2018_cox2_pge2_painrelief"
  units <- list(
    time = "h",
    dosing = "(not applicable; the PGE2 change is supplied as the time-varying covariate PGE2_DELTA)",
    concentration = "umol/L (PGE2 change via PGE2_DELTA; the output score is a pain relief score difference)"
  )

  covariateData <- list(
    PGE2_DELTA = list(
      description = "Drug-treated minus drug-free PGE2 concentration (negative = decrease); the driver of the Hill pain-relief model",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. model() derives the paper's positive-for-suppression",
        "axis as pge2_decr <- -PGE2_DELTA (Thiel 2018 Figure 5 x-axis",
        "'PGE2 decrease [umol/l]'; Methods 'difference of PGE2 concentrations",
        "between the drug-free treatment and the drug-treated case'). In the",
        "source the value is model-predicted, not measured: the Yang 2007",
        "arachidonic-acid network (human polymorphonuclear leukocytes) with",
        "competitive COX-2 inhibition, Ki(celecoxib, COX-2) = 0.017, driven",
        "by the PK-Sim 6.3 PBPK unbound celecoxib concentration in venous",
        "blood cells after a dissolved oral dose. After 400 mg it rose",
        "steadily to about 103 umol/L at 6 h (Figure 5b)."
      ),
      source_name = "PGE2 decrease (sign-reversed)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 57L,
    n_studies = 1L,
    disease_state = "Moderate-to-severe acute pain after oral (dental) surgery",
    dose_range = "Single oral dose of celecoxib 400 mg",
    notes = paste(
      "Pain relief data from Cheung 2007 (Thiel 2018 ref 26): n = 57 on",
      "celecoxib (Thiel 2018 Methods). Arm-mean pain relief scores at 0,",
      "0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5 and 6 h, minus the placebo arm",
      "mean (Figure 5b). The PGE2 driver was simulated for a PBPK reference",
      "individual (European man, 73 kg, 30 years), not for the trial",
      "patients. A 24 h celecoxib fit is evaluated in Supplementary Table S4",
      "(RMSE 0.22) but its parameters are not reported, so it is not",
      "included."
    )
  )

  ini({
    lemax <- log(2)
    label("Maximal placebo-adjusted mean pain relief score (score units)")
    # Thiel 2018 Results, Hill fit 'CEL: a = 2 CI (1.9, 2.1)'.
    lhill <- log(5.7)
    label("Hill coefficient of the PGE2 decrease effect (unitless)")
    # Thiel 2018 Results, 'CEL: ... b = 5.7 CI (3.7, 7.8)'.
    lec50 <- log(43.4)
    label("PGE2 decrease giving half the maximal pain relief (umol/L)")
    # Thiel 2018 Results, 'CEL: ... c = 43.4 CI (40.6, 46.1)'; units from Figure 5b x-axis.

    addSd <- 0.10
    label("Additive residual SD on the placebo-adjusted mean pain relief score (score units)")
    # Thiel 2018 Supplementary Table S4, row 'CEL (6h) Hill (3)' RMSE = 0.10 (SSE 0.08).
  })

  model({
    pge2_decr <- -PGE2_DELTA

    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    score <- emax * pge2_decr^hill / (ec50^hill + pge2_decr^hill)
    score ~ add(addSd)
  })
}
