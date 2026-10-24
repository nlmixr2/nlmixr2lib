Thiel_2018_diclofenac_painrelief <- function() {
  description <- paste(
    "Exposure-response model linking the predicted decrease in prostaglandin",
    "E2 (PGE2) to placebo-adjusted mean pain relief after oral diclofenac in",
    "adults with acute pain after oral (dental) surgery (Thiel 2018).",
    "Direct sigmoid Hill model without baseline:",
    "score = emax * x^hill / (ec50^hill + x^hill), where x = -PGE2_DELTA is",
    "the PGE2 decrease (umol/L) relative to the drug-free scenario and score",
    "is the mean pain relief score minus the placebo mean at the same time.",
    "Fit by least squares to arm-mean data pooled over 25, 50 and 100 mg and",
    "the 0-6 h time course (33 points); residual SD is the fit RMSE. The",
    "driver is supplied as a time-varying covariate: in the source it is the",
    "output of the Yang 2007 arachidonic-acid network model driven by a",
    "PK-Sim whole-body PBPK prediction of diclofenac in venous blood cells,",
    "and neither of those layers is part of this model. Arm-mean simulator;",
    "no inter-individual variability was estimated.",
    sep = " "
  )
  reference <- paste(
    "Thiel C, Smit I, Baier V, Cordes H, Fabry B, Blank LM, Kuepfer L. (2018).",
    "Using quantitative systems pharmacology to evaluate the drug efficacy of",
    "COX-2 and 5-LOX inhibitors in therapeutic situations.",
    "npj Systems Biology and Applications 4:28.",
    "doi:10.1038/s41540-018-0062-3.",
    "Observed pain relief: Hersh EV et al. Dose-ranging analgesic study of",
    "Prosorb diclofenac potassium in postsurgical dental pain.",
    "Clin Ther 2004;26:1215-1227.",
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
        "competitive COX-2 inhibition, Ki(diclofenac, COX-2) = 0.018, driven",
        "by the PK-Sim 6.3 PBPK unbound diclofenac concentration in venous",
        "blood cells after a dissolved oral dose. Over 0-6 h it rose to",
        "about 41, 48 and 55 umol/L at 25, 50 and 100 mg (Figure 5a)."
      ),
      source_name = "PGE2 decrease (sign-reversed)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 197L,
    n_studies = 1L,
    disease_state = "Moderate-to-severe acute pain after oral (dental) surgery",
    dose_range = "Single oral dose of diclofenac potassium 25, 50 or 100 mg (Prosorb dissolved formulation in the PBPK simulation)",
    notes = paste(
      "Pain relief data from Hersh 2004 (Thiel 2018 ref 27): 25 mg n = 63,",
      "50 mg n = 68, 100 mg n = 66 (Thiel 2018 Methods). Arm-mean pain",
      "relief scores at 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5 and 6 h,",
      "minus the placebo arm mean (Figure 5a). The PGE2 driver was simulated",
      "for a PBPK reference individual (European man, 73 kg, 30 years), not",
      "for the trial patients. The paper reports Lissy 2010 NCA as its check",
      "of the dissolved formulation: 25 mg Cmax 788 vs observed 1125 ng/mL,",
      "AUC 724 vs 603 ng*h/mL; 50 mg Cmax 1609 vs 2035 ng/mL, AUC 1491 vs",
      "1232 ng*h/mL."
    )
  )

  ini({
    lemax <- log(2.8)
    label("Maximal placebo-adjusted mean pain relief score (score units)")
    # Thiel 2018 Results, Hill fit 'DFN: a = 2.8 CI (2.0, 3.5)'.
    lhill <- log(4.4)
    label("Hill coefficient of the PGE2 decrease effect (unitless)")
    # Thiel 2018 Results, 'DFN: ... b = 4.4 CI (2.6, 6.3)'.
    lec50 <- log(37.2)
    label("PGE2 decrease giving half the maximal pain relief (umol/L)")
    # Thiel 2018 Results, 'DFN: ... c = 37.2 CI (31.9, 42.5)'; units from Figure 5a x-axis.

    addSd <- 0.21
    label("Additive residual SD on the placebo-adjusted mean pain relief score (score units)")
    # Thiel 2018 Supplementary Table S4, row 'DFN (6h) Hill (3)' RMSE = 0.21 (SSE 1.30).
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
