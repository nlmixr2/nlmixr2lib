Hayato_2018_lenvatinib_tumor <- function() {
  description <- paste0(
    "Exposure-response tumor-growth-inhibition model for lenvatinib in ",
    "radioiodine-refractory differentiated thyroid cancer (Hayato 2018; ",
    "phase 3 SELECT study, 248 lenvatinib and 126 placebo patients). The ",
    "RECIST 1.1 sum of target-lesion diameters grows first-order (rate KG) ",
    "and is shrunk by an Emax function of lenvatinib average daily AUC ",
    "whose EC50 rises exponentially with time on treatment (acquired ",
    "resistance, rate lambda). Additive IIV on KG, Emax and lambda; ",
    "exponential IIV on EC50; combined additive and proportional residual ",
    "error. Lenvatinib exposure enters as the AUC_LEN covariate (no PK ",
    "ODE); the source derived it from the Gupta 2016 population PK model ",
    "(modellib('Gupta_2016_lenvatinib')). Time is in weeks. Companion ",
    "model: modellib('Hayato_2018_lenvatinib_teae_dosemod')."
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
    time = "week",
    dosing = "n/a (no drug-dosing events; lenvatinib exposure enters as the AUC_LEN covariate, not via a PK ODE)",
    concentration = "mm (RECIST 1.1 sum of diameters of target lesions; not a drug concentration)"
  )

  covariateData <- list(
    AUC_LEN = list(
      description = "Lenvatinib steady-state daily AUC computed from the average dose between two consecutive tumor assessments (AVAUC), driving the Emax tumor-shrinkage term.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying, step-wise between tumor assessments. Hayato 2018",
        "Methods define AVAUC as 'the exposure of lenvatinib based on",
        "average dose between two tumor assessments', derived from the",
        "individual PK parameters of the Gupta 2016 population PK model",
        "(i.e. F * average daily dose / (CL/F)). The paper reports AVAUC",
        "and the EC50 in ug*h/mL; this model carries the canonical ng*h/mL",
        "scale and converts the EC50 by x 1000 instead, so the Emax term is",
        "numerically unchanged. Set to 0 for placebo subjects and for",
        "off-treatment intervals, which removes the drug term exactly and",
        "leaves pure first-order growth."
      ),
      source_name = "AVAUC"
    ),
    TUM_SLD = list(
      description = "Observed baseline RECIST 1.1 sum of diameters of target lesions at study entry; per-subject, time-fixed. Initial condition of the tumor ODE.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Hayato 2018 Table 1 contains no baseline-size parameter and no",
        "baseline IIV, so the initial condition is taken as the observed",
        "baseline sum of diameters (the usual NONMEM data-item",
        "construction). Because the ODE is first-order in tumor size, the",
        "RELATIVE change from baseline does not depend on this value; it",
        "matters only through the additive residual error."
      ),
      source_name = "y(0) (baseline sum of diameters of target lesions)"
    )
  )

  compartmentData <- list(
    tumor = list(
      analyte = "tumor burden (RECIST 1.1 sum of diameters of target lesions)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 374L,
    n_studies = 1L,
    n_observations = "2373 independent-radiologic-review tumor assessments from 374 patients (248 lenvatinib, 126 placebo) with a baseline and at least one post-dose assessment.",
    disease_state = "Radioiodine-refractory differentiated thyroid cancer (RR-DTC) with measurable target lesions by RECIST 1.1.",
    dose_range = "Lenvatinib 24 mg orally once daily starting dose with protocol-driven reductions (20, 14, 10 mg), or placebo.",
    regions = "Multinational phase 3 SELECT study (E7080-G000-303, NCT01321554).",
    notes = paste(
      "Of 392 SELECT patients, 17 had no post-dose tumor assessment and 1",
      "had no PK data (Hayato 2018 Methods). Demographics of the analysis",
      "set are not reported in Hayato 2018; they are those of SELECT",
      "(Schlumberger 2015) and of the Gupta 2016 PK analysis. No covariate",
      "analysis was performed. Estimation by FOCEI in NONMEM 7.2.0 with",
      "PDxPop 5.0."
    )
  )

  ini({
    # Hayato 2018 Methods equation and Table 1:
    #   dy/dt = KG*y - Emax*AVAUC/(AVAUC + EC50*R(t))*y,  R(t) = exp(lambda*t)
    # Time base is weeks. The IIV on KG, Emax and lambda is ADDITIVE
    # ('described by an additive error model') and Table 1 reports those
    # three rows as SD, so these typical values are not log-transformed.
    kgrow_pop <- 0.000932; label("First-order tumor growth rate constant KG (1/week)") # Hayato 2018 Table 1: KG = 0.000932 /week (%RSE 63.0)
    emax_pop <- 0.0834; label("Maximum lenvatinib effect on the tumor shrinkage rate Emax (1/week)") # Hayato 2018 Table 1: Emax = 0.0834 /week (%RSE 16.2)
    lec50 <- log(2030); label("Lenvatinib average daily AUC giving half of Emax at t = 0, EC50 (ng*h/mL)") # Hayato 2018 Table 1: EC50 = 2.03 ug*h/mL (%RSE 32.1), x 1000 to ng*h/mL
    lam_pop <- 0.170; label("Resistance appearance rate constant lambda (1/week)") # Hayato 2018 Table 1: lambda = 0.170 /week (%RSE 8.47)

    # Inter-individual variability. Additive-eta rows are SDs, so the
    # variances below are SD^2. The EC50 row is reported as 'CV %' for an
    # exponential eta; it is taken as sqrt(omega^2) x 100 (the PDxPop and
    # same-sponsor Majid 2024 convention), so omega^2 = 0.405^2. The
    # log(1 + CV^2) reading would give 0.152; see the vignette.
    etakgrow ~ 0.00858^2 # Hayato 2018 Table 1: KG IIV SD = 0.00858 /week -> 0.00858^2
    etaemax ~ 0.0508^2 # Hayato 2018 Table 1: Emax IIV SD = 0.0508 /week -> 0.0508^2
    etalec50 ~ 0.405^2 # Hayato 2018 Table 1: EC50 IIV = 40.5 CV% -> 0.405^2
    etalam ~ 0.0791^2 # Hayato 2018 Table 1: lambda IIV SD = 0.0791 /week -> 0.0791^2

    # Residual error: 'combined additive and proportional'. Table 1 prints
    # the proportional part as 3.89 CV%, i.e. 0.0389 as a fraction.
    propSd <- 0.0389; label("Proportional residual error on tumor size (fraction)") # Hayato 2018 Table 1: Proportional = 3.89 CV% (%RSE 4.79)
    addSd <- 1.74; label("Additive residual error on tumor size (mm)") # Hayato 2018 Table 1: Additive SD = 1.74 mm (%RSE 3.18)
  })

  model({
    kgrow <- kgrow_pop + etakgrow
    emax <- emax_pop + etaemax
    ec50 <- exp(lec50 + etalec50)
    lam <- lam_pop + etalam

    # Resistance: the effective EC50 rises as exp(lambda * t) with t the
    # time on study in weeks (Hayato 2018 Methods, R(t)).
    drv <- emax * AUC_LEN / (AUC_LEN + ec50 * exp(lam * t))

    tumor(0) <- TUM_SLD
    d/dt(tumor) <- kgrow * tumor - drv * tumor

    # The source does not say whether the two residual components are
    # summed as standard deviations (combined1) or as variances (combined2);
    # combined1 follows the same sponsor's deposited control stream for the
    # successor Majid 2024 lenvatinib tumor model.
    tumor ~ add(addSd) + prop(propSd) + combined1()
  })
}
