Ruiz_2025_thrombin_generation <- function() {
  description <- "In vitro (calibrated automated thrombogram, platelet-poor plasma, 5 pM recombinant tissue factor). Semi-mechanistic model of thrombin generation in normal subjects and trauma patients: tissue factor activates prothrombin directly (second order, with thrombin positive feedback) and through a seven-compartment transit chain (second order); thrombin is degraded first order. All rate constants and their IIV are estimated separately for normal and trauma subjects."
  reference <- paste(
    "Ruiz L, Jaramillo S, Calvo A, Torrente MA, Tassies D, Reverter JC,",
    "Blasi A, Troconiz I. Dynamics of thrombin generation: Filling the gap",
    "between the system pharmacology theory and clinical practice in clinical",
    "pharmacology and therapeutics. Pharmacol Res Perspect. 2025;13(1):e70014.",
    "doi:10.1002/prp2.70014",
    sep = " "
  )
  vignette <- "Ruiz_2025_thrombin_generation"
  units <- list(time = "min", dosing = "nmol", concentration = "nM")

  covariateData <- list(
    DIS_TRAUMA = list(
      description = "Subject condition: 1 = trauma patient, 0 = normal (control) subject",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the group-specific typical value and IIV variance of every rate constant (Table 2 reports each parameter separately for normal and trauma subjects; Results 'the only covariate included in the model was the individual condition').",
      source_name = "condition (normal / trauma)"
    ),
    PROTHROMBIN_BASE = list(
      description = "Baseline plasma factor II (total prothrombin) concentration, measured per subject",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Initial condition of the prothrombin state (the model 'requires the values of FII and TF to predict FIIa', Discussion). Table 1: normal median 1108.5 (725.1-1338.06) nM; trauma median 1052.8 (474.09-1784.83) nM; all subjects 1066.71 nM.",
      source_name = "FII"
    )
  )

  compartmentData <- list(
    tissue_factor = list(analyte = "recombinant tissue factor", units = "nM", specimen = "plasma", verified = TRUE),
    effect1 = list(
      analyte = "tissue-factor-driven activation signal (transit 1)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect2 = list(
      analyte = "tissue-factor-driven activation signal (transit 2)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect3 = list(
      analyte = "tissue-factor-driven activation signal (transit 3)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect4 = list(
      analyte = "tissue-factor-driven activation signal (transit 4)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect5 = list(
      analyte = "tissue-factor-driven activation signal (transit 5)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect6 = list(
      analyte = "tissue-factor-driven activation signal (transit 6)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect7 = list(
      analyte = "tissue-factor-driven activation signal (transit 7)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    prothrombin = list(analyte = "prothrombin (factor II)", units = "nM", specimen = "plasma", verified = TRUE),
    thrombin = list(analyte = "thrombin (factor IIa)", units = "nM", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "in vitro (human platelet-poor plasma)",
    n_subjects = 60L,
    n_studies = 1L,
    disease_state = "20 normal (control) subjects and 40 trauma patients",
    n_observations = "1978 thrombin observations (normal) and 3736 (trauma); data after 25 min excluded",
    assay = "Calibrated automated thrombogram (CAT, Thrombinoscope) in platelet-poor plasma with 5 pM recombinant tissue factor and 4 uM procoagulant phospholipids",
    prothrombin_baseline = "Table 1 FII: normal median 1108.5 nM (725.1-1338.06); trauma median 1052.8 nM (474.09-1784.83)",
    notes = "Thrombin generation and coagulation-factor data from Menezes et al. (Ruiz 2025 Methods 'Research data'); baseline factors in Table 1. No age, weight or sex data are reported."
  )

  ini({
    # Assay input. Methods: 'a final concentration of 5 pM of recombinant tissue
    # factor', i.e. 0.005 nM in the model's nM concentration units.
    tf0 <- fixed(0.005); label("Initial tissue factor concentration (nM)")
    # Results / Discussion: kdegTF 'resulted negligible (p > .05)' and was 'fixed to 0'.
    kdeg_tf <- fixed(0); label("Tissue factor first-order degradation rate constant (1/min)")

    # Table 2 typical values, normal subjects (DIS_TRAUMA = 0)
    lksyn_thrombin_tf <- log(4.7e-3); label("Direct TF-driven prothrombin activation, normal (1/nM/min)") # Table 2 ksynIIa Normal = 4.7 x 10^-3
    lktr <- log(1.4); label("Transit rate constant, normal (1/min)") # Table 2 ktran Normal = 1.4
    lksyn_thrombin_effect <- log(47); label("Transit-chain-driven prothrombin activation, normal (1/nM/min)") # Table 2 ksynII2a Normal = 47
    lkdeg_thrombin <- log(1.04); label("Thrombin first-order degradation rate constant, normal (1/min)") # Table 2 kdegIIa Normal = 1.04

    # Table 2 typical values, trauma patients (DIS_TRAUMA = 1)
    lksyn_thrombin_tf_trauma <- log(0.15); label("Direct TF-driven prothrombin activation, trauma (1/nM/min)") # Table 2 ksynIIa Trauma = 0.15
    lktr_trauma <- log(1.8); label("Transit rate constant, trauma (1/min)") # Table 2 ktran Trauma = 1.8
    lksyn_thrombin_effect_trauma <- log(45); label("Transit-chain-driven prothrombin activation, trauma (1/nM/min)") # Table 2 ksynII2a Trauma = 45
    lkdeg_thrombin_trauma <- log(0.8); label("Thrombin first-order degradation rate constant, trauma (1/min)") # Table 2 kdegIIa Trauma = 0.8

    # IIV: Table 2 'IIV' column is CV% = sqrt(omega^2) * 100 (table footnote),
    # so omega^2 = (CV/100)^2. Separate variances per condition; no covariances.
    etalksyn_thrombin_tf_normal ~ 2.40^2 # Table 2 ksynIIa IIV normal 240%
    etalksyn_thrombin_tf_trauma ~ 0.51^2 # Table 2 ksynIIa IIV trauma 51%
    etalkdeg_thrombin_normal ~ 0.20^2 # Table 2 kdegIIa IIV normal 20%
    etalkdeg_thrombin_trauma ~ 0.25^2 # Table 2 kdegIIa IIV trauma 25%
    etalktr_normal ~ 0.23^2 # Table 2 ktran IIV normal 23%
    etalktr_trauma ~ 0.42^2 # Table 2 ktran IIV trauma 42%
    etalksyn_thrombin_effect_normal ~ 0.14^2 # Table 2 ksynII2a IIV normal 14%
    etalksyn_thrombin_effect_trauma ~ 0.34^2 # Table 2 ksynII2a IIV trauma 34%

    # Residual error: combined before 15 min, additive from 15 min (Results).
    propSd <- 0.18; label("Proportional residual error before 15 min (fraction)") # Table 2 'Proportinal error' 0.18
    addSd <- 2.4; label("Additive residual error before 15 min (nM)") # Table 2 Additive error Time < 15 min = 2.4
    addSd_late <- 2.4; label("Additive residual error from 15 min (nM)") # Table 2 Additive error Time > 15 min = 2.4
  })

  model({
    nrm <- 1 - DIS_TRAUMA

    # Individual rate constants: the condition selects both the typical value
    # and the IIV variance (Table 2 reports each separately per group).
    ksyn_thrombin_tf <- exp(
      nrm * (lksyn_thrombin_tf + etalksyn_thrombin_tf_normal) +
        DIS_TRAUMA * (lksyn_thrombin_tf_trauma + etalksyn_thrombin_tf_trauma)
    )
    ktr <- exp(
      nrm * (lktr + etalktr_normal) +
        DIS_TRAUMA * (lktr_trauma + etalktr_trauma)
    )
    ksyn_thrombin_effect <- exp(
      nrm * (lksyn_thrombin_effect + etalksyn_thrombin_effect_normal) +
        DIS_TRAUMA * (lksyn_thrombin_effect_trauma + etalksyn_thrombin_effect_trauma)
    )
    kdeg_thrombin <- exp(
      nrm * (lkdeg_thrombin + etalkdeg_thrombin_normal) +
        DIS_TRAUMA * (lkdeg_thrombin_trauma + etalkdeg_thrombin_trauma)
    )

    # Figure 2 equations. Tissue factor feeds the transit chain without being
    # consumed by it. The direct pathway carries the thrombin positive feedback
    # as printed, TF * FII * (1 + FIIa), with FIIa entering as its nM value.
    d/dt(tissue_factor) <- -kdeg_tf * tissue_factor
    d/dt(effect1) <- ktr * tissue_factor - ktr * effect1
    d/dt(effect2) <- ktr * effect1 - ktr * effect2
    d/dt(effect3) <- ktr * effect2 - ktr * effect3
    d/dt(effect4) <- ktr * effect3 - ktr * effect4
    d/dt(effect5) <- ktr * effect4 - ktr * effect5
    d/dt(effect6) <- ktr * effect5 - ktr * effect6
    d/dt(effect7) <- ktr * effect6 - ktr * effect7

    act_tf <- ksyn_thrombin_tf * tissue_factor * prothrombin * (1 + thrombin)
    act_effect <- ksyn_thrombin_effect * effect7 * prothrombin
    d/dt(prothrombin) <- -act_tf - act_effect
    d/dt(thrombin) <- act_tf + act_effect - kdeg_thrombin * thrombin

    tissue_factor(0) <- tf0
    prothrombin(0) <- PROTHROMBIN_BASE

    # Combined error up to 15 min, additive only afterwards (Results).
    early <- 1 * (time < 15)
    addSd_t <- early * addSd + (1 - early) * addSd_late
    propSd_t <- early * propSd
    thrombin ~ add(addSd_t) + prop(propSd_t)
  })
}
