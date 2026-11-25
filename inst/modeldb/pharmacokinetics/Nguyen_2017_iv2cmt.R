Nguyen_2017_iv2cmt <- function() {
  description <- "Methodology reference. Two-compartment IV-bolus generative model for a hypothetical drug from the ISoP Model Evaluation Group tutorial on goodness-of-fit graphics (Nguyen 2017, pharmacokinetic case example). There is no molecule and no patients: the authors chose the parameter values and simulated 180 virtual subjects (10, 100 or 1000 mg) to illustrate how evaluation graphs respond to correct and deliberately misspecified models. Clearance carries a binary concomitant-treatment ('drug X') effect and central volume a log-weight effect; IIV is log-normal with a common SD of 0.5 and a 0.5 correlation between CL and V1; the residual error is combined additive plus proportional. Typical values are the as-run simulation values, which are one tenth of those printed in Supplementary Table S1 (see the vignette)."
  reference <- paste(
    "Nguyen THT, Mouksassi M-S, Holford N, Al-Huniti N, Freedman I, Hooker AC, John J,",
    "Karlsson MO, Mould DR, Perez Ruixo JJ, Plan EL, Savic R, van Hasselt JGC, Weber B,",
    "Zhou C, Comets E, Mentre F; Model Evaluation Group of the International Society of",
    "Pharmacometrics (ISoP) Best Practice Committee. Model Evaluation of Continuous Data",
    "Pharmacometric Models: Metrics and Graphics. CPT Pharmacometrics Syst Pharmacol.",
    "2017;6(2):87-109. doi:10.1002/psp4.12161. PMID 27884052. PMC5321813."
  )
  vignette <- "Nguyen_2017_model_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "hypothetical drug", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hypothetical drug", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters central volume through T_Weight = log(WT/70) (Nguyen 2017 Eq. 6), i.e. V1 = V1_typ * (WT/70)^0.75. Simulated as Normal(70, 15) kg (PK case example, Structural model).",
      source_name = "Weight (T_Weight = log(Weight/70) in the deposited dataset)"
    ),
    CONMED_DRUGX = list(
      description = "Concomitant treatment with the hypothetical interacting 'drug X' (1 = with, 0 = without)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant treatment)",
      notes = "CL = CL_typ * exp(0.5 * CONMED_DRUGX), a 1.65-fold increase with the concomitant treatment (Eq. 5; Supplementary Table S1 'beta concomitant treatment (Yes)' = 0.5). Half of the simulated subjects received it. Figure 1d labels the groups 'co-administration with drug X'.",
      source_name = "Sex (column name in the deposited PK dataset and individual-parameter files)"
    )
  )

  population <- list(
    species = "None (methodology paper; simulation-only hypothetical drug with no patients and no fitted estimates).",
    n_subjects = 180L,
    n_studies = 1L,
    weight_range = "Simulated Normal(70, 15) kg",
    disease_state = "N/A (Monte Carlo simulation for a model-evaluation tutorial; not a fit of any real molecule).",
    dose_range = "Single IV bolus of 10, 100 or 1000 mg (60 subjects per dose); samples at 0.5, 1, 4, 12 and 24 h (standard design), with 2- and 1-sample sparse designs drawn from those times.",
    regions = "N/A",
    scope_note = paste(
      "Filed under inst/modeldb/pharmacokinetics/ rather than specificDrugs/ because the paper's",
      "pharmacokinetic case example is a hypothetical drug. The values are the authors'",
      "data-generating (true) model; the misspecified models fitted to the simulated data have no",
      "printed estimates and are not packaged."
    ),
    notes = paste(
      "Supplementary Table S1 prints CL = 1 L/h, V1 = 5 L, Q = 1 L/h, V2 = 10 L. The deposited",
      "simulation (Supplementary R codes and data) was generated with values one tenth of these:",
      "the geometric means of the 180 deposited empirical-Bayes individual parameters, with the",
      "covariate effects removed, are 0.102 L/h, 0.528 L, 0.113 L/h and 0.998 L; the deposited script",
      "carries the true-model fit as 0.0977 / 0.518 / 0.113 / 1.01; and the published Figure 1",
      "shows observations up to about 4500 mg/L, which needs V1 near 0.5 L at 1000 mg. The packaged",
      "model uses the as-run values. Because all four disposition parameters scale together, the",
      "Table S1 values give the same curve shape at one tenth of the concentration."
    )
  )

  ini({
    lcl <- fixed(log(0.1)); label("Clearance at CONMED_DRUGX = 0 (L/h)") # Supp Table S1 'muCL' printed 1 L/h; as-run value 0.1 L/h (deposited simulation; see population$notes)
    lvc <- fixed(log(0.5)); label("Central volume V1 at 70 kg (L)") # Supp Table S1 'muV1' printed 5 L; as-run value 0.5 L
    lq <- fixed(log(0.1)); label("Intercompartmental clearance Q (L/h)") # Supp Table S1 'muQ' printed 1 L/h; as-run value 0.1 L/h
    lvp <- fixed(log(1)); label("Peripheral volume V2 (L)") # Supp Table S1 'muV2' printed 10 L; as-run value 1 L

    e_conmed_drugx_cl <- fixed(0.5); label("Log-scale effect of concomitant drug X on CL (unitless)") # Supp Table S1 'beta concomitant treatment (Yes)' = 0.5
    e_wt_vc <- fixed(0.75); label("Coefficient of log(WT/70) on V1 (unitless)") # Supp Table S1 'beta T_weight' = 0.75

    # Supp Table S1 omega row (symbol lost in the docx; Methods: 'all diagonal elements equal to
    # omega^2, where omega is the SD') = 0.5, so each variance is 0.25; corr(etaCL, etaV1) = 0.5
    # gives covariance 0.5 * 0.5 * 0.5 = 0.125.
    etalcl + etalvc ~ fixed(c(0.25, 0.125, 0.25))
    etalq ~ fixed(0.25) # Supp Table S1 omega = 0.5 (SD), variance 0.25
    etalvp ~ fixed(0.25) # Supp Table S1 omega = 0.5 (SD), variance 0.25

    addSd <- fixed(0.3); label("Additive residual SD (mg/L)") # Supp Table S1 'sigma inter' = 0.3
    propSd <- fixed(0.15); label("Proportional residual SD (fraction)") # Supp Table S1 'sigma slope' = 0.15
  })

  model({
    # Eq. 5-6: theta = mu * exp(beta * z) * exp(eta), with T_Weight = log(WT/70)
    cl <- exp(lcl + e_conmed_drugx_cl * CONMED_DRUGX + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / 70) + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Eq. 1
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    # Eq. 3: residual SD = sigma_add + sigma_prop * f
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
