Qin_2026_erlotinib <- function() {
  description <- paste(
    "PBPK (GastroPlus 9.8.2) pH-dependent acid-reducing-agent drug",
    "interaction paper, compartmental disposition layer.",
    "Two-compartment intravenous erlotinib model with body-weight-proportional",
    "central volume and clearance, fitted by the authors to published mean",
    "plasma concentrations after a 25 mg intravenous infusion over 0.5 h",
    "(Frohna et al. 2006). The paper's oral ACAT absorption, its gastric-pH",
    "interaction simulations and its XGBoost machine-learning model are not",
    "included (see the vignette). Deterministic typical-value model with no",
    "between-subject or residual variability.",
    sep = " "
  )
  reference <- paste(
    "Qin Y, Yu L, Chen T, Zuo Y, Wang H, Yang G, Pei Q. Machine Learning",
    "Enables Rapid Prediction of Acid-Reducing Agent Drug Interactions: A",
    "Streamlined Complement to PBPK Modeling. CPT Pharmacometrics Syst",
    "Pharmacol. 2026;15(8):e70314. doi:10.1002/psp4.70314.",
    "Disposition parameters are from Supplementary Table S1 (row",
    "'Disposition model parameters'); the intravenous checks are from",
    "Supplementary Tables S2 and S3.",
    "Source intravenous data: Frohna P, Lu J, Eppler S, et al. J Clin Pharmacol. 2006;46(3):282-290.",
    sep = " "
  )
  vignette <- "Qin_2026_acid_reducing_agent_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("central")

  compartmentData <- list(
    central = list(analyte = "erlotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "erlotinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "GastroPlus compartmental parameters are entered per kg and scale",
        "linearly with the simulated subject's weight (Supplementary Table S1",
        "prints Vc in L/kg and CL in L/h/kg)."
      ),
      source_name = "Weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "Healthy adult volunteers; mean age 33.7 years (SD 13.16) (Frohna 2006, Supplementary Table S2)",
    weight_range = "Not reported for the intravenous study (Supplementary Table S2)",
    disease_state = "Healthy volunteers (literature mean intravenous profiles)",
    dose_range = "25 mg IV over 0.5 h (fitted); 150 mg oral tablets (paper, not included)",
    regions = "Not reported",
    notes = paste(
      "No individual data were fitted. The authors fitted the compartmental",
      "disposition model to a published mean intravenous profile and used it",
      "as the systemic layer of a GastroPlus ACAT model of the oral drug,",
      "which was verified against oral profiles and then used to simulate the",
      "effect of raised gastric pH on oral exposure (Table 1). The simulated",
      "intravenous AUC and Cmax are compared with observed values in",
      "Supplementary Table S3. No other intravenous study was used.",
      sep = " "
    )
  )

  ini({
    lvc <- fixed(log(0.695)); label("Central volume per kg body weight (L/kg)") # Supplementary Table S1, Erlotinib, 'Vc = 0.695 L/kg', 'Estimated by fitting to IV' (Frohna et al. 2006)
    lcl <- fixed(log(0.059)); label("Systemic clearance per kg body weight (L/h/kg)") # Supplementary Table S1, Erlotinib, 'CL = 0.059 L/h/kg'
    lk12 <- fixed(log(0.215)); label("Central-to-peripheral1 rate constant K12 (1/h)") # Supplementary Table S1, Erlotinib, 'k12 = 0.215 1/h'
    lk21 <- fixed(log(0.141)); label("Peripheral1-to-central rate constant K21 (1/h)") # Supplementary Table S1, Erlotinib, 'k21 = 0.141 1/h'
  })

  model({
    vc <- exp(lvc) * WT
    cl <- exp(lcl) * WT
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma erlotinib concentration, mg/L (= ug/mL).
    Cc <- central / vc
  })
}
