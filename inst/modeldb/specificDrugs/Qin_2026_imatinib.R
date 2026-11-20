Qin_2026_imatinib <- function() {
  description <- paste(
    "PBPK (GastroPlus 9.8.2) pH-dependent acid-reducing-agent drug",
    "interaction paper, compartmental disposition layer.",
    "Three-compartment intravenous imatinib model with body-weight-proportional",
    "central volume and clearance, fitted by the authors to published mean",
    "plasma concentrations after a 100 mg intravenous infusion over 1 h",
    "(Peng et al. 2004). The paper's oral ACAT absorption, its gastric-pH",
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
    "Source intravenous data: Peng B, Dutreix C, Mehring G, et al. J Clin Pharmacol. 2004;44(2):158-162.",
    sep = " "
  )
  vignette <- "Qin_2026_acid_reducing_agent_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("central")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imatinib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "imatinib", units = "mg", specimen = "tissue", verified = TRUE)
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
    n_subjects = 17L,
    n_studies = 1L,
    age_range = "Adults aged 40 to 58 years (Peng 2004, Supplementary Table S2)",
    weight_range = "Mean 73.3 kg (SD 7) (Peng 2004, Supplementary Table S2)",
    disease_state = "Healthy volunteers (literature mean intravenous profiles)",
    dose_range = "100 mg IV over 1 h (fitted); 100 to 400 mg oral tablets and capsules (paper, not included)",
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
    lvc <- fixed(log(0.247)); label("Central volume per kg body weight (L/kg)") # Supplementary Table S1, Imatinib, 'Vc = 0.247 L/kg', 'Estimated by fitting to IV' (Peng et al. 2004)
    lcl <- fixed(log(0.189)); label("Systemic clearance per kg body weight (L/h/kg)") # Supplementary Table S1, Imatinib, 'CL = 0.189 L/h/kg'
    lk12 <- fixed(log(4.471)); label("Central-to-peripheral1 rate constant K12 (1/h)") # Supplementary Table S1, Imatinib, 'k12 = 4.471 1/h'
    lk21 <- fixed(log(1.226)); label("Peripheral1-to-central rate constant K21 (1/h)") # Supplementary Table S1, Imatinib, 'k21 = 1.226 1/h'
    lk13 <- fixed(log(1.706)); label("Central-to-peripheral2 rate constant K13 (1/h)") # Supplementary Table S1, Imatinib, 'k13 = 1.706 1/h'
    lk31 <- fixed(log(0.219)); label("Peripheral2-to-central rate constant K31 (1/h)") # Supplementary Table S1, Imatinib, 'k31 = 0.219 1/h'
  })

  model({
    vc <- exp(lvc) * WT
    cl <- exp(lcl) * WT
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Plasma imatinib concentration, mg/L (= ug/mL).
    Cc <- central / vc
  })
}
