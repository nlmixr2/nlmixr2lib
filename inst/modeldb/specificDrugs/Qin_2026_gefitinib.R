Qin_2026_gefitinib <- function() {
  description <- paste(
    "PBPK (GastroPlus 9.8.2) pH-dependent acid-reducing-agent drug",
    "interaction paper, compartmental disposition layer.",
    "Three-compartment intravenous gefitinib model with body-weight-proportional",
    "central volume and clearance, fitted by the authors to published mean",
    "plasma concentrations after a 50 and 100 mg intravenous infusions over 1 h",
    "(Swaisland et al. 2005). The paper's oral ACAT absorption, its gastric-pH",
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
    "Source intravenous data: Swaisland HC, Smith RP, Laight A, et al. Clin Pharmacokinet. 2005;44(11):1165-1177.",
    sep = " "
  )
  vignette <- "Qin_2026_acid_reducing_agent_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("central")

  compartmentData <- list(
    central = list(analyte = "gefitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gefitinib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "gefitinib", units = "mg", specimen = "tissue", verified = TRUE)
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
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "Healthy adult volunteers aged 18 to 65 years (Swaisland 2005, Supplementary Table S2)",
    weight_range = "Not reported for the intravenous studies (Supplementary Table S2)",
    disease_state = "Healthy volunteers (literature mean intravenous profiles)",
    dose_range = "50 and 100 mg IV over 1 h (fitted); 250 mg oral tablets (paper, not included)",
    regions = "Not reported",
    notes = paste(
      "No individual data were fitted. The authors fitted the compartmental",
      "disposition model to a published mean intravenous profile and used it",
      "as the systemic layer of a GastroPlus ACAT model of the oral drug,",
      "which was verified against oral profiles and then used to simulate the",
      "effect of raised gastric pH on oral exposure (Table 1). The simulated",
      "intravenous AUC and Cmax are compared with observed values in",
      "Supplementary Table S3. Both intravenous doses are also the intravenous checks in Supplementary Table S3.",
      sep = " "
    )
  )

  ini({
    lvc <- fixed(log(1.046)); label("Central volume per kg body weight (L/kg)") # Supplementary Table S1, Gefitinib, 'Vc = 1.046 L/kg', 'Estimated by fitting to IV' (Swaisland et al. 2005)
    lcl <- fixed(log(0.792)); label("Systemic clearance per kg body weight (L/h/kg)") # Supplementary Table S1, Gefitinib, 'CL = 0.792 L/h/kg'
    lk12 <- fixed(log(6.237)); label("Central-to-peripheral1 rate constant K12 (1/h)") # Supplementary Table S1, Gefitinib, 'k12 = 6.237 1/h'
    lk21 <- fixed(log(0.547)); label("Peripheral1-to-central rate constant K21 (1/h)") # Supplementary Table S1, Gefitinib, 'k21 = 0.547 1/h'
    lk13 <- fixed(log(0.237)); label("Central-to-peripheral2 rate constant K13 (1/h)") # Supplementary Table S1, Gefitinib, 'k13 = 0.237 1/h'
    lk31 <- fixed(log(0.031)); label("Peripheral2-to-central rate constant K31 (1/h)") # Supplementary Table S1, Gefitinib, 'k31 = 0.031 1/h'
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

    # Plasma gefitinib concentration, mg/L (= ug/mL).
    Cc <- central / vc
  })
}
