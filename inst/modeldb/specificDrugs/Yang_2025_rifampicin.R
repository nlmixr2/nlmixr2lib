Yang_2025_rifampicin <- function() {
  description <- paste(
    "PBPK-DDI perpetrator model (GastroPlus 9.7), compartmental disposition",
    "layer. Two-compartment intravenous rifampicin model in a 70 kg adult with",
    "body-weight-proportional volume and clearance, as deposited in the",
    "paper's GastroPlus compound record 'Rifampicin-ComPK'. Systemic",
    "clearance is the sum of three linear arms: CYP3A4-mediated hepatic",
    "metabolism, the remaining hepatic metabolism and renal clearance. The",
    "paper's oral ACAT absorption, its saturable CYP3A4 kinetics, the",
    "CYP3A4 induction (EC50/Emax) layer and the 28 victim-drug models are",
    "not included (see the vignette). Deterministic typical-value model with",
    "no between-subject or residual variability.",
    sep = " "
  )
  reference <- paste(
    "Yang CG, Chen T, Si WT, Wang AH, Ren HC, Wang L. High-performance PBPK",
    "model for predicting CYP3A4",
    "induction-mediated drug interactions: a refined and validated approach.",
    "Front Pharmacol. 2025;16:1521068. doi:10.3389/fphar.2025.1521068.",
    "Disposition parameters are from Table 1 and the deposited GastroPlus",
    "compound database (Supplementary Material 2, 'Rifampicin GastroPlus PK",
    "models.mdb', table Drugs, record 'Rifampicin-ComPK').",
    sep = " "
  )
  vignette <- "Yang_2025_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("central")

  compartmentData <- list(
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rifampicin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "GastroPlus compartmental parameters are entered per kg and scale",
        "linearly with the simulated patient's weight; the deposited",
        "'Rifampicin-ComPK' record uses PatientWeight = 70 kg."
      ),
      source_name = "PatientWeight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "Adults (literature mean profiles)",
    weight_range = "70 kg reference individual (deposited compound record)",
    disease_state = "Healthy volunteers",
    dose_range = "600 mg IV (this model); 200 and 400 mg single oral and 600 mg QD oral (paper, not included)",
    regions = "Not reported",
    notes = paste(
      "No individual data were fitted. The disposition parameters (Vc, K12,",
      "K21) were calculated with the GastroPlus PKPlus module from published",
      "mean plasma profiles after intravenous rifampicin, and the clearance",
      "arms come from the 'Rifampicin-PBPK' record of the GastroPlus DDI",
      "standards database (Table 1). The model was checked against mean",
      "profiles after 600 mg IV, 200 and 400 mg single oral doses and 600 mg",
      "QD oral dosing (Figure 2, Table 2; Acocella 1978, Loos 1985, Peloquin",
      "1997).",
      sep = " "
    )
  )

  ini({
    lvc <- fixed(log(0.15145)); label("Central volume per kg body weight (L/kg)") # Table 1 'Vc (L/Kg)' = 0.15145; mdb Drugs 'Rifampicin-ComPK' VolumeDistribution = 0.15145
    lk12 <- fixed(log(1.4219)); label("Central-to-peripheral rate constant K12 (1/h)") # Table 1 'K12 (1/h)' = 1.4219; mdb Drugs 'Rifampicin-ComPK' Kcp = 1.4219
    lk21 <- fixed(log(1.7949)); label("Peripheral-to-central rate constant K21 (1/h)") # Table 1 'K21 (1/h)' = 1.7949; mdb Drugs 'Rifampicin-ComPK' Kpc = 1.7949
    lcl <- fixed(log(0.0971)); label("Total systemic clearance per kg body weight (L/h/kg)") # Methods 'PBPK model of rifampicin': sum of the three arms 0.0179 + 0.0623 + 0.0169 = 0.0971
    lcl_3a4 <- fixed(log(0.0179)); label("CYP3A4-mediated hepatic clearance per kg body weight (L/h/kg)") # Methods 'PBPK model of rifampicin': 'hepatic non-saturated metabolism mediated by CYP3A4 (about 0.0179 L/h/kg)'
    lcl_renal <- fixed(log(0.0169)); label("Renal clearance per kg body weight (L/h/kg)") # Table 1 'Renal CL (L/h/kg)' = 0.0169; mdb Drugs 'Rifampicin-ComPK' CLr = 0.0169
  })

  model({
    vc <- exp(lvc) * WT
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Clearance arms, each per kg and scaled by body weight. The residual arm
    # is the paper's remaining hepatic metabolism, Table 1 'Other liver CL'
    # = 0.0623 L/h/kg (mdb CLperKg = 0.0623).
    cl <- exp(lcl) * WT
    cl_3a4 <- exp(lcl_3a4) * WT
    cl_renal <- exp(lcl_renal) * WT
    cl_other <- cl - cl_3a4 - cl_renal
    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma concentration, mg/L (= ug/mL).
    Cc <- central / vc
  })
}
