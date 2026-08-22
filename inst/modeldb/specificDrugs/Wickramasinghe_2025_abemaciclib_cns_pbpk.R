Wickramasinghe_2025_abemaciclib_cns_pbpk <- function() {
  description <- paste(
    "PBPK (9-compartment permeability-limited CNS model, SpatialCNS-PBPK",
    "R/Shiny app v1.0).",
    "Spatial pharmacokinetics of abemaciclib in the human central nervous",
    "system and brain tumors. Nine concentration states: brain blood",
    "(brain_vascular), brain parenchyma adjacent to the CSF tract",
    "(brain_csf_adjacent) and deep brain parenchyma (brain_deep), three",
    "CSF compartments (ventricular, cranial subarachnoid, spinal",
    "subarachnoid) and three tumor compartments (infiltrative rim, bulk",
    "tumor, core). Transport across the BBB, the blood-brain-tumor",
    "barrier and the blood-CSF barrier is the sum of a passive",
    "permeability-surface-area clearance acting on unbound AND unionized",
    "drug (the pH-dependent unionization fractions lam* encode the",
    "progressively acidic tumor interior) and transporter-mediated efflux",
    "and uptake clearances acting on unbound drug. CSF flows from the",
    "ventricles to the cranial and spinal subarachnoid spaces and drains",
    "to blood via arachnoid villi and nerve sheaths; drug also moves",
    "between the CSF tract and brain/tumor by paravascular (glymphatic)",
    "bulk flow and simple diffusion. The model is driven by the systemic",
    "plasma concentration supplied as the time-varying covariate",
    "CP_ABEMACICLIB_UM: plasma acts purely as a forcing function and CNS",
    "uptake does not deplete it, exactly as published. Every one of the 85",
    "system- and drug-specific parameters is fixed to the published value;",
    "none is estimated here."
  )
  reference <- paste(
    "Wickramasinghe CD, Kim S, Li J. SpatialCNS-PBPK: An R/Shiny",
    "Web-Based Application for Physiologically Based Pharmacokinetic",
    "Modeling of Spatial Pharmacokinetics in the Human Central Nervous",
    "System and Brain Tumors. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(5):864-875. doi:10.1002/psp4.70026.",
    "The nine differential equations are given in Data S1 of the",
    "Supporting Information; the system-specific parameters in Table 1;",
    "the drug-specific parameter definitions in Table 2 and their",
    "abemaciclib values, the interindividual variability and the driving",
    "plasma profile in Table S2. The 9-CNS model was originally developed",
    "and validated across six drugs in Li J, Wickramasinghe C, Jiang J,",
    "et al. Clin Pharmacol Ther. 2025;117(3):690-703 (reference 7 of the",
    "2025 tutorial), which is not on disk for this extraction."
  )
  vignette <- "Wickramasinghe_2025_spatial_cns_pbpk"

  units <- list(
    time          = "h",
    dosing        = "n/a (no dose events; plasma enters as the CP_ABEMACICLIB_UM covariate)",
    concentration = "umol/L"
  )

  covariateData <- list(
    CP_ABEMACICLIB_UM = list(
      description        = paste(
        "Instantaneous total abemaciclib plasma concentration, supplied as a",
        "time-varying covariate column. This is the input function of the",
        "9-CNS model: it drives cerebral delivery through the Qbrain term and",
        "is never depleted by CNS uptake."
      ),
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table S2 supplies a 1010-point population-mean profile over 0-168 h",
        "for glioblastoma patients receiving twice-daily oral abemaciclib",
        "(peak 0.4191 at 124.61 h; median inter-peak interval 11.4 h). The",
        "dose amount is not stated in this paper. Interpolate linearly",
        "between the supplied times, as the app does",
        "(rxSolve(covsInterpolation = \"linear\"))."
      ),
      source_name        = "Plasma1"
    )
  )

  compartmentData <- list(
    brain_vascular        = list(analyte = "abemaciclib", units = "umol/L", specimen = "plasma", verified = TRUE),
    brain_csf_adjacent    = list(analyte = "abemaciclib", units = "umol/L", specimen = "tissue", verified = TRUE),
    brain_deep            = list(analyte = "abemaciclib", units = "umol/L", specimen = "tissue", verified = TRUE),
    tumor_rim             = list(analyte = "abemaciclib", units = "umol/L", specimen = "tumor",  verified = TRUE),
    tumor_bulk            = list(analyte = "abemaciclib", units = "umol/L", specimen = "tumor",  verified = TRUE),
    tumor_core            = list(analyte = "abemaciclib", units = "umol/L", specimen = "tumor",  verified = TRUE),
    brain_csf_ventricular = list(analyte = "abemaciclib", units = "umol/L", specimen = "CSF",    verified = TRUE),
    brain_csf_sas_cranial = list(analyte = "abemaciclib", units = "umol/L", specimen = "CSF",    verified = TRUE),
    brain_csf_sas_spinal  = list(analyte = "abemaciclib", units = "umol/L", specimen = "CSF",    verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 39,
    disease_state = "glioblastoma (recurrent / newly diagnosed) undergoing tumor resection",
    dose_range    = "twice-daily oral abemaciclib; the dose amount is not stated in this paper",
    notes         = paste(
      "The driving plasma profile in Table S2 is the population-mean",
      "abemaciclib concentration-time profile determined in glioblastoma",
      "patients and is attributed to reference 7 (Li 2025 Clin Pharmacol",
      "Ther), which is not on disk. Table S2 also carries paired observed",
      "abemaciclib concentrations for 39 patients (subject IDs 100-0049 to",
      "300-0013): 39 non-enhancing tumor, 37 contrast-enhancing tumor and 34",
      "CSF samples (110 in total), collected 120.0-130.1 h after the start of",
      "dosing. No",
      "demographics (age, weight, sex, race) are reported in this paper.",
      "System-specific parameters in Table 1 describe a reference adult",
      "brain (1.2 L parenchyma, 0.15 L total CSF) carrying a 35 mL",
      "contrast-enhancing tumor."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All 85 parameters are FIXED to the published values: the 9-CNS model is
    # parameterised from physiology (Table 1) and from in vitro-in vivo
    # extrapolation of abemaciclib permeability / transporter / binding data
    # (Table 2, values in Table S2). Nothing is estimated in this paper.
    #
    # Parameters with a non-zero typical value are carried on the log scale so
    # the interindividual variability below is lognormal. The 13 parameters
    # that the paper sets to exactly zero (no measurable uptake transport, no
    # brain metabolism, no cranial-to-spinal CSF back-flow) cannot be
    # log-transformed and carry no eta -- matching the app, whose own virtual
    # patients (Table S3) return exactly zero for every one of them.
    # ---------------------------------------------------------------------

    lVbb <- fixed(log(0.063)); label("Volume of brain blood (L)")  # Table 1 + Table S2 'Sim1'
    lVbm1 <- fixed(log(0.12)); label("Volume of brain parenchyma adjacent to the CSF tract (L)")  # Table 1 + Table S2 'Sim1'
    lVbm2 <- fixed(log(1.08)); label("Volume of deep brain parenchyma (L)")  # Table 1 + Table S2 'Sim1'
    lVT1 <- fixed(log(0.07)); label("Volume of the infiltrative tumor rim (L)")  # Table 1 + Table S2 'Sim1'
    lVT2 <- fixed(log(0.035)); label("Volume of the bulky tumor region (L)")  # Table 1 + Table S2 'Sim1'
    lVT3 <- fixed(log(0.0035)); label("Volume of the tumor core (L)")  # Table 1 + Table S2 'Sim1'
    lVvcsf <- fixed(log(0.025)); label("Volume of ventricular CSF (L)")  # Table 1 + Table S2 'Sim1'
    lVccsf <- fixed(log(0.045)); label("Volume of cranial subarachnoid CSF (L)")  # Table 1 + Table S2 'Sim1'
    lVscsf <- fixed(log(0.08)); label("Volume of spinal subarachnoid CSF (L)")  # Table 1 + Table S2 'Sim1'
    lQbrain <- fixed(log(39)); label("Cerebral blood flow (L/h)")  # Table 1 + Table S2 'Sim1'
    lQCsink <- fixed(log(0.013)); label("Absorption rate of cranial CSF into blood via arachnoid villi (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSsink <- fixed(log(0.008)); label("Absorption rate of spinal CSF into blood via arachnoid villi (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkCB1 <- fixed(log(0.0013)); label("Paravascular bulk flow rate, cranial subarachnoid CSF to brain parenchyma 1 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB1C <- fixed(log(0.0016)); label("Paravascular bulk flow rate, brain parenchyma 1 to cranial subarachnoid CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkVB1 <- fixed(log(0.0001)); label("Paravascular bulk flow rate, ventricular CSF to brain parenchyma 1 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB1V <- fixed(log(0.0002)); label("Paravascular bulk flow rate, brain parenchyma 1 to ventricular CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB2B1 <- fixed(log(0.0005)); label("Convective bulk flow rate, brain parenchyma 2 to 1 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB1B2 <- fixed(log(0.0005)); label("Convective bulk flow rate, brain parenchyma 1 to 2 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB2T1 <- fixed(log(0.0005)); label("Convective bulk flow rate, brain parenchyma 2 to tumor rim (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT1B2 <- fixed(log(0.0005)); label("Convective bulk flow rate, tumor rim to brain parenchyma 2 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkCB2 <- fixed(log(0.0113)); label("Paravascular bulk flow rate, cranial subarachnoid CSF to brain parenchyma 2 (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkB2C <- fixed(log(0.0142)); label("Paravascular bulk flow rate, brain parenchyma 2 to cranial subarachnoid CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT2T1 <- fixed(log(0.0002)); label("Convective bulk flow rate, bulk tumor to tumor rim (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT1T2 <- fixed(log(0.0002)); label("Convective bulk flow rate, tumor rim to bulk tumor (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkCT1 <- fixed(log(0.0013)); label("Paravascular bulk flow rate, cranial subarachnoid CSF to tumor rim (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT1C <- fixed(log(0.0016)); label("Paravascular bulk flow rate, tumor rim to cranial subarachnoid CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT3T2 <- fixed(log(0.0002)); label("Convective bulk flow rate, tumor core to bulk tumor (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT2T3 <- fixed(log(0.0002)); label("Convective bulk flow rate, bulk tumor to tumor core (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkCT2 <- fixed(log(0.0002)); label("Paravascular bulk flow rate, cranial subarachnoid CSF to bulk tumor (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT2C <- fixed(log(0.0003)); label("Paravascular bulk flow rate, bulk tumor to cranial subarachnoid CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkCT3 <- fixed(log(0.0002)); label("Paravascular bulk flow rate, cranial subarachnoid CSF to tumor core (L/h)")  # Table 1 + Table S2 'Sim1'
    lQbulkT3C <- fixed(log(0.0002)); label("Paravascular bulk flow rate, tumor core to cranial subarachnoid CSF (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSin1r <- fixed(log(0.0013)); label("CSF back-flow rate from the cranial subarachnoid space to the ventricle (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSin1 <- fixed(log(0.0126)); label("CSF flow rate from the ventricle to the cranial subarachnoid space (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSin2r <- fixed(log(0.0008)); label("CSF back-flow rate from the spinal subarachnoid space to the ventricle (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSin2 <- fixed(log(0.0084)); label("CSF flow rate from the ventricle to the spinal subarachnoid space (L/h)")  # Table 1 + Table S2 'Sim1'
    lQSout <- fixed(log(0.0004)); label("CSF flow rate from the spinal to the cranial subarachnoid space (L/h)")  # Table 1 + Table S2 'Sim1'
    QSoutr <- fixed(0); label("CSF flow rate from the cranial to the spinal subarachnoid space (L/h)")  # Table 1 + Table S2 'Sim1'
    lPSB1 <- fixed(log(2.548)); label("Passive permeability clearance at the BBB, brain blood to adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSB2 <- fixed(log(25.48)); label("Passive permeability clearance at the BBB, brain blood to deep brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lPST1 <- fixed(log(2.0384)); label("Passive permeability clearance at the BBTB, brain blood to tumor rim (L/h)")  # Table 2 + Table S2 'Sim1'
    lPST2 <- fixed(log(5.096)); label("Passive permeability clearance at the BBTB, brain blood to bulk tumor (L/h)")  # Table 2 + Table S2 'Sim1'
    lPST3 <- fixed(log(0.5096)); label("Passive permeability clearance at the BBTB, brain blood to tumor core (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSV <- fixed(log(0.8408)); label("Passive permeability clearance at the blood-CSF barrier, brain blood to ventricular CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSC <- fixed(log(1.7072)); label("Passive permeability clearance at the blood-CSF barrier, brain blood to cranial subarachnoid CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSE1 <- fixed(log(50.96)); label("Simple diffusion rate between cranial subarachnoid CSF and adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSE2 <- fixed(log(25.48)); label("Simple diffusion rate between ventricular CSF and adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSB1B2 <- fixed(log(0.01)); label("Simple diffusion rate between the adjacent and deep brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lPSB2T1 <- fixed(log(0.01)); label("Simple diffusion rate between the deep brain parenchyma and tumor rim (L/h)")  # Table 2 + Table S2 'Sim1'
    lPST1T2 <- fixed(log(0.01)); label("Simple diffusion rate between tumor rim and bulk tumor (L/h)")  # Table 2 + Table S2 'Sim1'
    lPST2T3 <- fixed(log(0.01)); label("Simple diffusion rate between bulk tumor and tumor core (L/h)")  # Table 2 + Table S2 'Sim1'
    lCLeffbbb1 <- fixed(log(0.292)); label("Efflux transporter-mediated clearance at the BBB, brain blood to adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupbbb1 <- fixed(0); label("Uptake transporter-mediated influx clearance at the BBB, brain blood to adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lCLeffbbb2 <- fixed(log(2.92)); label("Efflux transporter-mediated clearance at the BBB, brain blood to deep brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupbbb2 <- fixed(0); label("Uptake transporter-mediated influx clearance at the BBB, brain blood to deep brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lCLeffT1 <- fixed(log(0.0584)); label("Efflux transporter-mediated clearance at the BBTB, brain blood to tumor rim (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupT1 <- fixed(0); label("Uptake transporter-mediated influx clearance at the BBB, brain blood to tumor rim (L/h)")  # Table 2 + Table S2 'Sim1'
    lCLeffT2 <- fixed(log(0.0088)); label("Efflux transporter-mediated clearance at the BBTB, brain blood to bulk tumor (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupT2 <- fixed(0); label("Uptake transporter-mediated influx clearance at the BBB, brain blood to bulk tumor (L/h)")  # Table 2 + Table S2 'Sim1'
    CLeffT3 <- fixed(0); label("Efflux transporter-mediated clearance at the BBTB, brain blood to tumor core (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupT3 <- fixed(0); label("Uptake transporter-mediated influx clearance at the BBB, brain blood to tumor core (L/h)")  # Table 2 + Table S2 'Sim1'
    CLeffvcsf <- fixed(0); label("Efflux transporter-mediated clearance at the blood-CSF barrier, brain blood to ventricular CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupvcsf <- fixed(0); label("Uptake transporter-mediated clearance at the blood-CSF barrier, brain blood to ventricular CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    CLeffccsf <- fixed(0); label("Efflux transporter-mediated clearance at the blood-CSF barrier, brain blood to cranial subarachnoid CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    CLupccsf <- fixed(0); label("Uptake transporter-mediated clearance at the blood-CSF barrier, brain blood to cranial subarachnoid CSF (L/h)")  # Table 2 + Table S2 'Sim1'
    CLmet1 <- fixed(0); label("Drug metabolism clearance in the adjacent brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    CLmet2 <- fixed(0); label("Drug metabolism clearance in the deep brain parenchyma (L/h)")  # Table 2 + Table S2 'Sim1'
    lfubb <- fixed(log(0.026)); label("Unbound fraction in brain blood (fraction)")  # Table 2 + Table S2 'Sim1'
    lfubm1 <- fixed(log(0.006)); label("Unbound fraction in the adjacent brain parenchyma (fraction)")  # Table 2 + Table S2 'Sim1'
    lfubm2 <- fixed(log(0.006)); label("Unbound fraction in the deep brain parenchyma (fraction)")  # Table 2 + Table S2 'Sim1'
    lfuT1 <- fixed(log(0.009)); label("Unbound fraction in the tumor rim (fraction)")  # Table 2 + Table S2 'Sim1'
    lfuT2 <- fixed(log(0.016)); label("Unbound fraction in the bulk tumor (fraction)")  # Table 2 + Table S2 'Sim1'
    lfuT3 <- fixed(log(0.016)); label("Unbound fraction in the tumor core (fraction)")  # Table 2 + Table S2 'Sim1'
    lfuvcsf <- fixed(log(0.289)); label("Unbound fraction in ventricular CSF (fraction)")  # Table 2 + Table S2 'Sim1'
    lfuccsf <- fixed(log(0.289)); label("Unbound fraction in cranial subarachnoid CSF (fraction)")  # Table 2 + Table S2 'Sim1'
    llambb <- fixed(log(0.2238)); label("Unionization fraction in brain blood (pH 7.4) (fraction)")  # Table 2 + Table S2 'Sim1'
    llambm1 <- fixed(log(0.154)); label("Unionization fraction in the adjacent brain parenchyma (pH 7.2) (fraction)")  # Table 2 + Table S2 'Sim1'
    llambm2 <- fixed(log(0.154)); label("Unionization fraction in the deep brain parenchyma (pH 7.2) (fraction)")  # Table 2 + Table S2 'Sim1'
    llamT1 <- fixed(log(0.0676)); label("Unionization fraction in the tumor rim (pH 6.8) (fraction)")  # Table 2 + Table S2 'Sim1'
    llamT2 <- fixed(log(0.035)); label("Unionization fraction in the bulk tumor (pH 6.5) (fraction)")  # Table 2 + Table S2 'Sim1'
    llamT3 <- fixed(log(0.0114)); label("Unionization fraction in the tumor core (pH 6.2) (fraction)")  # Table 2 + Table S2 'Sim1'
    llamvcsf <- fixed(log(0.1864)); label("Unionization fraction in ventricular CSF (pH 7.3) (fraction)")  # Table 2 + Table S2 'Sim1'
    llamccsf <- fixed(log(0.1864)); label("Unionization fraction in cranial subarachnoid CSF (pH 7.3) (fraction)")  # Table 2 + Table S2 'Sim1'
    lQglyccsf <- fixed(log(0.0065)); label("Absorption rate of cranial CSF via olfactory mucosa and cranial nerve sheaths (L/h)")  # Table 1 + Table S2 'Sim1'
    lQglyscsf <- fixed(log(0.004)); label("Absorption rate of spinal CSF via spinal nerve sheaths (L/h)")  # Table 1 + Table S2 'Sim1'

    # ---------------------------------------------------------------------
    # Interindividual variability. Table S2 assigns the SAME 20% coefficient
    # of variation ("IIV1") to every one of the 85 parameters, so
    # var = log(1 + 0.20^2) = 0.03922071 on the log scale, uncorrelated.
    #
    # The lognormal form is not stated in the paper; it was identified from
    # the app's own output. Table S3 lists five virtual patients that the
    # SpatialCNS-PBPK app generated from Sim1 + IIV1. Across the 355 non-zero
    # parameter draws, log(simulated / typical) is normal (Shapiro-Wilk
    # p = 0.15, skewness 0.12) with SD 0.194 -- against sqrt(log(1 + 0.20^2))
    # = 0.198 -- while the untransformed ratios are clearly right-skewed
    # (Shapiro-Wilk p = 1.6e-8, skewness 0.94). Lognormal at 20% CV therefore
    # reproduces the app; a normal distribution does not.
    #
    # Note that this uniform 20% is the input-file TEMPLATE default rather
    # than an estimated variance -- see the vignette Errata.
    # ---------------------------------------------------------------------

    etalVbb ~ 0.03922071
    etalVbm1 ~ 0.03922071
    etalVbm2 ~ 0.03922071
    etalVT1 ~ 0.03922071
    etalVT2 ~ 0.03922071
    etalVT3 ~ 0.03922071
    etalVvcsf ~ 0.03922071
    etalVccsf ~ 0.03922071
    etalVscsf ~ 0.03922071
    etalQbrain ~ 0.03922071
    etalQCsink ~ 0.03922071
    etalQSsink ~ 0.03922071
    etalQbulkCB1 ~ 0.03922071
    etalQbulkB1C ~ 0.03922071
    etalQbulkVB1 ~ 0.03922071
    etalQbulkB1V ~ 0.03922071
    etalQbulkB2B1 ~ 0.03922071
    etalQbulkB1B2 ~ 0.03922071
    etalQbulkB2T1 ~ 0.03922071
    etalQbulkT1B2 ~ 0.03922071
    etalQbulkCB2 ~ 0.03922071
    etalQbulkB2C ~ 0.03922071
    etalQbulkT2T1 ~ 0.03922071
    etalQbulkT1T2 ~ 0.03922071
    etalQbulkCT1 ~ 0.03922071
    etalQbulkT1C ~ 0.03922071
    etalQbulkT3T2 ~ 0.03922071
    etalQbulkT2T3 ~ 0.03922071
    etalQbulkCT2 ~ 0.03922071
    etalQbulkT2C ~ 0.03922071
    etalQbulkCT3 ~ 0.03922071
    etalQbulkT3C ~ 0.03922071
    etalQSin1r ~ 0.03922071
    etalQSin1 ~ 0.03922071
    etalQSin2r ~ 0.03922071
    etalQSin2 ~ 0.03922071
    etalQSout ~ 0.03922071
    etalPSB1 ~ 0.03922071
    etalPSB2 ~ 0.03922071
    etalPST1 ~ 0.03922071
    etalPST2 ~ 0.03922071
    etalPST3 ~ 0.03922071
    etalPSV ~ 0.03922071
    etalPSC ~ 0.03922071
    etalPSE1 ~ 0.03922071
    etalPSE2 ~ 0.03922071
    etalPSB1B2 ~ 0.03922071
    etalPSB2T1 ~ 0.03922071
    etalPST1T2 ~ 0.03922071
    etalPST2T3 ~ 0.03922071
    etalCLeffbbb1 ~ 0.03922071
    etalCLeffbbb2 ~ 0.03922071
    etalCLeffT1 ~ 0.03922071
    etalCLeffT2 ~ 0.03922071
    etalfubb ~ 0.03922071
    etalfubm1 ~ 0.03922071
    etalfubm2 ~ 0.03922071
    etalfuT1 ~ 0.03922071
    etalfuT2 ~ 0.03922071
    etalfuT3 ~ 0.03922071
    etalfuvcsf ~ 0.03922071
    etalfuccsf ~ 0.03922071
    etallambb ~ 0.03922071
    etallambm1 ~ 0.03922071
    etallambm2 ~ 0.03922071
    etallamT1 ~ 0.03922071
    etallamT2 ~ 0.03922071
    etallamT3 ~ 0.03922071
    etallamvcsf ~ 0.03922071
    etallamccsf ~ 0.03922071
    etalQglyccsf ~ 0.03922071
    etalQglyscsf ~ 0.03922071

  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual parameter values (lognormal IIV, 20% CV on each).
    # -------------------------------------------------------------------

    Vbb <- exp(lVbb + etalVbb)
    Vbm1 <- exp(lVbm1 + etalVbm1)
    Vbm2 <- exp(lVbm2 + etalVbm2)
    VT1 <- exp(lVT1 + etalVT1)
    VT2 <- exp(lVT2 + etalVT2)
    VT3 <- exp(lVT3 + etalVT3)
    Vvcsf <- exp(lVvcsf + etalVvcsf)
    Vccsf <- exp(lVccsf + etalVccsf)
    Vscsf <- exp(lVscsf + etalVscsf)
    Qbrain <- exp(lQbrain + etalQbrain)
    QCsink <- exp(lQCsink + etalQCsink)
    QSsink <- exp(lQSsink + etalQSsink)
    QbulkCB1 <- exp(lQbulkCB1 + etalQbulkCB1)
    QbulkB1C <- exp(lQbulkB1C + etalQbulkB1C)
    QbulkVB1 <- exp(lQbulkVB1 + etalQbulkVB1)
    QbulkB1V <- exp(lQbulkB1V + etalQbulkB1V)
    QbulkB2B1 <- exp(lQbulkB2B1 + etalQbulkB2B1)
    QbulkB1B2 <- exp(lQbulkB1B2 + etalQbulkB1B2)
    QbulkB2T1 <- exp(lQbulkB2T1 + etalQbulkB2T1)
    QbulkT1B2 <- exp(lQbulkT1B2 + etalQbulkT1B2)
    QbulkCB2 <- exp(lQbulkCB2 + etalQbulkCB2)
    QbulkB2C <- exp(lQbulkB2C + etalQbulkB2C)
    QbulkT2T1 <- exp(lQbulkT2T1 + etalQbulkT2T1)
    QbulkT1T2 <- exp(lQbulkT1T2 + etalQbulkT1T2)
    QbulkCT1 <- exp(lQbulkCT1 + etalQbulkCT1)
    QbulkT1C <- exp(lQbulkT1C + etalQbulkT1C)
    QbulkT3T2 <- exp(lQbulkT3T2 + etalQbulkT3T2)
    QbulkT2T3 <- exp(lQbulkT2T3 + etalQbulkT2T3)
    QbulkCT2 <- exp(lQbulkCT2 + etalQbulkCT2)
    QbulkT2C <- exp(lQbulkT2C + etalQbulkT2C)
    QbulkCT3 <- exp(lQbulkCT3 + etalQbulkCT3)
    QbulkT3C <- exp(lQbulkT3C + etalQbulkT3C)
    QSin1r <- exp(lQSin1r + etalQSin1r)
    QSin1 <- exp(lQSin1 + etalQSin1)
    QSin2r <- exp(lQSin2r + etalQSin2r)
    QSin2 <- exp(lQSin2 + etalQSin2)
    QSout <- exp(lQSout + etalQSout)
    PSB1 <- exp(lPSB1 + etalPSB1)
    PSB2 <- exp(lPSB2 + etalPSB2)
    PST1 <- exp(lPST1 + etalPST1)
    PST2 <- exp(lPST2 + etalPST2)
    PST3 <- exp(lPST3 + etalPST3)
    PSV <- exp(lPSV + etalPSV)
    PSC <- exp(lPSC + etalPSC)
    PSE1 <- exp(lPSE1 + etalPSE1)
    PSE2 <- exp(lPSE2 + etalPSE2)
    PSB1B2 <- exp(lPSB1B2 + etalPSB1B2)
    PSB2T1 <- exp(lPSB2T1 + etalPSB2T1)
    PST1T2 <- exp(lPST1T2 + etalPST1T2)
    PST2T3 <- exp(lPST2T3 + etalPST2T3)
    CLeffbbb1 <- exp(lCLeffbbb1 + etalCLeffbbb1)
    CLeffbbb2 <- exp(lCLeffbbb2 + etalCLeffbbb2)
    CLeffT1 <- exp(lCLeffT1 + etalCLeffT1)
    CLeffT2 <- exp(lCLeffT2 + etalCLeffT2)
    fubb <- exp(lfubb + etalfubb)
    fubm1 <- exp(lfubm1 + etalfubm1)
    fubm2 <- exp(lfubm2 + etalfubm2)
    fuT1 <- exp(lfuT1 + etalfuT1)
    fuT2 <- exp(lfuT2 + etalfuT2)
    fuT3 <- exp(lfuT3 + etalfuT3)
    fuvcsf <- exp(lfuvcsf + etalfuvcsf)
    fuccsf <- exp(lfuccsf + etalfuccsf)
    lambb <- exp(llambb + etallambb)
    lambm1 <- exp(llambm1 + etallambm1)
    lambm2 <- exp(llambm2 + etallambm2)
    lamT1 <- exp(llamT1 + etallamT1)
    lamT2 <- exp(llamT2 + etallamT2)
    lamT3 <- exp(llamT3 + etallamT3)
    lamvcsf <- exp(llamvcsf + etallamvcsf)
    lamccsf <- exp(llamccsf + etallamccsf)
    Qglyccsf <- exp(lQglyccsf + etalQglyccsf)
    Qglyscsf <- exp(lQglyscsf + etalQglyscsf)

    # -------------------------------------------------------------------
    # 2. Input function. The systemic plasma concentration is supplied
    #    exogenously and is NOT depleted by CNS uptake, exactly as
    #    published ("the drug plasma concentration-time profile serving as
    #    the input function of the 9-CNS model").
    # -------------------------------------------------------------------
    Cart <- CP_ABEMACICLIB_UM

    # -------------------------------------------------------------------
    # 3. Local aliases so the ODEs below read exactly as Data S1 prints
    #    them. Each state holds a TOTAL drug concentration; the unbound and
    #    unionized driving concentration on each side of a barrier is
    #    lam<region> * fu<region> * C<region> for passive permeability, and
    #    fu<region> * C<region> for transporter-mediated clearance (which
    #    acts on unbound drug whether ionized or not).
    # -------------------------------------------------------------------
    Cbb   <- brain_vascular
    Cbm1  <- brain_csf_adjacent
    Cbm2  <- brain_deep
    CT1   <- tumor_rim
    CT2   <- tumor_bulk
    CT3   <- tumor_core
    Cvcsf <- brain_csf_ventricular
    Cccsf <- brain_csf_sas_cranial
    Cscsf <- brain_csf_sas_spinal

    # -------------------------------------------------------------------
    # 4. The nine differential equations, Data S1. Each is written as
    #    V * dC/dt = (fluxes) in the source and divided by V here.
    # -------------------------------------------------------------------

    # Brain blood
    d/dt(brain_vascular) <- (
      Qbrain * (Cart - Cbb) +
      PSB1 * (lambm1 * fubm1 * Cbm1 - lambb * fubb * Cbb) +
        CLeffbbb1 * fubm1 * Cbm1 - CLupbbb1 * fubb * Cbb +
      PSB2 * (lambm2 * fubm2 * Cbm2 - lambb * fubb * Cbb) +
        CLeffbbb2 * fubm2 * Cbm2 - CLupbbb2 * fubb * Cbb +
      PST1 * (lamT1 * fuT1 * CT1 - lambb * fubb * Cbb) +
        CLeffT1 * fuT1 * CT1 - CLupT1 * fubb * Cbb +
      PST2 * (lamT2 * fuT2 * CT2 - lambb * fubb * Cbb) +
        CLeffT2 * fuT2 * CT2 - CLupT2 * fubb * Cbb +
      PST3 * (lamT3 * fuT3 * CT3 - lambb * fubb * Cbb) +
        CLeffT3 * fuT3 * CT3 - CLupT3 * fubb * Cbb +
      PSV * (lamvcsf * fuvcsf * Cvcsf - lambb * fubb * Cbb) +
        CLeffvcsf * fuvcsf * Cvcsf - CLupvcsf * fubb * Cbb +
      PSC * (lamccsf * fuccsf * Cccsf - lambb * fubb * Cbb) +
        CLeffccsf * fuccsf * Cccsf - CLupccsf * fubb * Cbb +
      QCsink * Cccsf + QSsink * Cscsf
    ) / Vbb

    # Brain parenchyma 1 (brain tissue adjacent to the CSF tract)
    d/dt(brain_csf_adjacent) <- (
      PSB1 * (lambb * fubb * Cbb - lambm1 * fubm1 * Cbm1) -
        CLeffbbb1 * fubm1 * Cbm1 + CLupbbb1 * fubb * Cbb +
      PSE1 * (lamccsf * fuccsf * Cccsf - lambm1 * fubm1 * Cbm1) +
        QbulkCB1 * Cccsf - QbulkB1C * fubm1 * Cbm1 +
      PSE2 * (lamvcsf * fuvcsf * Cvcsf - lambm1 * fubm1 * Cbm1) +
        QbulkVB1 * Cvcsf - QbulkB1V * fubm1 * Cbm1 +
      PSB1B2 * (lambm2 * fubm2 * Cbm2 - lambm1 * fubm1 * Cbm1) +
        QbulkB2B1 * fubm2 * Cbm2 - QbulkB1B2 * fubm1 * Cbm1 -
      CLmet1 * fubm1 * Cbm1
    ) / Vbm1

    # Brain parenchyma 2 (deep brain parenchyma)
    d/dt(brain_deep) <- (
      PSB2 * (lambb * fubb * Cbb - lambm2 * fubm2 * Cbm2) -
        CLeffbbb2 * fubm2 * Cbm2 + CLupbbb2 * fubb * Cbb +
      PSB1B2 * (lambm1 * fubm1 * Cbm1 - lambm2 * fubm2 * Cbm2) +
      PSB2T1 * (lamT1 * fuT1 * CT1 - lambm2 * fubm2 * Cbm2) -
        QbulkB2B1 * fubm2 * Cbm2 + QbulkB1B2 * fubm1 * Cbm1 -
        QbulkB2T1 * fubm2 * Cbm2 + QbulkT1B2 * fuT1 * CT1 +
        QbulkCB2 * Cccsf - QbulkB2C * fubm2 * Cbm2 -
      CLmet2 * fubm2 * Cbm2
    ) / Vbm2

    # Tumor mass 1 (infiltrative tumor rim)
    d/dt(tumor_rim) <- (
      PST1 * (lambb * fubb * Cbb - lamT1 * fuT1 * CT1) -
        CLeffT1 * fuT1 * CT1 + CLupT1 * fubb * Cbb +
      PSB2T1 * (lambm2 * fubm2 * Cbm2 - lamT1 * fuT1 * CT1) +
      PST1T2 * (lamT2 * fuT2 * CT2 - lamT1 * fuT1 * CT1) +
        QbulkB2T1 * fubm2 * Cbm2 - QbulkT1B2 * fuT1 * CT1 +
        QbulkT2T1 * fuT2 * CT2 - QbulkT1T2 * fuT1 * CT1 +
        QbulkCT1 * Cccsf - QbulkT1C * fuT1 * CT1
    ) / VT1

    # Tumor mass 2 (bulk tumor)
    d/dt(tumor_bulk) <- (
      PST2 * (lambb * fubb * Cbb - lamT2 * fuT2 * CT2) -
        CLeffT2 * fuT2 * CT2 + CLupT2 * fubb * Cbb +
      PST1T2 * (lamT1 * fuT1 * CT1 - lamT2 * fuT2 * CT2) +
      PST2T3 * (lamT3 * fuT3 * CT3 - lamT2 * fuT2 * CT2) +
        QbulkT1T2 * fuT1 * CT1 - QbulkT2T1 * fuT2 * CT2 +
        QbulkT3T2 * fuT3 * CT3 - QbulkT2T3 * fuT2 * CT2 +
        QbulkCT2 * Cccsf - QbulkT2C * fuT2 * CT2
    ) / VT2

    # Tumor mass 3 (tumor core)
    d/dt(tumor_core) <- (
      PST3 * (lambb * fubb * Cbb - lamT3 * fuT3 * CT3) -
        CLeffT3 * fuT3 * CT3 + CLupT3 * fubb * Cbb +
      PST2T3 * (lamT2 * fuT2 * CT2 - lamT3 * fuT3 * CT3) +
        QbulkT2T3 * fuT2 * CT2 - QbulkT3T2 * fuT3 * CT3 +
        QbulkCT3 * Cccsf - QbulkT3C * fuT3 * CT3
    ) / VT3

    # Ventricular CSF
    d/dt(brain_csf_ventricular) <- (
      PSV * (lambb * fubb * Cbb - lamvcsf * fuvcsf * Cvcsf) +
        CLupvcsf * fubb * Cbb - CLeffvcsf * fuvcsf * Cvcsf +
      PSE2 * (lambm1 * fubm1 * Cbm1 - lamvcsf * fuvcsf * Cvcsf) +
        QbulkB1V * fubm1 * Cbm1 - QbulkVB1 * Cvcsf +
        QSin1r * Cccsf - QSin1 * Cvcsf +
        QSin2r * Cscsf - QSin2 * Cvcsf
    ) / Vvcsf

    # Cranial subarachnoid CSF
    d/dt(brain_csf_sas_cranial) <- (
      PSE1 * (lambm1 * fubm1 * Cbm1 - lamccsf * fuccsf * Cccsf) -
      PSC * (lamccsf * fuccsf * Cccsf - lambb * fubb * Cbb) -
        CLeffccsf * fuccsf * Cccsf + CLupccsf * fubb * Cbb +
        QbulkB1C * fubm1 * Cbm1 - QbulkCB1 * Cccsf +
        QbulkB2C * fubm2 * Cbm2 - QbulkCB2 * Cccsf -
        QbulkCT1 * Cccsf + QbulkT1C * fuT1 * CT1 -
        QbulkCT2 * Cccsf + QbulkT2C * fuT2 * CT2 -
        QbulkCT3 * Cccsf + QbulkT3C * fuT3 * CT3 +
        QSin1 * Cvcsf + QSout * Cscsf - QSoutr * Cccsf -
        QSin1r * Cccsf - QCsink * Cccsf - Qglyccsf * Cccsf
    ) / Vccsf

    # Spinal subarachnoid CSF
    d/dt(brain_csf_sas_spinal) <- (
      QSin2 * Cvcsf + QSoutr * Cccsf - QSout * Cscsf -
        QSsink * Cscsf - QSin2r * Cscsf - Qglyscsf * Cscsf
    ) / Vscsf

    # -------------------------------------------------------------------
    # 5. Observations. Cbb ... Cscsf above are the nine total-concentration
    #    outputs and carry the paper's own names, so they line up
    #    one-for-one with the panels of Figure 2. The corresponding unbound
    #    concentrations are the app's "(Unbound)" output tabs.
    #
    #    No residual-error model is given: the paper reports no residual
    #    variability and the app is a deterministic simulator whose only
    #    stochastic element is the parameter-level IIV above.
    # -------------------------------------------------------------------
    Cbbu   <- fubb * Cbb
    Cbm1u  <- fubm1 * Cbm1
    Cbm2u  <- fubm2 * Cbm2
    CT1u   <- fuT1 * CT1
    CT2u   <- fuT2 * CT2
    CT3u   <- fuT3 * CT3
    Cvcsfu <- fuvcsf * Cvcsf
    Cccsfu <- fuccsf * Cccsf
  })
}

