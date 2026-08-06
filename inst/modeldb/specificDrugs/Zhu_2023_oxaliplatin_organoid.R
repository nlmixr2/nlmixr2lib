Zhu_2023_oxaliplatin_organoid <- function() {
  description <- paste(
    "In vitro (patient-derived colorectal-cancer tumour organoids, PDTOs).",
    "Exponential organoid-growth model with a sigmoidal Emax (Hill) killing",
    "term describing the 96 h concentration-viability relationship of",
    "oxaliplatin. A vehicle-control organoid and a treated organoid grow in",
    "parallel; the reported readout is cell viability, the ratio of treated to",
    "control organoid volume. This is the in vitro half of the Zhu 2023",
    "in vitro-to-in vivo translation: the same Emax / EC50 / Hill parameters",
    "are carried into the human oxaliplatin PK/PD model, see",
    "modellib('Zhu_2023_oxaliplatin')."
  )
  reference <- paste(
    "Zhu J, Zhang Y, Zhao Y, Zhang J, Hao K, He H. Translational",
    "Pharmacokinetic/Pharmacodynamic Modeling and Simulation of Oxaliplatin and",
    "Irinotecan in Colorectal Cancer. Pharmaceutics. 2023;15(9):2274.",
    "doi:10.3390/pharmaceutics15092274.",
    sep = " "
  )
  vignette <- "Zhu_2023_oxaliplatin_irinotecan_colorectal_cancer"

  paper_specific_compartments <- c("organoid_ctrl", "organoid_trt", "viability")

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    organoid_ctrl = list(
      analyte = "none", units = "arbitrary volume unit",
      specimen = "not applicable", verified = TRUE
    ),
    organoid_trt = list(
      analyte = "none", units = "arbitrary volume unit",
      specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    CONC_OXA_UM = list(
      description = paste(
        "Static oxaliplatin concentration applied to the culture medium of the",
        "96 h tumour-organoid drug-sensitivity assay."
      ),
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-invariant over the 96 h incubation (Zhu 2023 Section 2.1). Set to",
        "0 for the vehicle-control well. Drives the Hill killing term in",
        "Equation (2). Zhu 2023 does not tabulate the tested concentration",
        "grid; Figure 2A spans roughly 1-1000 umol/L."
      ),
      source_name = "C"
    )
  )

  population <- list(
    species = "in vitro (patient-derived colorectal-cancer tumour organoids)",
    n_subjects = 5L,
    n_studies = 1L,
    disease_state = "colorectal cancer",
    dose_range = "static medium concentrations over a 96 h incubation",
    notes = paste(
      "Organoid drug-sensitivity data were provided by Accurate International",
      "Biotechnology (Guangzhou, China); the assay is described in Zhu 2023",
      "Section 2.1 (96 h incubation, CellTiter-Glo 3D luminescence normalised to",
      "vehicle controls). Zhu 2023 Table 2 reports patient-specific parameters",
      "for 5 organoid lines for oxaliplatin (patients 1-5; patients 6 and 7 have",
      "no oxaliplatin data). No demographic characteristics of the donor",
      "patients are reported."
    )
  )

  ini({
    # Growth rate of the organoids; held constant by the authors ("The growth
    # rate of PDTO was fixed at 0.03, according to the natural growth of the
    # PDTOs", Zhu 2023 Section 3.1) and reported as 0.03 for every patient in
    # Table 2.
    lkgrow <- fixed(log(0.03)); label("Natural organoid growth rate kg (1/h)")            # Zhu 2023 Table 2, kg column (0.03 for patients 1-5)

    # Drug-effect parameters. Zhu 2023 Table 2 reports only patient-specific
    # estimates (no typical value); Table 4 gives the geometric means of those
    # same in vitro estimates, which the authors used as the median of the
    # virtual-trial sampling distribution. The Table 4 geometric means are used
    # here as the reference typical values.
    lemax <- log(0.073);      label("Maximum organoid killing rate Emax (1/h)")           # Zhu 2023 Table 4, Emax_OXA = 0.073 (geometric mean of Table 2 patients 1-5)
    lec50 <- log(358);        label("Oxaliplatin concentration at half-maximal killing EC50 (umol/L)")  # Zhu 2023 Table 4, EC50_OXA = 358 umol/L (geometric mean of Table 2)
    lhill <- log(0.61);       label("Hill coefficient of the killing term (unitless)")    # Zhu 2023 Table 4, hill_OXA = 0.61 (geometric mean of Table 2)

    # Zhu 2023 fitted patient-specific parameters and reports no random-effect
    # variances or residual-error estimates for the in vitro PD model, so the
    # residual SD is held at zero. See the vignette Errata.
    addSd_viability <- fixed(0); label("Additive residual SD on cell viability (fraction; not reported in the source)")  # Zhu 2023 reports no residual-error model for the in vitro PD fit
  })

  model({
    kgrow <- exp(lkgrow)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Applied medium concentration; max() guards the fractional power against a
    # negative covariate value.
    cdrug <- max(CONC_OXA_UM, 0)

    # Sigmoidal Emax killing rate (1/h), Zhu 2023 Equation (2) inner term.
    kkill <- emax * cdrug^hill / (ec50^hill + cdrug^hill)

    # Zhu 2023 Equation (1): vehicle-control organoid, exponential growth.
    d/dt(organoid_ctrl) <- kgrow * organoid_ctrl

    # Zhu 2023 Equation (2): treated organoid, exponential growth minus the
    # Hill killing term.
    d/dt(organoid_trt) <- kgrow * organoid_trt - kkill * organoid_trt

    # Both organoids start from the same (arbitrary) seeding volume, so
    # viability is scale-free.
    organoid_ctrl(0) <- 1
    organoid_trt(0) <- 1

    # Zhu 2023 Equation (3): cell viability normalised to vehicle control.
    viability <- organoid_trt / organoid_ctrl

    viability ~ add(addSd_viability)
  })
}
