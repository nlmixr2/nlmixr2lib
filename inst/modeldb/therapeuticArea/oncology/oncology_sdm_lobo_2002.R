oncology_sdm_lobo_2002 <- function() {
  description <- "In vitro (Ehrlich ascites cell line). Signal transduction (transit compartment) model for delayed concentration effects on cancer cell growth: a drug concentration Cc drives a saturable cell-kill rate K = kmax*Cc/(ec50 + Cc), which passes through four transit compartments of mean transit time tau before killing exponentially growing cells (Lobo 2002, methotrexate, Ehrlich ascites parameter set)."
  reference <- paste(
    "Lobo ED, Balthasar JP. Pharmacodynamic modeling of chemotherapeutic",
    "effects: application of a transit compartment model to characterize",
    "methotrexate effects in vitro. AAPS PharmSci. 2002;4(4):E42.",
    "doi:10.1208/ps040442.",
    sep = " "
  )
  vignette <- "oncology_sdm_lobo_2002"
  # Cc (the drug concentration, M in Lobo 2002) and tumor_vol0 (the cell
  # number at the start of exposure) are inputs supplied by the user or by an
  # upstream PK model; they are listed here rather than in covariateData.
  depends <- c("Cc", "tumor_vol0")
  covariateData <- list()
  # No dose records: drug exposure enters only through the Cc input.
  dosing <- character(0)
  units <- list(
    time = "h",
    dosing = "none (drug exposure is the Cc input, ug/mL)",
    concentration = "ug/mL (the Cc input); tumor_vol is a cell number (cells/well in Lobo 2002)"
  )

  compartmentData <- list(
    tumor_vol = list(
      analyte = "viable tumour cells (Ehrlich ascites cells in culture)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "drug-induced cell-kill rate signal K1",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "drug-induced cell-kill rate signal K2",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "drug-induced cell-kill rate signal K3",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "drug-induced cell-kill rate signal K4",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (Ehrlich ascites cell line)",
    n_subjects = NA_integer_,
    disease_state = "Ehrlich ascites carcinoma cells grown in 96-well plates (RPMI 1640, 10% fetal bovine serum, 37 C, 5% CO2)",
    dose_range = "Methotrexate 0.19, 2.0, 14.0 or 140.0 ug/mL in the medium for 24 h, then washed out; untreated control",
    notes = paste(
      "Lobo 2002 Methods (pages 2-3): 1000 cells seeded per well and allowed to",
      "attach for 48 h before methotrexate; 4 wells per concentration; viable",
      "cells counted by MTT assay on days 1-24 after exposure. The model was",
      "fitted in ADAPT II to the mean cell number of the 4 wells, all",
      "concentrations simultaneously. The paper fitted a second cell line",
      "(sarcoma 180, Table 3: kng 0.035 1/h, tau 30.0 h, Kmax 0.34 1/h,",
      "EC50 0.32 ug/mL) that is not packaged here; the article reproduces it",
      "through ini() overrides."
    )
  )

  ini({
    # Transit compartment model, Ehrlich ascites cells: Lobo 2002 Table 2 (page 9)
    # and Results text (page 7).
    lkng <- log(0.02) ; label("Cell net growth rate (growth minus death) (1/h)") # Table 2 transit compartment: kng = 0.02 1/h (%CV 4.9)
    ltau <- log(34.1) ; label("Mean transit time of each transit compartment (h)") # Table 2 transit compartment: tau = 34.1 h (%CV 3.4)
    lec50 <- log(0.1) ; label("Drug concentration giving half the maximal cell-kill rate (ug/mL)") # Table 2 transit compartment: EC50 = 0.1 ug/mL (%CV 48.7)
    kmax <- 0.29 ; label("Maximum drug-induced cell-kill rate (1/h)") # Table 2 transit compartment: Kmax = 0.29 1/h (%CV 5.8)

    # Not reported by Lobo 2002 (ADAPT II fit to mean cell counts; the
    # variance model is not given). Placeholder values for simulation.
    propSd_tumor_vol <- c(0, 0.3) ; label("Proportional residual error (fraction)")
    addSd_tumor_vol <- c(0, 50, 1000) ; label("Additive residual error (cells)")
  })
  model({
    # Lobo 2002 page 3, transit compartment model equations. Cc is the drug
    # concentration (M in the paper, methotrexate in the medium, ug/mL).
    kng <- exp(lkng)
    tau <- exp(ltau)
    ec50 <- exp(lec50)

    drugEffectTumorVol <- kmax*Cc/(ec50 + Cc)

    tumor_vol(0) <- tumor_vol0
    d/dt(tumor_vol) <- kng*tumor_vol - transit4*tumor_vol
    d/dt(transit1) <- (drugEffectTumorVol - transit1)/tau
    d/dt(transit2) <- (transit1 - transit2)/tau
    d/dt(transit3) <- (transit2 - transit3)/tau
    d/dt(transit4) <- (transit3 - transit4)/tau
    tumor_vol ~ prop(propSd_tumor_vol) + add(addSd_tumor_vol)
  })
}
