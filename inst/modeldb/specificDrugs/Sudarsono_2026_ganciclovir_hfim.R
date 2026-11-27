Sudarsono_2026_ganciclovir_hfim <- function() {
  description <- "In vitro (hollow-fiber infection model; MRC-5 human lung fibroblasts infected with cytomegalovirus). One-compartment model of ganciclovir in the hollow-fiber system (central reservoir plus hemodialyzer cartridge), fitted to concentrations measured in the intracapillary space during 1-h infusions every 12 h. Elimination is the dilution of the reservoir by drug-free medium, tuned to mimic a human-like plasma profile."
  reference <- "Sudarsono LM, Wenker SAM, Liu X, Brink J, van den Berg DJ, van Hasselt JGC, Martson AG. A novel hollow-fiber infection model (HFIM) for antiviral PK/PD studies of CMV infection. Antimicrob Agents Chemother. 2026;70(8):e00381-26. doi:10.1128/aac.00381-26. Fitted estimates: Table 1 'Estimates (%RSE)' column; design (input) values: Table 1 'Input' column and Methods Eqs 1-3."
  vignette <- "Sudarsono_2026_ganciclovir_hfim"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "ganciclovir",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (hollow-fiber infection model; MRC-5 cells, cytomegalovirus Towne and Merlin strains)",
    n_subjects = 2L,
    n_studies = 1L,
    disease_state = "Human cytomegalovirus infection of MRC-5 human lung fibroblasts (5.0e7 cells inoculated into the extracapillary space; MOI 0.1 for the lab-adapted Towne strain and 0.01 for the BAC-derived wild-type-like Merlin strain).",
    model_system = "Hollow-fiber infection model built from a low-cost FX paed hemodialyzer cartridge (polysulfone fibres, 3.3 nm pores) connected to a central reservoir held at 0.2 L; main circulation loop at 30 mL/min; drug-free medium pumped into the reservoir with excess volume to waste to mimic clearance; 37 degrees C, 5% CO2.",
    dose_range = "Ganciclovir 2 mg infused over 1 h (2 mg/h) into the central reservoir every 12 h for 3 days.",
    notes = paste(
      "n_subjects counts the two independent hollow-fiber experiments whose intracapillary-space concentrations were pooled for the fit (Figure 3).",
      "Pump settings were designed for a 10 mg/L peak and a 4-h half-life with the reservoir fixed at 0.2 L (Methods Eqs 1-3; Table 1 'Input' CL 0.03 L/h, V 0.2 L, Ke 0.15 1/h).",
      "The one-compartment fit (nlmixr2) estimated CL 0.033 L/h and V 0.3 L; the larger fitted volume presumably reflects the cartridge and tubing in addition to the reservoir.",
      "Concentrations in the cartridge and the central reservoir did not differ (Figure 4), supporting a single well-mixed compartment with no binding to the system.",
      "No between-experiment variability or residual-error terms were reported, so the model is encoded as typical values only.",
      "Viral replication under ganciclovir (Figure 5) was reported descriptively; no pharmacodynamic model was fitted."
    )
  )

  ini({
    lcl <- log(0.033)
    label("Clearance of the hollow-fiber system (L/h)") # Table 1 'Clearance (L/h)' estimate 0.033 (0.2 %RSE)
    lvc <- log(0.3)
    label("Volume of distribution of the hollow-fiber system (L)") # Table 1 'Volume of distribution (L)' estimate 0.3 (2.05 %RSE)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc
    d/dt(central) <- -kel * central
    Cc <- central / vc
  })
}
