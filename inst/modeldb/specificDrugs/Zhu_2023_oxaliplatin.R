Zhu_2023_oxaliplatin <- function() {
  description <- paste(
    "Two-compartment population PK model for ultrafilterable platinum in plasma",
    "after intravenous oxaliplatin in adults with metastatic colorectal cancer,",
    "estimated from human concentration-time profiles digitised from the",
    "literature. Zhu 2023 uses the central-compartment concentration as the",
    "driver of an in vitro-derived tumour-killing model; that in vivo",
    "tumour-growth layer is NOT included here because its growth equation and",
    "limiting parameter are not reported anywhere in the paper - see the",
    "vignette Errata. The in vitro drug-effect half is available as",
    "modellib('Zhu_2023_oxaliplatin_organoid')."
  )
  reference <- paste(
    "Zhu J, Zhang Y, Zhao Y, Zhang J, Hao K, He H. Translational",
    "Pharmacokinetic/Pharmacodynamic Modeling and Simulation of Oxaliplatin and",
    "Irinotecan in Colorectal Cancer. Pharmaceutics. 2023;15(9):2274.",
    "doi:10.3390/pharmaceutics15092274.",
    sep = " "
  )
  vignette <- "Zhu_2023_oxaliplatin_irinotecan_colorectal_cancer"

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    central = list(
      analyte = "ultrafilterable platinum", units = "umol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ultrafilterable platinum", units = "umol",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "metastatic colorectal cancer",
    dose_range = "130 mg/m2 IV every 3 weeks or 85 mg/m2 IV every 2 weeks (simulated regimens, Zhu 2023 Table 1)",
    notes = paste(
      "The plasma PK parameters were estimated from published human",
      "ultrafilterable-platinum concentration-time profiles digitised from the",
      "literature with GetData Graph Digitizer (Zhu 2023 Section 2.1); the",
      "individual source publications, subject counts and demographics are not",
      "reported. Studies in hepatic or renal impairment, studies with",
      "co-administered antitumour drugs, and whole-blood-only studies were",
      "excluded. The reported interindividual variability therefore mixes",
      "between-subject with between-study variability."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK of ultrafilterable platinum, Zhu 2023 Equations (4)-(5).
    # All values from Zhu 2023 Table 3, "Human Estimates (RSE%)" column.
    # ------------------------------------------------------------------
    lvc <- log(49.9);   label("Central volume of ultrafilterable platinum VC_OXA (L)")                # Zhu 2023 Table 3, VC_OXA human = 49.9 L (RSE 46.3%)
    lvp <- log(538);    label("Peripheral volume of ultrafilterable platinum VP_OXA (L)")             # Zhu 2023 Table 3, VP_OXA human = 538 L (RSE 29.3%)
    lcl <- log(5.96);   label("Systemic clearance of ultrafilterable platinum CL_OXA (L/h)")          # Zhu 2023 Table 3, CL_OXA human = 5.96 L/h (RSE 42.5%)
    lq  <- log(49.3);   label("Intercompartmental clearance Q_OXA (L/h)")                             # Zhu 2023 Table 3, Q_OXA human = 49.3 L/h (RSE 29.7%)

    # ------------------------------------------------------------------
    # Interindividual variability. Monolix reports omega as the standard
    # deviation of the log-scale random effect, so the variance entered
    # here is the square of the tabulated IIV. VP_OXA and Q_OXA carry no
    # IIV in Zhu 2023 Table 3.
    # ------------------------------------------------------------------
    etalvc ~ 1.1025;    # Zhu 2023 Table 3, VC_OXA IIV = 1.05 (RSE 30.9%); variance = 1.05^2
    etalcl ~ 0.356409;  # Zhu 2023 Table 3, CL_OXA IIV = 0.597 (RSE 63.6%); variance = 0.597^2

    # Zhu 2023 reports no residual-error estimates, so the residual SD is held
    # at zero. See the vignette Errata.
    propSd <- fixed(0); label("Proportional residual error on plasma ultrafilterable platinum (fraction; not reported in the source)")  # Zhu 2023 reports no residual-error model
  })

  model({
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    cl <- exp(lcl + etalcl)
    q  <- exp(lq)

    # States are amounts (umol); the paper writes its equations as
    # V * dC/dt, so the amount derivative is exactly the published RHS.
    Cc <- central / vc
    Cp <- peripheral1 / vp

    # Zhu 2023 Equation (4)
    d/dt(central) <- -(cl + q) * Cc + q * Cp
    # Zhu 2023 Equation (5)
    d/dt(peripheral1) <- q * (Cc - Cp)

    Cc ~ prop(propSd)
  })
}
