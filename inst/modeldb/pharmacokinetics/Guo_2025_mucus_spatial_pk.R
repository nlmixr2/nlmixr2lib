Guo_2025_mucus_spatial_pk <- function() {
  description <- paste(
    "Mechanistic spatial PK (method-of-lines discretised PDE). Diffusion of a",
    "hypothetical intravenously dosed drug across a 500 um layer of cystic",
    "fibrosis airway mucus (Guo 2025). A one-compartment plasma model sets the",
    "free concentration at the blood-mucus interface (instantaneous transfer,",
    "partition 1, no protein binding); within mucus the free drug diffuses by",
    "Fick's second law with a Stokes-Einstein diffusion coefficient computed",
    "from the molecule/particle radius and a mucin-dependent micro-viscosity,",
    "binds reversibly to immobile mucins (kon, koff) and is removed with",
    "muco-ciliary clearance (free and bound alike); the air-mucus interface is",
    "a zero-flux face. The layer is discretised on the authors' 25 um grid",
    "(21 nodes; 19 free and 20 bound drug states). Defaults are the paper's",
    "Table 1 small molecule (r = 1 nm, D = 78.6 um^2/s, t1/2 3 h); the radius,",
    "half-life, mucin, binding and muco-ciliary knobs reproduce the Figure 3-5",
    "sweeps. No inter-individual or residual variability (simulation framework)."
  )
  reference <- paste(
    "Guo Y, Yin J, Yonucu S, Knibbe CAJ, Guo T, van Hasselt JGC. Spatial",
    "Pharmacokinetic and Pharmacodynamic Modeling in Airway Mucus.",
    "Clin Pharmacokinet. 2025;64(12):1869-1883. doi:10.1007/s40262-025-01575-4.",
    "Equations 1-5 and Table 1 of the article, Equation S1 and Table S1 of the",
    "supplement; plasma volume, infusion duration and molecular weight from the",
    "authors' deposited code (https://github.com/vanhasseltlab/spatialPK).",
    sep = " "
  )
  vignette <- "Guo_2025_mucus_spatial_pkpd"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "hypothetical drug", units = "mg", specimen = "plasma", verified = TRUE),
    mucus_free_slab2 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab3 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab4 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab5 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab6 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab7 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab8 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab9 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab10 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab11 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab12 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab13 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab14 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab15 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab16 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab17 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab18 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab19 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab20 = list(
      analyte = "hypothetical drug (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab1 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab2 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab3 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab4 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab5 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab6 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab7 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab8 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab9 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab10 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab11 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab12 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab13 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab14 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab15 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab16 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab17 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab18 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab19 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab20 = list(
      analyte = "hypothetical drug (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "Not applicable (typical adult)",
    weight_range = "Not applicable",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not applicable",
    disease_state = "Cystic fibrosis airway mucus (physiological system parameters from CF sputum and mucus-plug literature)",
    dose_range = "100 mg IV (0.5 h infusion) every 8 h for 72 h; 300 mg/day as q8h, q24h or continuous infusion (Figure 5); 600 mg every 2 weeks for large molecules (Figure 4B)",
    regions = "Not applicable",
    notes = paste(
      "No patients were studied. Hypothetical drugs spanning molecule/particle",
      "radius 1, 5, 20 and 50 nm, plasma half-life 3-24 h (small) or 1-4 weeks",
      "(large), mucin association constant 10-10000 1/M, total mucin 1e-6 to",
      "5e-5 mol/L and muco-ciliary clearance 0-0.36 1/h were simulated against",
      "CF physiological mucus parameters (Table 1, Table S1)."
    )
  )

  ini({
    # Plasma PK (Eq. 3). The paper gives the half-life (3 h, Fig. 4 reference)
    # but not the volume; V = 15 L and CL = 0.05775 L/min come from the authors'
    # deposited code (scripts/01 and 02), and reproduce the 3 h half-life.
    lcl <- fixed(log(3.465)); label("Plasma clearance (L/h)") # deposited code CL = 0.05775 L/min x 60; t1/2 = ln(2) x 15 / 3.465 = 3.0 h (Fig. 4 caption reference t1/2 = 3 h)
    lvc <- fixed(log(15)); label("Plasma volume of distribution (L)") # deposited code V1 = 15 L (not printed in the article)

    # Diffusion coefficient (Eq. 2, Eq. S1, Table S1)
    r_hyd <- fixed(1); label("Hydrodynamic radius of the molecule/particle (nm)") # Table 1 title 'small molecule (r = 1 nm)'; Table S1 r = 1, 5, 20, 50 nm
    kboltz <- fixed(1.380649e-23); label("Boltzmann constant (N*m/K)") # Table S1 kB = 1.380649e-23
    temp_k <- fixed(310); label("Absolute temperature (K)") # Table S1 T = 310
    visc_slope <- fixed(128.87); label("Slope of mucus micro-viscosity on total mucin ((N*s/m^2)/(mol/L))") # Eq. S1 mu = 128.87 * Mucin + 0.00176
    visc_int <- fixed(0.00176); label("Intercept of mucus micro-viscosity (N*s/m^2)") # Eq. S1 mu = 128.87 * Mucin + 0.00176

    # Mucus system and mucin binding (Eqs. 4-5, Table 1)
    thick_mucus <- fixed(500); label("Thickness of the mucus layer (um)") # Table 1 L = 500 um
    mucin <- fixed(2e-5); label("Total mucin concentration in mucus (mol/L)") # Table 1 M = 2E-5 mol/L
    kon <- fixed(100); label("Drug-mucin association rate constant (1/(M*s))") # Table 1 kon = 100 (M*sec)^-1
    koff <- fixed(0.1); label("Drug-mucin dissociation rate constant (1/s)") # Table 1 koff = 0.1 sec^-1
    kel_mcc <- fixed(0.036); label("Muco-ciliary elimination rate constant (1/h)") # Table 1 ke,m = 0.04 h^-1, printed rounded; Fig. 4 caption reference Ke,m = 0.036 h^-1 = Table S1 1e-05 sec^-1
    mw <- fixed(500); label("Molecular weight used to express bound drug in mol/L (g/mol)") # deposited code pde() MW = 500 (not printed in the article)
  })

  model({
    # Plasma PK (Eq. 3; dosed as an IV infusion into central)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Boundary condition at the blood-mucus interface (x = 0): free mucus
    # concentration equals the plasma concentration (partition 1, fu = 1)
    Cmucus_bc <- Cc

    # Diffusion coefficient: modified Stokes-Einstein equation (Eq. 2) with the
    # mucin-dependent micro-viscosity of Eq. S1; c = 2/3 for r <= 1 nm, else 1.
    # The factor 1e21 converts kB*T/(nm * N*s/m^2) from m^2/s to um^2/s.
    visc <- visc_slope * mucin + visc_int
    c_se <- 1
    if (r_hyd <= 1) {
      c_se <- 2 / 3
    }
    diff_coef <- 1e21 * kboltz * temp_k / (6 * c_se * pi * visc * r_hyd)

    # Method of lines on the authors' 25 um grid: 20 intervals, nodes at
    # x = (n - 1) * dx for n = 1..21. Rates converted to 1/h.
    dx <- thick_mucus / 20
    kdiff <- diff_coef * 3600 / dx^2
    konh <- kon * 3600
    koffh <- koff * 3600

    # Free drug, nodes 2-20 (node 1 = boundary, node 21 = node 20)
    d/dt(mucus_free_slab2) <- kdiff * (Cmucus_bc - 2 * mucus_free_slab2 + mucus_free_slab3) - konh * (mucin - mucus_bound_slab2 / (1000 * mw)) * mucus_free_slab2 + koffh * mucus_bound_slab2 - kel_mcc * mucus_free_slab2
    d/dt(mucus_free_slab3) <- kdiff * (mucus_free_slab2 - 2 * mucus_free_slab3 + mucus_free_slab4) - konh * (mucin - mucus_bound_slab3 / (1000 * mw)) * mucus_free_slab3 + koffh * mucus_bound_slab3 - kel_mcc * mucus_free_slab3
    d/dt(mucus_free_slab4) <- kdiff * (mucus_free_slab3 - 2 * mucus_free_slab4 + mucus_free_slab5) - konh * (mucin - mucus_bound_slab4 / (1000 * mw)) * mucus_free_slab4 + koffh * mucus_bound_slab4 - kel_mcc * mucus_free_slab4
    d/dt(mucus_free_slab5) <- kdiff * (mucus_free_slab4 - 2 * mucus_free_slab5 + mucus_free_slab6) - konh * (mucin - mucus_bound_slab5 / (1000 * mw)) * mucus_free_slab5 + koffh * mucus_bound_slab5 - kel_mcc * mucus_free_slab5
    d/dt(mucus_free_slab6) <- kdiff * (mucus_free_slab5 - 2 * mucus_free_slab6 + mucus_free_slab7) - konh * (mucin - mucus_bound_slab6 / (1000 * mw)) * mucus_free_slab6 + koffh * mucus_bound_slab6 - kel_mcc * mucus_free_slab6
    d/dt(mucus_free_slab7) <- kdiff * (mucus_free_slab6 - 2 * mucus_free_slab7 + mucus_free_slab8) - konh * (mucin - mucus_bound_slab7 / (1000 * mw)) * mucus_free_slab7 + koffh * mucus_bound_slab7 - kel_mcc * mucus_free_slab7
    d/dt(mucus_free_slab8) <- kdiff * (mucus_free_slab7 - 2 * mucus_free_slab8 + mucus_free_slab9) - konh * (mucin - mucus_bound_slab8 / (1000 * mw)) * mucus_free_slab8 + koffh * mucus_bound_slab8 - kel_mcc * mucus_free_slab8
    d/dt(mucus_free_slab9) <- kdiff * (mucus_free_slab8 - 2 * mucus_free_slab9 + mucus_free_slab10) - konh * (mucin - mucus_bound_slab9 / (1000 * mw)) * mucus_free_slab9 + koffh * mucus_bound_slab9 - kel_mcc * mucus_free_slab9
    d/dt(mucus_free_slab10) <- kdiff * (mucus_free_slab9 - 2 * mucus_free_slab10 + mucus_free_slab11) - konh * (mucin - mucus_bound_slab10 / (1000 * mw)) * mucus_free_slab10 + koffh * mucus_bound_slab10 - kel_mcc * mucus_free_slab10
    d/dt(mucus_free_slab11) <- kdiff * (mucus_free_slab10 - 2 * mucus_free_slab11 + mucus_free_slab12) - konh * (mucin - mucus_bound_slab11 / (1000 * mw)) * mucus_free_slab11 + koffh * mucus_bound_slab11 - kel_mcc * mucus_free_slab11
    d/dt(mucus_free_slab12) <- kdiff * (mucus_free_slab11 - 2 * mucus_free_slab12 + mucus_free_slab13) - konh * (mucin - mucus_bound_slab12 / (1000 * mw)) * mucus_free_slab12 + koffh * mucus_bound_slab12 - kel_mcc * mucus_free_slab12
    d/dt(mucus_free_slab13) <- kdiff * (mucus_free_slab12 - 2 * mucus_free_slab13 + mucus_free_slab14) - konh * (mucin - mucus_bound_slab13 / (1000 * mw)) * mucus_free_slab13 + koffh * mucus_bound_slab13 - kel_mcc * mucus_free_slab13
    d/dt(mucus_free_slab14) <- kdiff * (mucus_free_slab13 - 2 * mucus_free_slab14 + mucus_free_slab15) - konh * (mucin - mucus_bound_slab14 / (1000 * mw)) * mucus_free_slab14 + koffh * mucus_bound_slab14 - kel_mcc * mucus_free_slab14
    d/dt(mucus_free_slab15) <- kdiff * (mucus_free_slab14 - 2 * mucus_free_slab15 + mucus_free_slab16) - konh * (mucin - mucus_bound_slab15 / (1000 * mw)) * mucus_free_slab15 + koffh * mucus_bound_slab15 - kel_mcc * mucus_free_slab15
    d/dt(mucus_free_slab16) <- kdiff * (mucus_free_slab15 - 2 * mucus_free_slab16 + mucus_free_slab17) - konh * (mucin - mucus_bound_slab16 / (1000 * mw)) * mucus_free_slab16 + koffh * mucus_bound_slab16 - kel_mcc * mucus_free_slab16
    d/dt(mucus_free_slab17) <- kdiff * (mucus_free_slab16 - 2 * mucus_free_slab17 + mucus_free_slab18) - konh * (mucin - mucus_bound_slab17 / (1000 * mw)) * mucus_free_slab17 + koffh * mucus_bound_slab17 - kel_mcc * mucus_free_slab17
    d/dt(mucus_free_slab18) <- kdiff * (mucus_free_slab17 - 2 * mucus_free_slab18 + mucus_free_slab19) - konh * (mucin - mucus_bound_slab18 / (1000 * mw)) * mucus_free_slab18 + koffh * mucus_bound_slab18 - kel_mcc * mucus_free_slab18
    d/dt(mucus_free_slab19) <- kdiff * (mucus_free_slab18 - 2 * mucus_free_slab19 + mucus_free_slab20) - konh * (mucin - mucus_bound_slab19 / (1000 * mw)) * mucus_free_slab19 + koffh * mucus_bound_slab19 - kel_mcc * mucus_free_slab19
    d/dt(mucus_free_slab20) <- kdiff * (mucus_free_slab19 - mucus_free_slab20) - konh * (mucin - mucus_bound_slab20 / (1000 * mw)) * mucus_free_slab20 + koffh * mucus_bound_slab20 - kel_mcc * mucus_free_slab20
    # Mucin-bound drug, nodes 1-20 (immobile)
    d/dt(mucus_bound_slab1) <- konh * (mucin - mucus_bound_slab1 / (1000 * mw)) * Cmucus_bc - koffh * mucus_bound_slab1 - kel_mcc * mucus_bound_slab1
    d/dt(mucus_bound_slab2) <- konh * (mucin - mucus_bound_slab2 / (1000 * mw)) * mucus_free_slab2 - koffh * mucus_bound_slab2 - kel_mcc * mucus_bound_slab2
    d/dt(mucus_bound_slab3) <- konh * (mucin - mucus_bound_slab3 / (1000 * mw)) * mucus_free_slab3 - koffh * mucus_bound_slab3 - kel_mcc * mucus_bound_slab3
    d/dt(mucus_bound_slab4) <- konh * (mucin - mucus_bound_slab4 / (1000 * mw)) * mucus_free_slab4 - koffh * mucus_bound_slab4 - kel_mcc * mucus_bound_slab4
    d/dt(mucus_bound_slab5) <- konh * (mucin - mucus_bound_slab5 / (1000 * mw)) * mucus_free_slab5 - koffh * mucus_bound_slab5 - kel_mcc * mucus_bound_slab5
    d/dt(mucus_bound_slab6) <- konh * (mucin - mucus_bound_slab6 / (1000 * mw)) * mucus_free_slab6 - koffh * mucus_bound_slab6 - kel_mcc * mucus_bound_slab6
    d/dt(mucus_bound_slab7) <- konh * (mucin - mucus_bound_slab7 / (1000 * mw)) * mucus_free_slab7 - koffh * mucus_bound_slab7 - kel_mcc * mucus_bound_slab7
    d/dt(mucus_bound_slab8) <- konh * (mucin - mucus_bound_slab8 / (1000 * mw)) * mucus_free_slab8 - koffh * mucus_bound_slab8 - kel_mcc * mucus_bound_slab8
    d/dt(mucus_bound_slab9) <- konh * (mucin - mucus_bound_slab9 / (1000 * mw)) * mucus_free_slab9 - koffh * mucus_bound_slab9 - kel_mcc * mucus_bound_slab9
    d/dt(mucus_bound_slab10) <- konh * (mucin - mucus_bound_slab10 / (1000 * mw)) * mucus_free_slab10 - koffh * mucus_bound_slab10 - kel_mcc * mucus_bound_slab10
    d/dt(mucus_bound_slab11) <- konh * (mucin - mucus_bound_slab11 / (1000 * mw)) * mucus_free_slab11 - koffh * mucus_bound_slab11 - kel_mcc * mucus_bound_slab11
    d/dt(mucus_bound_slab12) <- konh * (mucin - mucus_bound_slab12 / (1000 * mw)) * mucus_free_slab12 - koffh * mucus_bound_slab12 - kel_mcc * mucus_bound_slab12
    d/dt(mucus_bound_slab13) <- konh * (mucin - mucus_bound_slab13 / (1000 * mw)) * mucus_free_slab13 - koffh * mucus_bound_slab13 - kel_mcc * mucus_bound_slab13
    d/dt(mucus_bound_slab14) <- konh * (mucin - mucus_bound_slab14 / (1000 * mw)) * mucus_free_slab14 - koffh * mucus_bound_slab14 - kel_mcc * mucus_bound_slab14
    d/dt(mucus_bound_slab15) <- konh * (mucin - mucus_bound_slab15 / (1000 * mw)) * mucus_free_slab15 - koffh * mucus_bound_slab15 - kel_mcc * mucus_bound_slab15
    d/dt(mucus_bound_slab16) <- konh * (mucin - mucus_bound_slab16 / (1000 * mw)) * mucus_free_slab16 - koffh * mucus_bound_slab16 - kel_mcc * mucus_bound_slab16
    d/dt(mucus_bound_slab17) <- konh * (mucin - mucus_bound_slab17 / (1000 * mw)) * mucus_free_slab17 - koffh * mucus_bound_slab17 - kel_mcc * mucus_bound_slab17
    d/dt(mucus_bound_slab18) <- konh * (mucin - mucus_bound_slab18 / (1000 * mw)) * mucus_free_slab18 - koffh * mucus_bound_slab18 - kel_mcc * mucus_bound_slab18
    d/dt(mucus_bound_slab19) <- konh * (mucin - mucus_bound_slab19 / (1000 * mw)) * mucus_free_slab19 - koffh * mucus_bound_slab19 - kel_mcc * mucus_bound_slab19
    d/dt(mucus_bound_slab20) <- konh * (mucin - mucus_bound_slab20 / (1000 * mw)) * mucus_free_slab20 - koffh * mucus_bound_slab20 - kel_mcc * mucus_bound_slab20
    
    # Spatial summaries over the 21 grid nodes (node 21 = node 20)
    Cmucus_free_mean <- (Cmucus_bc + mucus_free_slab2 + mucus_free_slab3 + mucus_free_slab4 + mucus_free_slab5 + mucus_free_slab6 + mucus_free_slab7 + mucus_free_slab8 + mucus_free_slab9 + mucus_free_slab10 + mucus_free_slab11 + mucus_free_slab12 + mucus_free_slab13 + mucus_free_slab14 + mucus_free_slab15 + mucus_free_slab16 + mucus_free_slab17 + mucus_free_slab18 + mucus_free_slab19 + mucus_free_slab20 + mucus_free_slab20) / 21
    Cmucus_bound_mean <- (mucus_bound_slab1 + mucus_bound_slab2 + mucus_bound_slab3 + mucus_bound_slab4 + mucus_bound_slab5 + mucus_bound_slab6 + mucus_bound_slab7 + mucus_bound_slab8 + mucus_bound_slab9 + mucus_bound_slab10 + mucus_bound_slab11 + mucus_bound_slab12 + mucus_bound_slab13 + mucus_bound_slab14 + mucus_bound_slab15 + mucus_bound_slab16 + mucus_bound_slab17 + mucus_bound_slab18 + mucus_bound_slab19 + mucus_bound_slab20 + mucus_bound_slab20) / 21
    Cmucus_free_cv <- 100 * sqrt(((Cmucus_bc - Cmucus_free_mean)^2 + (mucus_free_slab2 - Cmucus_free_mean)^2 + (mucus_free_slab3 - Cmucus_free_mean)^2 + (mucus_free_slab4 - Cmucus_free_mean)^2 + (mucus_free_slab5 - Cmucus_free_mean)^2 + (mucus_free_slab6 - Cmucus_free_mean)^2 + (mucus_free_slab7 - Cmucus_free_mean)^2 + (mucus_free_slab8 - Cmucus_free_mean)^2 + (mucus_free_slab9 - Cmucus_free_mean)^2 + (mucus_free_slab10 - Cmucus_free_mean)^2 + (mucus_free_slab11 - Cmucus_free_mean)^2 + (mucus_free_slab12 - Cmucus_free_mean)^2 + (mucus_free_slab13 - Cmucus_free_mean)^2 + (mucus_free_slab14 - Cmucus_free_mean)^2 + (mucus_free_slab15 - Cmucus_free_mean)^2 + (mucus_free_slab16 - Cmucus_free_mean)^2 + (mucus_free_slab17 - Cmucus_free_mean)^2 + (mucus_free_slab18 - Cmucus_free_mean)^2 + (mucus_free_slab19 - Cmucus_free_mean)^2 + (mucus_free_slab20 - Cmucus_free_mean)^2 + (mucus_free_slab20 - Cmucus_free_mean)^2) / 20) / Cmucus_free_mean
  })
}
