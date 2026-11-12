Guo_2025_imipenem_mucus <- function() {
  description <- paste(
    "Mechanistic spatial PK/PD (method-of-lines discretised PDE). Imipenem",
    "killing of Pseudomonas aeruginosa inside a 500 um layer of cystic",
    "fibrosis airway mucus (Guo 2025 case study). Plasma PK is the",
    "two-compartment imipenem model of van Hasselt 2016 (CL 11.52 L/h, V1",
    "9.37 L, V2 6.41 L, Q 13.68 L/h); the free concentration at the",
    "blood-mucus interface is the plasma concentration times the ELF",
    "penetration coefficient 0.444. Within mucus imipenem diffuses (Fick's",
    "second law, Stokes-Einstein D for r = 1 nm, about 78.6 um^2/s) and binds",
    "reversibly to immobile mucins; muco-ciliary clearance is set to zero.",
    "Immobile bacteria at every grid node follow a two-subpopulation",
    "(susceptible S, resistant R) logistic-growth model with linear",
    "concentration-driven kill and drug-independent S-to-R transfer.",
    "The layer is discretised on the authors' 25 um grid (21 nodes; 19 free",
    "drug, 20 bound drug, 20 S and 20 R states). The baseline inoculum is",
    "homogeneous at 10^5.67 CFU/mL; per-node log10 offsets",
    "(dlog10b0_slab<n>) give gradient baselines. Deterministic: no",
    "inter-individual or residual variability."
  )
  reference <- paste(
    "Guo Y, Yin J, Yonucu S, Knibbe CAJ, Guo T, van Hasselt JGC. Spatial",
    "Pharmacokinetic and Pharmacodynamic Modeling in Airway Mucus.",
    "Clin Pharmacokinet. 2025;64(12):1869-1883. doi:10.1007/s40262-025-01575-4.",
    "Equations 1-5, Equations S1-S3 and Table S2 of the supplement. Plasma PK",
    "from van Hasselt JGC, Rizk ML, Lala M, et al. Br J Clin Pharmacol.",
    "2016;81:1113-1123 (doi:10.1111/bcp.12901); PD from Guo Y, Yin J, Aulin LBS,",
    "et al. bioRxiv 2024 (doi:10.1101/2024.05.28.596312). Molecular weight from",
    "the authors' deposited code (https://github.com/vanhasseltlab/spatialPK).",
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
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    mucus_free_slab2 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab3 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab4 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab5 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab6 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab7 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab8 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab9 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab10 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab11 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab12 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab13 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab14 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab15 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab16 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab17 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab18 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab19 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_free_slab20 = list(
      analyte = "imipenem (free)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab1 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab2 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab3 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab4 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab5 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab6 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab7 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab8 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab9 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab10 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab11 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab12 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab13 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab14 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab15 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab16 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab17 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab18 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab19 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    mucus_bound_slab20 = list(
      analyte = "imipenem (mucin-bound)",
      units = "mg/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab1 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab1 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab2 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab2 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab3 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab3 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab4 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab4 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab5 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab5 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab6 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab6 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab7 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab7 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab8 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab8 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab9 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab9 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab10 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab10 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab11 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab11 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab12 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab12 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab13 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab13 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab14 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab14 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab15 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab15 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab16 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab16 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab17 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab17 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab18 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab18 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab19 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab19 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_susceptible_slab20 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-susceptible subpopulation (CFU/mL)",
      units = "CFU/mL",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    bact_resistant_slab20 = list(
      analyte = "Pseudomonas aeruginosa, imipenem-resistant subpopulation (CFU/mL)",
      units = "CFU/mL",
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
    disease_state = "Cystic fibrosis airway mucus infected with Pseudomonas aeruginosa (simulated case study)",
    dose_range = "Imipenem 250, 500 or 1000 mg IV per administration, as 0.5 h intermittent infusions or as a continuous infusion, for 72 h",
    regions = "Not applicable",
    notes = paste(
      "No patients were studied. A typical-value plasma PK model (van Hasselt",
      "2016) drives the mucus boundary condition; bacterial PD parameters come",
      "from in vitro P. aeruginosa time-kill data (Guo 2024 bioRxiv). Three",
      "baseline bacterial spatial organisations (homogeneous, increasing and",
      "decreasing gradient) were simulated (Fig. S1, Table S3)."
    )
  )

  ini({
    # Plasma PK of imipenem (Table S2, from van Hasselt 2016); converted from L/min to L/h
    lcl <- fixed(log(11.52)); label("Imipenem plasma clearance (L/h)") # Table S2 CL = 0.192 L/min x 60
    lvc <- fixed(log(9.37)); label("Central volume of distribution (L)") # Table S2 V1 = 9.37 L
    lvp <- fixed(log(6.41)); label("Peripheral volume of distribution (L)") # Table S2 V2 = 6.41 L
    lq <- fixed(log(13.68)); label("Inter-compartmental clearance (L/h)") # Table S2 Q = 0.228 L/min x 60
    relf <- fixed(0.444); label("ELF:plasma penetration coefficient setting the mucus boundary concentration (unitless)") # Table S2 P_ELF = 0.444

    # Diffusion coefficient (Eq. 2, Eq. S1, Table S1); imipenem treated as a small molecule
    r_hyd <- fixed(1); label("Hydrodynamic radius of imipenem (nm)") # deposited code scripts/03 Diffusion_coeff(1, vis(2e-5)); gives D = 78.6 um^2/s (Table S2 D = 79 um^2/s, rounded)
    kboltz <- fixed(1.380649e-23); label("Boltzmann constant (N*m/K)") # Table S1 kB = 1.380649e-23
    temp_k <- fixed(310); label("Absolute temperature (K)") # Table S1 T = 310
    visc_slope <- fixed(128.87); label("Slope of mucus micro-viscosity on total mucin ((N*s/m^2)/(mol/L))") # Eq. S1 mu = 128.87 * Mucin + 0.00176
    visc_int <- fixed(0.00176); label("Intercept of mucus micro-viscosity (N*s/m^2)") # Eq. S1 mu = 128.87 * Mucin + 0.00176

    # Mucus system and mucin binding (Eqs. 4-5, Table S2)
    thick_mucus <- fixed(500); label("Thickness of the mucus layer (um)") # Table S2 L = 500 um
    mucin <- fixed(2e-5); label("Total mucin concentration in mucus (mol/L)") # Table S2 M = 2e-5 mol/L
    kon <- fixed(100); label("Drug-mucin association rate constant (1/(M*s))") # Table S2 kon = 100 (M*sec)^-1
    koff <- fixed(0.1); label("Drug-mucin dissociation rate constant (1/s)") # Table S2 koff = 0.1 sec^-1
    kel_mcc <- fixed(0); label("Muco-ciliary elimination rate constant (1/h)") # Table S2 ke = 0 hour^-1
    mw <- fixed(500); label("Molecular weight used to express bound drug in mol/L (g/mol)") # deposited code pde_pkpd_emax() uses 500 g/mol (not printed in the article)

    # Bacterial PD (Eqs. S2-S3, Table S2, from Guo 2024 bioRxiv); rates per second converted in model()
    kgs <- fixed(1.67e-4); label("Natural growth rate of the susceptible subpopulation (1/s)") # Table S2 kgs = 1.67e-4 sec^-1
    kgr <- fixed(0); label("Natural growth rate of the resistant subpopulation (1/s)") # Table S2 kgr = 0 sec^-1
    ksr <- fixed(9.86e-6); label("Transfer rate from susceptible to resistant bacteria (1/s)") # Table S2 ksr = 9.86e-6 sec^-1
    slope_s <- fixed(2.38e-4); label("Linear kill coefficient on the susceptible subpopulation (L/(mg*s))") # Table S2 Slope_s = 2.38e-4 L/mg*sec
    slope_r <- fixed(3.78e-5); label("Linear kill coefficient on the resistant subpopulation (L/(mg*s))") # Table S2 Slope_r = 3.78e-5 L/mg*sec
    log10b0 <- fixed(5.67); label("Baseline total bacterial density, all susceptible (log10 CFU/mL)") # Table S2 B0 = 5.67 log10 CFU/mL
    log10bmax <- fixed(9.52); label("Bacterial density in stationary phase (log10 CFU/mL)") # Table S2 Bmax = 9.52 log10 CFU/mL

    # Per-node offsets of the baseline log10 density; 0 gives the homogeneous baseline
    dlog10b0_slab1 <- fixed(0); label("Offset of the node-1 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab2 <- fixed(0); label("Offset of the node-2 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab3 <- fixed(0); label("Offset of the node-3 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab4 <- fixed(0); label("Offset of the node-4 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab5 <- fixed(0); label("Offset of the node-5 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab6 <- fixed(0); label("Offset of the node-6 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab7 <- fixed(0); label("Offset of the node-7 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab8 <- fixed(0); label("Offset of the node-8 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab9 <- fixed(0); label("Offset of the node-9 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab10 <- fixed(0); label("Offset of the node-10 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab11 <- fixed(0); label("Offset of the node-11 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab12 <- fixed(0); label("Offset of the node-12 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab13 <- fixed(0); label("Offset of the node-13 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab14 <- fixed(0); label("Offset of the node-14 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab15 <- fixed(0); label("Offset of the node-15 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab16 <- fixed(0); label("Offset of the node-16 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab17 <- fixed(0); label("Offset of the node-17 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab18 <- fixed(0); label("Offset of the node-18 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab19 <- fixed(0); label("Offset of the node-19 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
    dlog10b0_slab20 <- fixed(0); label("Offset of the node-20 baseline from log10b0 (log10 CFU/mL)") # 0 = homogeneous baseline (Fig. S1); gradient baselines set per node
  })

  model({
    # Plasma PK (two-compartment, IV infusion into central)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc

    # Boundary condition at the blood-mucus interface (x = 0): plasma
    # concentration corrected by the ELF penetration coefficient
    Cmucus_bc <- relf * Cc

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

    # Bacterial PD rates converted to 1/h; densities in CFU/mL
    kgsh <- kgs * 3600
    kgrh <- kgr * 3600
    ksrh <- ksr * 3600
    slopesh <- slope_s * 3600
    sloperh <- slope_r * 3600
    bmax <- 10^log10bmax

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
    # Bacteria, nodes 1-20 (immobile); node 21 = node 20
    d/dt(bact_susceptible_slab1) <- (1 - (bact_susceptible_slab1 + bact_resistant_slab1) / bmax) * kgsh * bact_susceptible_slab1 - slopesh * Cmucus_bc * bact_susceptible_slab1 - ksrh * bact_susceptible_slab1
    d/dt(bact_resistant_slab1) <- (1 - (bact_susceptible_slab1 + bact_resistant_slab1) / bmax) * kgrh * bact_resistant_slab1 - sloperh * Cmucus_bc * bact_resistant_slab1 + ksrh * bact_susceptible_slab1
    d/dt(bact_susceptible_slab2) <- (1 - (bact_susceptible_slab2 + bact_resistant_slab2) / bmax) * kgsh * bact_susceptible_slab2 - slopesh * mucus_free_slab2 * bact_susceptible_slab2 - ksrh * bact_susceptible_slab2
    d/dt(bact_resistant_slab2) <- (1 - (bact_susceptible_slab2 + bact_resistant_slab2) / bmax) * kgrh * bact_resistant_slab2 - sloperh * mucus_free_slab2 * bact_resistant_slab2 + ksrh * bact_susceptible_slab2
    d/dt(bact_susceptible_slab3) <- (1 - (bact_susceptible_slab3 + bact_resistant_slab3) / bmax) * kgsh * bact_susceptible_slab3 - slopesh * mucus_free_slab3 * bact_susceptible_slab3 - ksrh * bact_susceptible_slab3
    d/dt(bact_resistant_slab3) <- (1 - (bact_susceptible_slab3 + bact_resistant_slab3) / bmax) * kgrh * bact_resistant_slab3 - sloperh * mucus_free_slab3 * bact_resistant_slab3 + ksrh * bact_susceptible_slab3
    d/dt(bact_susceptible_slab4) <- (1 - (bact_susceptible_slab4 + bact_resistant_slab4) / bmax) * kgsh * bact_susceptible_slab4 - slopesh * mucus_free_slab4 * bact_susceptible_slab4 - ksrh * bact_susceptible_slab4
    d/dt(bact_resistant_slab4) <- (1 - (bact_susceptible_slab4 + bact_resistant_slab4) / bmax) * kgrh * bact_resistant_slab4 - sloperh * mucus_free_slab4 * bact_resistant_slab4 + ksrh * bact_susceptible_slab4
    d/dt(bact_susceptible_slab5) <- (1 - (bact_susceptible_slab5 + bact_resistant_slab5) / bmax) * kgsh * bact_susceptible_slab5 - slopesh * mucus_free_slab5 * bact_susceptible_slab5 - ksrh * bact_susceptible_slab5
    d/dt(bact_resistant_slab5) <- (1 - (bact_susceptible_slab5 + bact_resistant_slab5) / bmax) * kgrh * bact_resistant_slab5 - sloperh * mucus_free_slab5 * bact_resistant_slab5 + ksrh * bact_susceptible_slab5
    d/dt(bact_susceptible_slab6) <- (1 - (bact_susceptible_slab6 + bact_resistant_slab6) / bmax) * kgsh * bact_susceptible_slab6 - slopesh * mucus_free_slab6 * bact_susceptible_slab6 - ksrh * bact_susceptible_slab6
    d/dt(bact_resistant_slab6) <- (1 - (bact_susceptible_slab6 + bact_resistant_slab6) / bmax) * kgrh * bact_resistant_slab6 - sloperh * mucus_free_slab6 * bact_resistant_slab6 + ksrh * bact_susceptible_slab6
    d/dt(bact_susceptible_slab7) <- (1 - (bact_susceptible_slab7 + bact_resistant_slab7) / bmax) * kgsh * bact_susceptible_slab7 - slopesh * mucus_free_slab7 * bact_susceptible_slab7 - ksrh * bact_susceptible_slab7
    d/dt(bact_resistant_slab7) <- (1 - (bact_susceptible_slab7 + bact_resistant_slab7) / bmax) * kgrh * bact_resistant_slab7 - sloperh * mucus_free_slab7 * bact_resistant_slab7 + ksrh * bact_susceptible_slab7
    d/dt(bact_susceptible_slab8) <- (1 - (bact_susceptible_slab8 + bact_resistant_slab8) / bmax) * kgsh * bact_susceptible_slab8 - slopesh * mucus_free_slab8 * bact_susceptible_slab8 - ksrh * bact_susceptible_slab8
    d/dt(bact_resistant_slab8) <- (1 - (bact_susceptible_slab8 + bact_resistant_slab8) / bmax) * kgrh * bact_resistant_slab8 - sloperh * mucus_free_slab8 * bact_resistant_slab8 + ksrh * bact_susceptible_slab8
    d/dt(bact_susceptible_slab9) <- (1 - (bact_susceptible_slab9 + bact_resistant_slab9) / bmax) * kgsh * bact_susceptible_slab9 - slopesh * mucus_free_slab9 * bact_susceptible_slab9 - ksrh * bact_susceptible_slab9
    d/dt(bact_resistant_slab9) <- (1 - (bact_susceptible_slab9 + bact_resistant_slab9) / bmax) * kgrh * bact_resistant_slab9 - sloperh * mucus_free_slab9 * bact_resistant_slab9 + ksrh * bact_susceptible_slab9
    d/dt(bact_susceptible_slab10) <- (1 - (bact_susceptible_slab10 + bact_resistant_slab10) / bmax) * kgsh * bact_susceptible_slab10 - slopesh * mucus_free_slab10 * bact_susceptible_slab10 - ksrh * bact_susceptible_slab10
    d/dt(bact_resistant_slab10) <- (1 - (bact_susceptible_slab10 + bact_resistant_slab10) / bmax) * kgrh * bact_resistant_slab10 - sloperh * mucus_free_slab10 * bact_resistant_slab10 + ksrh * bact_susceptible_slab10
    d/dt(bact_susceptible_slab11) <- (1 - (bact_susceptible_slab11 + bact_resistant_slab11) / bmax) * kgsh * bact_susceptible_slab11 - slopesh * mucus_free_slab11 * bact_susceptible_slab11 - ksrh * bact_susceptible_slab11
    d/dt(bact_resistant_slab11) <- (1 - (bact_susceptible_slab11 + bact_resistant_slab11) / bmax) * kgrh * bact_resistant_slab11 - sloperh * mucus_free_slab11 * bact_resistant_slab11 + ksrh * bact_susceptible_slab11
    d/dt(bact_susceptible_slab12) <- (1 - (bact_susceptible_slab12 + bact_resistant_slab12) / bmax) * kgsh * bact_susceptible_slab12 - slopesh * mucus_free_slab12 * bact_susceptible_slab12 - ksrh * bact_susceptible_slab12
    d/dt(bact_resistant_slab12) <- (1 - (bact_susceptible_slab12 + bact_resistant_slab12) / bmax) * kgrh * bact_resistant_slab12 - sloperh * mucus_free_slab12 * bact_resistant_slab12 + ksrh * bact_susceptible_slab12
    d/dt(bact_susceptible_slab13) <- (1 - (bact_susceptible_slab13 + bact_resistant_slab13) / bmax) * kgsh * bact_susceptible_slab13 - slopesh * mucus_free_slab13 * bact_susceptible_slab13 - ksrh * bact_susceptible_slab13
    d/dt(bact_resistant_slab13) <- (1 - (bact_susceptible_slab13 + bact_resistant_slab13) / bmax) * kgrh * bact_resistant_slab13 - sloperh * mucus_free_slab13 * bact_resistant_slab13 + ksrh * bact_susceptible_slab13
    d/dt(bact_susceptible_slab14) <- (1 - (bact_susceptible_slab14 + bact_resistant_slab14) / bmax) * kgsh * bact_susceptible_slab14 - slopesh * mucus_free_slab14 * bact_susceptible_slab14 - ksrh * bact_susceptible_slab14
    d/dt(bact_resistant_slab14) <- (1 - (bact_susceptible_slab14 + bact_resistant_slab14) / bmax) * kgrh * bact_resistant_slab14 - sloperh * mucus_free_slab14 * bact_resistant_slab14 + ksrh * bact_susceptible_slab14
    d/dt(bact_susceptible_slab15) <- (1 - (bact_susceptible_slab15 + bact_resistant_slab15) / bmax) * kgsh * bact_susceptible_slab15 - slopesh * mucus_free_slab15 * bact_susceptible_slab15 - ksrh * bact_susceptible_slab15
    d/dt(bact_resistant_slab15) <- (1 - (bact_susceptible_slab15 + bact_resistant_slab15) / bmax) * kgrh * bact_resistant_slab15 - sloperh * mucus_free_slab15 * bact_resistant_slab15 + ksrh * bact_susceptible_slab15
    d/dt(bact_susceptible_slab16) <- (1 - (bact_susceptible_slab16 + bact_resistant_slab16) / bmax) * kgsh * bact_susceptible_slab16 - slopesh * mucus_free_slab16 * bact_susceptible_slab16 - ksrh * bact_susceptible_slab16
    d/dt(bact_resistant_slab16) <- (1 - (bact_susceptible_slab16 + bact_resistant_slab16) / bmax) * kgrh * bact_resistant_slab16 - sloperh * mucus_free_slab16 * bact_resistant_slab16 + ksrh * bact_susceptible_slab16
    d/dt(bact_susceptible_slab17) <- (1 - (bact_susceptible_slab17 + bact_resistant_slab17) / bmax) * kgsh * bact_susceptible_slab17 - slopesh * mucus_free_slab17 * bact_susceptible_slab17 - ksrh * bact_susceptible_slab17
    d/dt(bact_resistant_slab17) <- (1 - (bact_susceptible_slab17 + bact_resistant_slab17) / bmax) * kgrh * bact_resistant_slab17 - sloperh * mucus_free_slab17 * bact_resistant_slab17 + ksrh * bact_susceptible_slab17
    d/dt(bact_susceptible_slab18) <- (1 - (bact_susceptible_slab18 + bact_resistant_slab18) / bmax) * kgsh * bact_susceptible_slab18 - slopesh * mucus_free_slab18 * bact_susceptible_slab18 - ksrh * bact_susceptible_slab18
    d/dt(bact_resistant_slab18) <- (1 - (bact_susceptible_slab18 + bact_resistant_slab18) / bmax) * kgrh * bact_resistant_slab18 - sloperh * mucus_free_slab18 * bact_resistant_slab18 + ksrh * bact_susceptible_slab18
    d/dt(bact_susceptible_slab19) <- (1 - (bact_susceptible_slab19 + bact_resistant_slab19) / bmax) * kgsh * bact_susceptible_slab19 - slopesh * mucus_free_slab19 * bact_susceptible_slab19 - ksrh * bact_susceptible_slab19
    d/dt(bact_resistant_slab19) <- (1 - (bact_susceptible_slab19 + bact_resistant_slab19) / bmax) * kgrh * bact_resistant_slab19 - sloperh * mucus_free_slab19 * bact_resistant_slab19 + ksrh * bact_susceptible_slab19
    d/dt(bact_susceptible_slab20) <- (1 - (bact_susceptible_slab20 + bact_resistant_slab20) / bmax) * kgsh * bact_susceptible_slab20 - slopesh * mucus_free_slab20 * bact_susceptible_slab20 - ksrh * bact_susceptible_slab20
    d/dt(bact_resistant_slab20) <- (1 - (bact_susceptible_slab20 + bact_resistant_slab20) / bmax) * kgrh * bact_resistant_slab20 - sloperh * mucus_free_slab20 * bact_resistant_slab20 + ksrh * bact_susceptible_slab20
    # Baseline inoculum: 10^log10b0 CFU/mL in every node, shifted per node by dlog10b0_slab<n>
    bact_susceptible_slab1(0) <- 10^(log10b0 + dlog10b0_slab1)
    bact_susceptible_slab2(0) <- 10^(log10b0 + dlog10b0_slab2)
    bact_susceptible_slab3(0) <- 10^(log10b0 + dlog10b0_slab3)
    bact_susceptible_slab4(0) <- 10^(log10b0 + dlog10b0_slab4)
    bact_susceptible_slab5(0) <- 10^(log10b0 + dlog10b0_slab5)
    bact_susceptible_slab6(0) <- 10^(log10b0 + dlog10b0_slab6)
    bact_susceptible_slab7(0) <- 10^(log10b0 + dlog10b0_slab7)
    bact_susceptible_slab8(0) <- 10^(log10b0 + dlog10b0_slab8)
    bact_susceptible_slab9(0) <- 10^(log10b0 + dlog10b0_slab9)
    bact_susceptible_slab10(0) <- 10^(log10b0 + dlog10b0_slab10)
    bact_susceptible_slab11(0) <- 10^(log10b0 + dlog10b0_slab11)
    bact_susceptible_slab12(0) <- 10^(log10b0 + dlog10b0_slab12)
    bact_susceptible_slab13(0) <- 10^(log10b0 + dlog10b0_slab13)
    bact_susceptible_slab14(0) <- 10^(log10b0 + dlog10b0_slab14)
    bact_susceptible_slab15(0) <- 10^(log10b0 + dlog10b0_slab15)
    bact_susceptible_slab16(0) <- 10^(log10b0 + dlog10b0_slab16)
    bact_susceptible_slab17(0) <- 10^(log10b0 + dlog10b0_slab17)
    bact_susceptible_slab18(0) <- 10^(log10b0 + dlog10b0_slab18)
    bact_susceptible_slab19(0) <- 10^(log10b0 + dlog10b0_slab19)
    bact_susceptible_slab20(0) <- 10^(log10b0 + dlog10b0_slab20)
    
    # Spatial summaries over the 21 grid nodes (node 21 = node 20)
    Cmucus_free_mean <- (Cmucus_bc + mucus_free_slab2 + mucus_free_slab3 + mucus_free_slab4 + mucus_free_slab5 + mucus_free_slab6 + mucus_free_slab7 + mucus_free_slab8 + mucus_free_slab9 + mucus_free_slab10 + mucus_free_slab11 + mucus_free_slab12 + mucus_free_slab13 + mucus_free_slab14 + mucus_free_slab15 + mucus_free_slab16 + mucus_free_slab17 + mucus_free_slab18 + mucus_free_slab19 + mucus_free_slab20 + mucus_free_slab20) / 21
    Cmucus_bound_mean <- (mucus_bound_slab1 + mucus_bound_slab2 + mucus_bound_slab3 + mucus_bound_slab4 + mucus_bound_slab5 + mucus_bound_slab6 + mucus_bound_slab7 + mucus_bound_slab8 + mucus_bound_slab9 + mucus_bound_slab10 + mucus_bound_slab11 + mucus_bound_slab12 + mucus_bound_slab13 + mucus_bound_slab14 + mucus_bound_slab15 + mucus_bound_slab16 + mucus_bound_slab17 + mucus_bound_slab18 + mucus_bound_slab19 + mucus_bound_slab20 + mucus_bound_slab20) / 21
    Cmucus_free_cv <- 100 * sqrt(((Cmucus_bc - Cmucus_free_mean)^2 + (mucus_free_slab2 - Cmucus_free_mean)^2 + (mucus_free_slab3 - Cmucus_free_mean)^2 + (mucus_free_slab4 - Cmucus_free_mean)^2 + (mucus_free_slab5 - Cmucus_free_mean)^2 + (mucus_free_slab6 - Cmucus_free_mean)^2 + (mucus_free_slab7 - Cmucus_free_mean)^2 + (mucus_free_slab8 - Cmucus_free_mean)^2 + (mucus_free_slab9 - Cmucus_free_mean)^2 + (mucus_free_slab10 - Cmucus_free_mean)^2 + (mucus_free_slab11 - Cmucus_free_mean)^2 + (mucus_free_slab12 - Cmucus_free_mean)^2 + (mucus_free_slab13 - Cmucus_free_mean)^2 + (mucus_free_slab14 - Cmucus_free_mean)^2 + (mucus_free_slab15 - Cmucus_free_mean)^2 + (mucus_free_slab16 - Cmucus_free_mean)^2 + (mucus_free_slab17 - Cmucus_free_mean)^2 + (mucus_free_slab18 - Cmucus_free_mean)^2 + (mucus_free_slab19 - Cmucus_free_mean)^2 + (mucus_free_slab20 - Cmucus_free_mean)^2 + (mucus_free_slab20 - Cmucus_free_mean)^2) / 20) / Cmucus_free_mean
    bact_total_mean <- (bact_susceptible_slab1 + bact_resistant_slab1 + bact_susceptible_slab2 + bact_resistant_slab2 + bact_susceptible_slab3 + bact_resistant_slab3 + bact_susceptible_slab4 + bact_resistant_slab4 + bact_susceptible_slab5 + bact_resistant_slab5 + bact_susceptible_slab6 + bact_resistant_slab6 + bact_susceptible_slab7 + bact_resistant_slab7 + bact_susceptible_slab8 + bact_resistant_slab8 + bact_susceptible_slab9 + bact_resistant_slab9 + bact_susceptible_slab10 + bact_resistant_slab10 + bact_susceptible_slab11 + bact_resistant_slab11 + bact_susceptible_slab12 + bact_resistant_slab12 + bact_susceptible_slab13 + bact_resistant_slab13 + bact_susceptible_slab14 + bact_resistant_slab14 + bact_susceptible_slab15 + bact_resistant_slab15 + bact_susceptible_slab16 + bact_resistant_slab16 + bact_susceptible_slab17 + bact_resistant_slab17 + bact_susceptible_slab18 + bact_resistant_slab18 + bact_susceptible_slab19 + bact_resistant_slab19 + bact_susceptible_slab20 + bact_resistant_slab20 + bact_susceptible_slab20 + bact_resistant_slab20) / 21
  })
}
