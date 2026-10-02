ReigLopez_2020_mbq167_her2_mouse_tgi <- function() {
  description <- paste(
    "Preclinical (mouse bearing an orthotopic GFP-MDA-MB-435 HER2+ mammary",
    "fat pad tumor). Simeoni 2004 tumor-growth-inhibition model of the Rac/Cdc42",
    "inhibitor MBQ-167: generalized exponential-then-linear unperturbed growth,",
    "a sigmoid Emax (Kmax, IC50, Hill) kill rate driven by the total plasma",
    "MBQ-167 concentration, and three damaged-cell transit compartments before",
    "cell death. PD-only: the paper's Simcyp Animal V19 whole-body PBPK model is",
    "not portable, so the plasma concentration enters as the time-varying",
    "covariate CP_MBQ167_NGML."
  )
  reference <- paste(
    "Reig-Lopez J, Maldonado MdM, Merino-Sanjuan M, Cruz-Collazo AM,",
    "Ruiz-Calderon JF, Mangas-Sanjuan V, Dharmawardhane S, Duconge J.",
    "Physiologically-Based Pharmacokinetic/Pharmacodynamic Model of MBQ-167 to",
    "Predict Tumor Growth Inhibition in Mice. Pharmaceutics. 2020;12(10):975.",
    "doi:10.3390/pharmaceutics12100975. Parameter values from Table 2 and the",
    "drug-effect equation from Figure 1. Structural TGI formulation from Simeoni M,",
    "Magni P, Cammia C, De Nicolao G, Croci V, Pesenti E, Germani M, Poggesi I,",
    "Rocchetti M. Predictive pharmacokinetic-pharmacodynamic modeling of tumor",
    "growth kinetics in xenograft models after administration of anticancer",
    "agents. Cancer Res. 2004;64(3):1094-1101. doi:10.1158/0008-5472.CAN-03-2524."
  )
  vignette <- "ReigLopez_2020_mbq167"

  units <- list(
    time = "day",
    dosing = "none (PD-only; MBQ-167 exposure is supplied as the CP_MBQ167_NGML covariate)",
    concentration = "ng/mL (CP_MBQ167_NGML); tumor volume in mL"
  )

  compartmentData <- list(
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mL", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "drug-damaged tumor cells (transit state 1)",
      units = "mL",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "drug-damaged tumor cells (transit state 2)",
      units = "mL",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "drug-damaged tumor cells (transit state 3)",
      units = "mL",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_MBQ167_NGML = list(
      description = paste(
        "Total (bound + unbound) MBQ-167 plasma concentration. Reig-Lopez 2020 drives",
        "the tumor-growth-inhibition model with the total plasma concentration",
        "predicted by its Simcyp Animal V19 whole-body PBPK model (Table 2 row 'Drug",
        "input'; Figure 1 'Cp(t)'); that PBPK model is not reproducible outside the",
        "platform, so the user supplies the trajectory. Set to 0 for vehicle controls."
      ),
      units = "ng/mL",
      type = "continuous",
      source_name = "Cp(t)",
      reference_category = NA_character_,
      notes = paste(
        "Time-varying. Converted inside the model to umol/L with the Table 1",
        "molecular weight (338.414 g/mol) because IC50 is reported in uM. The",
        "paper's predicted typical plasma profile after a single 10 mg/kg IP dose in",
        "a 20 g BALB/c mouse has Cmax 833.31 ng/mL, AUC0-12h 1549.1 ng*h/mL (Table",
        "3), Tmax 0.26 h and terminal half-life 2.98 h (Discussion). Because the Hill",
        "coefficient is 0.5, the kill rate is sensitive to the low-concentration tail",
        "of each dosing interval, so the trajectory should be supplied on a fine",
        "time grid and with linear covariate interpolation."
      )
    )
  )

  population <- list(
    species = "mouse (female athymic nude nu/nu, orthotopic GFP-MDA-MB-435 HER2+ mammary fat pad tumor)",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "4-5 weeks at purchase",
    sex_female_pct = 100,
    disease_state = "orthotopic HER2+ (MDA-MB-435, described as HER2++) human breast cancer xenograft",
    dose_range = "vehicle, 1 or 10 mg/kg MBQ-167 IP every other day, three times a week, until sacrifice at day 65",
    regions = "University of Puerto Rico Medical Sciences Campus (preclinical)",
    notes = paste(
      "n = 10 mice per treatment group (Methods 2.6); the tumor-growth data come",
      "from the earlier Humphries-Bickley 2017 experiment (reference 18 of the",
      "paper). Tumor growth was quantified as GFP fluorescence integrated density",
      "and reported as tumor volume (mL). The fit is a typical-mouse fit (Simcyp",
      "parameter estimation, weighted least squares); no inter-individual or",
      "residual variability is reported. The 5 mg/kg arm described in Methods",
      "2.6 is not modelled."
    )
  )

  ini({
    # Table 2, HER2+ column. Superscripts: a = assumed, b = estimated,
    # c = optimized to best fit the observed data.
    lrbase_tumor <- fixed(log(0.1)); label("Initial tumor volume at the start of treatment (mL)") # Table 2 'Initial tumor volume (mL)' = 0.1 (a: assumed)
    ltumorExpGrowth <- log(0.2); label("Tumor growth rate during the initial exponential phase, lambda0 (1/day)") # Table 2 'lambda0 (day-1)' = 0.2 (c)
    ltumorLinGrowth <- log(0.12); label("Tumor growth rate during the later linear phase, lambda1 (mL/day)") # Table 2 'lambda1 (g/day)' = 0.12 (c); g read as mL
    psi <- 0.7; label("Shape factor of the exponential-to-linear growth switch, Psi (unitless)") # Table 2 'Psi' = 0.7 (c)
    ldamageTransit <- log(0.39); label("Transit rate of damaged cells through the death chain, k1 (1/day)") # Table 2 'k1 (day-1)' = 0.39 (c)
    lic50 <- log(0.0187); label("Total plasma MBQ-167 concentration giving half of Kmax, IC50 (umol/L)") # Table 2 'IC50 (uM)' = 0.0187 (b)
    lkmax <- log(0.3683); label("Maximum drug-induced kill rate of cycling cells, Kmax (1/day)") # Table 2 'Kmax (day-1)' = 0.3683 (b)
    lhill <- log(0.5); label("Hill coefficient of the drug effect (unitless)") # Table 2 'H' = 0.5 (c)

    # No residual error is reported: Simcyp parameter estimation used weighted
    # least squares on the mean profile (Methods 2.8.1). Carried as fixed(0).
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on tumor volume (fraction; not reported in the source)")
  })

  model({
    rbase_tumor <- exp(lrbase_tumor)
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    damageTransit <- exp(ldamageTransit)
    ic50 <- exp(lic50)
    kmax <- exp(lkmax)
    hill <- exp(lhill)

    # Total plasma concentration (ng/mL) to umol/L; MW 338.414 g/mol (Table 1).
    cp_um <- CP_MBQ167_NGML / 338.414

    # Drug effect from Figure 1: Kmax * Cp^H / (IC50^H + Cp^H), acting as the
    # first-order rate at which cycling cells enter the damaged-cell chain.
    killRate <- kmax * cp_um^hill / (ic50^hill + cp_um^hill)

    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    # Simeoni 2004 unperturbed growth with the shape factor psi estimated
    # (Table 2) rather than fixed at 20. HER2+ uses three transit compartments
    # (Table 2 'Number of transit compartments' = 3; Figure 1 TS1-TS3).
    d/dt(cycling_cells) <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth / tumorLinGrowth * tumor_vol)^psi)^(1 / psi) -
      killRate * cycling_cells
    d/dt(damaged_cells1) <- killRate * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    cycling_cells(0) <- rbase_tumor

    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
