Hosoi_2026_lascufloxacin_plasma_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, ddY female, neutropenic). One-compartment PK model of",
    "total lascufloxacin in PLASMA after subcutaneous injection, with first-order",
    "absorption and first-order elimination, in a murine Streptococcus",
    "pneumoniae ATCC 49619 lung infection model (Hosoi 2026). Each single-dose",
    "level (10, 20 and 30 mg/kg) was fitted separately by naive-pooled",
    "nonlinear least squares in Phoenix WinNonlin; the ini() values are the",
    "arithmetic means of the three fits that the paper tabulates (Table S3).",
    "The volume and clearance are apparent (V/F, CL/F) because subcutaneous",
    "bioavailability was not estimated. Free plasma concentration Cu = fu * Cc",
    "with the measured murine unbound fraction fu = 0.354 (64.6% bound). The",
    "plasma PK/PD index retained by the paper, free AUC over the 24 h",
    "treatment window divided by the MIC (fauc_mic; MIC 0.032 mg/L), drives a",
    "sigmoid inhibitory Imax model of the 24 h change in lung bacterial burden",
    "(dlog10cfu, log10 CFU/lung, negative = kill; Table S6). dlog10cfu is the",
    "paper's 24 h endpoint and is meant to be read at t = 24 h after the start",
    "of treatment; before 24 h it reflects only the exposure accumulated so",
    "far. Companion to the ELF model Hosoi_2026_lascufloxacin_elf_mouse, which",
    "was fitted separately to the epithelial lining fluid concentrations. No",
    "between-animal variability or residual error was reported."
  )
  reference <- paste(
    "Hosoi A, Duan S, Shimada S, Suzuki K, Mizukami Y, Shintani R, Takahashi M,",
    "Satake N, Muraishi T, Okamoto Y, Igarashi Y, Enoki Y, Taguchi K, Matsumoto K.",
    "Three fluoroquinolones with different epithelial lining fluid penetration",
    "in a murine lung infection model: an experimental site-specific PK/PD",
    "study. Antimicrob Agents Chemother. 2026;70(8):e0037526.",
    "doi:10.1128/aac.00375-26. PMID: 42363833. PMCID: PMC13435811. PK",
    "parameters: Supplementary Table S3. Sigmoid Imax parameters:",
    "Supplementary Table S6. MIC and protein binding: Results. Model",
    "equations: Materials and Methods."
  )
  vignette <- "Hosoi_2026_fluoroquinolones_elf_mouse"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L (Cc, Cu); log10 CFU/lung (dlog10cfu)"
  )

  paper_specific_compartments <- c("fauc_0_24")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "lascufloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lascufloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    fauc_0_24 = list(analyte = "lascufloxacin (unbound)", units = "mg*h/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (ddY, female, 6 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5 weeks at purchase, used after a 1-week acclimation",
    sex_female_pct = 100,
    organism = "Streptococcus pneumoniae ATCC 49619; lascufloxacin MIC 0.032 mg/L by CLSI broth microdilution",
    disease_state = paste(
      "Neutropenic (cyclophosphamide 150 mg/kg IP on day -4 and 100 mg/kg on",
      "day -1) mice with S. pneumoniae lung infection (intratracheal",
      "inoculation of 20 uL, about 1e7 CFU per mouse)"
    ),
    dose_range = paste(
      "PK: single subcutaneous doses of 10, 20 and 30 mg/kg. PD: various",
      "subcutaneous regimens given from 2 to 26 h after infection (24 h",
      "treatment window), lungs harvested at 26 h"
    ),
    regions = "Keio University, Tokyo, Japan",
    notes = paste(
      "Destructive sampling: three mice per group, blood from the inferior vena",
      "cava at 5, 15, 30, 60, 120, 240 and 300 min after dosing, then",
      "bronchoalveolar lavage. Plasma protein binding by ultrafiltration 64.6%",
      "(SD 1.06). The PK model was fitted to TOTAL plasma concentrations; the",
      "paper's fAUC/MIC uses fu * Cc. The sigmoid Imax fit of fAUC/MIC had",
      "R^2 = 0.85 (Figure 3A); fT>MIC and fCmax/MIC fits were also reported",
      "(Table S6) but fAUC/MIC was the index most strongly correlated with",
      "the antibacterial effect, and only it is packaged."
    )
  )

  ini({
    # Arithmetic means of the three single-dose fits (Table S3, 'Mean' row);
    # per-dose estimates are given on each line.
    lka <- log(3.55)
    label("Subcutaneous absorption rate constant ka (1/h)") # Table S3 Plasma ka Mean = 3.55 (3.17, 1.43, 6.06 at 10, 20, 30 mg/kg)
    lcl <- log(3.96)
    label("Apparent clearance CL/F (L/h/kg)") # Table S3 Plasma CL Mean = 3.96 (3.88, 4.65, 3.35)
    lvc <- log(6.56)
    label("Apparent volume of distribution Vd/F (L/kg)") # Table S3 Plasma Vd Mean = 6.56 (4.00, 5.13, 10.56)
    fu <- fixed(0.354)
    label("Unbound fraction of lascufloxacin in mouse plasma (unitless)") # Results: %PB = 64.6% (SD 1.06), so fu = 1 - 0.646

    mic <- fixed(0.032)
    label("Lascufloxacin MIC against S. pneumoniae ATCC 49619 (mg/L)") # Results 'MIC': 0.032 mg/L

    # Sigmoid inhibitory Imax model of the plasma fAUC/MIC (Materials and
    # Methods 'PK/PD analysis'): E = E0 - Imax * X^gamma / (X^gamma + IC50^gamma).
    e0 <- 0.27
    label("Change in lung burden at zero exposure, E0 (log10 CFU/lung)") # Table S6 Plasma fAUC/MIC E0 = 0.27
    limax <- log(3.27)
    label("Maximum antibacterial effect, Imax (log10 CFU/lung)") # Table S6 Plasma fAUC/MIC Imax = 3.27
    lic50 <- log(18.47)
    label("fAUC/MIC giving half the maximum effect, IC50 (h)") # Table S6 Plasma fAUC/MIC IC50 = 18.47
    lhill <- log(3.93)
    label("Hill coefficient of the fAUC/MIC effect curve (unitless)") # Table S6 Plasma fAUC/MIC gamma = 3.93

    propSd <- fixed(0)
    label("Proportional residual SD on Cc (fraction); not reported for the least-squares fit") # Materials and Methods 'Pharmacokinetic analysis': nonlinear least squares, no error model reported
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc

    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cu <- fu * Cc

    # Free plasma AUC over the 24 h treatment window; constant after 24 h.
    d/dt(fauc_0_24) <- Cu * (t < 24)
    fauc_mic <- fauc_0_24 / mic

    # 24 h change in lung bacterial burden (log10 CFU/lung): e0 with no drug,
    # approaching e0 - imax at saturating exposure.
    dlog10cfu <- e0 - imax * fauc_mic^hill / (fauc_mic^hill + ic50^hill)

    Cc ~ prop(propSd)
  })
}
