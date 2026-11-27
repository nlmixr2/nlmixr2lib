Hosoi_2026_levofloxacin_elf_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, ddY female, neutropenic). One-compartment PK model of",
    "levofloxacin in pulmonary EPITHELIAL LINING FLUID (ELF) after subcutaneous",
    "injection, with first-order absorption from the injection site directly",
    "into the ELF compartment and first-order elimination from it, in a murine",
    "Streptococcus pneumoniae ATCC 49619 lung infection model (Hosoi 2026). The",
    "ELF concentrations (bronchoalveolar lavage corrected by the plasma-to-BALF",
    "urea ratio) were fitted on their own, with no link to plasma, separately",
    "for each single-dose level (50, 100 and 200 mg/kg) by naive-pooled",
    "nonlinear least squares in Phoenix WinNonlin; the ini() values are the",
    "arithmetic means of the three fits that the paper tabulates (Table S1).",
    "The volume and clearance are apparent. The ELF PK/PD index retained by",
    "the paper, total ELF AUC over the 24 h treatment window divided by the",
    "MIC (auc_mic; MIC 0.512 mg/L), drives a sigmoid inhibitory Imax model of",
    "the 24 h change in lung bacterial burden (dlog10cfu, log10 CFU/lung,",
    "negative = kill; Table S4). dlog10cfu is the paper's 24 h endpoint and is",
    "meant to be read at t = 24 h after the start of treatment. Companion to",
    "the plasma model Hosoi_2026_levofloxacin_plasma_mouse. No between-animal",
    "variability or residual error was reported."
  )
  reference <- paste(
    "Hosoi A, Duan S, Shimada S, Suzuki K, Mizukami Y, Shintani R, Takahashi M,",
    "Satake N, Muraishi T, Okamoto Y, Igarashi Y, Enoki Y, Taguchi K, Matsumoto K.",
    "Three fluoroquinolones with different epithelial lining fluid penetration",
    "in a murine lung infection model: an experimental site-specific PK/PD",
    "study. Antimicrob Agents Chemother. 2026;70(8):e0037526.",
    "doi:10.1128/aac.00375-26. PMID: 42363833. PMCID: PMC13435811. PK",
    "parameters: Supplementary Table S1. Sigmoid Imax parameters:",
    "Supplementary Table S4. MIC: Results. ELF by urea dilution and model",
    "equations: Materials and Methods."
  )
  vignette <- "Hosoi_2026_fluoroquinolones_elf_mouse"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L (Celf); log10 CFU/lung (dlog10cfu)"
  )

  paper_specific_compartments <- c("auc_elf_0_24")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "levofloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    elf = list(analyte = "levofloxacin", units = "mg/kg", specimen = "epithelial lining fluid", verified = TRUE),
    auc_elf_0_24 = list(
      analyte = "levofloxacin",
      units = "mg*h/L",
      specimen = "epithelial lining fluid",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (ddY, female, 6 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5 weeks at purchase, used after a 1-week acclimation",
    sex_female_pct = 100,
    organism = "Streptococcus pneumoniae ATCC 49619; levofloxacin MIC 0.512 mg/L by CLSI broth microdilution",
    disease_state = paste(
      "Neutropenic (cyclophosphamide 150 mg/kg IP on day -4 and 100 mg/kg on",
      "day -1) mice with S. pneumoniae lung infection (intratracheal",
      "inoculation of 20 uL, about 1e7 CFU per mouse)"
    ),
    dose_range = paste(
      "PK: single subcutaneous doses of 50, 100 and 200 mg/kg. PD: various",
      "subcutaneous regimens given from 2 to 26 h after infection (24 h",
      "treatment window), lungs harvested at 26 h"
    ),
    regions = "Keio University, Tokyo, Japan",
    notes = paste(
      "Destructive sampling: three mice per group; after blood collection at",
      "5-300 min, bronchoalveolar lavage with 3 x 0.5 mL saline. ELF",
      "concentration = BALF concentration * (plasma urea / BALF urea). The ELF",
      "index is the TOTAL ELF AUC/MIC (no unbound-fraction correction). The",
      "sigmoid Imax fit of AUC/MIC had R^2 = 0.72 (Figure 1D); T>MIC and",
      "Cmax/MIC fits were also reported (Table S4) but AUC/MIC was the ELF",
      "index most strongly correlated with the antibacterial effect, and only",
      "it is packaged."
    )
  )

  ini({
    # Arithmetic means of the three single-dose fits (Table S1, 'Mean' row);
    # per-dose estimates are given on each line.
    lka <- log(8.83)
    label("Apparent absorption rate constant from the subcutaneous site into ELF, ka (1/h)") # Table S1 ELF ka Mean = 8.83 (7.27, 12.67, 6.56 at 50, 100, 200 mg/kg)
    lcl_elf <- log(2.24)
    label("Apparent clearance from the ELF compartment, CL (L/h/kg)") # Table S1 ELF CL Mean = 2.24 (2.95, 1.93, 1.82)
    lv_elf <- log(2.60)
    label("Apparent ELF volume of distribution, Vd (L/kg)") # Table S1 ELF Vd Mean = 2.60 (1.75, 2.65, 3.39)

    mic <- fixed(0.512)
    label("Levofloxacin MIC against S. pneumoniae ATCC 49619 (mg/L)") # Results 'MIC': 0.512 mg/L

    # Sigmoid inhibitory Imax model of the ELF AUC/MIC (Materials and
    # Methods 'PK/PD analysis'): E = E0 - Imax * X^gamma / (X^gamma + IC50^gamma).
    e0 <- 0.76
    label("Change in lung burden at zero exposure, E0 (log10 CFU/lung)") # Table S4 ELF AUC/MIC E0 = 0.76
    limax <- log(6.82)
    label("Maximum antibacterial effect, Imax (log10 CFU/lung)") # Table S4 ELF AUC/MIC Imax = 6.82
    lic50 <- log(125.08)
    label("ELF AUC/MIC giving half the maximum effect, IC50 (h)") # Table S4 ELF AUC/MIC IC50 = 125.08
    lhill <- log(1.12)
    label("Hill coefficient of the ELF AUC/MIC effect curve (unitless)") # Table S4 ELF AUC/MIC gamma = 1.12

    propSd_Celf <- fixed(0)
    label("Proportional residual SD on Celf (fraction); not reported for the least-squares fit") # Materials and Methods 'Pharmacokinetic analysis': nonlinear least squares, no error model reported
  })

  model({
    ka <- exp(lka)
    cl_elf <- exp(lcl_elf)
    v_elf <- exp(lv_elf)
    kel_elf <- cl_elf / v_elf

    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(elf) <- ka * depot - kel_elf * elf

    Celf <- elf / v_elf

    # Total ELF AUC over the 24 h treatment window; constant after 24 h.
    d/dt(auc_elf_0_24) <- Celf * (t < 24)
    auc_mic <- auc_elf_0_24 / mic

    # 24 h change in lung bacterial burden (log10 CFU/lung): e0 with no drug,
    # approaching e0 - imax at saturating exposure.
    dlog10cfu <- e0 - imax * auc_mic^hill / (auc_mic^hill + ic50^hill)

    Celf ~ prop(propSd_Celf)
  })
}
