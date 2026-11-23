vandeKreeke_2026_antibiotic_immune_bacteremia <- function() {
  description <- paste(
    "QSP. Mechanism-based host-pathogen model of bacteremia with a generic",
    "(hypothetical) antibiotic (van de Kreeke 2026). Extracellular bacteria",
    "(bact_susceptible, CFU/mL) grow logistically at kgrowth up to bmax and",
    "are phagocytosed by neutrophils and monocytes into intracellular pools",
    "(bact_phagocytosed_neut, bact_phagocytosed_mono) from which they are",
    "digested. Each phagocytosis and digestion rate is an Emax function of",
    "the immune-cell concentration times a capacity term that falls as the",
    "intracellular-bacteria-to-cell ratio rises (hyperbolic for neutrophils,",
    "sigmoid with a Hill exponent for monocytes). Neutrophil (NEUT) and",
    "monocyte (MONO_ABS) counts are constant per subject, as in the",
    "calibrated in vivo model. The antibiotic follows one-compartment PK",
    "(V 5 L, half-life 5 h, no protein binding) and lowers the bacterial",
    "growth rate through the Regoes 2004 function of Cc / MIC; kgmin and",
    "hill select the antibiotic type (bacteriostatic kgmin = -1 or",
    "bactericidal -3 1/h; time-dependent hill = 0.5 or concentration-",
    "dependent 3). The default is bacteriostatic time-dependent; override",
    "kgmin and hill through ini() for the other three types. The printed",
    "immune-capacity equations reproduce the paper's antibiotic-only and",
    "low-inoculum immune results but under-predict immune clearance at high",
    "bacterial load (see the vignette)."
  )
  reference <- paste(
    "van de Kreeke M, Pham AD, Mehciz M, Jordens T, Uyterlinde MG, van Wijk",
    "RC, Aulin LBS, Knibbe CAJ, Krekels EHJ, Zwep LB, van Hasselt JGC. A",
    "mechanism-based model of immune status effects on antibiotic PK/PD",
    "targets in bacteremia. J Pharmacokinet Pharmacodyn. 2026.",
    "doi:10.1007/s10928-026-10063-6. PMID: 42771197. PMCID: PMC13597578.",
    "Parameter values: Table 1; equations: Eqs. 1, 7, 8, 9, 11 and 12 and",
    "Fig. 1; in vitro values and calibration: Online Resources 1 and 4.",
    "Antibiotic effect: Regoes RR et al. Antimicrob Agents Chemother.",
    "2004;48:3670-3676.",
    sep = " "
  )
  vignette <- "vandeKreeke_2026_antibiotic_immune_bacteremia"
  # central: antibiotic dose; bact_susceptible: the bacterial inoculum.
  dosing <- c("central", "bact_susceptible")

  units <- list(
    time = "h",
    dosing = "mg (antibiotic); CFU/mL (bacterial inoculum)",
    concentration = "mg/L (antibiotic, equal to ug/mL); CFU/mL (bacteria)"
  )

  covariateData <- list(
    NEUT = list(
      description = "Absolute neutrophil count, held constant over the simulation",
      units = "cells/mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converted in model() to the paper's cells/mL (N = NEUT * 1000).",
        "Paper immune states (Methods, Evaluation of immune status): immune",
        "competent 5000; mild neutropenia 1250; moderate 750; severe 500;",
        "fully immune compromised 0 cells/mm^3. Mouse calibration (Online",
        "Resource 1, Table S3): 170 cells/mm^3."
      ),
      source_name = "N"
    ),
    MONO_ABS = list(
      description = "Absolute monocyte count, held constant over the simulation",
      units = "cells/mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converted in model() to the paper's cells/mL (M = MONO_ABS * 1000).",
        "Paper immune states: immune competent 600; monocytopenia 200; fully",
        "immune compromised 0 cells/mm^3. Mouse calibration (Table S3): 1400",
        "cells/mm^3."
      ),
      source_name = "M"
    ),
    MIC = list(
      description = "Minimum inhibitory concentration of the antibiotic against the infecting bacterium",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper sets MIC = 1 ug/mL (1 mg/L) for every simulation, with a",
        "homogeneous bacterial population and no resistance development",
        "(Methods, Eq. 7 and Evaluation of immune status)."
      ),
      source_name = "MIC"
    )
  )

  # The two intracellular pools are the bacteria held inside neutrophils and
  # inside monocytes (B_N and B_M of Eqs. 8-12), not drug-resistance
  # phenotypes of the registered bact_<phenotype> family.
  paper_specific_compartments <- c(
    "bact_phagocytosed_neut",
    "bact_phagocytosed_mono"
  )

  compartmentData <- list(
    central = list(analyte = "generic antibiotic", units = "mg", specimen = "plasma", verified = TRUE),
    bact_susceptible = list(
      analyte = "extracellular bacteria (B)",
      units = "CFU/mL",
      specimen = "whole blood",
      verified = TRUE
    ),
    bact_phagocytosed_neut = list(
      analyte = "bacteria phagocytosed by neutrophils and not yet digested (B_N)",
      units = "CFU/mL",
      specimen = "whole blood",
      verified = TRUE
    ),
    bact_phagocytosed_mono = list(
      analyte = "bacteria phagocytosed by monocytes and not yet digested (B_M)",
      units = "CFU/mL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human + mouse",
    n_subjects = 4L,
    n_studies = 1L,
    disease_state = paste(
      "Bacteremia with a generic pathogen. Simulated adults in seven immune",
      "states: immune competent, monocytopenia, mild, moderate and severe",
      "neutropenia, severe neutropenia with monocytopenia, and fully immune",
      "compromised."
    ),
    dose_range = paste(
      "Generic antibiotic, 0-100 mg total daily dose given as 1, 2, 3 or 6",
      "intravenous boluses over 24 h; initial bacterial density 10^9 CFU/mL."
    ),
    notes = paste(
      "No patient data were fitted. Phagocytosis and digestion parameters",
      "were estimated from digitised in vitro neutrophil and monocyte",
      "time-courses (Online Resource 1, Table S2; Online Resources 2 and 3),",
      "and one scaling factor on both maximum phagocytosis rates (0.227) was",
      "estimated from 24-h blood CFU in four immune-competent mice",
      "intravenously inoculated with E. coli (4 x 10^7 CFU/mouse; Table S3,",
      "Online Resource 4). n_subjects counts those four mice. The human",
      "simulations use kgrowth 0.9 1/h and a 10^9 CFU/mL inoculum selected",
      "in Online Resource 1, Fig. S1."
    )
  )

  ini({
    # ---- Bacteria (Table 1) --------------------------------------------------
    lkgrowth <- fixed(log(0.9)); label("Net bacterial growth rate without antibiotic, kgrw (1/h)") # Table 1: kgrw = 0.9 /h (clinical value, footnote a; mouse calibration used 0.5885)
    lbmax <- fixed(log(2.5e9)); label("Maximum bacterial density, Bmax (CFU/mL)") # Table 1: Bmax = 2.5 x 10^9 CFU/mL

    # ---- Neutrophils (Table 1) -----------------------------------------------
    lknp_max <- log(1.90); label("Maximum neutrophil phagocytosis rate constant, knp,max (1/h)") # Table 1: knp,max = 1.90 /h (in vitro 8.35 x scale 0.227, Online Resource 4: 1.895)
    ln50_phag <- log(6.92e5); label("Neutrophil concentration giving half-maximum phagocytosis, Np50 (cells/mL)") # Table 1: Np50 = 6.92 x 10^5 cell/mL
    lbn50_phag <- log(3.92); label("Intracellular-bacteria-to-neutrophil ratio halving phagocytosis, B_N/Np50 (CFU/cell)") # Table 1: B_N/Np50 = 3.92 CFU/cell
    lknd_max <- log(5.61); label("Maximum neutrophil digestion rate constant, knd,max (1/h)") # Table 1: knd,max = 5.61 /h
    ln50_dig <- log(3.98e6); label("Neutrophil concentration giving half-maximum digestion, Nd50 (cells/mL)") # Table 1: Nd50 = 3.98 x 10^6 cell/mL
    lbn50_dig <- log(2.84); label("Intracellular-bacteria-to-neutrophil ratio halving digestion, B_N/Nd50 (CFU/cell)") # Table 1: B_N/Nd50 = 2.84 CFU/cell

    # ---- Monocytes (Table 1) -------------------------------------------------
    lkmp_max <- log(1.77); label("Maximum monocyte phagocytosis rate constant, kmp,max (1/h)") # Table 1: kmp,max = 1.77 /h (in vitro 7.78 x scale 0.227, Online Resource 4: 1.766)
    lm50_phag <- log(2.73e5); label("Monocyte concentration giving half-maximum phagocytosis, Mp50 (cells/mL)") # Table 1: Mp50 = 2.73 x 10^5 cell/mL
    lbm50_phag <- log(1.25); label("Intracellular-bacteria-to-monocyte ratio halving phagocytosis, B_M/Mp50 (CFU/cell)") # Table 1: B_M/Mp50 = 1.25 CFU/cell
    lhill_mono_phag <- log(5.54); label("Hill exponent of the monocyte phagocytosis capacity term, gamma_p (unitless)") # Table 1: gamma_p = 5.54
    lkmd_max <- log(6.23); label("Maximum monocyte digestion rate constant, kmd,max (1/h)") # Table 1: kmd,max = 6.23 /h
    lm50_dig <- fixed(log(3.98e6)); label("Monocyte concentration giving half-maximum digestion, Md50 (cells/mL)") # Table 1: Md50 = 3.98 x 10^6 cell/mL; Results: set equal to Nd50 (no data)
    lbm50_dig <- log(1.05); label("Intracellular-bacteria-to-monocyte ratio halving digestion, B_M/Md50 (CFU/cell)") # Table 1: B_M/Md50 = 1.05 CFU/cell
    lhill_mono_dig <- log(14); label("Hill exponent of the monocyte digestion capacity term, gamma_d (unitless)") # Table 1: gamma_d = 14

    # ---- Hypothetical antibiotic (Methods, Mechanism-based model) ------------
    lvc <- fixed(log(5)); label("Antibiotic volume of distribution (L)") # Methods: one-compartment, distribution volume of 5 L
    lcl <- fixed(log(5 * log(2) / 5)); label("Antibiotic clearance (L/h)") # Methods: 5-h half-life, so CL = V * ln(2) / 5 h = 0.693 L/h
    kgmin <- fixed(-1); label("Bacterial growth rate at maximal antibiotic effect, kgmin (1/h)") # Methods, Eq. 7: kgmin = -1 /h bacteriostatic, -3 /h bactericidal
    lhill <- fixed(log(0.5)); label("Hill coefficient of the antibiotic effect on growth, H (unitless)") # Methods, Eq. 7: H = 0.5 time-dependent, 3 concentration-dependent
  })

  model({
    kgrowth <- exp(lkgrowth)
    bmax <- exp(lbmax)
    knp_max <- exp(lknp_max)
    n50_phag <- exp(ln50_phag)
    bn50_phag <- exp(lbn50_phag)
    knd_max <- exp(lknd_max)
    n50_dig <- exp(ln50_dig)
    bn50_dig <- exp(lbn50_dig)
    kmp_max <- exp(lkmp_max)
    m50_phag <- exp(lm50_phag)
    bm50_phag <- exp(lbm50_phag)
    hill_mono_phag <- exp(lhill_mono_phag)
    kmd_max <- exp(lkmd_max)
    m50_dig <- exp(lm50_dig)
    bm50_dig <- exp(lbm50_dig)
    hill_mono_dig <- exp(lhill_mono_dig)
    vc <- exp(lvc)
    cl <- exp(lcl)
    kel <- cl / vc
    hill <- exp(lhill)

    # Immune-cell concentrations in the paper's cells/mL, constant over time
    # (Calibration: in vivo the apoptosis, phagocytosis-induced death and
    # time-driven decay terms of Eqs. 5, 6, 10, 11 and 12 were removed).
    neut <- NEUT * 1000
    mono <- MONO_ABS * 1000

    # Intracellular-bacteria-to-cell ratios B_N/N and B_M/M (CFU/cell); zero
    # when the cell type is absent, where the Emax term already switches the
    # process off.
    ratio_neut <- 0
    if (neut > 0) ratio_neut <- bact_phagocytosed_neut / neut
    ratio_mono <- 0
    if (mono > 0) ratio_mono <- bact_phagocytosed_mono / mono

    # Eqs. 8 and 9. Eq. 9 prints B/N in the last denominator term; the
    # intracellular ratio B_N/N is used as in Eq. 8 and Eq. 4.
    knp <- knp_max * neut / (n50_phag + neut) * (1 - ratio_neut / (bn50_phag + ratio_neut))
    knd <- knd_max * neut / (n50_dig + neut) * (1 - ratio_neut / (bn50_dig + ratio_neut))

    # Eqs. 11 and 12 without the time-driven decay (removed at calibration).
    # Both print B/M in the last denominator term; B_M/M is used as in Eq. 4.
    kmp <- kmp_max * mono / (m50_phag + mono) *
      (1 - ratio_mono^hill_mono_phag / (bm50_phag^hill_mono_phag + ratio_mono^hill_mono_phag))
    kmd <- kmd_max * mono / (m50_dig + mono) *
      (1 - ratio_mono^hill_mono_dig / (bm50_dig^hill_mono_dig + ratio_mono^hill_mono_dig))

    # Antibiotic PK (one compartment, no protein binding).
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 7: Regoes 2004 net growth rate under antibiotic.
    cmic <- (Cc / MIC)^hill
    kg <- kgrowth - (kgrowth - kgmin) * cmic / (cmic - kgmin / kgrowth)

    # Eq. 1 and Fig. 1. The intracellular pools have no printed equation; each
    # gains bacteria by phagocytosis and loses them by digestion (Fig. 1).
    d/dt(bact_susceptible) <- kg * bact_susceptible * (1 - bact_susceptible / bmax) -
      knp * bact_susceptible - kmp * bact_susceptible
    d/dt(bact_phagocytosed_neut) <- knp * bact_susceptible - knd * bact_phagocytosed_neut
    d/dt(bact_phagocytosed_mono) <- kmp * bact_susceptible - kmd * bact_phagocytosed_mono

    # Bacterial density (CFU/mL): the extracellular bacteria B of Eq. 1.
    # cfu_total adds the viable bacteria still held inside phagocytes.
    cfu <- bact_susceptible
    cfu_total <- bact_susceptible + bact_phagocytosed_neut + bact_phagocytosed_mono
  })
}
