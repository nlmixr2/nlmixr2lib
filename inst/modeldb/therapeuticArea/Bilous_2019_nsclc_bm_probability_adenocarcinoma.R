Bilous_2019_nsclc_bm_probability_adenocarcinoma <- function() {
  description <- paste(
    "Population mechanistic model of the probability of brain metastasis at",
    "diagnosis as a function of primary tumor size in lung adenocarcinoma",
    "(non-small cell lung cancer; Bilous 2019, Figure 1A). Time 0 is the",
    "first cancer cell. The untreated primary tumor (in cells) grows by",
    "Gompertz kinetics calibrated to a 1e12-cell carrying capacity and the",
    "adenocarcinoma volume doubling time of 201 days at a 35-mm cohort",
    "reference diameter; brain metastases are born at rate mu * tumor^gamma",
    "with gamma = 1 (Iwata framework). The cellular metastatic potential mu",
    "is log-normal between patients. A patient has brain metastasis at",
    "diagnosis when the number of metastases born by the diagnosis time td",
    "(set by the primary tumor diameter TUMSZ) exceeds 1. No drug input.",
    sep = " "
  )
  reference <- paste(
    "Bilous M, Serdjebi C, Boyer A, Tomasini P, Pouypoudat C, Barbolosi D,",
    "Barlesi F, Chomy F, Benzekry S. Quantitative mathematical modeling of",
    "clinical brain metastasis dynamics in non-small cell lung cancer.",
    "Sci Rep. 2019;9:13018. doi:10.1038/s41598-019-49407-3.",
    "Equations are from Materials and Methods (equations 2 and 4 and",
    "'Population level: probability of BM apparition'); mu_pop and mu_sigma",
    "are from the Figure 1A legend; the doubling time is from Supplementary",
    "Table S1. The data fitted are from Mujoomdar A et al. Radiology",
    "2007;242:882-888 (doi:10.1148/radiol.2423051707).",
    sep = " "
  )
  vignette <- "Bilous_2019_nsclc_brain_metastasis"

  units <- list(
    time = "day",
    dosing = "none (no drug input)",
    concentration = "cells (tumor size) and count (number of brain metastases); there is no drug concentration"
  )

  paper_specific_compartments <- c("metastases")

  compartmentData <- list(
    tumor = list(
      analyte = "primary lung tumor size",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    metastases = list(
      analyte = "cumulative number of brain metastases born, N(t)",
      units = "count",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUMSZ = list(
      description = "Largest diameter of the primary lung tumor at diagnosis",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets only the age of the tumor at diagnosis td (the time at which",
        "the cohort Gompertz curve reaches this size, assuming a sphere and",
        "1e6 cells/mm^3). The Gompertz parameters themselves are fixed for",
        "the cohort at the 35-mm reference diameter diam_ref. Mujoomdar 2007",
        "reports brain-metastasis frequencies by tumor-size interval; Bilous",
        "2019 drew diameters uniformly within each interval."
      ),
      source_name = "S"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 136L,
    n_studies = 1L,
    disease_state = "lung adenocarcinoma (non-small cell lung cancer) at diagnosis",
    dose_range = "none (natural history before treatment)",
    regions = "not reported in Bilous 2019 (literature series of Mujoomdar 2007)",
    notes = paste(
      "Population-level data are the brain-metastasis frequencies by primary",
      "tumor diameter interval reported by Mujoomdar et al. 2007 (Radiology",
      "242:882) for 136 adenocarcinoma patients, taken from the literature",
      "(Bilous 2019 Materials and Methods, Patient data; Figure 1A)."
    )
  )

  ini({
    ltsmax <- fixed(log(1e12))
    label("Natural log of the Gompertz carrying capacity (log cells)") # Methods: carrying capacity exp(alpha0p/betap) = 10^12 cells
    tdouble <- fixed(201)
    label("Primary tumor volume doubling time at the reference diameter, adenocarcinoma (day)") # Supplementary Table S1 'Average' adenocarcinoma = 201 days; Methods text 201 days
    diam_ref <- fixed(35)
    label("Cohort reference primary tumor diameter for the Gompertz calibration (mm)") # Results: 'an adenocarcinoma with diameter 35 mm (median value from ref.50)'
    cell_density <- fixed(1e6)
    label("Tumor cell density (cells/mm^3)") # Methods: '1 mm3 = 10^6 cells'
    gamma <- fixed(1)
    label("Exponent of the primary tumor size in the dissemination rate gamma (unitless)") # Methods, Population level: 'set gamma = 1'

    # Figure 1A: mu_pop = 1.39e-11 and mu_sigma = 1.24e-9 cell-1 day-1.
    # mu_sigma has the units of mu and the text gives CV > 6000% (here
    # 1.24e-9 / 1.39e-11 = 8921%), so mu_pop and mu_sigma are the mean and
    # SD of the log-normal mu; only this reading reproduces Figure 1A.
    # omega^2 = log(1 + CV^2) = 8.982; median = mean * exp(-omega^2 / 2).
    lmu <- log(1.39e-11) - 0.5 * log(1 + (1.24e-9 / 1.39e-11)^2)
    label("Log of the median cellular metastatic potential mu (log 1/(cell * day))") # Figure 1A legend mu_pop = 1.39e-11 and mu_sigma = 1.24e-9 cell-1 day-1 (mean and SD)
    etalmu ~ 8.982 # log(1 + (1.24e-9 / 1.39e-11)^2) from the Figure 1A legend mu_sigma and mu_pop
  })
  model({
    # Cohort Gompertz parameters: closed-form solution of the Methods
    # doubling-time equation with alpha0p = betap * ltsmax, evaluated at
    # the cohort reference size.
    s_ref <- cell_density * pi / 6 * diam_ref^3
    beta_gomp <- -log((log(2 * s_ref) - ltsmax) / (log(s_ref) - ltsmax)) / tdouble
    alpha_gomp <- beta_gomp * ltsmax

    # Age of the primary tumor when it reaches the diagnosis diameter.
    s_diag <- cell_density * pi / 6 * TUMSZ^3
    td <- -log(1 - log(s_diag) / ltsmax) / beta_gomp

    mu <- exp(lmu + etalmu)

    tumor(0) <- 1
    d/dt(tumor) <- (alpha_gomp - beta_gomp * log(tumor)) * tumor
    # Equations 2 and 4: N(t) = integral of mu * Sp^gamma.
    d/dt(metastases) <- mu * tumor^gamma

    # Brain metastasis present when N > 1; read at t = td for the
    # probability at diagnosis (Methods, P(N(Td) > 1)).
    bm_present <- 0
    if (metastases > 1) {
      bm_present <- 1
    }
  })
}
