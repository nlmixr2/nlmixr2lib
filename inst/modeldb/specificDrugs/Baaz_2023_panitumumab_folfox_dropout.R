Baaz_2023_panitumumab_folfox_dropout <- function() {
  description <- paste(
    "Weibull time-to-dropout model for patients with RAS exon 2 wild-type",
    "metastatic colorectal cancer on FOLFOX with or without panitumumab",
    "(Baaz 2023, PRIME study, n = 248). Dropout (censoring without target or",
    "non-target progression) is independent of tumor response and of",
    "treatment; the hazard is h = (alfa_haz / te_haz) *",
    "(t / te_haz)^(alfa_haz - 1) with shape alfa_haz and characteristic",
    "time te_haz. Time is in months",
    "from the start of treatment. No random effects. Calibrated separately",
    "from the joint TGI / non-target-progression model",
    "modellib('Baaz_2023_panitumumab_folfox'); the paper combines the two to",
    "predict progression-free survival.",
    sep = " "
  )
  reference <- paste(
    "Baaz M, Cardilin T, Jirstrand M. Model-based prediction of",
    "progression-free survival for combination therapies in oncology.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1227-1237.",
    "doi:10.1002/psp4.13003.",
    "Parameter estimates from Table 1; model code from the Monolix project in",
    "Data S1.",
    sep = " "
  )
  vignette <- "Baaz_2023_panitumumab_folfox"

  units <- list(
    time = "month",
    dosing = "n/a (no dose events)",
    concentration = "n/a (the output sur is the probability of not having dropped out)"
  )

  covariateData <- list()

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of dropout (uninformative censoring)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 248L,
    n_studies = 1L,
    disease_state = "metastatic colorectal cancer, RAS exon 2 wild-type, first-line treatment",
    dose_range = "FOLFOX4 every 2 weeks, alone or with panitumumab 6 mg/kg every 2 weeks",
    notes = paste(
      "Same 127 FOLFOX and 121 FOLFOX + panitumumab PRIME patients as",
      "Baaz_2023_panitumumab_folfox. Dropout is a patient without target or",
      "non-target progression at the final check-up (Baaz 2023 Methods);",
      "the Weibull form fitted better than exponential and Gompertz."
    )
  )

  ini({
    # Baaz 2023 Table 1, TTE section; Data S1 Monolix code
    # 'h = p/Te * (t/Te)^(p-1)'. No BSV was estimated.
    lalfa_haz <- log(1.50)
    label("Weibull dropout shape k (unitless)") # Table 1: k = 1.50 (RSE 11%); Monolix p
    lte_haz <- log(36)
    label("Weibull dropout scale lambda, characteristic time (month)") # Table 1: lambda = 36 month (RSE 14%); Monolix Te
  })
  model({
    alfa_haz <- exp(lalfa_haz)
    te_haz <- exp(lte_haz)

    # Weibull hazard h2(t) = (k/lambda) * (t/lambda)^(k-1) (Methods).
    haz_drop <- alfa_haz / te_haz * (t / te_haz)^(alfa_haz - 1)
    d/dt(cumhaz) <- haz_drop
    sur <- exp(-cumhaz)
  })
}
