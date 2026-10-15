Agyeman_2022_covid19_si_gastine <- function() {
  description <- paste(
    "Slope-intercept exponential decay (SI) viral dynamic model of SARS-CoV-2",
    "upper-respiratory viral load after symptom onset, fitted by Agyeman et al.",
    "(2022) to dataset A: 252 untreated COVID-19 patients (747 samples within",
    "14 days of symptom onset) compiled by the Gastine et al. (2021)",
    "patient-level meta-analysis. One ODE state, free virus 'virus'",
    "(copies/mL), starting at the viral load at symptom onset and declining",
    "first-order at the overall viral elimination rate delta (1/day). Time is",
    "days since symptom onset. Typical values only: the paper estimated",
    "log-normal IIV on each estimated parameter and an additive residual error",
    "on the natural-log viral load but did not report their magnitudes, so the",
    "etas and residual SD are carried as fixed(0)."
  )
  reference <- paste(
    "Agyeman AA, You T, Chan PLS, Lonsdale DO, Hadjichrysanthou C, Mahungu T,",
    "Wey EQ, Lowe DM, Lipman MCI, Breuer J, Kloprogge F, Standing JF.",
    "Comparative assessment of viral dynamic models for SARS-CoV-2 for",
    "pharmacodynamic assessment in early treatment trials.",
    "Br J Clin Pharmacol. 2022;88(12):5428-5433. doi:10.1111/bcp.15518"
  )
  vignette <- "Agyeman_2022_covid19_viral_dynamics"

  units <- list(
    time = "day",
    dosing = "no dosing (untreated natural-history viral dynamics)",
    concentration = "copies/mL (observed as log10 copies/mL)"
  )

  covariateData <- list()

  compartmentData <- list(
    virus = list(
      analyte = "SARS-CoV-2 viral RNA (free virus)",
      units = "copies/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 252L,
    n_observations = 747L,
    n_studies = NA_integer_,
    age_range = "Not reported in Agyeman 2022 (see Gastine et al. 2021).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = "COVID-19 (SARS-CoV-2 infection), predominantly mild disease with one reported death; untreated patients only.",
    dose_range = "No antiviral treatment (untreated patients selected from the dataset).",
    regions = "Multinational (patient-level data from published case series).",
    notes = paste(
      "Dataset A of Agyeman 2022 (Section 2.1): patient-level viral load data",
      "compiled by the systematic review and meta-analysis of Gastine et al.",
      "(Clin Pharmacol Ther 2021;110:321-333). Restricted to untreated",
      "patients, upper-respiratory-tract sampling sites and 0-14 days after",
      "symptom onset. Viral loads below the limit of detection were censored",
      "in the fit. Patient demographics are given in the Gastine 2021 source",
      "publication, not in Agyeman 2022."
    )
  )

  ini({
    lrbase <- log(1.8e7)
    label("Viral load at symptom onset V(0) (copies/mL)") # Table S2, SI column: V(0) = 1.8 x 10^7 copies/ml (RSE 2.3%)
    ldelta <- log(0.56)
    label("Overall viral elimination rate delta (1/day)") # Table S2, SI column: delta = 0.56 1/day (RSE 13.6%)

    # Section 2.3: log-normal IIV on every estimated parameter; variances not
    # reported in the paper or Supporting Information, so held at zero.
    etalrbase ~ fixed(0)
    etaldelta ~ fixed(0)

    # Section 2.3: additive (normal) residual error on the log-transformed viral
    # load; magnitude not reported, so held at zero.
    addSd <- fixed(0)
    label("Additive residual SD on log10 viral load (log10 copies/mL); not reported")
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    delta <- exp(ldelta + etaldelta)

    # Equation 1
    d/dt(virus) <- -delta * virus
    virus(0) <- rbase

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
