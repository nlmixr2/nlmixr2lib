Agyeman_2022_covid19_rtcl_neant <- function() {
  description <- paste(
    "Reduced target-cell-limited (rTCL) viral dynamic model of SARS-CoV-2",
    "upper-respiratory viral load after symptom onset, fitted by Agyeman et al.",
    "(2022) to dataset B: 321 hospitalised COVID-19 patients (563",
    "nasopharyngeal samples within 14 days of symptom onset) of the French",
    "COVID cohort analysed by Neant et al. (2021). Two ODE states (Kim et al. 2021",
    "formulation): the fraction of target cells remaining 'target' (f,",
    "starting at 1) is infected by virus at rate beta, and free virus 'virus'",
    "(copies/mL), starting at the viral load at symptom onset, replicates at",
    "maximum rate gamma (fixed to 1/day) and is eliminated at the overall",
    "viral elimination rate delta. Time is days since symptom onset. Typical",
    "values only: the paper estimated log-normal IIV on each estimated",
    "parameter and an additive residual error on the natural-log viral load",
    "but did not report their magnitudes, so the etas and residual SD are",
    "carried as fixed(0)."
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
    target = list(
      analyte = "uninfected respiratory epithelial target cells, as the fraction f of the initial number",
      units = "fraction (unitless)",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "SARS-CoV-2 viral RNA (free virus)",
      units = "copies/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 321L,
    n_observations = 563L,
    n_studies = NA_integer_,
    age_range = "Not reported in Agyeman 2022 (see Neant et al. 2021).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = "COVID-19 (SARS-CoV-2 infection) requiring hospitalisation in conventional wards or intensive care units; 78 deaths during follow-up.",
    dose_range = "No study drug; patients received routine care including antiviral, antibiotic, antifungal or corticosteroid treatment, which the model does not represent.",
    regions = "France (French COVID cohort).",
    notes = paste(
      "Dataset B of Agyeman 2022 (Section 2.1): nasopharyngeal viral load",
      "data from the prospective French COVID cohort of hospitalised patients",
      "analysed by Neant et al. (Proc Natl Acad Sci U S A 2021;118(8)).",
      "Restricted to 0-14 days after symptom onset. Viral loads below the",
      "limit of detection were censored in the fit. Patient demographics are",
      "given in the Neant 2021 source publication, not in Agyeman 2022."
    )
  )

  ini({
    lrbase <- log(1.3e8)
    label("Viral load at symptom onset V(0) (copies/mL)") # Table S3, rTCL column: V(0) = 1.3 x 10^8 copies/ml (RSE 2%)
    lbeta <- log(1.27e-3)
    label("Rate constant for virus infection beta ((copies/mL)^-1 day^-1)") # Table S3, rTCL column: beta = 1.27 x 10^-3 (RSE 14.2%)
    ldelta <- log(0.80)
    label("Overall viral elimination rate delta (1/day)") # Table S3, rTCL column: delta = 0.80 1/day (RSE 16.1%)
    lgamma <- fixed(log(1))
    label("Maximum rate constant for viral replication gamma (1/day)") # Table S3, rTCL column: gamma = 1 (fixed)
    f0 <- fixed(1)
    label("Fraction of target cells remaining at time zero f(0) (unitless)") # Table S3 footnote: f(0) = 1 (for rTCL)

    # Section 2.3: log-normal IIV on every estimated parameter; variances not
    # reported in the paper or Supporting Information, so held at zero.
    etalrbase ~ fixed(0)
    etalbeta ~ fixed(0)
    etaldelta ~ fixed(0)

    # Section 2.3: additive (normal) residual error on the log-transformed viral
    # load; magnitude not reported, so held at zero.
    addSd <- fixed(0)
    label("Additive residual SD on log10 viral load (log10 copies/mL); not reported")
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    beta <- exp(lbeta + etalbeta)
    delta <- exp(ldelta + etaldelta)
    gamma <- exp(lgamma)

    # Equation 2
    d/dt(target) <- -beta * target * virus
    d/dt(virus) <- gamma * target * virus - delta * virus
    target(0) <- f0
    virus(0) <- rbase

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
