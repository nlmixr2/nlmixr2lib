Rao_2021_ethambutol <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption, an absorption lag time and linear elimination for oral ethambutol in Ugandan adults living with HIV and hospitalised with sepsis (with or without tuberculous meningitis) who were starting first-line anti-tuberculosis therapy (Rao 2021). No covariates were retained; the additive residual-error magnitude is not reported and is carried as fixed(0)."
  reference <- "Rao PS, Moore CC, Mbonde AA, Nuwagira E, Orikiriza P, Nyehangane D, Al-Shaer MH, Peloquin CA, Gratz J, Pholwat S, Arinaitwe R, Boum Y, Mwanga-Amumpaire J, Houpt ER, Kagan L, Heysell SK, Muzoora C. Population Pharmacokinetics and Significant Under-Dosing of Anti-Tuberculosis Medications in People with HIV and Critical Illness. Antibiotics (Basel). 2021;10(6):739. doi:10.3390/antibiotics10060739"
  vignette <- "Rao_2021_antituberculosis_critical_illness"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods 4.4) but not retained for ethambutol (Results 2.4: 'None of the covariates screened were significant').",
      source_name = "age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Screened in the stepwise covariate search (Methods 4.4) but not retained for ethambutol (Results 2.4).",
      source_name = "sex"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods 4.4) but not retained for ethambutol (Results 2.4).",
      source_name = "total body weight"
    ),
    MUAC = list(
      description = "Mid-upper arm circumference",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods 4.4) but not retained for ethambutol (Results 2.4). The paper reports MUAC in cm (Table 1: median 22 cm); the register stores it in mm.",
      source_name = "MUAC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "adults over 18 years (enrolled cohort of 81: median 36 years, IQR 30-43.5)",
    age_median = "36 years (enrolled cohort of 81)",
    weight_range = "enrolled cohort of 81: IQR 45.3-58 kg",
    weight_median = "52.5 kg (enrolled cohort of 81)",
    sex_female_pct = 36.3,
    hiv_positive_pct = 100,
    disease_state = "People living with HIV hospitalised with sepsis (two or more SIRS criteria plus suspected infection), with or without suspected tuberculous meningitis (27/81 = 33% concurrent meningitis), starting first-line therapy for presumed drug-susceptible tuberculosis. Of the 49 analysed, 26 had microbiologically confirmed TB and 23 clinical TB. Median CD4 169 cells/mm3; 76% on antiretroviral therapy before admission.",
    dose_range = "Weight-based fixed-dose combination rifampin/isoniazid/pyrazinamide/ethambutol once daily; ethambutol median 18.4 mg/kg (IQR 16.5-19.6; printed as mg/L in Results 2.4). PK sampled at steady state two weeks after treatment start.",
    regions = "Uganda (Mbarara Regional Referral Hospital).",
    notes = "81 enrolled; 49 completed PK testing (18 died, 13 lost to follow-up or withdrew, 1 incomplete). Sparse serum sampling at 1, 2, 4 and 6 h post-dose; 173 ethambutol concentrations. LC-MS/MS assay validated over 0.2-10 ug/mL. Demographics are reported only for the 81 enrolled (Table 1), not for the 49 analysed. The TB diagnostic group (microbiologically confirmed versus clinical TB) was also screened (Methods 4.4) and not retained. Phoenix NLME, FOCE-ELS. ClinicalTrials.gov NCT03559582."
  )

  ini({
    # Rao 2021 Table 3, ethambutol rows: 'Estimate (%CV)' column. All values
    # are apparent oral parameters (no IV reference; F not estimated).
    lka <- log(0.15); label("Absorption rate constant ka (1/h)") # Table 3: Ka = 0.15 1/h (21.8 %CV)
    lvc <- log(75.17); label("Apparent volume of distribution V/F (L)") # Table 3: V = 75.17 L (23.7 %CV)
    lcl <- log(51.6); label("Apparent clearance CL/F (L/h)") # Table 3: Cl = 51.6 L/h (12.2 %CV)
    ltlag <- log(0.4); label("Absorption lag time (h)") # Table 3: Tlag = 0.4 h (12.3 %CV)

    # IIV: Table 3 'IIV (Shrinkage)' column, read as Phoenix NLME Omega
    # diagonal VARIANCES of exponential etas (Methods 4.4: 'Interindividual
    # variability was included using exponential function'); the number in
    # parentheses is the eta shrinkage fraction. See the vignette for the
    # Table 4 spread check that supports the variance reading.
    etalka ~ 0.6 # Table 3: Ka IIV 0.6 (shrinkage 0.4)
    etalvc ~ 0.4 # Table 3: V IIV 0.4 (shrinkage 0.1)
    etalcl ~ 0.3 # Table 3: Cl IIV 0.3 (shrinkage 0.3)
    etaltlag ~ 0.6 # Table 3: Tlag IIV 0.6 (shrinkage 0.4)

    # Additive residual error (Results 2.4: 'additive residual error was
    # selected'); its magnitude is not reported anywhere in the paper, which
    # has no supplement. Carried as zero rather than invented.
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- magnitude not reported in the source)")
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
