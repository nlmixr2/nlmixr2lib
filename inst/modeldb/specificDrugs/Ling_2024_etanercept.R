Ling_2024_etanercept <- function() {
  description <- "One-compartment population PK model for the etanercept biosimilar Benepali with first-order subcutaneous absorption and linear elimination in bDMARD-naive adults with active rheumatoid arthritis (Ling 2024). Apparent (CL/F, V/F) parameterisation; the absorption rate constant was fixed from Korth-Bradley 2000 and carries no between-subject variability, between-subject variability on V/F was removed by the authors for stability, and the additive residual error was fixed at 0.0001. No covariates were retained: age, body weight, sex and concurrent csDMARD therapy were screened and rejected. Unlike the companion adalimumab model from the same paper, this parameter table is internally consistent and reproduces the paper's own Figures 2 and 4."
  reference <- paste(
    "Ling SF, Ogungbenro K, Darwich AS, Ariff ABM, Nair N, Bluett J,",
    "Morgan AW, Isaacs JD, Wilson AG, Hyrich KL, Barton A, Plant D.",
    "Population Pharmacokinetic Analysis and Simulation of Alternative",
    "Dosing Regimens for Biosimilars to Adalimumab and Etanercept in",
    "Patients with Rheumatoid Arthritis. Pharmaceutics. 2024;16(6):702.",
    "doi:10.3390/pharmaceutics16060702.",
    sep = " "
  )
  vignette <- "Ling_2024_tnfi_biosimilars"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "etanercept", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "etanercept", units = "mg", specimen = "serum", verified = TRUE)
  )

  # Screened during covariate model building but NOT retained. Ling 2024
  # Section 2.7(ii) lists the four screened covariates; Section 3.2 states
  # for the etanercept model that "No covariates were included due to the
  # same reasons as for the adalimumab biosimilar model". No point estimate
  # is published for any of them, so none can be carried into model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at treatment initiation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F; not retained. Cohort median 57.5 years (IQR 56-59), Ling 2024 Table 2.",
      source_name        = "age"
    ),
    WT = list(
      description        = "Body weight at treatment initiation",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F; not retained. Cohort median 70.5 kg (IQR 69-84), Ling 2024 Table 2.",
      source_name        = "body weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened; not retained. 4 of 6 subjects were female (Ling 2024 Table 2).",
      source_name        = "sex"
    ),
    CONMED_CSDMARD = list(
      description        = "Concurrent conventional synthetic DMARD therapy indicator, 1 = on a csDMARD, 0 = not",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concurrent csDMARD)",
      notes              = "Screened; not retained. Ling 2024 Table 2 reports 4 (100.00%) with 2 missing observations, i.e. all four subjects with a recorded value were on a csDMARD, leaving the covariate with no contrast in this cohort. Documentation-only label; see the note in Ling_2024_adalimumab.R for why no canonical was registered.",
      source_name        = "concurrent csDMARD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 6L,                                                 # Ling 2024 Section 3.1: six RA patients commencing Benepali
    n_observations = 40L,                                                # Ling 2024 Section 3.2: 40 serum samples available for analysis
    n_studies      = 1L,                                                 # BRAGGSS-PD sub-study, single prospective cohort
    age_range      = "median 57.5 years (IQR 56-59); inclusion required age >= 18 years",  # Ling 2024 Table 2
    age_median     = "57.5 years",                                       # Ling 2024 Table 2
    weight_range   = "median 70.5 kg (IQR 69-84)",                       # Ling 2024 Table 2
    weight_median  = "70.5 kg",                                          # Ling 2024 Table 2
    sex_female_pct = 66.67,                                              # Ling 2024 Table 2: 4 of 6 female
    race_ethnicity = "5 of 6 White, 1 of 6 West African",                # Ling 2024 Section 3.1
    disease_state  = "Rheumatoid arthritis by the 1987 American College of Rheumatology criteria; bDMARD-naive; pre-treatment DAS28 >= 5.1 required for entry (cohort median DAS28 5.33, IQR 4.96-5.58).",
    dose_range     = "50 mg subcutaneously every 7 days (licensed Benepali regimen), self-administered by pre-filled auto-injector; followed for 12 weeks.",
    regions        = "United Kingdom (three rheumatology centres in Greater Manchester)",
    notes          = "Baseline demographics are Ling 2024 Table 2. Real-world NHS patients recruited to the Personalised Dosing sub-study of BRAGGSS (BRAGGSS-PD) between January 2019 and August 2021; recruitment was curtailed by the COVID-19 pandemic below the planned sample size. Serum etanercept was measured by Promonitor-ETN-1DV ELISA (Grifols). Optimal sampling times (baseline, 1 h, 6 days, then 2, 4, 6 and 12 weeks) were designed in PopDes; the 6-day sample is the extra timepoint relative to the adalimumab schedule. The paper notes this is the first published popPK analysis of Benepali and the first at the licensed 50 mg weekly subcutaneous regimen."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, Ling 2024 Table 4. Apparent (CL/F, V/F)
    # values: dosing is subcutaneous and absolute bioavailability is not
    # identifiable.
    #
    # Unlike the companion adalimumab table, Table 4 is internally
    # consistent: these values reproduce the published population
    # simulation (Figure 2) to within 1.5% on the 5th, 50th and 95th
    # percentiles, and the published typical-individual first peak
    # (Figure 4) to within 0.2%.
    # ------------------------------------------------------------------

    # Fixed, with no BSV estimated: Ling 2024 Section 3.2, "Due to a large
    # RSE% from estimation of ka (minimum 1.070%), this value was fixed to
    # 0.0396 hour-1, as per Korth-Bradley et al. [24], and the BSV was not
    # estimated for this parameter."
    lka <- fixed(log(0.0396)); label("First-order SC absorption rate ka (1/h)")  # Ling 2024 Table 4 (Fixed) and Sect. 3.2 text, both 0.0396

    lcl <- log(0.0404); label("Apparent clearance CL/F (L/h)")                    # Ling 2024 Table 4
    lvc <- log(7.76);   label("Apparent central volume of distribution V/F (L)")  # Ling 2024 Table 4

    # ------------------------------------------------------------------
    # Between-subject variability. Table 4 heads this row "(%)" but prints
    # the Monolix log-scale SD un-multiplied (0.173, i.e. about 17.3% CV);
    # the companion adalimumab table prints the same quantity multiplied
    # by 100. Either reading gives essentially the same variance here
    # (0.173 vs 0.1717), and both reproduce the Figure 2 percentile bands.
    #
    # BSV on V/F was estimated and then REMOVED by the authors: Ling 2024
    # Section 3.2, "Additionally, the model had a large RSE% for estimated
    # VD, so the BSV estimation was removed from this parameter." Only a
    # typical value is therefore available for V/F, and no etalvc term
    # appears below. Independent confirmation: the thesis Discussion
    # (reference 23, p. 120) quotes "mean VD (7.76 L)", exactly equal to
    # the typical value, which can only hold if every individual shares it.
    # ------------------------------------------------------------------
    etalcl ~ 0.029929  # Ling 2024 Table 4: omega CL 0.173 log-scale SD; variance 0.173^2

    # ------------------------------------------------------------------
    # Residual unexplained variability, combined additive + proportional
    # per Ling 2024 Section 3.2. The proportional SD is Table 4's 0.46
    # (printed un-multiplied under a "(%)" heading, i.e. 46%); it
    # reproduces the Figure 2 percentile bands to 1.5%.
    #
    # The additive SD was deliberately fixed by the authors and is not
    # tabulated: Ling 2024 Section 3.2, "Finally, the additive error
    # standard deviation (SD) was fixed to 0.0001 to ensure model
    # stability." It is a numerical stabiliser, not an assay-error
    # estimate, and is encoded as fixed() to preserve that provenance.
    # ------------------------------------------------------------------
    propSd <- 0.46;         label("Proportional residual error (fraction)")  # Ling 2024 Table 4: sigma_prop 0.46
    addSd  <- fixed(0.0001); label("Additive residual error (mg/L)")          # Ling 2024 Sect. 3.2 text: "the additive error standard deviation (SD) was fixed to 0.0001"
  })

  model({
    # Individual PK parameters. ka is fixed and V/F carries no BSV, so
    # only cl takes an eta.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    # One-compartment elimination rate constant.
    kel <- cl / vc

    # First-order absorption from the subcutaneous depot into serum.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Doses in mg and vc in L give central / vc directly in mg/L (= ug/mL),
    # the unit used in Ling 2024 Table 4, its figures, and Supplementary
    # Table S2.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
