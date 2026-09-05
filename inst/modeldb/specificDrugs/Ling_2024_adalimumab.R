Ling_2024_adalimumab <- function() {
  description <- "One-compartment population PK model for the adalimumab biosimilar Amgevita with first-order subcutaneous absorption and linear elimination in bDMARD-naive adults with active rheumatoid arthritis (Ling 2024). Apparent (CL/F, V/F) parameterisation; the absorption rate constant was fixed from Ternant 2015 and carries no between-subject variability. No covariates were retained: age, body weight, sex and concurrent csDMARD therapy were screened and rejected. IMPORTANT: the source Table 3 misprints both ka and CL/F; the values used here come from the paper's own Methods text and from the first author's thesis (reference 23), and are confirmed against the paper's Figures 1 and 3 and Supplementary Figures S1 and S3. See the vignette Errata."
  reference <- paste(
    "Ling SF, Ogungbenro K, Darwich AS, Ariff ABM, Nair N, Bluett J,",
    "Morgan AW, Isaacs JD, Wilson AG, Hyrich KL, Barton A, Plant D.",
    "Population Pharmacokinetic Analysis and Simulation of Alternative",
    "Dosing Regimens for Biosimilars to Adalimumab and Etanercept in",
    "Patients with Rheumatoid Arthritis. Pharmaceutics. 2024;16(6):702.",
    "doi:10.3390/pharmaceutics16060702.",
    "Corrected CL/F sourced from the paper's own reference 23:",
    "Ling S. Personalising Dosing of Biologic Therapies in Inflammatory",
    "Arthritis to Maximise Cost-Benefit. PhD Thesis, University of",
    "Manchester, 2022, Section 4.4 (Discussion), p. 115.",
    sep = " "
  )
  vignette <- "Ling_2024_tnfi_biosimilars"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "adalimumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "adalimumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  # Screened during covariate model building but NOT retained in the final
  # model. Ling 2024 Section 2.7(ii): "The following four covariates were
  # investigated: age and body weight (continuous covariates) and sex and
  # concurrent conventional synthetic disease-modifying anti-rheumatic drug
  # (csDMARD) therapy (binary covariates)." Section 3.2: "No covariates were
  # included, as covariates tested only demonstrated a modest improvement in
  # -2LL and AIC whilst complicating the sparse-sampling model with redundant
  # variables." No point estimate is published for any of them, so none can be
  # carried into model(); they are recorded here for provenance only.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at treatment initiation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F; not retained. Cohort median 50.5 years (IQR 46-61), Ling 2024 Table 1.",
      source_name        = "age"
    ),
    WT = list(
      description        = "Body weight at treatment initiation",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and V/F; not retained. Cohort median 85.5 kg (IQR 66-111), Ling 2024 Table 1. The Discussion notes this adalimumab-biosimilar cohort was heavier than the etanercept-biosimilar cohort (85.5 vs 70.5 kg) and speculates that increased adiposity may have impaired subcutaneous absorption, but no weight effect was estimable in n = 10.",
      source_name        = "body weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened; not retained. 9 of 10 subjects were female (Ling 2024 Table 1). The Discussion states explicitly that 'due to the skew between female and male participants, we cannot distinguish between BSV and the effect of sex', so this covariate is unidentifiable in this cohort rather than merely non-significant.",
      source_name        = "sex"
    ),
    CONMED_CSDMARD = list(
      description        = "Concurrent conventional synthetic DMARD therapy indicator, 1 = on a csDMARD, 0 = not",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concurrent csDMARD)",
      notes              = "Screened; not retained. 8 of 10 subjects (80%) were on a concurrent csDMARD (Ling 2024 Table 1); the paper does not name the individual agents. Documentation-only label: because the covariate was rejected and no coefficient is published, no entry was added to inst/references/covariate-columns.md. A future extraction that RETAINS a csDMARD-class indicator should register the canonical then; it would be a class-level member of the existing CONMED_* family (cf. CONMED_ABX, CONMED_AED, CONMED_IMMUNOMOD).",
      source_name        = "concurrent csDMARD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,                                                # Ling 2024 Section 3.1: ten RA patients commencing Amgevita
    n_observations = 58L,                                                # Ling 2024 Section 3.2: 58 serum samples available for analysis
    n_studies      = 1L,                                                 # BRAGGSS-PD sub-study, single prospective cohort
    age_range      = "median 50.5 years (IQR 46-61); inclusion required age >= 18 years",  # Ling 2024 Table 1
    age_median     = "50.5 years",                                       # Ling 2024 Table 1
    weight_range   = "median 85.5 kg (IQR 66-111)",                      # Ling 2024 Table 1
    weight_median  = "85.5 kg",                                          # Ling 2024 Table 1
    sex_female_pct = 90,                                                 # Ling 2024 Table 1: 9 of 10 female
    race_ethnicity = "100% White",                                       # Ling 2024 Section 3.1: "All patients were white."
    disease_state  = "Rheumatoid arthritis by the 1987 American College of Rheumatology criteria; bDMARD-naive; pre-treatment DAS28 >= 5.1 required for entry (cohort median DAS28 5.71, IQR 5.20-6.09).",
    dose_range     = "40 mg subcutaneously every 14 days (licensed Amgevita regimen), self-administered by pre-filled auto-injector; followed for 12 weeks.",
    regions        = "United Kingdom (three rheumatology centres in Greater Manchester)",
    notes          = "Baseline demographics are Ling 2024 Table 1. Real-world NHS patients recruited to the Personalised Dosing sub-study of BRAGGSS (BRAGGSS-PD) between January 2019 and August 2021; recruitment was curtailed by the COVID-19 pandemic below the planned sample size. Serum adalimumab was measured by Promonitor-ADL-1DV ELISA (Grifols). Optimal sampling times (baseline, 1 h, then 2, 4, 6 and 12 weeks, all pre-dose troughs after the second sample) were designed in PopDes. Samples were NOT tested for anti-drug antibodies, which the paper lists as a limitation; the model therefore carries no immunogenicity term. Doses self-administered between study visits were recorded at their nominal times."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Apparent (CL/F, V/F) values: dosing is
    # subcutaneous and absolute bioavailability is not identifiable.
    #
    # TWO OF THE THREE STRUCTURAL VALUES IN LING 2024 TABLE 3 ARE
    # MISPRINTS. Both corrections are taken from printed prose (not from
    # figure-fitting) and are then confirmed numerically against four
    # independent outputs the paper itself published. See the vignette
    # "Source trace" and "Errata" sections for the full audit.
    # ------------------------------------------------------------------

    # Table 3 prints ka = 0.1167 /h. The paper's own Methods text prints
    # 0.01167 /h and attributes it to Ternant 2015 (reference 10); the
    # thesis (reference 23, p. 94) prints 0.01167 /h in text and 0.1167 /h
    # in its table, i.e. the same internal conflict. 0.01167 /h = 0.28 /day
    # is Ternant's published value. Fixed, and no BSV was estimated on it:
    # Ling 2024 Section 3.2, "Due to a large percentage relative standard
    # error (RSE%) from estimation of ka (minimum 622%), this value was
    # fixed to 0.01167 hour-1, as per Ternant et al. [10], and the random
    # effect (BSV) was not estimated for this parameter."
    lka <- fixed(log(0.01167)); label("First-order SC absorption rate ka (1/h)")  # Ling 2024 Sect. 3.2 text (NOT Table 3, which misprints 0.1167)

    # Table 3 prints CL/F = 0.00283 L/h. The thesis Discussion (reference
    # 23, p. 115) prints "CL was estimated at 0.0121 L/hr, which falls
    # within the range determined from Humira trial data of 0.00676 -
    # 0.0322 L/hr" -- and 0.00283 lies BELOW that range, contradicting the
    # sentence that contains it. 0.0121 L/h is the value used here.
    lcl <- log(0.0121); label("Apparent clearance CL/F (L/h)")  # Ling 2022 thesis (Ling 2024 ref. 23) Discussion p. 115 (NOT Ling 2024 Table 3, which misprints 0.00283)

    lvc <- log(9.19); label("Apparent central volume of distribution V/F (L)")  # Ling 2024 Table 3 (confirmed by thesis p. 115: "The typical individual VD was estimated at 9.19 L")

    # ------------------------------------------------------------------
    # Between-subject variability. Ling 2024 Table 3 heads these rows
    # "(%)" and calls them CVs, but they are Monolix's log-scale SDs
    # multiplied by 100: reproducing the Figure 1 5th/50th/95th percentile
    # bands requires omega_CL = 0.689 and omega_V = 0.156 read as log-scale
    # SDs (0.5% error on all three percentiles) rather than as CVs
    # (5.7% error). The companion etanercept table prints the same
    # quantities un-multiplied (0.173) under the same "(%)" heading,
    # which is the clerical tell.
    # ------------------------------------------------------------------
    etalcl ~ 0.474721  # Ling 2024 Table 3: omega CL 68.90 -> log-scale SD 0.689; variance 0.689^2
    etalvc ~ 0.024336  # Ling 2024 Table 3: omega VD 15.60 -> log-scale SD 0.156; variance 0.156^2

    # ------------------------------------------------------------------
    # Residual unexplained variability. Ling 2024 Section 3.2 states a
    # combined additive + proportional model was used.
    #
    # The proportional term is sound: Table 3's 26.00 is 0.26 on the same
    # x100 convention as the omegas above, and 0.26 reproduces the Figure 1
    # percentile bands.
    #
    # The additive term is NOT usable. Table 3 prints 10.80 mg/L, which is
    # larger than every adalimumab concentration in the study except one
    # (Supplementary Table S1 maximum 14.477 mg/L). Two published outputs
    # falsify it: (1) the IWRES distribution (Supplementary Figure S2)
    # spans about +/-2.5 with SD near 1, whereas an additive SD of 10.8
    # would compress every IWRES into about +/-0.3; and (2) re-simulating
    # Figure 1 with it drives the 5th percentile to about -13 mg/L, while
    # the published 5th percentile sits at 2.78 mg/L. Proportional-only
    # error reproduces the published bands to 0.5%.
    #
    # No corrected additive value is printed in the paper, its supplement,
    # or the thesis. Per the standing policy for unreported/unusable RUV
    # (fixed(0) plus erratum; never invent a variance) the additive
    # component is fixed at zero, which is also what the paper's own
    # Figure 1 behaves as if it used.
    # ------------------------------------------------------------------
    propSd <- 0.26;    label("Proportional residual error (fraction)")  # Ling 2024 Table 3: sigma_prop 26.00 -> 0.26
    addSd  <- fixed(0); label("Additive residual error (mg/L)")          # Ling 2024 Table 3 prints 10.80 mg/L; falsified by Suppl. Fig. S2 and Fig. 1 - see Errata
  })

  model({
    # Individual PK parameters. ka carries no BSV (fixed, not estimated).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # One-compartment elimination rate constant.
    kel <- cl / vc

    # First-order absorption from the subcutaneous depot into serum.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Doses in mg and vc in L give central / vc directly in mg/L (= ug/mL),
    # the unit used in Ling 2024 Table 3, its figures, and Supplementary
    # Table S1.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
