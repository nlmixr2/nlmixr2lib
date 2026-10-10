Lin_2021_asparaginaseErwiniaRecombinant <- function() {
  description <- "One-compartment population PK model for recombinant Erwinia chrysanthemi asparaginase (JZP-458, marketed as Rylaze) given intramuscularly or as a 2-hour intravenous infusion to healthy adults (Lin 2021, phase 1 study JZP458-101). The measured quantity is serum asparaginase activity (SAA), so all amounts are activity units (IU) rather than mass. Intravenous doses enter the central compartment directly. Intramuscular doses use sequential mixed-order absorption: the bioavailable fraction F enters the depot as a zero-order input at a constant rate R1 while first-order absorption ka drains the depot into the central compartment, which makes the terminal phase absorption rate limited (flip-flop) after intramuscular dosing. Body weight is an allometric (power) covariate on clearance; interindividual variability is exponential on clearance and volume, with a proportional residual error."
  reference <- "Lin T, Dumas T, Kaullen J, Berry NS, Choi MR, Zomorodi K, Silverman JA. Population pharmacokinetic model development and simulation for recombinant Erwinia asparaginase produced in Pseudomonas fluorescens (JZP-458). Clin Pharmacol Drug Dev. 2021;10(12):1503-1513. doi:10.1002/cpdd.1002"
  vignette <- "Lin_2021_asparaginaseErwiniaRecombinant"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  # Amounts are asparaginase ACTIVITY units (IU), not mass: the assay reads
  # serum asparaginase activity in IU/mL and the zero-order absorption rate is
  # reported in IU/h (Lin 2021 Table 2).
  compartmentData <- list(
    depot = list(
      analyte = "recombinant Erwinia chrysanthemi asparaginase (JZP-458)",
      units = "IU",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "recombinant Erwinia chrysanthemi asparaginase (JZP-458)",
      units = "IU",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric (power) covariate on clearance only, normalized to 70 kg: CL (mL/h) = 146 * (WT/70)^0.863 (Lin 2021 Results, Covariate analysis and final covariate population PK model selection; Table 2). Volume carries no body-size covariate. Single-dose study, so weight is time-fixed per subject. Analysis-set weight 78.3 +/- 9.6 kg (mean +/- SD; Lin 2021 Table 1). Weight was chosen over BSA, which was also statistically significant, because of established allometric scaling of CL on body weight.",
      source_name = "WT (WTKG in the supplementary figures)"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Statistically significant on CL (3.4% of CL variability explained) but not retained; body weight was selected instead (Lin 2021 Results, Covariate analysis). BSA is still needed to compute the mg/m^2 clinical dose. Analysis set 1.9 +/- 0.1 m^2 (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and Vd but not retained: no trend after including weight (Lin 2021 Results, Covariate analysis). Analysis set 38.3 +/- 8.6 years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and Vd but not retained: no trend after including weight (Lin 2021 Results, Covariate analysis). 7 of 24 participants female (Table 1: 17 male, 71%)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race listed among the screened intrinsic covariates, but 'Race and ethnicity could not be evaluated due to the small participant numbers in each of these subgroups' (Lin 2021 Results, Covariate analysis). 4 of 24 participants Black/African American (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 331L,
    age_range = "18-55 years (eligibility); mean 38.3 +/- 8.6 years",
    weight_range = "mean 78.3 +/- 9.6 kg",
    bsa_range = "mean 1.9 +/- 0.1 m^2",
    sex_female_pct = 29.2,
    race_ethnicity = c(
      `White` = 83,
      `Black/African American` = 17,
      `Hispanic/Latino` = 96
    ),
    disease_state = "Healthy adult volunteers (body mass index 19.0-30.0 kg/m^2).",
    dose_range = "Single dose of JZP-458: 12.5 or 25 mg/m^2 intramuscular (n = 6 each; dorsogluteal or deltoid, at most 2 mL per injection site) or 25 or 37.5 mg/m^2 as a 2-hour intravenous infusion (n = 6 each).",
    regions = "United States (single center, Miami, Florida)",
    notes = "Lin 2021 Methods (Study design) and Table 1. Phase 1, randomized, single-center, open-label study JZP458-101 (November 2018 to May 2019); six further participants received Erwinia chrysanthemi asparaginase (ERW) and were not part of the PopPK analysis. 331 quantifiable SAA observations, intensive sampling to 96 h post dose. Assay lower limit of quantitation 0.025 IU/mL. Race percentages are of 24; ethnicity was self-reported (Table 1)."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural PK -- Lin 2021 Table 2 (Population pharmacokinetic
    # parameters of JZP-458 following IV and IM administration).
    # Reference subject: 70 kg body weight.
    #
    # AMOUNT UNITS ARE ACTIVITY UNITS (IU), NOT MASS. The zero-order
    # absorption rate is reported in IU/h and the observation is serum
    # asparaginase ACTIVITY in IU/mL, so central/vc is IU/mL only if the dose
    # is entered in IU. The paper does not state the mg-to-IU specific
    # activity of JZP-458, so a clinical mg/m^2 dose cannot be converted to
    # model units from the source -- see the vignette Errata.
    # -----------------------------------------------------------------
    lcl <- log(146)          ; label("Clearance CL for a 70 kg adult (mL/h)")             # Lin 2021 Table 2 (CL 146 x (WT/70)^0.863 mL/h, 95% CI 128.4-163.6, RSE 6.15%); Results 'CL (mL/h) = 146 (mL/h) x (weight [kg]/70)^0.863'
    lvc <- log(3030)         ; label("Central volume of distribution Vd (mL)")            # Lin 2021 Table 2 (Vd 3030 mL, 95% CI 2655-3405, RSE 6.32%); footnote 'IV, Vd = 3.03 L'
    lka <- log(0.0348)       ; label("First-order absorption rate constant ka (1/h)")     # Lin 2021 Table 2 (ka 0.0348 1/h, 95% CI 0.02942-0.04018, RSE 7.89%)
    lr1 <- log(4000)         ; label("Zero-order input rate R1 into depot (IU/h)")        # Lin 2021 Table 2 (Zero-order absorption 4000 IU/h, 95% CI 1569-6431, RSE 31.01%)
    lfdepot <- log(0.365)    ; label("Intramuscular bioavailability F relative to IV (fraction)") # Lin 2021 Table 2 (F 0.365, 95% CI 0.3074-0.4226, RSE 8.05%); Results 'bioavailability at 36.5%'

    # Allometric weight exponent on CL. Table 2 prints it only inside the CL
    # row, with no separate RSE or CI, but it is not a conventional fixed value
    # (0.75 / 1) and Methods state body weight was screened and retained by
    # the OFV-drop criterion, so it is encoded as estimated.
    e_wt_cl <- 0.863         ; label("Allometric exponent of body weight on CL (unitless)") # Lin 2021 Table 2 (CL row '146 x (WT/70)^0.863'); Results final CL equation

    # -----------------------------------------------------------------
    # Interindividual variability -- Lin 2021 Table 2 'BSV%' column; footnote
    # 'BSV was modeled as exponential'. The column is read as a CV and
    # converted with omega^2 = log(1 + CV^2):
    #   CL: log(1 + 0.1888^2) = 0.03502
    #   Vd: log(1 + 0.3206^2) = 0.09783
    # (reading it as 100*omega instead gives 0.03565 and 0.1028; the two
    # readings differ by under 5% at these small magnitudes). No CL-Vd
    # covariance is reported, so the etas are uncorrelated.
    # -----------------------------------------------------------------
    etalcl ~ 0.03502         # Lin 2021 Table 2 (CL BSV 18.88%)
    etalvc ~ 0.09783         # Lin 2021 Table 2 (Vd BSV 32.06%)

    # -----------------------------------------------------------------
    # Residual error -- Lin 2021 Table 2 'Error model proportional 20.6%'.
    # -----------------------------------------------------------------
    propSd <- 0.206          ; label("Proportional residual error SD (fraction)")         # Lin 2021 Table 2 (Error model proportional 20.6%)
  })

  model({
    # 1. Individual parameters. Reference subject is 70 kg.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)
    r1 <- exp(lr1)
    fdepot <- exp(lfdepot)
    kel <- cl / vc

    # 2. One-compartment disposition. Intravenous doses (2-hour infusion in
    #    the source study) go straight into `central`. Intramuscular doses go
    #    into `depot`, from which first-order absorption ka feeds `central`
    #    (Lin 2021 Results, Base model: 'a sequential mixed order absorption
    #    function was used to estimate both zero- and first-order absorption
    #    rate constant parameters').
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 3. Intramuscular bioavailability and the zero-order depot input, giving
    #    the NONMEM R1 semantics: the input duration is fdepot * amt / r1.
    #    Intramuscular dose records MUST carry `rate = -1`; with the default
    #    `rate = 0` rxode2 silently ignores the modelled rate and gives an
    #    instantaneous bolus into the depot. Intravenous records into
    #    `central` carry their own infusion rate (or `dur = 2`) and are not
    #    affected by fdepot.
    f(depot) <- fdepot
    rate(depot) <- r1

    # 4. Observation: serum asparaginase activity (SAA) in IU/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
