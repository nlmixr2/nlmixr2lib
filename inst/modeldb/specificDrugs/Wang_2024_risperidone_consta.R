Wang_2024_risperidone_consta <- function() {
  description <- "Population PK model for the risperidone active moiety (risperidone + 9-OH-risperidone) after intramuscular gluteal injection of RISPERDAL CONSTA, the reference biweekly long-acting injectable risperidone formulation, in adults with schizophrenia or schizoaffective disorder (Wang 2024). Same structure as the companion Rykindo model: one-compartment disposition fed by three parallel release pathways out of the injection site, namely an immediate zero-order release of fraction F2 into the central compartment over duration D2, a middle release of fraction F3 = 1 - F1 - F2 entering a second depot as a zero-order input of duration D3 after lag ALAG3 and then absorbed first-order with rate K32, and a main first-order release of fraction F1 beginning ALAG1 = 27 days after injection. The long main-release lag reproduces the well-known 3-week delay that obliges 3 weeks of oral risperidone supplementation when Consta is started. Because the active moiety displays flip-flop kinetics, the elimination rate constant K was set equal to KA, so the apparent central volume is derived as V = CL/KA. Apparent clearance is a single value across trials and sexes; only the main-release absorption rate constant differs between the two trials. The companion Rykindo model is modellib('Wang_2024_risperidone_rykindo')."
  reference <- paste(
    "Wang W, Wang X, Dong Y, Walling DP, Liu P, Liu W, Shi Y, Sun K.",
    "Population Pharmacokinetic Analysis to Support and Facilitate Switching",
    "from Risperidone Formulations to Rykindo in Patients with Schizophrenia.",
    "Neurol Ther. 2024;13(2):355-372. doi:10.1007/s40120-024-00578-w.",
    sep = " "
  )
  vignette <- "Wang_2024_risperidone_lai_switching"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    STUDY_CT_USA_102 = list(
      description        = "1 = subject enrolled in trial CT-USA-102 (NCT02091388), the multiple-dose relative-bioavailability trial in which either Rykindo or Consta 25 mg was given every 2 weeks for five injections; 0 = subject enrolled in CT-USA-104 (NCT02186769), the single-dose relative-bioavailability trial",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the single-dose trial CT-USA-104; Consta was not studied in CT-1S01)",
      notes              = "Time-fixed; set from the trial identifier. Switches the first-order absorption rate constant KA of the main release from 0.179 to 0.271 /day (Wang 2024 Table 2, Consta column, rows 'KA_CT-1S01 and CT-USA-104' -- footnote: KA_CT-USA-104 for Consta -- and 'KA_CT-USA-102'). This is the only covariate retained in the Consta model: Wang 2024 Results states 'The only covariate that is statistically significant is found in the studies on the primary absorption rate constant', and no sex effect on clearance was significant. The authors attribute the faster apparent absorption in the multiple-dose trial to its sparser PK sampling scheme rather than to a formulation difference (Wang 2024 Discussion). Because the model sets the elimination rate constant equal to KA (flip-flop kinetics), this indicator also shifts the apparent central volume V = CL/KA. Shared with the companion Rykindo model. Follows the auto-approved STUDY_<id> canonical family.",
      source_name        = "study"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Wang 2024 Table 2 footnote states that the single Consta clearance estimate of 108.0 L/day applies to both males and females, and the Results section reports that no gender effect on clearance was statistically significant for Consta. A sex effect on clearance WAS retained for the companion Rykindo model."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as part of the demographic covariate analysis (Wang 2024 Methods, 'Population PK Model Development') but not retained in the final Consta model."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as part of the demographic covariate analysis (Wang 2024 Methods) but not retained in the final Consta model."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as part of the demographic covariate analysis (Wang 2024 Methods) but not retained in the final Consta model."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "administration site", verified = TRUE),
    depot2  = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 66L,
    n_studies      = 2L,
    n_observations = "1766 plasma concentration records of risperidone and 9-OH-risperidone (Wang 2024 Results, 'Population PK model of Consta active moiety').",
    age_range      = "medians 53.0 / 55.0 years and means 49.0 (SD 10.9) / 46.6 (9.6) years in CT-USA-104 / CT-USA-102 (Wang 2024 Table 1)",
    weight_range   = "medians 79.3 / 88.0 kg and means 81.5 (SD 20.1) / 89.7 (15.4) kg in CT-USA-104 / CT-USA-102 (Wang 2024 Table 1)",
    bmi_range      = "medians 25.35 / 29.2 kg/m^2 and means 26.5 (SD 5.8) / 29.3 (4.3) kg/m^2 (Wang 2024 Table 1)",
    sex_female_pct = 26.1,
    race_ethnicity = "Not reported in Wang 2024. Race was screened as a covariate but no race effect was retained in the final model.",
    disease_state  = "Stable adults with schizophrenia or schizoaffective disorder",
    dose_range     = "25 and 50 mg single intramuscular gluteal injection (CT-USA-104); 25 mg every 2 weeks for five injections (CT-USA-102)",
    regions        = "United States",
    cyp2d6_status  = "Across the two Consta cohorts (Wang 2024 Table 1): 5 poor, 24 intermediate and 37 extensive metabolizers, 3 unknown. CYP2D6 genotype was not retained as a covariate because the active moiety (risperidone + 9-OH-risperidone) is insensitive to CYP2D6 status.",
    notes          = "Pooled from the Consta arms of two phase 1 trials: CT-USA-104 (NCT02186769) and CT-USA-102 (NCT02091388). Wang 2024 Table 1 lists 15 + 54 = 69 enrolled subjects in the two Consta cohorts, while the Results text states that 66 subjects contributed the 1766 PK records; the difference is not explained in the paper. Age, weight, body mass index, sex and race were screened as covariates; only trial (on the main-release absorption rate constant) was retained."
  )

  ini({
    # ==================================================================
    # Wang 2024 Table 2, Consta column (final model estimates; standard
    # errors from 200 bootstrap replicates because the covariance step
    # failed). All values are for the active moiety.
    #
    # Structure (Wang 2024 Fig. 2b) is identical to the Rykindo model:
    # F1 is released first-order from `depot` after lag ALAG1; F2 enters
    # `central` directly as a zero-order input of duration D2; the
    # remainder F3 = 1 - F1 - F2 enters `depot2` as a zero-order input
    # of duration D3 starting at ALAG3 and is then absorbed first-order
    # with rate K32.
    # ==================================================================

    # --- Disposition -------------------------------------------------
    lcl <- log(108.0); label("Apparent clearance CL/F of the active moiety (L/day)")  # Wang 2024 Table 2 Consta: CL = 108.0 L/day (%RSE 4.15); table footnote: applies to both males and females

    # --- Main (third) release ----------------------------------------
    lka <- log(0.179); label("First-order absorption rate constant KA of the main release in trial CT-USA-104 (1/day)")  # Wang 2024 Table 2 Consta: KA = 0.179 /day (%RSE 3.45); table footnote: KA_CT-USA-104 for Consta
    e_study102_ka <- log(0.271 / 0.179); label("Log-multiplicative effect of trial CT-USA-102 on the main-release absorption rate constant (unitless)")  # Derived from Wang 2024 Table 2 Consta: KA_CT-USA-102 / KA_CT-USA-104 = 0.271 / 0.179 = 1.514 (%RSE 9.87 on the 0.271 estimate)
    ltlag_main <- log(27.0); label("Lag time ALAG1 of the main release (day)")  # Wang 2024 Table 2 Consta: ALAG1 = 27.0 day (%RSE 0.871)

    # --- Immediate release -------------------------------------------
    ld2 <- log(0.0467); label("Zero-order duration D2 of the immediate release into the central compartment (day)")  # Wang 2024 Table 2 Consta: D2 = 0.0467 day (%RSE 31.9)

    # --- Middle release ----------------------------------------------
    ltlag_mid <- log(0.0417); label("Lag time ALAG3 of the middle release (day)")  # Wang 2024 Table 2 Consta: ALAG3 = 0.0417 day (%RSE 119)
    ld3 <- log(23.9); label("Zero-order duration D3 of the middle release into the second depot (day)")  # Wang 2024 Table 2 Consta: D3 = 23.9 day (%RSE 1.93)
    lka_mid <- log(0.0830); label("First-order absorption rate constant K32 from the middle-release depot into the central compartment (1/day)")  # Wang 2024 Table 2 Consta: K32 = 0.0830 /day (%RSE 7.39)

    # --- Release fractions -------------------------------------------
    # F1 and F2 are estimated and F3 = 1 - F1 - F2 is derived (Wang 2024
    # Table 2 footnote). Both are carried on a logit scale, with F2
    # expressed relative to the non-main-release remainder, exactly as
    # in the companion Rykindo model. At eta = 0 this reproduces the
    # published point estimates: F1 = 0.755, F2 = (1 - 0.755) * 0.4857
    # = 0.119 and F3 = 0.126 (Table 2 reports F3 = 0.126).
    logitffo <- qlogis(0.755); label("Logit of the main-release dose fraction F1 (logit units)")  # Wang 2024 Table 2 Consta: F1 = 0.755 (%RSE 2.12)
    logitfburst <- qlogis(0.119 / (1 - 0.755)); label("Logit of the immediate-release dose fraction F2 expressed relative to the non-main-release remainder, so F2 = (1 - F1) * expit(logitfburst)")  # Wang 2024 Table 2 Consta: F2 = 0.119 (%RSE 8.74); 0.119 / (1 - 0.755) = 0.4857

    # --- Between-subject variability ---------------------------------
    # Wang 2024 Table 2 reports IIV for KA and CL in the Consta column
    # as CV% (26% and 34%); converted here with omega^2 = log(1 + CV^2).
    #
    # The Consta column also carries IIV entries for F2 (0.67), ALAG1
    # (0.87) and F1 (0.57) WITHOUT a percent sign, unlike every other
    # IIV entry in the table. Those three are set to fixed(0) here
    # because no stated scale for them survives a check against the
    # paper's own visual predictive check (Wang 2024 Fig. 4a): read as
    # CV%, as a log-scale variance, or as a log-scale standard
    # deviation, the implied variability on the 27-day main-release lag
    # smears the single-dose concentration peak so badly that the
    # simulated median peak falls to about a quarter of the published
    # one and the 97.5th-to-median ratio more than doubles. With these
    # three variances set to zero the simulated prediction interval
    # reproduces Fig. 4a closely. See the vignette 'Assumptions and
    # deviations' section for the full comparison and for how a user
    # can restore them. D2, ALAG3, K32 and D3 are reported as
    # 'Fix to 0' in Table 2.
    etalcl ~ log(1 + 0.34^2)  # Wang 2024 Table 2 Consta: IIV on CL = 34% CV (shrinkage 0.82%)
    etalka ~ log(1 + 0.26^2)  # Wang 2024 Table 2 Consta: IIV on KA = 26% CV (shrinkage 14.4%)
    etald2 ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on D2 'Fix to 0'
    etalogitffo ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on F1 reported as 0.57 with no scale; see the note above
    etalogitfburst ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on F2 reported as 0.67 with no scale; see the note above
    etaltlag_main ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on ALAG1 reported as 0.87 with no scale; see the note above
    etaltlag_mid ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on ALAG3 'Fix to 0'
    etald3 ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on D3 'Fix to 0'
    etalka_mid ~ fixed(0)  # Wang 2024 Table 2 Consta: IIV on K32 'Fix to 0'

    # --- Residual error ----------------------------------------------
    addSd <- 0.0546; label("Additive residual error on the active-moiety concentration (ng/mL)")  # Wang 2024 Table 2 Consta: add = 0.0546 (%RSE 74.1)
    propSd <- 0.32; label("Proportional residual error on the active-moiety concentration (fraction)")  # Wang 2024 Table 2 Consta: prop = 32% (%RSE 3.25)
  })

  model({
    # ==================================================================
    # 1. Individual parameters
    # ==================================================================
    cl <- exp(lcl + etalcl)
    ka <- exp(lka + etalka + e_study102_ka * STUDY_CT_USA_102)
    ka_mid <- exp(lka_mid + etalka_mid)
    d2 <- exp(ld2 + etald2)
    d3 <- exp(ld3 + etald3)
    tlag_main <- exp(ltlag_main + etaltlag_main)
    tlag_mid <- exp(ltlag_mid + etaltlag_mid)

    # Release fractions. Stick-breaking keeps F1 + F2 <= 1 so the
    # derived middle-release fraction F3 is never negative.
    ffo <- expit(logitffo + etalogitffo)
    fburst <- (1 - ffo) * expit(logitfburst + etalogitfburst)
    fmid <- 1 - ffo - fburst

    # ==================================================================
    # 2. Micro-constants. Wang 2024 Results: the initial fit gave an
    # elimination rate constant K = 0.379 /day essentially equal to the
    # absorption rate constant KA = 0.322 /day, indicating flip-flop
    # kinetics, so K was set to KA in the final model. The apparent
    # central volume is therefore derived rather than estimated.
    # ==================================================================
    kel <- ka
    vc <- cl / kel

    # ==================================================================
    # 3. ODE system (Wang 2024 Fig. 2b)
    #   depot   = main-release pool at the injection site (NONMEM cmt 1)
    #   central = plasma active moiety                    (NONMEM cmt 2)
    #   depot2  = middle-release pool                     (NONMEM cmt 3)
    # ==================================================================
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka_mid * depot2
    d/dt(central) <- ka * depot + ka_mid * depot2 - kel * central

    # ==================================================================
    # 4. Dose splitting, lag times and zero-order release durations.
    #
    # Each intramuscular injection is entered by the USER as THREE
    # simultaneous dose records, all carrying the full injected amount;
    # the f() multipliers below do the splitting:
    #   1. cmt = "depot",   rate =  0  (first-order main release)
    #   2. cmt = "central", rate = -2  (zero-order immediate release,
    #                                   modelled duration D2)
    #   3. cmt = "depot2",  rate = -2  (zero-order middle release,
    #                                   modelled duration D3)
    # rxode2 rejects two rate = -2 records that share the same time
    # stamp, so offset one of them by a negligible amount (for example
    # 1e-4 day, about 9 seconds) in the event table. See the vignette
    # for a worked event-table recipe.
    # ==================================================================
    f(depot) <- ffo
    alag(depot) <- tlag_main
    f(central) <- fburst
    dur(central) <- d2
    f(depot2) <- fmid
    alag(depot2) <- tlag_mid
    dur(depot2) <- d3

    # ==================================================================
    # 5. Observation. Dose is in mg and vc in L, so central / vc is
    # mg/L = ug/mL; multiply by 1000 for ng/mL.
    # ==================================================================
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
