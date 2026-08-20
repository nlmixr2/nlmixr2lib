Wang_2024_risperidone_rykindo <- function() {
  description <- "Population PK model for the risperidone active moiety (risperidone + 9-OH-risperidone) after intramuscular gluteal injection of RYKINDO (LY03004), a biweekly long-acting injectable risperidone formulation, in adults with schizophrenia or schizoaffective disorder (Wang 2024). One-compartment disposition fed by three parallel release pathways out of the injection site: (1) an immediate zero-order release of fraction F2 directly into the central compartment over duration D2, (2) a middle release of fraction F3 = 1 - F1 - F2 delivered as a zero-order input of duration D3 into a second depot beginning ALAG3 after injection and then absorbed first-order with rate K32, and (3) a main first-order release of fraction F1 from the primary depot beginning ALAG1 after injection with absorption rate constant KA. Because the active moiety displays flip-flop kinetics, the elimination rate constant K was set equal to KA, so the apparent central volume is derived as V = CL/KA. Apparent clearance carries a sex effect (females about 18 percent lower) and a between-study effect (0.675-fold in the two US relative-bioavailability trials relative to the single-ascending-dose trial); KA is higher in the multiple-dose trial. The companion RISPERDAL CONSTA model from the same paper is modellib('Wang_2024_risperidone_consta')."
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
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed. Wang 2024 Table 2 reports separate apparent-clearance estimates for males (186.7 L/day) and females (153.4 L/day) in trial CT-1S01; the model encodes the male value as the reference and the female/male ratio as a log-multiplicative coefficient. The paper's Results describe the effect as 'about 15% lower clearance' in females; the ratio implied by the Table 2 point estimates is 153.4/186.7 = 0.822, i.e. 17.8% lower. Table 2 is used because it carries the final model estimates. The authors judged the sex effect 'unlikely clinically relevant'. No sex effect was retained for Consta.",
      source_name        = "gender"
    ),
    STUDY_CT_USA = list(
      description        = "1 = subject enrolled in one of the two US relative-bioavailability trials of the Wang 2024 Rykindo programme (CT-USA-104, NCT02186769, single dose 25/50 mg; or CT-USA-102, NCT02091388, 25 mg every 2 weeks x 5); 0 = subject enrolled in the single-ascending-dose trial CT-1S01 (NCT02055287)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CT-1S01, the single-ascending-dose trial, which anchors the reference apparent clearance)",
      notes              = "Time-fixed; set from the trial identifier. Switches apparent clearance of the active moiety by a factor of 0.675 relative to CT-1S01 (Wang 2024 Table 2 row 'CL ratio for CT-USA-104 and 102 to CT-1S01'). Applies to the Rykindo model only; the Consta model has a single clearance across trials. Note that Wang 2024 Results paragraph 'Population Pharmacokinetic Modeling of Rykindo' describes the clearance reduction as applying to 'the multiple-dose study' only, which conflicts with the Table 2 row label that names both CT-USA-104 and CT-USA-102; the table label is used here because it is the final-model parameter definition and because the resulting AUC prediction matches the observed CT-USA-104 exposures far better (see the vignette source-trace section). Follows the auto-approved STUDY_<id> canonical family; a group-level indicator over two trials, directly analogous to STUDY_DORZA_EARLY.",
      source_name        = "study"
    ),
    STUDY_CT_USA_102 = list(
      description        = "1 = subject enrolled in trial CT-USA-102 (NCT02091388), the multiple-dose relative-bioavailability trial in which either Rykindo or Consta 25 mg was given every 2 weeks for five injections; 0 = subject enrolled in CT-1S01 or CT-USA-104",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the single-dose trials CT-1S01 and CT-USA-104)",
      notes              = "Time-fixed; set from the trial identifier. Switches the first-order absorption rate constant KA of the main release: 0.288 -> 0.380 /day for Rykindo and 0.179 -> 0.271 /day for Consta (Wang 2024 Table 2, rows 'KA_CT-1S01 and CT-USA-104' and 'KA_CT-USA-102'). The authors attribute the faster apparent absorption in the multiple-dose trial to its sparser PK sampling scheme rather than to a formulation difference (Wang 2024 Discussion). Because the model sets the elimination rate constant equal to KA (flip-flop kinetics), this indicator also shifts the apparent central volume V = CL/KA. Used by both Wang 2024 models. Follows the auto-approved STUDY_<id> canonical family.",
      source_name        = "study"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "administration site", verified = TRUE),
    depot2  = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 97L,
    n_studies      = 3L,
    n_observations = "2216 plasma concentration records of risperidone and 9-OH-risperidone (Wang 2024 Results, 'Population Pharmacokinetic Modeling of Rykindo').",
    age_range      = "medians 51.0 / 52.0 / 54.0 years and means 48.0 (SD 9.9) / 48.0 (9.2) / 45.6 (9.6) years in CT-1S01 / CT-USA-104 / CT-USA-102 (Wang 2024 Table 1)",
    weight_range   = "medians 84.4 / 89.7 / 88.9 kg and means 85.1 (SD 14.5) / 85.8 (14.8) / 89.4 (16.7) kg in CT-1S01 / CT-USA-104 / CT-USA-102 (Wang 2024 Table 1)",
    bmi_range      = "medians 27.9 / 29.6 / 29.3 kg/m^2 and means 28.0 (SD 4.1) / 27.9 (4.9) / 29.6 (5.1) kg/m^2 (Wang 2024 Table 1)",
    sex_female_pct = 24.5,
    race_ethnicity = "Not reported in Wang 2024. Race was screened as a covariate but no race effect was retained in the final model.",
    disease_state  = "Stable adults with schizophrenia or schizoaffective disorder",
    dose_range     = "12.5, 25, 37.5 and 50 mg single intramuscular gluteal injection (CT-1S01); 25 and 50 mg single injection (CT-USA-104); 25 mg every 2 weeks for five injections (CT-USA-102)",
    regions        = "United States (CT-USA-104, CT-USA-102) and the CT-1S01 single-ascending-dose trial",
    cyp2d6_status  = "Across the three Rykindo cohorts (Wang 2024 Table 1): 4 poor, 41 intermediate, 53 extensive and 1 ultra-rapid metabolizer, 3 unknown. CYP2D6 genotype was not retained as a covariate because the active moiety (risperidone + 9-OH-risperidone) is insensitive to CYP2D6 status.",
    notes          = "Pooled from three phase 1 trials: CT-1S01 (NCT02055287), CT-USA-104 (NCT02186769) and CT-USA-102 (NCT02091388). Wang 2024 Table 1 lists 32 + 16 + 54 = 102 enrolled subjects across the three Rykindo cohorts, while the Results text states that 97 subjects contributed the 2216 PK records; the difference is not explained in the paper. Age, weight, body mass index, sex and race were screened as covariates; only sex (on clearance) and trial (on clearance and on the main-release absorption rate constant) were retained."
  )

  ini({
    # ==================================================================
    # Wang 2024 Table 2, Rykindo column (final model estimates; standard
    # errors from 200 bootstrap replicates because the covariance step
    # failed). All values are for the active moiety.
    #
    # Structure (Wang 2024 Fig. 2b): a single dose is split three ways
    # out of the injection site. F1 is released first-order from `depot`
    # after a lag ALAG1; F2 enters `central` directly as a zero-order
    # input of duration D2; the remainder F3 = 1 - F1 - F2 enters
    # `depot2` as a zero-order input of duration D3 starting at ALAG3
    # and is then absorbed first-order with rate K32.
    # ==================================================================

    # --- Disposition -------------------------------------------------
    lcl <- log(186.7); label("Apparent clearance CL/F of the active moiety in males in trial CT-1S01 (L/day)")  # Wang 2024 Table 2 Rykindo: CL_male_CT-1S01 = 186.7 L/day (%RSE 6.30)
    e_sexf_cl <- log(153.4 / 186.7); label("Log-multiplicative effect of female sex on apparent clearance (unitless)")  # Derived from Wang 2024 Table 2 Rykindo: CL_female_CT-1S01 / CL_male_CT-1S01 = 153.4 / 186.7 = 0.822 (female %RSE 6.96)
    e_studyctusa_cl <- log(0.675); label("Log-multiplicative effect of the CT-USA-104 / CT-USA-102 trials on apparent clearance (unitless)")  # Wang 2024 Table 2 Rykindo: 'CL ratio for CT-USA-104 and 102 to CT-1S01' = 0.675 (%RSE 6.03)

    # --- Main (third) release ----------------------------------------
    lka <- log(0.288); label("First-order absorption rate constant KA of the main release in trials CT-1S01 and CT-USA-104 (1/day)")  # Wang 2024 Table 2 Rykindo: KA_CT-1S01 and CT-USA-104 = 0.288 /day (%RSE 6.61)
    e_study102_ka <- log(0.380 / 0.288); label("Log-multiplicative effect of trial CT-USA-102 on the main-release absorption rate constant (unitless)")  # Derived from Wang 2024 Table 2 Rykindo: KA_CT-USA-102 / KA_CT-1S01 and CT-USA-104 = 0.380 / 0.288 = 1.319 (%RSE 3.66 on the 0.380 estimate)
    ltlag_main <- log(13.3); label("Lag time ALAG1 of the main release (day)")  # Wang 2024 Table 2 Rykindo: ALAG1 = 13.3 day (%RSE 0.571)

    # --- Immediate release -------------------------------------------
    ld2 <- log(0.764); label("Zero-order duration D2 of the immediate release into the central compartment (day)")  # Wang 2024 Table 2 Rykindo: D2 = 0.764 day (%RSE 14.7)

    # --- Middle release ----------------------------------------------
    ltlag_mid <- log(3.47); label("Lag time ALAG3 of the middle release (day)")  # Wang 2024 Table 2 Rykindo: ALAG3 = 3.47 day (%RSE 4.63)
    ld3 <- log(2.18); label("Zero-order duration D3 of the middle release into the second depot (day)")  # Wang 2024 Table 2 Rykindo: D3 = 2.18 day (%RSE 19.34)
    lka_mid <- log(0.118); label("First-order absorption rate constant K32 from the middle-release depot into the central compartment (1/day)")  # Wang 2024 Table 2 Rykindo: K32 = 0.118 /day (%RSE 8.55)

    # --- Release fractions -------------------------------------------
    # F1 and F2 are estimated and F3 = 1 - F1 - F2 is derived (Wang 2024
    # Table 2 footnote). Both estimated fractions are carried on a logit
    # scale so each stays inside (0, 1), and F2 is expressed relative to
    # the non-main-release remainder (stick-breaking) so that
    # F1 + F2 <= 1 and therefore F3 >= 0 holds for every simulated
    # subject. At eta = 0 this reproduces the published point estimates
    # exactly: F1 = 0.430, F2 = (1 - 0.430) * 0.2596 = 0.148 and
    # F3 = 0.422 (Table 2 reports F3 = 0.421 from the rounded values).
    # See the vignette 'Assumptions and deviations' section.
    logitffo <- qlogis(0.430); label("Logit of the main-release dose fraction F1 (logit units)")  # Wang 2024 Table 2 Rykindo: F1 = 0.430 (%RSE 5.45)
    logitfburst <- qlogis(0.148 / (1 - 0.430)); label("Logit of the immediate-release dose fraction F2 expressed relative to the non-main-release remainder, so F2 = (1 - F1) * expit(logitfburst)")  # Wang 2024 Table 2 Rykindo: F2 = 0.148 (%RSE 7.59); 0.148 / (1 - 0.430) = 0.2596

    # --- Between-subject variability ---------------------------------
    # Wang 2024 Table 2 reports IIV for Rykindo as CV%; converted to
    # log-scale variance with omega^2 = log(1 + CV^2). The two release
    # fractions carry their variability on the logit scale instead; the
    # logit-scale variances below were chosen by moment matching so the
    # simulated coefficient of variation of F1 and of F2 reproduces the
    # published 52% and 57% (see the vignette 'Assumptions and
    # deviations' section).
    # ALAG1, ALAG3, K32 and D3 are reported as 'Fix to 0' in Table 2.
    etalcl ~ log(1 + 0.35^2)  # Wang 2024 Table 2 Rykindo: IIV on CL = 35% CV (shrinkage 5.38%)
    etalka ~ log(1 + 0.27^2)  # Wang 2024 Table 2 Rykindo: IIV on KA = 27% CV (shrinkage 7.84%)
    etald2 ~ log(1 + 0.52^2)  # Wang 2024 Table 2 Rykindo: IIV on D2 = 52% CV (shrinkage 19.8%)
    etalogitffo ~ 1.3958  # Wang 2024 Table 2 Rykindo: IIV on F1 = 52% CV (shrinkage 12.0%); logit-scale variance moment-matched to CV(F1) = 0.52
    etalogitfburst ~ 0.2605  # Wang 2024 Table 2 Rykindo: IIV on F2 = 57% CV (shrinkage 14.7%); logit-scale variance moment-matched to CV(F2) = 0.57
    etaltlag_main ~ fixed(0)  # Wang 2024 Table 2 Rykindo: IIV on ALAG1 'Fix to 0'
    etaltlag_mid ~ fixed(0)  # Wang 2024 Table 2 Rykindo: IIV on ALAG3 'Fix to 0'
    etald3 ~ fixed(0)  # Wang 2024 Table 2 Rykindo: IIV on D3 'Fix to 0'
    etalka_mid ~ fixed(0)  # Wang 2024 Table 2 Rykindo: IIV on K32 'Fix to 0'

    # --- Residual error ----------------------------------------------
    addSd <- 0.363; label("Additive residual error on the active-moiety concentration (ng/mL)")  # Wang 2024 Table 2 Rykindo: add = 0.363 (%RSE 6.50, shrinkage 8.15%)
    propSd <- 0.26; label("Proportional residual error on the active-moiety concentration (fraction)")  # Wang 2024 Table 2 Rykindo: prop = 26% (%RSE 3.67)
  })

  model({
    # ==================================================================
    # 1. Individual parameters
    # ==================================================================
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_studyctusa_cl * STUDY_CT_USA)
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
