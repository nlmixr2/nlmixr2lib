Stott_2022_amphotericinB_liposomal <- function() {
  description <- "Two-compartment IV-infusion population PK model for liposomal amphotericin B (AmBisome) in adults with HIV-associated cryptococcal meningoencephalitis given high-dose, short-course regimens (AMBITION-cm Phase II and III), with no covariates; disposition is written with explicit k12 / k21 micro-constants and no q / vp pair, so solve it with rxSolve(useLinCmt = FALSE) or the peripheral compartment is silently discarded (Stott 2022)"
  reference <- "Stott KE, Moyo M, Ahmadu A, Kajanga C, Gondwe E, Chimang'anga W, Chasweka M, Leeme TB, Molefi M, Chofle A, Bidwell G, Changalucha J, Unsworth J, Jimenez-Valverde A, Lawrence DS, Mwandumba HC, Lalloo DG, Harrison TS, Jarvis JN, Hope W, Martson AG. Population pharmacokinetics of liposomal amphotericin B in adults with HIV-associated cryptococcal meningoencephalitis. J Antimicrob Chemother. 2023;78(1):276-283. doi:10.1093/jac/dkac389"
  vignette <- "Stott_2022_amphotericinB_liposomal"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Stott 2022 Results, 'Population PK model': X(1) and X(2) are the amounts
  # of drug in mg in the central and peripheral compartments. Materials and
  # methods, 'Bioanalysis of PK samples': the UPLC-MS/MS assay measured TOTAL
  # (liposome-associated plus non-liposome-associated) amphotericin B in plasma.
  compartmentData <- list(
    central = list(
      analyte = "amphotericin B (total, liposome-associated plus free)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "amphotericin B (total, liposome-associated plus free)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # Stott 2022 Results, 'Population PK model': multivariate linear regression
  # of the Bayesian posterior clearance and volume found no significant
  # association with any screened covariate, so the base model was retained.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by bidirectional stepwise multivariate linear regression against the Bayesian posterior estimates of the base-model PK parameters; not retained (Stott 2022 Results, 'Population PK model'). Cohort median 52.0 kg (IQR 46.5-58.7; Table 1). Dosing was weight-based (mg/kg), so the dose in mg still depends on weight.",
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Stott 2022 Results). Cohort median 37 years (IQR 32-43; Table 1).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Named among the covariates with no significant association in Stott 2022 Results, although Methods lists only age, weight, CD4+ count and serum creatinine. 32 of 87 patients were female (Table 1).",
      source_name = "sex"
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Stott 2022 Results). Table 1 prints the cohort median as 64.0 'mmol/L' (IQR 58.0-87.2); the unit is a misprint for umol/L, which is the only physiologically possible reading and the unit the Toxicity section uses for its 207 umol/L threshold.",
      source_name = "baseline serum creatinine"
    ),
    CD4_ABS = list(
      description = "Baseline CD4+ T-cell count",
      units = "cells/mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Stott 2022 Results). Cohort median 30 cells/mm^3 (IQR 12-60; Table 1).",
      source_name = "CD4+ cell count"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 2L,
    age_range = "IQR 32-43 years",
    age_median = "37 years",
    weight_range = "IQR 46.5-58.7 kg",
    weight_median = "52.0 kg",
    sex_female_pct = 36.8,
    race_ethnicity = "Not reported; all patients were enrolled at sub-Saharan African sites (Botswana, Tanzania, Malawi)",
    disease_state = "HIV-associated cryptococcal meningoencephalitis. Baseline median (IQR): haemoglobin 11.0 (9.7-12.25) g/dL, WBC 5.0 (3.5-7.1) x10^9/L, platelets 255 (181-328.7) x10^9/L, serum creatinine 64.0 (58.0-87.2) umol/L, CD4+ count 30 (12-60) cells/mm^3 (Stott 2022 Table 1).",
    dose_range = "AmBisome by 2 h IV infusion after pre-hydration with 1 L 0.9% saline plus 20 mmol KCl. Phase II arms: 3 mg/kg/day for 14 days; 10 mg/kg on day 1 only; 10 mg/kg on day 1 plus 5 mg/kg on day 3; 10 mg/kg on day 1 plus 5 mg/kg on days 3 and 7 (all with fluconazole 1200 mg/day for 14 days). Phase III intervention arm: single 10 mg/kg dose on day 1 plus 14 days of flucytosine 100 mg/kg/day and fluconazole 1200 mg/day.",
    regions = "Botswana (Princess Marina Hospital, Gaborone) and Tanzania (Bugando Medical Centre and Sekou Toure Hospital, Mwanza) for Phase II; Malawi (Queen Elizabeth Central Hospital, Blantyre) for Phase III",
    notes = "AMBITION-cm PK substudies: 56 patients from all four Phase II arms (January 2015 to August 2016) and 31 from the Phase III single-dose intervention arm (November 2018 to October 2019). 565 plasma observations, mean 6.5 (range 2-12) per patient. Phase II sampling at the end of infusion, 6 h and 24 h; Phase III at 0, 2, 4, 7, 12 and 23 h after the start of the day-1 infusion and at 2, 4, 7, 12 and 23 h on day 7. Total amphotericin B assayed in plasma by UPLC-MS/MS, calibration range 0.25-50.0 mg/L, LLOQ 0.25 mg/L, CV < 9.0%. Fitting used the nonparametric adaptive grid (NPAG) algorithm of Pmetrics 1.9.7 (Stott 2022 Materials and methods and Results)."
  )

  ini({
    # Structural parameters: Stott 2022 Table 2. Pmetrics NPAG reports the
    # MEAN, MEDIAN and SD of each parameter's nonparametric marginal
    # distribution. The means are the typical values because the Abstract
    # reports them as the population PK parameter estimates ('Mean (SD)
    # population PK parameter estimates were: clearance 0.416 (0.363) L/h ...').
    lcl <- log(0.416)
    label("Clearance (L/h)") # Stott 2022 Table 2, 'Clearance (L/h)' mean = 0.416 (median 0.345)
    lvc <- log(4.566)
    label("Central volume of distribution (L)") # Stott 2022 Table 2, 'Volume (L)' mean = 4.566 (median 3.698)
    lk12 <- log(2.222)
    label("Central-to-peripheral first-order rate constant KCP (1/h)") # Stott 2022 Table 2, 'KCP (h-1)' mean = 2.222 (median 0.218)
    lk21 <- log(2.951)
    label("Peripheral-to-central first-order rate constant KPC (1/h)") # Stott 2022 Table 2, 'KPC (h-1)' mean = 2.951 (median 0.484)

    # Inter-individual variability. NPAG estimates a discrete, non-Gaussian
    # joint distribution and Table 2 reports only its marginal mean and SD, so
    # each marginal is approximated by a log-normal with the same coefficient
    # of variation, omega^2 = log(1 + (SD / mean)^2), with the reported mean
    # carried as the log-normal median. No correlations are published.
    etalcl ~ 0.566123 # Table 2, 'Clearance' SD 0.363 on mean 0.416 -> log(1 + (0.363/0.416)^2)
    etalvc ~ 0.682635 # Table 2, 'Volume' SD 4.518 on mean 4.566 -> log(1 + (4.518/4.566)^2)
    etalk12 ~ 1.18612 # Table 2, 'KCP' SD 3.351 on mean 2.222 -> log(1 + (3.351/2.222)^2)
    etalk21 ~ 1.06546 # Table 2, 'KPC' SD 4.070 on mean 2.951 -> log(1 + (4.070/2.951)^2)

    # Residual error. Pmetrics combines the assay error polynomial with the
    # selected 'additive' (lambda) term in quadrature, sqrt(assaySD^2 +
    # lambda^2), which is the nlmixr2 default combined2 form of add() + prop()
    # when the polynomial intercept is zero. Stott 2022 prints neither the
    # polynomial nor the fitted lambda, so lambda is fixed to zero and the
    # assay SD is represented by the reported assay coefficient of variation.
    addSd <- fixed(0)
    label("Additive (Pmetrics lambda) residual SD, not reported and set to zero (mg/L)") # Stott 2022 Results: 'additive error was selected for the final model'; the fitted value is not printed
    propSd <- fixed(0.09)
    label("Proportional residual SD, taken from the reported assay imprecision (fraction)") # Stott 2022 Materials and methods: 'The coefficient of variation was <9.0% over the concentration range 0.25-50.0 mg/L'
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    # Stott 2022 Equations 1 and 2. R(1), the 2 h zero-order IV infusion into
    # the central compartment, is encoded on the dosing rows.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Stott 2022 Equation 3: Y(1) = X(1) / V; mg / L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
