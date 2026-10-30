Cha_2024_zolpidem <- function() {
  description <- paste(
    "One-compartment population PK model for oral immediate-release",
    "zolpidem in 23 healthy Korean volunteers (12 male, 11 female) given a",
    "single 10 mg film-coated tablet (Cha 2024). First-order absorption",
    "after an absorption lag time, first-order elimination; apparent (oral)",
    "clearance and volume with correlated inter-individual variability, and",
    "proportional residual error. No covariate was retained (sex was",
    "tested and was not significant). The concentrations were digitised",
    "from the individual-profile figures of a regulatory study report, and",
    "the model was built to simulate plasma concentrations after zolpidem",
    "overdose."
  )
  reference <- paste(
    "Cha HJ, Han S, Pak KC, Kim H. Development of a Web Application for",
    "Simulating Plasma Drug Concentrations in Patients with Zolpidem",
    "Intoxication. Pharmaceutics. 2024;16(5):689.",
    "doi:10.3390/pharmaceutics16050689"
  )
  vignette <- "Cha_2024_zolpidem"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "zolpidem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "zolpidem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Cha 2024 Methods 2.2: 'Covariate analysis was performed,",
        "particularly in terms of sex'. Results 3.2: 'sex was not",
        "significant in the current study'. No point estimate is reported.",
        "Sex was the only covariate available because individual",
        "demographics could not be recovered from the digitised source."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23,
    n_studies = 1,
    age_range = "not reported for the analysed subset",
    age_mean = "30.2 (SD 5.9) years in males and 29.9 (SD 6.2) years in females of the original 30-subject study",
    weight_mean = "73.1 (SD 7.6) kg in males and 54.7 (SD 3.9) kg in females of the original 30-subject study",
    sex_female_pct = 47.8,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral 10 mg immediate-release film-coated zolpidem tablet",
    regions = "Korea",
    notes = paste(
      "Open-label, parallel-group, single-dose study (KCT0003934) of 30",
      "healthy Korean volunteers, 15 per sex, sampled at 0, 0.25, 0.5,",
      "0.75, 1, 1.5, 2, 3, 4, 6, 8 and 12 h (Cha 2024 Methods 2.1).",
      "Concentrations were digitised with WebPlotDigitizer from the",
      "individual-profile figures of the Korean Ministry of Food and Drug",
      "Safety report; 7 subjects were excluded because their figures",
      "duplicated those of other subjects, leaving 23 (12 male, 11 female;",
      "Results 3.1). Demographics are those of the original 30 subjects;",
      "the 23-subject subset could not be summarised. The digitised NCA",
      "agreed with the original NCA (Table 1)."
    )
  )

  ini({
    # Typical values: Cha 2024 Table 2, final-model 'Estimates (RSE %)'
    # column. CL/F and Vd/F are apparent parameters (bioavailability not
    # estimated). NONMEM 7.5 ADVAN14, FOCE-I (Methods 2.2).

    lka <- log(5.41)
    label("Absorption rate constant (1/h)") # Cha 2024 Table 2: Ka = 5.41 1/h (RSE 45.47%; bootstrap median 5.66, 95% CI 3.24-11.24).

    lcl <- log(16.9)
    label("Apparent clearance CL/F (L/h)") # Cha 2024 Table 2: CL/F = 16.90 L/h (RSE 8.40%; bootstrap median 16.87, 95% CI 14.29-19.50).

    lvc <- log(61.7)
    label("Apparent volume of distribution Vd/F (L)") # Cha 2024 Table 2: Vd/F = 61.70 L (RSE 5.82%; bootstrap median 61.33, 95% CI 54.89-68.64).

    ltlag <- log(0.394)
    label("Absorption lag time (h)") # Cha 2024 Table 2: ALAG = 0.394 h (RSE 6.24%; bootstrap median 0.394, 95% CI 0.385-0.493). Results 3.2: adding ALAG lowered the OFV by 25.3.

    # Inter-individual variability. Table 2 prints each IIV row as
    # 'omega2 <par> (CV %)' with values 158.91, 22.10 and 32.60, so the
    # printed numbers are CV% values and the omega^2 variances are
    # back-transformed here as log(1 + CV^2). That reading, rather than
    # omega x 100, reproduces the paper's own Monte Carlo summary in
    # Table 3 (for example 2.6% of 560 mg subjects never exceeding
    # 4000 ng/mL: 2.1% under log(1 + CV^2) versus 6.6% under omega x 100,
    # 4000 simulated subjects);
    # see the vignette for the full comparison. The Vd-CL correlation is
    # printed directly as rho = 0.853, so the covariance is
    # 0.853 * sqrt(0.101 * 0.04769) = 0.0592. No IIV on the lag time.

    etalka ~ 1.26 # Cha 2024 Table 2: IIV Ka = 158.91 CV% (RSE 43.65%; bootstrap median 157.81, 95% CI 69.76-444.98); log(1 + 1.5891^2) = 1.26.
    etalcl + etalvc ~ c(0.101, 0.0592, 0.04769) # Cha 2024 Table 2: IIV CL = 32.60 CV% (log(1 + 0.326^2) = 0.101), IIV Vd = 22.10 CV% (log(1 + 0.221^2) = 0.04769), correlation Vd-CL = 0.853 (bootstrap 95% CI 0.796-0.869).

    # Residual error. Methods 2.2: DV = IPRED * (1 + eps_prop). Table 2
    # prints 'Proportional error 0.284'; read as the SD of eps_prop (a 28.4%
    # proportional error) because the variance reading (SD 0.533) widens
    # the simulated between-subject SD of NCA lambda z from about 0.08 to
    # 0.12 1/h, further from the digitised NCA (Table 1 lambda z SD
    # 0.04-0.06 1/h); see the vignette.
    propSd <- 0.284
    label("Proportional residual error (fraction)") # Cha 2024 Table 2: proportional error = 0.284 (RSE 6.16%; bootstrap median 0.283, 95% CI 0.242-0.319).
  })

  model({
    # No covariate was retained (Results 3.2), so each individual
    # parameter is its typical value times exp(eta) (Methods 2.2:
    # Pi = theta_pop * exp(eta_i)).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    lag(depot) <- tlag

    # Dose in mg over a volume in L gives mg/L; x 1000 reports ng/mL, the
    # unit of the paper's figures and thresholds.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
