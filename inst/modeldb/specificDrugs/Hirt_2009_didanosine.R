Hirt_2009_didanosine <- function() {
  description <- "One-compartment population PK model for didanosine (ddI) administered once daily as buffered chewable Videx tablets in West African HIV-1-infected children; first-order absorption with ka fixed at 4 1/h, additive residual error, exponential IIV on CL/F and Vc/F with off-diagonal covariance"
  reference <- "Hirt D, Bardin C, Diagbouga S, Nacro B, Hien H, Zoure E, Rouet F, Ouiminga A, Urien S, Foulongne V, Van De Perre P, Treluyer JM, Msellati P. Didanosine population pharmacokinetics in West African human immunodeficiency virus-infected children administered once-daily tablets in relation to efficacy after one year of treatment. Antimicrob Agents Chemother. 2009;53(10):4399-4406. doi:10.1128/AAC.01187-08"
  vignette <- "Hirt_2009_didanosine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # The Hirt 2009 covariate search tested body weight, postnatal age, body
  # surface area, serum creatinine, amylase, ALT, AST, total bilirubin, and a
  # binary duration-of-treatment indicator (2 weeks vs 2-5 months) on both
  # CL/F and Vc/F. None met the OFV-decrease threshold (chi-squared 6.63
  # units, P < 0.01) and none were retained in the final model. These are
  # documented for provenance but never referenced in model(), so they live
  # in covariatesDataExcluded rather than covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 11-37 kg (Table 1). Screened on CL/F and Vc/F via linear (per kg), allometric, and generalized-additive forms; not retained.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (chronological time since birth)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 2.5-14 years, median 6.5 (Table 1). Screened (as 'postnatal age') on CL/F and Vc/F; not retained. None of the children were younger than 2.5 years, so ontogeny of CL could not be resolved in this dataset (Discussion).",
      source_name        = "AGE"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 0.50-1.22 m^2, median 0.69 (Table 1). Screened on CL/F and Vc/F via linear (per m^2) and allometric forms; not retained. BSA also drove the dosing rule (240 mg/m^2 QD), so per-patient administered doses were 50-300 mg.",
      source_name        = "BSA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 6.6-126.4 umol/L, median 63 (Table 1). Screened on CL/F and Vc/F; not retained.",
      source_name        = "creatinine"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 7-92 IU/L, median 26 (Table 1). Screened on CL/F and Vc/F; not retained.",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 12-200 IU/L, median 44 (Table 1). Screened on CL/F and Vc/F; not retained.",
      source_name        = "AST"
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort range 0.9-42.8 umol/L, median 6.5 (Table 1). Screened on CL/F and Vc/F; not retained. The paper also screened serum amylase (range 6-409 U/L, median 79) but amylase is not a canonical covariate column in nlmixr2lib.",
      source_name        = "bilirubin"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49,
    n_studies      = 1,
    n_observations = 183,
    age_range      = "2.5-14 years",
    age_median     = "6.5 years",
    weight_range   = "11-37 kg",
    sex_female_pct = 38.8,
    race_ethnicity = "West African paediatric cohort (race not stratified in source).",
    disease_state  = "HIV-1-infected children (CDC clinical category A, B, C, or N); antiretroviral-naive at enrolment (except for prophylaxis of mother-to-child transmission). Mean Z-scores -1.91 (weight-for-age) and -1.99 (height-for-age) at enrolment.",
    dose_range     = "Videx (didanosine) chewable / dispersible tablets, 240 mg/m^2 once daily, rounded to combinations of 25, 50, 100, or 200 mg tablets. Median actual administered dose 213 mg/m^2 (range 164-313). Tablets dissolved in water; children fasted for 2 h before and 1 h after the dose. Aluminium hydroxide antacid (Maalox) added when fewer than two pills were prescribed.",
    regions        = "Burkina Faso (Bobo-Dioulasso).",
    notes          = "BURKINAM-ANRS 12103 open phase II trial (NCT00122538), once-daily didanosine + lamivudine + efavirenz combination. Pharmacokinetic sampling at day 15 (40 children) or between months 2 and 5 (9 children) of treatment, pre-dose and 1, 2, 3, 6, 12, and 24 h post-dose (10 children) or pre-dose and 1 and 3 h (39 children). Baseline median viral load 5.5 log10 copies/mL; 18 stage A, 25 stage B, 5 stage C, 1 stage N. Body surface area median 0.69 m^2 (range 0.50-1.22)."
  )

  ini({
    # Structural parameters (Table 2 final model). CL/F and Vc/F are reported
    # in absolute units (L/h and L) rather than per body surface area, with no
    # significant covariate scaling.
    lka  <- fixed(log(4));   label("Absorption rate constant (ka, 1/h; fixed)") # Table 2 (ka fixed to give tmax ~ 0.5 h)
    lcl  <- log(208);        label("Apparent clearance (CL/F, L/h)")            # Table 2
    lvc  <- log(278);        label("Apparent central volume of distribution (Vc/F, L)") # Table 2

    # IIV on log-transformed structural parameters. The paper reports omega
    # as a CV%; the conventional log-normal mapping is omega^2 = log(CV^2 + 1).
    # var(CL/F) = log(1 + 1.27^2) = 0.9605
    # var(Vc/F) = log(1 + 0.83^2) = 0.5241
    # cor(CL/F, Vc/F) = 0.557 (Table 2), so
    # cov = 0.557 * sqrt(0.9605 * 0.5241) = 0.3952
    etalcl + etalvc ~ c(0.9605,
                        0.3952, 0.5241) # Table 2 (omega CL/F = 127%, omega Vc/F = 83%, r = 0.557)

    # Residual error: additive on the linear concentration scale (Methods,
    # 'Modeling strategy'; Table 2 sigma). Negative predictions in the VPC
    # arise from this additive structure (Fig. 2 caption).
    addSd <- 0.077; label("Additive residual error (mg/L)") # Table 2
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, Vc/F in L -> concentration in mg/L
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
