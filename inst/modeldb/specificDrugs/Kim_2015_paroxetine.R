Kim_2015_paroxetine <- function() {
  description <- "One-compartment population PK model with first-order absorption for paroxetine (SSRI antidepressant) in Korean adults with major depressive disorder or anxiety disorder receiving therapeutic drug monitoring (Kim 2015)."
  reference <- "Kim J-R, Woo HI, Chun M-R, Lim S-W, Kim HD, Na HS, Chung MW, Myung W, Lee S-Y, Kim DK. Exposure-outcome analysis in depressed patients treated with paroxetine using population pharmacokinetics. Drug Des Devel Ther. 2015;9:5247-5255. doi:10.2147/DDDT.S84718"
  vignette <- "Kim_2015_paroxetine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    DOSE = list(
      description        = "Daily paroxetine dose administered at steady state (per-subject)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Use case (a): per-subject assigned daily dose used as a power-form covariate on CL/F to capture the nonlinear (saturable) elimination of paroxetine without instantiating an explicit Michaelis-Menten model. Reference 25 mg/day (median of the Kim 2015 cohort). The negative exponent (-0.363) reflects decreasing CL/F with increasing dose. Per Kim 2015 Methods, subjects were dosed at 10-52.5 mg/day (median 25 mg) with the bulk of observations between 10 and 25 mg/day.",
      source_name        = "DOSE"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on CL/F with reference 71 years (median of the Kim 2015 cohort; range 24-90 years). The negative exponent (-0.702) reflects decreasing CL/F with age, consistent with reduced hepatic CYP2D6 activity in elderly subjects.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 127L,
    n_studies      = 1L,
    age_range      = "24-90 years",
    age_median     = "71 years",
    weight_range   = "37.0-86.9 kg",
    weight_median  = "58 kg",
    sex_female_pct = 70.1,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Adults with major depressive disorder (65.4%) or anxiety disorder (34.6%; comprising generalized anxiety disorder, panic disorder, or social phobia) per DSM-IV criteria.",
    dose_range     = "10-52.5 mg PO once daily (one subject on twice-daily dosing); 23.2% on immediate-release tablets and 76.8% on controlled-release tablets.",
    regions        = "Republic of Korea (Samsung Medical Center, Seoul; single-center retrospective cohort 2005-2011).",
    n_observations = 271L,
    formulation_mix = "23.2% immediate-release; 76.8% controlled-release (relative bioavailability of CR 0.67 per the Paxil-CR product monograph).",
    notes          = "Baseline demographics per Kim 2015 Table 1. Retrospective TDM dataset (271 steady-state trough concentrations from 127 outpatients). No subjects on co-medication known to alter paroxetine PK. CYP2D6 genotype not characterised. Both serum albumin and weight data were missing for a minority of subjects; the published model carried these subjects with population-typical values."
  )

  ini({
    # Structural PK parameters - Kim 2015 Table 2.
    # The typical CL/F is estimated; Ka and Vd/F are fixed to literature values
    # (Venkatakrishnan & Obach 2005 Drug Metab Dispos 33(6):845-852 for Ka;
    # 17 L/kg x 60 kg = 1020 L from Findling 1999 and product-label data for
    # Vd/F) because the sparse trough-only TDM dataset (median 2 observations
    # per subject) cannot independently estimate Ka and V.
    lka <- fixed(log(0.908))
    label("Absorption rate constant Ka (1/h)")                              # Kim 2015 Table 2, fixed to literature (Venkatakrishnan 2005)
    lvc <- fixed(log(1020))
    label("Apparent volume of distribution Vd/F (L)")                       # Kim 2015 Table 2, fixed to literature (17 L/kg x 60 kg)
    lcl <- log(13.1)
    label("Apparent clearance CL/F at DOSE=25 mg/day, AGE=71 y (L/h)")      # Kim 2015 Table 2, RSE 5.9%

    # Covariate effects on CL/F (power form, centred at median values).
    # CL/F = exp(lcl + etalcl) * (DOSE/25)^e_dose_cl * (AGE/71)^e_age_cl
    e_dose_cl <- -0.363
    label("Power exponent of daily dose on CL/F (unitless, reference 25 mg/day)")  # Kim 2015 Table 2, RSE 34.4% (negative -> saturable elimination)
    e_age_cl  <- -0.702
    label("Power exponent of age on CL/F (unitless, reference 71 years)")          # Kim 2015 Table 2, RSE 30.1% (negative -> lower CL in elderly)

    # IIV on CL/F only - the published model reports CV% on the linear scale,
    # converted to log-normal variance via omega^2 = log(CV^2 + 1):
    #   CL/F BSV 40.2% -> omega^2 = log(0.402^2 + 1) = 0.1497839
    # No IIV on Ka or Vd (both fixed; the paper notes inclusion of further
    # IIVs did not improve fit).
    etalcl ~ 0.1497839                                                       # Kim 2015 Table 2, CV 40.2%, RSE 13.4%

    # Residual error - Kim 2015 Methods states an "additive error model with
    # log-transformed data" with sigma = 0.642 in log(ug/L). The NONMEM
    # log-additive (LTBS) pattern maps to nlmixr2's lnorm() residual:
    #   log(Cobs) = log(Cpred) + eps, eps ~ N(0, expSd^2).
    # For small expSd the linear-scale CV approximates expSd (here ~64%).
    expSd <- 0.642
    label("Log-normal residual SD on log(concentration) (unitless)")        # Kim 2015 Table 2, RSE 7.6%
  })

  model({
    # Individual PK parameters with power-form covariate effects on CL/F.
    # Both effects are centred at the cohort medians (DOSE = 25 mg/day,
    # AGE = 71 years), matching the typical-value estimate of lcl.
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl) * (DOSE / 25)^e_dose_cl * (AGE / 71)^e_age_cl

    # Micro-constant for the one-compartment system.
    kel <- cl / vc

    # ODE system - oral first-order absorption into the central compartment,
    # first-order elimination from central.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, vc in L -> mg/L; multiply by 1000 to convert
    # to ug/L (= ng/mL) to match the paper's reported observation unit.
    Cc <- 1000 * central / vc

    # Log-normal residual error (additive on log scale, per Kim 2015 Methods).
    Cc ~ lnorm(expSd)
  })
}
