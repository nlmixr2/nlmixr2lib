Chang_2023_rifapentine_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment PK model with first-order oral absorption for",
    "rifapentine in female BALB/c mice of the paucibacillary tuberculosis-preventive-therapy",
    "model, dosed orally twice daily on a tapered schedule designed to mimic the exposure",
    "profile of a long-acting injectable (Chang 2023, second PK/PD study). Apparent clearance",
    "increases geometrically by a fixed fraction per dosing day (autoinduction,",
    "CL/F on dosing day n = CL/F on day 1 x (1 + Ind_Factor)^(n - 1)). Parameters were",
    "fitted per regimen to naive-pooled mouse plasma concentrations by weighted least squares",
    "(no between-subject variability); the default values are the 'LAI target 0.6 ug/mL,",
    "2 doses' regimen fit, and the other regimen fits and the pilot / first-study fits",
    "(no induction term) are tabulated in the ini() comments and the vignette."
  )
  reference <- paste(
    "Chang YS, Li SY, Pertinez H, Betoudji F, Lee J, Rannard SP, Owen A, Nuermberger EL,",
    "Ammerman NC. Using Dynamic Oral Dosing of Rifapentine and Rifabutin to Simulate",
    "Exposure Profiles of Long-Acting Formulations in a Mouse Model of Tuberculosis",
    "Preventive Therapy. Antimicrob Agents Chemother. 2023;67(7):e00481-23.",
    "doi:10.1128/aac.00481-23."
  )
  vignette <- "Chang_2023_rifapentine_rifabutin_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "rifapentine",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "rifapentine",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (BALB/c)",
    sex = "female (male mice in the pilot PK study only)",
    n_subjects = paste(
      "Pilot PK: 84 female and 48 male uninfected mice (12 per dose group; 4 rifapentine",
      "groups in females, 2 in males). PK/PD study 1: 125 infected female mice in total",
      "(10-11 per rifapentine regimen). PK/PD study 2: 147 infected female mice in total",
      "(15 per single-LAI rifapentine regimen and 5 per two-LAI regimen)."
    ),
    n_studies = 3L,
    age_range = "8-10 weeks (pilot PK); 6-8 weeks at immunization (PK/PD studies)",
    weight_range = "doses in mg/kg computed for an average body mass of 20 g",
    disease_state = paste(
      "Pilot PK: uninfected. PK/PD studies: paucibacillary latent-tuberculosis model --",
      "aerosol immunization with M. bovis rBCG30, low-dose aerosol M. tuberculosis H37Rv",
      "challenge 6 weeks later, treatment started 6-7 weeks after challenge."
    ),
    dose_range = paste(
      "Pilot: 0.1, 0.3, 1 and 3 mg/kg once daily, 5 days per week (females; 0.3 and",
      "3 mg/kg in males). Study 1: 0.075-1.25 mg/kg twice daily with a 0.2-3.125 mg/kg",
      "loading dose, 7 days per week for 3 weeks. Study 2: 0.1-2.4 mg/kg twice daily,",
      "dose changed every 4 days over 8 weeks (Table S5); total 35.4-129.2 mg/kg."
    ),
    regions = "United States (Johns Hopkins University); analysis at the University of Liverpool",
    notes = paste(
      "Rifapentine tablets (Priftin) suspended in 0.05% agarose, 0.2 mL by oral gavage.",
      "Destructive / staggered mandibular-bleed sampling; LC-MS/MS LLOQ 0.10 ug/mL. All",
      "fits naive-pooled the individual mouse concentrations and minimised least squares",
      "weighted by 1/(predicted)^2 with pracma::lsqnonlin in R 4.0.3."
    )
  )

  ini({
    # Default values: Table 5, regimen 'RPT LAI target 0.6 ug/mL, 2 doses' (second PK/PD
    # study), the regimen fit with the lowest RSE on every parameter. The other Table 5
    # regimen fits (CL/F, V/F, Ka, Ind_Factor):
    #   LAI 0.6 ug/mL, 1 dose  : 0.02 [9.6], 0.37 [122.0], 0.23 [128.9], 1.1 [17.8]
    #   LAI 2 ug/mL, 1 dose    : 0.02 [9.4], 0.64 [29.9],  5.09 [0.004], 1.1 [23.4]
    #   LAI 2 ug/mL, 2 doses   : 0.02 [7.8], 0.48 [25.5],  0.62 [46.3],  1.2 [19.5]
    #   LAI 3.5 ug/mL, 1 dose  : 0.02 [12.6], 0.45 [70.4], 0.57 [101.4], 1.1 [24.7]
    # Pilot PK (Table 1, female, no induction term; CL/F, V/F, Ka):
    #   0.1 mg/kg 0.024/0.743/2.57; 0.3 mg/kg 0.029/0.736/1.75;
    #   1 mg/kg 0.021/0.643/1.64; 3 mg/kg 0.016/0.418/1.33
    # Pilot PK (Table S2, male): 0.3 mg/kg 0.020/0.449/1.05; 3 mg/kg 0.015/0.287/0.75
    # PK/PD study 1 (Fig. S3 panels A-D, no induction term):
    #   0.075 mg/kg BID 0.03/1.47/0.89; 0.25 mg/kg BID 0.03/1.10/1.01;
    #   0.75 mg/kg BID 0.02/0.88/1.33; 1.25 mg/kg BID 0.02/0.68/0.81
    lka <- log(0.42); label("First-order absorption rate constant (Ka, 1/h)") # Table 5, RPT LAI target 0.6 ug/mL 2 doses: Ka = 0.42 1/h [23.3% RSE]
    lcl <- log(0.02); label("Apparent clearance on dosing day 1 (CL/F, L/h/kg)") # Table 5, RPT LAI target 0.6 ug/mL 2 doses: CL/F = 0.02 L/h/kg [4.7% RSE]
    lvc <- log(0.35); label("Apparent volume of distribution (V/F, L/kg)") # Table 5, RPT LAI target 0.6 ug/mL 2 doses: V/F = 0.35 L/kg [16.8% RSE]
    cl_time_fday <- 0.011; label("Fractional cumulative increase in CL/F per dosing day (Ind_Factor, fraction)") # Table 5, RPT LAI target 0.6 ug/mL 2 doses: Ind_Factor = 1.1 percent [11.9% RSE]

    # Residual error: fits minimised least squares weighted by 1/(predicted)^2 (Methods,
    # PK and PK/PD analyses), i.e. a proportional error model, but no residual magnitude is
    # reported; fixed at zero so simulations reproduce the published fitted curves.
    propSd <- fixed(0); label("Proportional residual error (fraction; magnitude not reported)") # Methods: 1/pred^2-weighted least squares, residual SD not reported
  })

  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    # Methods equation: CL/F(day n) = CL/F(day 1) x (1 + Ind_Factor)^(n - 1). Dosing was
    # 7 days per week from t = 0 at the first dose, so dosing day n - 1 = floor(t / 24).
    cl <- exp(lcl) * (1 + cl_time_fday)^floor(t / 24)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
