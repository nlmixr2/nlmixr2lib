Chang_2023_rifabutin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment PK model with first-order oral absorption for",
    "rifabutin in female BALB/c mice of the paucibacillary tuberculosis-preventive-therapy",
    "model, dosed orally twice daily on a tapered schedule designed to mimic the exposure",
    "profile of a long-acting injectable (Chang 2023, second PK/PD study). Apparent clearance",
    "increases geometrically by a fixed fraction per dosing day (autoinduction,",
    "CL/F on dosing day n = CL/F on day 1 x (1 + Ind_Factor)^(n - 1)). Parameters were",
    "fitted to naive-pooled mouse plasma concentrations by weighted least squares (no",
    "between-subject variability); the default values are the 'LAI target 0.15 ug/mL,",
    "1 dose' regimen fit, the only rifabutin regimen of that study that could be fitted, and",
    "the pilot / first-study fits (no induction term) are tabulated in the ini() comments",
    "and the vignette."
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
      analyte = "rifabutin",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "rifabutin",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (BALB/c)",
    sex = "female (male mice in the pilot PK study only)",
    n_subjects = paste(
      "Pilot PK: 84 female and 48 male uninfected mice (12 per dose group; 3 rifabutin",
      "groups in females, 2 in males). PK/PD study 1: 125 infected female mice in total",
      "(10-11 per rifabutin regimen). PK/PD study 2: 147 infected female mice in total",
      "(15 per rifabutin regimen)."
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
      "Pilot: 0.3, 1 and 3 mg/kg once daily, 5 days per week (females; 0.3 and 3 mg/kg in",
      "males). Study 1: 0.75-5.5 mg/kg twice daily, 7 days per week for 3 weeks. Study 2:",
      "0.2-5 mg/kg twice daily, dose changed every 4 days over 8 weeks (Table S5); total",
      "70.4 and 254.8 mg/kg."
    ),
    regions = "United States (Johns Hopkins University); analysis at the University of Liverpool",
    notes = paste(
      "Rifabutin powder suspended in 0.05% agarose, 0.2 mL by oral gavage. Destructive /",
      "staggered mandibular-bleed sampling; LC-MS/MS LLOQ 0.05 ug/mL. All fits naive-pooled",
      "the individual mouse concentrations and minimised least squares weighted by",
      "1/(predicted)^2 with pracma::lsqnonlin in R 4.0.3. The study-2 'LAI target",
      "0.045 ug/mL' regimen was not fitted (most 9-h and 24-h samples below the LLOQ)."
    )
  )

  ini({
    # Default values: Table 5, regimen 'RFB LAI target 0.15 ug/mL, 1 dose' (second PK/PD
    # study); the 'LAI target 0.045 ug/mL' regimen was 'Not fitted' (Table 5 footnote f).
    # Pilot PK (Table 1, female, no induction term; CL/F, V/F, Ka):
    #   0.3 mg/kg 0.674/7.360/49.10; 1 mg/kg 1.040/7.490/16.40; 3 mg/kg 1.170/8.670/8.48
    # Pilot PK (Table S2, male): 0.3 mg/kg 1.39/9.87/48.9; 3 mg/kg 1.79/11.0/47.7
    # PK/PD study 1 (Fig. S5 panels A-D, no induction term):
    #   1.5 mg/kg BID 0.88/8.74/23.40; 5.5 mg/kg BID 0.86/8.56/24.09;
    #   0.75 mg/kg BID 0.77/8.12/25.33; 2.5 mg/kg BID 0.93/10.97/25.23
    # Several rifabutin Ka estimates carry an RSE of about 0 percent; the paper states no
    # parameter was fixed, so they are estimates in an uninformative region (absorption
    # complete before the first 0.5-1 h sample), not fixed values.
    lka <- log(10.03); label("First-order absorption rate constant (Ka, 1/h)") # Table 5, RFB LAI target 0.15 ug/mL 1 dose: Ka = 10.03 1/h [0.5% RSE]
    lcl <- log(1.06); label("Apparent clearance on dosing day 1 (CL/F, L/h/kg)") # Table 5, RFB LAI target 0.15 ug/mL 1 dose: CL/F = 1.06 L/h/kg [15.6% RSE]
    lvc <- log(13.23); label("Apparent volume of distribution (V/F, L/kg)") # Table 5, RFB LAI target 0.15 ug/mL 1 dose: V/F = 13.23 L/kg [19.8% RSE]
    cl_time_fday <- 0.009; label("Fractional cumulative increase in CL/F per dosing day (Ind_Factor, fraction)") # Table 5, RFB LAI target 0.15 ug/mL 1 dose: Ind_Factor = 0.9 percent [46.4% RSE]

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
