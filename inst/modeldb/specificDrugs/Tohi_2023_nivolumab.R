Tohi_2023_nivolumab <- function() {
  description <- "Two-compartment population PK model for intravenous nivolumab (anti-PD-1 IgG4) in Japanese adults with non-small cell lung cancer, with power effects of serum albumin and eGFR on clearance and a sigmoid Emax decrease of clearance with time since the first dose (Tohi 2023)"
  reference <- "Tohi M, Irie K, Mizuno T, Okuyoshi H, Hirabatake M, Ikesue H, Muroi N, Eto M, Fukushima S, Tomii K, Hashida T. Population Pharmacokinetics of Nivolumab in Japanese Patients with Nonsmall Cell Lung Cancer. Ther Drug Monit. 2023;45(1):110-116. doi:10.1097/FTD.0000000000000996"
  vignette <- "Tohi_2023_nivolumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Nivolumab was measured in remnant serum by LC-MS/MS (Methods, 'Sample
  # Collection and Pharmacokinetic Measurement'); the central compartment is
  # the sampled serum pool and the peripheral compartment its distribution
  # partner in a standard two-compartment model.
  compartmentData <- list(
    central = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (measured approximately every 2 weeks; Methods 'Covariate Analysis on Nivolumab Clearance'). Power effect on CL with reference 3.6 g/dL (the Table 1 median and the value printed in the final CL equation). The source reports albumin in g/dL; the canonical column is g/L, so model() converts with alb_gdL <- ALB * 0.1 before applying the published coefficient.",
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying eGFR calculated with the Japanese-population equation (Methods 'Patients', reference 14 of the source). Power effect on CL with reference 70 mL/min/1.73 m^2 (the Table 1 median and the value printed in the final CL equation).",
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL in the forward-inclusion step (Table S1, dOFV -0.41, p = 0.522) and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL as a time-varying covariate (Table S1, dOFV -0.365, p = 0.546) and not retained. Body weight enters only the 3 mg/kg dose in the published Monte Carlo simulation."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL as a time-varying covariate (Table S1, dOFV -0.944, p = 0.331) and not retained."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on CL (Table S1, dOFV -0.011, p = 0.916) and not retained."
    ),
    ECOG = list(
      description = "ECOG performance status",
      units = "(score)",
      type = "categorical",
      reference_category = NULL,
      notes = "Screened on CL (Table S1, dOFV -0.118, p = 0.731) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_observations = 223L,
    age_range = "38-83 years",
    age_median = "69 years",
    weight_range = "36.8-80.5 kg",
    weight_median = "62.7 kg",
    sex_female_pct = 26.5,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Non-small cell lung cancer (22 adenocarcinoma, 11 squamous cell carcinoma, 1 not otherwise specified) treated with nivolumab monotherapy",
    dose_range = "240 mg flat dose IV every 2 weeks (33 patients) or every 3 weeks (1 patient); median 16.5 doses (range 1-75), median treatment period 296 days (range 14-1217)",
    regions = "Japan (single centre, Kobe City Medical Center General Hospital)",
    renal_function = "eGFR median 70 (range 29-144) mL/min/1.73 m^2",
    albumin = "Serum albumin median 3.6 (range 2.5-4.8) g/dL",
    performance_status = "ECOG PS 0: 10, PS 1: 22, PS 2: 1, unknown: 1",
    notes = "Real-world opportunistic sampling of remnant serum (1-12 samples per patient). Baseline demographics from Tohi 2023 Table 1; longitudinal albumin and eGFR summaries in Supplemental Table S3."
  )

  ini({
    # Structural parameters at the reference patient (ALB 3.6 g/dL, eGFR
    # 70 mL/min/1.73 m^2, t = 0). Units as reported: L and L/h.
    lcl <- log(0.0064); label("Clearance at the first dose for the reference patient (L/h)") # Table 2: CL = 0.0064 L/h (RSE 5%)
    lvc <- log(2.28); label("Central volume of distribution (L)") # Table 2: V1 = 2.28 L (RSE 17%)
    lvp <- log(1.81); label("Peripheral volume of distribution (L)") # Table 2: V2 = 1.81 L (RSE 11%)
    lq <- log(0.018); label("Intercompartmental clearance (L/h)") # Table 2: Q = 0.018 L/h (RSE 26%)

    # Covariate effects on CL, power form (cov / median)^theta (Methods
    # continuous-covariate equation and the Results final CL equation).
    e_alb_cl <- -1.48; label("Power exponent of serum albumin on CL (unitless)") # Table 2: theta for albumin = -1.48 (RSE 33%)
    e_crcl_cl <- 0.411; label("Power exponent of eGFR on CL (unitless)") # Table 2: theta for eGFR = 0.411 (RSE 23%)

    # Time-varying CL, sigmoid Emax in time since the first dose. All three
    # parameters were taken from Osawa 2019 (source reference 12) and held
    # constant; Methods and Table 2 both print T50 = 1510 h, while the Results
    # final CL equation prints 1501 h (see the vignette Errata).
    cl_time_max <- fixed(-0.285); label("Maximal change in CL on the log scale (unitless)") # Table 2: theta for EMAX = -0.285 FIX; Methods: Emax = -0.285 from Osawa et al
    lcl_t50 <- fixed(log(1510)); label("Log of the time after the first dose at which the change in CL is half of its maximum (log h)") # Table 2: theta for TM50 = 1510 FIX; Methods: T50 = 1510 h
    lcl_time_hill <- fixed(log(2.02)); label("Log of the Hill coefficient of time on CL (log unitless)") # Table 2: theta for HILL = 2.02 FIX; Methods: gamma = 2.020

    # IIV, exponential (Methods). Table 2 reports CV%; converted with
    # omega^2 = log(1 + CV^2). IIV on Q and V2 could not be estimated and was
    # held at zero (Results 'Base Model'), so those parameters carry no eta.
    etalcl ~ 0.03695 # Table 2: IIV for CL = 19.4 CV%; log(1 + 0.194^2) = 0.03695
    etalvc ~ 0.19425 # Table 2: IIV for V1 = 46.3 CV%; log(1 + 0.463^2) = 0.19425

    propSd <- 0.160; label("Proportional residual error (fraction)") # Table 2: proportional error = 16.0 CV%
  })
  model({
    # Albumin is supplied in canonical g/L; the published coefficient and
    # reference (3.6) are in g/dL.
    alb_gdL <- ALB * 0.1

    # Time-varying clearance; t is time since the first dose in hours, so the
    # first dose must be at t = 0.
    cl_t50 <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl_time <- exp(cl_time_max * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))
    cl <- exp(lcl + etalcl) *
      (alb_gdL / 3.6)^e_alb_cl *
      (CRCL / 70)^e_crcl_cl *
      cl_time
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
