Prasanchaimontri_2024_favipiravir <- function() {
  description <- paste(
    "One-compartment first-order-absorption pharmacokinetic model for oral",
    "favipiravir in Thai adults with mild COVID-19 without pneumonia",
    "(Prasanchaimontri 2024). Apparent clearance changes linearly with each",
    "elapsed dosing day, CL/F(day n) = CL/F(day 1) * (1 + Factor * (n - 1)),",
    "describing the apparent change in favipiravir clearance with repeated",
    "dosing. The paper used a two-stage approach: each of the eight patients",
    "was fitted individually and the summary statistics of the individual",
    "estimates define the population. Typical values are the geometric means",
    "and the between-subject variances are the variances of the",
    "log-transformed individual estimates of CL/F, V/F and ka, the values the",
    "paper used for its log-normal population simulations. The daily change",
    "is carried at its arithmetic mean (2.4 percent per day) without",
    "between-subject variability. The paper fitted with 1/predicted^2",
    "weighting but reported no residual-error magnitude, so the proportional",
    "residual SD is fixed at zero."
  )
  reference <- paste(
    "Prasanchaimontri IO, Manosuthi W, Pertinez H, Owen A, Niyomnaitham S,",
    "Sirijatuphat R, Charoenpong L, Cressey TR, Copeland K, Mokmued P,",
    "Chokephaibulkit K (2024). Favipiravir pharmacokinetics in Thai adults with",
    "mild COVID-19: A sub-study of interpatient variability and ethnic",
    "differences in exposure. Pharmacol Res Perspect 12(6):e1233.",
    "doi:10.1002/prp2.1233"
  )
  vignette <- "Prasanchaimontri_2024_favipiravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "favipiravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "favipiravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "19-53 years",
    age_median = "39 years",
    weight_mean = "67.9 kg (SD 15.45 kg)",
    bmi_median = "27.9 kg/m^2 (range 18.0-33.6)",
    bsa_mean = "1.70 m^2 (SD 0.19 m^2)",
    sex_female_pct = 87.5,
    race_ethnicity = "Thai",
    disease_state = "Symptomatic, PCR-confirmed mild-to-moderate COVID-19 without pneumonia",
    dose_range = "Oral favipiravir 1800 mg twice daily on day 1, then 800 mg twice daily for 5 to 14 days (200 mg tablets)",
    regions = "Thailand",
    notes = paste(
      "PK sub-study of eight of the 62 patients in the favipiravir arm of a",
      "prospective randomized trial (Prasanchaimontri 2024 Methods and",
      "Results). Seven of eight were female; three had underlying diseases.",
      "Four patients were treated for 14 days, three for 7 days and one",
      "(Patient 2) stopped on day 4 because of a rash. Plasma was sampled",
      "pre-dose and 1 h post-dose on days 1-5, 7 and 14 and assayed by",
      "HPLC-UV (calibration range 0.1-50 mg/L). Patients 4 and 5 had",
      "measurable pre-dose concentrations, indicating favipiravir before",
      "study entry. Mean (SD) weight and BSA are from the Supplementary",
      "Methods; age and BMI from the Results."
    )
  )

  ini({
    # Two-stage analysis: the paper reports one individual fit per patient
    # (Table 1) and the summary statistics across the eight patients. No
    # population parameter was estimated, so every value is wrapped in
    # fixed(). The Methods state the population simulations drew the
    # parameters from log-normal distributions using these summary statistics.
    lka <- fixed(log(3.4)); label("First-order absorption rate constant ka (1/h)") # Table 1 'Geometric mean' K_A = 3.4 1/h
    lcl <- fixed(log(1.0)); label("Apparent clearance CL/F on dosing day 1 (L/h)") # Table 1 'Geometric mean' CL/F = 1.0 L/h
    lvc <- fixed(log(19.5)); label("Apparent volume of distribution V/F (L)") # Table 1 'Geometric mean' V/F = 19.5 L

    cl_time_slope <- fixed(0.024); label("Fractional change in CL/F per elapsed dosing day (1/day)") # Table 1 'Arithmetic mean' daily percent change in CL/F = 2.4 -> 0.024; equation e001

    etalka ~ fixed(1.8) # Table 1 'Variance of log-transformed parameter' K_A = 1.8
    etalcl ~ fixed(0.3) # Table 1 'Variance of log-transformed parameter' CL/F = 0.3
    etalvc ~ fixed(0.1) # Table 1 'Variance of log-transformed parameter' V/F = 0.1

    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)") # Methods: lsqnonlin with 1/predicted^2 weighting; no magnitude reported
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)

    # Methods equation e001: (CL/F)_day n = (CL/F)_day 1 * (1 + Factor * (n - 1)).
    # n is the dosing day counted from the first dose at time 0, so n - 1 is
    # the number of whole days elapsed (model time in hours).
    day_elapsed <- floor(t / 24)
    cl <- exp(lcl + etalcl) * (1 + cl_time_slope * day_elapsed)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
