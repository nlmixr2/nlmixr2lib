Wu_2022_meropenem <- function() {
  description <- "One-compartment IV-infusion population PK model for meropenem in Chinese neonates and young infants (Wu 2022). V scales linearly with current body weight (reference 2.68 kg); CL scales allometrically with current body weight (fixed exponent 0.75) and with a maturation factor that is a power function of postnatal age (reference 12 days) times a power function of gestational age (reference 36.5 weeks). IIV on V and CL; proportional residual error."
  reference <- "Wu YE, Kou C, Li X, Tang BH, Yao BF, Hao GX, Zheng Y, van den Anker J, You DP, Shen AD, Zhao W. Developmental Population Pharmacokinetics-Pharmacodynamics of Meropenem in Chinese Neonates and Young Infants: Dosing Recommendations for Late-Onset Sepsis. Children (Basel). 2022;9(12):1998. doi:10.3390/children9121998"
  vignette <- "Wu_2022_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Wu 2022 reports current weight (CW) in grams; the Table 2 equations normalise by the cohort median 2680 g, encoded here as WT / 2.68 with WT in kg. Table 1: median 2680 g (range 980-5310).",
      source_name = "CW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL through F_age, centred on the cohort median 36.5 weeks (Table 2). Table 1: median 36.5 weeks (range 26.3-42.1).",
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Wu 2022 F_age uses PNA/12 with 12 = the cohort median PNA of 12 days (Table 2 footnote; Table 1 PNA in days, median 12.0, range 1-113). model() converts the canonical months to days with PNA * 30.4375 before forming the ratio. The Table 2 footnote's 'postnatal age in weeks' is inconsistent with its own '12 days' median and with Table 1; days is the only reading under which the centring value is the cohort median.",
      source_name = "PNA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    n_concentrations = 110L,
    age_range = "PNA 1-113 days (median 12); PMA 27.4-46.1 weeks (median 38.2)",
    age_median = "PNA 12 days; PMA 38.2 weeks",
    ga_range = "26.3-42.1 weeks (median 36.5); 43 with GA <= 37 weeks and 35 with GA > 37 weeks",
    weight_range = "Current weight 980-5310 g (median 2680); birth weight 480-4100 g (median 2500)",
    weight_median = "2.68 kg",
    race_ethnicity = "Chinese",
    disease_state = "Neonates and young infants (PMA < 48 weeks) with bacterial infection treated with meropenem",
    dose_range = "Intravenous meropenem 20-160 mg per dose (median 55 mg), 9.3-40.4 mg/kg per dose (median 19.7 mg/kg)",
    regions = "China (Children's Hospital of Hebei Province, Shijiazhuang; Beijing Children's Hospital, Beijing)",
    renal_function = "Serum creatinine median 23.8 umol/L (range 1.1-145.2); BUN median 3.79 mmol/L (range 0.5-19.6); albumin median 25.6 g/L (range 16.7-44.6)",
    notes = "Baseline characteristics from Wu 2022 Table 1. Opportunistic sampling from scavenged routine blood samples; NONMEM 7.4 FOCE-I. External validation in 16 newborns from a third site. Sex distribution not reported."
  )

  ini({
    # Structural parameters (Wu 2022 Table 2, 'Full Dataset Final Estimate').
    # Reference subject: CW = 2680 g, PNA = 12 days, GA = 36.5 weeks.
    lvc <- log(1.63)
    label("Typical V at CW 2.68 kg (L)") # Table 2: theta1 = 1.63 (RSE 11.4%)
    lcl <- log(0.503)
    label("Typical CL at CW 2.68 kg, PNA 12 days, GA 36.5 weeks (L/h)") # Table 2: theta2 = 0.503 (RSE 7.90%)

    # Allometric exponents fixed a priori (Results 3.2: 'allometric
    # coefficients of 1 for V and 0.75 for CL').
    e_wt_vc <- fixed(1)
    label("Allometric exponent of current weight on V (unitless)") # Results 3.2; Table 2 V equation
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of current weight on CL (unitless)") # Results 3.2; Table 2 CL equation

    # Maturation factor F_age = (PNA/12)^theta3 * (GA/36.5)^theta4 (Table 2)
    e_pna_cl <- 0.209
    label("Power exponent of postnatal age on CL (unitless)") # Table 2: theta3 = 0.209 (RSE 21.4%)
    e_ga_cl <- 2.14
    label("Power exponent of gestational age on CL (unitless)") # Table 2: theta4 = 2.14 (RSE 22.9%)

    # IIV (Table 2 'Inter-individual variability (%)', exponential model;
    # omega^2 = log(CV^2 + 1)).
    etalvc ~ 0.03017 # Table 2: V IIV 17.5% -> log(0.175^2 + 1)
    etalcl ~ 0.22957 # Table 2: CL IIV 50.8% -> log(0.508^2 + 1)

    # Residual error (Results 3.2: proportional model; Table 2 28.8%).
    propSd <- 0.288
    label("Proportional residual error (fraction)") # Table 2: residual variability 28.8% (RSE 15.5%)
  })
  model({
    # Postnatal age in days (canonical PNA is in months) for the Table 2
    # F_age term, which is centred on the 12-day cohort median.
    pna_days <- PNA * 30.4375
    f_age <- (pna_days / 12)^e_pna_cl * (GA / 36.5)^e_ga_cl

    vc <- exp(lvc + etalvc) * (WT / 2.68)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 2.68)^e_wt_cl * f_age

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
