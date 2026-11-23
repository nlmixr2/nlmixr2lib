Suzuoka_2026_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (2-5 g/m^2 over 24 h) in Japanese children with acute",
    "lymphoblastic leukaemia (Suzuoka 2026; n = 27 patients in the",
    "model-development cohort, each treatment course modelled as a separate",
    "individual). Linear first-order elimination from the central",
    "compartment. Body surface area and the delivered intravenous hydration",
    "rate per body surface area (time-varying) raise clearance through power",
    "terms; age and the hydration rate raise intercompartmental clearance.",
    "Exponential between-course variability on CL and the peripheral volume,",
    "and a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Suzuoka T, Noda T, Tsuge K, Ishiwata Y, Hirai T, Isoda T, Takagi M,",
    "Nagata M (2026). Development and Evaluation of a Population",
    "Pharmacokinetic Model for High-Dose Methotrexate in Japanese Pediatric",
    "Patients with Acute Lymphoblastic Leukemia.",
    "J Clin Pharmacol 66(10):e70310. doi:10.1002/jcph.70310.",
    sep = " "
  )
  vignette <- "Suzuoka_2026_methotrexate"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Two-compartment ADVAN3 TRANS4 model (Methods, 'Development of a PPK
  # Model'). Amounts are in umol so that Cc is in umol/L, the unit in which
  # the paper reports every methotrexate concentration and target threshold;
  # convert a dose in mg with the methotrexate molecular weight 454.44 g/mol
  # (umol = mg / 454.44 * 1000). Methods: 'Plasma methotrexate concentrations
  # were measured using a homogeneous enzyme immunoassay', so the specimen is
  # plasma.
  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area (DuBois)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column BSA, calculated with the DuBois equation (Methods:",
        "'BSA calculated using the Dubois equation', chosen because dose and",
        "hydration rate are prescribed on DuBois BSA at the institution).",
        "Power-form effect on clearance normalised to the development-cohort",
        "median 0.74 m^2 (Table 1; range 0.52-1.6), Table 2 / Eq 11:",
        "'CL (L/h) = 3.49 x (BSA/0.74)^0.357 x (HYDR/129)^0.107'. Collected",
        "per treatment course; each course is a separate individual in the",
        "analysis."
      ),
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column Age (years). Power-form effect on intercompartmental",
        "clearance only, normalised to the development-cohort median 6 years",
        "(Table 1; range 2-13), Table 2 / Eq 12:",
        "'Q (L/h) = 0.0714 x (Age/6)^0.753 x (HYDR/129)^0.176'."
      ),
      source_name = "Age"
    ),
    HYDRATION_RATE_BSA = list(
      description = "Delivered intravenous hydration fluid rate per body surface area",
      units = "mL/m^2/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column HYDR. Methods Eq 1: 'HYDR (mL/m2/h) = Hydration volume",
        "rate (mL/h) / BSA (m2)', where the rate is the delivered rate",
        "recorded in the medical records. Started at 125 mL/m^2/h with the",
        "methotrexate infusion and raised by the physician when concentrations",
        "exceeded institutional thresholds. TIME-VARYING: 'When the hydration",
        "rate was changed during a course, the updated HYDR value was applied",
        "to all records from the time of the change onward'; supply it on",
        "every record of the dataset, with the new value from the time of the",
        "change. Power-form effects on clearance (exponent 0.107) and on",
        "intercompartmental clearance (exponent 0.176), both normalised to the",
        "development-cohort median 129 mL/m^2/h (Table 1; range 104-187). The",
        "effect on CL was pre-specified on clinical grounds and kept although",
        "not statistically significant (Table 3; RSE 125%)."
      ),
      source_name = "HYDR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    age_range = "2 to 13 years",
    age_median = "6 years",
    weight_range = "11.4 to 57.9 kg",
    weight_median = "18.9 kg",
    height_range = "84.4 to 159.9 cm (median 108.2)",
    bsa_range = "0.52 to 1.6 m^2 (median 0.74)",
    sex_female_pct = 48.1,
    race_ethnicity = "Japanese (single-centre cohort in Tokyo).",
    disease_state = "Paediatric acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate with hyperhydration and urine alkalinisation.",
    renal_function = "Serum creatinine 0.15-0.51 mg/dL (median 0.24); Uemura eGFR 117-209 mL/min/1.73 m^2 (median 163).",
    hepatic_function = "AST 12-118 U/L (median 23); ALT 6-96 U/L (median 13).",
    dose_range = "1.88-5.08 g/m^2 (median 3.04; 1100-7800 mg): 10% of the dose over the first 30 min, the remaining 90% over 23.5 h.",
    hydration = "Delivered IV hydration rate 104-187 mL/m^2/h (median 129); initial rate 125 mL/m^2/h.",
    regions = "Japan (Tokyo Medical and Dental University Hospital, 2013-2024).",
    notes = paste(
      "Demographics from Table 1, model-development group (14 male / 13",
      "female). A further 11 patients formed the validation group (age 2-15",
      "years, BSA 0.53-1.68 m^2, HYDR 114-191 mL/m^2/h); 504 concentrations",
      "were measured across all 38 patients (4-6 per course, sampled 24-72 h",
      "after the start of the infusion). Every treatment course (1-4 per",
      "patient) was treated as an independent individual, so the random",
      "effects are between-course rather than between-patient. Concentrations",
      "below the 0.04 umol/L LLOQ were excluded (M1)."
    )
  )

  ini({
    # Structural parameters -- Table 2 'Final' column (also Table S1 and
    # Eqs 11-12).
    lcl <- log(3.49); label("Clearance at BSA 0.74 m^2 and hydration 129 mL/m^2/h (L/h)") # Table 2 theta1 = 3.49 (95% CI 3.04-3.94)
    lvc <- log(10.4); label("Central volume of distribution (L)") # Table 2 theta2 = 10.4 (95% CI 8.82-12.0)
    lq <- log(0.0714); label("Intercompartmental clearance at age 6 y and hydration 129 mL/m^2/h (L/h)") # Table 2 theta3 = 0.0714 (95% CI 0.0494-0.0934)
    lvp <- log(1.07); label("Peripheral volume of distribution (L)") # Table 2 theta4 = 1.07 (95% CI 0.778-1.36)

    # Covariate exponents (power models centred on the cohort medians).
    e_bsa_cl <- 0.357; label("Power exponent of BSA on CL (unitless)") # Table 2 theta5 = 0.357 (95% CI 0.229-0.485)
    e_hydration_rate_bsa_cl <- 0.107; label("Power exponent of hydration rate on CL (unitless)") # Table 2 theta6 = 0.107 (95% CI 0-0.370)
    e_age_q <- 0.753; label("Power exponent of age on Q (unitless)") # Table 2 theta7 = 0.753 (95% CI 0.397-1.11)
    e_hydration_rate_bsa_q <- 0.176; label("Power exponent of hydration rate on Q (unitless)") # Table 2 theta8 = 0.176 (95% CI 0.104-0.248)

    # Between-course variability, exponential. Table 2 reports CV%;
    # converted with omega^2 = log(1 + CV^2). Table 4 records the IIV on Vc
    # and Q as 'Fixed at 0', so neither carries an eta here.
    etalcl ~ 0.026541 # Table 2 omega CL = 16.4% CV -> log(1 + 0.164^2)
    etalvp ~ 0.019410 # Table 2 omega Vp = 14.0% CV -> log(1 + 0.140^2)

    # Residual error -- proportional (Results: 'A two-compartment open model
    # was developed using a proportional error model').
    propSd <- 0.312; label("Proportional residual error (fraction)") # Table 2 sigma = 31.2% CV
  })

  model({
    # Individual PK parameters, Table 2 / Eqs 11-12. HYDRATION_RATE_BSA is
    # time-varying, so CL and Q change from the time the hydration rate
    # changes.
    cl <- exp(lcl + etalcl) * (BSA / 0.74)^e_bsa_cl * (HYDRATION_RATE_BSA / 129)^e_hydration_rate_bsa_cl
    vc <- exp(lvc)
    q <- exp(lq) * (AGE / 6)^e_age_q * (HYDRATION_RATE_BSA / 129)^e_hydration_rate_bsa_q
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion straight into the central compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma methotrexate concentration in umol/L (amounts umol, volumes L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
