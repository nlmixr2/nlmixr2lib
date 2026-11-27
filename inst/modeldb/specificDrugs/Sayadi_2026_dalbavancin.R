Sayadi_2026_dalbavancin <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous dalbavancin in adults with Gram-positive infections monitored by routine therapeutic drug monitoring (Limoges University Hospital, France). Clearance scales with CKD-EPI eGFR (power 0.22) and total body weight (power 0.61), both centred on the cohort medians (75 mL/min/1.73 m^2, 79 kg); the central and peripheral volumes scale linearly with body weight (exponent fixed at 1). Unbound concentration is returned as Ccu using the 7 percent free fraction the paper applied for its fAUC24/MIC >= 111.1 target-attainment simulations."
  reference <- paste(
    "Sayadi H, Gregoire M, Fromage Y, Labriffe M, Maizaud F, Monchaud C,",
    "Codde C, Boutoille D, Marquet P, Woillard J-B. Population",
    "pharmacokinetics of dalbavancin: external validation, model averaging,",
    "and implications for precision dosing in prolonged therapy. Antimicrob",
    "Agents Chemother. 2026;70(8):e00665-26. doi:10.1128/aac.00665-26.",
    sep = " "
  )
  vignette <- "Sayadi_2026_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (estimated exponent 0.61) and on V1 and V2",
        "(exponents fixed to 1), centred on the cohort median of 79 kg",
        "(Sayadi 2026 Table 2; Methods 'Continuous covariates were centered on",
        "their population median'). Cohort median 79 kg, IQR 65-88, range",
        "32-138 kg. Missing weights (27 of 183 patients) were imputed by",
        "k-nearest neighbours (k = 5) before fitting."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the creatinine CKD-EPI equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL (exponent 0.22) centred on the cohort median of",
        "75 mL/min/1.73 m^2 (Sayadi 2026 Table 2). The paper's eGFR is the",
        "CKD-EPI estimate (Methods, Study design and data source), which is",
        "BSA-normalised; it is not a Cockcroft-Gault creatinine clearance.",
        "Cohort median 75, IQR 49-87, range 14-135 mL/min/1.73 m^2. Missing",
        "values (29 of 183 patients) were imputed by k-nearest neighbours."
      ),
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Candidate covariate (Sayadi 2026 Methods, Covariate analysis) not retained. Cohort median 77 years, IQR 67-83, range 23-99 (Table 2)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Candidate covariate not retained. 37.2% female (Table 2)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Candidate covariate not retained. Cohort median 25 kg/m^2, IQR 23-30 (Table 2)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Candidate covariate not retained (eGFR was retained instead). Cohort median 87 umol/L, IQR 65-110 (Table 2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 1L,
    n_observations = 406L,
    age_range = "23-99 years",
    age_median = "77 years",
    weight_range = "32-138 kg",
    weight_median = "79 kg",
    sex_female_pct = 37.2,
    race_ethnicity = "Not reported; single-centre French cohort.",
    disease_state = paste(
      "Adults treated with dalbavancin for Gram-positive infections:",
      "osteoarticular infection 43.2%, infective endocarditis 12.6%, vascular",
      "prosthetic infection 6.0%, bacteraemia 2.2%, ABSSSI 1.1%, other or not",
      "specified 35.0% (Table 2). Pathogens: coagulase-negative staphylococci",
      "30.1%, Staphylococcus aureus 19.7%, enterococci 7.1%."
    ),
    dose_range = "Routine-care intravenous dalbavancin (typically 1500 mg); 1 to 23 doses per patient, 24.6% received 3 or more doses (Table 2). Infusion duration as recorded, or 30 min when missing.",
    regions = "France (Limoges University Hospital); external validation cohort from Nantes University Hospital (n = 30).",
    renal_function = "CKD-EPI eGFR median 75 mL/min/1.73 m^2, IQR 49-87, range 14-135 (Table 2).",
    notes = paste(
      "Retrospective observational TDM study, March 2021 to March 2026; 1-25",
      "samples per patient, mostly troughs and late elimination phase. Total",
      "plasma dalbavancin by LC-MS/MS (LLOQ 0.5 mg/L), all concentrations",
      "above LLOQ. Estimation by SAEM in Monolix 2024R1. Externally",
      "validated in 30 Nantes patients (92 concentrations)."
    )
  )

  ini({
    # Structural parameters -- Sayadi 2026 Table 4 (final model,
    # stochastic approximation), typical values for a 79 kg patient with
    # eGFR 75 mL/min/1.73 m^2.
    lcl <- log(0.035); label("Clearance CL (L/h)") # Sayadi 2026 Table 4: TV_Cl = 0.035 L/h (RSE 2.99%)
    lvc <- log(4.74); label("Central volume of distribution V1 (L)") # Sayadi 2026 Table 4: TV_V1 = 4.74 L (RSE 4.86%)
    lq <- log(0.033); label("Intercompartmental clearance Q (L/h)") # Sayadi 2026 Table 4: TV_Q = 0.033 L/h (RSE 13.66%)
    lvp <- log(8.74); label("Peripheral volume of distribution V2 (L)") # Sayadi 2026 Table 4: TV_V2 = 8.74 L (RSE 7.12%)

    # Covariate power exponents -- Sayadi 2026 Table 4. The volume
    # exponents were fixed to 1 (Results, Structural model and covariates;
    # Table 4 footnote).
    e_crcl_cl <- 0.22; label("Power exponent of eGFR on CL (unitless)") # Sayadi 2026 Table 4: beta eGFR,CL = 0.22 (RSE 19.92%)
    e_wt_cl <- 0.61; label("Power exponent of body weight on CL (unitless)") # Sayadi 2026 Table 4: beta WT,CL = 0.61 (RSE 12.67%)
    e_wt_vc <- fixed(1); label("Power exponent of body weight on V1 (unitless)") # Sayadi 2026 Table 4: beta WT,V1 = 1, footnote 'fixed to 1.0'
    e_wt_vp <- fixed(1); label("Power exponent of body weight on V2 (unitless)") # Sayadi 2026 Table 4: beta WT,V2 = 1, footnote 'fixed to 1.0'

    # Plasma unbound fraction, not estimated: the paper converts total to
    # unbound concentration with fu = 0.07 for its fAUC24/MIC target
    # (Methods, Monte Carlo simulations and probability of target attainment).
    fu <- fixed(0.07); label("Dalbavancin plasma unbound fraction (unitless)") # Sayadi 2026 Methods: 'fraction unbound (fu) of 0.07'

    # Inter-individual variability. Table 4 reports omega as the standard
    # deviation of the log-normal random effect; variance = omega^2.
    etalcl ~ 0.0289 # Sayadi 2026 Table 4: omega CL = 0.17 (CV 17.5%); 0.17^2
    etalvc ~ 0.04 # Sayadi 2026 Table 4: omega V1 = 0.2 (CV 20.02%); 0.2^2
    etalq ~ 0.2401 # Sayadi 2026 Table 4: omega Q = 0.49 (CV 52.39%); 0.49^2
    etalvp ~ 0.1156 # Sayadi 2026 Table 4: omega V2 = 0.34 (CV 35.22%); 0.34^2

    # Residual error -- proportional only.
    propSd <- 0.21; label("Proportional residual error (fraction)") # Sayadi 2026 Table 4: sigma prop = 0.21 (RSE 5.31%)
  })

  model({
    # Individual parameters. Continuous covariates enter as power functions
    # centred on the cohort medians (Table 2: WT 79 kg, eGFR 75).
    cl <- exp(lcl + etalcl) * (CRCL / 75)^e_crcl_cl * (WT / 79)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 79)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 79)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion into the central compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Cc is total plasma dalbavancin (the assayed and fitted quantity);
    # Ccu is unbound dalbavancin, which drives the paper's fAUC24/MIC target.
    Cc <- central / vc
    Ccu <- Cc * fu

    Cc ~ prop(propSd)
  })
}
