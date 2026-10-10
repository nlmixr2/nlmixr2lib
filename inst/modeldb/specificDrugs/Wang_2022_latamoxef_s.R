Wang_2022_latamoxef_s <- function() {
  description <- "Two-compartment IV population PK model for the S-epimer of latamoxef (moxalactam) in 145 Chinese children aged 0.08-10.58 years with bacterial infection (Wang 2022). The dose is the TOTAL latamoxef dose: because the R-epimer fraction r of the administered product is unknown (the paper's stated range is 0.4444-0.5833, from the Chinese Pharmacopoeia limits on the R:S ratio), the paper estimates apparent parameters V1/(1 - r), V2/(1 - r), CL/(1 - r) and Q/(1 - r), so Cc is the S-epimer serum concentration and the epimer-specific V1, V2, CL and Q are (1 - r) times the model values. BSA power scaling normalised to 0.39 m^2: exponent 1 (fixed) on V1/(1 - r) and V2/(1 - r), an estimated 1.33 on CL/(1 - r) and 0.75 (fixed) on Q/(1 - r). Exponential IIV on V1/(1 - r) and CL/(1 - r), additive residual error. Total latamoxef and the R-epimer are fitted as separate models in the same paper (modellib('Wang_2022_latamoxef'), modellib('Wang_2022_latamoxef_r'))."
  reference <- "Wang Y, Sun D, Mei Y, Wu S, Li X, Li S, Wang J, Gao L, Xu H, Tuo Y. Population Pharmacokinetics and Dosing Regimen Optimization of Latamoxef in Chinese Children. Pharmaceutics. 2022;14(5):1033. doi:10.3390/pharmaceutics14051033"
  vignette <- "Wang_2022_latamoxef"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "latamoxef S-epimer (state carries S-epimer amount / (1 - r), i.e. total-latamoxef-dose equivalents)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "latamoxef S-epimer (state carries S-epimer amount / (1 - r), i.e. total-latamoxef-dose equivalents)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Wang 2022 Table 1: mean 0.41 (SD 0.14), median 0.39 (range 0.20-1.03) m^2. Reference value 0.39 m^2 (the cohort median) in every power term of the Results S-epimer final-model equations. The BSA formula is not stated, but the cohort median height 68 cm and weight 8 kg give 0.389 m^2 by Mosteller, matching the tabulated median.",
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 145L,
    n_studies = 1L,
    n_concentrations = 165L,
    age_range = "0.08-10.58 years",
    age_median = "0.60 years (mean 1.08, SD 1.63)",
    weight_range = "2.9-27.5 kg",
    weight_median = "8 kg (mean 8.68, SD 4.11)",
    height_range = "49-140 cm (median 68)",
    bsa_range = "0.20-1.03 m^2 (median 0.39)",
    sex_female_pct = 37.2,
    race_ethnicity = "Chinese (single-centre cohort, Wuhan Children's Hospital)",
    disease_state = "Hospitalised children with bacterial infection treated with latamoxef (July-November 2021).",
    dose_range = "Latamoxef sodium 40-80 mg/kg/day by intravenous injection, divided into two or three doses.",
    regions = "China (Wuhan Children's Hospital, Tongji Medical College, Huazhong University of Science and Technology)",
    renal_function = "Modified-Schwartz eGFR median 123.76 (range 63.61-267.24) mL/min/1.73 m^2; serum creatinine median 20.90 (9.70-48.10) umol/L. No child had eGFR < 60.",
    notes = "Demographics from Wang 2022 Table 1 (91 male / 54 female). Sparse TDM sampling: 1-3 residual serum samples per child; S-epimer 0.06-62.97 ug/mL by chiral HPLC-UV (Chinese Pharmacopoeia 2020 method, LLOQ 1.5 ug/mL for latamoxef). Model fitted in Phoenix NLME 8.2; only BSA was retained as a covariate."
  )

  ini({
    # Apparent structural parameters relative to the TOTAL latamoxef dose
    # (Wang 2022 Table 3, group S, 'Final Model Estimate' column; Results
    # S-epimer final-model equations). Reference BSA 0.39 m^2. The
    # epimer-specific values are (1 - r) times these: Table 3 ranges
    # theta_V1-S 3.38-4.51 L, theta_V2-S 7.97-10.63 L, theta_CL-S 0.98-1.31 L/h
    # and theta_Q-S 0.79-1.05 L/h for r = 0.4444-0.5833.
    lvc <- log(8.12)
    label("Apparent central volume V1/(1 - r) at BSA = 0.39 m^2 (L)") # Table 3: theta_V1/(1-r) = 8.12 L (SE 18.13%)
    lvp <- log(19.13)
    label("Apparent peripheral volume V2/(1 - r) at BSA = 0.39 m^2 (L)") # Table 3: theta_V2/(1-r) = 19.13 L (SE 48.39%)
    lcl <- log(2.36)
    label("Apparent clearance CL/(1 - r) at BSA = 0.39 m^2 (L/h)") # Table 3: theta_CL/(1-r) = 2.36 L/h (SE 9.65%)
    lq <- log(1.89)
    label("Apparent inter-compartmental clearance Q/(1 - r) at BSA = 0.39 m^2 (L/h)") # Table 3: theta_Q/(1-r) = 1.89 L/h (SE 20.26%)

    # BSA power exponents (Table 3 theta9-theta12; Results equations
    # V1/(1-r) = 8.12 * (BSA/0.39), V2/(1-r) = 19.13 * (BSA/0.39),
    # CL/(1-r) = 2.36 * (BSA/0.39)^1.33, Q/(1-r) = 1.89 * (BSA/0.39)^0.75).
    e_bsa_vc <- fixed(1)
    label("BSA power exponent on V1/(1 - r) (unitless)") # Table 3: theta9 = 1.00 (fixed)
    e_bsa_vp <- fixed(1)
    label("BSA power exponent on V2/(1 - r) (unitless)") # Table 3: theta10 = 1.00 (fixed)
    e_bsa_cl <- 1.33
    label("BSA power exponent on CL/(1 - r) (unitless)") # Table 3: theta11 = 1.33 (SE 20.48%)
    e_bsa_q <- fixed(0.75)
    label("BSA power exponent on Q/(1 - r) (unitless)") # Table 3: theta12 = 0.75 (fixed)

    # IIV, exponential model (Methods Eq. 1). Table 3 prints omega in percent
    # and defines it as the 'square root of inter-individual variance', so
    # variance = (omega / 100)^2.
    etalvc ~ 1.35024 # Table 3: omega_V1/(1-r) = 116.20 %; 1.1620^2
    etalcl ~ 0.18662 # Table 3: omega_CL/(1-r) = 43.20 %; 0.4320^2

    # Additive residual error (Methods Eq. 2).
    addSd <- 3.81
    label("Additive residual error (mg/L)") # Table 3: sigma_S = 3.81 mg/L (SE 9.24%)
  })
  model({
    vc <- exp(lvc + etalvc) * (BSA / 0.39)^e_bsa_vc
    vp <- exp(lvp) * (BSA / 0.39)^e_bsa_vp
    cl <- exp(lcl + etalcl) * (BSA / 0.39)^e_bsa_cl
    q <- exp(lq) * (BSA / 0.39)^e_bsa_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose = total latamoxef (mg); vc = V1/(1 - r) (L) -> S-epimer concentration (mg/L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
