Wang_2022_latamoxef_r <- function() {
  description <- "Two-compartment IV population PK model for the R-epimer of latamoxef (moxalactam) in 145 Chinese children aged 0.08-10.58 years with bacterial infection (Wang 2022). The dose is the TOTAL latamoxef dose: because the R-epimer fraction r of the administered product is unknown (the paper's stated range is 0.4444-0.5833, from the Chinese Pharmacopoeia limits on the R:S ratio), the paper estimates apparent parameters V1/r, V2/r, CL/r and Q/r, so Cc is the R-epimer serum concentration and the epimer-specific V1, V2, CL and Q are r times the model values. BSA power scaling normalised to 0.39 m^2: exponent 1 (fixed) on V1/r and V2/r, an estimated 1.42 on CL/r and 0.75 (fixed) on Q/r. Exponential IIV on V1/r and CL/r, additive residual error. Total latamoxef and the S-epimer are fitted as separate models in the same paper (modellib('Wang_2022_latamoxef'), modellib('Wang_2022_latamoxef_s'))."
  reference <- "Wang Y, Sun D, Mei Y, Wu S, Li X, Li S, Wang J, Gao L, Xu H, Tuo Y. Population Pharmacokinetics and Dosing Regimen Optimization of Latamoxef in Chinese Children. Pharmaceutics. 2022;14(5):1033. doi:10.3390/pharmaceutics14051033"
  vignette <- "Wang_2022_latamoxef"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "latamoxef R-epimer (state carries R-epimer amount / r, i.e. total-latamoxef-dose equivalents)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "latamoxef R-epimer (state carries R-epimer amount / r, i.e. total-latamoxef-dose equivalents)",
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
      notes = "Wang 2022 Table 1: mean 0.41 (SD 0.14), median 0.39 (range 0.20-1.03) m^2. Reference value 0.39 m^2 (the cohort median) in every power term of the Results R-epimer final-model equations. The BSA formula is not stated, but the cohort median height 68 cm and weight 8 kg give 0.389 m^2 by Mosteller, matching the tabulated median.",
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
    notes = "Demographics from Wang 2022 Table 1 (91 male / 54 female). Sparse TDM sampling: 1-3 residual serum samples per child; R-epimer 0.19-60.92 ug/mL by chiral HPLC-UV (Chinese Pharmacopoeia 2020 method, LLOQ 1.5 ug/mL for latamoxef). Model fitted in Phoenix NLME 8.2; only BSA was retained as a covariate."
  )

  ini({
    # Apparent structural parameters relative to the TOTAL latamoxef dose
    # (Wang 2022 Table 3, group R, 'Final Model Estimate' column; Results
    # R-epimer final-model equations). Reference BSA 0.39 m^2. The
    # epimer-specific values are r times these: Table 3 ranges theta_V1-R
    # 4.31-5.65 L, theta_V2-R 14.67-19.25 L, theta_CL-R 0.75-0.98 L/h and
    # theta_Q-R 1.40-1.84 L/h for r = 0.4444-0.5833.
    lvc <- log(9.69)
    label("Apparent central volume V1/r at BSA = 0.39 m^2 (L)") # Table 3: theta_V1/r = 9.69 L (SE 16.00%)
    lvp <- log(33.00)
    label("Apparent peripheral volume V2/r at BSA = 0.39 m^2 (L)") # Table 3: theta_V2/r = 33.00 L (SE 33.75%)
    lcl <- log(1.68)
    label("Apparent clearance CL/r at BSA = 0.39 m^2 (L/h)") # Table 3: theta_CL/r = 1.68 L/h (SE 9.71%)
    lq <- log(3.15)
    label("Apparent inter-compartmental clearance Q/r at BSA = 0.39 m^2 (L/h)") # Table 3: theta_Q/r = 3.15 L/h (SE 21.70%)

    # BSA power exponents (Table 3 theta5-theta8; Results equations
    # V1/r = 9.69 * (BSA/0.39), V2/r = 33.00 * (BSA/0.39),
    # CL/r = 1.68 * (BSA/0.39)^1.42, Q/r = 3.15 * (BSA/0.39)^0.75).
    e_bsa_vc <- fixed(1)
    label("BSA power exponent on V1/r (unitless)") # Table 3: theta5 = 1.00 (fixed)
    e_bsa_vp <- fixed(1)
    label("BSA power exponent on V2/r (unitless)") # Table 3: theta6 = 1.00 (fixed)
    e_bsa_cl <- 1.42
    label("BSA power exponent on CL/r (unitless)") # Table 3: theta7 = 1.42 (SE 18.07%)
    e_bsa_q <- fixed(0.75)
    label("BSA power exponent on Q/r (unitless)") # Table 3: theta8 = 0.75 (fixed)

    # IIV, exponential model (Methods Eq. 1). Table 3 prints omega in percent
    # and defines it as the 'square root of inter-individual variance', so
    # variance = (omega / 100)^2.
    etalvc ~ 0.42393 # Table 3: omega_V1/r = 65.11 %; 0.6511^2
    etalcl ~ 0.12510 # Table 3: omega_CL/r = 35.37 %; 0.3537^2

    # Additive residual error (Methods Eq. 2).
    addSd <- 5.33
    label("Additive residual error (mg/L)") # Table 3: sigma_R = 5.33 mg/L (SE 13.07%)
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

    # Dose = total latamoxef (mg); vc = V1/r (L) -> R-epimer concentration (mg/L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
