Wang_2022_latamoxef <- function() {
  description <- "Two-compartment IV population PK model for total (R + S) latamoxef (moxalactam) in 145 Chinese children aged 0.08-10.58 years with bacterial infection (Wang 2022). All four structural parameters scale with body surface area normalised to the cohort median 0.39 m^2: exponent 1 (fixed) on V1 and V2, an estimated 1.49 on CL and 0.75 (fixed) on Q. Exponential IIV on V1 and CL, additive residual error. The R- and S-epimers were fitted as separate models in the same paper (modellib('Wang_2022_latamoxef_r'), modellib('Wang_2022_latamoxef_s'))."
  reference <- "Wang Y, Sun D, Mei Y, Wu S, Li X, Li S, Wang J, Gao L, Xu H, Tuo Y. Population Pharmacokinetics and Dosing Regimen Optimization of Latamoxef in Chinese Children. Pharmaceutics. 2022;14(5):1033. doi:10.3390/pharmaceutics14051033"
  vignette <- "Wang_2022_latamoxef"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "latamoxef (R + S)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "latamoxef (R + S)", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Wang 2022 Table 1: mean 0.41 (SD 0.14), median 0.39 (range 0.20-1.03) m^2. Reference value 0.39 m^2 (the cohort median; Table 2 Model II 'BSA/BSA_median') in every power term of the Results final-model equations. The BSA formula is not stated, but the cohort median height 68 cm and weight 8 kg give 0.389 m^2 by Mosteller, matching the tabulated median.",
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
    notes = "Demographics from Wang 2022 Table 1 (91 male / 54 female). Sparse TDM sampling: 1-3 residual serum samples per child, total latamoxef 1.84-117.88 ug/mL by chiral HPLC-UV (LLOQ 1.5 ug/mL), which also resolves the R- and S-epimers. Model fitted in Phoenix NLME 8.2. Covariates screened by stepwise forward inclusion / backward elimination (sex, age, weight, height, BSA, haematology, renal and hepatic markers, hs-CRP, procalcitonin); only BSA was retained."
  )

  ini({
    # Structural parameters (Wang 2022 Table 3, group R + S, 'Final Model
    # Estimate' column; Results final-model equations). Reference BSA 0.39 m^2.
    lvc <- log(4.84)
    label("Central volume of distribution V1 at BSA = 0.39 m^2 (L)") # Table 3: theta_V1 = 4.84 L (SE 15.85%)
    lvp <- log(16.18)
    label("Peripheral volume of distribution V2 at BSA = 0.39 m^2 (L)") # Table 3: theta_V2 = 16.18 L (SE 47.35%)
    lcl <- log(1.00)
    label("Clearance CL at BSA = 0.39 m^2 (L/h)") # Table 3: theta_CL = 1.00 L/h (SE 9.05%)
    lq <- log(0.97)
    label("Inter-compartmental clearance Q at BSA = 0.39 m^2 (L/h)") # Table 3: theta_Q = 0.97 L/h (SE 15.93%)

    # BSA power exponents (Table 3 theta1-theta4; Results equations
    # V1 = 4.84 * (BSA/0.39), V2 = 16.18 * (BSA/0.39), CL = 1.00 * (BSA/0.39)^1.49,
    # Q = 0.97 * (BSA/0.39)^0.75).
    e_bsa_vc <- fixed(1)
    label("BSA power exponent on V1 (unitless)") # Table 3: theta1 = 1.00 (fixed)
    e_bsa_vp <- fixed(1)
    label("BSA power exponent on V2 (unitless)") # Table 3: theta2 = 1.00 (fixed)
    e_bsa_cl <- 1.49
    label("BSA power exponent on CL (unitless)") # Table 3: theta3 = 1.49 (SE 14.69%)
    e_bsa_q <- fixed(0.75)
    label("BSA power exponent on Q (unitless)") # Table 3: theta4 = 0.75 (fixed)

    # IIV, exponential model P_i = theta * exp(eta_i) (Methods Eq. 1). Table 3
    # prints omega in percent and its footnote defines omega as the 'square root
    # of inter-individual variance', so variance = (omega / 100)^2.
    etalvc ~ 1.10334 # Table 3: omega_V1 = 105.04 %; 1.0504^2
    etalcl ~ 0.08317 # Table 3: omega_CL = 28.84 %; 0.2884^2

    # Additive residual error (Methods Eq. 2).
    addSd <- 7.29
    label("Additive residual error (mg/L)") # Table 3: sigma = 7.29 mg/L (SE 11.49%)
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

    # Dose in mg, vc in L -> mg/L (= ug/mL, the paper's concentration unit).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
