Li_2026_atorvastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with zero-order absorption for",
    "oral microdose atorvastatin (OATP1B / BCRP / P-gp / CYP3A probe) in",
    "Chinese healthy volunteers and type 2 diabetes patients (Li 2026). The",
    "dose enters the central compartment zero-order over 0.5 h, and the",
    "relative bioavailability is 1.73-fold higher in type 2 diabetes. Dose",
    "records go to central with rate = -2."
  )
  reference <- paste(
    "Li Y, Yang J, Zhou E, Geng K, Lei Z, Yao X, Feng F, Jin Y, Kong W, Hu Y,",
    "Liang H, Li H, Yang G-P, Hong T, Liu D.",
    "A microdose cocktail and population pharmacokinetic study suggests",
    "potentially reduced CYP3A and P-gp activities in type 2 diabetes.",
    "Front Pharmacol. 2026;17:1924307. doi:10.3389/fphar.2026.1924307.",
    sep = " "
  )
  vignette <- "Li_2026_microdose_cocktail_t2d"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus indicator (1 = type 2 diabetes, 0 = healthy volunteer)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "All patients in Li 2026 had type 2 diabetes. Multiplies the",
        "relative bioavailability by 1.73 (Eq. 6)."
      ),
      source_name = "type 2 diabetes"
    )
  )

  compartmentData <- list(
    central = list(analyte = "atorvastatin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "atorvastatin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 2L,
    age_range = "18-65 years (inclusion criterion)",
    age_mean = "35.6 years (healthy volunteers); 49.2 years (type 2 diabetes)",
    weight_mean = "58.8 kg (healthy volunteers); 67.1 kg (type 2 diabetes)",
    sex_female_pct = 66.7,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "14 healthy volunteers and 13 type 2 diabetes patients (HbA1c 7.62 +/-",
      "0.77%, disease duration 1.5-8.5 years, normal hepatic and renal",
      "function; 7 on metformin, 2 acarbose, 2 glyburide, 1 dapagliflozin)."
    ),
    dose_range = paste(
      "Single oral dose of atorvastatin 100 ug in a fasted microdose cocktail",
      "with midazolam 10 ug, dabigatran etexilate 375 ug, pitavastatin 10 ug",
      "and rosuvastatin 50 ug."
    ),
    regions = "China (Peking University Third Hospital, Beijing)",
    renal_function = "Mean eGFR 106.6 (healthy volunteers) and 103.6 (type 2 diabetes) mL/min/1.73 m^2",
    notes = paste(
      "Li 2026 Table 1. Healthy-volunteer data come from the earlier Kong",
      "2025 cocktail study and were sampled at 0.5-24 h; type 2 diabetes",
      "patients were sampled at 1, 4 and 24 h. NCT06273215."
    )
  )

  ini({
    lcl <- log(416); label("Apparent clearance CL/F (L/h)") # Table S4 ATV CL/F 416 (RSE 11.0%)
    lvc <- log(4250); label("Apparent central volume Vc/F (L)") # Table S4 ATV VC/F 4250 (RSE 17.6%)
    lvp <- log(6050); label("Apparent peripheral volume Vp/F (L)") # Table S4 ATV VT/F 6050 (RSE 10.5%)
    lq <- log(3490); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 ATV Q/F 3490 (RSE 12.0%)
    ld2 <- log(0.500); label("Duration of the zero-order input into central, D2 (h)") # Table S4 ATV D2 0.500 (RSE 1.80%)

    e_dis_diab_fcentral <- 1.73; label("Relative bioavailability in type 2 diabetes vs healthy volunteers (ratio)") # Table S4 'Type 2 diabetes on F' 1.73 (RSE 9.00%); Eq. 6

    etalcl ~ 0.0502 # Table S4 ATV omega^2 (CL/F) 0.0502 (RSE 80.5%)
    etalvc ~ 0.541 # Table S4 ATV omega^2 (VC/F) 0.541 (RSE 36.6%)

    propSd <- sqrt(0.0878); label("Proportional residual error (fraction)") # Table S4 ATV sigma^2 0.0878 (RSE 17.0%); SD = sqrt(variance)
  })
  model({
    # Eq. 6: F = 1 in healthy volunteers, 1.73 in type 2 diabetes
    fcentral <- e_dis_diab_fcentral^DIS_DIAB

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    d2 <- exp(ld2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Figure 2d: the dose enters central zero-order over D2 (dose records
    # to central with rate = -2).
    f(central) <- fcentral
    dur(central) <- d2

    # ug / L = ng/mL; x 1000 -> pg/mL (the paper's concentration unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
