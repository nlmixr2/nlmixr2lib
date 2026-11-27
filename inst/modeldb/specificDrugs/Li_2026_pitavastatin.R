Li_2026_pitavastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with parallel first-order and",
    "zero-order absorption for oral microdose pitavastatin (OATP1B probe) in",
    "Chinese healthy volunteers and type 2 diabetes patients (Li 2026). The",
    "dose splits into a first-order fraction (41%) and a zero-order fraction",
    "infused into the central compartment over 0.654 h; CL/F is 1.42-fold",
    "higher in men. Type 2 diabetes is not a covariate. Each administration",
    "needs two dose records with the same amount: one to depot and one to",
    "central with rate = -2."
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
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = paste(
        "Li 2026 Eq. 4 gives CL/F = 11.7 L/h in women and 11.7 x 1.42 in",
        "men, so the male factor is applied as 1.42^(1 - SEXF)."
      ),
      source_name = "Male (Eq. 4); Table S4 'Male on CL/F'"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pitavastatin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pitavastatin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pitavastatin", units = "ug", specimen = "plasma", verified = TRUE)
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
      "Single oral dose of pitavastatin 10 ug in a fasted microdose cocktail",
      "with midazolam 10 ug, dabigatran etexilate 375 ug, rosuvastatin 50 ug",
      "and atorvastatin 100 ug."
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
    lka <- log(0.230); label("First-order absorption rate constant of the depot fraction (1/h)") # Table S4 PTV Ka 0.230 (RSE 13.1%); ALAG1 0 FIX, so no lag
    lcl <- log(11.7); label("Apparent clearance CL/F in women (L/h)") # Table S4 PTV CL/F 11.7 (RSE 23.5%); Eq. 4
    lvc <- log(12.2); label("Apparent central volume Vc/F (L)") # Table S4 PTV VC/F 12.2 (RSE 16.1%)
    lvp <- log(378); label("Apparent peripheral volume Vp/F (L)") # Table S4 PTV VT/F 378 (RSE 60.8%)
    lq <- log(10.9); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 PTV Q/F 10.9 (RSE 18.2%)
    lfdepot <- log(0.410); label("Fraction of the dose absorbed first-order through the depot, F1 (fraction)") # Table S4 PTV F1 0.410 (RSE 6.60%)
    ld2 <- log(0.654); label("Duration of the zero-order input into central, D2 (h)") # Table S4 PTV D2 0.654 (RSE 6.90%)

    e_male_cl <- 1.42; label("CL/F ratio of men to women (ratio)") # Table S4 'Male on CL/F' 1.42 (RSE 14.5%); Eq. 4

    etalcl ~ 0.078 # Table S4 PTV omega^2 (CL/F) 0.078 (RSE 54.9%)
    etalvc ~ 0.412 # Table S4 PTV omega^2 (VC/F) 0.412 (RSE 33.7%)

    propSd <- sqrt(0.0366); label("Proportional residual error (fraction)") # Table S4 PTV sigma^2 0.0366 (RSE 16.6%); SD = sqrt(variance)
  })
  model({
    # Eq. 4: CL/F = 11.7 in women, 11.7 x 1.42 in men
    sex_cl <- e_male_cl^(1 - SEXF)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * sex_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    fdepot <- exp(lfdepot)
    d2 <- exp(ld2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Figure 2b: DOSE x F1 is absorbed first-order from depot and
    # DOSE x (1 - F1) enters central zero-order over D2, in parallel. Each
    # administration is two dose records of the full amount (depot, and
    # central with rate = -2); f() splits them.
    f(depot) <- fdepot
    f(central) <- 1 - fdepot
    dur(central) <- d2

    # ug / L = ng/mL; x 1000 -> pg/mL (the paper's concentration unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
