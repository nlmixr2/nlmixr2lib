Li_2026_dabigatran <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "dabigatran after an oral microdose of dabigatran etexilate (intestinal",
    "P-gp probe) in Chinese healthy volunteers and type 2 diabetes patients",
    "(Li 2026). CL/F falls with serum urea (power model centred on 4.5",
    "mmol/L) and the relative bioavailability is 1.23-fold higher in type 2",
    "diabetes."
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
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Li 2026 reports serum urea in mmol/L (Table 1: 4.63 +/- 0.99 in",
        "healthy volunteers, 4.72 +/- 1.60 in type 2 diabetes). Enters as",
        "(Urea / 4.5)^-0.749 on CL/F (Eq. 2)."
      ),
      source_name = "Urea"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus indicator (1 = type 2 diabetes, 0 = healthy volunteer)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "All patients in Li 2026 had type 2 diabetes. Multiplies the",
        "relative bioavailability by 1.23 (Eq. 3)."
      ),
      source_name = "type 2 diabetes"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dabigatran etexilate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dabigatran", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dabigatran", units = "ug", specimen = "plasma", verified = TRUE)
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
      "Single oral dose of dabigatran etexilate 375 ug in a fasted microdose",
      "cocktail with midazolam 10 ug, pitavastatin 10 ug, rosuvastatin 50 ug",
      "and atorvastatin 100 ug. Dabigatran (the active moiety) was assayed."
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
    lka <- log(1.37); label("Absorption rate constant (1/h)") # Table S4 DAB Ka 1.37 (RSE 20.6%)
    lcl <- log(150); label("Apparent clearance CL/F at urea 4.5 mmol/L (L/h)") # Table S4 DAB CL/F 150 (RSE 6.10%); Eq. 2
    lvc <- log(865); label("Apparent central volume Vc/F (L)") # Table S4 DAB VC/F 865 (RSE 14.1%)
    lvp <- log(535); label("Apparent peripheral volume Vp/F (L)") # Table S4 DAB VT/F 535 (RSE 25.4%)
    lq <- log(66.2); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 DAB Q/F 66.2 (RSE 37.9%)

    e_bun_cl <- -0.749; label("Power exponent of urea/4.5 on CL/F (unitless)") # Table S4 'Urea on CL/F' -0.749 (RSE 21.5%); Eq. 2
    e_dis_diab_fdepot <- 1.23; label("Relative bioavailability in type 2 diabetes vs healthy volunteers (ratio)") # Table S4 'Type 2 diabetes on F' 1.23 (RSE 7.2%); Eq. 3

    etalcl ~ 0.0459 # Table S4 DAB omega^2 (CL/F) 0.0459 (RSE 35.1%)
    etalvc ~ 0.0556 # Table S4 DAB omega^2 (VC/F) 0.0556 (RSE 40.8%)

    propSd <- sqrt(0.0391); label("Proportional residual error (fraction)") # Table S4 DAB sigma^2 0.0391 (RSE 16.7%); SD = sqrt(variance)
  })
  model({
    # Eq. 2: CL/F = 150 x (Urea / 4.5)^-0.749
    urea_cl <- (BUN / 4.5)^e_bun_cl
    # Eq. 3: F = 1 in healthy volunteers, 1.23 in type 2 diabetes
    fdepot <- e_dis_diab_fdepot^DIS_DIAB

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * urea_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Dose is ug of dabigatran etexilate; CL/F, Vc/F, Vp/F and Q/F are
    # apparent values on that dose scale. ug / L = ng/mL; x 1000 -> pg/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
