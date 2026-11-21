Kong_2025_atorvastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral microdose atorvastatin (mixed CYP3A4/OATP1B/P-gp/BCRP substrate)",
    "in Chinese healthy volunteers and haemodialysis-dependent end-stage",
    "renal disease patients (Kong 2025). No covariates were retained."
  )
  reference <- paste(
    "Kong W, Pan Y, Wu Y, Hu Y, Jiang Z, Tian X, Bi S, Wang S, Feng F, Jin Y,",
    "Li J, Li H, Wang Y, Liang H, Tang W, Liu D.",
    "Microdose Cocktail Study Reveals the Activity and Key Influencing",
    "Factors of OATP1B, P-Gp, BCRP, and CYP3A in End-Stage Renal Disease",
    "Patients. Clin Pharmacol Ther. 2025;117(5):1303-.",
    "doi:10.1002/cpt.3546.",
    sep = " "
  )
  vignette <- "Kong_2025_microdose_cocktail_esrd"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  compartmentData <- list(
    depot = list(analyte = "atorvastatin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "atorvastatin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "atorvastatin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "18-65 years (inclusion criterion)",
    age_mean = "35.6 years (healthy volunteers); 43.4 years (ESRD)",
    sex_female_pct = 50.0,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "14 healthy volunteers and 10 ESRD patients on haemodialysis three",
      "times a week (chronic glomerulonephritis in all, no diabetes)."
    ),
    dose_range = paste(
      "Single oral dose of atorvastatin 100 ug in a microdose cocktail with",
      "dabigatran etexilate 375 ug, pitavastatin 10 ug, rosuvastatin 50 ug",
      "and midazolam 10 ug, fasted, 24 h before dialysis."
    ),
    regions = "China (Peking University Third Hospital, Beijing)",
    renal_function = "Healthy volunteers: mean eGFR 106.6 mL/min/1.73 m^2; ESRD: mean eGFR 4 mL/min/1.73 m^2",
    notes = paste(
      "Kong 2025 Table 1. BMI 22.4 (HV) and 21.9 (ESRD) kg/m^2; serum",
      "albumin 44.9 and 42.1 g/L. Female 9/14 HV and 3/10 ESRD. Samples at",
      "0.5-24 h (HV) and 1, 4 and 24 h pre-dialysis (ESRD) were pooled.",
      "One ESRD patient was excluded from the atorvastatin analysis because",
      "she had taken atorvastatin within the preceding 14 days. NCT05747768."
    )
  )

  ini({
    lka <- log(7.53); label("Absorption rate constant (1/h)") # Table 3 ATV ka 7.53 (RSE 16.7%)
    lcl <- log(363); label("Apparent clearance CL/F (L/h)") # Table 3 ATV CL/F 363 (RSE 9.20%)
    lvc <- log(835); label("Apparent central volume Vc/F (L)") # Table 3 ATV Vc/F 835 (RSE 37.2%)
    lvp <- log(5830); label("Apparent peripheral volume Vp/F (L)") # Table 3 ATV Vp/F 5,830 (RSE 2.10%)
    lq <- log(3920); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 ATV Q/F 3,920 (RSE 3.80%)

    etalcl ~ log(1 + 0.480^2) # Table 3 ATV CL/F IIV 48.0% (RSE 44.6%); omega^2 = log(1 + CV^2)
    etalvc ~ log(1 + 0.993^2) # Table 3 ATV Vc/F IIV 99.3% (RSE 63.3%); omega^2 = log(1 + CV^2)
    etalvp ~ log(1 + 0.626^2) # Table 3 ATV Vp/F IIV 62.6% (RSE 19.1%); omega^2 = log(1 + CV^2)
    etalq ~ log(1 + 0.302^2) # Table 3 ATV Q/F IIV 30.2% (RSE 84.7%); omega^2 = log(1 + CV^2)

    propSd <- fixed(0); label("Proportional residual error (fraction); not reported in the source") # residual error model and magnitude not reported in Kong 2025
  })
  model({
    # Eq. 5: CL/F = TVCL x e^eta (no covariates)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ug / L = ng/mL; x 1000 -> pg/mL (the paper's concentration unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
