Kong_2025_midazolam <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral microdose midazolam (CYP3A probe) in Chinese healthy volunteers",
    "and haemodialysis-dependent end-stage renal disease (ESRD) patients",
    "(Kong 2025). CL/F is 94% higher in ESRD."
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

  covariateData <- list(
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator (1 = ESRD on maintenance haemodialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer, eGFR >= 90 mL/min/1.73 m^2)",
      notes = paste(
        "Kong 2025 'disease state (with or without ESRD)'. ESRD = eGFR <= 15",
        "mL/min/1.73 m^2 (CKD-EPI) on haemodialysis three times a week; the",
        "comparator group was healthy volunteers with eGFR >= 90. Enters as",
        "(1 + DIS_MDZ) on CL/F (Eq. 3). Only pre-dialysis samples were",
        "modelled, so the model carries no intradialytic clearance."
      ),
      source_name = "DIS"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
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
      "Single oral dose of midazolam 10 ug in a microdose cocktail with",
      "dabigatran etexilate 375 ug, pitavastatin 10 ug, rosuvastatin 50 ug",
      "and atorvastatin 100 ug, fasted, 24 h before dialysis."
    ),
    regions = "China (Peking University Third Hospital, Beijing)",
    renal_function = "Healthy volunteers: mean eGFR 106.6 mL/min/1.73 m^2; ESRD: mean eGFR 4 mL/min/1.73 m^2",
    notes = paste(
      "Kong 2025 Table 1. BMI 22.4 (HV) and 21.9 (ESRD) kg/m^2; serum",
      "albumin 44.9 and 42.1 g/L. Female 9/14 HV and 3/10 ESRD. Samples at",
      "0.5-24 h (HV) and 1, 4 and 24 h pre-dialysis (ESRD) were pooled.",
      "NCT05747768."
    )
  )

  ini({
    lka <- fixed(log(13)); label("Absorption rate constant (1/h)") # Table 3 MDZ ka 13 (RSE NA; bootstrap 13 (13, 13)) -> fixed
    lcl <- log(67.6); label("Apparent clearance CL/F in healthy volunteers (L/h)") # Table 3 MDZ CL/F 67.6 (RSE 24.6%)
    lvc <- log(195); label("Apparent central volume Vc/F (L)") # Table 3 MDZ Vc/F 195 (RSE 18.8%)
    lvp <- log(333); label("Apparent peripheral volume Vp/F (L)") # Table 3 MDZ Vp/F 333 (RSE 25.6%)
    lq <- log(34.7); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 MDZ Q/F 34.7 (RSE 15.5%)

    e_renalimp_esrd_cl <- 0.94; label("Fractional change in CL/F in ESRD (fraction)") # Table 3 'DISMDZ on CL/F' 0.94 (RSE 52.1%); Eq. 3 (1 + DIS_MDZ)

    etalcl ~ log(1 + 0.379^2) # Table 3 MDZ CL/F IIV 37.9% (RSE 27.3%); omega^2 = log(1 + CV^2)
    etalvc ~ log(1 + 0.421^2) # Table 3 MDZ Vc/F IIV 42.1% (RSE 43.2%); omega^2 = log(1 + CV^2)

    propSd <- fixed(0); label("Proportional residual error (fraction); not reported in the source") # residual error model and magnitude not reported in Kong 2025
  })
  model({
    # Eq. 3: CL/F = TVCL x e^eta x (1 + DIS_MDZ)
    esrd_cl <- 1 + e_renalimp_esrd_cl * RENALIMP_ESRD

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * esrd_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

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
