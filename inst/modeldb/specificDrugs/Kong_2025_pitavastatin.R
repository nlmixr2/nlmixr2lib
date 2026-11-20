Kong_2025_pitavastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral microdose pitavastatin (OATP1B probe) in Chinese healthy",
    "volunteers and haemodialysis-dependent end-stage renal disease (ESRD)",
    "patients (Kong 2025). CL/F is reduced by 38.1% in ESRD and falls further",
    "with the stool relative abundance of the gut genus Veillonella."
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
        "(1 + DIS_PTV) on CL/F (Eq. 1). Only pre-dialysis samples were",
        "modelled, so the model carries no intradialytic clearance."
      ),
      source_name = "DIS"
    ),
    GUT_VEILLONELLA_PCT = list(
      description = "Stool relative abundance of the gut genus Veillonella",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "16S rRNA V3-V4 sequencing (primers 338F/806R, MGI-2000). Enters as",
        "(exp(Veil) / 1.08)^-0.677 on CL/F (Eq. 1, Table 3), i.e. the effect",
        "is 1 at Veil = log(1.08) = 0.077%. Observed range 0.015-1.42% in the",
        "ESRD patients (Discussion); the healthy-volunteer range is not",
        "reported."
      ),
      source_name = "Veil"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pitavastatin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pitavastatin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pitavastatin", units = "ug", specimen = "plasma", verified = TRUE)
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
      "Single oral dose of pitavastatin 10 ug in a microdose cocktail with",
      "dabigatran etexilate 375 ug, rosuvastatin 50 ug, midazolam 10 ug and",
      "atorvastatin 100 ug, fasted, 24 h before dialysis."
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
    lka <- fixed(log(18)); label("Absorption rate constant (1/h)") # Table 3 PTV ka 18 (RSE NA; bootstrap 18 (18, 18)) -> fixed
    lcl <- log(16.5); label("Apparent clearance CL/F in healthy volunteers at Veil = 0.077% (L/h)") # Table 3 PTV CL/F 16.5 (RSE 4.50%)
    lvc <- log(25.2); label("Apparent central volume Vc/F (L)") # Table 3 PTV Vc/F 25.2 (RSE 7.00%)
    lvp <- log(76.1); label("Apparent peripheral volume Vp/F (L)") # Table 3 PTV Vp/F 76.1 (RSE 9.90%)
    lq <- log(8.01); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 PTV Q/F 8.01 (RSE printed 0; bootstrap 8.1 (6.6, 9.4))

    e_renalimp_esrd_cl <- -0.381; label("Fractional change in CL/F in ESRD (fraction)") # Table 3 'DISPTV on CL/F' -0.381 (RSE 10.8%); Eq. 1 (1 + DIS_PTV)
    e_gut_veillonella_pct_cl <- -0.677; label("Exponent of exp(Veil)/1.08 on CL/F (unitless)") # Table 3 'Veil on CL/F' -0.677 (RSE 21.7%); Eq. 1

    etalcl ~ log(1 + 0.096^2) # Table 3 PTV CL/F IIV 9.60% (RSE 37.6%); omega^2 = log(1 + CV^2)
    etalq ~ log(1 + 0.076^2) # Table 3 PTV Q/F IIV 7.60% (RSE not reported); omega^2 = log(1 + CV^2)

    propSd <- fixed(0); label("Proportional residual error (fraction); not reported in the source") # residual error model and magnitude not reported in Kong 2025
  })
  model({
    # Eq. 1: CL/F = TVCL x e^eta x (e^Veil / 1.08)^theta_Veil x (1 + DIS_PTV)
    veil_cl <- (exp(GUT_VEILLONELLA_PCT) / 1.08)^e_gut_veillonella_pct_cl
    esrd_cl <- 1 + e_renalimp_esrd_cl * RENALIMP_ESRD

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * veil_cl * esrd_cl
    vc <- exp(lvc)
    vp <- exp(lvp)
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
