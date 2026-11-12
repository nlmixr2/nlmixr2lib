Kong_2025_rosuvastatin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral microdose rosuvastatin (OATP1B/BCRP probe) in Chinese healthy",
    "volunteers and haemodialysis-dependent end-stage renal disease patients",
    "(Kong 2025). CL/F falls with the stool relative abundance of the gut",
    "taxon Clostridium cluster XIVb and Vc/F falls with serum albumin; ESRD",
    "itself is not a covariate."
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
    GUT_CLOSTRIDIUM_XIVB_PCT = list(
      description = "Stool relative abundance of the gut taxon Clostridium cluster XIVb",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "16S rRNA V3-V4 sequencing (primers 338F/806R, MGI-2000); the RDP",
        "genus-level taxon Clostridium_XlVb. Enters as (exp(Clos) /",
        "1.08)^-3.26 on CL/F (Eq. 4, Table 3), i.e. the effect is 1 at Clos",
        "= log(1.08) = 0.077%. The observed range is not reported."
      ),
      source_name = "Clos"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 3 'ALB on Vc/F' -6.42 without a printed functional form; the",
        "maintainers encode it as the power form (ALB / 44)^-6.42, with 44",
        "g/L the rounded pooled mean of Table 1 (44.9 g/L in 14 HV, 42.1 g/L",
        "in 10 ESRD). Observed means span roughly 40-48 g/L, over which the",
        "factor runs from about 1.8 to 0.6."
      ),
      source_name = "ALB"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rosuvastatin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rosuvastatin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rosuvastatin", units = "ug", specimen = "plasma", verified = TRUE)
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
      "Single oral dose of rosuvastatin 50 ug in a microdose cocktail with",
      "dabigatran etexilate 375 ug, pitavastatin 10 ug, midazolam 10 ug and",
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
    lka <- log(0.30); label("Absorption rate constant (1/h)") # Table 3 RSV ka 0.30 (RSE 17.4%)
    lcl <- log(43.1); label("Apparent clearance CL/F at Clos = 0.077% (L/h)") # Table 3 RSV CL/F 43.1 (RSE 8.60%)
    lvc <- log(289); label("Apparent central volume Vc/F at ALB 44 g/L (L)") # Table 3 RSV Vc/F 289 (RSE 17.2%)
    lvp <- log(2740); label("Apparent peripheral volume Vp/F (L)") # Table 3 RSV Vp/F 2,740 (RSE 22.5%)
    lq <- log(48.6); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 RSV Q/F 48.6 (RSE 11.7%)

    e_gut_clostridium_xivb_pct_cl <- -3.26; label("Exponent of exp(Clos)/1.08 on CL/F (unitless)") # Table 3 'Clos on CL/F' -3.26 (RSE 54.3%); Eq. 4
    e_alb_vc <- -6.42; label("Power exponent of ALB/44 on Vc/F (unitless)") # Table 3 'ALB on Vc/F' -6.42 (RSE 24.3%); functional form not printed

    etalcl ~ log(1 + 0.508^2) # Table 3 RSV CL/F IIV 50.8% (RSE 30.4%); omega^2 = log(1 + CV^2)
    etalvc ~ log(1 + 0.224^2) # Table 3 RSV Vc/F IIV 22.4% (RSE 59.4%); omega^2 = log(1 + CV^2)

    propSd <- fixed(0); label("Proportional residual error (fraction); not reported in the source") # residual error model and magnitude not reported in Kong 2025
  })
  model({
    # Eq. 4: CL/F = TVCL x e^eta x (e^Clos / 1.08)^theta_Clos
    clos_cl <- (exp(GUT_CLOSTRIDIUM_XIVB_PCT) / 1.08)^e_gut_clostridium_xivb_pct_cl
    # ALB on Vc/F: power form centred on the rounded pooled mean 44 g/L
    alb_vc <- (ALB / 44)^e_alb_vc

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * clos_cl
    vc <- exp(lvc + etalvc) * alb_vc
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
