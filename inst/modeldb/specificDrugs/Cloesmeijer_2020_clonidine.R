Cloesmeijer_2020_clonidine <- function() {
  description <- paste0(
    "Two-compartment population PK model for intravenous clonidine in ",
    "critically ill, intubated and sedated adult intensive care unit ",
    "patients receiving continuous IV infusion (600-1800 ug/day, with or ",
    "without a 4-h loading infusion). Central volume is allometrically ",
    "scaled to body weight with a fixed exponent of 1 (reference 70 kg); ",
    "clearance increases linearly with time after the start of the ",
    "clonidine infusion (0.213 percent per hour). IIV on CL and V1 only; ",
    "combined additive and proportional residual error."
  )
  reference <- paste(
    "Cloesmeijer ME, van den Oever HLA, Mathot RAA, Zeeman M,",
    "Kruisdijk-Gerritsen A, Bles CMA, Nassikovker P, de Meijer AR,",
    "van Steveninck FL, Arbouw MEL.",
    "Optimising the dose of clonidine to achieve sedation in intensive care",
    "unit patients with population pharmacokinetics.",
    "Br J Clin Pharmacol. 2020;86(8):1620-1631.",
    "doi:10.1111/bcp.14273.",
    sep = " "
  )
  vignette <- "Cloesmeijer_2020_clonidine"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    central = list(analyte = "clonidine", units = "ug", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "clonidine", units = "ug", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of V1 only, (WT/70)^1 with the exponent fixed at 1 (Methods Eq 4; Results 3.1.3 'Allometric scaling based on body weight was applied to V1 with an exponent of 1'). CL, Q and V2 are not weight-scaled in the final model (Table 2 model 9; Table 3 units L/h and L). Time-fixed per subject.",
      source_name = "Bodyweight"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on V1 in forward addition (Table 2 model 3) but removed in backward elimination; not in the final model. Methods 2.3.2: Cockcroft-Gault, CKD-EPI or 24-h urine creatinine clearance, adjusted for the patient's BSA.",
      source_name = "CLcr"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on V1 in forward addition (Table 2 model 4) but removed in backward elimination; not in the final model.",
      source_name = "albumin"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on V1 in forward addition (Table 2 model 5) but removed in backward elimination; not in the final model.",
      source_name = "bilirubin"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous veno-venous haemofiltration (1 = on CVVH, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Tested on CL and V1 (Results 3.1.2; 3 of 24 patients) and not significant; not in the final model.",
      source_name = "CVVH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "25-83 years",
    age_median = "67 years",
    weight_range = "53-113 kg",
    weight_median = "84 kg",
    sex_female_pct = 33.3,
    race_ethnicity = NULL,
    disease_state = "Intubated and sedated critically ill adults in the ICU with an expected stay of at least 3 days; 3 of 24 received continuous veno-venous haemofiltration (CVVH). Clonidine was added to standard sedation with morphine plus midazolam or propofol.",
    dose_range = "Continuous IV infusion of 600, 1200 or 1800 ug/day (25, 50 or 75 ug/h; 8 patients per group); 4 patients per group also received a loading dose of 50 percent of the daily dose over 4 h.",
    regions = "Single centre, Deventer Hospital ICU, The Netherlands (NCT02466373).",
    notes = "Table 1: height 173 (155-189) cm, BMI 27 (20-44) kg/m^2, BSA 2.0 (1.5-2.4) m^2, serum creatinine 74 (32-441) umol/L, albumin 20 (12-32) g/L, bilirubin 7 (<3-59) umol/L; treatment duration 96 (25-171) h. 275 plasma concentrations (5-16 per patient); 11 samples below the 0.1 ug/L LLOQ were discarded in the final model. Screened covariates not retained: BSA, BMI, height, age, sex, CLcr, albumin, bilirubin, CVVH."
  )

  ini({
    # Structural parameters -- Table 3 'Final parameter values'. Dose in ug,
    # volumes in L, clearances in L/h, so ug / L = ug/L.
    lcl <- log(17.0); label("Clearance at the start of the clonidine infusion (L/h)") # Table 3 'CL (L/h) 17.0 (10)'
    lvc <- log(124); label("Central volume of distribution for a 70 kg patient (L)") # Table 3 'V1 (L/70 kg) 124 (36)'
    lq <- log(83.7); label("Intercompartmental clearance (L/h)") # Table 3 'Q(L/h) 83.7 (35)'
    lvp <- log(178); label("Peripheral volume of distribution (L)") # Table 3 'V2 (L) 178 (35)'

    # Allometric exponent on V1, fixed at 1 (Methods Eq 4 and following text;
    # Results 3.1.3).
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V1 (unitless)") # Methods 2.3.2 'the power exponent was fixed at 1 for central volume of distribution (V1)'

    # Linear increase of CL with time after the start of the clonidine
    # infusion. Table 3 prints 'Increase CL per hour 0.213' and the Discussion
    # states 'CL increased linearly with 0.213%/h from baseline ... 17 L/h at
    # the start of the treatment and increased to 20.4 L/h after 4 days':
    # 17 * (1 + 0.00213 * 96) = 20.5 L/h, so the slope is 0.00213 per hour.
    cl_time_slope <- 0.00213; label("Fractional linear increase in CL per hour after the start of infusion (1/h)") # Table 3 'Increase CL per hour 0.213 (19)' as percent per hour (Discussion)

    # IIV -- Table 3 reports %CV; omega^2 = log(CV^2 + 1).
    #   log(0.333^2 + 1) = 0.105189
    #   log(0.668^2 + 1) = 0.369689
    etalcl ~ 0.105189 # Table 3 IIV CL 33.3 percent CV
    etalvc ~ 0.369689 # Table 3 IIV V1 66.8 percent CV

    # Residual error -- combined additive + proportional (Results 3.1.1;
    # Methods Eqs 2-3). Table 3 values taken as standard deviations.
    propSd <- 0.141; label("Proportional residual error (fraction)") # Table 3 'Proportional error 0.141 (4)'
    addSd <- 0.0532; label("Additive residual error (ug/L)") # Table 3 'Additive error (ug/L) 0.0532 (14)'
  })

  model({
    # Individual parameters. t is time in hours since the start of the
    # clonidine infusion (the first dose must be at t = 0).
    cl <- exp(lcl + etalcl) * (1 + cl_time_slope * t)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
