Chen_2024_hyperoside_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous pharmacokinetic model for",
    "hyperoside (HYP), one of eight constituents of guhong injection (GHI)",
    "quantified in plasma, in male Sprague-Dawley rats subjected to 30 min",
    "left-anterior-descending ligation followed by 1 h reperfusion (myocardial",
    "ischemia/reperfusion, MI/R) (Chen 2024). GHI was given as a single",
    "tail-vein injection of 2.5, 5 or 10 mL/kg; the hyperoside dose is the GHI",
    "volume dose times its content in GHI (3.4 ug/mL), i.e. 8.5 ug/kg, 17",
    "ug/kg, 34 ug/kg. Disposition was fitted separately in each dose group with",
    "Drug and Statistics (DAS) v3.2.6, so V1, V2, CL1 and Q are selected from",
    "the covariate DOSE_GHI_MLKG rather than through a dose-covariate function",
    "the authors did not fit. Chen 2024 fitted no concentration-effect model",
    "for hyperoside: its weighted PLSR regression coefficients against all four",
    "myocardial biomarkers were positive (Table 10), so the paper's PK/PD",
    "tables (Tables 12-15) contain no entry for it and this file carries the PK",
    "layer only. No between-subject variability or residual error was reported;",
    "every parameter is fixed at the published mean and the residual SDs are",
    "fixed at zero."
  )
  reference <- paste(
    "Chen HY, Li C, Shao CY, Wu YJ, Wan HT, He Y (2024). An auxiliary strategy",
    "of partial least squares regression in pharmacokinetic/pharmacodynamic",
    "studies: A case of application of guhong injection in myocardial",
    "ischemia/reperfusion rats. J Food Drug Anal 32(1):79-102.",
    "doi:10.38212/2224-6614.3492"
  )
  vignette <- "Chen_2024_guhongInjection"
  units <- list(time = "h", dosing = "ug/kg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_GHI_MLKG = list(
      description        = "Administered guhong injection volume dose; selects the dose-group-specific disposition and sigmoid-Emax parameter set",
      units              = "mL/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per rat. Chen 2024 Section 2.3 gives three GHI dose groups -",
        "GHI-L 2.5 mL/kg, GHI-M 5 mL/kg and GHI-H 10 mL/kg - each fitted as an",
        "independent two-compartment model, so this column is a group selector",
        "rather than a continuous covariate with an estimated exponent. The model",
        "assigns GHI-L below 3.75 mL/kg, GHI-M from 3.75 to below 7.5 mL/kg and",
        "GHI-H at or above 7.5 mL/kg (midpoints of the studied levels); values",
        "outside 2.5-10 mL/kg are extrapolations of the nearest fitted group."
      ),
      source_name        = "GHI dose (mL/kg)"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "hyperoside", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hyperoside", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "rat (Sprague-Dawley)",
    n_subjects    = 18,
    n_studies     = 1,
    sex           = "male",
    weight_range  = "320-340 g",
    disease_state = "myocardial ischemia/reperfusion injury (30 min LAD ligation, 1 h reperfusion)",
    dose_range    = "guhong injection 2.5, 5 or 10 mL/kg as a single tail-vein injection",
    regions       = "China",
    notes         = paste(
      "Chen 2024 Sections 2.2-2.4: adult male Sprague-Dawley rats, 320-340 g,",
      "randomised to six groups (sham, model, GHI-L, GHI-M, GHI-H, verapamil",
      "hydrochloride 0.95 mg/kg positive control). The PK parameters in Table 4",
      "come from the three GHI-treated groups with n = 6 rats per group (18 rats",
      "contribute to this model); the four PD biomarker series were measured with",
      "n = 3 per group. Plasma was sampled from the retro-orbital venous plexus",
      "at 2, 5, 10, 20, 40, 60, 90, 120, 240 and 360 min. No age or",
      "race/ethnicity data apply."
    )
  )

  ini({
    # Two-compartment disposition, fitted independently per GHI dose group
    # (Chen 2024 Table 4; V1 = central volume, V2 = peripheral volume,
    #  CL1 = elimination clearance, Q = intercompartmental clearance).
    # --- GHI-L (2.5 mL/kg GHI; hyperoside 8.5 ug/kg) ---
    lvc_ghil <- fixed(log(0.493)); label("Central volume V1, GHI-L (L/kg)")  # Table 4, V1, GHI-L
    lvp_ghil <- fixed(log(2.979)); label("Peripheral volume V2, GHI-L (L/kg)")  # Table 4, V2, GHI-L
    lcl_ghil <- fixed(log(5.336)); label("Elimination clearance CL1, GHI-L (L/h/kg)")  # Table 4, CL1, GHI-L
    lq_ghil <- fixed(log(8.69)); label("Intercompartmental clearance Q, GHI-L (L/h/kg)")  # Table 4, Q, GHI-L
    # --- GHI-M (5 mL/kg GHI; hyperoside 17 ug/kg) ---
    lvc_ghim <- fixed(log(0.34)); label("Central volume V1, GHI-M (L/kg)")  # Table 4, V1, GHI-M
    lvp_ghim <- fixed(log(0.935)); label("Peripheral volume V2, GHI-M (L/kg)")  # Table 4, V2, GHI-M
    lcl_ghim <- fixed(log(4.595)); label("Elimination clearance CL1, GHI-M (L/h/kg)")  # Table 4, CL1, GHI-M
    lq_ghim <- fixed(log(4.572)); label("Intercompartmental clearance Q, GHI-M (L/h/kg)")  # Table 4, Q, GHI-M
    # --- GHI-H (10 mL/kg GHI; hyperoside 34 ug/kg) ---
    lvc_ghih <- fixed(log(0.309)); label("Central volume V1, GHI-H (L/kg)")  # Table 4, V1, GHI-H
    lvp_ghih <- fixed(log(1.207)); label("Peripheral volume V2, GHI-H (L/kg)")  # Table 4, V2, GHI-H
    lcl_ghih <- fixed(log(4.915)); label("Elimination clearance CL1, GHI-H (L/h/kg)")  # Table 4, CL1, GHI-H
    lq_ghih <- fixed(log(7.26)); label("Intercompartmental clearance Q, GHI-H (L/h/kg)")  # Table 4, Q, GHI-H

    addSd <- fixed(0); label("Additive residual SD on hyperoside plasma concentration (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
  })

  model({
    # Dose-group selectors; exactly one is 1 for any DOSE_GHI_MLKG.
    isGhiL <- (DOSE_GHI_MLKG < 3.75)
    isGhiM <- (DOSE_GHI_MLKG >= 3.75) * (DOSE_GHI_MLKG < 7.5)
    isGhiH <- (DOSE_GHI_MLKG >= 7.5)

    # Disposition parameters for the rat's dose group.
    vc <- isGhiL * exp(lvc_ghil) + isGhiM * exp(lvc_ghim) + isGhiH * exp(lvc_ghih)
    vp <- isGhiL * exp(lvp_ghil) + isGhiM * exp(lvp_ghim) + isGhiH * exp(lvp_ghih)
    cl <- isGhiL * exp(lcl_ghil) + isGhiM * exp(lcl_ghim) + isGhiH * exp(lcl_ghih)
    q  <- isGhiL * exp(lq_ghil)  + isGhiM * exp(lq_ghim)  + isGhiH * exp(lq_ghih)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
