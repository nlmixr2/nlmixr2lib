Chen_2024_kaempferol3ORutinoside_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous pharmacokinetic model for",
    "kaempferol-3-O-rutinoside (K-3-R), one of eight constituents of guhong",
    "injection (GHI) quantified in plasma, in male Sprague-Dawley rats",
    "subjected to 30 min left-anterior-descending ligation followed by 1 h",
    "reperfusion (myocardial ischemia/reperfusion, MI/R) (Chen 2024). GHI was",
    "given as a single tail-vein injection of 2.5, 5 or 10 mL/kg; the",
    "kaempferol-3-O-rutinoside dose is the GHI volume dose times its content in",
    "GHI (65.7 ug/mL), i.e. 164.25 ug/kg, 328.5 ug/kg, 657 ug/kg. Disposition",
    "was fitted separately in each dose group with Drug and Statistics (DAS)",
    "v3.2.6, so V1, V2, CL1 and Q are selected from the covariate DOSE_GHI_MLKG",
    "rather than through a dose-covariate function the authors did not fit.",
    "Direct-effect sigmoid-Emax models link the kaempferol-3-O-rutinoside",
    "plasma concentration to the GHI-minus-model-group difference in creatine",
    "kinase-MB (CK-MB) (E = Emax * C^gamma / (ED50^gamma + C^gamma); Tables",
    "12). Chen 2024 fitted a PK/PD model only for the analyte/biomarker/dose",
    "combinations whose PLSR coefficient was negative, so the effect of an",
    "unfitted combination is returned as zero rather than extrapolated. No",
    "between-subject variability or residual error was reported; every",
    "parameter is fixed at the published mean and the residual SDs are fixed at",
    "zero."
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
    central     = list(analyte = "kaempferol-3-O-rutinoside", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "kaempferol-3-O-rutinoside", units = "ug", specimen = "plasma", verified = TRUE)
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
      "hydrochloride 0.95 mg/kg positive control). The PK parameters in Table 6",
      "come from the three GHI-treated groups with n = 6 rats per group (18 rats",
      "contribute to this model); the four PD biomarker series were measured with",
      "n = 3 per group. Plasma was sampled from the retro-orbital venous plexus",
      "at 2, 5, 10, 20, 40, 60, 90, 120, 240 and 360 min. No age or",
      "race/ethnicity data apply."
    )
  )

  ini({
    # Two-compartment disposition, fitted independently per GHI dose group
    # (Chen 2024 Table 6; V1 = central volume, V2 = peripheral volume,
    #  CL1 = elimination clearance, Q = intercompartmental clearance).
    # --- GHI-L (2.5 mL/kg GHI; kaempferol-3-O-rutinoside 164.25 ug/kg) ---
    lvc_ghil <- fixed(log(0.546)); label("Central volume V1, GHI-L (L/kg)")  # Table 6, V1, GHI-L
    lvp_ghil <- fixed(log(2.437)); label("Peripheral volume V2, GHI-L (L/kg)")  # Table 6, V2, GHI-L
    lcl_ghil <- fixed(log(5.158)); label("Elimination clearance CL1, GHI-L (L/h/kg)")  # Table 6, CL1, GHI-L
    lq_ghil <- fixed(log(9.311)); label("Intercompartmental clearance Q, GHI-L (L/h/kg)")  # Table 6, Q, GHI-L
    # --- GHI-M (5 mL/kg GHI; kaempferol-3-O-rutinoside 328.5 ug/kg) ---
    lvc_ghim <- fixed(log(0.603)); label("Central volume V1, GHI-M (L/kg)")  # Table 6, V1, GHI-M
    lvp_ghim <- fixed(log(1.669)); label("Peripheral volume V2, GHI-M (L/kg)")  # Table 6, V2, GHI-M
    lcl_ghim <- fixed(log(4.134)); label("Elimination clearance CL1, GHI-M (L/h/kg)")  # Table 6, CL1, GHI-M
    lq_ghim <- fixed(log(4.919)); label("Intercompartmental clearance Q, GHI-M (L/h/kg)")  # Table 6, Q, GHI-M
    # --- GHI-H (10 mL/kg GHI; kaempferol-3-O-rutinoside 657 ug/kg) ---
    lvc_ghih <- fixed(log(0.49)); label("Central volume V1, GHI-H (L/kg)")  # Table 6, V1, GHI-H
    lvp_ghih <- fixed(log(1.244)); label("Peripheral volume V2, GHI-H (L/kg)")  # Table 6, V2, GHI-H
    lcl_ghih <- fixed(log(4.481)); label("Elimination clearance CL1, GHI-H (L/h/kg)")  # Table 6, CL1, GHI-H
    lq_ghih <- fixed(log(8.485)); label("Intercompartmental clearance Q, GHI-H (L/h/kg)")  # Table 6, Q, GHI-H

    # Sigmoid-Emax effect on creatine kinase-MB (CK-MB) (Chen 2024 Table 12;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ng/mL).
    lemax_ckmb_ghil <- fixed(log(11.04)); label("Maximum effect Emax on dCK-MB, GHI-L (ng/mL)")  # Table 12, GHI-L, kaempferol-3-O-rutinoside
    lec50_ckmb_ghil <- fixed(log(14.175)); label("Half-maximal concentration ED50 on dCK-MB, GHI-L (ng/mL)")  # Table 12, GHI-L, kaempferol-3-O-rutinoside
    lhill_ckmb_ghil <- fixed(log(0.081)); label("Sigmoidicity gamma on dCK-MB, GHI-L (unitless)")  # Table 12, GHI-L, kaempferol-3-O-rutinoside
    # GHI-M: no K-3-R/CKMB equation reported (PLSR coefficient not negative).
    # GHI-H: no K-3-R/CKMB equation reported (PLSR coefficient not negative).

    addSd <- fixed(0); label("Additive residual SD on kaempferol-3-O-rutinoside plasma concentration (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dCKMB <- fixed(0); label("Additive residual SD on dCK-MB (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
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

    emax_ckmb_ghil <- exp(lemax_ckmb_ghil)
    ec50_ckmb_ghil <- exp(lec50_ckmb_ghil)
    hill_ckmb_ghil <- exp(lhill_ckmb_ghil)
    eff_ckmb_ghil <- emax_ckmb_ghil * Cc^hill_ckmb_ghil / (ec50_ckmb_ghil^hill_ckmb_ghil + Cc^hill_ckmb_ghil)
    # dCKMB is 0 in dose groups for which Chen 2024 fitted no K-3-R equation.
    dCKMB <- isGhiL * eff_ckmb_ghil

    Cc ~ add(addSd)
    dCKMB ~ add(addSd_dCKMB)
  })
}
