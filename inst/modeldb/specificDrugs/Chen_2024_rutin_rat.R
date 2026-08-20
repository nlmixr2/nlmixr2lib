Chen_2024_rutin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous pharmacokinetic model for",
    "rutin (RT), one of eight constituents of guhong injection (GHI) quantified",
    "in plasma, in male Sprague-Dawley rats subjected to 30 min",
    "left-anterior-descending ligation followed by 1 h reperfusion (myocardial",
    "ischemia/reperfusion, MI/R) (Chen 2024). GHI was given as a single",
    "tail-vein injection of 2.5, 5 or 10 mL/kg; the rutin dose is the GHI",
    "volume dose times its content in GHI (11.7 ug/mL), i.e. 29.25 ug/kg, 58.5",
    "ug/kg, 117 ug/kg. Disposition was fitted separately in each dose group",
    "with Drug and Statistics (DAS) v3.2.6, so V1, V2, CL1 and Q are selected",
    "from the covariate DOSE_GHI_MLKG rather than through a dose-covariate",
    "function the authors did not fit. Direct-effect sigmoid-Emax models link",
    "the rutin plasma concentration to the GHI-minus-model-group difference in",
    "creatine kinase-MB (CK-MB) and cardiac troponin I (cTn I) (E = Emax *",
    "C^gamma / (ED50^gamma + C^gamma); Tables 12, 14). Chen 2024 fitted a PK/PD",
    "model only for the analyte/biomarker/dose combinations whose PLSR",
    "coefficient was negative, so the effect of an unfitted combination is",
    "returned as zero rather than extrapolated. No between-subject variability",
    "or residual error was reported; every parameter is fixed at the published",
    "mean and the residual SDs are fixed at zero."
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
    central     = list(analyte = "rutin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rutin", units = "ug", specimen = "plasma", verified = TRUE)
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
      "hydrochloride 0.95 mg/kg positive control). The PK parameters in Table 9",
      "come from the three GHI-treated groups with n = 6 rats per group (18 rats",
      "contribute to this model); the four PD biomarker series were measured with",
      "n = 3 per group. Plasma was sampled from the retro-orbital venous plexus",
      "at 2, 5, 10, 20, 40, 60, 90, 120, 240 and 360 min. No age or",
      "race/ethnicity data apply."
    )
  )

  ini({
    # Two-compartment disposition, fitted independently per GHI dose group
    # (Chen 2024 Table 9; V1 = central volume, V2 = peripheral volume,
    #  CL1 = elimination clearance, Q = intercompartmental clearance).
    # --- GHI-L (2.5 mL/kg GHI; rutin 29.25 ug/kg) ---
    lvc_ghil <- fixed(log(0.379)); label("Central volume V1, GHI-L (L/kg)")  # Table 9, V1, GHI-L
    lvp_ghil <- fixed(log(1.671)); label("Peripheral volume V2, GHI-L (L/kg)")  # Table 9, V2, GHI-L
    lcl_ghil <- fixed(log(3.795)); label("Elimination clearance CL1, GHI-L (L/h/kg)")  # Table 9, CL1, GHI-L
    lq_ghil <- fixed(log(6.267)); label("Intercompartmental clearance Q, GHI-L (L/h/kg)")  # Table 9, Q, GHI-L
    # --- GHI-M (5 mL/kg GHI; rutin 58.5 ug/kg) ---
    lvc_ghim <- fixed(log(0.466)); label("Central volume V1, GHI-M (L/kg)")  # Table 9, V1, GHI-M
    lvp_ghim <- fixed(log(1.075)); label("Peripheral volume V2, GHI-M (L/kg)")  # Table 9, V2, GHI-M
    lcl_ghim <- fixed(log(2.872)); label("Elimination clearance CL1, GHI-M (L/h/kg)")  # Table 9, CL1, GHI-M
    lq_ghim <- fixed(log(5.018)); label("Intercompartmental clearance Q, GHI-M (L/h/kg)")  # Table 9, Q, GHI-M
    # --- GHI-H (10 mL/kg GHI; rutin 117 ug/kg) ---
    lvc_ghih <- fixed(log(0.239)); label("Central volume V1, GHI-H (L/kg)")  # Table 9, V1, GHI-H
    lvp_ghih <- fixed(log(0.855)); label("Peripheral volume V2, GHI-H (L/kg)")  # Table 9, V2, GHI-H
    lcl_ghih <- fixed(log(2.813)); label("Elimination clearance CL1, GHI-H (L/h/kg)")  # Table 9, CL1, GHI-H
    lq_ghih <- fixed(log(6.706)); label("Intercompartmental clearance Q, GHI-H (L/h/kg)")  # Table 9, Q, GHI-H

    # Sigmoid-Emax effect on creatine kinase-MB (CK-MB) (Chen 2024 Table 12;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ng/mL).
    # GHI-L: no RT/CKMB equation reported (PLSR coefficient not negative).
    lemax_ckmb_ghim <- fixed(log(11.991)); label("Maximum effect Emax on dCK-MB, GHI-M (ng/mL)")  # Table 12, GHI-M, rutin
    lec50_ckmb_ghim <- fixed(log(4.448)); label("Half-maximal concentration ED50 on dCK-MB, GHI-M (ng/mL)")  # Table 12, GHI-M, rutin
    lhill_ckmb_ghim <- fixed(log(1.357)); label("Sigmoidicity gamma on dCK-MB, GHI-M (unitless)")  # Table 12, GHI-M, rutin
    lemax_ckmb_ghih <- fixed(log(15.908)); label("Maximum effect Emax on dCK-MB, GHI-H (ng/mL)")  # Table 12, GHI-H, rutin
    lec50_ckmb_ghih <- fixed(log(17.239)); label("Half-maximal concentration ED50 on dCK-MB, GHI-H (ng/mL)")  # Table 12, GHI-H, rutin
    lhill_ckmb_ghih <- fixed(log(1.459)); label("Sigmoidicity gamma on dCK-MB, GHI-H (unitless)")  # Table 12, GHI-H, rutin

    # Sigmoid-Emax effect on cardiac troponin I (cTn I) (Chen 2024 Table 14;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ng/mL).
    lemax_ctni_ghil <- fixed(log(4.632)); label("Maximum effect Emax on dcTn I, GHI-L (ng/mL)")  # Table 14, GHI-L, rutin
    lec50_ctni_ghil <- fixed(log(2.595)); label("Half-maximal concentration ED50 on dcTn I, GHI-L (ng/mL)")  # Table 14, GHI-L, rutin
    lhill_ctni_ghil <- fixed(log(0.01)); label("Sigmoidicity gamma on dcTn I, GHI-L (unitless)")  # Table 14, GHI-L, rutin
    # GHI-M: no RT/CTNI equation reported (PLSR coefficient not negative).
    # GHI-H: no RT/CTNI equation reported (PLSR coefficient not negative).

    addSd <- fixed(0); label("Additive residual SD on rutin plasma concentration (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dCKMB <- fixed(0); label("Additive residual SD on dCK-MB (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dCTNI <- fixed(0); label("Additive residual SD on dcTn I (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
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

    emax_ckmb_ghim <- exp(lemax_ckmb_ghim)
    ec50_ckmb_ghim <- exp(lec50_ckmb_ghim)
    hill_ckmb_ghim <- exp(lhill_ckmb_ghim)
    eff_ckmb_ghim <- emax_ckmb_ghim * Cc^hill_ckmb_ghim / (ec50_ckmb_ghim^hill_ckmb_ghim + Cc^hill_ckmb_ghim)
    emax_ckmb_ghih <- exp(lemax_ckmb_ghih)
    ec50_ckmb_ghih <- exp(lec50_ckmb_ghih)
    hill_ckmb_ghih <- exp(lhill_ckmb_ghih)
    eff_ckmb_ghih <- emax_ckmb_ghih * Cc^hill_ckmb_ghih / (ec50_ckmb_ghih^hill_ckmb_ghih + Cc^hill_ckmb_ghih)
    # dCKMB is 0 in dose groups for which Chen 2024 fitted no RT equation.
    dCKMB <- isGhiM * eff_ckmb_ghim +
      isGhiH * eff_ckmb_ghih

    emax_ctni_ghil <- exp(lemax_ctni_ghil)
    ec50_ctni_ghil <- exp(lec50_ctni_ghil)
    hill_ctni_ghil <- exp(lhill_ctni_ghil)
    eff_ctni_ghil <- emax_ctni_ghil * Cc^hill_ctni_ghil / (ec50_ctni_ghil^hill_ctni_ghil + Cc^hill_ctni_ghil)
    # dCTNI is 0 in dose groups for which Chen 2024 fitted no RT equation.
    dCTNI <- isGhiL * eff_ctni_ghil

    Cc ~ add(addSd)
    dCKMB ~ add(addSd_dCKMB)
    dCTNI ~ add(addSd_dCTNI)
  })
}
