Chen_2024_nAcetylglutamine_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous pharmacokinetic model for",
    "N-acetyl-L-glutamine (NAG), one of eight constituents of guhong injection",
    "(GHI) quantified in plasma, in male Sprague-Dawley rats subjected to 30",
    "min left-anterior-descending ligation followed by 1 h reperfusion",
    "(myocardial ischemia/reperfusion, MI/R) (Chen 2024). GHI was given as a",
    "single tail-vein injection of 2.5, 5 or 10 mL/kg; the N-acetyl-L-glutamine",
    "dose is the GHI volume dose times its content in GHI (30 mg/mL), i.e. 75",
    "mg/kg, 150 mg/kg, 300 mg/kg. Disposition was fitted separately in each",
    "dose group with Drug and Statistics (DAS) v3.2.6, so V1, V2, CL1 and Q are",
    "selected from the covariate DOSE_GHI_MLKG rather than through a",
    "dose-covariate function the authors did not fit. Direct-effect",
    "sigmoid-Emax models link the N-acetyl-L-glutamine plasma concentration to",
    "the GHI-minus-model-group difference in creatine kinase-MB (CK-MB),",
    "ischemia-modified albumin (IMA), cardiac troponin I (cTn I) and",
    "alpha-hydroxybutyrate dehydrogenase (alpha-HBDH) (E = Emax * C^gamma /",
    "(ED50^gamma + C^gamma); Tables 12, 13, 14, 15). Chen 2024 fitted a PK/PD",
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
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

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
    central     = list(analyte = "N-acetyl-L-glutamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "N-acetyl-L-glutamine", units = "mg", specimen = "plasma", verified = TRUE)
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
      "hydrochloride 0.95 mg/kg positive control). The PK parameters in Table 7",
      "come from the three GHI-treated groups with n = 6 rats per group (18 rats",
      "contribute to this model); the four PD biomarker series were measured with",
      "n = 3 per group. Plasma was sampled from the retro-orbital venous plexus",
      "at 2, 5, 10, 20, 40, 60, 90, 120, 240 and 360 min. No age or",
      "race/ethnicity data apply."
    )
  )

  ini({
    # Two-compartment disposition, fitted independently per GHI dose group
    # (Chen 2024 Table 7; V1 = central volume, V2 = peripheral volume,
    #  CL1 = elimination clearance, Q = intercompartmental clearance).
    # --- GHI-L (2.5 mL/kg GHI; N-acetyl-L-glutamine 75 mg/kg) ---
    lvc_ghil <- fixed(log(0.223)); label("Central volume V1, GHI-L (L/kg)")  # Table 7, V1, GHI-L
    lvp_ghil <- fixed(log(1.381)); label("Peripheral volume V2, GHI-L (L/kg)")  # Table 7, V2, GHI-L
    lcl_ghil <- fixed(log(1.577)); label("Elimination clearance CL1, GHI-L (L/h/kg)")  # Table 7, CL1, GHI-L
    lq_ghil <- fixed(log(8.185)); label("Intercompartmental clearance Q, GHI-L (L/h/kg)")  # Table 7, Q, GHI-L
    # --- GHI-M (5 mL/kg GHI; N-acetyl-L-glutamine 150 mg/kg) ---
    lvc_ghim <- fixed(log(0.024)); label("Central volume V1, GHI-M (L/kg)")  # Table 7, V1, GHI-M
    lvp_ghim <- fixed(log(0.445)); label("Peripheral volume V2, GHI-M (L/kg)")  # Table 7, V2, GHI-M
    lcl_ghim <- fixed(log(1.133)); label("Elimination clearance CL1, GHI-M (L/h/kg)")  # Table 7, CL1, GHI-M
    lq_ghim <- fixed(log(1.771)); label("Intercompartmental clearance Q, GHI-M (L/h/kg)")  # Table 7, Q, GHI-M
    # --- GHI-H (10 mL/kg GHI; N-acetyl-L-glutamine 300 mg/kg) ---
    lvc_ghih <- fixed(log(0.339)); label("Central volume V1, GHI-H (L/kg)")  # Table 7, V1, GHI-H
    lvp_ghih <- fixed(log(1.144)); label("Peripheral volume V2, GHI-H (L/kg)")  # Table 7, V2, GHI-H
    lcl_ghih <- fixed(log(1.633)); label("Elimination clearance CL1, GHI-H (L/h/kg)")  # Table 7, CL1, GHI-H
    lq_ghih <- fixed(log(8.631)); label("Intercompartmental clearance Q, GHI-H (L/h/kg)")  # Table 7, Q, GHI-H

    # Sigmoid-Emax effect on creatine kinase-MB (CK-MB) (Chen 2024 Table 12;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ug/mL).
    lemax_ckmb_ghil <- fixed(log(11.04)); label("Maximum effect Emax on dCK-MB, GHI-L (ng/mL)")  # Table 12, GHI-L, N-acetyl-L-glutamine
    lec50_ckmb_ghil <- fixed(log(23.588)); label("Half-maximal concentration ED50 on dCK-MB, GHI-L (ug/mL)")  # Table 12, GHI-L, N-acetyl-L-glutamine
    lhill_ckmb_ghil <- fixed(log(0.077)); label("Sigmoidicity gamma on dCK-MB, GHI-L (unitless)")  # Table 12, GHI-L, N-acetyl-L-glutamine
    lemax_ckmb_ghim <- fixed(log(13.089)); label("Maximum effect Emax on dCK-MB, GHI-M (ng/mL)")  # Table 12, GHI-M, N-acetyl-L-glutamine
    lec50_ckmb_ghim <- fixed(log(87.945)); label("Half-maximal concentration ED50 on dCK-MB, GHI-M (ug/mL)")  # Table 12, GHI-M, N-acetyl-L-glutamine
    lhill_ckmb_ghim <- fixed(log(0.096)); label("Sigmoidicity gamma on dCK-MB, GHI-M (unitless)")  # Table 12, GHI-M, N-acetyl-L-glutamine
    # GHI-H: no NAG/CKMB equation reported (PLSR coefficient not negative).

    # Sigmoid-Emax effect on ischemia-modified albumin (IMA) (Chen 2024 Table 13;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in U/mL, C in ug/mL).
    lemax_ima_ghil <- fixed(log(31.136)); label("Maximum effect Emax on dIMA, GHI-L (U/mL)")  # Table 13, GHI-L, N-acetyl-L-glutamine
    lec50_ima_ghil <- fixed(log(34.028)); label("Half-maximal concentration ED50 on dIMA, GHI-L (ug/mL)")  # Table 13, GHI-L, N-acetyl-L-glutamine
    lhill_ima_ghil <- fixed(log(0.055)); label("Sigmoidicity gamma on dIMA, GHI-L (unitless)")  # Table 13, GHI-L, N-acetyl-L-glutamine
    lemax_ima_ghim <- fixed(log(42.733)); label("Maximum effect Emax on dIMA, GHI-M (U/mL)")  # Table 13, GHI-M, N-acetyl-L-glutamine
    lec50_ima_ghim <- fixed(log(257.817)); label("Half-maximal concentration ED50 on dIMA, GHI-M (ug/mL)")  # Table 13, GHI-M, N-acetyl-L-glutamine
    lhill_ima_ghim <- fixed(log(0.01)); label("Sigmoidicity gamma on dIMA, GHI-M (unitless)")  # Table 13, GHI-M, N-acetyl-L-glutamine
    lemax_ima_ghih <- fixed(log(77.916)); label("Maximum effect Emax on dIMA, GHI-H (U/mL)")  # Table 13, GHI-H, N-acetyl-L-glutamine
    lec50_ima_ghih <- fixed(log(278.429)); label("Half-maximal concentration ED50 on dIMA, GHI-H (ug/mL)")  # Table 13, GHI-H, N-acetyl-L-glutamine
    lhill_ima_ghih <- fixed(log(0.398)); label("Sigmoidicity gamma on dIMA, GHI-H (unitless)")  # Table 13, GHI-H, N-acetyl-L-glutamine

    # Sigmoid-Emax effect on cardiac troponin I (cTn I) (Chen 2024 Table 14;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ug/mL).
    lemax_ctni_ghil <- fixed(log(4.608)); label("Maximum effect Emax on dcTn I, GHI-L (ng/mL)")  # Table 14, GHI-L, N-acetyl-L-glutamine
    lec50_ctni_ghil <- fixed(log(39.987)); label("Half-maximal concentration ED50 on dcTn I, GHI-L (ug/mL)")  # Table 14, GHI-L, N-acetyl-L-glutamine
    lhill_ctni_ghil <- fixed(log(0.037)); label("Sigmoidicity gamma on dcTn I, GHI-L (unitless)")  # Table 14, GHI-L, N-acetyl-L-glutamine
    lemax_ctni_ghim <- fixed(log(6.714)); label("Maximum effect Emax on dcTn I, GHI-M (ng/mL)")  # Table 14, GHI-M, N-acetyl-L-glutamine
    lec50_ctni_ghim <- fixed(log(79.994)); label("Half-maximal concentration ED50 on dcTn I, GHI-M (ug/mL)")  # Table 14, GHI-M, N-acetyl-L-glutamine
    lhill_ctni_ghim <- fixed(log(0.035)); label("Sigmoidicity gamma on dcTn I, GHI-M (unitless)")  # Table 14, GHI-M, N-acetyl-L-glutamine
    lemax_ctni_ghih <- fixed(log(10.158)); label("Maximum effect Emax on dcTn I, GHI-H (ng/mL)")  # Table 14, GHI-H, N-acetyl-L-glutamine
    lec50_ctni_ghih <- fixed(log(273.62)); label("Half-maximal concentration ED50 on dcTn I, GHI-H (ug/mL)")  # Table 14, GHI-H, N-acetyl-L-glutamine
    lhill_ctni_ghih <- fixed(log(0.228)); label("Sigmoidicity gamma on dcTn I, GHI-H (unitless)")  # Table 14, GHI-H, N-acetyl-L-glutamine

    # Sigmoid-Emax effect on alpha-hydroxybutyrate dehydrogenase (alpha-HBDH) (Chen 2024 Table 15;
    #  E = Emax * C^gamma / (ED50^gamma + C^gamma), E in ng/mL, C in ug/mL).
    lemax_hbdh_ghil <- fixed(log(48.85)); label("Maximum effect Emax on dalpha-HBDH, GHI-L (ng/mL)")  # Table 15, GHI-L, N-acetyl-L-glutamine
    lec50_hbdh_ghil <- fixed(log(39.96)); label("Half-maximal concentration ED50 on dalpha-HBDH, GHI-L (ug/mL)")  # Table 15, GHI-L, N-acetyl-L-glutamine
    lhill_hbdh_ghil <- fixed(log(0.023)); label("Sigmoidicity gamma on dalpha-HBDH, GHI-L (unitless)")  # Table 15, GHI-L, N-acetyl-L-glutamine
    lemax_hbdh_ghim <- fixed(log(72.312)); label("Maximum effect Emax on dalpha-HBDH, GHI-M (ng/mL)")  # Table 15, GHI-M, N-acetyl-L-glutamine
    lec50_hbdh_ghim <- fixed(log(211.779)); label("Half-maximal concentration ED50 on dalpha-HBDH, GHI-M (ug/mL)")  # Table 15, GHI-M, N-acetyl-L-glutamine
    lhill_hbdh_ghim <- fixed(log(0.01)); label("Sigmoidicity gamma on dalpha-HBDH, GHI-M (unitless)")  # Table 15, GHI-M, N-acetyl-L-glutamine
    lemax_hbdh_ghih <- fixed(log(91.925)); label("Maximum effect Emax on dalpha-HBDH, GHI-H (ng/mL)")  # Table 15, GHI-H, N-acetyl-L-glutamine
    lec50_hbdh_ghih <- fixed(log(154.598)); label("Half-maximal concentration ED50 on dalpha-HBDH, GHI-H (ug/mL)")  # Table 15, GHI-H, N-acetyl-L-glutamine
    lhill_hbdh_ghih <- fixed(log(0.136)); label("Sigmoidicity gamma on dalpha-HBDH, GHI-H (unitless)")  # Table 15, GHI-H, N-acetyl-L-glutamine

    addSd <- fixed(0); label("Additive residual SD on N-acetyl-L-glutamine plasma concentration (ug/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dCKMB <- fixed(0); label("Additive residual SD on dCK-MB (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dIMA <- fixed(0); label("Additive residual SD on dIMA (U/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dCTNI <- fixed(0); label("Additive residual SD on dcTn I (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
    addSd_dHBDH <- fixed(0); label("Additive residual SD on dalpha-HBDH (ng/mL; not reported)")  # Chen 2024 reports no residual-error model
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
    emax_ckmb_ghim <- exp(lemax_ckmb_ghim)
    ec50_ckmb_ghim <- exp(lec50_ckmb_ghim)
    hill_ckmb_ghim <- exp(lhill_ckmb_ghim)
    eff_ckmb_ghim <- emax_ckmb_ghim * Cc^hill_ckmb_ghim / (ec50_ckmb_ghim^hill_ckmb_ghim + Cc^hill_ckmb_ghim)
    # dCKMB is 0 in dose groups for which Chen 2024 fitted no NAG equation.
    dCKMB <- isGhiL * eff_ckmb_ghil +
      isGhiM * eff_ckmb_ghim

    emax_ima_ghil <- exp(lemax_ima_ghil)
    ec50_ima_ghil <- exp(lec50_ima_ghil)
    hill_ima_ghil <- exp(lhill_ima_ghil)
    eff_ima_ghil <- emax_ima_ghil * Cc^hill_ima_ghil / (ec50_ima_ghil^hill_ima_ghil + Cc^hill_ima_ghil)
    emax_ima_ghim <- exp(lemax_ima_ghim)
    ec50_ima_ghim <- exp(lec50_ima_ghim)
    hill_ima_ghim <- exp(lhill_ima_ghim)
    eff_ima_ghim <- emax_ima_ghim * Cc^hill_ima_ghim / (ec50_ima_ghim^hill_ima_ghim + Cc^hill_ima_ghim)
    emax_ima_ghih <- exp(lemax_ima_ghih)
    ec50_ima_ghih <- exp(lec50_ima_ghih)
    hill_ima_ghih <- exp(lhill_ima_ghih)
    eff_ima_ghih <- emax_ima_ghih * Cc^hill_ima_ghih / (ec50_ima_ghih^hill_ima_ghih + Cc^hill_ima_ghih)
    # dIMA is 0 in dose groups for which Chen 2024 fitted no NAG equation.
    dIMA <- isGhiL * eff_ima_ghil +
      isGhiM * eff_ima_ghim +
      isGhiH * eff_ima_ghih

    emax_ctni_ghil <- exp(lemax_ctni_ghil)
    ec50_ctni_ghil <- exp(lec50_ctni_ghil)
    hill_ctni_ghil <- exp(lhill_ctni_ghil)
    eff_ctni_ghil <- emax_ctni_ghil * Cc^hill_ctni_ghil / (ec50_ctni_ghil^hill_ctni_ghil + Cc^hill_ctni_ghil)
    emax_ctni_ghim <- exp(lemax_ctni_ghim)
    ec50_ctni_ghim <- exp(lec50_ctni_ghim)
    hill_ctni_ghim <- exp(lhill_ctni_ghim)
    eff_ctni_ghim <- emax_ctni_ghim * Cc^hill_ctni_ghim / (ec50_ctni_ghim^hill_ctni_ghim + Cc^hill_ctni_ghim)
    emax_ctni_ghih <- exp(lemax_ctni_ghih)
    ec50_ctni_ghih <- exp(lec50_ctni_ghih)
    hill_ctni_ghih <- exp(lhill_ctni_ghih)
    eff_ctni_ghih <- emax_ctni_ghih * Cc^hill_ctni_ghih / (ec50_ctni_ghih^hill_ctni_ghih + Cc^hill_ctni_ghih)
    # dCTNI is 0 in dose groups for which Chen 2024 fitted no NAG equation.
    dCTNI <- isGhiL * eff_ctni_ghil +
      isGhiM * eff_ctni_ghim +
      isGhiH * eff_ctni_ghih

    emax_hbdh_ghil <- exp(lemax_hbdh_ghil)
    ec50_hbdh_ghil <- exp(lec50_hbdh_ghil)
    hill_hbdh_ghil <- exp(lhill_hbdh_ghil)
    eff_hbdh_ghil <- emax_hbdh_ghil * Cc^hill_hbdh_ghil / (ec50_hbdh_ghil^hill_hbdh_ghil + Cc^hill_hbdh_ghil)
    emax_hbdh_ghim <- exp(lemax_hbdh_ghim)
    ec50_hbdh_ghim <- exp(lec50_hbdh_ghim)
    hill_hbdh_ghim <- exp(lhill_hbdh_ghim)
    eff_hbdh_ghim <- emax_hbdh_ghim * Cc^hill_hbdh_ghim / (ec50_hbdh_ghim^hill_hbdh_ghim + Cc^hill_hbdh_ghim)
    emax_hbdh_ghih <- exp(lemax_hbdh_ghih)
    ec50_hbdh_ghih <- exp(lec50_hbdh_ghih)
    hill_hbdh_ghih <- exp(lhill_hbdh_ghih)
    eff_hbdh_ghih <- emax_hbdh_ghih * Cc^hill_hbdh_ghih / (ec50_hbdh_ghih^hill_hbdh_ghih + Cc^hill_hbdh_ghih)
    # dHBDH is 0 in dose groups for which Chen 2024 fitted no NAG equation.
    dHBDH <- isGhiL * eff_hbdh_ghil +
      isGhiM * eff_hbdh_ghim +
      isGhiH * eff_hbdh_ghih

    Cc ~ add(addSd)
    dCKMB ~ add(addSd_dCKMB)
    dIMA ~ add(addSd_dIMA)
    dCTNI ~ add(addSd_dCTNI)
    dHBDH ~ add(addSd_dHBDH)
  })
}
