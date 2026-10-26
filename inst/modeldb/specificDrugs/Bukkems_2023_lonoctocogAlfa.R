Bukkems_2023_lonoctocogAlfa <- function() {
  description <- "Two-compartment population PK model for factor VIII (FVIII) activity (IU/dL) after intravenous rFVIII-SingleChain (lonoctocog alfa, Afstyla) prophylaxis in children, adolescents and adults with severe hemophilia A, fitted simultaneously to chromogenic (CSA) and one-stage (OSA) assay FVIII activity (Bukkems 2023). Central volume and clearance scale with body weight (power exponents 0.834 and 0.702, centred at 68 kg) and clearance falls with von Willebrand factor antigen (exponent -0.652, centred at 113 IU/dL). An estimated endogenous FVIII level (0.265 IU/dL, IIV 166% CV, correlated with CL) is added to the prediction, and a residual FVIII level from an unknown pre-study dose (1.41 IU/dL, decaying at the terminal rate) is added for subjects whose pre-first-dose FVIII exceeded 1 IU/dL. OSA activity reads 0.526 times the CSA activity of the drug. Residual error is combined proportional and additive for CSA and power plus proportional for OSA, all scaled by an individual random effect (IIV 36.8% CV). Sister repeated time-to-event bleeding model from the same paper: modellib('Bukkems_2023_lonoctocogAlfa_bleed')."
  reference <- paste(
    "Bukkems LH, Jonsson S, Cnossen MH, Karlsson MO, Mathot RAA; OPTI-CLOT",
    "studies and SYMPHONY consortium. Relationship between factor VIII levels",
    "and bleeding for rFVIII-SingleChain in severe hemophilia A: A repeated",
    "time-to-event analysis. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(5):706-718. doi:10.1002/psp4.12938. PMID:36965157.",
    "PMCID:PMC10196421. Open Access under CC BY-NC 4.0.",
    "Final estimates are in Table 2; the NONMEM control stream of the final",
    "PK model is Appendix S2 of the supplement, and its $THETA / $OMEGA /",
    "$SIGMA values (the unrounded final estimates) are used here.",
    "The starting model was Zhang Y et al. J Thromb Haemost.",
    "2017;15(6):1106-1114 (doi:10.1111/jth.13662), re-estimated on the larger",
    "pooled dataset.",
    sep = " "
  )
  vignette <- "Bukkems_2023_lonoctocogAlfa"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  # Bookkeeping state for the residual FVIII from an unknown pre-study dose
  # (Appendix S2: IPRED = BASE + CSA + EXTR * TERM, TERM = EXP(-BETA * TIME)
  # when the pre-first-dose FVIII exceeded 1 IU/dL). A unit dose into
  # `prior_dose` at the first record, scaled by f(prior_dose) = EXTR, gives
  # exactly EXTR * exp(-beta * t). Same construction as the `virtual` state of
  # deJager_2020_factorviii.R.
  paper_specific_compartments <- c("prior_dose")
  # One random effect multiplies every residual-error component
  # (Appendix S2: Y = IPRED + (...) * EXP(ETA(6))); it has no ini() typical value.
  paper_specific_etas <- "etaresSd"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling of V1 (exponent 0.834) and CL (exponent 0.702), both centred at 68 kg (Appendix S2 $PK 'TVV1 = THETA(1) * (WT/68)**THETA(5)'; Table 2 footnote equations '(WT/68)'). The Table 2 footnote text 'typical values are obtained for a patient with a bodyweight of 70 kg' disagrees with its own equations; 68 kg, the value in the equations and the code, is used. Cohort median 60 kg (range 10-120; Table 1).",
      source_name = "WT"
    ),
    VWF = list(
      description = "Von Willebrand factor antigen (VWF:Ag)",
      units = "IU/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, exponent -0.652, centred at 113 IU/dL (Table 2 footnote 'CL = 2.08 * (WT/68)**0.70 * (VWF/113)**-0.65'). VWF antigen, one value per subject. VWF was missing for 154 of 241 subjects (Table 1 footnote); for those the source imputed it inside the model from age, VWF = 79.367 + 1.6672 * AGE + eta with eta variance 1931 (IU/dL)^2 fixed (Appendix S2 $PK; Figure S1 prints 79.58 + 1.66 * age). That imputation is a missing-data device and is not encoded here: supply VWF for every simulated subject (the vignette draws it from the same age regression). Cohort median 113 IU/dL (range 43-296; Table 1).",
      source_name = "VWF"
    ),
    ASSAY_OSA = list(
      description = "FVIII activity assay: 1 = one-stage clotting assay (OSA), 0 = chromogenic substrate assay (CSA)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CSA; the structural parameters are on the CSA scale)",
      notes = "Per-observation. The drug-derived FVIII reads 0.526 times lower by OSA than by CSA (Table 2 'OSA instead of CSA' 0.53; Appendix S2 'OSA = CSA * THETA(9)'); the endogenous and residual-previous-dose terms are not scaled. The residual-error model also switches with the assay. Source column ASSAY coded 1 = OSA, 2 = CSA, so ASSAY_OSA = 2 - ASSAY.",
      source_name = "ASSAY"
    )
  )

  compartmentData <- list(
    central = list(analyte = "lonoctocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "lonoctocog alfa (FVIII activity)",
      units = "IU",
      specimen = "plasma",
      verified = TRUE
    ),
    prior_dose = list(
      analyte = "FVIII activity remaining from an unknown pre-study dose",
      units = "IU/dL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 241L,
    n_studies = 3L,
    age_range = "1-58 years",
    age_median = "19 years",
    weight_range = "10-120 kg",
    weight_median = "60 kg",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Severe hemophilia A (endogenous FVIII < 1 IU/dL) on rFVIII-SingleChain prophylaxis; persons treated on demand were excluded.",
    dose_range = "Intravenous rFVIII-SingleChain prophylaxis, 15-50 IU/kg every second day or two to three times per week (median weekly dose 91.6 IU/kg, range 19.6-246).",
    regions = "Europe 56%, Japan 4%, United States 9%, rest of the world 31% (Australia, Canada, Lebanon, Malaysia, Philippines, Russian Federation, South Africa, Ukraine).",
    biomarkers = "Median VWF antigen 113 IU/dL (range 43-296), missing in 154 subjects.",
    notes = "Pooled from CSL627_1001 (phase I-III, >= 12 years; NCT01486927), CSL627_3002 (children < 12 years; NCT02093897) and the CSL627_3001 extension (all ages; NCT02172950); 124 subjects who continued into the extension were analysed as one person (Table 1). 3052 OSA and 2545 CSA FVIII levels; 6.9% were below the LLOQ (0.1 IU/dL OSA, 1 IU/dL CSA) and were set to half the LLOQ. Estimated with FOCE-I in NONMEM 7.4.1."
  )

  ini({
    # Structural parameters (CSA scale). Appendix S2 $THETA final estimates;
    # Table 2 prints them rounded.
    lcl <- log(2.07866); label("Clearance at 68 kg and VWF 113 IU/dL (dL/h)") # Appendix S2 THETA(3) = 2.07866; Table 2 CL 2.08 dL/h (95% CI 1.99-2.18)
    lvc <- log(34.3587); label("Central volume of distribution at 68 kg (dL)") # Appendix S2 THETA(1) = 34.3587; Table 2 V1 34.4 dL (95% CI 33.3-35.4)
    lq <- log(1.52841); label("Intercompartmental clearance (dL/h)") # Appendix S2 THETA(4) = 1.52841; Table 2 Q 1.53 dL/h (95% CI 0.96-2.43)
    lvp <- log(2.40027); label("Peripheral volume of distribution (dL)") # Appendix S2 THETA(2) = 2.40027; Table 2 V2 2.40 dL (95% CI 1.71-3.27)

    lrbase <- log(0.265462); label("Endogenous FVIII activity (IU/dL)") # Appendix S2 THETA(7) = 0.265462; Table 2 endogenous FVIII level 0.27 IU/dL (95% CI 0.23-0.31)
    lfprior_dose <- log(1.40884); label("Residual FVIII activity from an unknown pre-study dose at the first record (IU/dL per unit dose into prior_dose)") # Appendix S2 THETA(8) = 1.40884; Table 2 residual FVIII level previous dose 1.41 IU/dL (95% CI 0.69-2.57)

    e_wt_vc <- 0.833885; label("Power exponent of body weight on V1, centred at 68 kg (unitless)") # Appendix S2 THETA(5) = 0.833885; Table 2 0.84 (95% CI 0.79-0.89)
    e_wt_cl <- 0.702087; label("Power exponent of body weight on CL, centred at 68 kg (unitless)") # Appendix S2 THETA(6) = 0.702087; Table 2 0.70 (95% CI 0.61-0.79)
    e_vwf_cl <- -0.652459; label("Power exponent of VWF antigen on CL, centred at 113 IU/dL (unitless)") # Appendix S2 THETA(10) = -0.652459; Table 2 -0.65 (95% CI -0.83 to -0.49)

    e_assay_osa_cc <- 0.52613; label("Ratio of OSA to CSA activity for the drug-derived FVIII (unitless)") # Appendix S2 THETA(9) = 0.52613; Table 2 'OSA instead of CSA' 0.53 (95% CI 0.51-0.54)

    # IIV, log-normal. Table 2 reports CV% = sqrt(exp(omega^2) - 1).
    etalvc ~ 0.0238443 # Appendix S2 $OMEGA 0.0238443 (IIV V1); Table 2 15.5 CV%
    etalcl + etalrbase ~ c(
      0.0611255,
      0.0427114, 1.31873
    ) # Appendix S2 $OMEGA BLOCK(2): IIV CL 0.0611255 (Table 2 25.1 CV%), CL-BASE covariance 0.0427114 (not tabulated), IIV BASE 1.31873 (Table 2 165.6 CV%)
    etalfprior_dose ~ 3.505 # Appendix S2 $OMEGA 3.505 (IIV residual previous dose); Table 2 prints 596.6 CV%, sqrt(exp(3.505) - 1) = 568 CV%
    etaresSd ~ 0.12698 # Appendix S2 $OMEGA 0.12698 (IIV on the residual error, ETA(6)); Table 2 36.8 CV%

    # Residual error. Appendix S2 $SIGMA holds variances; the SDs below are
    # their square roots. CSA rows: Y = IPRED + (IPRED * EPS(1) + IPRED *
    # EPS(5|6) + EPS(3)) * EXP(ETA(6)). OSA rows: Y = IPRED + (IPRED**THETA(11)
    # * EPS(2) + IPRED * EPS(7|8|9) + EPS(4)) * EXP(ETA(6)), with EPS(4) 0 FIX.
    propSd_csa <- 0.0717229; label("Proportional residual SD, CSA (fraction)") # Appendix S2 $SIGMA 0.00514417 (EPS(1)); sqrt = 0.0717; Table 2 'Proportional error CSA' 7.17 CV%
    addSd_csa <- 0.739159; label("Additive residual SD, CSA (IU/dL)") # Appendix S2 $SIGMA 0.546356 (EPS(3)); sqrt = 0.739; Table 2 'Additive error CSA' prints the variance, 0.55
    propSd_rep_csa <- 0.133594; label("Proportional residual SD of repeated measurements of the same sample, CSA (fraction)") # Appendix S2 $SIGMA 0.0178472 (EPS(5), EPS(6) SAME); sqrt = 0.1336; Table 2 prints 13.3 against the OSA row
    powSd_osa <- 0.415944; label("Residual SD scale of the power error term, OSA (IU/dL^(1 - powExp_osa))") # Appendix S2 $SIGMA 0.173009 (EPS(2)); sqrt = 0.416; Table 2 'Scale proportional error OSA' 41.6 CV%
    powExp_osa <- 0.726637; label("Exponent of the power error term, OSA (unitless)") # Appendix S2 THETA(11) = 0.726637; Table 2 'Exponent proportional error OSA' 0.73
    propSd_rep_osa <- 0.213505; label("Proportional residual SD of repeated measurements of the same sample, OSA (fraction)") # Appendix S2 $SIGMA 0.0455845 (EPS(7), EPS(8), EPS(9) SAME); sqrt = 0.2135; Table 2 prints 21.4 against the CSA row
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 68)^e_wt_cl * (VWF / 113)^e_vwf_cl
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Terminal (beta) disposition rate constant (Appendix S2 $PK).
    beta_term <- (kel + k12 + k21 - sqrt((kel + k12 + k21)^2 - 4 * kel * k21)) / 2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(prior_dose) <- -beta_term * prior_dose

    f(prior_dose) <- exp(lfprior_dose + etalfprior_dose)

    # Drug-derived FVIII on the CSA scale; OSA reads e_assay_osa_cc times it.
    fviii_csa <- central / vc
    Cc <- rbase + fviii_csa * e_assay_osa_cc^ASSAY_OSA + prior_dose

    # Total residual SD per row: the components are independent, so their
    # variances add; ETA(6) then scales the whole error.
    assay_csa <- 1 - ASSAY_OSA
    var_csa <- (propSd_csa^2 + propSd_rep_csa^2) * Cc^2 + addSd_csa^2
    var_osa <- powSd_osa^2 * Cc^(2 * powExp_osa) + propSd_rep_osa^2 * Cc^2
    resSd <- exp(etaresSd) * sqrt(assay_csa * var_csa + ASSAY_OSA * var_osa)
    Cc ~ add(resSd)
  })
}
