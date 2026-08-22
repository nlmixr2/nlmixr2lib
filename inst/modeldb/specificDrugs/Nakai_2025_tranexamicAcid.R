Nakai_2025_tranexamicAcid <- function() {
  description <- "Two-compartment population PK model for intravenous tranexamic acid (TXA) with first-order elimination, in adults undergoing cardiac surgery with cardiopulmonary bypass; body weight on central volume and Cockcroft-Gault creatinine clearance on clearance (Nakai 2025)."
  reference   <- "Nakai T, Tamura T, Miyagawa Y, Inagaki T, Mutsuga M, Yamada S, Yamada K, Nishiwaki K, Mizoguchi H. Population pharmacokinetic model of tranexamic acid in patients who undergo cardiac surgery with cardiopulmonary bypass. Eur J Clin Pharmacol. 2025;81(3):441-449. doi:10.1007/s00228-025-03802-0"
  vignette    <- "Nakai_2025_tranexamicAcid"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (WT/61.4)^0.911 on the central volume V1 only; 61.4 kg is the cohort median body weight (Nakai 2025 Table 1). Body weight was also tested on CL2 (Table 2 model 3) but the final model (model 7) retains it only on V1. Time-fixed at the pre-operative value.",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation using ACTUAL body weight, NOT body-surface-area normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (CRCL/61.0)^0.752 on CL1; 61.0 mL/min is the cohort median (Nakai 2025 Table 1). Cockcroft-Gault as printed in Nakai 2025 Methods: CLcr (mL/min) = ((140 - age) * BW) / (Scr * 72), x 0.85 if female, with Scr in mg/dL. Serum creatinine was measured pre-operatively, so CRCL is time-fixed at baseline. The summary equation in the Abstract/Results writes 'CLcr (L/h)' but the normalizing constant 61.0 is the mL/min median from Table 1; the Discussion worked example (BW 61.4 kg, CLcr 60.0 mL/min -> CL1 = 3.223 L/h) confirms the ratio is formed in mL/min. Raw (non-BSA-normalized) Cockcroft-Gault CRCL follows the precedent of Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R and Wada_2023_sparsentan.R.",
      source_name        = "CLcr"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tranexamic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 1L,
    n_observations = 453L,
    age_range      = "26-84 years",
    age_median     = "69 years (IQR 60-75)",
    weight_range   = "38.0-97.6 kg",
    weight_median  = "61.4 kg (IQR 54.6-75.2)",
    bmi_range      = "15.1-35.0 kg/m^2",
    sex_female_pct = 33.8,
    disease_state  = "Adults undergoing cardiac surgery (CABG, single or complex valve, aortic, or combined procedures) with cardiopulmonary bypass",
    renal_function = "Cockcroft-Gault creatinine clearance 21.8-147.5 mL/min (median 61.0, IQR 48.1-76.3); serum creatinine 38.0-271.4 umol/L (median 85.7). Patients on haemodialysis were excluded.",
    dose_range     = "Double-bolus regimen: 1000 mg intravenous TXA at the start of the operative procedure and a further 1000 mg after cardiopulmonary bypass was discontinued (Nakai 2025 Methods, Study design)",
    regions        = "Single centre, Nagoya University Hospital, Japan (August 2021 to August 2022)",
    notes          = "Prospective observational study; 88 patients screened, 11 excluded (8 on haemodialysis, 2 without CPB, 1 with a different dosing regimen). Baseline demographics in Nakai 2025 Table 1. Sampling at approximately 0.5, 1, 2 and 5 h after the first dose and 1, 6 and 16 h after re-administration. Intraoperative CPB duration 76-424 min (median 170). Estimation in Phoenix NLME 8.4.0 with FOCE-ELS."
  )

  ini({
    # Structural parameters. Nakai 2025 Table 3 "Final model" Estimate column,
    # reported for the reference subject at the cohort median WT 61.4 kg and
    # median CRCL 61.0 mL/min. Reproduced verbatim in the Results and
    # Conclusions summary equation:
    #   V1 (L) = 12.77 x (BW/61.4)^0.911,  V2 (L) = 6.857,
    #   CL1 (L/h) = 3.263 x [CLcr/61.0]^0.752,  CL2 (L/h) = 2.859
    lvc <- log(12.77); label("Central volume of distribution V1 at WT 61.4 kg (L)")            # Nakai 2025 Table 3 (%RSE 3.582)
    lvp <- log(6.857); label("Peripheral volume of distribution V2 (L)")                       # Nakai 2025 Table 3 (%RSE 13.82)
    lcl <- log(3.263); label("Clearance CL1 at CRCL 61.0 mL/min (L/h)")                        # Nakai 2025 Table 3 (%RSE 4.198)
    lq  <- log(2.859); label("Inter-compartmental clearance CL2 (L/h)")                        # Nakai 2025 Table 3 (%RSE 13.93)

    # Covariate effects. Both were estimated (not fixed): Table 3 reports a
    # %RSE for each, and Table 2 shows the stepwise OFV drops that retained
    # them (model 4 added BW on V1, model 7 added CLcr on CL1).
    e_wt_vc   <- 0.911; label("Power exponent on (WT/61.4) for V1 (unitless)")                 # Nakai 2025 Table 3 "BW on V1" (%RSE 16.51)
    e_crcl_cl <- 0.752; label("Power exponent on (CRCL/61.0) for CL1 (unitless)")              # Nakai 2025 Table 3 "CLcr on CL1" (%RSE 16.55)

    # Inter-individual variability. Nakai 2025 Methods, Base model, declares
    # exponential IIV on "the parameters in the PK model" -- P_i = TV(P) x
    # exp(eta_i), eta ~ N(0, omega^2) -- but NO omega estimate is reported
    # anywhere in the paper. Table 3 lists six rows only (V1, V2, CL1, CL2,
    # BW on V1, CLcr on CL1), and the Results state that "the RSE of all the
    # parameters was 3.582-16.55%", a range whose endpoints are exactly the
    # smallest and largest of those six %RSE values -- so the six structural
    # and covariate estimates are the complete published parameter set. The
    # only supplementary item the paper cites (Supplemental Table S1) is the
    # characteristics of three low-concentration patients, not a parameter
    # table. The etas are therefore encoded as fixed(0): the structure the
    # authors declared is recorded, and the zero states that no variance was
    # published rather than that variability was estimated to be absent.
    # Never substitute an invented or class-typical omega.
    etalvc ~ fixed(0)  # Nakai 2025 Methods, Base model (exponential IIV declared; magnitude not reported)
    etalvp ~ fixed(0)  # Nakai 2025 Methods, Base model (exponential IIV declared; magnitude not reported)
    etalcl ~ fixed(0)  # Nakai 2025 Methods, Base model (exponential IIV declared; magnitude not reported)
    etalq  ~ fixed(0)  # Nakai 2025 Methods, Base model (exponential IIV declared; magnitude not reported)

    # Residual error. Nakai 2025 Methods, Base model, specifies a combined
    # additive-plus-multiplicative error model with SDs sigma1 (additive) and
    # sigma2 (multiplicative), and the Results state the combined error model
    # was selected. Neither sigma is reported in Table 3 (see the note above),
    # so both are encoded as fixed(0) to record the declared structure without
    # inventing a magnitude.
    addSd  <- fixed(0); label("Additive residual error (ug/mL; ZERO - declared by the source but not reported)")        # Nakai 2025 Methods, Base model (sigma1 not reported)
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - declared by the source but not reported)") # Nakai 2025 Methods, Base model (sigma2 not reported)
  })

  model({
    # Individual parameters. Body weight acts on V1 only and creatinine
    # clearance on CL1 only; V2 and CL2 carry no covariate in the final model
    # (Nakai 2025 Table 2 model 7, and the summary equation in Results).
    vc <- exp(lvc + etalvc) * (WT / 61.4)^e_wt_vc
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl) * (CRCL / 61.0)^e_crcl_cl
    q  <- exp(lq  + etalq)

    # Micro-constants for the explicit two-compartment system.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous administration: doses enter the central compartment directly
    # (Nakai 2025 Methods -- 1000 mg IV bolus at the start of surgery and again
    # after CPB was discontinued). No depot, no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give mg/L, which is numerically identical to the
    # ug/mL used throughout Nakai 2025.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
