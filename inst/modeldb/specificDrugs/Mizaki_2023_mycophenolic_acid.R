Mizaki_2023_mycophenolic_acid <- function() {
  description <- paste0(
    "Two-compartment population PK model for mycophenolic acid (MPA) after ",
    "oral mycophenolate mofetil (MMF) twice daily in adult Japanese patients ",
    "with lupus nephritis (Mizaki 2023). Absorption through a chain of six ",
    "transit compartments (Savic form, ktr = (6 + 1) / MTT) into an ",
    "absorption compartment drained by first-order ka; first-order ",
    "elimination from the central compartment; no enterohepatic circulation ",
    "term. Cockcroft-Gault creatinine clearance and serum albumin act on CL/F ",
    "as power functions centred on the cohort medians (6.2 L/h and 3.6 g/dL). ",
    "Doses are mg of MMF (no molecular-weight conversion; parameters are ",
    "apparent with respect to the MMF dose) and Cc is MPA in mg/L. IIV is ",
    "exponential on V1/F, CL/F and MTT; residual error is proportional ",
    "(34.09%)."
  )
  reference <- paste(
    "Mizaki T, Nobata H, Banno S, Yamaguchi M, Kinashi H, Iwagaitsu S,",
    "Ishimoto T, Kuru Y, Ohnishi M, Sako K, Ito Y. (2023). Population",
    "pharmacokinetics and limited sampling strategy for therapeutic drug",
    "monitoring of mycophenolate mofetil in Japanese patients with lupus",
    "nephritis. J Pharm Health Care Sci 9:1.",
    "doi:10.1186/s40780-022-00271-w",
    sep = " "
  )
  vignette <- "Mizaki_2023_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Doses are mg of MMF and concentrations ug/mL (= mg/L) of MPA; the paper
  # gives no molecular-weight conversion, so every state holds an amount in
  # MMF-dose-equivalent mg and all volumes and clearances are apparent with
  # respect to the MMF dose (Methods: V1/F, CL/F, V2/F and Q/F are reported
  # as V1, CL, V2 and Q).
  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "mycophenolate mofetil (MMF), absorption compartment",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "mycophenolic acid (MMF-dose-equivalent mg)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "mycophenolic acid (MMF-dose-equivalent mg)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, NOT BSA-normalized).",
      units = "mL/min",
      type = "continuous",
      reference_category = "6.2 L/h (= 103.3 mL/min), the cohort median; Table 1 reports the median as 103.6 mL/min (IQR 71.1-125.0).",
      notes = "Methods: CLcr calculated with the Cockcroft-Gault equation. The model was fitted with CLcr in L/h; Additional file 2 prints the power form CL = tvCL x (CLcr / CLcr_median)^theta with 'CLcr_median was 6.2 L/h'. Converted inside model() as CRCL x 60 / 1000.",
      source_name = "CLcr"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = "36 g/L (3.6 g/dL), the cohort median (Table 1; Additional file 2 'Alb_median was 3.6 g/dL').",
      notes = "The paper reports albumin in g/dL; model() converts the canonical g/L column back to g/dL (ALB x 0.1) so the published exponent applies unchanged.",
      source_name = "Alb"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    CRP = list(
      description = "C-reactive protein.",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    PRED_DOSE = list(
      description = "Concomitant prednisolone dose.",
      units = "mg/day",
      type = "continuous",
      notes = "Tested as a covariate; not retained (Methods)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor (1 = yes).",
      units = "(binary)",
      type = "binary",
      notes = "Exponential effect on V1/F (-1.07) in the four-covariate candidate model of Additional file 2; dropped from the final model for overparameterization (RSE of V1 32.71%, condition number 5128.3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_observations = "302 plasma MPA concentrations (Results): one 12-h steady-state profile per patient at 0, 0.5, 1, 2, 3, 4, 6, 8 and 12 h after the morning dose (33 full profiles; one patient with 0-4 h only).",
    age_range = "Median 39.0 years (IQR 26.3-51.8); adults 18 years or older (Table 1).",
    weight_range = "Median 54.7 kg (IQR 47.2-58.6) (Table 1).",
    sex_female_pct = 82.4,
    race_ethnicity = c(Asian = 100),
    disease_state = "Lupus nephritis treated with MMF (CellCept) at Aichi Medical University Hospital, March 2015 - June 2022; retrospective chart review.",
    dose_range = "Oral MMF 250-1000 mg every 12 h; median 1500 mg/day (IQR 1500-1688), 27.9 mg/kg/day (Table 1).",
    co_medication = "Prednisolone in all 34 (median 10 mg); tacrolimus 8, proton-pump inhibitor 22, iron/magnesium oxide 5, NSAIDs 3 (Table 1).",
    renal_function = "Cockcroft-Gault creatinine clearance median 103.6 mL/min (IQR 71.1-125.0); eGFR median 77.4 mL/min (Table 1).",
    albumin = "Serum albumin median 3.6 g/dL (IQR 2.8-3.9) (Table 1).",
    assay = "Enzyme immunoassay (cobas MPA kit on cobas 6000 c501), LLOQ 0.40 ug/mL.",
    regions = "Japan (single centre).",
    notes = "Male 6, female 28. Estimation by FOCE-ELS in Phoenix NLME 8.1."
  )

  ini({
    lka <- log(2.98); label("Absorption rate constant ka (1/h)") # Table 2 'Ka (h-1)' 2.98
    lvc <- log(22.95); label("Apparent central volume V1/F (L)") # Table 2 'V1 (L)' 22.95
    lvp <- log(336.03); label("Apparent peripheral volume V2/F (L)") # Table 2 'V2 (L)' 336.03
    lcl <- log(13.15); label("Apparent clearance CL/F at CLcr 6.2 L/h and albumin 3.6 g/dL (L/h)") # Table 2 'CL (L/h)' 13.15
    lq <- log(26.44); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q(L/h)' 26.44
    lmtt <- log(0.45); label("Mean transit time of the absorption chain MTT (h)") # Table 2 'MTT (h)' 0.45
    ntr <- fixed(6); label("Number of transit compartments (count)") # Discussion 'fixed 6 transit compartments'; Results 'The model with 6 compartments'; Fig. 3 caption ktr = (n + 1) / MTT

    e_crcl_cl <- 0.78; label("Power exponent of CLcr / 6.2 L/h on CL/F (unitless)") # Table 2 'Effect CLcr on CL' 0.78; power form from Additional file 2
    e_alb_cl <- -0.88; label("Power exponent of albumin / 3.6 g/dL on CL/F (unitless)") # Table 2 'Effect Alb on CL' -0.88; power form from Additional file 2

    # IIV: Table 2 reports CV% of exponential random effects;
    # omega^2 = log(1 + CV^2). No covariances are reported.
    etalvc ~ 0.56315 # Table 2 'IIV V1 (CV%)' 86.96; log(1 + 0.8696^2)
    etalcl ~ 0.09883 # Table 2 'IIV CL (CV%)' 32.23; log(1 + 0.3223^2)
    etalmtt ~ 0.29433 # Table 2 'IIV MTT (CV%)' 58.50; log(1 + 0.5850^2)

    propSd <- 0.3409; label("Proportional residual error (fraction)") # Table 2 'Proportional error (CV%)' 34.09
  })

  model({
    # Additional file 2 power covariate form, centred on the cohort medians:
    # CL = tvCL x (CLcr / 6.2 L/h)^theta1 x (Alb / 3.6 g/dL)^theta2 x exp(eta).
    crcl_lh <- CRCL * 60 / 1000
    alb_gdl <- ALB * 0.1

    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    cl <- exp(lcl + etalcl) * (crcl_lh / 6.2)^e_crcl_cl * (alb_gdl / 3.6)^e_alb_cl
    q <- exp(lq)
    mtt <- exp(lmtt + etalmtt)

    # Savic transit form, ktr = (n + 1) / MTT: the dose compartment (depot)
    # and the six transit compartments all drain at ktr, so the mean arrival
    # time into the absorption compartment (transit7) equals MTT; transit7
    # then empties into central at ka.
    ktr <- (ntr + 1) / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ka * transit7
    d/dt(central) <- ka * transit7 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
