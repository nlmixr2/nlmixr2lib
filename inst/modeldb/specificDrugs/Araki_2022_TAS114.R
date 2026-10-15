Araki_2022_TAS114 <- function() {
  description <- "Semimechanistic population PK/PD model of oral TAS-114 (dual dUTPase / dihydropyrimidine dehydrogenase inhibitor) in healthy adult men and adults with advanced solid tumours (Araki 2022). TAS-114 follows a two-compartment model with first-order absorption and an absorption lag; clearance is scaled by the relative amount of a metabolising enzyme (CYP3A) whose zero-order synthesis is stimulated by the central TAS-114 concentration through an Emax function (enzyme-turnover autoinduction). Age (power) and AST (exponential) act on CL/F and body surface area (exponential) on Vc/F. Plasma uracil, the endogenous DPD substrate, follows an indirect-response model in which TAS-114 inhibits the first-order uracil elimination (Imax model) with the elimination rate fixed from the literature uracil half-life."
  reference <- "Araki H, Takenaka T, Takahashi K, Yamashita F, Matsuoka K, Yoshisue K, Ieiri I. A semimechanistic population pharmacokinetic and pharmacodynamic model incorporating autoinduction for the dose justification of TAS-114. CPT Pharmacometrics Syst Pharmacol. 2022;11(5):604-615. doi:10.1002/psp4.12747"
  vignette <- "Araki_2022_TAS114"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F, (AGE / 59)^e_age_cl, centred on the data-set median of 59 years (Table 2; Text S1 'CLAGE = ((AGE/59)**THETA(3))').",
      source_name = "AGE"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on CL/F, exp(e_ast_cl * (AST - 24)), centred on the data-set median of 24 U/L (Table 2; Text S1 'CLAST = EXP(THETA(2)*(AST - 24))').",
      source_name = "AST"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on Vc/F, exp(e_bsa_vc * (BSA - 1.7)), centred on the data-set median of 1.70 m^2 (Table 2; Text S1 'V2BSA = EXP(THETA(1)*(BSA - 1.7))'). The BSA formula is not stated.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "TAS-114", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "TAS-114", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "TAS-114", units = "mg", specimen = "tissue", verified = TRUE),
    enzyme = list(
      analyte = "CYP3A (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    uracil = list(analyte = "uracil", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  # Plasma uracil (the endogenous DPD substrate) is carried as its own
  # indirect-response state; it is the only non-canonical compartment.
  paper_specific_compartments <- c("uracil")

  population <- list(
    species = "human",
    n_subjects = 185L,
    n_studies = 4L,
    age_range = "20-81 years",
    age_median = "59 years",
    weight_range = "36-119 kg",
    weight_median = "62 kg",
    bsa_range = "1.25-2.35 m^2",
    bsa_median = "1.70 m^2",
    sex_female_pct = 42.7,
    race_ethnicity = c(Japanese = 51.9, Caucasian = 41.6, `African American` = 3.8, Other = 2.7),
    disease_state = "Healthy adult men (study 10057010, n = 28) and adults with advanced solid tumours (studies 10057020, TPU-TAS-114-101 and TPU-TAS-114-102, n = 157)",
    dose_range = "6-800 mg TAS-114 orally, single dose or twice daily for 14 days (with S-1 or capecitabine in the patient studies)",
    co_medication = "S-1 (studies 10057020 and TPU-TAS-114-102) or capecitabine (TPU-TAS-114-101); none in study 10057010",
    pd_subset = "Uracil PD model fit to 240 plasma uracil concentrations from the 24 healthy men of study 10057010",
    notes = "Demographics from Table 2 (sex 106 male / 79 female; race 96 Japanese / 77 Caucasian / 7 African American / 5 other; ALT 19 (6-151) U/L; AST 24 (7-140) U/L; BUN 16.0 (5.7-354.2) mg/dL). Study designs in Table 1. 2661 plasma TAS-114 concentrations were modelled."
  )

  ini({
    # ------------------------------------------------------------------
    # TAS-114 PK (Table 3, final population PK model; NONMEM code in
    # supplementary Text S1). Values are the Table 3 final estimates; the
    # Text S1 $THETA block carries a mix of final and initial values.
    # ------------------------------------------------------------------
    lka <- log(0.508); label("First-order absorption rate constant (1/h)") # Table 3 ka = 0.508 1/h (RSE 3.6%)
    lvc <- log(17.0); label("Apparent central volume of distribution Vc/F (L)") # Table 3 Vc/F = 17.0 L (RSE 7.0%)
    lvp <- log(10.5); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 3 Vp/F = 10.5 L (RSE 5.7%)
    lq <- log(4.02); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 Q/F = 4.02 L/h (RSE 13.1%)
    lcl <- log(8.74); label("Apparent clearance CL/F at baseline enzyme amount (L/h)") # Table 3 CL/F = 8.74 L/h (RSE 6.2%)
    ltlag <- log(0.217); label("Absorption lag time (h)") # Table 3 tlag = 0.217 h (RSE 8.4%)

    # Enzyme-turnover autoinduction (Methods equations; Text S1 $DES).
    lemax <- log(4.69); label("Maximum stimulation of enzyme synthesis Emax (unitless)") # Table 3 Emax = 4.69 (RSE 24.0%)
    lec50 <- log(5870); label("TAS-114 concentration at half-maximal enzyme induction EC50 (ng/mL)") # Table 3 EC50 = 5870 ng/mL (RSE 28.2%)
    lkenz <- log(0.0230); label("First-order enzyme degradation rate constant kenz,deg (1/h)") # Table 3 kenz,deg = 0.0230 1/h (RSE 16.4%; t1/2 30.1 h)

    # Covariate effects (Table 3; functional forms from Text S1 $PK).
    e_bsa_vc <- 1.52; label("Exponential BSA effect on Vc/F (1/m^2), centred at 1.70 m^2") # Table 3 'Effects of BSA on Vc/F' = 1.52 (RSE 21.4%)
    e_ast_cl <- -0.00753; label("Exponential AST effect on CL/F (1/(U/L)), centred at 24 U/L") # Table 3 'Effects of AST on CL/F' = -0.00753 (RSE 19.2%)
    e_age_cl <- -0.983; label("Power exponent of age on CL/F (unitless), reference 59 years") # Table 3 'Effects of AGE on CL/F' = -0.983 (RSE 15.6%)

    # IIV (exponential). Table 3 reports CV% = sqrt(exp(omega^2) - 1); the
    # variances below are log(1 + CV^2). They agree with the Text S1 $OMEGA
    # values (0.0432399, 0.39638, 0.357514) to the printed precision.
    etalka ~ 0.04316 # Table 3 IIV ka = 21.0 CV%
    etalvc ~ 0.39692 # Table 3 IIV Vc/F = 69.8 CV%
    etalcl ~ 0.35882 # Table 3 IIV CL/F = 65.7 CV%

    # Proportional residual error. Text S1 codes W = THETA(13)*IPRED with
    # $SIGMA 1 FIX, so THETA(13) is the SD; its value 0.37177 is the square
    # of the 0.610 that Table 3 prints as 'Proportional error (CV%) 61.0',
    # i.e. the table reports sqrt(THETA). The Figure 2 pcVPC bands support
    # the control-stream reading (see the vignette).
    propSd <- 0.37177; label("Proportional residual error, TAS-114 (fraction)") # Text S1 $THETA 13 = 0.37177 (Table 3 shows sqrt = 61.0%)

    # ------------------------------------------------------------------
    # Uracil indirect-response PD (Table 4, final population PD model;
    # NONMEM code in supplementary Text S2).
    # ------------------------------------------------------------------
    lkout <- fixed(log(2.67)); label("First-order uracil elimination rate constant kout (1/h)") # Table 4 kout = 2.67 1/h, not estimated (literature uracil t1/2 0.260 h)
    lic50 <- log(1046); label("TAS-114 concentration at half-maximal inhibition of uracil elimination IC50 (ng/mL)") # Table 4 IC50 = 1046 ng/mL (RSE 21.7%)
    limax <- log(0.888); label("Maximum fractional inhibition of uracil elimination Imax (unitless)") # Table 4 Imax = 0.888 (RSE 8.5%)
    lrbase <- log(9.27); label("Baseline plasma uracil concentration (ng/mL)") # Table 4 Baseline = 9.27 ng/mL (RSE 3.5%)

    etalrbase ~ 0.02344 # Table 4 IIV baseline = 15.4 CV%; log(1 + 0.154^2)

    # Text S2 codes the uracil error the same way as Text S1 (W = THETA*IPRED,
    # $SIGMA 1 FIX). Read the same way as the PK error, Table 4's
    # 'Proportional error (CV%) 39.9' is sqrt(THETA), so the SD is 0.399^2.
    # The Figure 2e pcVPC baseline band supports this reading (vignette).
    propSd_uracil <- 0.159; label("Proportional residual error, plasma uracil (fraction)") # Table 4 'Proportional error (CV%) 39.9' squared: 0.399^2 = 0.159
  })

  model({
    # Individual PK parameters (Text S1 $PK)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * exp(e_bsa_vc * (BSA - 1.7))
    vp <- exp(lvp)
    q <- exp(lq)
    cl <- exp(lcl + etalcl) * (AGE / 59)^e_age_cl * exp(e_ast_cl * (AST - 24))
    tlag <- exp(ltlag)

    emax <- exp(lemax)
    ec50 <- exp(lec50)
    kenz <- exp(lkenz)

    kout <- exp(lkout)
    ic50 <- exp(lic50)
    imax <- exp(limax)
    rbase <- exp(lrbase + etalrbase)
    kin <- kout * rbase

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL (the NONMEM
    # data set dosed in ug, so A(2)/V2 was already in ng/mL).
    Cc <- 1000 * central / vc

    # Emax stimulation of enzyme synthesis (Methods Eq. for f(Cp))
    ind <- emax * Cc / (ec50 + Cc)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * enzyme * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(enzyme) <- kenz * (1 + ind) - kenz * enzyme
    d/dt(uracil) <- kin - kout * (1 - imax * Cc / (ic50 + Cc)) * uracil

    enzyme(0) <- 1
    uracil(0) <- rbase
    alag(depot) <- tlag

    Cc ~ prop(propSd)
    uracil ~ prop(propSd_uracil)
  })
}
