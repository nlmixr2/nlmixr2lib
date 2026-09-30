Hanke_2021_rosuvastatin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous and fasted oral",
    "rosuvastatin in healthy adult men (Hanke 2021, Electronic Supplementary",
    "Material section 2.4), fitted in NONMEM to individual profiles from two",
    "10 mg oral tablet studies plus the digitised mean profile of one 8 mg",
    "4-h intravenous infusion study. The slow absorption and late Cmax are",
    "described by splitting each oral dose into two portions that share one",
    "first-order absorption rate constant: the first portion (fraction",
    "1 - VF2, typical 63.4 percent) is absorbed from depot without delay and",
    "the second (fraction VF2, typical 36.6 percent) from depot2 after a lag",
    "time of 2.30 h. Total oral bioavailability is 7.93 percent and",
    "elimination is first order from the central compartment. Random effects",
    "are between-subject variability on VF2 (logit scale), CL and the",
    "bioavailability of the first oral portion; there are no covariates. This",
    "is the rosuvastatin-monotherapy population-PK analysis the paper used to",
    "derive the split-dose input for its PK-Sim whole-body PBPK model; the",
    "PBPK layer and the drug-drug-interaction extension of the population-PK",
    "model (whose interaction factors are not reported) are not reproduced here."
  )
  reference <- paste(
    "Hanke N, Gomez-Mantilla JD, Ishiguro N, Stopfer P, Nock V.",
    "Physiologically Based Pharmacokinetic Modeling of Rosuvastatin to",
    "Predict Transporter-Mediated Drug-Drug Interactions.",
    "Pharm Res. 2021;38(10):1645-1661. doi:10.1007/s11095-021-03109-6.",
    "Population-PK parameters from Electronic Supplementary Material Table",
    "S2.4.1 ('without DDI' column); structure and residual-error model from",
    "the NONMEM control stream in ESM section 2.4.2.",
    sep = " "
  )
  vignette <- "Hanke_2021_rosuvastatin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # An oral dose is given as TWO dose records of the full total rosuvastatin
  # dose, one to depot and one to depot2; f(depot) and f(depot2) then route
  # (1 - VF2) * FTOT * exp(ETA(3)) and VF2 * FTOT of it, as in the NONMEM $PK
  # block (F1 = VF1 * FTOT * EXP(ETA(3)), F2 = VF2_2 * FTOT). Intravenous
  # doses go to central.
  dosing <- c("depot", "depot2", "central")

  compartmentData <- list(
    depot = list(analyte = "rosuvastatin", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "rosuvastatin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rosuvastatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rosuvastatin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_studies = 3,
    n_subjects = 54,
    n_profiles = 45,
    age_range = "20-55 years",
    weight_range = "67-105 kg",
    sex_female_pct = 0,
    race_ethnicity = "European",
    disease_state = "healthy volunteers",
    dose_range = "10 mg oral tablet (fasted, single dose); 8 mg intravenous 4-h infusion",
    regions = "Europe",
    unit_of_analysis = "individual concentration-time profiles (oral studies) plus one digitised study-mean profile (intravenous study)",
    notes = paste(
      "Hanke 2021 ESM Table S2.3.1: Martin et al. 2003c, 8 mg intravenous",
      "4-h infusion, n = 10, entered as 1 mean profile (age 21-51 years,",
      "weight 68-85 kg); Stopfer et al. 2016, 10 mg oral tablet fasted, 19",
      "individual profiles (23-49 years, 68-99 kg); Stopfer et al. 2018b, 10",
      "mg oral tablet fasted, 25 individual profiles (20-55 years, 67-105",
      "kg). ESM Table S3.2.1 lists all three studies as 100 percent male",
      "European cohorts. The same ESM section later adds three",
      "drug-drug-interaction studies (rifampicin, probenecid, gemfibrozil)",
      "for a 'with DDI' refit; that refit is not the model encoded here.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # Absorption and bioavailability (ESM Table S2.4.1 'without DDI' column;
    # control stream section 2.4.2 $PK). KA1 is shared by both oral
    # portions; ALAG2 delays only the second portion; FTOT is the absolute
    # bioavailability applied to the whole oral dose.
    # ========================================================================
    lka <- log(0.464); label("First-order absorption rate constant, both oral portions (1/h)") # Table S2.4.1 'Ka' = 0.464 1/h (RSE 4.8%)
    lfdepot <- log(0.0793); label("Absolute oral bioavailability FTOT (fraction)") # Table S2.4.1 'Ftot' = 7.93% (RSE 10.8%)
    logitfrac <- log(0.366 / (1 - 0.366)); label("Logit of the fraction of the total oral dose attributed to the second portion VF2 (unitless)") # Table S2.4.1 'VF2' = 0.366 (RSE 7.0%); $PK PHI_2 = LOG(VF2/(1-VF2))
    ltlag2 <- log(2.30); label("Lag time of the second oral portion ALAG2 (h)") # Table S2.4.1 'ALAG2' = 2.30 h (RSE 1.5%)

    # ========================================================================
    # Disposition: two-compartment model with first-order elimination from
    # central (ESM section 2.4.1; $DES K30 = CL/V3, K34 = Q/V3, K43 = Q/V4).
    # ========================================================================
    lcl <- log(19.1); label("Clearance (L/h)") # Table S2.4.1 'CL' = 19.1 l/h (RSE 4.2%)
    lvc <- log(79.4); label("Central volume of distribution (L)") # Table S2.4.1 'V3' = 79.4 l (RSE 1.7%)
    lq <- log(12.0); label("Intercompartmental clearance (L/h)") # Table S2.4.1 'Q' = 12.0 l/h (RSE 7.5%)
    lvp <- log(199); label("Peripheral volume of distribution (L)") # Table S2.4.1 'V4' = 199 l (RSE 11.8%)

    # ========================================================================
    # Between-subject variability. ESM section 2.3.2: 'IIVs were modeled
    # exponentially'; the control stream puts ETA(2) on CL as
    # THETA * EXP(ETA), ETA(1) on the logit of VF2 (PHI_2 + ETA(1)) and ETA(3)
    # on F1 only (F1 = VF1 * FTOT * EXP(ETA(3)); F2 = VF2_2 * FTOT carries no
    # ETA(3)), although Table S2.4.1 labels ETA(3) 'IIV Ftot'. Table S2.4.1
    # reports each IIV as a %CV; the ESM uses the log-normal relation
    # CV = sqrt(exp(omega^2) - 1) elsewhere (Table S7.0.1 footnote: '35 % CV
    # was assumed (= 1.40 GeoSD)'), so omega^2 = log(1 + CV^2). The same
    # conversion is applied to the logit-scale VF2 row. The $OMEGA block is
    # diagonal and its printed numbers are initial estimates, not finals.
    # ========================================================================
    etalogitfrac ~ log(1 + 0.778^2) # Table S2.4.1 'IIV VF2' = 77.8 %CV (RSE 14.4%)
    etalcl ~ log(1 + 0.261^2) # Table S2.4.1 'IIV CL' = 26.1 %CV (RSE 22.7%)
    etalfdepot ~ log(1 + 0.801^2) # Table S2.4.1 'IIV Ftot' = 80.1 %CV (RSE 9.0%); acts on F1 only per $PK

    # ========================================================================
    # Residual error: $ERROR Y = IPRED + IPRED * EPS(1) + EPS(2), two
    # independent epsilons (combined proportional plus additive). The
    # additive $SIGMA 0.00006 carries the FIX flag, so its control-stream
    # value is the final value; SD = sqrt(0.00006) = 0.0077460 ng/mL, which
    # Table S2.4.1 prints as 'Add RE 0.00775 +- ng/ml, fixed'.
    # ========================================================================
    propSd <- 0.222; label("Proportional residual SD (fraction)") # Table S2.4.1 'Prop RE' = 22.2% (RSE 5.8%)
    addSd <- fixed(0.0077460); label("Additive residual SD (ng/mL)") # $SIGMA EPS(2) = 0.00006 with FIX flag; sqrt(0.00006) = 0.0077460; Table S2.4.1 'Add RE' = 0.00775 ng/ml
  })

  model({
    # ---- Individual parameters (control stream $PK) ------------------------
    ka <- exp(lka)
    tlag2 <- exp(ltlag2)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)

    # VF2_2 = EXP(PHI_2 + ETA(1)) / (1 + EXP(PHI_2 + ETA(1)))
    frac <- expit(logitfrac + etalogitfrac)

    # ---- Micro-constants ($DES) --------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system ($DES DADT(1)-DADT(4)) ---------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka * depot2
    d/dt(central) <- ka * depot + ka * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Bioavailability and lag ($PK F1, F2, ALAG2) -----------------------
    f(depot) <- (1 - frac) * fdepot * exp(etalfdepot)
    f(depot2) <- frac * fdepot
    alag(depot2) <- tlag2

    # ---- Observation and residual error ($ERROR; S3 = V3/1000) -------------
    # Dose in mg and volume in L give mg/L; S3 = V3/1000 converts to ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
