Carral_2023_amoxicillin <- function() {
  description <- "One-compartment population PK simulation model for oral amoxicillin in virtual adult male outpatients with community-acquired pneumonia (Carral 2023). Parallel first-order (depot) and lagged zero-order (central) absorption, with the zero-order arm's bioavailability falling with dose as KD / (KD + DOSE); apparent clearance proportional to Cockcroft-Gault creatinine clearance (CL/F = CL_slope * CLcr / 102) and apparent volume proportional to body weight. The paper's own Monte Carlo simulation model for adherence / relative-forgiveness analysis, built from the Sjovall data."
  reference <- "Carral N, Lukas JC, Estrade O, Jauregizar N, Morillas H, Suarez E. Non-Adherence in Adult Male Patients with Community-Acquired Pneumonia: Relative Forgiveness of Amoxicillin versus Respiratory Fluoroquinolones. Antibiotics (Basel). 2023;12(5):838. doi:10.3390/antibiotics12050838"
  vignette <- "Carral_2023_cap_nonadherence"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the Cockcroft-Gault creatinine clearance. In the simulation control stream age is generated per virtual patient as SIMAGE = 45 * exp(ETA(6)) with omega^2 = 0.025 (Appendix A); here it is a data column so a user can supply real ages, and the vignette reproduces the paper's draw.",
      source_name = "SIMAGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Used in the Cockcroft-Gault creatinine clearance and as a linear multiplier on V/F (V/F = 0.33 L/kg * WT). Generated in the control stream as SIMWT = 70 * exp(ETA(5)) with omega^2 = 0.001 (Appendix A).",
      source_name = "SIMWT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the Cockcroft-Gault creatinine clearance CLcr = (140 - AGE) * WT / (72 * CREAT) (male form, no 0.85 factor; Appendix A). The paper's three renal subgroups use 0.7, 0.9 and 1.3 mg/dL (Methods 4.1; Supplementary Table S1), giving CLcr of about 131, 100 and 70 mL/min.",
      source_name = "Cr"
    ),
    DOSE = list(
      description = "Administered amoxicillin dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-dependent bioavailability of the zero-order absorption arm, DF = 1 - DOSE / (KD + DOSE) (Appendix A $INPUT column DOSE). Supply it on every record; the paper simulates 1000 mg every 8 h.",
      source_name = "DOSE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "amoxicillin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "31-65 years (virtual; typical 45 years)",
    weight_range = "typical 70 kg (virtual)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "Virtual adult male outpatients with community-acquired pneumonia (CAP) caused by Streptococcus pneumoniae",
    dose_range = "Oral amoxicillin 1000 mg every 8 h for 10 days",
    renal_function = "Three virtual subgroups by serum creatinine 0.7 / 0.9 / 1.3 mg/dL: Cockcroft-Gault CLcr 131 +/- 10, 100 +/- 8 and 70 +/- 6 mL/min (Supplementary Table S1)",
    notes = "Simulation-only study: no patients were enrolled. The PK model was developed by the authors from the Sjovall et al. amoxicillin data (Methods 4.2) and used in NONMEM 7.4.3 ONLYSIM mode with 10,000 virtual patients per renal subgroup. Virtual patients were adult Caucasian men aged 31-65 years weighing 70 kg (Methods 4.1)."
  )

  ini({
    # Appendix A control stream ($THETA). The stream is run with $SIMUL
    # ONLYSIM; its $THETA values are the simulation model's typical values.
    lka <- log(0.635); label("First-order absorption rate constant KA (1/h)") # Appendix A THETA(1) = 0.635 '1.KA (1/hr)'
    lcl <- log(15.5); label("Apparent clearance CL/F at CLcr = 102 mL/min (L/h)") # Appendix A THETA(2) = 15.5; TVCL = THETA(2) * CLcr / 102 (Methods 4.2 'CL/F = CL_slope x ClCR/102')
    lvc <- log(0.33); label("Apparent volume of distribution V/F per kg body weight (L/kg)") # Appendix A THETA(3) = 0.33; TVV = THETA(3) * SIMWT
    logitffo <- log(0.177); label("Logit of the first-order absorbed fraction F1 (unitless)") # Appendix A THETA(5) = 0.177; F1 = 1 - 1/(1 + THETA(5)) = expit(log(THETA(5))) = 0.150
    ld2 <- log(1.44); label("Duration of the zero-order absorption arm D2 (h)") # Appendix A THETA(6) = 1.44; D2 = THETA(6)
    ltlag <- log(0.190); label("Absorption lag of the zero-order arm ALAG2 (h)") # Appendix A THETA(7) = 0.190; ALAG2 = ALAG1 + THETA(7) with ALAG1 = THETA(4) = 0 FIXED
    led50 <- log(1300); label("KD: dose halving the zero-order-arm bioavailability (mg)") # Appendix A THETA(8) = 1300; DF = 1 - DOSE/(KD + DOSE)

    # Appendix A $OMEGA, in ETA order. The first entry is printed as
    # '; 004'; read as 0.04 (20% CV), the CV Methods 4.2 states for the PK
    # parameters. ETA(5) (weight) and ETA(6) (age) generate the virtual
    # covariates and are reproduced in the vignette cohort, not here.
    etalka ~ 0.04 # Appendix A $OMEGA ETA(1), printed '; 004'
    etalcl ~ 0.036 # Appendix A $OMEGA ETA(2)
    etalvc ~ 0.01 # Appendix A $OMEGA ETA(3)
    etaled50 ~ 0.04 # Appendix A $OMEGA ETA(4)

    # Appendix A $ERROR: W = IPRE + 0.001; Y = F + W * EPS(1), $SIGMA 0.4
    # (variance). SD = sqrt(0.4) * (IPRE + 0.001): combined1 with
    # propSd = sqrt(0.4) and addSd = 0.001 * sqrt(0.4).
    propSd <- sqrt(0.4); label("Proportional residual error (fraction)") # Appendix A $SIGMA 0.4 (variance)
    addSd <- 0.001 * sqrt(0.4); label("Additive residual error (mg/L)") # Appendix A W = IPRE + .001 with $SIGMA 0.4
  })

  model({
    # Cockcroft-Gault creatinine clearance, male form (Appendix A).
    crcl_cg <- (140 - AGE) * WT / (CREAT * 72)

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * crcl_cg / 102
    vc <- exp(lvc + etalvc) * WT
    kel <- cl / vc

    d2 <- exp(ld2)
    tlag <- exp(ltlag)
    kd <- exp(led50 + etaled50)
    ffo <- expit(logitffo)
    # Bioavailability of the zero-order arm: the remaining (1 - F1),
    # reduced by the dose-dependent saturation factor DF (fsat).
    fsat <- 1 - DOSE / (kd + DOSE)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Each administration is two dose records with the same amt: one to
    # depot (first-order arm, no lag: ALAG1 = 0 FIXED) and one to central
    # with rate = -2 (zero-order arm of modelled duration D2).
    f(depot) <- ffo
    f(central) <- (1 - ffo) * fsat
    dur(central) <- d2
    alag(central) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
