Wu_2022_caspofungin <- function() {
  description <- "Two-compartment population PK model for intravenous caspofungin in critically ill adults after cardiac surgery, pooling heart transplant recipients and non-transplant cardiac-surgery patients (Wu 2022). Serum albumin enters as a power covariate (centred at 37.42 g/L) on clearance; heart transplantation, ECMO and CRRT were screened and not retained. Exponential IIV on CL, Vc and Vp (Q IIV fixed to 0); residual error is the linear sum of an additive and a proportional SD. NONMEM 7.3 FOCE-I; 58 patients, 414 samples."
  reference <- "Wu Z, Lan J, Wang X, Wu Y, Yao F, Wang Y, Zhao B, Wang Y, Chen J, Chen C. Population pharmacokinetics of caspofungin and dose simulations in heart transplant recipients. Antimicrob Agents Chemother. 2022;66(5):e02249-21. doi:10.1128/aac.02249-21. PMC9116478."
  vignette <- "Wu_2022_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Caspofungin was given as a 1-h intravenous infusion and total plasma
  # caspofungin was assayed by LC-MS/MS with caspofungin-d4 as the internal
  # standard (Methods, 'Analytical procedures'; linear range 0.4-25 mg/L).
  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Evaluated on the day of sample collection (Methods, 'Demographic",
        "characteristics and data collection'), so it may be time-varying.",
        "Enters as (ALB / 37.42)^-1.01 on clearance, the divisor being the",
        "population average used in the supplementary NONMEM control stream",
        "($PK: CL = THETA(1) * EXP(ETA(1)) * (ALB/37.42)**THETA(7)). Table 1",
        "medians: 40.3 g/L (29.7-49.3) in the heart-transplant group and",
        "33.7 g/L (26.41-48.40) in the control group. Must be strictly",
        "positive (power term)."
      ),
      source_name = "ALB"
    )
  )

  # Screened but not retained (Methods, 'Population PK model'; Results;
  # Discussion; Figure S2). Documentation only; none is referenced in model().
  covariatesDataExcluded <- list(
    TX_ANY = list(
      description = "Transplant-recipient indicator (1 = heart transplant recipient, 0 = non-transplant cardiac-surgery patient)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'HTx status', not retained -- the paper's headline negative finding. Every transplant in this cohort is a heart transplant (27 of 58 patients)."
    ),
    ECMO_STATUS = list(
      description = "ECMO treatment-status indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. 8/27 heart-transplant and 8/31 control patients were on ECMO (Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy status indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (with CRRT blood flow rate and dialysate), not retained. 12/27 and 16/31 patients on CRRT (Table 1)."
    ),
    AGE = list(description = "Age", units = "years", type = "continuous", notes = "Screened, not retained (Table 1)."),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened, not retained; the Discussion notes the narrow weight range (median 60 kg, 43.5-100)."
    ),
    HT = list(description = "Height", units = "cm", type = "continuous", notes = "Screened, not retained (Table 1)."),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; dropped OFV by 11.93 alone but added nothing after ALB (collinearity, Discussion)."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    PROCALCITONIN = list(
      description = "Procalcitonin",
      units = "ug/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained (Table 1, 'Serum creatine')."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened, not retained. Raw mL/min (Table 1)."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score at baseline",
      units = "(score)",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    n_samples = 414L,
    age_range = "Median 50 years (20-73) in the heart-transplant group and 58 (22-78) in the control group (Table 1)",
    weight_range = "Median 59.5 kg (43.5-76) in the heart-transplant group and 62.0 kg (48.0-100.0) in the control group (Table 1)",
    sex_female_pct = 25.9,
    race_ethnicity = "Not reported; single centre in Guangzhou, China",
    disease_state = "Critically ill adults in the surgical ICU after cardiac surgery who received caspofungin: 27 heart transplant recipients and 31 non-transplant cardiac-surgery controls. 28 were on CRRT and 16 on ECMO; median SOFA 9-10.",
    dose_range = "70 mg loading dose then 50 mg every 24 h, as a 1-h intravenous infusion (Methods, 'Blood sample collection').",
    regions = "China (Guangdong Provincial People's Hospital, Guangzhou; July 2019 - July 2021)",
    notes = paste(
      "Sampling: pre-dose and 1, 2, 6, 10, 16 and 24 h after the start of the",
      "infusion, on at least 3 days. One-, two- and three-compartment models",
      "were tested; two compartments were selected. Final estimates come from",
      "Table 2 and the supplementary NONMEM control stream, which agree."
    )
  )

  ini({
    # Table 2 'Estimates' column (final model); identical to the $THETA
    # values of the supplementary NONMEM control stream. Typical values apply
    # at ALB = 37.42 g/L.
    lcl <- log(0.385)
    label("Clearance at ALB = 37.42 g/L (L/h)") # Table 2: CL = 0.385 L/h (RSE 5%); control stream THETA(1)
    lvc <- log(4.27)
    label("Central volume (L)") # Table 2: Vc = 4.27 L (RSE 12%); control stream THETA(2)
    lq <- log(2.85)
    label("Intercompartmental clearance (L/h)") # Table 2: Q = 2.85 L/h (RSE 11%); control stream THETA(3)
    lvp <- log(6.01)
    label("Peripheral volume (L)") # Table 2: Vp = 6.01 L (RSE 13%); control stream THETA(4)

    e_alb_cl <- -1.01
    label("Power exponent of (ALB / 37.42) on CL (unitless)") # Table 2: theta ALB-CL = -1.01 (RSE 14%); control stream THETA(7)

    # IIV: exponential model. Table 2 'INTER VAR (%)' = 100 * sqrt(omega^2);
    # the control stream $OMEGA block prints the variances (0.112, 0.455,
    # 0 FIX, 0.22).
    etalcl ~ 0.112 # Table 2: CL IIV 33.5% (0.335^2 = 0.112); control stream $OMEGA 0.112
    etalvc ~ 0.455 # Table 2: Vc IIV 67.5% (0.675^2 = 0.456); control stream $OMEGA 0.455
    etalq ~ fixed(0) # Table 2: Q IIV '0 FIX'; control stream $OMEGA '0 FIX'
    etalvp ~ 0.2275 # Table 2: Vp IIV 47.7% (0.477^2 = 0.2275); control stream $OMEGA prints it rounded as 0.22

    # $ERROR: W = IPRED * THETA(5) + THETA(6); Y = IPRED + W * ERR(1) with
    # $SIGMA 1 FIX -- the SD is the linear sum of the two terms (combined1).
    propSd <- 0.134
    label("Proportional residual error (fraction)") # Table 2: proportional error 13.4% (RSE 11%); control stream THETA(5)
    addSd <- 0.213
    label("Additive residual error (mg/L)") # Table 2: additive error 0.213 mg/L (RSE 29%); control stream THETA(6)
  })

  model({
    # Supplementary control stream $PK:
    #   CL = THETA(1) * EXP(ETA(1)) * (ALB/37.42)**THETA(7)
    #   V1 = THETA(2) * EXP(ETA(2)); Q = THETA(3) * EXP(ETA(3))
    #   V2 = THETA(4) * EXP(ETA(4)); ADVAN3 TRANS4
    cl <- exp(lcl + etalcl) * (ALB / 37.42)^e_alb_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
