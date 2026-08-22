Kim_2025_infliximab_aubourg <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in adults with Crohn's disease (Aubourg 2015 model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Aubourg A, Picon L, Lecomte T, Bejan-Angoulvant T, Paintaud G, Ternant D.",
    "A robust estimation of infliximab pharmacokinetic parameters in Crohn's",
    "disease. Eur J Clin Pharmacol. 2015;71(12):1541-1542.",
    "doi:10.1007/s00228-015-1942-8 (reference 16 of Kim 2025).",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #1) with MAXEVAL=0,",
    "i.e. every THETA / OMEGA / SIGMA is held at the originally published value.",
    "Values here are transcribed from that control stream and cross-checked",
    "against Kim 2025 Table S2.",
    sep = " "
  )
  vignette <- "Kim_2025_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "no reference level: separate typical CL and Vc are estimated for each sex",
      notes              = "Drives separate typical values of both CL and Vc rather than a multiplicative offset: CL = 0.336 L/day (female) or 0.456 L/day (male); Vc = 2.6 L (female) or 3.2 L (male) before the weight effect. The source control stream codes SEX = 1 for female and SEX = 0 for male, which matches the canonical SEXF polarity directly (no value inversion needed). Time-invariant.",
      source_name        = "SEX"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vc only, normalized to a 60 kg reference: Vc = Vc_sex * (WT/60)^0.22. CL carries no weight effect in this model. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement).",
      source_name        = "WGT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 133L,
    n_studies      = 1L,
    age_range      = "Adults; age not specified in Kim 2025 Table S1 for this model.",
    weight_range   = "Not reported; median body weight 60 kg (Kim 2025 Table S1), which is the reference weight used by the model.",
    sex_female_pct = 59,
    race_ethnicity = "Not specified; developed in a French Crohn's disease cohort.",
    disease_state  = "Crohn's disease (n = 133).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "France.",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Aubourg model: Crohn's disease only (n = 133), adults,",
      "induction and maintenance phases, peak and trough sampling, 0% ATI",
      "positive, 59% female, median body weight 60 kg. Age, albumin, ESR and",
      "immunomodulator use were not specified.",
      "The external-validation cohort in Kim 2025 itself (142 Korean adult IBD",
      "patients with 294 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model gave MPE 29.0%, MAE 1.1 mg/L and rRMSE 200.0%",
      "(Kim 2025 Table 2).",
      sep = " "
    )
  )

  ini({
    # Structural parameters. This model estimates a separate typical clearance
    # and central volume in each sex rather than a reference value plus an
    # offset, so both strata carry an explicit suffix (no bare lcl / lvc).
    lcl_female <- log(0.336); label("Typical clearance in females (CL, L/day)")   # Data S1 model #1 $THETA1 = 0.336 L/day (reported 0.014 L/h); Table S2 "F: 0.336"
    lcl_male   <- log(0.456); label("Typical clearance in males (CL, L/day)")     # Data S1 model #1 $THETA2 = 0.456 L/day (reported 0.019 L/h); Table S2 "M: 0.456"
    lvc_female <- log(2.6);   label("Typical central volume in females at 60 kg (Vc, L)")  # Data S1 model #1 $THETA3 = 2.6 L; Table S2 "F: 2.6"
    lvc_male   <- log(3.2);   label("Typical central volume in males at 60 kg (Vc, L)")    # Data S1 model #1 $THETA4 = 3.2 L; Table S2 "M: 3.2"
    lvp        <- log(4.5);   label("Peripheral volume of distribution (Vp, L)")           # Data S1 model #1 $THETA5 = 4.5 L; Table S2 "4.5"
    lq         <- log(1.992); label("Inter-compartmental clearance (Q, L/day)")            # Data S1 model #1 $THETA6 = 1.992 L/day; Table S2 "1.992"

    # Covariate effect
    e_wt_vc <- 0.22; label("Power exponent of body weight on Vc ((WT/60)^e_wt_vc)")  # Data S1 model #1 $THETA7 = 0.22

    # Inter-individual variability. The source control stream states
    # "$OMEGA ; Interindividual standard deviation OMEGA = (%CV/100)**2", i.e.
    # the variance is the squared coefficient of variation rather than the exact
    # log-normal log(1 + CV^2). Values are used exactly as published:
    #   0.2209 = 0.47^2 (47% CV on CL); 0.0729 = 0.27^2 (27% CV on Vc).
    # A single CL eta and a single Vc eta are shared across both sex strata.
    etalcl ~ 0.2209   # Data S1 model #1 $OMEGA ETA1 = 0.2209 (47% CV); Table S2 '(47%)'
    etalvc ~ 0.0729   # Data S1 model #1 $OMEGA ETA2 = 0.0729 (27% CV); Table S2 '(27%)'

    # Residual error -- combined proportional + additive. The control stream
    # reports variances with the corresponding SDs in its comments; the SDs are
    # what nlmixr2 expects.
    propSd <- 0.21; label("Proportional residual error (fraction)")  # Data S1 model #1 $SIGMA 0.0441, comment "proportional error (SD = 0.21)"
    addSd  <- 2.3;  label("Additive residual error (mg/L)")         # Data S1 model #1 $SIGMA 5.29, comment "additive error (SD = 2.3 ug/ml)"
  })
  model({
    # Individual PK parameters. Sex selects the typical value for CL and Vc;
    # SEXF = 1 is female and SEXF = 0 is male, matching the source SEX coding.
    #   CL (L/day) = 0.336 (female) or 0.456 (male)
    #   Vc (L)     = [2.6 (female) or 3.2 (male)] * (WT/60)^0.22
    #   Vp (L)     = 4.5
    #   Q  (L/day) = 1.992
    cl <- exp(lcl_female * SEXF + lcl_male * (1 - SEXF) + etalcl)
    vc <- exp(lvc_female * SEXF + lvc_male * (1 - SEXF) + etalvc) * (WT / 60)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model; infliximab is given as an IV infusion into central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L, numerically identical to ug/mL.
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
