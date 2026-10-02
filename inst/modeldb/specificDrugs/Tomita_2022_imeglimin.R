Tomita_2022_imeglimin <- function() {
  description <- "Two-compartment population PK model for oral imeglimin in Japanese and Western healthy volunteers and patients with type 2 diabetes, including patients with chronic kidney disease down to eGFR 14 mL/min/1.73 m^2 (Tomita 2022). First-order absorption with a lag time; dose-dependent relative bioavailability (inhibitory Emax on dose, normalised to F = 1 at 1000 mg) and dose-dependent ka (power on dose); formulation and fasting effects on lag time and ka; linear eGFR effect on CL/F capped at 120 mL/min/1.73 m^2 plus power effects of body weight and age on CL/F, body weight on Vc/F, and age and Japanese ethnicity on Q/F. Log-scale residual error with separate magnitudes for phase I studies and for pre-dose samples. Doses must be supplied as imeglimin FREE BASE (labelled imeglimin hydrochloride mass x 0.810, e.g. 1000 mg tablet = 810 mg); the dose nonlinearity on F and ka is evaluated on the labelled dose, recovered in model() as podo(depot) / 0.810."
  reference <- "Tomita Y, Hansson E, Mazuir F, Wellhagen GJ, Ooi QX, Mezzalana E, Kitamura A, Nemoto D, Bolze S. Imeglimin population pharmacokinetics and dose adjustment predictions for renal impairment in Japanese and Western patients with type 2 diabetes. Clin Transl Sci. 2022;15(4):1014-1026. doi:10.1111/cts.13221"
  vignette <- "Tomita_2022_imeglimin"
  paper_specific_residual_sds <- c("expSdPhase23", "expSdPhase1", "expSdPredose")
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate normalised to 1.73 m^2 body surface area, at baseline",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Creatinine-based eGFR: the three-variable Japanese equation for Japanese subjects and the CKD-EPI equation for Western subjects (Tomita 2022 Methods 'Clinical studies'). Enters CL/F linearly, centred at the population median 81.4 and capped at 120 (Equation 8): (1 + (min(CRCL, 120) - 81.4) * e_crcl_cl). Observed range 14.1-152.",
      source_name = "eGFR"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (Equation 9) and Vc/F (Equation 1), normalised to the median 77.35 kg (Table 3 footnote).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (Equation 10) and Q/F (Equation 1), normalised to the median 59 years (Table 3 footnote).",
      source_name = "age"
    ),
    RACE_JAPANESE = list(
      description = "Japanese ethnicity indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = Japanese subject, 0 = Western subject (mostly White; Table 2 note). Fractional effect on Q/F (Equation 2): (1 + e_japanese_q * RACE_JAPANESE).",
      source_name = "Japanese"
    ),
    FORM_CAPSULE = list(
      description = "Imeglimin capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = capsule (Western renal-impairment study and Western phase IIa studies 003 / 004, Table 1); 0 = tablet. The reference formulation is the optimised tablet used in the Japanese phase IIb / III studies (Table 3 footnote a). Mutually exclusive with FORM_IMEGLIMIN_CONVENTIONAL_TABLET.",
      source_name = "Formulation"
    ),
    FORM_IMEGLIMIN_CONVENTIONAL_TABLET = list(
      description = "Imeglimin conventional (non-optimised) tablet formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = conventional tablet (single/multiple ascending dose study and Western phase IIb, Table 1); 0 = capsule or optimised tablet. Reference = optimised tablet (Table 3 footnote a).",
      source_name = "Formulation"
    ),
    FED = list(
      description = "Non-fasting condition at the time of dosing",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "1 = non-fasting (regular meal, high-fat meal, or no specific instruction on meal times); 0 = fasted or semi-fasted (Table 3 footnote b). The paper's reference is the non-fasting condition, so the fasting effects enter as (1 + e_fasted_tlag * (1 - FED)) and (1 + e_fasted_ka * (1 - FED)). Per dose record.",
      source_name = "Food"
    ),
    STUDY_IMEGLIMIN_PHASE1 = list(
      description = "Phase I (clinical pharmacology) study indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = subject from one of the four clinical pharmacology studies (single/multiple ascending dose, Western renal impairment, Japanese renal impairment, T2DM with renal impairment; Table 1, Table 2 'Phase I' column, 160 subjects); 0 = phase II / III study. Selects the phase I residual-error magnitude (Table 3 'Proportional phase I RUV').",
      source_name = "Phase I"
    ),
    SAMPLE_PREDOSE = list(
      description = "Pre-dose plasma sample indicator (per observation)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = observation drawn before the dose at that visit; 0 = post-dose sample. Selects the pre-dose residual-error magnitude (Table 3 'Proportional predose RUV'). The paper does not state which magnitude applies to a pre-dose sample in a phase I study; the model lets the pre-dose magnitude take precedence over the phase I one (an assumption, see the vignette).",
      source_name = "predose"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "imeglimin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "imeglimin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imeglimin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 867L,
    n_studies = 9L,
    n_observations = 8256L,
    age_range = "20-80 years",
    weight_range = "35.6-148 kg",
    sex_female_pct = 43.1,
    race_ethnicity = c(Japanese = 44.4, Western = 55.6),
    disease_state = "Healthy volunteers and patients with type 2 diabetes mellitus (745 of 867), including patients with chronic kidney disease stages G3a-G5",
    renal_function = "eGFR 14.1-152 mL/min/1.73 m^2 (CKD stages G1-G5)",
    dose_range = "500-8000 mg single oral doses; 500-2000 mg b.i.d. and 1000-2000 mg q.d. multiple oral doses (imeglimin hydrochloride)",
    regions = "Japan, Europe, United States",
    notes = "Nine studies used for model development (Table 1): single/multiple ascending dose, Western and Japanese renal impairment, T2DM with renal impairment, Western phase IIa (003, 004), Western and Japanese phase IIb, and Japanese phase III TIMES 1. TIMES 2 (n = 690) and TIMES 3 (n = 106) were used only for external evaluation. Demographics from Table 2 (phase I n = 160, phase II n = 604, TIMES 1 n = 103); sex and Japanese percentages are computed from the Table 2 counts. Western subjects were mostly White (Table 2 note)."
  )

  ini({
    # Structural parameters (Table 3) for the reference patient: Western, optimised
    # tablet, non-fasting, 1000 mg, WT 77.35 kg, eGFR 81.4 mL/min/1.73 m^2, age 59 y.
    lka <- log(0.144); label("Absorption rate constant at the 1000 mg reference dose (1/h)") # Table 3: ka = 0.144 h-1 (RSE 4.42%)
    ltlag <- log(0.229); label("Absorption lag time (h)") # Table 3: ALAG = 0.229 h (RSE 11.0%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability at the 1000 mg reference dose (unitless)") # Table 3: F = 1 (FIX)
    lcl <- log(66.9); label("Apparent clearance CL/F (L/h)") # Table 3: CL/F = 66.9 L/h (RSE 1.70%)
    lvc <- log(142); label("Apparent central volume Vc/F (L)") # Table 3: Vc/F = 142 L (RSE 5.65%)
    lq <- log(15.9); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 3: Q/F = 15.9 L/h (RSE 7.27%)
    lvp <- log(374); label("Apparent peripheral volume Vp/F (L)") # Table 3: Vp/F = 374 L (RSE 6.56%)

    # Dose nonlinearity (labelled imeglimin hydrochloride dose, mg)
    ld50_fdepot <- log(2410); label("Labelled dose at 50% of the maximal decrease in F, D50 (mg)") # Table 3: D50 = 2410 mg (RSE 13.0%); Equation 7
    e_dose_ka <- -0.138; label("Power exponent on (labelled dose / 1000 mg) for ka (unitless)") # Table 3: theta ka,Dose = -0.138 (RSE 20.5%); Equation 1

    # Renal function, body size and age
    e_crcl_cl <- 0.00951; label("Linear slope of eGFR (capped at 120) on CL/F, centred at 81.4 (per mL/min/1.73 m^2)") # Table 3: theta CL/F,eGFR = 0.00951 (RSE 3.60%); Equation 8
    e_wt_cl <- 0.388; label("Power exponent on (WT / 77.35 kg) for CL/F (unitless)") # Table 3: theta CL/F,WT = 0.388 (RSE 13.4%); Equation 9
    e_age_cl <- -0.343; label("Power exponent on (AGE / 59 y) for CL/F (unitless)") # Table 3: theta CL/F,Age = -0.343 (RSE 12.6%); Equation 10
    e_wt_vc <- 0.802; label("Power exponent on (WT / 77.35 kg) for Vc/F (unitless)") # Table 3: theta Vc/F,WT = 0.802 (RSE 22.6%); Equation 1
    e_age_q <- -0.859; label("Power exponent on (AGE / 59 y) for Q/F (unitless)") # Table 3: 'Age on Q/F' = -0.859 (RSE 16.0%); Equation 1
    e_japanese_q <- -0.291; label("Fractional change in Q/F for Japanese subjects (unitless)") # Table 3: 'Japanese on Q/F' = -0.291 (RSE 17.7%); Equation 2

    # Formulation and food (fractional effects, Equation 2)
    e_capsule_tlag <- 2.45; label("Fractional change in lag time for the capsule vs optimised tablet (unitless)") # Table 3: 'Capsule formulation on ALAG' = 2.45 (RSE 16.0%)
    e_convtab_tlag <- 0.719; label("Fractional change in lag time for the conventional vs optimised tablet (unitless)") # Table 3: 'Conventional tablet formulation on ALAG' = 0.719 (RSE 28.6%)
    e_capsconv_ka <- 0.298; label("Fractional change in ka for the capsule or conventional tablet vs optimised tablet (unitless)") # Table 3: 'Capsule or conventional tablet formulation on ka' = 0.298 (RSE 15.4%)
    e_fasted_tlag <- -0.449; label("Fractional change in lag time when fasted or semi-fasted vs non-fasting (unitless)") # Table 3: 'Fasted or semi-fasted on ALAG' = -0.449 (RSE 12.4%)
    e_fasted_ka <- 0.389; label("Fractional change in ka when fasted or semi-fasted vs non-fasting (unitless)") # Table 3: 'Fasted or semi-fasted on ka' = 0.389 (RSE 27.2%)

    # IIV. Table 3 reports each IIV as a 'CV' with RSEs 'on approximate standard
    # deviation scale'; the residual rows in the same column are log-scale SDs, so the
    # IIV entries are read as sqrt(omega) and squared here.
    etalka ~ 0.051984 # Table 3: IIV ka (CV) 0.228 -> 0.228^2
    etalcl + etalfdepot ~ c(0.212521, 0.200216, 0.279841) # Table 3: IIV CL 0.461, IIV F 0.529, correlation CL-F 0.821 -> 0.461^2, 0.821*0.461*0.529, 0.529^2
    etalvc ~ 0.416025 # Table 3: IIV Vc (CV) 0.645 -> 0.645^2
    etalvp ~ 0.366025 # Table 3: IIV Vp (CV) 0.605 -> 0.605^2

    # Residual error: additive on the log scale (Methods), three magnitudes
    expSdPhase23 <- 0.359; label("Log-scale residual SD, post-dose samples in phase II / III studies") # Table 3: 'Proportional RUV (CV)' = 0.359 (RSE 7.05%)
    expSdPhase1 <- 0.190; label("Log-scale residual SD, post-dose samples in phase I studies") # Table 3: 'Proportional phase I RUV (CV)' = 0.190 (RSE 4.08%)
    expSdPredose <- 0.505; label("Log-scale residual SD, pre-dose samples") # Table 3: 'Proportional predose RUV (CV)' = 0.505 (RSE 3.49%)
  })

  model({
    # 1. Dose nonlinearity on the labelled (hydrochloride) dose. amt is free base
    #    (labelled x 0.810, Methods 'Clinical studies'), so the labelled dose is
    #    podo(depot) / 0.810. Before the first dose arrives in the depot (including the
    #    first lag interval) podo(depot) is NA or 0, which would poison the ODE through
    #    (dose / 1000)^e_dose_ka; the depot is empty then, so the neutral reference dose
    #    is substituted.
    dose_label <- podo(depot) / 0.810
    if (is.na(dose_label) || dose_label <= 0) {
      dose_label <- 1000
    }
    d50_fdepot <- exp(ld50_fdepot)
    # Equation 7: inhibitory Emax normalised to F = 1 at 1000 mg
    f_dose <- 1 - (dose_label / (dose_label + d50_fdepot) - 1000 / (1000 + d50_fdepot))

    # 2. Covariate terms
    # Equation 8: linear eGFR effect capped at 120, centred at 81.4
    crcl_cap <- min(CRCL, 120)
    f_crcl_cl <- 1 + (crcl_cap - 81.4) * e_crcl_cl
    fasted <- 1 - FED

    # 3. Individual parameters
    ka <- exp(lka + etalka) * (dose_label / 1000)^e_dose_ka *
      (1 + e_capsconv_ka * (FORM_CAPSULE + FORM_IMEGLIMIN_CONVENTIONAL_TABLET)) *
      (1 + e_fasted_ka * fasted)
    tlag <- exp(ltlag) * (1 + e_capsule_tlag * FORM_CAPSULE) *
      (1 + e_convtab_tlag * FORM_IMEGLIMIN_CONVENTIONAL_TABLET) *
      (1 + e_fasted_tlag * fasted)
    fdepot <- exp(lfdepot + etalfdepot) * f_dose
    cl <- exp(lcl + etalcl) * f_crcl_cl * (WT / 77.35)^e_wt_cl * (AGE / 59)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 77.35)^e_wt_vc
    q <- exp(lq) * (AGE / 59)^e_age_q * (1 + e_japanese_q * RACE_JAPANESE)
    vp <- exp(lvp + etalvp)

    # 4. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 6. Bioavailability and lag time
    f(depot) <- fdepot
    alag(depot) <- tlag

    # 7. Observation (mg/L = ug/mL of imeglimin free base) and residual error.
    #    Pre-dose samples take the pre-dose magnitude; otherwise phase I vs phase II/III.
    Cc <- central / vc
    sd_ruv <- expSdPredose * SAMPLE_PREDOSE +
      (1 - SAMPLE_PREDOSE) * (expSdPhase1 * STUDY_IMEGLIMIN_PHASE1 + expSdPhase23 * (1 - STUDY_IMEGLIMIN_PHASE1))
    Cc ~ lnorm(sd_ruv)
  })
}
