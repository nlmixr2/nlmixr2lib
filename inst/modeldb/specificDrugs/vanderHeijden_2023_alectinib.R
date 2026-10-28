vanderHeijden_2023_alectinib <- function() {
  description <- paste(
    "One-compartment population PK model for oral alectinib (an ALK",
    "inhibitor) in adults with ALK-positive non-small cell lung cancer,",
    "used by van der Heijden 2023 to simulate a stable-isotope-labelled",
    "microtracer (alectinib-d6) food-effect study during steady-state",
    "alectinib treatment. The structural model and its estimates are the",
    "European Medicines Agency Alecensa assessment-report model (138",
    "patients, 600 mg twice daily with a meal): sequential zero-order input",
    "into the depot over D1 = 3.44 h followed by first-order absorption at",
    "ka = 1.30 1/h, CL/F = 81.9 L/h and V/F = 4016 L, each allometrically",
    "scaled by body weight (exponents 0.75 and 1, reference 70 kg).",
    "Relative bioavailability is 1 for a dose taken with food and 0.714",
    "fasted (the assumed 40% food effect), with an assumed 8% (CV)",
    "between-occasion variability on bioavailability. Correlated",
    "between-subject variability on CL/F and V/F, independent variability",
    "on ka and D1, and a proportional residual error of 19%; the",
    "assessment-report model's additive error is dropped for the microtracer",
    "concentration range.",
    sep = " "
  )
  reference <- paste(
    "van der Heijden LT, Steeghs N, Beijnen JH, Huitema ADR, Dorlo TPC.",
    "The use of microtracers in food-effect trials: An alternative study",
    "design for toxic drugs with long half-lives exemplified by the case for",
    "alectinib. Clin Transl Sci. 2023;16:2557-2564.",
    "doi:10.1111/cts.13647. Structural model and estimates from the European",
    "Medicines Agency (EMA) Alecensa Assessment Report (2016), as reproduced",
    "in the paper's Data S1 (Table S2 and NONMEM code).",
    sep = " "
  )
  vignette <- "vanderHeijden_2023_alectinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "alectinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alectinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1)",
        "with reference weight 70 kg (Data S1 NONMEM code: TVCL =",
        "THETA(1)*(WT/70)**0.75, TVV = THETA(2)*(WT/70)**1). The trial",
        "simulations resampled body weight from a previous alectinib NSCLC",
        "cohort (median 77.5 kg, range 49-123 kg; Methods 'Clinical trial",
        "simulation')."
      ),
      source_name = "WT"
    ),
    FED = list(
      description = "Fed-state dose-record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Per-dose-record indicator: 1 = dose taken with a (standardized)",
        "meal, 0 = fasted. Selects between the two relative-bioavailability",
        "thetas of the Data S1 NONMEM code: F1 = THETA(5) = 1 FIX with food",
        "and, when FLAG = 1, F1 = THETA(6) = 0.714 under fasting conditions",
        "(Table S2). The paper's FLAG column is therefore FED = 1 - FLAG.",
        "The reference state of the model (F = 1) is the fed state, which",
        "is how alectinib was dosed in the source population; FED = 0 is",
        "this column's reference category only by register convention.",
        "In the proposed design the fed microtracer dose is study day 1",
        "(OCC = 1) and the fasted one study day 9 (OCC = 2)."
      ),
      source_name = "FLAG (inverted)"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on relative bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "1 = the first microtracer administration (fed, study day 1); 2 =",
        "the second (fasted, study day 9). Data S1 NONMEM code: REPLACE",
        "ETA(OCC_F1)=ETA(5,6) with $OMEGA BLOCK(1) 0.0064 and $OMEGA",
        "BLOCK(1) SAME, i.e. one eta per occasion with a shared variance.",
        "Records with OCC outside 1:2 (e.g. OCC = 0 for routine",
        "therapeutic doses) carry no between-occasion variability."
      ),
      source_name = "OCC_F1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 138L,
    n_studies = 2L,
    age_range = "22-79 years",
    age_median = "52 years",
    sex_female_pct = 56,
    race_ethnicity = c(
      White = 67,
      Asian = 26,
      Other = 3,
      `Black or African American` = 0.7,
      `American Indian or Alaska Native` = 0.7,
      Unknown = 2
    ),
    disease_state = paste(
      "Locally advanced or metastatic ALK-positive non-small cell lung",
      "cancer after crizotinib; CNS metastases 61%, adenocarcinoma 95%,",
      "prior chemotherapy 80%, ECOG PS 0/1/2 32%/59%/9% (Table S1)."
    ),
    dose_range = "Alectinib 600 mg orally twice daily with a meal.",
    regions = "Global (pivotal phase II studies NP28673 and NP28761).",
    notes = paste(
      "Demographics of the model-development dataset are from Data S1",
      "Table S1 (138 patients from the pooled phase II studies NP28673 and",
      "NP28761, reference 14). Weight was not tabulated for this cohort.",
      "The paper itself fitted no model to data: it re-used the EMA",
      "assessment-report population PK model (reference 13) to simulate a",
      "single 100 ug alectinib-d6 microtracer given with food and, 8 days",
      "later, fasted, in patients at steady state on 600 mg twice daily,",
      "assuming the microtracer and alectinib share the same PK. The",
      "food-effect bioavailability (0.714) and the between-occasion",
      "variability (0.0064) were added for the simulation (Table S2",
      "footnote)."
    )
  )

  ini({
    # Structural PK -- Data S1 Table S2 'Fixed effects' (EMA assessment-report
    # estimates), identical to the $THETA initial values of the Data S1 code.
    # Reference weight 70 kg.
    lcl <- log(81.9) ; label("Apparent clearance CL/F at 70 kg (L/h)") # Table S2 'Cl/F (L/h)' = 81.9; $THETA(1)
    lvc <- log(4016) ; label("Apparent volume of distribution V/F at 70 kg (L)") # Table S2 'V/F (L)' = 4016; $THETA(2)
    lka <- log(1.30) ; label("First-order absorption rate constant ka (1/h)") # Table S2 'Ka (h-1)' = 1.30; $THETA(3)
    ld1 <- log(3.44) ; label("Duration of the zero-order input into the depot D1 (h)") # Table S2 'D1 (h)' = 3.44; $THETA(4)

    # Allometric exponents, fixed in the Data S1 code (Methods: 'power
    # coefficients of 0.75 for clearance and 1.0 for volume of distribution').
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL/F (unitless)") # Data S1 code TVCL = THETA(1)*(WT/70)**0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of body weight on V/F (unitless)") # Data S1 code TVV = THETA(2)*(WT/70)**1

    # Relative bioavailability. The two values are separate thetas in the Data
    # S1 code, selected by FLAG (FLAG = 1 is the fasted occasion). Both are
    # fixed: 1 for the fed reference, 0.714 the assumed 40% food effect
    # (Equation 1: FE = (1 - 0.714)/0.714 = 40%).
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability with food (fraction)") # Table S2 'F (fed conditions)' = 1; $THETA(5) 1 FIX
    lfdepot_fasted <- fixed(log(0.714)) ; label("Relative bioavailability fasted (fraction)") # Table S2 'F (fasted conditions)' = 0.714 (added for the simulation); Methods 'fixed to 0.714 for microtracer administration under fasting conditions'

    # Between-subject variability (Table S2; variances on the log scale, all
    # parameters EXP(ETA)). The CL/F-V/F covariance is 0.0787 in Table S2 and
    # 0.07865012 in the $OMEGA BLOCK(2) of the Data S1 code; the unrounded
    # code value is used.
    etalcl + etalvc ~ c(0.162, 0.07865012, 0.161) # Table S2 'IIV Cl/F' = 0.162, 'Covariance Cl/F ~ Vd/F' = 0.0787, 'IIV Vd/F' = 0.161
    etalka ~ 0.398 # Table S2 'IIV Ka' = 0.398
    etald1 ~ 0.178 # Table S2 'IIV D1' = 0.178

    # Between-occasion variability on bioavailability, one eta per occasion with
    # a shared variance (Data S1 $OMEGA BLOCK(1) 0.0064 and BLOCK(1) SAME).
    # The value is an assumption of the simulation (Methods: 'a relatively
    # small inter-occasion variability of 8.0% was assumed'; 0.08^2 = 0.0064),
    # not an estimate, so both are fixed.
    etaiov_fdepot_1 ~ fixed(0.0064) # Table S2 'IOV' = 0.0064 (added for the simulation)
    etaiov_fdepot_2 ~ fixed(0.0064) # $OMEGA BLOCK(1) SAME

    # Residual error. Data S1 $ERROR Y = F*(1+EPS(1))+EPS(2) with $SIGMA 0.0361
    # (proportional variance) and 0 FIX (additive); the additive term was
    # excluded for the microtracer concentration range (Methods), so only the
    # proportional term is carried.
    propSd <- 0.19 ; label("Proportional residual error (fraction)") # Table S2 'Proportional error' = 0.0361 (variance); sqrt(0.0361) = 0.19, the 'proportional error of 19%' of the Figure 3 legend
  })

  model({
    # Occasion indicators; OCC outside 1:2 removes the between-occasion term.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    # F1 = THETA(5)*EXP(ETA(OCC_F1)) fed; THETA(6)*EXP(ETA(OCC_F1)) fasted.
    fdepot <- exp((lfdepot * FED + lfdepot_fasted * (1 - FED)) + iov_fdepot)

    kel <- cl / vc

    # Sequential zero-order then first-order absorption (ADVAN2 with a
    # modelled duration D1 on the depot). Dose records need rate = -2 so the
    # dose enters the depot over dur(depot).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    dur(depot) <- d1
    f(depot) <- fdepot

    # Dose in mg and V/F in L give mg/L; x1000 gives ng/mL (= ug/L, the unit
    # of the Data S1 additive-error comment).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
