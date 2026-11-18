Goeyvaerts_2026_canagliflozin <- function() {
  description <- paste(
    "Population PK and PK/PD model for oral canagliflozin and HbA1c in",
    "paediatric patients aged 10-17 years with type 2 diabetes mellitus.",
    "Two-compartment PK parameterised in rate constants, with a lag time",
    "followed by zero-order input into the depot and first-order absorption;",
    "body weight and sex on Vc/F, total daily dose on ke and",
    "overencapsulation on the lag time, all fixed to the adult model; log-scale",
    "residual SD re-estimated separately for rich (phase 1) and sparse",
    "(phase 3) sampling. HbA1c turnover with an additive placebo effect and a",
    "canagliflozin Emax effect scaled by (baseline - 5)/3; background insulin",
    "and randomised arm on baseline HbA1c, background insulin on the placebo",
    "effect."
  )
  reference <- paste(
    "Goeyvaerts N, Gisleskog PO, Neyens M, Perez Ruixo JJ, Ali SR, Cirillo I,",
    "Ackaert O.",
    "Population pharmacokinetics and pharmacodynamics of canagliflozin in",
    "paediatric patients with type 2 diabetes mellitus.",
    "Br J Clin Pharmacol. 2026;92(10):3622-3633. doi:10.1002/bcp.70644."
  )
  vignette <- "Goeyvaerts_2026_canagliflozin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (canagliflozin, Cc); % NGSP (HbA1c, Hba1c)"
  )

  # depot: the canagliflozin dose, which must carry rate = -2 so that the
  # modelled zero-order input duration dur(depot) applies.
  dosing <- "depot"

  # Two observation outputs with their own log-scale residual SDs.
  paper_specific_residual_sds <- c("expSdPhase1", "expSdPhase3")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters Vc/F as (WT/75.8)^0.583 (Table 2 footnote b); because CL/F = Vc/F * ke and ke carries no weight term, CL/F scales with the same exponent (Methods 'popPK and PK/PD models').",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex, female indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the typical Vc/F at 75.8 kg: 99.3 L for males, 82.6 L for females (Table 2).",
      source_name = "Sex"
    ),
    DOSE_CANA_MGD = list(
      description = "Total daily canagliflozin dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters ke as (TDD/100)^-0.0631 (Table 2 footnote c). 100 mg for",
        "100 mg QD or 50 mg BID; 300 mg for 300 mg QD or 150 mg BID. For",
        "placebo subjects (no canagliflozin given) any value may be supplied;",
        "the model substitutes the 100 mg reference when the value is 0 so that",
        "ke stays finite."
      ),
      source_name = "TDD"
    ),
    FORM_CANAGLIFLOZIN_OVERENCAP = list(
      description = "Overencapsulated canagliflozin tablet indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nonencapsulated tablet)",
      notes = "Selects the lag time: 0.147 h for nonencapsulated, 0.262 h for overencapsulated tablets (Table 2). Every paediatric patient received overencapsulated tablets, so set 1 to reproduce the paediatric analysis (Table 2 footnote a).",
      source_name = "overencapsulation"
    ),
    STUDY_PHASE3 = list(
      description = "Sparse-sampling (phase 2/3) record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1, rich sampling)",
      notes = paste(
        "Selects the canagliflozin log-scale residual SD only: 0.227 for the",
        "phase 1 rich data, 1.14 for the 'Phase 2 and 3' sparse data (Table 2;",
        "the paediatric data set held only phase 3 sparse data in that stratum).",
        "Does not affect the typical-value prediction; set 0 to simulate the",
        "intensive-sampling residual error."
      ),
      source_name = "Study phase"
    ),
    CONMED_INSULIN = list(
      description = "Background insulin therapy (alone or with metformin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no background insulin)",
      notes = "Selects the baseline HbA1c (Table 4, 'Insulin' versus 'No insulin' rows) and the placebo effect Efp (0.171 % with insulin, -0.238 % without). 64 of 160 PK/PD patients (40%) were on insulin or metformin plus insulin (Table S1).",
      source_name = "background insulin"
    ),
    TRT_CANAGLIFLOZIN = list(
      description = "Randomised canagliflozin treatment-arm indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (randomised to placebo)",
      notes = "Selects the baseline HbA1c only (Table 4, 'canagliflozin' versus 'placebo' rows; the arms were imbalanced at baseline, median 7.50 % versus 8.00 %, Results 'PK/PD analysis'). The drug effect itself is concentration-driven and vanishes when no canagliflozin is dosed.",
      source_name = "treatment group"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Age on Vc/F in the adult model (Hoeben 2016) was removed during the paediatric model update (Results 'PopPK analysis')."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "eGFR on ke in the adult model was removed during the paediatric model update (Results 'PopPK analysis')."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "BMI on ka and on the lag time in the adult model were removed during the paediatric model update (Results 'PopPK analysis')."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "canagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "canagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    Hba1c = list(analyte = "HbA1c", units = "%", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 160L,
    n_studies = 2L,
    n_observations = "330 canagliflozin concentrations (90 patients); 702 HbA1c values (160 patients)",
    studies = "Phase 1 open-label multiple-dose study (n = 17, 14 days) and phase 3 randomised placebo-controlled study (n = 73 canagliflozin, 87 placebo; 52 weeks).",
    age_range = "10-17 years",
    age_median = "15 years",
    weight_range = "44.8-170 kg (PK analysis); 40.6-162 kg (PK/PD analysis)",
    weight_median = "79.5 kg (PK analysis); 75.2 kg (PK/PD analysis)",
    sex_female_pct = 67.8,
    race_ethnicity = c(
      White_not_Hispanic = 16.7,
      Black = 20.0,
      White_Hispanic = 22.2,
      Asian = 33.3,
      American_Indian_Alaska_Native = 4.4,
      Other = 3.3
    ),
    disease_state = "Type 2 diabetes mellitus (paediatric, 10-17 years)",
    dose_range = "Canagliflozin 100 or 300 mg orally once daily (phase 3: 100 mg, with uptitration to 300 mg at week 13 for 12 patients) or placebo.",
    co_medication = "PK/PD cohort: metformin 46.2%, metformin + insulin 28.7%, insulin 11.2%, diet and exercise only 13.8% (Table S1).",
    renal_function = "Median eGFR 160 mL/min/1.73 m^2 (range 70.3-277; Schwartz formula).",
    regions = "Multinational",
    notes = "Demographics of the 90 PK patients from Table 1 and of the 160 PK/PD patients from Table S1. All PK parameters except the residual SDs, and the HbA1c half-life, Emax, EC50 and IIV, are fixed to the adult canagliflozin models (Hoeben 2016 Clin Pharmacokinet 55:209; de Winter 2017 Br J Clin Pharmacol 83:1072)."
  )

  ini({
    # PK: Table 2. Every structural, covariate and IIV value is FIX to the
    # adult popPK model (Hoeben 2016); only the residual SDs were estimated.
    lvc_male <- fixed(log(99.3)); label("Apparent central volume Vc/F for males at 75.8 kg (L)") # Table 2 'Vc/F Males' = 99.3 FIX
    lvc_female <- fixed(log(82.6)); label("Apparent central volume Vc/F for females at 75.8 kg (L)") # Table 2 'Vc/F Females' = 82.6 FIX
    lkel <- fixed(log(0.145)); label("Elimination rate constant ke at 100 mg/day (1/h)") # Table 2 'ke' = 0.145 FIX
    lka <- fixed(log(3.68)); label("First-order absorption rate constant ka (1/h)") # Table 2 'ka' = 3.68 FIX
    lk12 <- fixed(log(0.101)); label("Distribution rate constant central to peripheral k23 (1/h)") # Table 2 'k23' = 0.101 FIX
    lk21 <- fixed(log(0.0856)); label("Distribution rate constant peripheral to central k32 (1/h)") # Table 2 'k32' = 0.0856 FIX
    ltlag_nonencap <- fixed(log(0.147)); label("Lag time for nonencapsulated tablets (h)") # Table 2 'Tlag Nonencapsulated tablets' = 0.147 FIX
    ltlag_overencap <- fixed(log(0.262)); label("Lag time for overencapsulated tablets (h)") # Table 2 'Tlag Overencapsulated tablets' = 0.262 FIX
    ld1 <- fixed(log(0.604)); label("Duration of zero-order input into the depot D1 (h)") # Table 2 'D1' = 0.604 FIX

    e_wt_vc <- fixed(0.583); label("Power exponent of (WT/75.8) on Vc/F (unitless)") # Table 2 'Body weight on Vc/F' = 0.583 FIX; footnote b reference 75.8 kg
    e_dose_cana_mgd_kel <- fixed(-0.0631); label("Power exponent of (TDD/100) on ke (unitless)") # Table 2 'Total daily dose on ke' = -0.0631 FIX; footnote c reference 100 mg

    # IIV: Table 2 CV% (FIX), converted with omega^2 = log(1 + CV^2).
    etalvc ~ fixed(0.021092) # Table 2 Vc/F IIV 14.6 CV% FIX -> log(1 + 0.146^2)
    etalkel ~ fixed(0.051986) # Table 2 ke IIV 23.1 CV% FIX -> log(1 + 0.231^2)
    etalka ~ fixed(1.5117) # Table 2 ka IIV 188 CV% FIX -> log(1 + 1.88^2)
    etalk21 ~ fixed(0.12314) # Table 2 k32 IIV 36.2 CV% FIX -> log(1 + 0.362^2)
    etaltlag ~ fixed(0.62919) # Table 2 Tlag IIV 93.6 CV% FIX -> log(1 + 0.936^2)

    expSdPhase1 <- 0.227; label("Log-scale residual SD, phase 1 rich PK data (unitless)") # Table 2 'Residual error (SD) Phase 1 paediatric' = 0.227 (RSE 13.4%)
    expSdPhase3 <- 1.14; label("Log-scale residual SD, phase 2/3 sparse PK data (unitless)") # Table 2 'Residual error (SD) Phase 2 and 3 paediatric' = 1.14 (RSE 9.12%)

    # PK/PD: Table 4. Baseline, placebo effect and residual SD estimated on the
    # paediatric data; half-life, Emax, EC50 and IIV FIX to the adult model
    # (de Winter 2017).
    lrbase_ins_pbo <- log(8.50); label("Baseline HbA1c, background insulin, placebo arm (%)") # Table 4 'Baseline HbA1c Insulin placebo' = 8.50 (SE 0.213)
    lrbase_ins_cana <- log(7.82); label("Baseline HbA1c, background insulin, canagliflozin arm (%)") # Table 4 'Baseline HbA1c Insulin canagliflozin' = 7.82 (SE 0.270)
    lrbase_noins_pbo <- log(7.91); label("Baseline HbA1c, no background insulin, placebo arm (%)") # Table 4 'Baseline HbA1c No insulin placebo' = 7.91 (SE 0.190)
    lrbase_noins_cana <- log(7.22); label("Baseline HbA1c, no background insulin, canagliflozin arm (%)") # Table 4 'Baseline HbA1c No insulin canagliflozin' = 7.22 (SE 0.154)
    lkout <- fixed(log(log(2) / 28.2)); label("HbA1c first-order loss rate constant kout (1/day; divided by 24 in model())") # Table 4 'HbA1c half-life' = 28.2 days FIX; footnote a t1/2 = log(2)/kout
    pbo_ins <- 0.171; label("Steady-state placebo effect Efp with background insulin at baseline HbA1c 8 % (% HbA1c)") # Table 4 'Efp insulin' = 0.171 (SE 0.148)
    pbo_noins <- -0.238; label("Steady-state placebo effect Efp without background insulin at baseline HbA1c 8 % (% HbA1c)") # Table 4 'Efp no insulin' = -0.238 (SE 0.107)
    emax <- fixed(-0.738); label("Maximum canagliflozin effect Emax at baseline HbA1c 8 % (% HbA1c)") # Table 4 'Emax' = -0.738 FIX
    lec50 <- fixed(4.12); label("Log of the canagliflozin EC50 (log ng/mL)") # Table 4 'Log EC50' = 4.12 FIX; footnote b EC50 = 61.6 ng/mL

    etalrbase ~ fixed(0.011) # Table 4 'IIV of baseline HbA1c (variance)' = 0.011 FIX (lognormal)
    etapbo ~ fixed(0.369) # Table 4 'IIV of Efp (variance)' = 0.369 FIX (normal, additive on Efp)

    expSd_Hba1c <- 0.0773; label("Log-scale residual SD for HbA1c (unitless)") # Table 4 'Residual error (SD)' = 0.0773 (SE 0.00522)
  })
  model({
    # PK (Methods 'popPK and PK/PD models'; Table 2 and footnotes)
    vc_typ <- exp(lvc_male) * (1 - SEXF) + exp(lvc_female) * SEXF
    vc <- vc_typ * exp(etalvc) * (WT / 75.8)^e_wt_vc

    # Placebo subjects carry no drug; the 100 mg reference keeps ke finite.
    tdd <- DOSE_CANA_MGD
    if (DOSE_CANA_MGD <= 0) tdd <- 100
    kel <- exp(lkel + etalkel) * (tdd / 100)^e_dose_cana_mgd_kel
    ka <- exp(lka + etalka)
    k12 <- exp(lk12)
    k21 <- exp(lk21 + etalk21)
    tlag_typ <- exp(ltlag_nonencap) * (1 - FORM_CANAGLIFLOZIN_OVERENCAP) +
      exp(ltlag_overencap) * FORM_CANAGLIFLOZIN_OVERENCAP
    tlag <- tlag_typ * exp(etaltlag)
    d1 <- exp(ld1)
    cl <- kel * vc

    # PK/PD (Supplement Equations 1-3; Table 4)
    rbase_typ <- exp(lrbase_ins_pbo) * CONMED_INSULIN * (1 - TRT_CANAGLIFLOZIN) +
      exp(lrbase_ins_cana) * CONMED_INSULIN * TRT_CANAGLIFLOZIN +
      exp(lrbase_noins_pbo) * (1 - CONMED_INSULIN) * (1 - TRT_CANAGLIFLOZIN) +
      exp(lrbase_noins_cana) * (1 - CONMED_INSULIN) * TRT_CANAGLIFLOZIN
    rbase <- rbase_typ * exp(etalrbase)
    pbo <- pbo_ins * CONMED_INSULIN + pbo_noins * (1 - CONMED_INSULIN) + etapbo
    kout <- exp(lkout) / 24
    kin <- kout * rbase
    ec50 <- exp(lec50)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Lag, then zero-order input into the depot over D1 (dose rate = -2).
    alag(depot) <- tlag
    dur(depot) <- d1

    # mg / L * 1000 = ng/mL
    Cc <- central / vc * 1000

    # Equation 3: canagliflozin effect; Equation 2: combined effect scaled by
    # (H(0) - 5)/(8 - 5); Equation 1: turnover with kin = kout * H(0).
    efc <- emax * Cc / (Cc + ec50)
    ef <- kout * (efc + pbo) * (rbase - 5) / (8 - 5)
    d / dt(Hba1c) <- ef + kin - kout * Hba1c
    Hba1c(0) <- rbase

    expSdCc <- expSdPhase1 * (1 - STUDY_PHASE3) + expSdPhase3 * STUDY_PHASE3
    Cc ~ lnorm(expSdCc)
    Hba1c ~ lnorm(expSd_Hba1c)
  })
}
