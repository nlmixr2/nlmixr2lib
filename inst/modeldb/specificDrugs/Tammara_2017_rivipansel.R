Tammara_2017_rivipansel <- function() {
  description <- paste(
    "Three-compartment IV population PK model for rivipansel in adults",
    "and adolescents with sickle cell disease (SCD) and in healthy adult",
    "volunteers (Tammara 2017). Rivipansel is a pan-selectin antagonist",
    "given as a 20-minute IV infusion; renal excretion of unchanged drug",
    "is the primary clearance mechanism. The integrated population PK",
    "model pools 109 subjects across three phase I studies (rivipansel",
    "studies 101, 102, 103) and one phase II SCD study (NCT01119833,",
    "Telen 2015). Clearance is a power function of creatinine clearance",
    "(CRCL, raw Cockcroft-Gault mL/min reference 150) with an additive",
    "23.4% shift in the phase II SCD cohort (STUDY_RIV201) attributed to",
    "glomerular hyperfiltration. The central, first peripheral, and",
    "second peripheral volumes share a single estimated body-weight",
    "exponent (0.569, reference 70 kg). The additive and proportional",
    "residual error magnitudes differ between the phase I and phase II",
    "cohorts and are selected per observation via STUDY_RIV201."
  )
  reference <- paste(
    "Tammara BK, Harnisch LO.",
    "Dose Selection Based on Modeling and Simulation for Rivipansel in",
    "Pediatric Patients Aged 6 to 11 Years With Sickle Cell Disease.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(12):845-854.",
    "doi:10.1002/psp4.12263."
  )
  vignette <- "Tammara_2017_rivipansel"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed body weight in kg. Shared power-form covariate on",
        "central volume V1, first peripheral V2, and second peripheral",
        "V3, all normalised to a reference of 70 kg",
        "(Table 1 footnote c: V1,2,3 = theta_V1,2,3 * (WT/70)^0.569)."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance in raw (NOT BSA-normalised) mL/min.",
        "Estimated in adults by the Cockcroft-Gault formula and in",
        "children by the Schwartz formula with body-surface-area",
        "adjustment so that both yield absolute mL/min on a common",
        "scale across the age range."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Reference value 150 mL/min",
        "(Table 1 footnote b: CL = 1.25 * (CRCL/150)^0.468 * (1 + 0.234 * STUD)).",
        "NOT BSA-normalised: the per-model unit is raw mL/min, matching",
        "Delattre 2010 amikacin rather than the BSA-normalised CRCL",
        "default of the canonical covariate. Document the unit",
        "explicitly in covariateData[[CRCL]]$units so downstream users",
        "do not silently feed mL/min/1.73 m^2 into the power form."
      ),
      source_name        = "CRCL"
    ),
    STUDY_RIV201 = list(
      description        = paste(
        "Phase II rivipansel SCD study (NCT01119833) indicator.",
        "1 = subject is from the phase II study in SCD patients with",
        "vaso-occlusive crisis; 0 = subject is from one of the three",
        "phase I studies (studies 101, 102, 103) pooled into the",
        "integrated population PK analysis."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rivipansel phase I studies 101, 102, 103)",
      notes              = paste(
        "Subject-level / time-fixed. Used both as an additive shift on",
        "typical clearance (1 + 0.234 * STUDY_RIV201) -- interpreted by",
        "the authors as the SCD-hyperfiltration component -- and to",
        "select between phase I and study 201 cohort-specific residual",
        "error magnitudes in model() (additive 0.29 vs 0.165 ug/mL,",
        "proportional 9.52% vs 23.2% CV). For simulation use cases",
        "targeting the SCD population (the paper's stated goal), set",
        "STUDY_RIV201 = 1 for every subject."
      ),
      source_name        = "STUD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 109,
    n_studies      = 4,
    age_range      = "12-51 years (pooled phase I + phase II)",
    age_median     = "not reported in main text",
    weight_range   = "not reported in main text",
    weight_median  = "70 kg (reference value used in the allometric covariate form, Table 1 footnote c)",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Pooled cohort across three rivipansel phase I studies (studies",
      "101 and 102 in healthy adult volunteers, study 103 NCT00911495",
      "in adults with SCD not in vaso-occlusive crisis) and the",
      "rivipansel phase II study (NCT01119833; Telen 2015 Blood",
      "125:2656-2664) in adolescents and adults with SCD hospitalized",
      "for vaso-occlusive crisis. 12 children aged 12-17 years were",
      "enrolled in the phase II study."
    ),
    dose_range     = paste(
      "Rivipansel 2-40 mg/kg IV (20-minute infusion). Phase I single",
      "dose 2-40 mg/kg (study 101); phase I multiple dose 2-20 mg/kg",
      "q8h with a 40 mg/kg loading dose followed by 20 mg/kg q8h",
      "(study 102); phase I in SCD 20 mg/kg loading + 10 mg/kg at 10 h",
      "(study 103). Phase II 20 mg/kg loading + 10 mg/kg q12h x 7 days",
      "(low dose) or 40 mg/kg loading + 20 mg/kg q12h x 7 days (high",
      "dose). Renal clearance is the primary elimination mechanism",
      "(> 90% of dose recovered unchanged in urine)."
    ),
    regions        = "United States (Pfizer-sponsored development program)",
    notes          = paste(
      "Demographics summarized in Tammara 2017 Methods / Study design",
      "and data sources. NONMEM 7.3.0 with FOCE-INTER; precision of",
      "parameter estimates obtained via Sampling Importance Resampling",
      "(SIR; Dosne 2016). The structural three-compartmental model and",
      "variance components were carried from the earlier per-study",
      "analyses (unpublished Pfizer reports) and all parameters were",
      "reestimated against the integrated dataset; the paper's stated",
      "intent is dose extrapolation to pediatric SCD patients aged",
      "6-11 years for the phase III study (NCT02187003), not de novo",
      "popPK model building."
    )
  )

  ini({
    # ============================================================
    # Structural parameters -- Tammara 2017 Table 1 "Structural
    # model" rows. Reference CRCL = 150 mL/min; reference WT = 70 kg
    # (Table 1 footnotes b and c).
    # ============================================================
    lcl  <- log(1.25)
    label("Clearance at CRCL = 150 mL/min, phase I cohort (L/h)")           # Table 1: CL = 1.25 L/h (SIR %RSE 2.12, 95% CI 1.2-1.3)
    lvc  <- log(6.24)
    label("Central volume V1 at WT = 70 kg (L)")                            # Table 1: V1 = 6.24 L (SIR %RSE 2.83, 95% CI 5.91-6.6)
    lq   <- log(2.62)
    label("Rapid inter-compartmental clearance Q_rapid (L/h)")              # Table 1: Q_rapid = 2.62 L/h (SIR %RSE 5.86, 95% CI 2.34-2.96)
    lvp  <- log(4.02)
    label("First peripheral volume V2 at WT = 70 kg (L)")                   # Table 1: V2 = 4.02 L (SIR %RSE 2.72, 95% CI 3.8-4.24)
    lq2  <- log(0.0316)
    label("Slow inter-compartmental clearance Q_slow (L/h)")                # Table 1: Q_slow = 0.0316 L/h (SIR %RSE 17.5, 95% CI 0.0226-0.0441)
    lvp2 <- log(0.656)
    label("Second peripheral volume V3 at WT = 70 kg (L)")                  # Table 1: V3 = 0.656 L (SIR %RSE 7.85, 95% CI 0.563-0.76)

    # ============================================================
    # Covariate effects -- Tammara 2017 Table 1 "Covariate model"
    # rows.
    # ============================================================
    e_study_riv201_cl <- 0.234
    label("Additive phase II SCD study effect on CL (fraction; 1 + e * STUDY_RIV201)")  # Table 1: study effect on CL = 23.4% (SIR %RSE 21.2, 95% CI 14.2-33.9); paper interprets as SCD-hyperfiltration component
    e_crcl_cl         <- 0.468
    label("Power exponent for CRCL on CL (unitless)")                       # Table 1: exponent for CRCL on CL = 0.468 (SIR %RSE 15.9, 95% CI 0.325-0.619)
    e_wt_vc_vp_vp2    <- 0.569
    label("Body-weight exponent shared across V1, V2, V3 (unitless)")       # Table 1: exponent for WT on V1, V2, V3 = 0.569 (SIR %RSE 15.1, 95% CI 0.406-0.740)

    # ============================================================
    # Inter-individual variability -- Tammara 2017 Table 1 column
    # "%CV" with omega^2 = log(1 + CV^2). No correlation structure
    # is reported; all six etas are independent.
    # ============================================================
    etalcl  ~ log(1 + 0.188^2)   # Table 1 IIV: CL %CV = 18.8 (SIR %RSE 15.4, 95% CI 16.3-21.8). g-shrinkage 4.5%
    etalvc  ~ log(1 + 0.241^2)   # Table 1 IIV: V1 %CV = 24.1 (SIR %RSE 17.8, 95% CI 20.6-28.9)
    etalq   ~ log(1 + 0.198^2)   # Table 1 IIV: Q_rapid %CV = 19.8 (SIR %RSE 65.7, 95% CI 5.84-32.4). g-shrinkage 79%
    etalvp  ~ log(1 + 0.155^2)   # Table 1 IIV: V2 %CV = 15.5 (SIR %RSE 26.8, 95% CI 11.9-19.7)
    etalq2  ~ log(1 + 0.155^2)   # Table 1 IIV: Q_slow %CV = 15.5 (SIR %RSE 191, 95% CI 1.47-41.4)
    etalvp2 ~ log(1 + 0.153^2)   # Table 1 IIV: V3 %CV = 15.3 (SIR %RSE 63.5, 95% CI 5.61-24.4)

    # ============================================================
    # Residual error -- Tammara 2017 Table 1 reports separate
    # additive (ug/mL) and proportional (%CV) magnitudes for the
    # phase I studies versus study 201 (the phase II SCD study).
    # The cohort-specific magnitudes below are combined in model()
    # via the STUDY_RIV201 indicator into the canonical propSd and
    # addSd. e-shrinkage 7.8% for phase I, 11.9% for study 201.
    # ============================================================
    addSd_phase1   <- 0.29
    label("Additive residual SD, phase I cohort (ug/mL)")                   # Table 1: additive residual, phase I = 0.29 ug/mL (SIR %RSE 26.7, 95% CI 0.172-0.479)
    addSd_study201 <- 0.165
    label("Additive residual SD, study 201 / phase II SCD cohort (ug/mL)")  # Table 1: additive residual, study 201 = 0.165 ug/mL (SIR %RSE 70.7, 95% CI 0.0144-0.457)
    propSd_phase1   <- 0.0952
    label("Proportional residual SD, phase I cohort (fraction)")            # Table 1: proportional residual, phase I = 9.52% CV (SIR %RSE 4.33, 95% CI 9.14-9.95)
    propSd_study201 <- 0.232
    label("Proportional residual SD, study 201 / phase II SCD cohort (fraction)")  # Table 1: proportional residual, study 201 = 23.2% CV (SIR %RSE 12, 95% CI 20.9-26.2)
  })

  model({
    # ------------------------------------------------------------
    # Individual structural parameters. CL is a power function of
    # CRCL with an additive phase II SCD shift; the three volumes
    # share a single body-weight exponent.
    # ------------------------------------------------------------
    cl  <- exp(lcl  + etalcl)  * (CRCL / 150)^e_crcl_cl *
           (1 + e_study_riv201_cl * STUDY_RIV201)
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc_vp_vp2
    q   <- exp(lq   + etalq)
    vp  <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vc_vp_vp2
    q2  <- exp(lq2  + etalq2)
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc_vp_vp2

    # ------------------------------------------------------------
    # Three-compartment IV disposition. Dose enters central (IV
    # infusion); two peripheral compartments exchange with central
    # via the rapid (Q) and slow (Q2) inter-compartmental
    # clearances. Concentration in central is in mg/L = ug/mL with
    # mg dosing and L volumes.
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central)     <-  k21 * peripheral1 + k31 * peripheral2 -
                          (kel + k12 + k13) * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    Cc <- central / vc

    # ------------------------------------------------------------
    # Cohort-conditional residual error: STUDY_RIV201 = 1 selects
    # the phase II SCD magnitudes (Table 1 "study 201" rows);
    # STUDY_RIV201 = 0 selects the phase I magnitudes. Same idiom
    # as deHoogd_2017_morphine.
    # ------------------------------------------------------------
    addSd  <-  addSd_phase1  * (1 - STUDY_RIV201) +  addSd_study201 * STUDY_RIV201
    propSd <- propSd_phase1  * (1 - STUDY_RIV201) + propSd_study201 * STUDY_RIV201

    Cc ~ add(addSd) + prop(propSd)
  })
}
