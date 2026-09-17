Lee_2026_oxacillin <- function() {
  description <- paste0(
    "One-compartment intravenous population pharmacokinetic model with first-order ",
    "elimination for oxacillin in 22 preterm and term neonates and infants aged 4 to 82 ",
    "days postnatal (gestational age 23.9 to 40.3 weeks), each given a 25 mg/kg loading ",
    "dose over 30 min followed immediately by a continuous infusion of 160 mg/kg/day (120 ",
    "mg/kg/day for the single infant born at <32 weeks gestation and enrolled at <14 days ",
    "of life), with 79 evaluable plasma concentrations fitted in NONMEM 7.5.0 by FOCE-I ",
    "(ADVAN1 TRANS2). Clearance carries a fixed allometric weight exponent of 0.75 and an ",
    "estimated power effect of postnatal age (exponent 0.433), both normalised to the ",
    "typical infant of 3.4 kg and 36 days; central volume is linear in weight over the ",
    "same 3.4 kg reference. Postnatal age was the only covariate retained: gestational ",
    "age, postmenstrual age, height, body surface area, albumin, serum creatinine and ",
    "baseline transaminases were screened and rejected, and they are recorded in ",
    "covariatesDataExcluded rather than in the model. Between-subject variability could ",
    "only be estimated for clearance. Oxacillin is dosed directly into `central`. The ",
    "unbound concentration Cu that drives the paper's fT > MIC target-attainment analysis ",
    "is derived algebraically as 10% of the total concentration, an unbound fraction the ",
    "authors fixed from the literature rather than measuring. The magnitudes of the ",
    "combined additive-plus-proportional residual error are not reported anywhere in the ",
    "paper or its supplements, so both residual standard deviations are encoded as ",
    "fixed(0); see the vignette Errata."
  )
  reference <- paste(
    "Lee A, Liu C, Tran MT, Phal S, Peloquin CA, Nieves D, Capparelli E, Arrieta AC.",
    "(2026). Population pharmacokinetics and safety of continuous oxacillin in preterm",
    "and term neonates and infants.",
    "Antimicrobial Agents and Chemotherapy 70(6).",
    "doi:10.1128/aac.01777-25.",
    sep = " "
  )
  vignette <- "Lee_2026_oxacillin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "oxacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Cohort median 3.29 kg, range 0.605-7.025 kg (Supplemental Table S2, column ",
        "'Current weight (kg)'). The published equations normalise weight to 3.4 kg, the ",
        "weight of the 'typical infant 36 days postnatal age and weighing 3.4 kg' named in ",
        "Lee 2026 Results, 'Pharmacokinetic model'. Clearance uses the fixed allometric ",
        "exponent 0.75 and central volume is linear in weight (exponent 1). Weight was ",
        "carried into the model before any other covariate was screened."
      ),
      source_name = "WT_Kg"
    ),
    PNA = list(
      description = "Postnatal (chronological) age since birth",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Lee 2026 reports postnatal age in DAYS (cohort median 32 days, IQR 12-46, Table ",
        "1; raw range 4-82 days in Supplemental Table S2) and normalises the clearance ",
        "covariate to 36 days. The canonical PNA column carries MONTHS ",
        "(inst/references/covariate-columns.md), so model() converts back with ",
        "1 month = 30.4375 days before forming the age ratio -- the same reparameterisation ",
        "used by Zhao_2018_omeprazole.R (days) and Bardhi_2026_ampicillin_foal.R (hours). ",
        "The reference 36 days therefore corresponds to 1.18275 months. Treated as ",
        "time-fixed at enrolment in the source analysis, whose PK sampling spans at most ",
        "96 h."
      ),
      source_name = "PNA_days"
    )
  )

  # Screened by Lee 2026 (Materials and methods, 'Pharmacokinetic model') but NOT retained
  # in the final model, so they must not appear in model(). Recorded here so the paper's
  # covariate screen is preserved without a 'declared but not referenced' convention
  # warning. None of these has a published point estimate -- the paper reports only that
  # the OFV drop failed the 3.84 threshold and/or that the bootstrap CI spanned no effect.
  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Cohort median 38 weeks (IQR 30-39), Lee 2026 Table 1. Used to define the four ",
        "enrolment cohorts (<32 vs >=32 weeks) but not retained as a model covariate. The ",
        "Discussion is explicit about why: only one premature neonate under 14 days of age ",
        "could be enrolled (Table 1, Cohort 4, n = 1), 'thus, we had a limited ability to ",
        "assess the impact of GA at birth'."
      ),
      source_name = "GA (weeks)"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Screened and not retained. Lee 2026 Discussion: 'postmenstrual age was also not ",
        "included in the final model as it is driven by the PNA component of postmenstrual ",
        "age' -- i.e. it was collinear with the retained PNA term. Note that the canonical ",
        "PAGE column carries months; this entry records the source's weeks convention only ",
        "because the covariate is documentation, never referenced in model()."
      ),
      source_name = "postmenstrual age"
    ),
    HT = list(
      description = "Body length (recumbent height)",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Raw values 32-56.5 cm (Supplemental Table S2, column 'Length (cm)'). Screened as ",
        "'height (cm)' in Lee 2026 Materials and methods and not retained."
      ),
      source_name = "Length (cm)"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in Lee 2026 Materials and methods and not retained; no derivation formula is given.",
      source_name = "body surface area"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Cohort median 3.2 g/dL (IQR 3.0-3.6), Lee 2026 Table 1. Note the unit: the ",
        "canonical ALB column is g/L, and this entry preserves the source's US convention ",
        "because the covariate is documentation only. Screened and not retained. The ",
        "Conclusion flags albumin as a plausible driver of the unmeasured protein binding ",
        "rather than of clearance."
      ),
      source_name = "Serum albumin (g/dL)"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Cohort median 0.3 mg/dL (IQR 0.2-0.36), Lee 2026 Table 1. Screened and not ",
        "retained. The Discussion gives three reasons, all about the covariate's ",
        "resolution rather than its biology: values were tightly distributed in infants ",
        "over 14 days, were highly correlated with the retained PNA term, and were ",
        "'only reported to one significant digit with a lower limit of detection of ",
        "0.2 mg/dL'. It also warns that creatinine above 0.4 mg/dL in the under-14-day ",
        "infants 'may partially reflect maternal creatinine rather than neonatal ",
        "clearance'."
      ),
      source_name = "Serum creatinine (mg/dL)"
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as 'baseline serum transaminases' in Lee 2026 Materials and methods and not retained.",
      source_name = "Baseline AST"
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as 'baseline serum transaminases' in Lee 2026 Materials and methods and not retained.",
      source_name = "Baseline ALT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    n_observations = 79L,
    age_range = "4-82 days postnatal (inclusion criteria >3 to <=90 days)",
    age_median = "32 days postnatal (IQR 12-46)",
    ga_range = "23.9-40.3 weeks gestational age at birth",
    ga_median = "38 weeks (IQR 30-39)",
    weight_range = "0.605-7.025 kg",
    weight_median = "3.29 kg",
    sex_female_pct = 27,
    race_ethnicity = c(Hispanic = 73, White = 14, Black = 5, `Hispanic/Black` = 5, `Hispanic/White` = 5),
    disease_state = paste0(
      "Hospitalized neonates and young infants receiving oxacillin as standard of care: ",
      "11 empiric therapy for suspected infection, 9 confirmed methicillin-susceptible ",
      "Staphylococcus aureus infection (osteomyelitis n = 1, pneumonia n = 3, bacteraemia ",
      "n = 3, skin/soft tissue n = 2), 1 methicillin-susceptible S. epidermidis ",
      "bacteraemia, 1 perioperative prophylaxis. Critically ill infants could not be ",
      "enrolled, and renal dysfunction (dialysis, urine output <0.5 mL/kg/h, or serum ",
      "creatinine >1.7 mg/dL), transaminases >5x ULN, therapeutic hypothermia within 24 h ",
      "and ECMO were exclusion criteria, so the model should not be extrapolated to those ",
      "states."
    ),
    dose_range = paste0(
      "25 mg/kg intravenous loading dose over 30 min, immediately followed by a continuous ",
      "intravenous infusion of 160 mg/kg/day (21 of 22 infants) or 120 mg/kg/day (1 infant ",
      "in Cohort 4: gestational age <32 weeks and postnatal age <14 days)"
    ),
    sampling = paste0(
      "Convenience sampling, maximum five samples per infant (minimum 25 uL each): a ",
      "baseline sample before the loading dose and then 30-120 min, 8-16 h and 16-96 h ",
      "after the start of the continuous infusion, plus one within 1 h after the end of ",
      "infusion. Group 1 and four further infants had at most three samples to limit ",
      "phlebotomy. Of 89 plasma samples, 10 were below the 2 mcg/mL lower limit of ",
      "quantification (9 of them baseline) and were excluded, leaving 79. Plasma assayed ",
      "by validated LC-MS/MS at the University of Florida Infectious Disease ",
      "Pharmacokinetics Laboratory over a 2-100 mcg/mL calibration range with within- and ",
      "between-day CV <10%. One cerebrospinal-fluid concentration was collected and was ",
      "NOT included in the final model."
    ),
    regions = "United States (single centre, Children's Hospital of Orange County, California)",
    notes = paste0(
      "Prospective, phase 1, open-label, single-centre study. Demographics are Lee 2026 ",
      "Table 1; per-patient raw values are Supplemental Table S2. Enrolment was stratified ",
      "into four cohorts by gestational age (<32 vs >=32 weeks) and postnatal age (<14 vs ",
      ">=14 days) targeting 6 infants each, but Cohort 4 (GA <32 weeks, PNA <14 days) ",
      "recruited only 1 of 6. Model reliability was assessed by 1,000-set bootstrapping in ",
      "Wings for NONMEM; every final estimate fell inside its bootstrap 95% confidence ",
      "interval. eta-shrinkage on clearance was 8.27%. Table 1's total male count (15, ",
      "68%) disagrees with the sum of its own per-cohort entries (8 + 3 + 4 + 1 = 16) and ",
      "with Supplemental Table S2, which lists 16 male and 6 female infants; the 27% ",
      "female recorded here follows the raw data. Median measured total oxacillin ",
      "concentration 20.2 mcg/mL (IQR 11.8-36.2). Safety: 45 adverse events in 18 infants, ",
      "12 possibly oxacillin-related, 3 serious events none of them related, no deaths."
    )
  )

  ini({
    # =================================================================
    # Disposition -- Lee 2026 Table 2 ('Final pharmacokinetic model
    # parameters') and the final-model equations printed in Results,
    # 'Pharmacokinetic model' and repeated in the Table 2 footnote b:
    #
    #   CL (L/h) = 1.01 * (WT_Kg / 3.4)^0.75 * (PNA_days / 36)^0.433
    #   Vd (L)   = 1.87 * WT_Kg / 3.4
    #
    # Both are absolute (L/h and L) at the reference infant of 3.4 kg
    # and 36 days postnatal age.
    # =================================================================
    lcl <- log(1.01)
    label("Clearance for a 3.4 kg, 36-day-old infant (L/h)")  # Lee 2026 Table 2, CL (Theta1) = 1.01 (bootstrap 95% CI 0.791-1.3, bootstrap median 1.01)

    lvc <- log(1.87)
    label("Central compartment volume of distribution for a 3.4 kg infant (L)")  # Lee 2026 Table 2, V (Theta2) = 1.87 (bootstrap 95% CI 1.07-2.61, bootstrap median 1.78)

    # -----------------------------------------------------------------
    # Size scaling. The 0.75 exponent on clearance is the theoretical
    # allometric value: Lee 2026 Materials and methods states clearance
    # 'was allometrically scaled by weight (WT_Kg) and included before
    # assessment of other covariates', and no exponent appears in Table
    # 2 with an estimate, a shrinkage or a bootstrap interval. The
    # volume exponent is 1 because the printed Vd equation is linear in
    # weight, not a power function. Both are therefore fixed, not
    # estimated.
    # -----------------------------------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on clearance (unitless)")  # Lee 2026 Table 2 footnote b: (WT_Kg / 3.4)^0.75; a-priori allometric value, not estimated

    e_wt_vc <- fixed(1)
    label("Exponent of body weight on central volume (unitless)")  # Lee 2026 Table 2 footnote b: Vd = 1.87 * WT_Kg / 3.4, i.e. linear in weight

    # -----------------------------------------------------------------
    # Postnatal-age maturation of clearance. This is the only covariate
    # the paper retained; Results quantifies it as 'a 2.2-fold increase
    # in CL (weight adjusted) from 14 to 90 days of age', which the
    # exponent reproduces: (90/14)^0.433 = 2.22.
    # -----------------------------------------------------------------
    e_pna_cl <- 0.433
    label("Power exponent of postnatal age on clearance (unitless)")  # Lee 2026 Table 2, 'Effect of PNA on CL (Theta3)' = 0.433 (bootstrap 95% CI 0.133-0.709, bootstrap median 0.422)

    # -----------------------------------------------------------------
    # Unbound fraction. NOT measured in this study. Lee 2026 Materials
    # and methods, 'Simulations and pharmacodynamic assessments':
    # 'Free-drug concentrations of oxacillin were estimated at 10% of
    # total concentrations as reported in the literature (9)'. The
    # Discussion warns this is likely an UNDER-estimate of the true free
    # fraction in preterm infants. It is a literature constant, so it is
    # fixed and carries no variability.
    # -----------------------------------------------------------------
    fu <- fixed(0.1)
    label("Fraction of oxacillin unbound in plasma (unitless)")  # Lee 2026 Materials and methods, 'Simulations and pharmacodynamic assessments'; Introduction cites ~90% protein binding from reference 9

    # =================================================================
    # Between-subject variability. Lee 2026 Results: 'Between-subject
    # variability (BSV) could only be estimated for CL', so there is no
    # eta on volume.
    #
    # SCALE OF THE PUBLISHED 0.568 -- this is the one reading the paper
    # does not state, and it is resolved here from the paper's own data
    # rather than assumed. Table 2 prints 'IIV CL (etaCL) 0.568' with no
    # unit and no 'variance' / 'SD' / 'CV' qualifier. Read as a NONMEM
    # OMEGA variance it implies an SD of sqrt(0.568) = 0.754 on the log
    # scale; read as a standard deviation it is 0.568 itself. Two
    # independent checks against on-disk sources both exclude the
    # variance reading:
    #
    #   1. Supplemental Table S2 gives each infant's dose rate, weight,
    #      postnatal age and steady-state concentrations. Under a
    #      continuous infusion CL_i = rate / Css_i, so the individual
    #      clearances are recoverable directly. The spread of
    #      log(CL_i / CL_predicted) about the published covariate model
    #      is SD 0.571 using the 16-96 h samples and 0.589 using the
    #      per-infant mean of the 8-16 h and 16-96 h samples. That
    #      spread contains BOTH between-subject variability AND residual
    #      error, so it is an UPPER bound on the BSV SD. 0.754 exceeds
    #      it; 0.568 sits just underneath it, exactly as it should.
    #
    #   2. The Figure 3 target-attainment curves pin the SD without
    #      needing the level. For a continuous infusion the paper's
    #      fT 100% > MIC criterion is simply Css_unbound > MIC, so
    #      probit(PTA) is linear in log(MIC) with slope -1 / SD. Reading
    #      the resolved mid-range segments of the red 160 mg/kg/day
    #      curves gives SD 0.61 and 0.55 (PNAD 90), 0.61 (PNAD 28), 0.53
    #      (PNAD 14) and 0.56 (PNAD 7) -- mean 0.57, and no panel
    #      approaches 0.754. The same reading also reproduces the level:
    #      the PNAD 90 curve crosses 50% at MIC 1.45 mg/L against 1.50
    #      predicted for a cohort-median 3.29 kg infant.
    #
    # 0.568 is therefore the standard deviation of eta_CL on the log
    # scale. nlmixr2's ini() takes VARIANCES, so the entry below is
    # 0.568^2. See the vignette section 'Reading the published
    # between-subject variability', which reproduces both checks as
    # executable code.
    #
    # Table 2 source row: 'IIV CL (etaCL)' final estimate 0.568,
    # bootstrap 95% CI 0.373-0.683, bootstrap median 0.54, eta-shrinkage
    # 8.27 percent. (Kept off the trailing comment below: rxode2 promotes
    # a trailing comment on an ini() line that carries no label() into
    # label(), and a bare percentage inside a label reads as a claimed
    # back-transform.)
    # =================================================================
    etalcl ~ 0.322624   # Lee 2026 Table 2, IIV CL read as a log-scale SD of 0.568; variance = 0.568^2

    # =================================================================
    # Residual error. Lee 2026 Materials and methods, 'Pharmacokinetic
    # model': 'A combined additive and proportional within-subject error
    # model was chosen to characterize the residual error.' The STRUCTURE
    # is therefore known, but neither magnitude is reported: Table 2 has
    # no sigma row, the Results text gives none, and none of the three
    # supplemental workbooks (S1 adverse events, S2 raw data, S3 raw
    # concentrations) contains one. Per the standing policy on
    # unreported variability, both are encoded as fixed(0) rather than
    # invented; simulations from this model are therefore free of
    # residual noise and reproduce IPRED, not DV.
    # =================================================================
    addSd <- fixed(0)
    label("Additive residual error SD (ug/mL; ZERO - magnitude not reported in source)")  # Lee 2026 declares a combined error model but publishes no sigma estimates

    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - magnitude not reported in source)")  # Lee 2026 declares a combined error model but publishes no sigma estimates
  })

  model({
    # 1. Derived covariate term. The canonical PNA covariate column is in
    #    MONTHS; Lee 2026 reports postnatal age in DAYS and normalises to
    #    36 days, so convert here. 1 month = 30.4375 days.
    pnaDays  <- PNA * 30.4375
    ageRatio <- pnaDays / 36

    # 2. Individual parameters.
    #      CL = 1.01 * (WT/3.4)^0.75 * (PNA_days/36)^0.433
    #      Vd = 1.87 * (WT/3.4)
    cl <- exp(lcl + etalcl) * (WT / 3.4)^e_wt_cl * ageRatio^e_pna_cl
    vc <- exp(lvc) * (WT / 3.4)^e_wt_vc

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system. One compartment with first-order elimination
    #    (NONMEM ADVAN1 TRANS2). Oxacillin is given intravenously -- a
    #    30 min loading infusion followed by a continuous infusion -- so
    #    it is dosed straight into `central` with no absorption step.
    d/dt(central) <- -kel * central

    # 5. Observation and error. Cc is the TOTAL plasma concentration the
    #    assay measured; Cu is the unbound concentration that drives the
    #    paper's fT > MIC target, taken as a fixed 10% of total.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd)
  })
}
