Han_2024_ainuovirine <- function() {
  description <- paste0(
    "Two-compartment population PK model with first-order absorption, an ",
    "absorption lag time, dose-level-dependent relative bioavailability ",
    "and a step increase in apparent clearance between the first dose and ",
    "steady state, for the third-generation non-nucleoside reverse ",
    "transcriptase inhibitor ainuovirine (ANV) in treatment-naive adults ",
    "living with HIV-1 (Han 2024, n = 337, 1947 plasma concentrations ",
    "pooled from the phase 1 multiple-ascending-dose trial ",
    "ADYY-ACC007-103 (ChiCTR1800018022; 75, 150 or 300 mg orally once ",
    "daily for 10 days) and the phase 3 trial ADYY-ACC007-301 ",
    "(ChiCTR1800019041; 150 mg once daily at bedtime for 48 weeks)). ANV ",
    "is nonlinear in dose: the elimination half-life is essentially the ",
    "same at 75, 150 and 300 mg, so Han 2024 attributes the ",
    "less-than-proportional exposure to ABSORPTION and carries the ",
    "nonlinearity entirely in relative bioavailability -- F = 1 at 75 mg ",
    "(reference), 0.716 at 150 mg and 0.410 at 300 mg. Apparent clearance ",
    "carries a second, orthogonal step: CL/F at steady state is 2.47-fold ",
    "the first-dose value (a 147% increase, Table 3 'Drugno on CL' = ",
    "1.47), which the paper attributes chiefly to auto-induction of ",
    "CYP3A4 by ANV. No covariate reached significance -- age, sex, weight, ",
    "BMI, ALT, AST, total bilirubin, albumin, creatinine, creatinine ",
    "clearance, total cholesterol and concomitant medication were all ",
    "screened and rejected, so the base model IS the final model. ",
    "Companion exposure-response models in the Han_2024_ainuovirine_* ",
    "family cover the virologic and adverse-drug-reaction endpoints."
  )
  reference <- paste(
    "Han X, Sun J, Zhang Y, Jiang T, Zheng Q, Peng H, Wang Y, Xia W,",
    "Zhang T, Sun L, Yun X, Qin H, Wu H, Su B.",
    "Population pharmacokinetics of Ainuovirine and exposure-response",
    "analysis in human immunodeficiency virus-infected individuals.",
    "Chin Med J (Engl). 2024;137(20):2474-2482.",
    "doi:10.1097/CM9.0000000000002917.",
    sep = " "
  )
  vignette <- "Han_2024_ainuovirine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Han 2024 Methods (oral once-daily
  # dosing; venous plasma ANV concentrations by the study bioanalytical
  # assay) and Table 3 (Vc/F and Vp/F reported in L).
  compartmentData <- list(
    depot       = list(analyte = "ainuovirine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ainuovirine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ainuovirine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_ANV_MG = list(
      description        = "Administered ainuovirine dose level for the record's dosing event, in mg.",
      units              = "mg",
      type               = "continuous",
      reference_category = "75 mg (the dose level at which relative bioavailability is anchored to 1)",
      notes              = paste(
        "Used ONLY to select which of the three studied dose levels the",
        "record belongs to; it is not a continuous regressor. Han 2024",
        "estimates relative bioavailability as a piecewise (categorical)",
        "function of dose level -- Methods state 'categorical covariates",
        "were described using a piecewise model' -- with F = 1 at 75 mg,",
        "0.716 at 150 mg and 0.410 at 300 mg (Table 3, rows 'F 150 mg' and",
        "'F 300 mg', whose footnote reads 'F 150 mg is the bioavailability",
        "of 150 mg relative to 75 mg'). model() therefore bands this column",
        "at the arithmetic midpoints of the studied levels (<= 112.5 mg ->",
        "75 mg; 112.5-225 mg -> 150 mg; > 225 mg -> 300 mg) rather than",
        "interpolating. The model is calibrated ONLY at 75, 150 and 300 mg;",
        "a dose outside those levels is snapped to the nearest studied one",
        "and the resulting F carries no support from the source data. This",
        "matters more than usual here because the F contrast IS the",
        "nonlinearity -- doubling 150 mg to 300 mg raises steady-state",
        "exposure by only 0.410/0.716 = 1.15-fold, not 2-fold.",
        "Time-varying in principle (a subject whose dose level changed",
        "would carry different values), but constant within subject in both",
        "source trials.",
        "Distinct from MULTI_DOSE_PT below: DOSE_ANV_MG selects the",
        "bioavailability level (a between-dose-group contrast), while",
        "MULTI_DOSE_PT selects the clearance level (a within-subject",
        "first-dose-vs-steady-state contrast). The two are orthogonal and",
        "both apply to the phase 1 records."
      ),
      source_name        = "Dose group (75 mg / 150 mg / 300 mg)"
    ),
    MULTI_DOSE_PT = list(
      description        = "Multiple-dose (steady-state) record indicator: 1 = the record was collected after repeated once-daily ainuovirine dosing had reached steady state; 0 = the record was collected after the first dose.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (first-dose records)",
      notes              = paste(
        "Per-record indicator switching apparent clearance via the",
        "multiplicative form cl * (1 + 1.47 * MULTI_DOSE_PT), i.e. CL/F at",
        "steady state is 2.47-fold the first-dose value. Han 2024 Table 3",
        "reports the coefficient on the row 'Drugno on CL' = 1.47 (RSE",
        "21.40%, 95% CI 0.85-2.09) and states the resulting multiplier",
        "outright in the table footnote: 'CL (steady-state) pop = CL",
        "typical x 2.47'; the Results and Discussion both restate it as",
        "'a 147.0% increase in steady-state CL/F over the first dose'.",
        "The paper's own source column name is Drugno (dose number).",
        "Mechanism per the Discussion: auto-induction of CYP3A4 -- ANV",
        "induces CYP3A enzyme activity and CYP3A4 mRNA in primary human",
        "hepatocytes -- with disease progression and reduced renal plasma",
        "clearance offered as partial contributors.",
        "The switch is a STEP, not a time course: Han 2024 estimates no",
        "induction rate constant and no half-life of onset, because the",
        "phase 1 design samples only two occasions (after the first dose on",
        "day 1, and after the last dose on day 10). Records between those",
        "landmarks are not informed by the source data. Assign 0 to day-1",
        "records and 1 to records at or after steady state; all phase 3",
        "ADYY-ACC007-301 samples (weeks 12, 24, 36 and 48) are steady-state",
        "records and take 1.",
        "Record-level, not subject-level: a phase 1 subject contributes",
        "both levels."
      ),
      source_name        = "Drugno (dose number; first dose vs steady state)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "The ONLY screened covariate that produced a significant OFV drop,",
        "and it was still rejected. Han 2024 Results: 'The inclusion of age",
        "as a covariate resulted in a significant decrease in the OFV of",
        "the PopPK model, but the estimated value was close to 0 (-0.00625)",
        "and the inter-individual random effect of clearance (CL) decreased",
        "by only 0.7. Therefore, age was not included in the final model.'",
        "The Discussion attributes the null result to the narrow age range",
        "studied (18-61 years) and contrasts it with doravirine, another",
        "NNRTI, for which age IS a significant covariate. Cohort mean 31.18",
        "years (SD 9.84) (Table 2)."
      )
    ),
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and rejected (Han 2024 Results: 'sex, weight, BMI, ALT, AST, TB, ALB, CRE, CrCL, CHOL, and COMB were not influential covariates'). Cohort mean 66.86 kg (SD 11.17) (Table 2). Note that no allometric scaling is applied in the final model -- CL/F and the volumes are population typical values, not weight-normalised."
    ),
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and rejected. The Discussion flags the limitation directly: only 18 of 337 subjects (5.34%) were female, so the analysis had little power to detect a sex effect (Table 2)."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 22.56 kg/m^2 (SD 3.19) (Table 2)."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 26.50 IU/L (SD 16.59) (Table 2). The Discussion cautions that few participants with hepatic impairment were enrolled, so the null result does not license extrapolation to hepatic impairment."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 23.61 IU/L (SD 8.46) (Table 2)."
    ),
    TBIL = list(
      description = "Total bilirubin at baseline.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 12.34 umol/L (SD 4.99) (Table 2). Source abbreviation TB."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 47.08 g/L (SD 2.96) (Table 2)."
    ),
    CREAT = list(
      description = "Serum creatinine at baseline.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 68.11 umol/L (SD 10.49) (Table 2). Source abbreviation CRE."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 78.12 mL/min (SD 21.86) (Table 2). As with ALT/AST, the Discussion notes that few participants with renal impairment were enrolled, so this null result does not license extrapolation to renal impairment -- especially given that the Discussion elsewhere offers reduced renal plasma clearance as a partial explanation for the steady-state CL/F increase."
    ),
    TCHOL = list(
      description = "Total cholesterol at baseline.",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened and rejected. Cohort mean 4.20 mmol/L (SD 0.74) (Table 2). Source abbreviation CHOL."
    ),
    CONMED_ANY = list(
      description = "Any concomitant medication indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and rejected. Han 2024 lists 'combination medication (COMB)' among the tested covariates but does not report which drugs it covered, how many subjects carried it, or its coding, so it cannot be reconstructed as a usable model term."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 337L,
    n_studies      = 2L,
    n_observations = "1947 plasma ainuovirine concentrations. 341 participants took ANV orally; 4 whose plasma concentrations were all below the lower limit of quantification were treated as missing, leaving 337 in the PopPK dataset (28 from ADYY-ACC007-103, 309 from ADYY-ACC007-301)",
    age_range      = "18-61 years (range quoted in the Discussion); mean 31.18 years (SD 9.84) overall, by group 29.63 (6.37) at 75 mg, 34.70 (8.43) at phase 1 150 mg, 33.40 (9.90) at 300 mg, 31.03 (9.95) in phase 3 (Han 2024 Table 2)",
    weight_range   = "mean 66.86 kg (SD 11.17) overall; by group 62.13 (9.61), 64.96 (6.97), 66.20 (7.42) and 67.06 (11.42) kg (Han 2024 Table 2). Range not reported",
    sex_female_pct = 5.34,
    race_ethnicity = "Not reported. Both trials were conducted at a single Chinese centre (Beijing Youan Hospital, Capital Medical University), so the cohort is effectively Chinese Han-majority, but Han 2024 tabulates no race or ethnicity data",
    disease_state  = "Antiretroviral-therapy-naive people living with HIV-1 (PLWH). Baseline HIV-RNA and CD4 counts are not tabulated in Han 2024; the phase 3 efficacy endpoint was HIV-RNA < 50 copies/mL at week 48",
    dose_range     = "Phase 1 (ADYY-ACC007-103): 75 mg (n = 8), 150 mg (n = 10) or 300 mg (n = 10) orally once daily for 10 days. Phase 3 (ADYY-ACC007-301): 150 mg once daily at bedtime on an empty stomach for 48 weeks (n = 309), as part of an ANV-based regimen compared against an efavirenz-based regimen",
    regions        = "Single centre, Beijing, China",
    baseline_labs  = "ALT 26.50 IU/L (SD 16.59); AST 23.61 (8.46) IU/L; total bilirubin 12.34 (4.99) umol/L; albumin 47.08 (2.96) g/L; creatinine 68.11 (10.49) umol/L; creatinine clearance 78.12 (21.86) mL/min; total cholesterol 4.20 (0.74) mmol/L; BMI 22.56 (3.19) kg/m^2 (Han 2024 Table 2)",
    sampling       = "Phase 1: pre-dose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 168 and 192 h after the FIRST dose, then pre-dose and 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36, 72 and 120 h after the LAST (day 10) dose. Phase 3: sparse sampling at weeks 12, 24, 36 and 48 (Han 2024 Table 1)",
    notes          = paste0(
      "The Discussion names two limitations that bound reuse: (1) the ",
      "narrow age range and the very low proportion of female ",
      "participants limited the power of the covariate screen, so the ",
      "all-null covariate result is a statement about this cohort rather ",
      "than about ANV pharmacology generally; (2) all data come from a ",
      "single institution. Only a few participants with hepatic or renal ",
      "impairment were enrolled, so the model should not be used to ",
      "predict exposure in organ impairment."
    )
  )

  ini({
    # ==================================================================
    # Han 2024 Table 3, "Final model / Estimates (RSE [%])" column. Every
    # value below is that column; the 95% CI and the 1000-replicate
    # bootstrap median + 95% CI in the adjacent columns agree throughout
    # (the Results note that they "largely overlapped").
    #
    # SD-versus-VARIANCE. The Results prose and Table 3 report the random
    # effects on DIFFERENT scales, and Table 3 is the one to encode:
    #
    #   quantity          Table 3 (SD)     Results prose      SD^2
    #   omega(CL)         30.9%            "9.54%"            0.0954
    #   sigma(prop)       27.9%            "7.8%"             0.0778
    #   sigma(add)        8.89 ng/mL       "78.99 ng/mL"      79.03
    #
    # All three prose figures are the squares of the Table 3 figures, so
    # the prose is quoting variances. Three independent facts settle it
    # in favour of Table 3: (a) the RSEs printed alongside each pair are
    # identical (3.70 and 27.80), so these are the same rows and not
    # different quantities; (b) Table 3's own 95% CIs bracket the Table 3
    # values (30.9 in 28.4-33.4; 27.9 in 25.9-29.9; 8.89 in 4.05-13.70)
    # and would not bracket the prose values; (c) the bootstrap medians
    # (30.9, 27.7, 8.24) track the Table 3 values. See the vignette
    # Errata section.
    #
    # NO ALLOMETRY AND NO COVARIATES. Han 2024's covariate screen
    # rejected everything (Results: "Consequently, the basic model was
    # considered the final model"), so CL/F, Vc/F, Vp/F and Q/F below are
    # population typical values for the whole cohort, NOT values at a
    # reference weight. Do not add a (WT/70)^0.75 term when reusing this
    # model -- it is not in the source.
    # ==================================================================

    # ----- Structural disposition -----
    lcl <- log(6.46);   label("Apparent clearance after the FIRST dose (CL/F, L/h)")     # Han 2024 Table 3, CL/F = 6.46 (RSE 15.00%), 95% CI 4.56-8.36; bootstrap median 6.35 (4.17-8.36). Steady-state CL/F is this value x 2.47 -- see e_md_cl
    lvc <- log(11.5);   label("Apparent central volume of distribution (Vc/F, L)")       # Han 2024 Table 3, Vc/F = 11.5 (RSE 13.7%), 95% CI 8.4-14.6; bootstrap median 11.4 (8.86-15.90)
    lvp <- log(293.0);  label("Apparent peripheral volume of distribution (Vp/F, L)")    # Han 2024 Table 3, Vp/F = 293.0 (RSE 10.5%), 95% CI 233.0-353.0; bootstrap median 297.0 (240.0-377.0)
    lq  <- log(17.6);   label("Apparent intercompartmental clearance (Q/F, L/h)")        # Han 2024 Table 3, Q/F = 17.6 (RSE 11.0%), 95% CI 13.8-21.4; bootstrap median 17.7 (14.5-22.8)

    # ----- Absorption -----
    # KA is the one structural parameter Han 2024 did NOT estimate in the
    # pooled fit. Table 3 marks the row with an asterisk and footnotes it
    # "KA fixed as the typical value of KA for the PopPK model
    # constructed from ADYY-ACC007-103 study data", prints no RSE and no
    # confidence interval, and the bootstrap column simply repeats the
    # same number. The Discussion gives the reason: pooling the sparse
    # phase 3 data (sampled in the distribution and elimination phases
    # only) with the phase 1 data left several parameters poorly
    # estimated, and "the accuracy of the model parameters was
    # significantly improved by directly adopting the population typical
    # values of KA of the 103 study modeling".
    #
    # Table 3's unit for this row reads "L/h", which is a typographical
    # error -- KA is a first-order rate constant. The Results text gives
    # it correctly as "0.0985 h -1".
    lka   <- fixed(log(0.0985)); label("First-order absorption rate constant (KA, 1/h)")  # Han 2024 Table 3, KA = 0.0985, carried over from the ADYY-ACC007-103-only model
    ltlag <- log(0.208);         label("Absorption lag time (ALAG, h)")                   # Han 2024 Table 3, ALAG = 0.208 (RSE 26.700%), 95% CI 0.099-0.317; bootstrap median 0.216 (0.083-0.309)

    # ----- Relative bioavailability by dose level -----
    # The 75 mg group is the reference: Table 3's footnote defines both
    # estimated rows against it ("F 150 mg is the bioavailability of
    # 150 mg relative to 75 mg; F 300 mg is the bioavailability of 300 mg
    # relative to 75 mg"), so F(75 mg) is a structural anchor of 1 rather
    # than an estimate.
    lfdepot           <- fixed(log(1)); label("Relative bioavailability at the 75 mg reference dose level (F, unitless)")  # Han 2024 Table 3 footnote: the 150 mg and 300 mg rows are expressed relative to 75 mg, so 75 mg anchors F at 1
    e_dose150_fdepot  <- 0.716;         label("Relative bioavailability at 150 mg versus the 75 mg reference (F, unitless)")  # Han 2024 Table 3, F 150 mg = 0.716 (RSE 10.600%), 95% CI 0.567-0.865; bootstrap median 0.715 (0.570-0.900)
    e_dose300_fdepot  <- 0.410;         label("Relative bioavailability at 300 mg versus the 75 mg reference (F, unitless)")  # Han 2024 Table 3, F 300 mg = 0.410 (RSE 12.900%), 95% CI 0.306-0.514; bootstrap median 0.411 (0.311-0.549)

    # ----- First-dose versus steady-state clearance -----
    e_md_cl <- 1.47; label("Multiplicative effect of MULTI_DOSE_PT on apparent clearance (unitless; CL ~ cl * (1 + e_md_cl * MULTI_DOSE_PT))")  # Han 2024 Table 3, "Drugno on CL" = 1.47 (RSE 21.40%), 95% CI 0.85-2.09; bootstrap median 1.51 (1.02-2.83). Table 3 footnote: "CL (steady-state) pop = CL typical x 2.47", i.e. 1 + 1.47

    # ----- Inter-individual variability -----
    # Han 2024 estimated IIV on clearance ONLY. The Results state that
    # "the distribution of individual random effects was close to normal
    # basically, and there were no significant correlated pharmacokinetic
    # parameters", and Table 3 carries a single "Inter-individual
    # omega (CL)" row. Shrinkage 5.0% (Table 3 footnote).
    #
    # Table 3's 30.9% is the SD on the log scale (see the SD-versus-
    # variance note above), so the variance rxode2 wants is 0.309^2 =
    # 0.0954 -- which is exactly the "9.54%" the Results prose quotes.
    # Note this is omega^2 read straight off the paper, NOT a
    # log(CV^2 + 1) back-transform: the two agree to 5% here
    # (log(1 + 0.309^2) = 0.0911) but the paper's own arithmetic is the
    # authority.
    etalcl ~ 0.0954  # Han 2024 Table 3, omega(CL) = 30.9% (RSE 4.1%), 95% CI 28.4-33.4%; squared to the variance the Results prose reports as 9.54%

    # ----- Residual error -----
    # Combined proportional + additive. Shrinkage 7.0% for both
    # (Table 3 footnote).
    propSd <- 0.279; label("Proportional residual error (fraction)")   # Han 2024 Table 3, sigma(Prop) = 27.9% (RSE 3.7%), 95% CI 25.9-29.9%; bootstrap median 27.7 (25.8-29.7)
    addSd  <- 8.89;  label("Additive residual error (ng/mL)")          # Han 2024 Table 3, sigma(Add) = 8.89 ng/mL (RSE 27.80%), 95% CI 4.05-13.70; bootstrap median 8.24 (4.00-13.10)
  })

  model({
    # ==================================================================
    # 1. Derived covariate terms
    # ==================================================================
    # Dose-level indicators. Han 2024 fits relative bioavailability as a
    # piecewise function of the three studied dose levels, so DOSE_ANV_MG
    # is banded at the arithmetic midpoints (112.5 and 225 mg) rather
    # than interpolated. Exactly one indicator is 1 for any dose.
    is_dose150 <- (DOSE_ANV_MG > 112.5) * (DOSE_ANV_MG <= 225)
    is_dose300 <- (DOSE_ANV_MG > 225)
    is_dose75  <- 1 - is_dose150 - is_dose300

    fdepot <- exp(lfdepot) *
      (is_dose75 + is_dose150 * e_dose150_fdepot + is_dose300 * e_dose300_fdepot)

    # ==================================================================
    # 2. Individual PK parameters
    # ==================================================================
    # The (1 + e_md_cl * MULTI_DOSE_PT) factor is a STEP between the
    # first dose (MULTI_DOSE_PT = 0, CL/F = 6.46 L/h) and steady state
    # (MULTI_DOSE_PT = 1, CL/F = 6.46 * 2.47 = 15.96 L/h). Han 2024
    # estimates no time course for the transition.
    cl <- exp(lcl + etalcl) * (1 + e_md_cl * MULTI_DOSE_PT)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    tlag <- exp(ltlag)

    # ==================================================================
    # 3. Micro-constants
    # ==================================================================
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ==================================================================
    # 4. ODE system
    # ==================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ==================================================================
    # 5. Bioavailability and lag time
    # ==================================================================
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ==================================================================
    # 6. Observation and error
    # ==================================================================
    # Dose is carried in mg and vc in L, so central/vc is mg/L; the
    # factor of 1000 converts to the ng/mL that Han 2024 reports (and in
    # which the additive residual 8.89 is expressed).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
