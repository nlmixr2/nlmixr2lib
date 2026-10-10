Takada_2022_roxadustat <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for oral roxadustat (HIF-prolyl hydroxylase inhibitor) in Japanese dialysis-dependent chronic kidney disease patients with anaemia, on haemodialysis or peritoneal dialysis (Takada 2022). Age >= 65 years lowers CL/F; five concomitant phosphate binders (sevelamer hydrochloride/bixalomer, calcium carbonate, lanthanum carbonate, ferric citrate, sucroferric oxyhydroxide) each lower the relative bioavailability, with an additional decrease when phosphate binders were taken without time separation (phase II study)."
  reference <- "Takada A, Shibata T, Shiga T, Groenendaal-van de Meent D, Komatsu K. Population pharmacokinetics of roxadustat in Japanese dialysis-dependent chronic kidney disease patients with anaemia. Br J Clin Pharmacol. 2022;88(2):787-797. doi:10.1111/bcp.15023"
  vignette <- "Takada_2022_roxadustat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "roxadustat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "roxadustat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "roxadustat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only as the categorical threshold age >= 65 years (Methods 2.5: 'Categorical age (>= 65 or <65 y) was tested in addition to continuous age'; Results 3.4 and Table 4 row 'Age >= 65 y for CL'). A subject aged exactly 65 years is in the >= 65 group. Baseline value per Methods 2.5.",
      source_name = "Age"
    ),
    CONMED_SEVELAMER_BIXALOMER = list(
      description = "Concomitant sevelamer hydrochloride or bixalomer (polymer phosphate binders, pooled), 1 = coadministered / 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant sevelamer hydrochloride or bixalomer",
      notes = "Time-varying: derived from daily administration records (Methods 2.5). Sevelamer hydrochloride and bixalomer were grouped 'due to the same mechanism of action and similar effects on bioavailability'. Multiplicative factor 0.744 on F1 (Table 4 row 'SBUSE on F1'). Concomitant-medication doses were not considered.",
      source_name = "SBUSE"
    ),
    CONMED_CALCIUM_CARBONATE = list(
      description = "Concomitant calcium carbonate used as a phosphate binder, 1 = coadministered / 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant calcium carbonate",
      notes = "Time-varying from daily administration records. Multiplicative factor 0.931 on F1 (Table 4 row 'CUSE on F1'). Used here as a phosphate binder (chelation of roxadustat in the gut), not as an antacid; see CONMED_ANTACID for the gastric-pH covariate.",
      source_name = "CUSE"
    ),
    CONMED_LANTHANUM_CARBONATE = list(
      description = "Concomitant lanthanum carbonate (phosphate binder), 1 = coadministered / 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant lanthanum carbonate",
      notes = "Time-varying from daily administration records. Multiplicative factor 0.969 on F1 (Table 4 row 'LUSE on F1').",
      source_name = "LUSE"
    ),
    CONMED_FERRIC_CITRATE = list(
      description = "Concomitant ferric citrate hydrate (phosphate binder), 1 = coadministered / 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant ferric citrate",
      notes = "Time-varying from daily administration records. Multiplicative factor 0.744 on F1 (Table 4 row 'FUSE on F1'). Oral iron given as an iron supplement was a separate covariate that was screened and not retained.",
      source_name = "FUSE"
    ),
    CONMED_SUCROFERRIC_OXYHYDROXIDE = list(
      description = "Concomitant sucroferric oxyhydroxide (phosphate binder), 1 = coadministered / 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant sucroferric oxyhydroxide",
      notes = "Time-varying from daily administration records. Multiplicative factor 0.837 on F1 (Table 4 row 'SCUSE on F1').",
      source_name = "SCUSE"
    ),
    STUDY_PHASE2 = list(
      description = "Record from the phase II study 1517-CL-0304, in which roxadustat was taken without the 1-hour time separation from phosphate binders required in the phase III studies; 1 = phase II / 0 = phase III",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = phase III studies (1517-CL-0308, -0307, -0302), with phosphate binders taken at least 1 hour before or after roxadustat",
      notes = "Multiplies F1 by 0.935 (Table 4 row 'No time separation of PBs in phase 2 on F1'; Table 4 footnote 'the additional effect of the phosphate binders due to the uncontrolled intake') only while at least one of the five phosphate binders is coadministered. To simulate a patient who takes phosphate binders together with roxadustat, set STUDY_PHASE2 = 1; set 0 for the staggered (>= 1 h separation) dosing used in the phase III studies and the label.",
      source_name = "Study 1517-CL-0304"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female / 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Screened on CL/F (Table 2); not retained."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F (Table 2, power model normalised to the mean); not retained."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F (Table 2); not retained."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F (Table 2) because roxadustat is ~99% albumin-bound; not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Hepatic parameter screened on CL/F (Table 2 footnote b); not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Hepatic parameter screened on CL/F (Table 2 footnote b); not retained."
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Hepatic parameter screened on CL/F (Table 2 footnote b); not retained."
    ),
    TPRO = list(
      description = "Baseline total protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Hepatic parameter screened on CL/F (Table 2 footnote b); not retained."
    ),
    PERIT_DIAL = list(
      description = "Dialysis modality, 1 = peritoneal dialysis / 0 = haemodialysis",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = haemodialysis",
      notes = "Type of dialysis screened on CL/F (Table 2, Results 3.2); not significant, so haemodialysis and peritoneal dialysis patients share one model."
    ),
    CONMED_CLOPIDOGREL = list(
      description = "Concomitant clopidogrel (CYP2C8 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concomitant clopidogrel",
      notes = "Screened on CL/F (Table 2); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 367L,
    n_studies = 4L,
    n_observations = 1285L,
    age_mean = "64.4 years (SD 10.9)",
    weight_mean = "59.5 kg (SD 11.8)",
    bmi_mean = "22.8 kg/m^2 (SD 3.5)",
    sex_female_pct = 30.2,
    race_ethnicity = "Japanese",
    disease_state = "Dialysis-dependent chronic kidney disease with renal anaemia: haemodialysis (studies 1517-CL-0304, -0308, -0307; n = 311) or peritoneal dialysis (1517-CL-0302; n = 56)",
    dose_range = "Oral roxadustat three times weekly; initial 50, 70 or 100 mg titrated on haemoglobin to 20-250 mg",
    regions = "Japan",
    co_medication = "At least one phosphate binder in 301 patients (82.0%): sevelamer hydrochloride/bixalomer 79, calcium carbonate 175, lanthanum carbonate 139 (Table 3; the Results text says 131), ferric citrate 59, sucroferric oxyhydroxide 17",
    notes = "Pooled phase II (1517-CL-0304, n = 92) and three phase III (1517-CL-0308 n = 74, 1517-CL-0307 n = 145, 1517-CL-0302 n = 56) studies; Table 1 and Table 3. Sparse sampling at study visits at any time relative to dose and dialysis. LLOQ 1.00 ng/mL; 15 BLQ samples (1.2%) excluded."
  )

  ini({
    # Structural parameters -- Takada 2022 Table 4 (final model). Reference
    # patient: age < 65 years, no concomitant phosphate binder.
    lcl <- log(0.923); label("Apparent clearance CL/F (L/h)") # Table 4 'CL/F (L/h)' 0.923 (95% CI 0.796-1.05; RSE 7.0%)
    lvc <- log(14.6); label("Apparent central volume Vc/F (L)") # Table 4 'Vc/F (L)' 14.6 (95% CI 11.9-17.3; RSE 9.5%)
    lka <- log(0.63); label("First-order absorption rate constant ka (1/h)") # Table 4 'Ka (1/h)' 0.63 (95% CI 0.405-0.855; RSE 18.3%)
    lq <- log(0.134); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 4 'Q/F (L/h)' 0.134 (95% CI 0.00601-0.262; RSE 48.7%)
    lvp <- log(2.89); label("Apparent peripheral volume Vp/F (L)") # Table 4 'Vp/F (L)' 2.89 (95% CI 1.59-4.19; RSE 23.0%)
    ltlag <- log(0.287); label("Absorption lag time ALAG1 (h)") # Table 4 'ALAG1 (h)' 0.287 (95% CI 0.273-0.301; RSE 2.6%)

    # Categorical covariate effects, proportional-change model P = theta1 * theta2^(0 or 1)
    # (Methods 2.5 equation). Each is a multiplicative factor applied when the indicator is 1.
    e_age_ge65_cl <- 0.792; label("Multiplicative factor on CL/F for age >= 65 years (unitless)") # Table 4 'Age >= 65 y for CL' 0.792 (95% CI 0.716-0.868; RSE 4.9%)
    e_conmed_sevelamer_bixalomer_f <- 0.744; label("Multiplicative factor on F1 for sevelamer hydrochloride/bixalomer (unitless)") # Table 4 'SBUSE on F1' 0.744 (95% CI 0.606-0.882; RSE 9.5%)
    e_conmed_calcium_carbonate_f <- 0.931; label("Multiplicative factor on F1 for calcium carbonate (unitless)") # Table 4 'CUSE on F1' 0.931 (95% CI 0.799-1.063; RSE 7.2%)
    e_conmed_lanthanum_carbonate_f <- 0.969; label("Multiplicative factor on F1 for lanthanum carbonate (unitless)") # Table 4 'LUSE on F1' 0.969 (95% CI 0.837-1.101; RSE 7.0%)
    e_conmed_ferric_citrate_f <- 0.744; label("Multiplicative factor on F1 for ferric citrate (unitless)") # Table 4 'FUSE on F1' 0.744 (95% CI 0.602-0.886; RSE 9.7%)
    e_conmed_sucroferric_oxyhydroxide_f <- 0.837; label("Multiplicative factor on F1 for sucroferric oxyhydroxide (unitless)") # Table 4 'SCUSE on F1' 0.837 (95% CI 0.490-1.184; RSE 21.1%)
    e_study_phase2_f <- 0.935; label("Multiplicative factor on F1 for phosphate binders taken without time separation, phase II (unitless)") # Table 4 'No time separation of PBs in phase 2 on F1' 0.935 (95% CI 0.888-0.982; RSE 2.6%)

    # IIV: Table 4 prints %CV = 100 * sqrt(omega^2) (the printed 95% CIs are
    # omega^2 +/- 1.96 * RSE * omega^2 back-transformed by sqrt, e.g. ka
    # 3.389 * (1 -/+ 1.96 * 0.152) -> 154.3-209.7%). Diagonal OMEGA.
    etalcl ~ 0.173889 # Table 4 'IIV CL/F (%CV)' 41.7% -> 0.417^2 (shrinkage 6.9%)
    etalvc ~ 0.0484 # Table 4 'IIV Vc/F (%CV)' 22.0% -> 0.220^2 (shrinkage 69.2%)
    etalka ~ 3.389281 # Table 4 'IIV ka (%CV)' 184.1% -> 1.841^2 (shrinkage 50.7%)

    # Residual error: Yij = IPRED + W * eps, W = sqrt((IPRED * theta_prop)^2 + theta_add^2),
    # EPS fixed to variance 1 (Methods 2.4), so both thetas are SDs.
    propSd <- 0.439; label("Proportional residual error SD (fraction)") # Table 4 'Proportional error (CV%)' 43.9% (95% CI 41.3-46.5%; RSE 3.1%)
    addSd <- 1.88; label("Additive residual error SD (ng/mL)") # Table 4 'Additive error (ng/mL)' 1.88 (95% CI 0.892-2.87; RSE 26.8%)
  })
  model({
    # Age >= 65 years (categorical; Results 3.4).
    age_ge65 <- 0
    if (AGE >= 65) age_ge65 <- 1

    # Any of the five phosphate binders coadministered at this time.
    pb_any <- 1 -
      (1 - CONMED_SEVELAMER_BIXALOMER) *
        (1 - CONMED_CALCIUM_CARBONATE) *
        (1 - CONMED_LANTHANUM_CARBONATE) *
        (1 - CONMED_FERRIC_CITRATE) *
        (1 - CONMED_SUCROFERRIC_OXYHYDROXIDE)

    # Relative bioavailability (reference F1 = 1, no phosphate binder).
    fdepot <- e_conmed_sevelamer_bixalomer_f^CONMED_SEVELAMER_BIXALOMER *
      e_conmed_calcium_carbonate_f^CONMED_CALCIUM_CARBONATE *
      e_conmed_lanthanum_carbonate_f^CONMED_LANTHANUM_CARBONATE *
      e_conmed_ferric_citrate_f^CONMED_FERRIC_CITRATE *
      e_conmed_sucroferric_oxyhydroxide_f^CONMED_SUCROFERRIC_OXYHYDROXIDE *
      e_study_phase2_f^(STUDY_PHASE2 * pb_any)

    cl <- exp(lcl + etalcl) * e_age_ge65_cl^age_ge65
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    q <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
