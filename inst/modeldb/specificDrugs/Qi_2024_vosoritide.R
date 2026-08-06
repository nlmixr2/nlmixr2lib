Qi_2024_vosoritide <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination and",
    "change-point first-order absorption (the absorption rate constant",
    "switches from Ka1 to Ka2 at an estimated time after each dose) for",
    "subcutaneous vosoritide (BMN 111, a C-type natriuretic peptide analog)",
    "in children with achondroplasia aged 0.95-15 years, pooled from five",
    "clinical trials (Qi 2024). Body weight is a power covariate on CL/F",
    "(exponent 0.356) and on V/F (exponent 1.09), both referenced to 20 kg.",
    "Relative bioavailability rises exponentially with time on treatment and",
    "is 56% higher for the 0.2 mg/mL dosing solution used only in study",
    "111-202. Study-level random effects nested inside the subject-level IIV",
    "on CL/F and V/F reproduce the paper's secondary study identity number",
    "(SIDN) hierarchy, and separate log-scale residual errors are carried for",
    "the ELISA and the electrochemiluminescence assays. The model was used to",
    "derive the eight-band weight-band dosing regimen of Qi 2024 Table 6.")
  reference   <- "Qi Y, Chan ML, Mould DR, Larimore K, Fisheleva E, Cherukuri A, Day J, Savarirayan R, Irving M, Bacino CA, Hoover-Fong J, Ozono K, Mohnike K, Wilcox WR, Bober MB, Henshaw J. Development of a weight-band dosing approach for vosoritide in children with achondroplasia using a population pharmacokinetic model. Clin Pharmacokinet. 2024;63(5):707-719. doi:10.1007/s40262-024-01371-6"
  vignette    <- "Qi_2024_vosoritide"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    depot   = list(analyte = "vosoritide", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vosoritide", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on both CL/F and V/F with reference weight 20 kg (Qi 2024 Sect. 2.4, 'The reference body weight was established as 20 kg'). Applied as (WT/20)^0.356 on CL/F and (WT/20)^1.09 on V/F, reproducing Qi 2024 Table 2 exactly at 9, 20, 40, 60 and 74.5 kg. Qi 2024 Discussion reports that time-varying weight was tested alongside baseline weight and alternative size metrics (LBW, FFM, BMI, BSA); the published final model carries body weight, and the paper does not state whether the retained column is baseline or time-varying, so either may be supplied.",
      source_name        = "WT"
    ),
    FORM_VOSO_SOLN02 = list(
      description        = "Indicator that the dose was prepared from the 0.2 mg/mL vosoritide dosing solution, 1 = 0.2 mg/mL, 0 = 0.8 or 2 mg/mL.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the 0.8 mg/mL and 2 mg/mL dosing solutions, which Qi 2024 Discussion states 'did not have any effect on F'; these are the phase III and commercial strengths).",
      notes              = "Per-dose-record indicator. The 0.2 mg/mL solution was used only in study 111-202 (Qi 2024 Sect. 2.3). Multiplies relative bioavailability by 1.56 (Qi 2024 Table 5 'Effect of SOLNC (0.2 mg/mL)' = 1.56; Table 3 shows the resulting relative bioavailability at 1.56 at time 0 and a constant 56% increase over the reference at every time point). Encoded on the log scale as e_form_voso_soln02_fdepot = log(1.56).",
      source_name        = "SOLNC"
    ),
    ELISA = list(
      description        = "Bioanalytical assay indicator, 1 = plasma vosoritide measured by the validated enzyme-linked immunosorbent assay (ELISA), 0 = measured by the optimized electrochemiluminescence (ECL) assay.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECL assay; studies 111-206, 111-301 and 111-302, LLOQ 0.137 ug/L).",
      notes              = "Per-observation indicator. ELISA = 1 for studies 111-202 and 111-205 (LLOQ 0.391 ug/L); ELISA = 0 for studies 111-206, 111-301 and 111-302 (Qi 2024 Sect. 2.3). Selects the log-scale residual SD: 0.665 for ELISA (Qi 2024 Table 5 'Residual error 1') and 0.610 for ECL ('Residual error 2'). Assay choice is study-fixed, but the indicator is carried per observation so pooled datasets apply the correct residual to each row.",
      source_name        = "assay type (Qi 2024 Sect. 2.3 / Table 5 footnotes a and b)"
    ),
    SIDN = list(
      description        = "Secondary study identity number: an integer grouping column identifying the clinical trial a subject's records belong to, used as the nesting level for the study-level random effects on CL/F and V/F.",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- SIDN is a nesting level, not a covariate with a reference category. It indexes the second level of random effects (eta6 on CL/F, eta7 on V/F) rather than entering any typical-value equation.",
      notes              = "Qi 2024 Sect. 3.2: 'An additional secondary study identity number (SIDN) was used to represent the clinical trial in which each patient was enrolled, but allowed for the fact that patients may have been enrolled in > 1 study. The effects of SIDN on the IIV of CL/F and V/F were modeled by an additional hierarchical level of effect (StudyCL and StudyV).' The paper states there were only three SIDN values in the analysis database but does not enumerate which trials map to which value, so the column is deliberately opaque: supply 1, 2 or 3. Implemented with rxode2/nlmixr2 native nested random effects (`etalcl_study ~ ... | SIDN`), which places a separate draw per SIDN level on top of the per-subject etas.",
      source_name        = "SIDN"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.", units = "years", type = "continuous",
      notes = "Screened graphically and in single-covariate models (Qi 2024 Sect. 2.4, Table S3) but not retained in the final model; Qi 2024 Key Points state 'no other covariates were identified as being predictive'. Baseline range 0.95-15 years, mean 8.43 (Table 1)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.", units = "(binary)", type = "binary",
      notes = "Screened but not retained (Qi 2024 Key Points). Table 1: 84 male, 74 female."
    ),
    BMI = list(
      description = "Body mass index.", units = "kg/m^2", type = "continuous",
      notes = "Screened as an alternative body-size metric (Qi 2024 Discussion) but not retained; body weight was the retained size descriptor."
    ),
    BSA = list(
      description = "Body surface area (Mosteller).", units = "m^2", type = "continuous",
      notes = "Screened as an alternative body-size metric, including a paediatric BSA variant computed from lean body weight (Qi 2024 Discussion, Table 1 'BSA2'); not retained."
    ),
    FFM = list(
      description = "Fat-free mass.", units = "kg", type = "continuous",
      notes = "Screened as an alternative body-size metric (Qi 2024 Discussion); not retained."
    ),
    LBW = list(
      description = "Lean body weight.", units = "kg", type = "continuous",
      notes = "Screened as an alternative body-size metric, both directly and as the input to BMI and BSA variants (Qi 2024 Discussion); not retained."
    ),
    ALB = list(
      description = "Serum albumin.", units = "g/dL (values in Qi 2024 Table 1 are tabulated in g/L)", type = "continuous",
      notes = "Screened in the covariate evaluation (Qi 2024 Table 1 baseline summary, Table S3) but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase.", units = "IU/L", type = "continuous",
      notes = "Screened (Qi 2024 Table 1) but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase.", units = "IU/L", type = "continuous",
      notes = "Screened (Qi 2024 Table 1) but not retained."
    ),
    BILI = list(
      description = "Total bilirubin.", units = "umol/L", type = "continuous",
      notes = "Screened (Qi 2024 Table 1) but not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance.", units = "mL/min", type = "continuous",
      notes = "Screened (Qi 2024 Table 1) but not retained."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate.", units = "mL/min/1.73m^2", type = "continuous",
      notes = "Screened (Qi 2024 Table 1) but not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 158L,
    n_studies      = 5L,
    age_range      = "0.95-15 years (mean 8.43, median 8; Qi 2024 Table 1)",
    weight_range   = "9-74.5 kg (mean 23.8, median 22.2; Qi 2024 Table 1)",
    sex_female_pct = 46.8,
    race_ethnicity = c(White = 72.2, Black = 3.8, Asian = 17.7, Other = 6.3),
    disease_state  = "Children with achondroplasia (autosomal dominant FGFR3 gain-of-function skeletal dysplasia) with open epiphyses.",
    dose_range     = "2.5 ug/kg/day (6 patients), 7.5 ug/kg/day (12 patients), 15 ug/kg/day (151 patients) and 30 ug/kg/day (11 patients) subcutaneous, once daily (Qi 2024 Sect. 3.1).",
    regions        = "Multi-national BioMarin development program: NCT02055157 (111-202), NCT02724228 (111-205), NCT03197766 (111-301), NCT03424018 (111-302), NCT03583697 (111-206).",
    n_observations = "4741 concentrations retained from an initial 6181 (23.3% excluded, mostly pre-dose samples expected to be non-measurable given vosoritide's short half-life; Qi 2024 Sect. 3.1, Table S4).",
    notes          = "One patient was excluded for a substantial number of oscillations in weight over time. Study 111-206 (children <= 5 years) contributed interim data from sentinel patients only (8 subjects in Table 1), so the low-body-weight end of the covariate range is sparsely informed; Qi 2024 flags this as a study limitation and plans a model update. Table 1 lists 29 / 8 / 60 / 61 subjects for studies 111-202 (with its 111-205 extension) / 111-206 / 111-301 / 111-302."
  )

  ini({
    # ==================================================================
    # Final model parameter estimates from Qi 2024 Table 5 (typical
    # values with %SE and bootstrap 2.5th / median / 97.5th percentiles).
    # Table 5 reports back-transformed typical values; the paper's
    # Sect. 2.4 equations show how the thetas enter:
    #
    #   CL/F      = exp(theta1 + eta1 + eta6) * (WT(Kg)/20)^theta14
    #   V/F       = exp(theta2 + eta2 + eta7) * (WT(Kg)/20)^theta15
    #   Ka1       = exp(theta5 + eta5)
    #   Ka2       = exp(theta6 + eta8)
    #   Change-point = theta7 + eta9
    #   SD1 = theta8 ; SD2 = theta9
    #   F         = exp(theta16 * Time/10000)
    #
    # eta6 / eta7 are the additional hierarchical (study) level indexed
    # by SIDN. eta5 and eta8 (on Ka1 and Ka2) were not estimated: Qi
    # 2024 Sect. 3.2 states the model "accounted for the IIV in CL/F,
    # V/F, and change-point" and Table 5 has no Ka IIV rows, so Ka1 and
    # Ka2 are typical-value only here.
    #
    # Scale of the "(CV, %)" column in Table 5.  Every variability row
    # is reported under a single "(CV, %)" heading.  The change-point
    # IIV is the row that settles the scale: it is the only *fixed*
    # variance in the table (22.4, %SE 0, bootstrap NE), and
    # 100 * sqrt(0.05) = 22.36 -> 22.4, i.e. the reported number is
    # 100 * sqrt(omega^2) for a round fixed omega^2 of 0.05.  The
    # alternative log-normal back-transform 100 * sqrt(exp(0.05) - 1)
    # gives 22.6, which does not match.  So throughout this file the
    # Table 5 percentages are read as 100 x SD on the parameter's own
    # scale: omega^2 = (CV/100)^2 and the residual SDs are CV/100.
    # This also agrees with the residual rows, where the paper's own
    # equations name the estimated quantity a standard deviation
    # ("SD1 = theta8", "SD2 = theta9") rather than a CV.
    # ==================================================================

    # ---------------- Structural (typical-value) parameters ----------------
    lcl <- log(47.47)
    label("Log apparent clearance CL/F at the 20 kg reference weight (log(L/h))")     # Qi 2024 Table 5: CL/F = 47.47 L/h (%SE 1.8; bootstrap 40.45 / 47.7 / 57.97). Corroborated by Table 2 row 'WT 20 (reference)' = 47.47 L/h = 100.00% of reference.
    lvc <- log(17.99)
    label("Log apparent central volume of distribution V/F at the 20 kg reference weight (log(L))")  # Qi 2024 Table 5: V/F = 17.99 L (%SE 3.3; bootstrap 15.18 / 18.36 / 23.09). Corroborated by Table 2 row 'WT 20 (reference)' = 17.99 L = 100.00% of reference.
    lka1 <- log(2.21)
    label("Log first-order absorption rate constant Ka1, active before the change-point (log(1/h))")  # Qi 2024 Table 5: Ka1 = 2.21 1/h (%SE 14.7; bootstrap 1.8 / 2.3 / 3.19)
    lka2 <- log(0.06)
    label("Log first-order absorption rate constant Ka2, active from the change-point onward (log(1/h))")  # Qi 2024 Table 5: Ka2 = 0.06 1/h (%SE 3.8; bootstrap 0.04 / 0.06 / 0.09)
    lmtime <- log(0.31)
    label("Log time after dose at which the absorption rate constant changes from Ka1 to Ka2 (log(h))")  # Qi 2024 Table 5: Change-point = 0.31 h (%SE 7.3; bootstrap 0.27 / 0.3 / 0.33). Estimated in NONMEM with the MTIME feature (Sect. 2.4).

    # ---------------- Covariate effects (Qi 2024 Table 5 covariate rows) ----------------
    e_wt_cl <- 0.356
    label("Power exponent of (WT/20) on CL/F (unitless)")                              # Qi 2024 Table 5 'Effect of weight on CL/F' = 0.356 (%SE 25.1; bootstrap 0.02 / 0.35 / 0.68). Verified against Table 2: (9/20)^0.356 = 0.7526 -> 75.26% and (74.5/20)^0.356 = 1.5971 -> 159.71%, matching the tabulated percentages exactly.
    e_wt_vc <- 1.09
    label("Power exponent of (WT/20) on V/F (unitless)")                               # Qi 2024 Table 5 'Effect of weight on V/F' = 1.09 (%SE 8.1; bootstrap 0.85 / 1.1 / 1.31). Verified against Table 2: (9/20)^1.09 = 0.4188 -> 41.88% and (74.5/20)^1.09 = 4.1930 -> 419.30%, matching the tabulated percentages exactly.
    e_form_voso_soln02_fdepot <- log(1.56)
    label("Log relative bioavailability multiplier for the 0.2 mg/mL dosing solution (unitless)")  # Qi 2024 Table 5 'Effect of SOLNC (0.2 mg/mL)' = 1.56 (%SE 15.3; bootstrap 1.18 / 1.54 / 2.03). Table 3 lists this as a constant 56% increase over the 0.8 / 2.0 mg/mL reference at every time point, i.e. a multiplicative factor of 1.56 on F.
    e_time_fdepot <- 0.21
    label("Slope of log relative bioavailability on time after starting dose, per 10000 h (unitless)")  # Qi 2024 Table 5 'Effect of time on F' = 0.21 (%SE 8.1; bootstrap 0.16 / 0.21 / 0.26), entering as F = exp(0.21 * Time/10000) per the Sect. 2.4 equation. Verified against Table 4: exp(0.21 * 1000/10000) = 1.021 -> 2.12% and exp(0.21 * 25000/10000) = 1.690 -> 69.05%, matching the tabulated values.

    # ---------------- Subject-level between-subject variability (Table 5 IIV rows) ----------------
    # omega^2 = (CV/100)^2 per the scale note above.
    etalcl ~ 0.112896                                                                 # Qi 2024 Table 5 'IIV CL (CV, %)' = 33.6 (%SE 8.3; bootstrap 27.3 / 33 / 40.49); 0.336^2 = 0.112896
    etalvc ~ 0.058564                                                                 # Qi 2024 Table 5 'IIV V (CV, %)' = 24.2 (%SE 13.6; bootstrap 17.7 / 23.9 / 29.99); 0.242^2 = 0.058564
    etalmtime ~ fixed(0.05)                                                           # Qi 2024 Table 5 'IIV change-point (fixed)' = 22.4, %SE 0, bootstrap NE; 100*sqrt(0.05) = 22.36 -> 22.4

    # ---------------- Study-level (nested) variability (Table 5 study IIV rows) ----------------
    # The additional hierarchical level of Qi 2024 Sect. 2.4 (eta6 on
    # CL/F, eta7 on V/F), indexed by the SIDN grouping column. Precision
    # is poor for these two terms because only three SIDN values were
    # present in the database (Qi 2024 Sect. 3.2).
    etalcl_study ~ 0.066049 | SIDN                                                    # Qi 2024 Table 5 'IIV study CL (CV, %)' = 25.7 (%SE 47.4; bootstrap 18.51 / 26.15 / 34.6); 0.257^2 = 0.066049
    etalvc_study ~ 0.000144 | SIDN                                                    # Qi 2024 Table 5 'IIV study V (CV, %)' = 1.2 (%SE 51.8; bootstrap 0.4 / 2.05 / 5.7); 0.012^2 = 0.000144

    # ---------------- Residual variability (Table 5 residual rows) ----------------
    # A log-transform-both-sides approach was used (Qi 2024 Sect. 2.4),
    # so the residual is additive on the natural-log concentration
    # scale, i.e. log-normal in linear space -> nlmixr2 lnorm().
    expSdElisa <- 0.665
    label("Log-scale residual SD for concentrations measured by the ELISA assay (log(ug/L))")  # Qi 2024 Table 5 'Residual error 1' = 66.5 (%SE 1.6; bootstrap 61.71 / 66.3 / 72); footnote a: ELISA assay (studies 111-202, 111-205)
    expSdEcl <- 0.610
    label("Log-scale residual SD for concentrations measured by the electrochemiluminescence assay (log(ug/L))")  # Qi 2024 Table 5 'Residual error 2' = 61 (%SE 1.5; bootstrap 56.61 / 60.85 / 65.4); footnote b: ECL assay (studies 111-206, 111-301, 111-302)
  })

  model({
    # ------------------------------------------------------------------
    # Reference body weight for the two power covariates (Qi 2024
    # Sect. 2.4: "The reference body weight was established as 20 kg").
    # ------------------------------------------------------------------
    ref_wt <- 20

    # ------------------------------------------------------------------
    # Individual PK parameters. Qi 2024 Sect. 2.4 equations, with the
    # subject-level eta (eta1 / eta2) and the study-level eta (eta6 /
    # eta7) both entering additively on the log scale.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl + etalcl_study) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc + etalvc_study) * (WT / ref_wt)^e_wt_vc

    ka1 <- exp(lka1)
    ka2 <- exp(lka2)

    # Change-point time after dose. Qi 2024 prints "Change-point =
    # theta7 + eta9" (additive), but the fixed variability on that eta
    # is 22.4% and an additive eta with SD 0.224 h on a 0.31 h
    # change-point would be a 72% relative spread and would go negative
    # for roughly 8% of subjects. Read as a log-normal parameter
    # (omega^2 = 0.05 fixed, i.e. 22.4% CV) the spread is 22.4% -- which
    # is what the table reports -- and positivity is guaranteed. The two
    # readings give practically identical SDs on the change-point
    # (0.0693 h log-normal vs 0.0694 h additive).
    t_change <- exp(lmtime + etalmtime)

    kel <- cl / vc

    # ------------------------------------------------------------------
    # Relative bioavailability (Qi 2024 Sect. 2.4 and Tables 3-4):
    #   F = 1.56^(SOLNC == 0.2 mg/mL) * exp(theta16 * Time / 10000)
    # "Time" is time after starting dose in hours (Table 4 cross-tabs
    # 1000 h against 0.11 years), so tafd() is the right clock. The
    # printed F equation carries only the time term; the SOLNC factor
    # comes from Table 5's "Effect of SOLNC (0.2 mg/mL)" = 1.56 and is
    # confirmed by Table 3, where the 0.2 mg/mL column is exactly 1.56x
    # the Table 4 time-only column at every tabulated time.
    # ------------------------------------------------------------------
    f(depot) <- exp(e_form_voso_soln02_fdepot * FORM_VOSO_SOLN02 +
                      e_time_fdepot * tafd() / 10000)

    # ------------------------------------------------------------------
    # Change-point first-order absorption (Qi 2024 Sect. 2.4): "If the
    # condition 'time after dose < change-point' is true, then the value
    # of absorption rate (Ka) is set to Ka1. Otherwise, the value is set
    # to Ka2." NONMEM used MTIME to place the switch exactly.
    #
    # The switch is written INLINE in the d/dt right-hand sides rather
    # than assigned to a local `ka`: rxode2 evaluates standalone model()
    # assignments once per output record and holds them constant across
    # the integration interval, which would make the absorption profile
    # depend on the observation grid. Inline terms are re-evaluated at
    # every solver step. `t - tlast()` is time after the most recent
    # dose (equivalent to tad()).
    # ------------------------------------------------------------------
    d/dt(depot) <- -(ka1 * ((t - tlast()) < t_change) +
                       ka2 * ((t - tlast()) >= t_change)) * depot
    d/dt(central) <- (ka1 * ((t - tlast()) < t_change) +
                        ka2 * ((t - tlast()) >= t_change)) * depot -
      kel * central

    # NOTE: deliberately no record-level `ka <- ka1 * (tad() < t_change) +
    # ka2 * (tad() >= t_change)` reporting copy. rxode2 5.1.3 normalises
    # `t - tlast()` to `tad()`, so such a line becomes syntactically
    # identical to the inline d/dt switch above; common-subexpression
    # elimination then folds the ODE right-hand sides onto that single
    # `ka` symbol, the system matches the one-compartment-with-depot
    # pattern, and rxSolve's automatic ODE -> linCmt conversion replaces
    # the whole ODE system with an analytic solution driven by a
    # record-level `ka`. That silently reintroduces the grid dependence
    # this inline form exists to avoid. Simulations should additionally
    # pass useLinCmt = FALSE to rxSolve().

    # ------------------------------------------------------------------
    # Observation. Doses are in ug and volumes in L, so central / vc is
    # in ug/L, the unit used for concentrations throughout Qi 2024
    # (Fig. 1 caption: "Concentration is given in ug/L").
    # ------------------------------------------------------------------
    Cc <- central / vc

    # Assay-dependent log-scale residual SD, reproducing the paper's two
    # separate residual errors as a linear switch on the ELISA
    # indicator (same encoding pattern as Cirincione_2017_exenatide).
    expSd <- expSdElisa * ELISA + expSdEcl * (1 - ELISA)

    Cc ~ lnorm(expSd)
  })
}
