Valadez_2025_cefepime <- function() {
  description <- "Two-compartment population PK model for intravenous cefepime in mechanically ventilated adults admitted to the medical ICU with suspected hospital-acquired pneumonia, with and without extracorporeal membrane oxygenation (ECMO) and excluding concurrent renal replacement therapy. Total clearance scales with raw Cockcroft-Gault creatinine clearance through a non-linear power term standardised to 120 mL/min; central volume scales linearly with total body weight standardised to 70 kg and is multiplied by exp(beta2) in patients cannulated onto ECMO, a 2.8-fold expansion. Inter-compartmental distribution is parameterised directly as the micro-constants KCP and KPC rather than as Q and Vp. Estimated with the Pmetrics non-parametric adaptive grid (NPAG); the source reports only the weighted median of the non-parametric population distribution with 95% credible intervals, so no between-subject variance and no residual-error model are available (see the vignette Assumptions and deviations). Valadez 2025, n = 70 patients (9 on ECMO), 114 plasma samples."
  reference <- "Valadez A, Zurawska M, Harlan E, Scheetz MH, Neely MN, Yarnold PR, Kang M, Korth E, Martinez F, Giblin B, Donnelly HK, Dedicatoria K, Medernach R, Nozick S, Hauser AR, Ozer EA, Diaz E, Misharin AV, Wunderink RG, Rhodes NJ. Individual target pharmacokinetic/pharmacodynamic attainment rates among cefepime-treated patients admitted to the ICU with hospital-acquired pneumonia with and without ECMO. Antimicrob Agents Chemother. 2025;69(6):e0010225. doi:10.1128/aac.00102-25. PMID 40372025; PMC12135513."
  vignette <- "Valadez_2025_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefepime is administered as an intravenous infusion
  # and total (not unbound) plasma cefepime was assayed by LC-MS/MS
  # (Valadez 2025, "Cefepime dosing, sample collection, and bioanalysis").
  compartmentData <- list(
    central     = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation from serum creatinine, age, sex and total body weight",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT BSA-normalized to mL/min/1.73 m^2 (Valadez 2025 Methods, 'Design and patients': 'Creatinine clearance (CrCl) was calculated using the method of Cockcroft and Gault equation'). Enters clearance as the non-linear power term (CRCL / 120)^beta1 with beta1 = 0.76 (Eq 1); the 120 mL/min standardisation is stated in the Table 2 footnote ('creatinine clearance normalized clearance ... (CrCl/120 mL/min)^beta1'). Cohort baseline mean +/- SD 115.8 +/- 88.6 mL/min (Table 1); the very large SD reflects an ICU population spanning renal impairment to augmented renal clearance. Patients on hemodialysis, peritoneal dialysis or CRRT were EXCLUDED from the analysis, so the covariate is never paired with a renal-replacement indicator in this model. Covariates were carried at multiple time points in the source analysis (Fig. 1 legend: 'the Bayesian posterior distribution based on covariate values across all time points'), so CRCL is naturally time-varying; the Results text quotes posterior estimates 'based on covariate values at time zero'. Must be strictly positive -- the covariate enters a power term.",
      source_name        = "CrCl"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Total body weight (TBW), not ideal or adjusted body weight (Valadez 2025 Table 1 footnote: 'Total body weight (TBW)'). Enters central volume LINEARLY as the ratio (WT / 70), i.e. an allometric exponent fixed at exactly 1 rather than an estimated or 0.75 exponent (Eq 2, and the Table 2 footnote 'weight-normalized (WT/70 kg) volume of distribution (V1 in L)'). Cohort baseline mean +/- SD 83.3 +/- 26.5 kg (Table 1). Weight was extracted from the electronic health record alongside the other covariates and is treated as available at each covariate time point.",
      source_name        = "TBW"
    ),
    ECMO_STATUS = list(
      description        = "Extracorporeal membrane oxygenation treatment-status indicator (1 = cannulated onto ECMO, 0 = not on ECMO)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not receiving ECMO)",
      notes              = "Binary classifier, stated explicitly in the source: 'ECMO is treated as a binary classifier (0 = no ECMO, 1 = ECMO)' (Valadez 2025 Results, 'Model development and selection'). 9 of 70 patients (12.9%) required ECMO during treatment (Table 1). Modality was predominantly veno-venous (VV) ECMO for respiratory support (Discussion). Enters central volume as the EXPONENTIAL multiplier exp(beta2 * ECMO_STATUS) -- see Eq 2 and Table S1 run 8 ('Model 6 with ECMO on V: V1 * e(beta2*ECMO)'). Retained on Vd only: an ECMO effect on CL improved fit individually (dOFV = 10.2) but was rejected for low precision (CV > 200%) with point estimates near zero, and the combined CL + Vd model was worse than Vd alone (Table S1 runs 7 and 9). Patients on concurrent renal replacement therapy were excluded, so this indicator is not confounded by CRRT in this cohort. Treated as a subject-level indicator here, matching how the source reports it; the source does not resolve within-subject cannulation / decannulation timing (the Discussion lists 'the timing of ECMO initiation ... was not standardized across patients' as a limitation).",
      source_name        = "ECMO"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 70L,
    n_studies        = 1L,
    n_samples        = 114L,
    age_mean_sd      = "62.1 +/- 14.4 years (mean +/- SD; Table 1). Full range not reported.",
    weight_mean_sd   = "83.3 +/- 26.5 kg total body weight (mean +/- SD; Table 1). Full range not reported.",
    bsa_mean_sd      = "1.96 +/- 0.33 m^2 (mean +/- SD; Table 1)",
    sex_female_pct   = 40,
    race_ethnicity   = "Not reported. The source tabulates only age, total body weight, BSA, serum creatinine, creatinine clearance, ECMO status and sex (Table 1).",
    disease_state    = "Mechanically ventilated adults admitted to the Medical Intensive Care Unit at Northwestern Memorial Hospital with suspected pneumonia, treated with cefepime. 9 of 70 (12.9%) required ECMO during treatment, predominantly veno-venous for severe respiratory failure. Baseline serum creatinine 1.12 +/- 0.66 mg/dL and Cockcroft-Gault creatinine clearance 115.8 +/- 88.6 mL/min (Table 1).",
    renal_function   = "Patients requiring concurrent hemodialysis, peritoneal dialysis or continuous renal replacement therapy were EXCLUDED, which the authors identify as the design feature that let them isolate an ECMO-specific effect uncontaminated by RRT. Renal function spans impairment to augmented renal clearance (CrCl mean 115.8, SD 88.6 mL/min).",
    dose_range       = "Institutional renal-function-based cefepime protocols; mean +/- SD initial 24 h dose 4.2 +/- 1.7 g/day. Monte Carlo simulations evaluated 2 g IV every 8 h as a 0.5 h intermittent infusion or a 4 h extended infusion, with and without a 2 g or 3 g loading dose given over 0.5 h.",
    regions          = "United States (single centre, Chicago, Illinois)",
    protein_binding  = "Not fitted. Total cefepime was quantified in plasma; unbound concentrations for the PK/PD target analysis were derived post hoc by applying a published protein binding of 20% (unbound fraction 0.80) (Valadez 2025 Methods, 'Individual PK/PD target attainment'). The model therefore outputs TOTAL plasma cefepime as Cc; multiply by 0.80 to obtain the free concentration the paper's fT>MIC targets are evaluated against.",
    sampling         = "Opportunistic residual plasma sampling, 1-14 samples per patient (median 2 in both the ECMO and non-ECMO groups). Measured concentrations 1.7-142.9 mg/L; mean +/- SD time after dose 7.8 +/- 6 h. Assay linear 1-100 mg/L by LC-MS/MS.",
    notes            = "Retrospective PK/PD study nested within the prospective SCRIPT study (Successful Clinical Response in Pneumonia Treatment), samples collected June 2018 - March 2024. Model estimated with the NPAG non-parametric algorithm in Pmetrics 2.1.1 for R 4.1.2. Final model diagnostics: population predictions R^2 = 0.611, slope 0.901, intercept 6.97, bias 0.213 mg/L, imprecision 8.42 mg^2/L^2; individual posterior predictions R^2 = 0.896, slope 1.02, intercept 1.05, bias -0.0736 mg/L, imprecision 0.989 mg^2/L^2."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All six values are the "Median (weighted)"
    # column of Valadez 2025 Table 2, which summarises the Pmetrics NPAG
    # non-parametric joint population distribution. The paired
    # "95% Credible interval" column is NOT carried into an omega -- see
    # the Inter-individual variability note below for why.
    # ---------------------------------------------------------------------

    lcl <- log(4.45)
    label("Total clearance at CRCL = 120 mL/min (L/h)")
    # Valadez 2025 Table 2: CL1 = 4.45 L/h (95% CrI 3.22-5.37). Table 2
    # footnote defines CL1 as "creatinine clearance normalized clearance".

    e_crcl_cl <- 0.76
    label("Power exponent on (CRCL / 120) for CL (unitless)")
    # Valadez 2025 Table 2: beta1 = 0.76 (95% CrI 0.52-0.92). Table 2
    # footnote: "non-linear scaling effect of clearance [i.e.,
    # (CrCl/120 mL/min)^beta1]". Estimated, not fixed: Table S1 run 11
    # shows that eliminating beta1 from the final model worsened fit by
    # dOFV > 3.84, which is why it was retained.

    lvc <- log(16.3)
    label("Central volume of distribution at WT = 70 kg without ECMO (L)")
    # Valadez 2025 Table 2: V1 = 16.3 L (95% CrI 6.54-33.5). Table 2
    # footnote: "weight-normalized (WT/70 kg) volume of distribution
    # (V1 in L)". Cross-checked against the Discussion, which quotes
    # "the median population parameter values for CL and Vd in non-ECMO
    # patients were 4.45 L/h and 16.3 L".

    e_ecmo_status_vc <- 1.043
    label("Exponential coefficient of ECMO_STATUS on central volume (unitless)")
    # Valadez 2025 Table 2: beta2 = 1.043 (95% CrI -0.40 to 1.59).
    # Enters as exp(beta2 * ECMO), NOT as (1 + beta2 * ECMO), per Eq 2
    # and Table S1 run 8 ("V1 * e(beta2*ECMO)"). Two independent
    # arithmetic checks in the source confirm the exponential form:
    #   exp(1.043)        = 2.837  -> Abstract / Conclusion "2.8-fold
    #                                 increase in Vd"
    #   16.3 * exp(1.043) = 46.25  -> Discussion "ECMO patients had a
    #                                 2.8-fold-greater Vd (i.e., 46.25 L
    #                                 vs 16.3 L)"
    # The linear reading 16.3 * (1 + 1.043) = 33.3 L matches neither.

    lk12 <- log(3.23)
    label("Central-to-peripheral rate constant KCP (1/h)")
    # Valadez 2025 Table 2: KCP = 3.23 1/h (95% CrI 0.49-8.61). Pmetrics
    # names the two-compartment micro-constants by their direction, so
    # KCP is central -> peripheral (k12) and KPC is peripheral ->
    # central (k21); Table S1 run 2 lists the base model's primary
    # variables as "CL, V, KPC, KCP".

    lk21 <- log(2.68)
    label("Peripheral-to-central rate constant KPC (1/h)")
    # Valadez 2025 Table 2: KPC = 2.68 1/h (95% CrI 2.01-9.12).

    # ---------------------------------------------------------------------
    # Inter-individual variability: NONE is encoded.
    #
    # NPAG estimates a discrete non-parametric joint density over support
    # points rather than a parametric omega matrix. Valadez 2025 Table 2
    # summarises that density with a weighted median and a "95% Credible
    # interval" only -- it reports no SD, no CV% and no omega, so there is
    # no published between-subject variance to carry. (Contrast
    # Tsai_2023_ceftriaxone.R, whose Table 2 does report per-parameter CV%
    # and which therefore CAN reconstruct a log-normal omega.)
    #
    # The Table 2 credible intervals are parameter PRECISION, not
    # between-subject spread, so converting them to omegas would fabricate
    # IIV. Three checks in the source agree on that reading:
    #   1. beta2's interval is -0.40 to 1.59, i.e. it straddles zero. beta2
    #      is a population covariate coefficient; a negative ECMO effect on
    #      Vd is not a subject-level phenotype, it is an imprecisely
    #      estimated fixed effect.
    #   2. The authors themselves read these spreads as precision: the
    #      ECMO-on-CL model was "rejected due to low precision (i.e.,
    #      CV > 200%) and point estimates near zero (i.e., a null effect)"
    #      (Results, "Model development and selection").
    #   3. The intervals are far too narrow to be between-subject spread.
    #      The CL1 interval implies a log-scale SD of about 0.13 (CV ~13%),
    #      whereas the Bayesian posterior CL of the non-ECMO patients spans
    #      an IQR of 2.0-4.1 L/h (Results, "Patient characteristics"),
    #      which implies a log-scale SD of about 0.53 (CV ~57%).
    #
    # The model is therefore a TYPICAL-VALUE model. `fixed(0)` etas are
    # deliberately NOT used: a zero-variance omega is singular and breaks
    # rxode2's Cholesky sampler. Users who need a stochastic simulation
    # must supply their own omega.
    #
    # Residual error: also NOT reported. Pmetrics parameterises residual
    # error as an assay-error polynomial SD(c) = C0 + C1*c + C2*c^2 + C3*c^3
    # scaled by an estimated gamma; the paper prints none of those
    # coefficients. The additive SD below is the ONLY residual statistic
    # the source publishes, carried as fixed() so it can never be mistaken
    # for an estimated $SIGMA. See the vignette Assumptions and deviations.
    # ---------------------------------------------------------------------

    addSd <- fixed(0.994)
    label("Additive residual error on Cc (mg/L), derived from the reported individual-prediction imprecision")
    # Valadez 2025 Results, "Model development and selection": the final
    # individual (Bayesian posterior) predictions had imprecision
    # 0.989 mg^2/L^2, which is Pmetrics' bias-adjusted mean squared
    # prediction error. sqrt(0.989) = 0.9945 mg/L. This is a derived
    # value, NOT a published residual-error parameter: the paper reports
    # no assay-error polynomial, no gamma, and no proportional term, so
    # the purely additive form is the minimum-assumption reading of the
    # only residual number available. The companion POPULATION-prediction
    # imprecision is 8.42 mg^2/L^2 (sqrt = 2.90 mg/L), which additionally
    # carries the unmodelled between-subject variability and is therefore
    # not the residual.
  })

  model({
    # Covariate reference values, exactly as printed in Valadez 2025
    # Eq 1 and Eq 2 and repeated in the Table 2 footnote.
    crcl_ref <- 120 # mL/min, raw Cockcroft-Gault
    wt_ref   <- 70  # kg, total body weight

    # Eq 1: CL = CL1 * (CrCL / 120)^beta1
    cl <- exp(lcl) * (CRCL / crcl_ref)^e_crcl_cl

    # Eq 2: Vd = V1 * (WT / 70) * e^(beta2 * ECMO)
    # The weight term is a plain ratio (exponent exactly 1), and the ECMO
    # term is exponential rather than proportional -- see the ini() note.
    vc <- exp(lvc) * (WT / wt_ref) * exp(e_ecmo_status_vc * ECMO_STATUS)

    # Distribution is parameterised as the micro-constants themselves, so
    # KCP and KPC are covariate-independent; the equivalent Q and Vp move
    # with vc (q = k12 * vc, vp = k12 * vc / k21).
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # TOTAL plasma cefepime. The paper's PK/PD targets are evaluated on
    # unbound drug, obtained by multiplying by the published unbound
    # fraction of 0.80 (protein binding 20%); that factor is a literature
    # constant applied outside the PK model, not a fitted parameter, so it
    # is recorded in population$protein_binding rather than in ini().
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
