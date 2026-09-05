Sokolov_2023_dapagliflozin <- function() {
  description <- paste(
    "Semi-mechanistic exposure-response (PD-only) model for dapagliflozin as",
    "an insulin adjunct in adults with type 1 diabetes (T1D). Per-subject",
    "steady-state 24 h dapagliflozin AUC (AUC_DAPA, supplied as a covariate",
    "column from an upstream population PK analysis) drives a three-stage",
    "cascade of ratio-to-pre-treatment-baseline endpoints: (1) an Imax",
    "reduction in total daily basal insulin dose (rins), (2) average daily",
    "plasma glucose by continuous glucose monitoring (rglu) as a power",
    "function of rins times a second Imax term plus a linear",
    "treatment-independent glucose drift, and (3) HbA1c (rhba1c) as a power",
    "function of rglu times a third, direct Imax term capturing the",
    "glucose-independent dapagliflozin benefit. Absolute HbA1c is recovered",
    "from the per-subject baseline covariate HBA1C. No ODEs and no dosing",
    "events; model time is in weeks."
  )
  reference <- paste(
    "Sokolov V, Yakovleva T, Penland RC, Boulton DW, Tang W.",
    "Effectiveness of dapagliflozin as an insulin adjunct in type 1",
    "diabetes: a semi-mechanistic exposure-response model.",
    "Front Pharmacol. 2023;14:1229255. doi:10.3389/fphar.2023.1229255."
  )
  vignette <- "Sokolov_2023_dapagliflozin"
  units <- list(
    time          = "week",
    dosing        = "n/a (no drug dosing events; dapagliflozin exposure enters through the per-subject AUC_DAPA covariate)",
    concentration = "n/a (multi-output PD-only model; rins, rglu and rhba1c are unitless ratios to the pre-treatment baseline, and hba1c is in % NGSP)",
    AUC_DAPA      = "ng*h/mL"
  )
  covariateData <- list(
    AUC_DAPA = list(
      description        = "Steady-state 24 h dapagliflozin AUC supplied as a per-subject (time-fixed) drug-exposure covariate.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject steady-state 24 h dapagliflozin exposure. Sokolov 2023",
        "Table 1 states the values were taken from a previously performed",
        "population PK analysis (Melin et al. 2022) and supplied to the",
        "exposure-response model as a regressor; the PD model contains no PK",
        "sub-model. AUC_DAPA = 0 recovers the placebo arm exactly (all three",
        "Imax terms vanish, so rins = 1 and the only remaining dynamics are",
        "the linear glucose drift and its propagation to HbA1c). Mean AUC by",
        "dose reported in the Figure 10 caption: 51.4 (1 mg), 130.6",
        "(2.5 mg), 294.5 (5 mg) and 594.3 (10 mg) ng*h/mL once daily."
      ),
      source_name        = "AUC"
    ),
    HBA1C = list(
      description        = "Baseline (pre-treatment) HbA1c, per-subject and time-fixed",
      units              = "% (NGSP)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to convert the fitted ratio-to-baseline HbA1c (rhba1c)",
        "into an absolute HbA1c trajectory: hba1c = HBA1C * rhba1c. The",
        "Sokolov 2023 model itself was fitted entirely on the ratio scale",
        "(Supplementary Equations 1-3 contain no baseline terms), so HBA1C",
        "does not enter any estimated relationship. The forward simulations",
        "in Figure 10 use the pooled mean baseline HbA1c = 8.48 %; Table 4",
        "reports study medians of 8.4, 8.4 and 8.3 %."
      ),
      source_name        = "HbA1c"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 1661L,
    n_subjects_estimation = 883L,
    n_subjects_validation = 778L,
    n_measurements       = 12460L,
    n_studies            = 3L,
    studies              = paste(
      "Estimation: NCT01498185 (phase 2 dose-ranging, N = 70, first 7 inpatient",
      "days only) pooled with NCT02460978 (DEPICT-2, phase 3, 24-week",
      "double-blind period). External validation: NCT02268214 (DEPICT-1,",
      "phase 3, 24-week double-blind period), not used in model development."
    ),
    age_range            = "18-75 years (study medians 30, 43 and 43 years; Table 4)",
    weight_range         = "44.6-184.8 kg (study medians 74.8, 80.8 and 76.8 kg; Table 4)",
    bmi_range            = "18.2-65.8 kg/m^2 (study medians 23.9, 27.8 and 26.9 kg/m^2; Table 4)",
    sex_female_pct       = c(NCT01498185 = 42.9, NCT02268214 = 52.1, NCT02460978 = 56.0),
    race_ethnicity       = paste(
      "Predominantly White (88.6 %, 95.6 % and 78.4 % by study); NCT02460978",
      "enrolled 19.7 % Asian patients (Table 4)."
    ),
    disease_state        = paste(
      "Adults with inadequately controlled type 1 diabetes on background",
      "basal-bolus insulin (multiple daily injections or continuous",
      "subcutaneous insulin infusion). Phase 3 randomisation required HbA1c",
      "7.5-10.5 % (58-91 mmol/mol). Baseline medians: HbA1c 8.3-8.4 %,",
      "24 h mean CGM glucose 170-190 mg/dL, total daily insulin dose",
      "48-54 units, eGFR 89-91 mL/min/1.73 m^2, diabetes duration 17-19 years."
    ),
    dose_range           = "Dapagliflozin 1, 2.5, 5 or 10 mg once daily, or placebo; exposure enters through AUC_DAPA.",
    regions              = "Multi-national (NCT01498185, DEPICT-1 and DEPICT-2 trial programmes).",
    notes                = paste(
      "Non-linear mixed-effects model fitted in Monolix 2020R1 using a",
      "three-step sequential strategy (Supplementary 'Structural model' and",
      "'Model development'). Step 1 fits Equation 1 (basal insulin dose ratio",
      "vs AUC) on the active arms only; step 2 fits Equation 2 (glucose ratio)",
      "with the observed rINS and AUC as regressors; step 3 fits Equation 3",
      "(HbA1c ratio) with the observed rGLU and AUC as regressors. Each step's",
      "estimates were fixed before the next. This file encodes the composed",
      "forward-simulation form used for Figures 7, 9 and 10, in which",
      "Equation 1 feeds Equation 2 and Equation 2 feeds Equation 3."
    )
  )

  ini({
    # ---- Step 1: basal daily insulin dose ratio (Supplementary Equation 1; Table 3 step 1) ----
    # Imax_ins carries logit-scale IIV (Supplementary Equation 6), so the
    # typical value is stored on the logit scale here.
    logitimax_ins  <- qlogis(0.0941); label("Logit of Imax_ins -- maximum fractional dapagliflozin-mediated reduction in total daily basal insulin dose (unitless)")  # Table 3 step 1: Imax_ins = 0.0941, RSE 8.14 %; qlogis(0.0941) = -2.2648
    lauc50_ins     <- log(38.8);      label("IAUC50_ins -- dapagliflozin AUC giving half the maximal basal insulin dose reduction (ng*h/mL; log-scale)")            # Table 3 step 1: IAUC50_ins = 38.8 ng*h/mL, RSE 35.2 %

    # ---- Step 2: average daily CGM glucose ratio (Supplementary Equation 2; Table 3 step 2) ----
    # SIGN NOTE (see the vignette Assumptions and deviations section): Table 3
    # prints k1 = 0.0674 without a sign, but Supplementary Equation 2 applies
    # it as the bare exponent (rINS)^k1. A positive exponent makes glucose FALL
    # when insulin is withdrawn, which contradicts the paper's own reported
    # results. The negative value is required to reproduce, simultaneously:
    #   - the Abstract / Figure 10A statement that a 50 % basal insulin
    #     reduction gives "~5 % increase in glucose exposure"
    #     (0.5^-0.0674 = 1.0478, i.e. +4.8 %);
    #   - the Figure 10A pair -0.50 % (no insulin adjustment) vs -0.42 %
    #     (50 % basal insulin reduction) at 10 mg;
    #   - the negative slope (-0.113) of the observed glucose-vs-basal-insulin
    #     regression in Figure 4A.
    # The magnitude 0.0674 is used exactly as printed. The sibling M-EASE-1
    # empagliflozin T1D model (Johnston_2021_empagliflozin_MEASE1) likewise
    # carries a negative insulin-to-glucose power exponent (-0.261).
    ins_glu_eff    <- -0.0674;        label("k1 -- power exponent of the basal insulin dose ratio on the glucose ratio (unitless)")                                  # Table 3 step 2: k1 = 0.0674, RSE 29.7 %; sign recovered from the paper's own simulation results (see note above)
    pbo_glu_rate   <- 0.0015;         label("keff -- linear treatment-independent drift in the glucose ratio (per week)")                                            # Table 3 step 2: keff = 0.0015, RSE 21.5 %. Table 3 labels the units "mg/dL/week", but Equation 2 adds this term to a unitless ratio, so it is a fraction-of-baseline per week: 0.0015 * 24 = +3.6 % at week 24, matching the "2 %-4 % increase in glucose at weeks 12 and 24" in the Discussion.
    imax_glu       <- 0.15;           label("Imax_glu -- maximum fractional dapagliflozin-mediated reduction in average daily glucose (unitless)")                   # Table 3 step 2: Imax_glu = 0.15, RSE 7.35 %
    lauc50_glu     <- log(67.4);      label("IAUC50_glu -- dapagliflozin AUC giving half the maximal glucose reduction (ng*h/mL; log-scale)")                        # Table 3 step 2: IAUC50_glu = 67.4 ng*h/mL, RSE 28.2 %

    # ---- Step 3: HbA1c ratio (Supplementary Equation 3; Table 3 step 3) ----
    lgamma_glueff  <- log(0.165);     label("k2 -- power exponent of the glucose ratio on the HbA1c ratio (unitless; log-scale)")                                    # Table 3 step 3: k2 = 0.165, RSE 8.03 %; log-scale because the IIV is log-transformed (Supplementary Equation 5)
    imax_hba1c     <- 0.0421;         label("Imax_hba1c -- maximum fractional direct (glucose-independent) dapagliflozin-mediated reduction in HbA1c (unitless)")    # Table 3 step 3: Imax_hba1c = 0.0421, RSE 6.69 %
    lauc50_hba1c   <- fixed(log(67.4)); label("IAUC50_hba1c -- dapagliflozin AUC giving half the maximal direct HbA1c reduction (ng*h/mL; log-scale)")               # Table 3 step 3: "same as IAUC50_glu", no RSE reported. Supplementary "Model development": IAUC50_hba1c "was not identifiable (RSE > 50%) and had to be fixed at the previously defined value of IAUC50_glu = 67.4 ng/mL*h".

    # ---- Inter-individual variability (Table 3 parenthesised omega column) ----
    # Monolix reports omega as the standard deviation of the random effect;
    # nlmixr2 `ini()` takes the variance, so each entry below is omega^2.
    # Only three of the nine structural parameters carry IIV, matching the
    # Supplementary "Model development" narrative (one random effect per step).
    etalogitimax_ins ~ 1.1881    # omega_Imax_ins = 1.09 (RSE 6.43 %) on the logit scale -> variance 1.09^2
    etapbo_glu_rate  ~ 9e-06     # omega_keff = 0.003 (RSE 9.02 %), additive on the untransformed parameter (Supplementary Equation 4) -> variance 0.003^2
    etalgamma_glueff ~ 0.546121  # omega_k2 = 0.739 (RSE 8.08 %) on the log scale -> variance 0.739^2

    # ---- Residual error (Table 3; one error model per modelling step) ----
    propSd_rins   <- 0.185; label("Proportional residual error on the basal insulin dose ratio (fraction)")  # Table 3 step 1: b_ins = 0.185, RSE 1.93 %
    propSd_rglu   <- 0.162; label("Proportional residual error on the glucose ratio (fraction)")             # Table 3 step 2: b_glu = 0.162, RSE 1.99 %
    addSd_rhba1c  <- 0.061; label("Constant residual error on the HbA1c ratio (fraction of baseline)")       # Table 3 step 3: a_hba1c = 0.061, RSE 2.15 %; constant on the ratio scale, i.e. about 0.5 % HbA1c at a baseline of 8.48 %
  })

  model({
    # ---- Step 1 -- Supplementary Equation 1 ----
    #   rINS_ik = 1 - Imax_ins,i * AUC_i / (IAUC50_ins + AUC_i)
    # Imax_ins is logit-normally distributed (Supplementary Equation 6), which
    # keeps every individual's maximum insulin-dose reduction inside [0, 1].
    imax_ins  <- expit(logitimax_ins + etalogitimax_ins)
    auc50_ins <- exp(lauc50_ins)
    rins      <- 1 - imax_ins * AUC_DAPA / (auc50_ins + AUC_DAPA)

    # ---- Step 2 -- Supplementary Equation 2 ----
    #   rGLU_ik = (rINS_ik)^k1 * (1 - Imax_glu * AUC_i / (IAUC50_glu + AUC_i))
    #             + keff_i * TIME_k
    # keff carries additive (untransformed) IIV per Supplementary Equation 4.
    # `time` is model time in WEEKS (see the file `units` field): TIME_k in
    # Equation 2 is "the time of the kth visit (weeks)".
    keff_i    <- pbo_glu_rate + etapbo_glu_rate
    auc50_glu <- exp(lauc50_glu)
    rglu      <- rins^ins_glu_eff *
                 (1 - imax_glu * AUC_DAPA / (auc50_glu + AUC_DAPA)) +
                 keff_i * time

    # ---- Step 3 -- Supplementary Equation 3 ----
    #   rHBA1C_ik = (rGLU_ik)^k2_i
    #               * (1 - Imax_hba1c * AUC_i / (IAUC50_hba1c + AUC_i))
    # The second factor is the additional direct dapagliflozin effect on HbA1c
    # that the glucose pathway alone could not explain (Results, Figure 7C).
    k2_i        <- exp(lgamma_glueff + etalgamma_glueff)
    auc50_hba1c <- exp(lauc50_hba1c)
    rhba1c      <- rglu^k2_i *
                   (1 - imax_hba1c * AUC_DAPA / (auc50_hba1c + AUC_DAPA))

    # ---- Absolute HbA1c ----
    # The model was fitted on the ratio scale; the absolute trajectory is
    # recovered from the per-subject pre-treatment baseline. Figure 10 uses the
    # pooled mean baseline HbA1c = 8.48 %.
    hba1c <- HBA1C * rhba1c

    # ---- Multi-output residual error (Table 3) ----
    rins   ~ prop(propSd_rins)
    rglu   ~ prop(propSd_rglu)
    rhba1c ~ add(addSd_rhba1c)
  })
}
