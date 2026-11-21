Cerou_2024_amcenestrant <- function() {
  description <- paste(
    "Joint tumor growth inhibition (TGI) and progression-free survival (PFS)",
    "model for oral amcenestrant monotherapy in women with ER+/HER2- advanced",
    "breast cancer (Cerou 2024, final model calibrated on AMEERA-1-2-3,",
    "214 patients). Target-lesion tumor size (sum of longest diameters, mm) is",
    "the sum of resistant and treatment-sensitive cell populations sharing one",
    "proliferation rate Kg = Ks / RE; amcenestrant inhibits proliferation of",
    "the sensitive cells only (Imax model, Imax = 1) through two effect-delay",
    "transit compartments fed by the plasma concentration. Amcenestrant PK",
    "(depot plus three transit states into a two-compartment linear model) is",
    "driven by per-subject empirical-Bayes PK parameters supplied as data",
    "columns (KTR_INDIV, CL_INDIV, VC_INDIV, Q_INDIV, VP_INDIV), exactly as",
    "the authors' Monolix regressors; the authors' population PK model is",
    "unpublished. PFS follows a Weibull hazard multiplied by",
    "exp(betaL * dTS/dt)."
  )
  reference <- paste(
    "Cerou M, Thai HT, Deyme L, Fliscounakis-Huynh S, Cartot-Cotton S,",
    "Veyrat-Follet C. Joint modeling of tumor dynamics and progression-free",
    "survival in advanced breast cancer: Leveraging data from amcenestrant",
    "early phase I-II trials. CPT Pharmacometrics Syst Pharmacol.",
    "2024;13(6):941-953. doi:10.1002/psp4.13128"
  )
  vignette <- "Cerou_2024_amcenestrant"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L (plasma Cc); umol/L (effect-site Ce, IC50); mm (tumor size TS)"
  )
  paper_specific_compartments <- c("tumor_resistant", "tumor_sensitive")

  covariateData <- list(
    KTR_INDIV = list(
      description = "Per-subject empirical-Bayes amcenestrant transit / absorption rate constant (shared by the depot and the three transit states)",
      units = "1/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Monolix regressor 'KTR' (Data S1 model file). Individual estimate from",
        "the authors' amcenestrant population PK model (Cerou 2022 PAGE abstract",
        "10015), which is not published with parameter values. For simulation,",
        "supply each subject's value; the vignette uses an illustrative typical PK",
        "built by the maintainers from published noncompartmental results."
      ),
      source_name = "KTR"
    ),
    CL_INDIV = list(
      description = "Per-subject empirical-Bayes apparent amcenestrant clearance (CL/F)",
      units = "L/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Monolix regressor 'CL' (Data S1). Same upstream source and simulation guidance as KTR_INDIV. Units are L/day because the model time unit is day.",
      source_name = "CL"
    ),
    VC_INDIV = list(
      description = "Per-subject empirical-Bayes apparent amcenestrant central volume (V1/F)",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = "Monolix regressor 'V1' (Data S1). Same upstream source and simulation guidance as KTR_INDIV.",
      source_name = "V1"
    ),
    Q_INDIV = list(
      description = "Per-subject empirical-Bayes apparent amcenestrant intercompartmental clearance (Q/F)",
      units = "L/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Monolix regressor 'Q' (Data S1). Same upstream source and simulation guidance as KTR_INDIV.",
      source_name = "Q"
    ),
    VP_INDIV = list(
      description = "Per-subject empirical-Bayes apparent amcenestrant peripheral volume (V2/F)",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = "Monolix regressor 'V2' (Data S1). Same upstream source and simulation guidance as KTR_INDIV.",
      source_name = "V2"
    ),
    MET_GE2 = list(
      description = "Two or more organs with metastasis at baseline (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Derived from the paper's categorical NMET (number of organs with",
        "metastasis; categories 1, 2 and >= 3, reference >= 3; Table 3).",
        "The model forms the paper's indicators as NMET = 1 <=> 1 - MET_GE2 and",
        "NMET = 2 <=> MET_GE2 - MET_GE3."
      ),
      source_name = "NMET"
    ),
    MET_GE3 = list(
      description = "Three or more organs with metastasis at baseline (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Derived from NMET (Table 2 row 'Number of organs with metastasis >= 3': 45.3% of the calibration set). See MET_GE2.",
      source_name = "NMET"
    ),
    LMET = list(
      description = "Liver metastasis at baseline (1 = present)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper LIVMET; 47.2% of the calibration set (Table 2). Acts on TS0, on RE and on the PFS baseline hazard.",
      source_name = "LIVMET"
    ),
    PREMENOPAUSAL = list(
      description = "Pre- or peri-menopausal status at baseline (1 = pre/peri-menopausal, 0 = postmenopausal)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper MENOPAUSAL (pre/peri vs post). AMEERA-1 and AMEERA-2 enrolled postmenopausal women only; AMEERA-3 also enrolled premenopausal women on a GnRH analogue (Data S1).",
      source_name = "MENOPAUSAL"
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper enters AST normalised to the upper limit of normal (ASTN) as",
        "log(ASTN / median ASTN). With a common upper limit this equals",
        "log(AST / median AST); the model centres raw AST on the calibration-set",
        "median 28.4 IU/L (Table 2) because the median of ASTN is not reported.",
        "The common-limit assumption is reasonable for AST but not for every",
        "laboratory: the authors' 2022 PAGE poster shows AMEERA-2 used a different",
        "alkaline-phosphatase upper limit, so users with mixed-laboratory data",
        "should convert AST to a common reference range first."
      ),
      source_name = "ASTN"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as log(AGE / 60); 60 years is the calibration-set median (Table 2).",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "amcenestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "amcenestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "amcenestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "amcenestrant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amcenestrant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amcenestrant", units = "mg", specimen = "plasma", verified = TRUE),
    effect1 = list(analyte = "amcenestrant", units = "umol/L", specimen = "not applicable", verified = TRUE),
    effect2 = list(analyte = "amcenestrant", units = "umol/L", specimen = "not applicable", verified = TRUE),
    tumor_resistant = list(
      analyte = "resistant tumor cells (target-lesion diameter)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_sensitive = list(
      analyte = "treatment-sensitive tumor cells (target-lesion diameter)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "progression-free-survival hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 214L,
    n_studies = 3L,
    age_range = "29-88 years",
    age_median = "60 years",
    sex_female_pct = 100,
    disease_state = "ER+/HER2- advanced (locally advanced or metastatic) breast cancer",
    dose_range = "amcenestrant 20-600 mg QD or 300 mg BID orally in 28-day cycles; 194/214 (91%) at 400 mg QD",
    regions = "International (AMEERA-1 and AMEERA-3 non-Japanese; AMEERA-2 Japan)",
    n_measurable_ts = 178L,
    n_nonmeasurable_ts = 36L,
    liver_metastasis_pct = 47.2,
    nmet_ge3_pct = 45.3,
    prior_cdk46_pct = 74.8,
    notes = paste(
      "Calibration dataset AMEERA-1 (n = 68), AMEERA-2 (Japanese, n = 8) and",
      "the amcenestrant arm of AMEERA-3 (n = 138); 178 patients had measurable",
      "target lesions and 36 had non-measurable disease only (PFS data only).",
      "Baseline covariates: Cerou 2024 Table 2 'All (n = 214)' column; study",
      "design: Table 1. Monolix 2020R1 SAEM."
    )
  )

  ini({
    # Tumor size at baseline and its covariates (Table 3, AMEERA-1-2-3 column)
    lrbase <- log(43.17); label("Baseline target-lesion tumor size TS0 (mm)") # Table 3: TS0 = 43.17 mm
    e_nmet1_rbase <- -0.38; label("Effect of NMET = 1 (vs >= 3) on log TS0 (unitless)") # Table 3: beta TS0 NMET=1 = -0.38
    e_nmet2_rbase <- -0.40; label("Effect of NMET = 2 (vs >= 3) on log TS0 (unitless)") # Table 3: beta TS0 NMET=2 = -0.40
    e_lmet_rbase <- 0.25; label("Effect of liver metastasis on log TS0 (unitless)") # Table 3: beta TS0 LIVMET yes = 0.25
    e_premenopausal_rbase <- 0.34; label("Effect of pre/peri-menopausal status on log TS0 (unitless)") # Table 3: beta TS0 MENOPAUSAL pre/peri = 0.34

    # Tumor cell kinetics
    lks <- fixed(log(0.09)); label("Tumor cell apoptosis rate Ks (1/day)") # Table 3: Ks = 0.09 (fixed)
    logitre <- logit(0.992); label("Logit of RE = Ks/Kg (unitless; RE in (0, 1])") # Table 3: RE = 0.992
    e_ast_re <- -0.57; label("Effect of log(AST/median) on logit RE (unitless)") # Table 3: beta RE logtASTN = -0.57
    e_lmet_re <- -0.78; label("Effect of liver metastasis on logit RE (unitless)") # Table 3: beta RE LIVMET yes = -0.78
    logitpr <- logit(0.660); label("Logit of baseline fraction of resistant cells pR (unitless)") # Table 3: pR = 66.0%
    e_age_pr <- 1.35; label("Effect of log(AGE/60) on logit pR (unitless)") # Table 3: beta pR logtAGE = 1.35

    # Amcenestrant effect on proliferation of sensitive cells
    limax <- fixed(log(1)); label("Maximum inhibition of sensitive-cell proliferation Imax (fraction)") # Methods: Imax was fixed to 1
    lic50 <- log(3.74); label("Effect-site concentration giving half-maximal inhibition IC50 (umol/L)") # Table 3: IC50 = 3.74 (printed mol/L; 3.74 x 554.48 = 2073 ng/mL per Results)
    lke0 <- fixed(log(0.005)); label("Effect-delay transit rate constant ke0 (1/day)") # Table 3: ke0 = 0.005 (fixed by grid search, Data S1)

    # PFS: Weibull baseline hazard and TS-slope link
    llam_haz <- log(297.73); label("Weibull scale Te of the PFS baseline hazard (day)") # Table 3: Te = 297.73 day
    lk_haz <- log(0.84); label("Weibull shape s of the PFS baseline hazard (unitless)") # Table 3: s = 0.84
    e_tsslope_haz <- 10.73; label("Log-hazard coefficient on the instantaneous TS slope dTS/dt (day/mm)") # Table 3: beta Slope = 10.73
    e_ast_haz <- 0.68; label("Effect of log(AST/median) on log PFS hazard (unitless)") # Table 3: beta E logtASTN = 0.68
    e_lmet_haz <- 0.52; label("Effect of liver metastasis on log PFS hazard (unitless)") # Table 3: beta E LIVMET yes = 0.52

    # IIV: Table 3 prints omega SD x 100 (Table S2 lists the same quantities as SDs);
    # variance = (omega / 100)^2
    etalrbase ~ 0.3249 # Table 3: omega TS0 = 57 -> 0.57^2
    etalogitre ~ 1.5625 # Table 3: omega RE = 125 -> 1.25^2 (logit-normal)
    etalogitpr ~ 0.7056 # Table 3: omega PR = 84 -> 0.84^2
    etalke0 ~ 12.3201 # Table 3: omega ke0 = 351 -> 3.51^2

    # Residual error on TS (combined additive + proportional)
    addSd <- 0.75; label("Additive residual error on TS (mm)") # Table 3: sigma additive = 0.75 mm
    propSd <- 0.06; label("Proportional residual error on TS (fraction)") # Table 3: sigma proportional = 6%
  })
  model({
    # Covariate terms (Methods: continuous covariates enter as log(cov / median))
    nmet1 <- 1 - MET_GE2
    nmet2 <- MET_GE2 - MET_GE3
    logt_ast <- log(AST / 28.4)
    logt_age <- log(AGE / 60)

    # Individual parameters
    ts0 <- exp(
      lrbase + e_nmet1_rbase * nmet1 + e_nmet2_rbase * nmet2 +
        e_lmet_rbase * LMET + e_premenopausal_rbase * PREMENOPAUSAL + etalrbase
    )
    ks <- exp(lks)
    re <- expit(logitre + e_ast_re * logt_ast + e_lmet_re * LMET + etalogitre)
    pr <- expit(logitpr + e_age_pr * logt_age + etalogitpr)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    ke0 <- exp(lke0 + etalke0)
    lam_haz <- exp(llam_haz)
    k_haz <- exp(lk_haz)
    e_cov_haz <- exp(e_ast_haz * logt_ast + e_lmet_haz * LMET)

    # Amcenestrant PK from per-subject inputs (Data S1 Monolix model)
    ktr <- KTR_INDIV
    kel <- CL_INDIV / VC_INDIV
    k12 <- Q_INDIV / VC_INDIV
    k21 <- Q_INDIV / VP_INDIV

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / VC_INDIV # mg/L
    Cmolar <- Cc / 554.48 * 1000 # umol/L; MW 554.48 g/mol (Data S1)

    # Two-compartment effect-delay chain; the second state drives the effect
    d/dt(effect1) <- (Cmolar - effect1) * ke0
    d/dt(effect2) <- (effect1 - effect2) * ke0
    inhibition <- imax * effect2 / (effect2 + ic50)

    # Resistant / sensitive tumor cells (sum of longest diameters, mm)
    kg <- ks / re
    tumor_resistant(0) <- ts0 * pr
    tumor_sensitive(0) <- ts0 * (1 - pr)
    dres <- (kg - ks) * tumor_resistant
    dsens <- (kg * (1 - inhibition) - ks) * tumor_sensitive
    d/dt(tumor_resistant) <- dres
    d/dt(tumor_sensitive) <- dsens
    TS <- tumor_resistant + tumor_sensitive
    ts_slope <- dres + dsens # mm/day

    # PFS hazard: Weibull x covariates x exp(betaL * dTS/dt). Time 0 is the
    # start of treatment (Monolix Tstart). The 1e-6 day offset keeps the
    # s < 1 Weibull term finite at t = 0; the link exponent is capped at 50
    # because an ODE keeps integrating past progression, and exp() of a fast
    # regrowth slope overflows while survival is already zero.
    link <- min(e_tsslope_haz * ts_slope, 50)
    hazard <- (k_haz / lam_haz) * ((t + 1e-6) / lam_haz)^(k_haz - 1) * e_cov_haz * exp(link)
    d/dt(cumhaz) <- hazard
    surv <- exp(-cumhaz)

    TS ~ add(addSd) + prop(propSd)
  })
}
