Cerou_2024_amcenestrant_ameera12 <- function() {
  description <- paste(
    "Joint tumor growth inhibition (TGI) and progression-free survival (PFS)",
    "model for oral amcenestrant monotherapy in women with ER+/HER2- advanced",
    "breast cancer (Cerou 2024, development model fitted to AMEERA-1 and",
    "AMEERA-2, 75 patients, and externally evaluated on AMEERA-3). Same",
    "structure as Cerou_2024_amcenestrant (the final model recalibrated on",
    "AMEERA-1-2-3): resistant and treatment-sensitive tumor cells sharing one",
    "proliferation rate Kg = Ks / RE, Imax inhibition (Imax = 1) of",
    "sensitive-cell proliferation through two effect-delay transit",
    "compartments, amcenestrant PK driven by per-subject empirical-Bayes PK",
    "parameters supplied as data columns (the authors' population PK model is",
    "unpublished), and a Weibull PFS hazard multiplied by exp(betaL * dTS/dt)."
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
    MET_GE3 = list(
      description = "Three or more organs with metastasis at baseline (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Derived from the paper's NMET (number of organs with metastasis); the paper's indicator NMET < 3 is 1 - MET_GE3 (Table 3).",
      source_name = "NMET"
    ),
    LMET = list(
      description = "Liver metastasis at baseline (1 = present)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper LIVMET (Table 3). In this model it acts on TS0 only.",
      source_name = "LIVMET"
    ),
    LNMET = list(
      description = "Lymph-node metastasis at baseline (1 = present)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper LYMMET (Table S1 'Presence of lymph node metastasis at baseline'; Table 3).",
      source_name = "LYMMET"
    ),
    PRIOR_CDK46I = list(
      description = "Prior CDK4/6-inhibitor therapy (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Paper PRIORCDK4_6 (Table S1). Table 3 reports the effect for 'PRIOR CDK4/6 = no', so the model applies the coefficient to (1 - PRIOR_CDK46I).",
      source_name = "PRIORCDK4_6"
    ),
    ALBR = list(
      description = "Baseline serum albumin normalised to the laboratory upper limit of normal",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Paper ALBN. Enters as log(ALBR / 0.74). The paper does not print the",
        "learning-set median; 0.74 is the median in AMEERA-1 Part A and Part B",
        "(0.72 in AMEERA-2) in the authors' PAGE 2022 poster (abstract 10015)",
        "on the same 75 patients, which also reports a median albumin of 41 g/L."
      ),
      source_name = "ALBN"
    ),
    ALPR = list(
      description = "Baseline serum alkaline phosphatase normalised to the laboratory upper limit of normal",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Paper ALKN. Enters as log(ALPR / 0.89). The paper does not print the",
        "learning-set median; the authors' PAGE 2022 poster reports per-part",
        "medians 0.94 (AMEERA-1 Part A, n = 21), 0.89 (Part B, n = 46) and 0.59",
        "(AMEERA-2, n = 8), so the pooled median is close to 0.89. Do not",
        "substitute raw ALP: AMEERA-2 used a different upper limit (raw median",
        "202.5 IU/L but normalised median 0.59)."
      ),
      source_name = "ALKN"
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
    n_subjects = 75L,
    n_studies = 2L,
    age_range = "37-88 years",
    age_median = "59 (AMEERA-1 Part A), 64 (Part B), 66 (AMEERA-2) years",
    sex_female_pct = 100,
    disease_state = "ER+/HER2- advanced breast cancer, postmenopausal, heavily pretreated",
    dose_range = "amcenestrant 20-600 mg QD or 300 mg BID orally in 28-day cycles; 52/75 (70%) at 400 mg QD",
    regions = "International (AMEERA-1) and Japan (AMEERA-2)",
    liver_metastasis_pct = 61.3,
    notes = paste(
      "Learning dataset: AMEERA-1 Part A (n = 21), Part B (n = 46) and",
      "AMEERA-2 (n = 8), patients with measurable target lesions at baseline",
      "and at least one post-treatment measurement (Cerou 2024 Tables 1-2).",
      "The pooled liver-metastasis percentage (46/75) is computed from the",
      "per-part counts of Table 2. External evaluation: 98 patients of the",
      "AMEERA-3 amcenestrant arm (Figure 2)."
    )
  )

  ini({
    # Tumor size at baseline and its covariates (Table 3, AMEERA-1-2 column)
    lrbase <- log(37.67); label("Baseline target-lesion tumor size TS0 (mm)") # Table 3: TS0 = 37.67 mm
    e_nmetlt3_rbase <- -0.36; label("Effect of NMET < 3 (vs >= 3) on log TS0 (unitless)") # Table 3: beta TS0 NMET<3 = -0.36
    e_lmet_rbase <- 0.41; label("Effect of liver metastasis on log TS0 (unitless)") # Table 3: beta TS0 LIVMET yes = 0.41

    # Tumor cell kinetics
    lks <- fixed(log(0.09)); label("Tumor cell apoptosis rate Ks (1/day)") # Table 3: Ks = 0.09 (fixed)
    logitre <- logit(0.988); label("Logit of RE = Ks/Kg (unitless; RE in (0, 1])") # Table 3: RE = 0.988
    e_albr_re <- 5.77; label("Effect of log(ALBR/0.74) on logit RE (unitless)") # Table 3: beta RE logtALBN = 5.77
    e_lnmet_re <- -0.90; label("Effect of lymph-node metastasis on logit RE (unitless)") # Table 3: beta RE LYMMET yes = -0.90
    e_noprior_cdk46i_re <- 1.73; label("Effect of no prior CDK4/6 inhibitor on logit RE (unitless)") # Table 3: beta RE PRIOR CDK4/6=no = 1.73
    logitpr <- logit(0.714); label("Logit of baseline fraction of resistant cells pR (unitless)") # Table 3: pR = 71.4%

    # Amcenestrant effect on proliferation of sensitive cells
    limax <- fixed(log(1)); label("Maximum inhibition of sensitive-cell proliferation Imax (fraction)") # Methods: Imax was fixed to 1
    lic50 <- log(1.83); label("Effect-site concentration giving half-maximal inhibition IC50 (umol/L)") # Table 3: IC50 = 1.83 (printed mol/L; 1014 ng/mL per Results)
    lke0 <- fixed(log(0.005)); label("Effect-delay transit rate constant ke0 (1/day)") # Table 3: ke0 = 0.005 (fixed by grid search, Data S1)

    # PFS: Weibull baseline hazard and TS-slope link
    llam_haz <- log(165); label("Weibull scale Te of the PFS baseline hazard (day)") # Table 3: Te = 165 day
    lk_haz <- log(0.79); label("Weibull shape s of the PFS baseline hazard (unitless)") # Table 3: s = 0.79
    e_tsslope_haz <- 10.20; label("Log-hazard coefficient on the instantaneous TS slope dTS/dt (day/mm)") # Table 3: beta Slope = 10.20
    e_alpr_haz <- 0.93; label("Effect of log(ALPR/0.89) on log PFS hazard (unitless)") # Table 3: beta E logtALKN = 0.93

    # IIV: Table 3 prints omega SD x 100 (Table S2 lists the same quantities as SDs);
    # variance = (omega / 100)^2
    etalrbase ~ 0.3249 # Table 3: omega TS0 = 57 -> 0.57^2
    etalogitre ~ 0.7744 # Table 3: omega RE = 88 -> 0.88^2 (logit-normal)
    etalogitpr ~ 1.96 # Table 3: omega PR = 140 -> 1.40^2
    etalke0 ~ 6.76 # Table 3: omega ke0 = 260 -> 2.60^2

    # Residual error on TS (combined additive + proportional)
    addSd <- 0.53; label("Additive residual error on TS (mm)") # Table 3: sigma additive = 0.53 mm
    propSd <- 0.0822; label("Proportional residual error on TS (fraction)") # Table 3: sigma proportional = 8.22%
  })
  model({
    # Covariate terms (Methods: continuous covariates enter as log(cov / median))
    nmetlt3 <- 1 - MET_GE3
    logt_albr <- log(ALBR / 0.74)
    logt_alpr <- log(ALPR / 0.89)

    # Individual parameters
    ts0 <- exp(lrbase + e_nmetlt3_rbase * nmetlt3 + e_lmet_rbase * LMET + etalrbase)
    ks <- exp(lks)
    re <- expit(
      logitre + e_albr_re * logt_albr + e_lnmet_re * LNMET +
        e_noprior_cdk46i_re * (1 - PRIOR_CDK46I) + etalogitre
    )
    pr <- expit(logitpr + etalogitpr)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    ke0 <- exp(lke0 + etalke0)
    lam_haz <- exp(llam_haz)
    k_haz <- exp(lk_haz)
    e_cov_haz <- exp(e_alpr_haz * logt_alpr)

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
