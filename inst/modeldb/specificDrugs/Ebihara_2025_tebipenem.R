Ebihara_2025_tebipenem <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "tebipenem (given as the prodrug tebipenem pivoxil) in Japanese adults,",
    "stratified by renal function (Ebihara 2025). Apparent clearance CL/F,",
    "apparent volume of distribution Vd/F, and the absorption rate constant ka",
    "each take a separate value in each of four Cockcroft-Gault creatinine",
    "clearance strata (CRCL >= 80, 50 to 80, 30 to 50 and < 30 mL/min); the",
    "stratum is selected inside model() from the continuous CRCL covariate. The",
    "stratum-specific values are the cohort means of a 17-subject Japanese phase",
    "1 study in healthy volunteers and subjects with reduced renal function",
    "(Nakashima 2009), and the inter-individual variability is borrowed from a",
    "faropenem population PK analysis (Hirooka 2005) because no tebipenem",
    "population PK model existed for Japanese adults. Every parameter is",
    "therefore fixed: the source estimates nothing and performs Monte Carlo",
    "simulation only, so propSd is also fixed at 0. The plasma unbound fraction",
    "fu = 0.33 is carried as a parameter and the model reports unbound",
    "concentration Cu alongside total Cc, because the source's target attainment",
    "criterion is on unbound drug: fAUC(0-24)/MIC/tau >= 34.58 for urinary tract",
    "infections caused by ESBL-producing Enterobacterales."
  )
  reference <- paste(
    "Ebihara F, Maruyama T, Kasai H, Shiokawa M, Matsunaga N, Hamada Y (2025).",
    "Preclinical Study of Pharmacokinetic/Pharmacodynamic Analysis of Tebipenem",
    "Using Monte Carlo Simulation for Extended-Spectrum beta-Lactamase-Producing",
    "Bacterial Urinary Tract Infections in Japanese Patients According to Renal",
    "Function. Antibiotics 14(7):648. doi:10.3390/antibiotics14070648.",
    "The structural parameter values originate in Nakashima M, Morita J,",
    "Aizawa K (2009) Jpn J Chemother 57:109-114 and the inter-individual",
    "variability in Hirooka H, Sou M, Miyata K, Shiba K, Tanigawara Y (2005)",
    "Jpn J Clin Pharmacol Ther 36:197-207 (faropenem); both are reproduced in",
    "Ebihara 2025, which is the source-trace of record for every value here.",
    sep = " "
  )
  vignette <- "Ebihara_2025_tebipenem"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: Ebihara 2025 Methods 4.1 states a
  # one-compartment model with first-order absorption for orally administered
  # tebipenem, and Table 3 reports the apparent (oral) parameters CL/F and
  # Vd/F, so the central state is plasma tebipenem and the depot is the
  # absorption site.
  compartmentData <- list(
    depot   = list(analyte = "tebipenem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tebipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance, reported by the source as raw (NOT",
        "BSA-normalized) mL/min. Ebihara 2025 abbreviates it CCR and uses it",
        "only to allocate a subject to one of four renal-function strata; no",
        "continuous within-stratum relationship is estimated."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Stratum boundaries are the Ebihara 2025 Table 3 column headers and",
        "Methods 4.1: 'Renal function was evaluated in the following four",
        "groups: CCR >= 80 mL/min, 50 <= CCR < 80 mL/min, 30 <= CCR < 50",
        "mL/min, and CCR < 30 mL/min.' The four strata are exhaustive and",
        "mutually exclusive, so model() selects exactly one value of CL/F,",
        "Vd/F and ka per subject; the model is therefore a step function of",
        "CRCL and is piecewise-constant within a stratum. Values are raw",
        "Cockcroft-Gault mL/min and are NOT BSA-normalized, following the",
        "Delattre 2010 amikacin, Georges 2009 ceftazidime, Chen 2023",
        "nemonoxacin, Wada 2023 sparsentan and Shu 2024 posaconazole",
        "precedents in this register; the source underlying cohort mean CCR",
        "per stratum is 106.3 mL/min (n = 6), 65.8 mL/min (n = 6), 40.2",
        "mL/min (n = 2) and 9.1 mL/min (n = 3) from highest to lowest,",
        "reported in Ebihara 2025 Methods 4.4 'CCR simulation approach'.",
        "Binary stratification of CRCL at a threshold has precedent in",
        "NA_NA_lidocaine.R; this model generalises it to four bands."
      ),
      source_name        = "CCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 17,
    n_studies      = 1,
    weight_median  = "61.9 kg (cohort mean)",
    height_median  = "164.3 cm (cohort mean)",
    bsa_median     = "1.67 m^2 (cohort mean)",
    regions        = "Japan",
    disease_state  = paste(
      "Underlying PK cohort: healthy Japanese volunteers plus subjects with",
      "reduced renal function (Nakashima 2009 phase 1). Simulation target",
      "population: Japanese adults with urinary tract infection caused by",
      "extended-spectrum beta-lactamase (ESBL)-producing Enterobacterales",
      "(Escherichia coli and Klebsiella pneumoniae, MIC90 0.03 mg/L)."
    ),
    renal_function = paste(
      "Four Cockcroft-Gault creatinine clearance strata with cohort mean CCR",
      "and n: CCR >= 80 mL/min (106.3 mL/min, n = 6); 50 <= CCR < 80 mL/min",
      "(65.8 mL/min, n = 6); 30 <= CCR < 50 mL/min (40.2 mL/min, n = 2);",
      "CCR < 30 mL/min (9.1 mL/min, n = 3)."
    ),
    dose_range     = paste(
      "Simulated oral regimens: 150 mg q12h, 250 mg q12h, 300 mg q8h and",
      "600 mg q8h. Tebipenem pivoxil is approved in Japan for pediatric use",
      "only (4 mg/kg twice daily, up to 6 mg/kg); no adult indication is",
      "established, which is what this analysis was designed to support."
    ),
    notes          = paste(
      "Demographics from Ebihara 2025 Methods 4.1 ('The Japanese study",
      "population (n = 17) had the following demographic characteristics:",
      "mean height, 164.3 cm; mean body weight, 61.9 kg; and mean BSA, 1.67",
      "m^2') and Methods 4.4. Per-stratum n sums to 17 (6 + 6 + 2 + 3). The",
      "source reports no age or sex distribution. The word 'Preclinical' in",
      "the source title denotes pre-approval simulation evidence supporting a",
      "Japanese adult regulatory application, NOT an animal study: the model",
      "is a human model throughout. Monte Carlo simulations used 1000 virtual",
      "patients per renal-function stratum in Phoenix 8.3."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed(): Ebihara 2025 estimates nothing. Methods
    # 4.1: 'To our knowledge, there is no population pharmacokinetic model
    # of TBPM for Japanese subjects; therefore, we used mean pharmacokinetic
    # parameter values obtained from clinical trials [15]' (= Nakashima
    # 2009). The paper is a Monte Carlo simulation of these fixed inputs.
    #
    # Stratum suffixes follow parameter-names.md section 'Stratum-suffixed
    # parameters': every stratum carries a suffix, none keeps the bare
    # canonical, and the suffix names the stratum rather than the paper's
    # symbol. Only CL/F, Vd/F and ka are stratum-specific; fu and all three
    # IIV terms are shared across strata and keep their bare canonicals.
    # ------------------------------------------------------------------

    # Apparent clearance CL/F by creatinine-clearance stratum (L/h).
    lcl_crclge80   <- fixed(log(21.738)); label("Apparent clearance CL/F, CRCL >= 80 mL/min (L/h)")      # Ebihara 2025 Table 3, 'CCR >= 80' column, CL/F row = 21.738
    lcl_crcl50to80 <- fixed(log(16.176)); label("Apparent clearance CL/F, 50 <= CRCL < 80 mL/min (L/h)")  # Ebihara 2025 Table 3, '50 <= CCR < 80' column, CL/F row = 16.176
    lcl_crcl30to50 <- fixed(log(8.61));   label("Apparent clearance CL/F, 30 <= CRCL < 50 mL/min (L/h)")  # Ebihara 2025 Table 3, '30 <= CCR < 50' column, CL/F row = 8.61
    lcl_crcllt30   <- fixed(log(2.724));  label("Apparent clearance CL/F, CRCL < 30 mL/min (L/h)")        # Ebihara 2025 Table 3, 'CCR < 30' column, CL/F row = 2.724

    # Apparent volume of distribution Vd/F by creatinine-clearance stratum (L).
    # Note the source's Vd/F is NOT monotone in renal function (34.352 L in
    # the 50-80 stratum exceeds the 26.461 L of the >= 80 stratum, and the
    # 30-50 stratum is the smallest at 17.906 L). Transcribed as printed;
    # the source limitations paragraph attributes the irregularity to the
    # very small renal-impairment cohorts (n = 2 and n = 3).
    lvc_crclge80   <- fixed(log(26.461)); label("Apparent volume of distribution Vd/F, CRCL >= 80 mL/min (L)")      # Ebihara 2025 Table 3, 'CCR >= 80' column, Vd/F row = 26.461
    lvc_crcl50to80 <- fixed(log(34.352)); label("Apparent volume of distribution Vd/F, 50 <= CRCL < 80 mL/min (L)")  # Ebihara 2025 Table 3, '50 <= CCR < 80' column, Vd/F row = 34.352
    lvc_crcl30to50 <- fixed(log(17.906)); label("Apparent volume of distribution Vd/F, 30 <= CRCL < 50 mL/min (L)")  # Ebihara 2025 Table 3, '30 <= CCR < 50' column, Vd/F row = 17.906
    lvc_crcllt30   <- fixed(log(15.754)); label("Apparent volume of distribution Vd/F, CRCL < 30 mL/min (L)")        # Ebihara 2025 Table 3, 'CCR < 30' column, Vd/F row = 15.754

    # First-order absorption rate constant ka by creatinine-clearance
    # stratum (1/h). Also non-monotone as printed (2.5, 1.1, 2.8, 1.7);
    # transcribed as reported. ka does not affect AUC and therefore does not
    # affect the source's target-attainment metric.
    lka_crclge80   <- fixed(log(2.5)); label("Absorption rate constant ka, CRCL >= 80 mL/min (1/h)")      # Ebihara 2025 Table 3, 'CCR >= 80' column, ka row = 2.5
    lka_crcl50to80 <- fixed(log(1.1)); label("Absorption rate constant ka, 50 <= CRCL < 80 mL/min (1/h)")  # Ebihara 2025 Table 3, '50 <= CCR < 80' column, ka row = 1.1
    lka_crcl30to50 <- fixed(log(2.8)); label("Absorption rate constant ka, 30 <= CRCL < 50 mL/min (1/h)")  # Ebihara 2025 Table 3, '30 <= CCR < 50' column, ka row = 2.8
    lka_crcllt30   <- fixed(log(1.7)); label("Absorption rate constant ka, CRCL < 30 mL/min (1/h)")        # Ebihara 2025 Table 3, 'CCR < 30' column, ka row = 1.7

    # Plasma unbound fraction, identical in all four strata.
    fu <- fixed(0.33); label("Fraction unbound in plasma (unitless)")  # Ebihara 2025 Table 3 'Fu' row (0.33 in every column) and Methods 4.1: 'The unbound fraction of TBPM was set at 0.33 [28]' (= Kijima 2009)

    # ------------------------------------------------------------------
    # Inter-individual variability, shared across all four strata and fixed
    # from a prior publication. Ebihara 2025 Discussion: 'The inter-individual
    # variability values used from FRPM (wCL/F = 53%; wVd/F = 37%; wka = 67%)'
    # -- FRPM is faropenem, from Hirooka 2005 [17]. Methods 4.1: the values
    # were 'applied to the apparent total body clearance (CL/F), apparent
    # volume of distribution (Vd/F), and absorption rate constant (ka)'.
    #
    # Scale: the source prints these as bare percentages with no OMEGA block
    # and no exponential-IIV equation, so the omega-SD reading (omega = pct /
    # 100) and the %CV reading (omega = sqrt(log(1 + CV^2))) are both
    # admissible. The choice was resolved against the source's own published
    # PTA grid, which for this linear model reduces to a closed form in the
    # CL/F distribution alone (see the vignette 'Resolving the IIV scale'
    # section): RMSE over the 18 published PTA cells is 1.75 points for
    # omega = 0.53, 2.12 points for the %CV reading, and 4.74 points for a
    # variance reading (refuted). The omega-SD reading is encoded here; the
    # two admissible readings differ by only 6% in omega.
    # ------------------------------------------------------------------
    etalcl ~ fixed(0.2809)  # Ebihara 2025 Discussion wCL/F = 53% read as an omega-SD: 0.53^2 = 0.2809
    etalvc ~ fixed(0.1369)  # Ebihara 2025 Discussion wVd/F = 37% read as an omega-SD: 0.37^2 = 0.1369
    etalka ~ fixed(0.4489)  # Ebihara 2025 Discussion wka   = 67% read as an omega-SD: 0.67^2 = 0.4489

    # Residual error: the source fits no data (Monte Carlo simulation of
    # fixed literature parameters) and reports no residual error model.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # not reported: Ebihara 2025 performs simulation only
  })

  model({
    # 1. Renal-function stratum indicators. Boundaries are the Ebihara 2025
    #    Table 3 column headers / Methods 4.1 group definitions. The four
    #    strata are exhaustive and mutually exclusive, so exactly one
    #    indicator equals 1 for any CRCL >= 0.
    crclGe80   <- (CRCL >= 80)
    crcl50to80 <- (CRCL >= 50) * (CRCL < 80)
    crcl30to50 <- (CRCL >= 30) * (CRCL < 50)
    crclLt30   <- (CRCL < 30)

    # 2. Individual parameters. Each is the indicator-weighted sum of the
    #    four stratum-specific typical values times the shared exponential
    #    eta. The exp() is applied to each theta BEFORE summing rather than
    #    to a sum of thetas, because rxode2 rejects a bare (theta + theta)
    #    subexpression as a mu-reference block.
    cl <- (exp(lcl_crclge80)   * crclGe80   +
           exp(lcl_crcl50to80) * crcl50to80 +
           exp(lcl_crcl30to50) * crcl30to50 +
           exp(lcl_crcllt30)   * crclLt30) * exp(etalcl)
    vc <- (exp(lvc_crclge80)   * crclGe80   +
           exp(lvc_crcl50to80) * crcl50to80 +
           exp(lvc_crcl30to50) * crcl30to50 +
           exp(lvc_crcllt30)   * crclLt30) * exp(etalvc)
    ka <- (exp(lka_crclge80)   * crclGe80   +
           exp(lka_crcl50to80) * crcl50to80 +
           exp(lka_crcl30to50) * crcl30to50 +
           exp(lka_crcllt30)   * crclLt30) * exp(etalka)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system: one compartment with first-order absorption
    #    (Ebihara 2025 Table 3 'Pharmacokinetic Model' row).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Bioavailability is not identifiable from the source: CL/F and Vd/F
    #    are apparent (oral) parameters, so F is folded into them and no
    #    f(depot) term is applied.

    # 6. Observation. Cc is total plasma tebipenem; Cu is the unbound
    #    concentration that drives the source's PK/PD target
    #    fAUC(0-24)/MIC/tau >= 34.58 (Methods 4.3).
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
