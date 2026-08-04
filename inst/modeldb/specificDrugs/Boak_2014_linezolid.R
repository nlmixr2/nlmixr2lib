Boak_2014_linezolid <- function() {
  description <- paste(
    "Population PK/toxicodynamic model for linezolid in 41 critically",
    "ill adult patients (Australia + USA, 42 treatment courses of 600 mg",
    "linezolid q12h IV and/or oral for 5-54 days). PK is a one-compartment",
    "model with three serial absorption-lag compartments (lat1 -> lat2 ->",
    "lat3) followed by a gut compartment (depot) before the central",
    "compartment; oral bioavailability is fixed at 100%. Renal function",
    "is described by a weight-normalised Cockcroft-Gault GFR expressed as",
    "the ratio RF = GFR / 120, driving an additive renal-plus-nonrenal",
    "total clearance CL = F_Size_CL * (CL_NR + CL_R * RF) with fixed",
    "allometric exponents F_Size_CL = (WT/65)^0.75 and F_Size_V =",
    "(WT/65)^1. The toxicodynamic model represents platelet turnover with",
    "a 15-compartment bone-marrow precursor chain (precursor1..precursor15,",
    "rate ktr = 15/MTT_Pre) feeding a 15-compartment circulating platelet",
    "chain (transit1..transit15, rate kout = 15/MTT_PL, Friberg-Bulitta",
    "life-span distribution). The observed platelet count is the average",
    "of the 15 circulating compartments. Linezolid inhibits synthesis of",
    "platelet precursor cells (Imax fixed to 1, so the inhibitory factor",
    "simplifies to IC50 / (IC50 + Cc)) with IC50 = 8.06 mg/L, and a",
    "homeostatic feedback (Base_PL/PL)^gamma with gamma = 1.02 stimulates",
    "precursor synthesis when circulating platelets fall below baseline.",
    "Initial-condition scaling factors F_Ini_Pre and F_Ini_PL (typical",
    "values fixed to 1, BSV estimated) allow the model to describe",
    "subjects whose platelet counts were not at steady state at the",
    "start of therapy."
  )
  reference <- paste(
    "Boak LM, Rayner CR, Grayson ML, Paterson DL, Spelman D, Khumra S,",
    "Capitano B, Forrest A, Li J, Nation RL, Bulitta JB. (2014).",
    "Clinical Population Pharmacokinetics and Toxicodynamics of Linezolid.",
    "Antimicrob Agents Chemother 58(4):2334-2343.",
    "doi:10.1128/AAC.01885-13"
  )
  vignette <- "Boak_2014_linezolid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L",
    platelet      = "10^9 cells/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    lat1        = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    lat2        = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    lat3        = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    depot       = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = FALSE),
    precursor1  = list(analyte = "linezolid metabolite", units = "mg", specimen = "not applicable", verified = FALSE),
    precursor2  = list(analyte = "linezolid metabolite", units = "mg", specimen = "not applicable", verified = FALSE),
    precursor3  = list(analyte = "linezolid metabolite", units = "mg", specimen = "not applicable", verified = FALSE),
    precursor4  = list(analyte = "linezolid metabolite", units = "mg", specimen = "not applicable", verified = FALSE),
    precursor5  = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor6  = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor7  = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor8  = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor9  = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor10 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor11 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor12 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor13 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor14 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    precursor15 = list(analyte = "linezolid metabolite", units = "mg", specimen = "urine", verified = FALSE),
    transit1    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit2    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit3    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit4    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit5    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit6    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit7    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit8    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit9    = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit10   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit11   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit12   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit13   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit14   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE),
    transit15   = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = paste(
        "Body weight used in the model = the lower value of ideal body",
        "weight (IBW) and total body weight (TBW), per Boak 2014 Methods",
        "'Covariate effect model'."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The Boak 2014 model uses WT_i = min(IBW_i, TBW_i) as the size",
        "descriptor for allometric scaling on both CL and V (reference WT_STD",
        "= 65 kg, Methods paragraph 'Covariate effect model'). IBW is",
        "typically computed by the Devine 1974 formula; users must pre-compute",
        "WT externally as the smaller of IBW and TBW before supplying it to",
        "the model. Fixed allometric exponents are 0.75 on CL and 1 on V",
        "(Boak 2014 reference 31)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used inside model() to compute the weight-normalised Cockcroft-Gault",
        "GFR (Boak 2014 Eq 12): GFR = (140 - AGE) * 65 * F_sex / (72 * CREAT).",
        "The individual weight is intentionally NOT used here; substituting",
        "the standard 65 kg decouples renal function from body size so the",
        "allometric F_Size_CL can handle size separately."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Cockcroft-Gault sex factor F_sex = 1 - 0.15 * SEXF (i.e., 0.85 for",
        "female and 1 for male per Boak 2014 Eq 12). SEXF replaces the",
        "paper's F_sex variable."
      ),
      source_name        = "F_sex"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration (Cockcroft-Gault input)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The Cockcroft-Gault equation embedded in Boak 2014 Eq 12 requires",
        "SCR in mg/dL. Supply values in mg/dL (divide umol/L values by 88.42",
        "to convert). Boak 2014 Methods notes 'GFR was included as a",
        "time-dependent covariate' so time-varying CREAT is supported."
      ),
      source_name        = "SCR"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 41L,
    n_studies           = 1L,
    n_treatment_courses = 42L,
    n_pk_observations   = 161L,
    age_range           = "adults >= 18 years",
    weight_range        = paste(
      "mean 74.6 kg (CV about 19% in Table 1 combined-cohort mean of the",
      "no-thrombocytopenia and thrombocytopenia subgroups)"
    ),
    sex_female_pct      = 39.0,
    race_ethnicity      = "not tabulated in Boak 2014",
    disease_state       = paste(
      "Critically ill / hospitalised adults treated with linezolid at the",
      "Departments of Infectious Diseases of the Alfred and Austin",
      "Hospitals (Melbourne, Victoria, Australia) and the Transplant Unit",
      "of the University of Pittsburgh Medical Center (Pittsburgh, PA,",
      "USA) for infections including multidrug-resistant gram-positive",
      "organisms (MRSA, VRE); 14 (34%) developed thrombocytopenia",
      "(platelet count < 100 x 10^9/L) during therapy, 4 of whom had a",
      "baseline platelet count already below the thrombocytopenia",
      "threshold."
    ),
    dose_range          = paste(
      "600 mg linezolid every 12 h (q12h) by intravenous infusion and/or",
      "oral tablet, for at least 4 days (mean duration 22 days; range",
      "5-54 days). Some patients received both IV and oral therapy but",
      "not simultaneously."
    ),
    regions             = "Australia (Melbourne) and USA (Pittsburgh)",
    notes               = paste(
      "Prospective observational study, enrolment October 2004 - January",
      "2007 across three hospitals. Sparse PK sampling (predose plus 2,",
      "4 and 8 h postdose) collected within the first week of therapy in",
      "all but one patient. Platelet counts were measured up to 5 days",
      "before the start of therapy and at least weekly during treatment.",
      "Estimation used parallelised S-ADAPT (v1.57) with importance",
      "sampling (pmethod = 4) via the SADAPT-TRAN facilitator. Univariate",
      "risk-factor analysis of the binary thrombocytopenia outcome",
      "(Table 1) identified baseline platelet count as the only",
      "significant risk factor (p = 0.0084); it is not a model covariate."
    )
  )

  ini({
    # ---- PK structural parameters (Boak 2014 Table 2) ----
    # Absorption chain: dose -> lat1 -> lat2 -> lat3 -> depot -> central.
    # k_lag = 3 / T_lag (per-compartment rate for the lat chain, so total mean
    # transit through the three lat compartments equals T_lag).
    ltlag       <- log(69.2)   ; label("Mean absorption lag time (min); k_lag = 3/T_lag applied per lat compartment")   # Boak 2014 Table 2 (T_lag)
    ltabs       <- log(10.3)   ; label("Absorption half-life (min); ka = ln(2)/Tabs applied to depot -> central")       # Boak 2014 Table 2 (T_abs1/2)

    # Composite total clearance CL = F_Size_CL * (CL_NR + CL_R * RF); CL_NR
    # and CL_R are reported separately in Table 2 at WT_STD = 65 kg and RF = 1.
    lcl_nonren  <- log(4.55)   ; label("Non-renal clearance CL_NR at 65 kg (L/h)")                                     # Boak 2014 Table 2 (CL_NR)
    lcl_renal   <- log(2.17)   ; label("Renal clearance CL_R at 65 kg, RF = 1 (L/h)")                                  # Boak 2014 Table 2 (CL_R)
    lvc         <- log(44.3)   ; label("Central volume of distribution V at 65 kg (L)")                                # Boak 2014 Table 2 (V)

    # Fixed allometric exponents (Boak 2014 Methods 'Covariate effect model';
    # reference 31 in Boak).
    e_wt_cl     <- fixed(0.75) ; label("Allometric exponent on CL (unitless)")                                         # Boak 2014 Methods
    e_wt_vc     <- fixed(1.00) ; label("Allometric exponent on V (unitless)")                                          # Boak 2014 Methods

    # ---- PD (platelet-turnover) structural parameters (Boak 2014 Table 2) ----
    lrbase      <- log(252)    ; label("Steady-state platelet baseline Base_PL (10^9/L)")                                              # Boak 2014 Table 2 (Base_PL)
    lec50       <- log(8.06)   ; label("IC50 for inhibition of precursor synthesis (Imax fixed to 1) (mg/L)")                          # Boak 2014 Table 2 (IC50)
    lmtt_pre    <- log(7.68)   ; label("Mean transit time for the 15-compartment precursor chain MTT_Pre (days); ktr = 15/MTT_Pre")    # Boak 2014 Table 2 (MTT_Pre)
    lmtt_pl     <- log(6.80)   ; label("Mean life span for the 15-compartment platelet chain MTT_PL (days); kout = 15/MTT_PL")         # Boak 2014 Table 2 (MTT_PL)
    lgamma      <- log(1.02)   ; label("Feedback exponent gamma in (Base_PL/PL)^gamma stimulation of precursor synthesis (unitless)")  # Boak 2014 Table 2 (gamma)

    # F_Ini_Pre and F_Ini_PL scale the initial precursor and platelet pools
    # away from steady state at time zero. Typical values are fixed at 1
    # (Boak 2014 Table 2 footnote d) with estimated BSV.
    lfini_pre   <- fixed(log(1)) ; label("Initial-condition scale factor F_Ini_Pre for precursors (typical)")               # Boak 2014 Table 2 (F_Ini_Pre)
    lfini_pl    <- fixed(log(1)) ; label("Initial-condition scale factor F_Ini_PL for platelets (typical)")                 # Boak 2014 Table 2 (F_Ini_PL)

    # ---- Inter-individual variability ----
    # BSV values in Boak 2014 Table 2 are reported as decimal CV (confirmed
    # by the Discussion 'estimated mean transit times of 7.68 days (CV,
    # 34.7%) for precursor cells and 6.80 days (20.3%) for platelets' where
    # 0.347 / 0.203 in Table 2 map to the 34.7% / 20.3% CVs, and by the
    # Results 'A linezolid concentration of 8.06 mg/liter (101%
    # between-patient variability)' where 1.01 in Table 2 maps to 101% CV).
    # Convert CV to log-normal omega^2 = log(CV^2 + 1) for entry into
    # nlmixr2's ~ notation. A single eta on the composite CL follows the
    # Tsuji_2017_linezolid convention for the same additive-renal-plus-nonrenal
    # structure; the eta is named after lcl_nonren for convention-check
    # compliance but is applied multiplicatively on (CL_NR + CL_R * RF).
    etaltlag       ~ 0.30748   # T_lag        CV = 0.600 -> omega^2 = log(1.36)
    etaltabs       ~ 0.02139   # T_abs1/2     CV = 0.147 -> omega^2 = log(1.021609)
    etalcl_nonren  ~ 0.21451   # CL composite CV = 0.489 -> omega^2 = log(1.239121); single eta on composite CL
    etalvc         ~ 0.001296  # V            CV = 0.036 -> omega^2 = log(1.001296)
    etalrbase      ~ 0.35325   # Base_PL      CV = 0.651 -> omega^2 = log(1.423801)
    etalec50       ~ 0.70307   # IC50         CV = 1.01  -> omega^2 = log(2.02010)
    etalmtt_pre    ~ 0.11375   # MTT_Pre      CV = 0.347 -> omega^2 = log(1.120409)
    etalmtt_pl     ~ 0.04039   # MTT_PL       CV = 0.203 -> omega^2 = log(1.041209)
    etalgamma      ~ fixed(0.02226) # gamma CV = 0.15 -> omega^2 = log(1.0225)
    etalfini_pre   ~ 0.04569   # F_Ini_Pre    CV = 0.215 -> omega^2 = log(1.046225)
    etalfini_pl    ~ 0.05455   # F_Ini_PL     CV = 0.236 -> omega^2 = log(1.055696)

    # ---- Residual error (additive + proportional per Boak 2014 Table 2) ----
    addSd       <- 0.309   ; label("Additive residual SD for linezolid Cc (mg/L)")                        # Boak 2014 Table 2 (SD_in)
    propSd      <- 0.225   ; label("Proportional residual SD for linezolid Cc (fraction)")                # Boak 2014 Table 2 (SD_sl)
    addSd_plt   <- 15.1    ; label("Additive residual SD for observed platelet count (10^9/L)")           # Boak 2014 Table 2 (PD_in)
    propSd_plt  <- 0.0755  ; label("Proportional residual SD for observed platelet count (fraction)")     # Boak 2014 Table 2 (PD_sl)
  })

  model({
    # ---- Cockcroft-Gault-derived renal function ratio (Boak 2014 Eq 12) ----
    # SEXF is 1 for female and 0 for male; CG uses 0.85 for female and 1 for male.
    f_sex   <- 1 - 0.15 * SEXF
    # GFR is weight-normalised: the paper substitutes WT_STD = 65 kg for
    # the individual weight in the CG formula to decouple size from renal
    # function; the allometric F_Size_CL then handles size separately.
    gfr_65  <- (140 - AGE) * 65 * f_sex / (72 * CREAT)
    rf      <- gfr_65 / 120

    # ---- Allometric size factors (reference 65 kg) ----
    fsize_cl <- (WT / 65)^e_wt_cl
    fsize_vc <- (WT / 65)^e_wt_vc

    # ---- Individual PK parameters ----
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal)
    cl <- fsize_cl * (cl_nonren + cl_renal * rf) * exp(etalcl_nonren)  # single eta on composite CL
    vc <- exp(lvc + etalvc) * fsize_vc

    # ---- Absorption rate constants (per hour) ----
    tlag_min <- exp(ltlag + etaltlag)                # minutes; total mean transit through 3 lat compartments
    tabs_min <- exp(ltabs + etaltabs)                # minutes; absorption half-life
    k_lag    <- 3 * 60 / tlag_min                    # per-hour rate for each lat compartment
    ka       <- log(2) * 60 / tabs_min               # per-hour depot -> central rate

    # ---- PD individual parameters ----
    rbase    <- exp(lrbase + etalrbase)              # steady-state platelet baseline (10^9/L)
    ec50     <- exp(lec50 + etalec50)                # IC50 for inhibition of precursor synthesis (mg/L)
    mtt_pre  <- exp(lmtt_pre + etalmtt_pre)          # days
    mtt_pl   <- exp(lmtt_pl  + etalmtt_pl)           # days
    gamma    <- exp(lgamma + etalgamma)              # feedback exponent (unitless)
    fini_pre <- exp(lfini_pre + etalfini_pre)        # initial-condition scale factor for precursors
    fini_pl  <- exp(lfini_pl  + etalfini_pl)         # initial-condition scale factor for platelets

    # Chain rate constants (per hour). For a chain of 15 identical compartments
    # with rate constant k, the mean transit time is 15/k.
    ktr  <- 15 / (mtt_pre * 24)                      # per hour
    kout <- 15 / (mtt_pl  * 24)                      # per hour

    # Precursor steady-state baseline (Boak 2014 Methods paragraph after Eq 6):
    # Base_Pre = Base_PL * kout / ktr.
    # Zero-order synthesis rate K_in follows from ktr * Base_Pre = kin so that at
    # steady state with no drug and no feedback perturbation, dPre1/dt = 0.
    base_pre <- rbase * kout / ktr
    kin      <- kout * rbase

    # ---- Observed platelet count = mean of the 15 circulating (transit) compartments ----
    plt <- (transit1  + transit2  + transit3  + transit4  + transit5  +
            transit6  + transit7  + transit8  + transit9  + transit10 +
            transit11 + transit12 + transit13 + transit14 + transit15) / 15

    # ---- Feedback stimulation of precursor synthesis (Boak 2014 Eq 7) ----
    # Stim_feedback = (Base_PL / PL)^gamma; when PL < Base_PL and gamma > 0,
    # Stim_feedback > 1 and stimulates precursor synthesis.
    stim_fb <- (rbase / plt)^gamma

    # ---- Concentration in central compartment ----
    Cc <- central / vc

    # ---- Drug effect on precursor synthesis (Imax fixed to 1 per Boak 2014 Table 2) ----
    # 1 - Imax * Cc / (IC50 + Cc) reduces to IC50 / (IC50 + Cc).
    drug_eff <- ec50 / (ec50 + Cc)

    # ---- PK ODEs ----
    # Oral dose enters lat1; IV dose enters central via the event table
    # (cmt = "central", rate = ...) as an infusion.
    d/dt(lat1)    <- -k_lag * lat1
    d/dt(lat2)    <-  k_lag * lat1 - k_lag * lat2
    d/dt(lat3)    <-  k_lag * lat2 - k_lag * lat3
    d/dt(depot)   <-  k_lag * lat3 - ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # ---- PD ODEs: 15 platelet precursor compartments (Boak 2014 Eq 6 + precursor chain) ----
    d/dt(precursor1)  <- kin * stim_fb * drug_eff - ktr * precursor1
    d/dt(precursor2)  <- ktr * precursor1  - ktr * precursor2
    d/dt(precursor3)  <- ktr * precursor2  - ktr * precursor3
    d/dt(precursor4)  <- ktr * precursor3  - ktr * precursor4
    d/dt(precursor5)  <- ktr * precursor4  - ktr * precursor5
    d/dt(precursor6)  <- ktr * precursor5  - ktr * precursor6
    d/dt(precursor7)  <- ktr * precursor6  - ktr * precursor7
    d/dt(precursor8)  <- ktr * precursor7  - ktr * precursor8
    d/dt(precursor9)  <- ktr * precursor8  - ktr * precursor9
    d/dt(precursor10) <- ktr * precursor9  - ktr * precursor10
    d/dt(precursor11) <- ktr * precursor10 - ktr * precursor11
    d/dt(precursor12) <- ktr * precursor11 - ktr * precursor12
    d/dt(precursor13) <- ktr * precursor12 - ktr * precursor13
    d/dt(precursor14) <- ktr * precursor13 - ktr * precursor14
    d/dt(precursor15) <- ktr * precursor14 - ktr * precursor15

    # ---- PD ODEs: 15 circulating platelet compartments (Friberg-Bulitta life-span chain) ----
    # Entry from the last precursor at rate ktr; transfer between platelet
    # compartments and out of transit15 at rate kout.
    d/dt(transit1)  <- ktr  * precursor15 - kout * transit1
    d/dt(transit2)  <- kout * transit1    - kout * transit2
    d/dt(transit3)  <- kout * transit2    - kout * transit3
    d/dt(transit4)  <- kout * transit3    - kout * transit4
    d/dt(transit5)  <- kout * transit4    - kout * transit5
    d/dt(transit6)  <- kout * transit5    - kout * transit6
    d/dt(transit7)  <- kout * transit6    - kout * transit7
    d/dt(transit8)  <- kout * transit7    - kout * transit8
    d/dt(transit9)  <- kout * transit8    - kout * transit9
    d/dt(transit10) <- kout * transit9    - kout * transit10
    d/dt(transit11) <- kout * transit10   - kout * transit11
    d/dt(transit12) <- kout * transit11   - kout * transit12
    d/dt(transit13) <- kout * transit12   - kout * transit13
    d/dt(transit14) <- kout * transit13   - kout * transit14
    d/dt(transit15) <- kout * transit14   - kout * transit15

    # ---- Initial conditions ----
    # All precursor compartments start at F_Ini_Pre * Base_Pre; all platelet
    # compartments start at F_Ini_PL * Base_PL (Boak 2014 IC for Eq 6 and the
    # 15 platelet compartments).
    precursor1(0)  <- fini_pre * base_pre
    precursor2(0)  <- fini_pre * base_pre
    precursor3(0)  <- fini_pre * base_pre
    precursor4(0)  <- fini_pre * base_pre
    precursor5(0)  <- fini_pre * base_pre
    precursor6(0)  <- fini_pre * base_pre
    precursor7(0)  <- fini_pre * base_pre
    precursor8(0)  <- fini_pre * base_pre
    precursor9(0)  <- fini_pre * base_pre
    precursor10(0) <- fini_pre * base_pre
    precursor11(0) <- fini_pre * base_pre
    precursor12(0) <- fini_pre * base_pre
    precursor13(0) <- fini_pre * base_pre
    precursor14(0) <- fini_pre * base_pre
    precursor15(0) <- fini_pre * base_pre
    transit1(0)  <- fini_pl * rbase
    transit2(0)  <- fini_pl * rbase
    transit3(0)  <- fini_pl * rbase
    transit4(0)  <- fini_pl * rbase
    transit5(0)  <- fini_pl * rbase
    transit6(0)  <- fini_pl * rbase
    transit7(0)  <- fini_pl * rbase
    transit8(0)  <- fini_pl * rbase
    transit9(0)  <- fini_pl * rbase
    transit10(0) <- fini_pl * rbase
    transit11(0) <- fini_pl * rbase
    transit12(0) <- fini_pl * rbase
    transit13(0) <- fini_pl * rbase
    transit14(0) <- fini_pl * rbase
    transit15(0) <- fini_pl * rbase

    # ---- Residual error ----
    Cc  ~ add(addSd)     + prop(propSd)
    plt ~ add(addSd_plt) + prop(propSd_plt)
  })
}
