Riggs_2014_empagliflozin <- function() {
  description <- paste(
    "Exposure-response (E-R) model for empagliflozin in patients with type 2",
    "diabetes mellitus (T2DM), linking steady-state exposure to 24 h urinary",
    "glucose excretion (UGE), fasting plasma glucose (FPG) and HbA1c (Riggs",
    "2014 Br J Clin Pharmacol). The PK layer is the two-compartment,",
    "lagged-first-order-absorption population PK model whose parameter values",
    "Riggs 2014 restates in Results, with an imposed allometric body-weight",
    "effect. The dosing-interval AUC (AUCss = DOSE_EMPA_MGD * 1e6 / MW / CL)",
    "drives a single stimulation function STIM that simultaneously lowers FPG",
    "and raises UGE. UGE is an algebraic sum of an exponential baseline in FPG",
    "and a saturable drug term (Umax, Ustim50); FPG is a turnover state whose",
    "removal is enhanced by STIM; FPG in turn drives first-order HbA1c",
    "production against a physiologic-limit boundary. Study-specific baselines",
    "and study-specific UGE / potency factors are carried for all five",
    "contributing studies (A-E). Data: 974 T2DM patients, 1-100 mg once daily,",
    "up to 12 weeks, from five randomised placebo-controlled trials."
  )
  reference <- paste(
    "Riggs MM, Seman LJ, Staab A, MacGregor TR, Gillespie W, Gastonguay MR,",
    "Woerle HJ, Macha S. Exposure-response modelling for empagliflozin, a",
    "sodium glucose cotransporter 2 (SGLT2) inhibitor, in patients with type 2",
    "diabetes. Br J Clin Pharmacol. 2014;78(6):1407-1418.",
    "doi:10.1111/bcp.12453. The allometric exponents and the 70 kg reference",
    "weight of the PK layer are not printed in Riggs 2014 (which states only",
    "that 'weight was included allometrically on each of these parameters');",
    "they are taken from the companion population PK paper whose model Riggs",
    "2014 uses: Riggs MM, Staab A, Seman L, MacGregor TR, Bergsma TT,",
    "Gastonguay MR, Macha S. Population pharmacokinetics of empagliflozin, a",
    "sodium glucose cotransporter 2 inhibitor, in patients with type 2",
    "diabetes. J Clin Pharmacol. 2013;53(10):1028-1038. doi:10.1002/jcph.147",
    "(Table 3: CL/F and Q/F scale as (WT/70)^0.75 FIXED; V2/F and V3/F as",
    "(WT/70)^1 FIXED)."
  )
  vignette <- "Riggs_2014_empagliflozin"
  units <- list(
    time          = "h",
    dosing        = "mg empagliflozin (oral, once daily)",
    concentration = paste(
      "Cc in nmol/L (= nM; converted from mg/L via MW 450.9 g/mol as stated in",
      "Riggs 2014 Results); AUC in nmol*h/L; FPG in mmol/L; HbA1c in % (NGSP);",
      "UGE in g glucose per 24 h"
    )
  )
  paper_specific_compartments <- c("hba1c")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    glucose     = list(analyte = "glucose", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hba1c       = list(analyte = "HbA1c", units = "% (NGSP)", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Imposed (not estimated) allometric effect on all four disposition",
        "parameters, normalised to a 70 kg reference: (WT/70)^0.75 on CL/F and",
        "Q/F, (WT/70)^1 on V2/F (vc) and V3/F (vp). Riggs 2014 Results states",
        "only that 'weight was included allometrically on each of these",
        "parameters'; the exponents and the 70 kg reference are read from the",
        "companion population PK paper Riggs 2013 J Clin Pharmacol Table 3,",
        "where all four rows are flagged FIXED. Cohort mean weight 85 kg",
        "(range 44 - 152) across the five studies."
      ),
      source_name        = "WT"
    ),
    DOSE_EMPA_MGD = list(
      description        = "Patient's own once-daily empagliflozin dose at the current dosing record",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the dosing-interval exposure that feeds the STIM function:",
        "AUC [nmol*h/L] = DOSE_EMPA_MGD * 1e6 / MW_empa / cl with MW_empa =",
        "450.9 g/mol (Riggs 2014 Results). For a linear disposition model and",
        "once-daily dosing this identity is exact for the steady-state dosing",
        "interval and reproduces the four typical AUC values Riggs 2014 prints",
        "in Results to within rounding: 1 mg -> 225, 3 mg -> 674, 10 mg -> 2250",
        "and 25 mg -> 5620 nmol*h/L. Set to 0 for placebo arms (AUC = 0, STIM =",
        "0). Studied doses 1, 2.5, 5, 10, 25, 50 and 100 mg once daily. Supply",
        "per dose record, constant within an inter-dose interval, and keep",
        "consistent with the `amt` administered into `depot` (also in mg)."
      ),
      source_name        = "DOSE"
    ),
    STUDY_EMPA_A = list(
      description        = "Study A cohort indicator (EudraCT 2007-000654-32; phase I, 8 days, Germany)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is not in study A)",
      notes              = paste(
        "Selects the study A typical baseline FPG (7.85 mmol/L, Table 2",
        "theta_1) and is the reference level for the UGE baseline, gamma_base,",
        "Umax, Ustim50 and C*50 parameters. Study A contributed no HbA1c to the",
        "E-R analysis (Table 1 footnote), so hba1c_base is not defined for this",
        "study. Exactly one of STUDY_EMPA_A..E must be 1."
      ),
      source_name        = "Study A"
    ),
    STUDY_EMPA_B = list(
      description        = "Study B cohort indicator (NCT00558571; phase I, 4 weeks, Germany)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is not in study B)",
      notes              = paste(
        "Selects the study B typical baseline FPG (8.50 mmol/L, Table 2",
        "theta_2) and baseline HbA1c (7.18 pct, theta_23), and scales the UGE",
        "baseline by theta_20 = 0.320. Shares the reference gamma_base, Umax,",
        "Ustim50 and C*50 with study A. Exactly one of STUDY_EMPA_A..E must be 1."
      ),
      source_name        = "Study B"
    ),
    STUDY_EMPA_C = list(
      description        = "Study C cohort indicator (NCT00885118; phase II, 4 weeks, Japan, Japanese patients only)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is not in study C)",
      notes              = paste(
        "The only study with its own value for every UGE / potency parameter:",
        "UGE baseline x theta_16 = 0.632, gamma_base x theta_17 = 1.16, Umax x",
        "theta_14 = 1.11, Ustim50 x theta_15 = 1.58 and C*50 x theta_19 = 0.169",
        "(Table 2). Also selects baseline FPG 8.76 mmol/L (theta_18) and",
        "baseline HbA1c 7.85 pct (theta_28). Riggs 2014 Discussion attributes",
        "the greater response at low doses in study C partly to the lower body",
        "weight of the Japanese cohort (mean 67.9 kg). Exactly one of",
        "STUDY_EMPA_A..E must be 1."
      ),
      source_name        = "Study C"
    ),
    STUDY_EMPA_D = list(
      description        = "Study D cohort indicator (NCT00789035; phase IIb, 12 weeks, multinational, empagliflozin monotherapy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is not in study D)",
      notes              = paste(
        "Selects baseline FPG 9.30 mmol/L (Table 2 theta_3) and baseline HbA1c",
        "7.85 pct (theta_24). Shares the reference C*50 with studies A and B",
        "(the pooled A + B + D estimate is the primary AUC50 = 626 nmol*h/L).",
        "Study D contributed no UGE observations, so the UGE parameters are not",
        "identified for this study and fall back to the study A reference.",
        "Exactly one of STUDY_EMPA_A..E must be 1."
      ),
      source_name        = "Study D"
    ),
    STUDY_EMPA_E = list(
      description        = "Study E cohort indicator (NCT00749190; phase IIb, 12 weeks, multinational, on background metformin)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject is not in study E)",
      notes              = paste(
        "Selects baseline FPG 9.49 mmol/L (Table 2 theta_4), baseline HbA1c",
        "7.89 pct (theta_25) and scales C*50 by theta_21 = 1.93, giving the",
        "study E AUC50 of 1210 nmol*h/L (Table 2 'Calculated parameters').",
        "Study E contributed no UGE observations, so the UGE parameters fall",
        "back to the study A reference. Exactly one of STUDY_EMPA_A..E must be 1."
      ),
      source_name        = "Study E"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description        = "Cockcroft-Gault estimated creatinine clearance (screened graphically, not retained)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Riggs 2014 Results investigated renal function on the efficacy E-R",
        "graphically (Supplementary Figures S1 and S2) and reported 'no",
        "apparent influence of creatinine clearance ... on either FPG or HbA1c",
        "response down to the minimum observed CLcr of approximately",
        "50 mL/min'. No coefficient is reported, so no term is carried.",
        "Cohort means 94 - 117 mL/min by study (Table 1)."
      ),
      source_name        = "CLcr"
    ),
    SEXF = list(
      description        = "Sex (female indicator; screened on the tolerability E-R only, no effect retained)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Riggs 2014 Methods states 'Gender was also considered as a covariate'",
        "for the exposure-tolerability relationships, which were described by",
        "non-parametric GAM smooths rather than by a parametric model. Figure 6",
        "shows the smooths separately for females and males and Results reports",
        "no increase in adverse-event probability with exposure up to 50 mg.",
        "No parametric coefficient exists to carry. Cohort 574 male / 400",
        "female (Table 1)."
      ),
      source_name        = "Gender"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 974L,
    n_subjects_efficacy  = 974L,
    n_subjects_tolerability = 748L,
    n_studies            = 5L,
    studies              = paste(
      "Study A (EudraCT 2007-000654-32; phase I, n = 48, 8 days, Germany);",
      "study B (NCT00558571; phase I, n = 78, 4 weeks, Germany); study C",
      "(NCT00885118; phase II, n = 100, 4 weeks, Japan); study D",
      "(NCT00789035; phase IIb, n = 324 analysed, 12 weeks, multinational);",
      "study E (NCT00749190; phase IIb, n = 424, 12 weeks, multinational, on",
      "background metformin). Open-label metformin (study D) and open-label",
      "sitagliptin (study E) arms were excluded from the E-R analyses."
    ),
    age_range            = "28 - 80 years (study means 56.7 - 58.4)",
    weight_range         = "44 - 152 kg (study means 67.9 - 94.6; study C, the Japanese cohort, is the lightest)",
    sex_female_pct       = 41.1,
    race_ethnicity       = c(White = 77.1, Asian = 21.8, Black = 0.8, HawaiianPacific = 0.2),
    disease_state        = paste(
      "Type 2 diabetes mellitus. Baseline FPG study means 8.3 - 9.7 mmol/L",
      "(overall range 2.8 - 21.0); baseline HbA1c study means 7.1 - 8.1 pct",
      "(overall range 5.6 - 10.4). Baseline serum creatinine 0.8 - 0.9 mg/dL",
      "and Cockcroft-Gault creatinine clearance 94 - 117 mL/min by study;",
      "fewer than 1.5 pct of patients had CLcr < 50 mL/min."
    ),
    dose_range           = "1 - 100 mg empagliflozin orally once daily for up to 12 weeks (1, 2.5, 5, 10, 25, 50 and 100 mg arms across the five studies), plus placebo",
    regions              = "Germany (studies A and B), Japan (study C), multinational (studies D and E)",
    endpoints            = paste(
      "24 h urinary glucose excretion (studies A, B, C only), fasting plasma",
      "glucose (all five studies) and HbA1c (studies B, C, D, E; study A HbA1c",
      "was not included in the E-R analyses per the Table 1 footnote).",
      "Tolerability endpoints (hypoglycaemia n = 4, urinary tract infection",
      "n = 17, genital / vulvovaginal events n = 16) were evaluated in studies",
      "D and E by non-parametric GAM smoothing and produced no parametric",
      "exposure-response model."
    ),
    notes                = paste(
      "Baseline demographics are tabulated per study in Riggs 2014 Table 1.",
      "The percentages above are pooled across the five studies from the Table",
      "1 counts. Exposures were individual dosing-interval AUCs generated by",
      "the companion population PK model (Riggs 2013 J Clin Pharmacol), whose",
      "structural parameter values Riggs 2014 restates in Results."
    )
  )

  ini({
    # ================================================================
    # PK layer. Riggs 2014 Results, 'Efficacy markers ... relationships
    # with empagliflozin exposure': "Population estimates (interindividual
    # variance [IIV] estimates, coefficient of variation [CV%]) of CL/F,
    # central and peripheral volumes of distribution and inter-compartmental
    # clearance were 9.87 l h-1 (26.9%), 3.02 l, 60.4 l (30.8%) and
    # 5.16 l h-1"; "Following a 0.5 h lag, the typical oral absorption rate
    # constant was estimated to be 0.224 h-1 (IIV, 15.2%)".
    # Independent check: these four disposition values give a terminal
    # half-life of 12.5 h, matching the "approximately 12 h" the same
    # paragraph reports.
    # ================================================================
    lcl <- log(9.87)  ; label("Apparent oral clearance CL/F (L/h) at the 70 kg reference weight")            # Riggs 2014 Results: CL/F = 9.87 L/h
    lvc <- log(3.02)  ; label("Apparent central volume V2/F (L) at the 70 kg reference weight; V2 -> vc per canonical compartment names")   # Riggs 2014 Results: central volume = 3.02 L
    lvp <- log(60.4)  ; label("Apparent peripheral volume V3/F (L) at the 70 kg reference weight; V3 -> vp per canonical compartment names") # Riggs 2014 Results: peripheral volume = 60.4 L
    lq  <- log(5.16)  ; label("Apparent inter-compartmental clearance Q/F (L/h) at the 70 kg reference weight")  # Riggs 2014 Results: inter-compartmental clearance = 5.16 L/h
    lka <- log(0.224) ; label("First-order oral absorption rate constant ka (1/h)")                            # Riggs 2014 Results: ka = 0.224 1/h
    lalag <- fixed(log(0.5)) ; label("Absorption lag time ALAG1 (h)")                                          # Riggs 2014 Results: "Following a 0.5 h lag"; FIXED in the companion popPK model (Riggs 2013 Table 3)

    # Imposed allometric weight exponents. Riggs 2014 says only that "weight
    # was included allometrically on each of these parameters"; the exponents
    # and the 70 kg reference come from Riggs 2013 J Clin Pharmacol Table 3,
    # where every row is flagged FIXED (non-paper provenance: companion
    # population PK publication, on disk).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL/F")   # Riggs 2013 Table 3 row 'CL/F (WT/70)^theta_20' = 0.75 FIXED
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent of (WT/70) on Q/F")    # Riggs 2013 Table 3 row 'Q/F (WT/70)^theta_22' = 0.75 FIXED
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT/70) on V2/F (vc)")  # Riggs 2013 Table 3 row 'V2/F (WT/70)^theta_21' = 1 FIXED
    e_wt_vp <- fixed(1)    ; label("Allometric exponent of (WT/70) on V3/F (vp)")  # Riggs 2013 Table 3 row 'V3/F (WT/70)^theta_23' = 1 FIXED

    # PK IIV, back-transformed from the CV% Riggs 2014 reports in Results
    # via omega^2 = log(1 + CV^2). No IIV was reported on V2/F or Q/F, and
    # Riggs 2014 reports no correlations among the PK etas.
    etalcl ~ 0.0698695   # CV 26.9% -> log(1 + 0.269^2)
    etalvp ~ 0.0906272   # CV 30.8% -> log(1 + 0.308^2)
    etalka ~ 0.0228412   # CV 15.2% -> log(1 + 0.152^2)

    # ================================================================
    # UGE model (Equation 1). Study A is the reference; studies B and C
    # scale the baseline, and study C additionally scales gamma_base,
    # Umax and Ustim50. Studies D and E contributed no UGE data.
    # ================================================================
    luge_base    <- log(3.71)  ; label("Baseline 24 h UGE at FPG = 8 mmol/L, study A (g per 24 h)")   # Table 2 theta_12 = 3.71 (1.79 pct RSE)
    lf_uge_base_b <- log(0.320); label("Multiplicative factor on baseline UGE for study B (log scale)")  # Table 2 theta_20 = 0.320 (18.6 pct RSE)
    lf_uge_base_c <- log(0.632); label("Multiplicative factor on baseline UGE for study C (log scale)")  # Table 2 theta_16 = 0.632 (51.3 pct RSE)

    lgamma_base    <- log(5.31); label("Exponent gamma_base of (FPG/8) in the UGE baseline term")        # Table 2 theta_13 = 5.31 (5.22 pct RSE)
    lf_gamma_base_c <- log(1.16); label("Multiplicative factor on gamma_base for study C (log scale)")   # Table 2 theta_17 = 1.16 (118 pct RSE)

    lumax    <- log(121)  ; label("Maximal drug-related 24 h UGE increase Umax (g per 24 h)")            # Table 2 theta_10 = 121 (1.04 pct RSE)
    lf_umax_c <- log(1.11); label("Multiplicative factor on Umax for study C (log scale)")               # Table 2 theta_14 = 1.11 (177 pct RSE)

    lustim50    <- log(0.590); label("Ustim50: value of FPG * STIM giving half of Umax (mmol/L)")        # Table 2 theta_11 = 0.590 (69.8 pct RSE)
    lf_ustim50_c <- log(1.58); label("Multiplicative factor on Ustim50 for study C (log scale)")         # Table 2 theta_15 = 1.58 (105 pct RSE)

    # IIV on baseline UGE, from the reported CV% (Riggs 2014 Results:
    # "IIV for baseline UGE (CV%) was estimated to be 158.4%").
    etaluge_base ~ 1.2553471   # CV 158.4% -> log(1 + 1.584^2)

    # ================================================================
    # STIM function (Equation 2), in the paper's estimation-stable
    # reparameterisation (Riggs 2014 Results cites refs 16-17 for it):
    #   Emax  = (beta + 1) / beta * Emax_truncated
    #   AUC50 = C*50 / beta
    # Both identities are checkable against the Table 2 'Calculated
    # parameters' block: 1.795/0.795 * 0.0701 = 0.1583 vs a reported
    # 0.158, and 498/0.795 = 626.4 vs a reported 626.
    # ================================================================
    emax_trunc <- 0.0701 ; label("Truncated maximal fractional FPG reduction Emax,truncated")     # Table 2 theta_6 = 0.0701 (18.2 pct RSE)
    beta_fpg   <- 0.795  ; label("Reparameterisation constant beta_FPG,stimulation")              # Table 2 theta_8 = 0.795 (30.4 pct RSE)
    gamma_fpg  <- 1.47   ; label("Exponent gamma_FPG,stimulation of (baseline FPG / 8) on the maximal stimulation")  # Table 2 theta_9 = 1.47 (30.8 pct RSE)

    lc50star     <- fixed(log(498)); label("C*50: reparameterised half-maximal exposure constant (nmol*h/L), studies A + B + D") # Table 2 theta_7 = 498, RSE reported as FIXED
    lf_c50star_c <- log(0.169)     ; label("Multiplicative factor on C*50 for study C (log scale)")  # Table 2 theta_19 = 0.169 (22.8 pct RSE)
    lf_c50star_e <- log(1.93)      ; label("Multiplicative factor on C*50 for study E (log scale)")  # Table 2 theta_21 = 1.93 (RSE not reported)

    # ================================================================
    # FPG turnover (Equation 3) and its study-specific typical baselines.
    # ================================================================
    lkfpgout <- log(0.0407) ; label("First-order FPG removal rate constant kFPG,out (1/h)")   # Table 2 theta_5 = 0.0407 1/h (18 pct RSE)

    fpg_base_a <- 7.85 ; label("Typical baseline FPG, study A (mmol/L)")   # Table 2 theta_1  = 7.85 (1.29 pct RSE)
    fpg_base_b <- 8.50 ; label("Typical baseline FPG, study B (mmol/L)")   # Table 2 theta_2  = 8.50 (1.38 pct RSE)
    fpg_base_c <- 8.76 ; label("Typical baseline FPG, study C (mmol/L)")   # Table 2 theta_18 = 8.76 (0.81 pct RSE)
    fpg_base_d <- 9.30 ; label("Typical baseline FPG, study D (mmol/L)")   # Table 2 theta_3  = 9.30 (0.493 pct RSE)
    fpg_base_e <- 9.49 ; label("Typical baseline FPG, study E (mmol/L)")   # Table 2 theta_4  = 9.49 (0.44 pct RSE)

    # ================================================================
    # HbA1c turnover (Equation 4). Riggs 2014 reports kHbA1c,out in
    # week^-1 and kHbA1c,in in % week^-1 mmol/L^-1; both are converted to
    # the model's hour time base by dividing by 168 h/week.
    # ================================================================
    lkhba1cout <- log(0.167 / 168) ; label("First-order HbA1c removal rate constant kHbA1c,out (1/h)")   # Table 2 theta_26 = 0.167 week^-1 -> /168 h per week
    lkhba1cin  <- log(0.078 / 168) ; label("HbA1c production rate constant kHbA1c,in (pct per h per mmol/L FPG)")  # Table 2 theta_27 = 0.078 pct week^-1 mmol/L^-1 -> /168
    hba1c_limit <- 3.34 ; label("HbA1c physiologic limit parameter HbA1c,limit (pct)")   # Table 2 theta_22 = 3.34 pct (bootstrap median 3.52)
    eta_share   <- 2.7  ; label("Shared-eta scaling: eta on kHbA1c,out equals theta_29 times the eta on kHbA1c,in")  # Table 2 theta_29 = 2.7 (bootstrap median 2.59)

    # Baseline HbA1c is study-specific. Study B carries the un-suffixed
    # name so that the IIV below pairs with a fixed effect of the same
    # stem; each study's value is an absolute Table 2 theta in its own
    # right, written on the log scale so the log-normal IIV is additive.
    lhba1c_base   <- log(7.18) ; label("Typical baseline HbA1c, study B (pct)")   # Table 2 theta_23 = 7.18 (bootstrap median 7.16)
    lhba1c_base_c <- log(7.85) ; label("Typical baseline HbA1c, study C (pct)")   # Table 2 theta_28 = 7.85 (bootstrap median 7.85)
    lhba1c_base_d <- log(7.85) ; label("Typical baseline HbA1c, study D (pct)")   # Table 2 theta_24 = 7.85 (bootstrap median 7.85)
    lhba1c_base_e <- log(7.89) ; label("Typical baseline HbA1c, study E (pct)")   # Table 2 theta_25 = 7.89 (bootstrap median 7.89)

    # Riggs 2014 Results: "IIV estimates (CV%) for baseline HbA1c and
    # kHbA1c,in were 9.53% and 8.23%, respectively, with a correlation
    # estimate of -0.310." Back-transformed as omega^2 = log(1 + CV^2)
    # with covariance = rho * omega_1 * omega_2.
    etalhba1c_base + etalkhba1cin ~ c(0.00904107,
                                      -0.00242173, 0.00675043)

    # ================================================================
    # Residual error. Table 2 'Residual variance' reports a variance in
    # the point-estimate column and a CV% in the right-hand column; the
    # two are related by CV = sqrt(exp(variance) - 1) for all three
    # endpoints (0.380 -> 68.0 vs 67.9 pct; 0.01461 -> 12.13 vs 12.1 pct;
    # 0.001287 -> 3.589 vs 3.59 pct), i.e. the residuals are log-normal
    # (exponential), which maps to `~ lnorm(...)` rather than `~ prop(...)`.
    # ================================================================
    expSd_uge     <- sqrt(0.380)    ; label("Log-normal residual SD for 24 h UGE (67.9 pct CV)")   # Table 2 residual variance UGE = 0.380 (11.9 pct RSE)
    expSd_glucose <- sqrt(0.01461)  ; label("Log-normal residual SD for FPG (12.1 pct CV)")        # Table 2 residual variance FPG = 0.01461
    expSd_hba1c   <- sqrt(0.001287) ; label("Log-normal residual SD for HbA1c (3.59 pct CV)")      # Table 2 residual variance HbA1c = 0.001287
  })

  model({
    # ---- Constants ----
    # Empagliflozin molecular weight, stated in Riggs 2014 Results
    # ("molecular weight = 450.9 g mol-1").
    mw_empa     <- 450.9
    nmol_per_mg <- 1e6 / mw_empa

    ref_wt  <- 70   # kg (allometric reference weight; Riggs 2013 Table 3)
    ref_fpg <- 8    # mmol/L (144 mg/dL); the FPG normalisation constant in Equations 1 and 2

    # ---- Individual PK parameters (allometric weight effect only) ----
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc)          * (WT / ref_wt)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / ref_wt)^e_wt_vp
    q  <- exp(lq)           * (WT / ref_wt)^e_wt_q
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PK ODEs (two-compartment oral with an absorption lag) ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    alag(depot)       <- exp(lalag)

    # Plasma empagliflozin in nmol/L; the central state holds mg.
    Cc <- central / vc * nmol_per_mg

    # ---- Exposure driver ----
    # AUC_ijk in Equations 2 is the individual dosing-interval AUC. For a
    # linear disposition model dosed once daily, the steady-state
    # dosing-interval AUC is exactly Dose / (CL/F); this is the same
    # identity the successor analysis (Baron 2016) uses, and it reproduces
    # the four typical AUC values printed in Riggs 2014 Results.
    auc <- DOSE_EMPA_MGD * nmol_per_mg / cl

    # ---- Study selection (exactly one indicator is 1) ----
    fpg_base <- STUDY_EMPA_A * fpg_base_a +
                STUDY_EMPA_B * fpg_base_b +
                STUDY_EMPA_C * fpg_base_c +
                STUDY_EMPA_D * fpg_base_d +
                STUDY_EMPA_E * fpg_base_e

    # Study A contributed no HbA1c to the E-R analysis, so it has no
    # theta of its own; a study A subject is given the study B baseline
    # (the other short phase I cohort) purely so that the HbA1c state has
    # a defined initial condition. HbA1c CHANGES are unaffected by this
    # choice because the drug effect enters only through FPG.
    lhba1c_base_i <- (STUDY_EMPA_A + STUDY_EMPA_B) * lhba1c_base +
                     STUDY_EMPA_C * lhba1c_base_c +
                     STUDY_EMPA_D * lhba1c_base_d +
                     STUDY_EMPA_E * lhba1c_base_e
    hba1c_base <- exp(lhba1c_base_i + etalhba1c_base)

    # ---- UGE parameters (Equation 1); study C scales all four ----
    uge_base   <- exp(luge_base +
                      lf_uge_base_b * STUDY_EMPA_B +
                      lf_uge_base_c * STUDY_EMPA_C +
                      etaluge_base)
    gamma_base <- exp(lgamma_base + lf_gamma_base_c * STUDY_EMPA_C)
    umax       <- exp(lumax       + lf_umax_c       * STUDY_EMPA_C)
    ustim50    <- exp(lustim50    + lf_ustim50_c    * STUDY_EMPA_C)

    # ---- STIM function (Equation 2) ----
    c50star <- exp(lc50star +
                   lf_c50star_c * STUDY_EMPA_C +
                   lf_c50star_e * STUDY_EMPA_E)

    # Back out the interpretable Emax / AUC50 pair from the reported
    # reparameterisation. For the reference studies these evaluate to
    # 0.1583 and 626.4, matching the Table 2 'Calculated parameters'
    # rows (0.158 and 626); for study C 106 and for study E 1210.
    emax_fpg <- (beta_fpg + 1) / beta_fpg * emax_trunc
    auc50    <- c50star / beta_fpg

    stim <- (fpg_base / ref_fpg)^gamma_fpg * emax_fpg * auc / (auc50 + auc)

    # ---- FPG turnover (Equation 3) ----
    # STIM sits on the kFPG,out arrow in the Figure 3 schematic and the
    # Results text says it "directly affected the removal rate of FPG",
    # so the drug enhances removal. The steady state is
    # FPG_ss = fpg_base * (1 - STIM), i.e. the maximal fractional FPG
    # decrease equals STIM itself -- which is exactly how Table 2 labels
    # the calculated 0.158 ("FPG maximal decrease (proportional)") and
    # what reproduces all four maximal decreases Riggs 2014 reports
    # (1.3, 1.0, 1.7 and 2.2 mmol/L at baseline FPG 8, 7.4, 9.1 and
    # 10 mmol/L). See the vignette Errata: the printed Equation 3 has
    # the (1 - STIM) factor MULTIPLYING the removal term, which would
    # make empagliflozin raise FPG and contradicts the schematic, the
    # narrative and every number the paper reports.
    kfpgout <- exp(lkfpgout)
    kfpgin  <- kfpgout * fpg_base
    d/dt(glucose) <- kfpgin - kfpgout * glucose / (1 - stim)
    glucose(0)    <- fpg_base

    # ---- UGE (Equation 1), evaluated at the current FPG ----
    # An exponential baseline term in FPG plus a saturable drug term in
    # (FPG * STIM). The two are kept as named intermediates because Riggs
    # 2014 quotes them separately: the drug term alone is the "UGE increase"
    # of the Results narrative (72 and 75 g per 24 h at 10 and 25 mg for a
    # baseline FPG of 8 mmol/L), while the baseline term alone is the
    # placebo relationship of Figure 4.
    uge_baseline <- uge_base * (glucose / ref_fpg)^gamma_base
    uge_drug     <- glucose * umax * stim / (ustim50 + glucose * stim)
    uge          <- uge_baseline + uge_drug

    # ---- HbA1c turnover (Equation 4) ----
    # d(HbA1c)/dt = kHbA1c,in * FPG - kHbA1c,out * HbA1c * (1 - limit/HbA1c)
    # The eta on kHbA1c,out is theta_29 times the eta on kHbA1c,in
    # (Table 2 theta_29 = 2.7), not an independent random effect.
    khba1cin  <- exp(lkhba1cin  + etalkhba1cin)
    khba1cout <- exp(lkhba1cout + eta_share * etalkhba1cin)
    d/dt(hba1c) <- khba1cin * glucose - khba1cout * hba1c * (1 - hba1c_limit / hba1c)
    hba1c(0)    <- hba1c_base

    # ---- Observations ----
    # Riggs 2014 reports no residual error for the plasma concentration
    # (the PK model is reproduced from the companion popPK paper), so Cc
    # is returned as a derived quantity without an error model.
    uge     ~ lnorm(expSd_uge)
    glucose ~ lnorm(expSd_glucose)
    hba1c   ~ lnorm(expSd_hba1c)
  })
}
