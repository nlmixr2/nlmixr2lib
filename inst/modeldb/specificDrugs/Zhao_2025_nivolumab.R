Zhao_2025_nivolumab <- function() {
  description <- "Two-compartment population PK model for subcutaneous and intravenous nivolumab (anti-PD-1 IgG4) with first-order subcutaneous absorption, logit-scale bioavailability, and time-varying clearance (sigmoid Emax); the pre-specified model for the CheckMate 67T PK non-inferiority analysis (Zhao 2025)"
  reference <- "Zhao Y, Vezina H, Hu Z, Kondic A, Zhu L, Roy A. Model-informed drug development of subcutaneous nivolumab: comparison of pharmacokinetic analysis methodologies using clinical trial simulation. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2107-2117. doi:10.1002/psp4.70120"
  vignette <- "Zhao_2025_nivolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE because the source control streams
  # (File S1 / File S2, ADVAN4 TRANS4 with S2 = V2) identify the depot as the
  # subcutaneous injection site and compartment 2 as the sampled serum pool.
  compartmentData <- list(
    depot       = list(analyte = "nivolumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 80 kg. The reference value is taken from the source control streams (File S1 and File S2, 'BBWT_R = 80 ; refernce value, kg' [sic]), which are the code that produced the Table S2 estimates. Table S2's explanatory note instead describes the reference subject as 'weighing 75 kg'; the control-stream value governs because it is what the estimation actually used. See the vignette Errata.",
      source_name        = "BBWT"
    ),
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 90 mL/min (File S1 / File S2 'BGFR_R = 90 ; reference value, mL/min'). Source column name is BGFR (baseline eGFR); stored under the canonical CRCL. Table 1 reports eGFR in mL/min/1.73 m^2, so the column is BSA-normalized even though the control-stream comment abbreviates the units to mL/min. The source paper does not name the estimating equation; the sibling nivolumab analysis (Bajaj 2017) used CKD-EPI.",
      source_name        = "BGFR"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "The source encodes sex as SEXN with 1 = male and 2 = female, and applies each effect under 'IF (SEX_I .EQ. 2)', so the paper's reference category is male. That matches the canonical SEXF orientation directly (1 = female, 0 = male) with no sign inversion: the effect is applied as SEXF = 1. Female sex carries an exponential effect on CL and Vc and a multiplicative effect on subcutaneous bioavailability.",
      source_name        = "SEXN"
    ),
    ECOG_GE1 = list(
      description        = "Baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator (1 if ECOG >= 1, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status = 0, i.e., fully active)",
      notes              = "Exponential effect on CL and a multiplicative effect on subcutaneous bioavailability for patients with ECOG >= 1 (File S1 / File S2 'CL_PS_1 = THETA(14); effect of PS 1+:0' applied under 'IF (PS_I .GE. 1)'). Renamed from the source column PS to the canonical ECOG_GE1. The source imputes missing PS to 1 rather than to the reference level 0.",
      source_name        = "PS"
    ),
    TUMTP_GASTRIC = list(
      description        = "Indicator for gastric cancer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types; the source's reference tumor type is second-line NSCLC)",
      notes              = "Exponential effect on CL (File S1 / File S2 'IF (TTYPE_I .EQ. 7) TVCL = TVCL * EXP(CL_GC); reference is TTYPE=1 (NSCLC_2L)'). Decomposed from the source's categorical TTYPEN2 column, in which levels 2, 4, 5, 6, 8 and 9 are collapsed to level 3 ('other') before the effect is applied, so only gastric cancer (level 7) and classical Hodgkin lymphoma (level 10) retain distinct effects. Gastric cancer was 1.6 percent of the pooled dataset (Table 1).",
      source_name        = "TTYPEN2"
    ),
    TUMTP_HODGKIN_CLASSICAL = list(
      description        = "Indicator for classical Hodgkin lymphoma",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types; the source's reference tumor type is second-line NSCLC)",
      notes              = "Exponential effect on CL (File S1 / File S2 'IF (TTYPE_I .EQ. 10) TVCL = TVCL * EXP(CL_CHL); reference is TTYPE=1 (NSCLC_2L)'). Decomposed from the source's categorical TTYPEN2 column, level 10. Classical Hodgkin lymphoma was 7.2 percent of the pooled dataset (Table 1) and carries the largest single covariate effect in the model (CL about 28 percent lower).",
      source_name        = "TTYPEN2"
    )
  )

  # Covariates the source screened but did not retain in the pre-specified
  # model. The control-stream header records the screen explicitly
  # ("COVMODEL: CL ~ BBWT+GFR+SEX+PS+TTYPE+RACE [No LDH/Albumin/Hepa/Age]")
  # and the corresponding $THETA lines are commented out, so no point
  # estimate is available for any of these.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Summarized in Table 1 (overall median 61.0 years, range 18.0-90.0) but explicitly excluded from the covariate model by the control-stream header '[No LDH/Albumin/Hepa/Age]'. No point estimate reported."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Summarized in Table 1 (overall median 3.90 g/dL, 35.4 percent missing) and read into the full-model dataset as column BALB, but excluded from the covariate model by the control-stream header '[No LDH/Albumin/Hepa/Age]'. No point estimate reported."
    ),
    RACE_ASIAN = list(
      description = "Indicator for Asian race",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as CL_RAAS ('effect of Asian:White/Other'). The $THETA line is commented out in both File S1 (';(0.0587) ; CL_RAAA' / ';(-0.0758) ; CL_RAAS') and File S2, and no value appears in Table S2, so the effect was dropped from the pre-specified model. Retained in Bajaj 2017 (see Bajaj_2017_nivolumab.R) but not here."
    ),
    RACE_BLACK = list(
      description = "Indicator for Black / African American race",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as CL_RAAA ('effect of African American: White/Other'); the $THETA line is commented out in both control streams and no value appears in Table S2."
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Excluded from the covariate model by the control-stream header '[No LDH/Albumin/Hepa/Age]'. Not summarized in Table 1 and no point estimate reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3554L,
    n_studies      = 20L,
    age_range      = "median 61.0 years (range 18.0 - 90.0); mean 59.6 (SD 13.0)",
    age_median     = "61.0 years",
    weight_range   = "median 76.7 kg (range 34.1 - 180); mean 78.3 (SD 19.2)",
    weight_median  = "76.7 kg",
    sex_female_pct = 32.8,
    race_ethnicity = c(
      White = 88.6, Asian = 6.2, `Black/African American` = 2.9,
      `Other/unknown` = 2.1, `American Indian/Alaska Native` = 0.1,
      `Native Hawaiian/other Pacific Islander` = 0.0
    ),
    disease_state  = "Advanced / metastatic solid tumors and classical Hodgkin lymphoma (NSCLC 25.9%, melanoma 24.7%, RCC 20.1%, bladder 8.7%, classical Hodgkin lymphoma 7.2%, SCCHN 4.8%, SCLC 2.7%, CRC 1.7%, gastric 1.6%, other 2.4%)",
    dose_range     = "Intravenous 0.1 - 20 mg/kg (and 240 mg flat) 1-hour infusions Q2W or Q3W across 19 studies; subcutaneous 720, 960 and 1200 mg Q4W and 600 mg Q2W co-formulated with rHuPH20 in CheckMate 8KX (CA2098KX)",
    regions        = "Global (phase I / II / III studies)",
    ecog_distribution = "ECOG 0 43.6%, ECOG 1 54.7%, ECOG 2 1.7%, missing 0.1%",
    renal_function = "Baseline eGFR median 84.3 (range 18.1 - 191) mL/min/1.73 m^2; mean 81.9 (SD 23.0)",
    albumin        = "Baseline albumin median 3.90 (range 1.40 - 5.30) g/dL; 35.4% missing",
    n_subjects_sc  = 66L,
    n_subjects_iv  = 3488L,
    notes          = "Baseline demographics per Table 1 (N = 3554: 66 subcutaneous from CheckMate 8KX and 3488 intravenous). Studies pooled per Table S1: MDX1106-01 (CA209001), MDX1106-03 (CA209003), ONO-4538-01 (CA209005), CA209009, CA209010, CA209017, CA209025, CA209026, CA209032, CA209037, CA209039, ONO-4538-02 (CA209-051), CA209057, CA209063, CA209066, CA209067, CA209141, CA209205, CA209275 (19 intravenous studies) plus CA2098KX / CheckMate 8KX (subcutaneous). Table 1 footnote b states 18 intravenous studies; Table S1 lists 19 (see the vignette Errata). Subjects in CheckMate 8KX who received subcutaneous nivolumab without rHuPH20 were excluded from the analysis. The model was pre-specified to analyze CheckMate 67T (NCT04810078), a phase III PK non-inferiority trial of 1200 mg subcutaneous Q4W versus 3 mg/kg intravenous Q2W in previously-treated advanced or metastatic clear-cell RCC."
  )

  ini({
    # Structural parameters (Table S2, 'Parameter estimates of nivolumab
    # prespecified model'). Reference subject: 80 kg, eGFR 90 mL/min, male,
    # ECOG 0, second-line NSCLC.
    #
    # UNITS: Table S2 labels CL0_REF and Q_REF '[mL/hr]' but the values 0.0109
    # and 0.0327 are L/h. Three independent proofs: (1) both control streams
    # annotate the same THETAs 'CL [L/h]' and 'Q [L/h]'; (2) Table S3 reports
    # the same parameters as 10.6-11.0 and 32.4-32.5 mL/hr, i.e. 1000-fold
    # larger under the mL/hr label; (3) 10.9 mL/h is the physiologically
    # expected nivolumab clearance. Converted here to L/day (x 24) because
    # this model keeps time in days.
    lcl <- log(0.0109 * 24); label("Baseline clearance CL0_REF at the reference covariates (L/day)")     # Table S2: CL0_REF = 0.0109 L/h (mislabeled mL/hr); x 24 = 0.2616 L/day
    lvc <- log(4.25);        label("Central volume of distribution Vc_REF (L)")                          # Table S2: Vc_REF = 4.25 L
    lq  <- log(0.0327 * 24); label("Intercompartmental clearance Q_REF (L/day)")                         # Table S2: Q_REF = 0.0327 L/h (mislabeled mL/hr); x 24 = 0.7848 L/day
    lvp <- log(2.63);        label("Peripheral volume of distribution Vp_REF (L)")                       # Table S2: Vp_REF = 2.63 L

    # Subcutaneous absorption. Table S2 reports Ka_REF in day^-1, which this
    # model uses directly; the $THETAP prior in File S2 gives the same value
    # on the model's native hourly scale (0.307 / 24 = 0.0128 ~ 0.013 /h),
    # which confirms the time unit of the printed estimate.
    lka <- log(0.307);                    label("First-order subcutaneous absorption rate constant Ka_REF (1/day)")     # Table S2: Ka_REF = 0.307 /day; File S2 $THETAP KA = 0.013 /h = 0.312 /day
    logitfdepot <- log(0.752 / (1 - 0.752)); label("Logit of subcutaneous bioavailability F_REF (fraction)")             # Table S2: F_REF = 0.752; logit(0.752) = 1.1093

    # Covariate effects on CL (File S1 / File S2 $PK). Power on WT and eGFR;
    # exponential on female sex, ECOG >= 1 and the two retained tumor types.
    e_wt_cl                      <-  0.622; label("Power exponent of WT on CL (unitless)")                              # Table S2: CL_WTB  = 0.622
    e_crcl_cl                    <-  0.139; label("Power exponent of CRCL (eGFR) on CL (unitless)")                     # Table S2: CL_eGFR = 0.139
    e_sex_cl                     <- -0.158; label("Exponential coefficient of female sex on CL (unitless)")             # Table S2: CL_SEX  = -0.158
    e_ecog_ge1_cl                <-  0.174; label("Exponential coefficient of ECOG_GE1 on CL (unitless)")               # Table S2: CL_PS   = 0.174
    e_tumtp_gastric_cl           <-  0.180; label("Exponential coefficient of gastric cancer on CL (unitless)")         # Table S2: CL_GC   = 0.180
    e_tumtp_hodgkin_classical_cl <- -0.330; label("Exponential coefficient of classical Hodgkin lymphoma on CL (unitless)") # Table S2: CL_CHL  = -0.330

    # Covariate effects on Vc (File S1 / File S2 $PK). Power on WT;
    # exponential on female sex.
    e_wt_vc  <-  0.630; label("Power exponent of WT on Vc (unitless)")                                                  # Table S2: Vc_WTB = 0.630
    e_sex_vc <- -0.134; label("Exponential coefficient of female sex on Vc (unitless)")                                 # Table S2: Vc_SEX = -0.134

    # Covariate effects on subcutaneous bioavailability. These are applied
    # MULTIPLICATIVELY on the natural fraction scale, not exponentially and
    # not on the logit scale: File S1 / File S2 read
    # 'TVF1 = TVF1 * (F1_FEMALE)' and 'TVF1 = TVF1 * (F1_PS_1)' with no EXP().
    e_sex_fdepot      <- 0.859; label("Multiplicative factor of female sex on subcutaneous bioavailability (unitless)")  # Table S2: F_SEX = 0.859
    e_ecog_ge1_fdepot <- 1.07;  label("Multiplicative factor of ECOG_GE1 on subcutaneous bioavailability (unitless)")    # Table S2: F_PS  = 1.07

    # Time-varying clearance (sigmoid Emax in time since first dose;
    # File S1 / File S2 'CL_TIME = EXP(EMAX*TIME**HILL/(T50**HILL+TIME**HILL))').
    # Emax is the maximal fractional change in CL on the log scale: at
    # t >> T50, CL approaches CL_base * exp(-0.303) = 0.739 x baseline, a
    # 26.1 percent reduction. T50 is reported in hours and converted to days
    # (/ 24); the ratio t^HILL / (T50^HILL + t^HILL) is unit-invariant provided
    # t and T50 share a unit.
    cl_time_max  <- -0.303;     label("Maximal fractional change in CL, Emax_REF (unitless, log scale)")                # Table S2: Emax_REF = -0.303
    cl_t50       <-  1.40e3 / 24; label("Time at which the change in CL is 50%% of Emax (days)")                        # Table S2: T50 = 1.40e3 h; / 24 = 58.33 days
    cl_time_hill <-  2.82;      label("Hill / sigmoidicity exponent of time on CL (unitless)")                          # Table S2: HILL = 2.82

    # Inter-individual variability (Table S2 'Random effects'). Reported as
    # variance with the standard deviation in parentheses for diagonal
    # elements and as covariance with the correlation in parentheses for
    # off-diagonal elements; every parenthetical reproduces exactly
    # (sqrt(0.114) = 0.338; 0.0377 / (0.338 x 0.355) = 0.314), which confirms
    # the printed numbers are VARIANCES.
    #
    # CL and Vc are log-normal with one covariance; Vp is an independent
    # log-normal eta; Emax carries an INDEPENDENT ADDITIVE eta
    # (File S1 / File S2 'EMAX = AEMAX + ZEMAX'); Ka and the logit of F share
    # a second correlated block.
    etalcl + etalvc ~ c(0.114,
                        0.0377, 0.126)                                                                                  # Table S2: omega^2_CL = 0.114 (0.338), cov CL:Vc = 0.0377 (0.314), omega^2_Vc = 0.126 (0.355)
    etalvp ~ 0.235                                                                                                      # Table S2: omega^2_Vp = 0.235 (0.485)
    etacl_time_max ~ 0.0519                                                                                             # Table S2: omega^2_Emax = 0.0519 (0.228); additive on Emax
    etalka + etalogitfdepot ~ c(0.0955,
                                0.213, 0.862)                                                                           # Table S2: omega^2_Ka = 0.0955 (0.309), cov Ka:F = 0.213 (0.742), omega^2_F = 0.862 (0.928)

    # Residual error. Proportional only. The magnitude 0.204 is a standard
    # deviation, not a variance: both control streams fix the residual
    # variance with '$SIGMA 1 FIX' and build the error as
    # 'REWT = F*PERR + AERR' with 'AERR = 0' and 'Y = IPRED + REWT*EPS(1)',
    # so THETA(5) = PERR multiplies the prediction directly.
    propSd <- 0.204; label("Proportional residual error (fraction)")                                                    # Table S2: proportional = 0.204; File S2 $THETAP PERR = 0.204 FIX with $SIGMA 1 FIX
  })
  model({
    # Individual baseline CL and Vc with covariate adjustments (File S1 /
    # File S2 $PK). Reference subject is 80 kg, eGFR 90 mL/min, male, ECOG 0,
    # second-line NSCLC, so every covariate term equals 1 at the reference.
    cl_base <- exp(lcl + etalcl) *
      (WT   / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_sex_cl                     * SEXF) *
      exp(e_ecog_ge1_cl                * ECOG_GE1) *
      exp(e_tumtp_gastric_cl           * TUMTP_GASTRIC) *
      exp(e_tumtp_hodgkin_classical_cl * TUMTP_HODGKIN_CLASSICAL)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sex_vc * SEXF)

    # Time-varying clearance with additive IIV on Emax. t is time since the
    # first dose, matching NONMEM's TIME in the source control streams.
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    # Subcutaneous bioavailability. The covariate factors act on the natural
    # fraction scale and the IIV acts on the logit scale, reproducing the
    # source's two-step construction:
    #   TVF1   = AF1 * F1_FEMALE^SEXF * F1_PS_1^ECOG_GE1
    #   LOGITF1 = LOG(TVF1/(1-TVF1)) + ZF1
    #   F1     = EXP(LOGITF1)/(1+EXP(LOGITF1))
    # which is the paper's printed equation for F_i. Across the four
    # covariate cells TVF1 spans 0.646 - 0.805, so the logit is always finite.
    fdepot_tv <- expit(logitfdepot) *
      e_sex_fdepot^SEXF *
      e_ecog_ge1_fdepot^ECOG_GE1
    fdepot <- expit(log(fdepot_tv / (1 - fdepot_tv)) + etalogitfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Bioavailability applies only to the depot, so the dosing route is
    # expressed structurally: subcutaneous doses go to depot (absorbed with
    # ka and fraction fdepot), intravenous doses go to central (F = 1). This
    # reproduces the source's ADVAN4 setup, in which F1 governs compartment 1
    # only and the 'IF (ROUTEN.EQ.0)' branches that set TVKA = 1 and
    # TVF1 = 0.4 are inert placeholders for intravenous records dosed into
    # compartment 2.
    f(depot) <- fdepot

    # Dose in mg and volumes in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
