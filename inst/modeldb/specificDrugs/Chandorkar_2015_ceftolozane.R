Chandorkar_2015_ceftolozane <- function() {
  description <- paste(
    "Two-compartment population PK model for ceftolozane given as a 1-hour",
    "intravenous infusion, fitted to 5,048 plasma concentrations from 376",
    "adults pooled across 10 studies: five Phase 1 studies in healthy",
    "volunteers, three studies spanning mild to severe renal impairment, and",
    "two Phase 2 studies in patients with complicated urinary tract infection",
    "(cUTI) or complicated intra-abdominal infection (cIAI) (Chandorkar 2015).",
    "Elimination is first order. Baseline Cockcroft-Gault creatinine clearance",
    "acts on clearance through a power function, and body weight acts",
    "proportionally on central volume. Infection type shifts both clearance",
    "and central volume, and in cIAI patients the body-weight effect on",
    "central volume is switched off because the paper found no significant",
    "weight-volume correlation in that group. Between-subject variability is",
    "diagonal on CL and Vc only; the peripheral parameters carry none.",
    "Companion model to Chandorkar_2015_tazobactam, which the same paper",
    "fitted separately because the two drugs do not interact."
  )
  reference <- paste(
    "Chandorkar G, Xiao A, Mouksassi MS, Hershberger E, Krishna G.",
    "Population pharmacokinetics of ceftolozane/tazobactam in healthy",
    "volunteers, subjects with varying degrees of renal function and patients",
    "with bacterial infections. J Clin Pharmacol. 2015;55(2):230-239.",
    "doi:10.1002/jcph.395.",
    "All fixed-effect, random-effect and residual-error estimates are Table 3",
    "panel [A]. No supplement was deposited with the article (EuropePMC",
    "hasSuppl 'N'); Table 3A, the Results narrative and Figure 1A together",
    "report the complete final model.",
    sep = " "
  )
  vignette <- "Chandorkar_2015_ceftolozane_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are in mg and the paper reports plasma
  # concentrations in ug/mL, which is numerically identical to mg/L, so
  # amount(mg) / volume(L) lands directly in the reported concentration unit.
  compartmentData <- list(
    central     = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Baseline creatinine clearance estimated by the Cockcroft-Gault",
        "formula; raw mL/min, NOT body-surface-area normalized"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chandorkar 2015 Table 3A: 'CL (L/h), No infection, 5.11",
        "(2.15)*(CrCL/109)^0.715 (6.14)'. Reference value 109 mL/min; power",
        "exponent 0.715 (RSE 6.14%). The Discussion restates the same",
        "equation: 'CL of ceftolozane = 5.11 L/h * (CrCL/109)^0.715'.",
        "ASSAY FORM. Methods, 'Sources of Variability and Covariate",
        "Analysis': 'The CrCL was estimated using the Cockcroft-Gault",
        "formula ... where CrCL is creatinine clearance (mL/min), age is in",
        "years, WT is actual body weight (kg), and SCr is serum creatinine",
        "(mg/dL); for female subjects the value was multiplied by a factor of",
        "0.85.' Actual (not ideal or lean) body weight, and NOT normalized to",
        "1.73 m^2 -- do not interchange this column with a BSA-normalized",
        "eGFR.",
        "BASELINE (time-fixed) per subject. Results, 'Datasets': 'Baseline",
        "CrCL was used to describe renal function since serum creatinine was",
        "stable across the short treatment duration, with a median value of",
        "the actual changes (increase or decrease) of approximately 5% and a",
        "median value of the absolute changes of <15%.'",
        "RANGE FITTED. Table 2 gives cohort means 101.0 mL/min (range 19-215)",
        "without infection and 97.4 (range 41-309) with infection; Figure 1A",
        "labels the pooled observed range as 19.1-308.5 mL/min. The renal",
        "impairment CATEGORIES quoted throughout the paper (normal >=90; mild",
        ">=50 to <90; moderate >=30 to <50; severe 15 to <30 mL/min) are",
        "descriptive strata only -- the model itself uses the continuous",
        "column, with no category indicators."
      ),
      source_name        = "CrCL"
    ),
    WT = list(
      description        = "Actual total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chandorkar 2015 Table 3A: 'Vc (L), No infection, 11.4",
        "(2.70)*(weight/74)'. Reference value 74 kg. The exponent is",
        "STRUCTURALLY 1, not estimated: no RSE is attached to a weight",
        "exponent anywhere in Table 3A (contrast the CrCL row, where the",
        "6.14% RSE is attached to the 0.715), and the Results state 'the Vc",
        "changed proportionally (linearly) with body weight in subjects",
        "without cIAI'. The Discussion confirms the slope: 'Vc would be",
        "expected to change by about 20% for every 20% change (increase or",
        "decrease) in body weight'. Encoded as e_wt_vc <- fixed(1).",
        "SWITCHED OFF IN cIAI PATIENTS. Table 3A gives the cIAI Vc row as a",
        "bare 'x1.59 (12.3)' with no '*(weight/74)' term, unlike the",
        "no-infection and cUTI rows which both carry it. Results: 'in cIAI",
        "patients, there was no significant correlation between Vc and body",
        "weight given the large observed variability.' Discussion: 'except in",
        "cIAI patients, where volume was not a function of weight'. The",
        "weight exponent is therefore multiplied by (1 - DIS_CIAI) in",
        "model().",
        "Weight also entered the clearance screen -- Results: 'Covariate",
        "analysis showed that both CL and Vc increased with body weight' --",
        "but Table 3A retains it on Vc only, so no CL term is carried.",
        "Cohort means 73.5 kg (range 49-106) without infection and 79.6",
        "(range 43-173) with infection (Table 2). Baseline (time-fixed)."
      ),
      source_name        = "weight"
    ),
    DIS_CUTI = list(
      description        = "Complicated urinary tract infection cohort indicator (1 = cUTI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer or renal-impairment subject with no infection; DIS_CIAI also 0)",
      notes              = paste(
        "Chandorkar 2015 Table 3A rows 'CL (L/h), With cUTI, x1.21 (24.6)'",
        "and 'Vc (L), With cUTI, x1.21 (30.1)*(weight/74)'. Both printed",
        "factors are exp(beta) under the paper's own stated parameterisation",
        "(Methods: categorical covariates entered as 'a linear model with an",
        "exponentiated factor relative to the reference'), so the ini()",
        "coefficients are log(1.21).",
        "The cUTI factor on CL is independently confirmed inside the Figure",
        "1A tornado panel, which prints the cUTI bar as spanning relative CL",
        "1 to 1.21.",
        "Mutually exclusive with DIS_CIAI; both 0 is the uninfected",
        "reference. 73 of the 376 subjects had cUTI (Table 1, Umeh 2010 /",
        "NCT00921024), and those patients received ceftolozane 1000 mg q8h",
        "WITHOUT tazobactam -- which is why the companion",
        "Chandorkar_2015_tazobactam model carries no cUTI term at all rather",
        "than having omitted one.",
        "Unlike Vc, the cUTI Vc row RETAINS the body-weight term; only cIAI",
        "switches it off. See the WT notes."
      ),
      source_name        = "cUTI"
    ),
    DIS_CIAI = list(
      description        = "Complicated intra-abdominal infection cohort indicator (1 = cIAI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer or renal-impairment subject with no infection; DIS_CUTI also 0)",
      notes              = paste(
        "Chandorkar 2015 Table 3A rows 'CL (L/h), With cIAI, x1.22 (22.5)'",
        "and 'Vc (L), With cIAI, x1.59 (12.3)'. Printed as exp(beta) factors,",
        "so the ini() coefficients are log(1.22) and log(1.59).",
        "The cIAI factor on CL is independently confirmed inside the Figure",
        "1A tornado panel, which prints the cIAI bar as spanning relative CL",
        "1 to 1.22. The Vc factor is confirmed by the Discussion's worked",
        "example: 'with a typical body weight of 74 kg, Vc would be about 11,",
        "14, and 18 L in healthy subjects, cUTI patients, and cIAI patients'",
        "-- 11.4 * 1.59 = 18.1 L, and 'Vc was about 30% different between",
        "these two patient groups (13.8 L at 74 kg body weight for cUTI and",
        "18.2 L for cIAI)'.",
        "BEWARE the Discussion's competing sentence 'the cIAI effect was",
        "increasing it by factor of 1.50-fold'. That 1.50 is inconsistent",
        "with every other statement in the paper -- Table 3A's 1.59, the",
        "'about 18 L' worked example and the '18.2 L' figure both require",
        "1.59 (11.4 * 1.50 = 17.1 L, not 18.2). The table value is used.",
        "THIS INDICATOR ALSO GATES THE BODY-WEIGHT EFFECT on Vc: in cIAI",
        "patients Vc is not a function of weight. See the WT notes for the",
        "textual evidence and the model() encoding.",
        "Mutually exclusive with DIS_CUTI. 77 of the 376 subjects had cIAI",
        "(Table 1, Lucasti 2014 / NCT01147640), receiving",
        "ceftolozane/tazobactam 1000/500 mg q8h."
      ),
      source_name        = "cIAI"
    )
  )

  # Covariates the paper screened but did not retain in the final ceftolozane
  # model. Documented so the provenance of the covariate screen survives,
  # without declaring covariateData entries that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as an intrinsic covariate (Methods) and found directional",
        "but not retained. Results: 'A small negative trend between age and",
        "CL was also observed but it was not clinically meaningful.' No age",
        "coefficient appears in Table 3A. Cohort range 18-86 years."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and rejected. Results: 'Other covariates such as race, sex,",
        "dose level, and drug-drug interaction did not significantly affect",
        "CL or Vc of ceftolozane.' Note that sex nonetheless enters the model",
        "INDIRECTLY through the Cockcroft-Gault CRCL column, which multiplies",
        "by 0.85 for female subjects."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and rejected (Results: race 'did not significantly affect",
        "CL or Vc'). The cohort was 82.7% white without infection and 96.7%",
        "with infection (Table 2), so the analysis had little power to",
        "resolve a race effect."
      )
    ),
    DOSE_CEFTOLOZANE_MG = list(
      description = "Administered ceftolozane dose",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "Screened as an extrinsic covariate to test for dose-dependent",
        "(nonlinear) PK and rejected, supporting the linear structural model.",
        "Results: dose level 'did not significantly affect CL or Vc'.",
        "Discussion: 'the PK of ceftolozane/tazobactam is dose-proportional",
        "and linear'. Doses spanned 250-3000 mg (Table 1)."
      )
    ),
    CONMED_TAZOBACTAM = list(
      description = "Co-administration of tazobactam indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as the drug-drug-interaction covariate and rejected.",
        "Discussion: 'similar to previous observations, no drug-drug",
        "interaction was observed between ceftolozane and tazobactam and the",
        "PK profile of ceftolozane was unaffected by administration of",
        "tazobactam.' This is the finding that licenses fitting -- and",
        "extracting -- ceftolozane and tazobactam as two independent models."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 376L,
    n_studies        = 10L,
    n_concentrations = 5048L,
    age_range        = "18-86 years",
    age_median       = "means 44.7 years (no infection) and 53.5 years (infection); medians not reported",
    weight_range     = "43-173 kg",
    weight_median    = "means 73.5 kg (no infection) and 79.6 kg (infection); medians not reported",
    sex_female_pct   = 43.6,
    race_ethnicity   = c(White = 88.3, Other = 11.7),
    disease_state    = paste(
      "Pooled healthy adults, adults with mild to severe renal impairment,",
      "and hospitalized patients with complicated urinary tract infection or",
      "complicated intra-abdominal infection. Table 2: 226 of 376 subjects",
      "had no infection and 150 (39.9%) had cUTI or cIAI; 121 (32.2%) were",
      "renally impaired."
    ),
    dose_range       = paste(
      "All doses given as 1-hour intravenous infusions, ceftolozane alone or",
      "with tazobactam at a fixed 2:1 ratio. Healthy volunteers received",
      "single or multiple (q8h or q12h) doses of ceftolozane 250, 500, 1000,",
      "1500 or 2000 mg, or ceftolozane/tazobactam 500/250, 1000/500, 1500/750,",
      "2000/1000 or 3000/1500 mg. Renal-impairment subjects received a single",
      "1000 mg (mild/moderate) or 500 mg (severe) ceftolozane dose. cUTI",
      "patients received ceftolozane 1000 mg q8h alone; cIAI patients received",
      "ceftolozane/tazobactam 1000/500 mg q8h (Table 1)."
    ),
    renal_function   = paste(
      "Spans severe renal impairment to augmented clearance. Estimated CrCL",
      "means 101.0 mL/min (range 19-215) without infection and 97.4 mL/min",
      "(range 41-309) with infection; Figure 1A gives the pooled observed",
      "range as 19.1-308.5 mL/min. Category counts (Table 2, no infection /",
      "infection): normal 186/69, mild 28/78, moderate 6/3, severe 6/0. No",
      "end-stage-renal-disease or dialysis subjects were enrolled, so the",
      "model carries no information below about 15 mL/min."
    ),
    bmi_range        = "17-56 kg/m^2 (means 25.8 without infection, 27.3 with infection)",
    pkpd_target      = paste(
      "Not fitted in this paper. Discussion: 'the therapeutic efficacy of",
      "ceftolozane is best correlated with the percentage of time the plasma",
      "drug concentration exceeds the MIC for the target organism (%T>MIC)'.",
      "The paper positions this PK model as the input to a later probability",
      "-of-target-attainment analysis rather than performing one."
    ),
    notes            = paste(
      "sex_female_pct is the pooled ceftolozane data set: (97 + 67) / 376 =",
      "43.6% (Table 2). race_ethnicity is likewise pooled white (187 + 145) /",
      "376 = 88.3%; the paper reports only the white percentage, so the",
      "remainder is recorded as Other rather than broken out.",
      "Estimation was first-order conditional estimation-extended least",
      "squares (FOCE-ELS) in Phoenix NLME 1.2 (Certara), NOT NONMEM.",
      "Model evaluation used goodness-of-fit diagnostics, a visual predictive",
      "check with 1,000 replicates, and 1,000-sample nonparametric bootstrap",
      "resampling, with bootstrap-versus-final differences under 5%.",
      "Concentrations below the limit of quantitation were treated as missing",
      "with no substitution. Twenty-five samples from 20 subjects had",
      "CWRES > 4 and were retained; excluding them shifted PK parameters by",
      "-0.2% to 6.7%."
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Chandorkar 2015 Table 3A, "Population
    # estimates (RSE %)" column. Values refer to the reference subject:
    # CrCL = 109 mL/min, weight = 74 kg, no infection.
    # =====================================================================
    lcl <- log(5.11); label("Ceftolozane clearance at CrCL = 109 mL/min, no infection (L/h)")            # Table 3A: CL, No infection = 5.11 L/h (RSE 2.15%)
    lvc <- log(11.4); label("Ceftolozane central volume of distribution at weight = 74 kg, no infection (L)") # Table 3A: Vc, No infection = 11.4 L (RSE 2.70%)
    lq  <- log(1.19); label("Ceftolozane inter-compartmental clearance (L/h)")                           # Table 3A: CL2 = 1.19 L/h (RSE 2.24%)

    lvp <- fixed(log(2.88)); label("Ceftolozane peripheral volume of distribution (L)")
    # Table 3A: Vp = "2.88 (fixed)" -- the only parameter in the table whose
    # parenthetical is the word "fixed" rather than an RSE, so it is held
    # constant rather than estimated. The Discussion rounds it: "volume of
    # distribution in the peripheral compartment was about 3 L".

    # =====================================================================
    # COVARIATE EFFECTS.
    #
    # Methods, "Sources of Variability and Covariate Analysis": "Covariates
    # were introduced in a multiplicative order using a power model
    # standardized by the median for continuous covariates and a linear
    # model with an exponentiated factor relative to the reference for
    # categorical covariates."
    #
    # So the continuous covariates enter as (cov / ref)^exponent and the
    # categorical ones as exp(beta), which is what Table 3A prints as
    # "x1.21" etc. The ini() coefficients below are therefore the logs of
    # the printed factors.
    # =====================================================================
    e_crcl_cl <- 0.715; label("Power exponent of baseline CrCL on ceftolozane CL, CRCL/109 (unitless)")
    # Table 3A: 5.11*(CrCL/109)^0.715, RSE 6.14% on the exponent.
    # Independently confirmed against the Figure 1A tornado panel, which
    # prints relative CL at six CrCL values; (CrCL/109)^0.715 reproduces all
    # six exactly: 19.1 -> 0.29, 15 -> 0.24, 30 -> 0.40, 50 -> 0.57,
    # 90 -> 0.87, 308.5 -> 2.10.

    e_wt_vc <- fixed(1); label("Power exponent of body weight on ceftolozane Vc, WT/74 (unitless)")
    # Table 3A prints the Vc rows as "*(weight/74)" with NO exponent and no
    # RSE attached, and the Results say Vc "changed proportionally
    # (linearly) with body weight". Structurally 1, not estimated -- hence
    # fixed(). See covariateData$WT. In cIAI patients this exponent is
    # switched off entirely in model().

    e_cuti_cl <- log(1.21); label("Effect of cUTI on ceftolozane CL (log multiplicative factor)")  # Table 3A: CL, With cUTI = x1.21 (RSE 24.6%); Figure 1A tornado bar 1 to 1.21
    e_ciai_cl <- log(1.22); label("Effect of cIAI on ceftolozane CL (log multiplicative factor)")  # Table 3A: CL, With cIAI = x1.22 (RSE 22.5%); Figure 1A tornado bar 1 to 1.22
    e_cuti_vc <- log(1.21); label("Effect of cUTI on ceftolozane Vc (log multiplicative factor)")  # Table 3A: Vc, With cUTI = x1.21 (RSE 30.1%)
    e_ciai_vc <- log(1.59); label("Effect of cIAI on ceftolozane Vc (log multiplicative factor)")  # Table 3A: Vc, With cIAI = x1.59 (RSE 12.3%)

    # =====================================================================
    # BETWEEN-SUBJECT VARIABILITY -- Table 3A, "BSV % (RSE %)" column.
    #
    # Results: "A two-compartmental structural model with a diagonal
    # variance (omega) of CL, Vc, peripheral volume of distribution (Vp),
    # and peripheral clearance (CL2) fixed to a value of 0 provided the best
    # data fit." Diagonal, and the CL2 / Vp rows read "Fixed at 0", so only
    # CL and Vc carry IIV.
    #
    # SCALE. Methods: "A variance component, which assumed a log-normal
    # distribution of PK parameters, was used to characterize the
    # between-subject variability". The column is a %CV on that log-normal,
    # so omega^2 = log(1 + CV^2). The reported RSEs corroborate that the
    # column is an SD-like quantity rather than a variance: for a variance
    # estimated on 376 subjects the asymptotic RSE is about sqrt(2/376) =
    # 7.3%, against about 3.6% for the corresponding SD, and the printed
    # RSEs are 3.94% and 4.50%.
    # =====================================================================
    etalcl ~ 0.103369  # BSV on CL 33.0% CV (RSE 3.94%, shrinkage 3.5%): log(1 + 0.330^2)
    etalvc ~ 0.147043  # BSV on Vc 39.8% CV (RSE 4.50%, shrinkage 8.3%): log(1 + 0.398^2)

    # =====================================================================
    # RESIDUAL VARIABILITY -- Table 3A. Results: "The residual variability
    # was found to be composite (both proportional and additive)."
    #
    # The paper supplies its own arithmetic check on both magnitudes at
    # once: "For a fitted ceftolozane concentration of 100 ug/mL, the total
    # residual error would be 16.85 ug/mL" -- i.e. 0.168 * 100 + 0.0524 =
    # 16.85, which pins the proportional term as a fraction and the additive
    # term as ug/mL.
    # =====================================================================
    propSd <- 0.168;  label("Proportional residual error (fraction)")  # Table 3A: Proportional error = 16.8% (RSE 11.8%)
    addSd  <- 0.0524; label("Additive residual error (ug/mL)")         # Table 3A: Additional error = 0.0524 ug/mL (RSE 8.07%)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # Clearance (Table 3A / Discussion "CL of ceftolozane = 5.11 L/h *
    # (CrCL/109)^0.715 with a multiplicative factor of 1.21 for patients
    # with cUTI and 1.22 for patients with cIAI"):
    #
    #   CL = 5.11 * (CrCL/109)^0.715 * 1.21^cUTI * 1.22^cIAI
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl + e_cuti_cl * DIS_CUTI + e_ciai_cl * DIS_CIAI) *
      (CRCL / 109)^e_crcl_cl

    # Central volume (Table 3A). The three rows of the Vc block are
    #
    #   no infection : 11.4 * (weight/74)
    #   cUTI         : 11.4 * 1.21 * (weight/74)
    #   cIAI         : 11.4 * 1.59              <- no weight term
    #
    # so the weight exponent is gated by (1 - DIS_CIAI). Multiplying the
    # exponent rather than branching keeps the expression continuous and
    # finite for every positive weight, and reduces to exactly 11.4 * 1.59
    # when DIS_CIAI is 1 because (WT/74)^0 is 1.
    vc <- exp(lvc + etalvc + e_cuti_vc * DIS_CUTI + e_ciai_vc * DIS_CIAI) *
      (WT / 74)^(e_wt_vc * (1 - DIS_CIAI))

    q  <- exp(lq)
    vp <- exp(lvp)

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 3. Two-compartment disposition. Ceftolozane is given as a zero-order
    #    1-hour intravenous infusion straight into the central compartment;
    #    the infusion is expressed in the event table (amt + dur or rate),
    #    not here.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Observation. Dose in mg, volumes in L -> mg/L, numerically
    #    identical to the ug/mL the paper reports. These are TOTAL plasma
    #    concentrations; the paper's %T>MIC efficacy target is defined on
    #    free concentrations, but no ceftolozane unbound fraction is
    #    reported here, so no conversion is carried.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
