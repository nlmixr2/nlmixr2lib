Ousey_2026_plozasiran <- function() {
  description <- paste(
    "Cascading kinetic-pharmacodynamic (K-PD) indirect-response population",
    "PD model for plozasiran (an APOC3-targeting GalNAc-conjugated siRNA) in",
    "adult patients with familial chylomicronemia syndrome (Ousey 2026,",
    "Phase 3 PALISADE). No plasma PK is used: the subcutaneous dose enters a",
    "volume-less biophase (liver) compartment `depot_kpd` in mg and is",
    "eliminated first-order at kel (t1/2 = 48.1 days), which is the",
    "rate-limiting step sustaining the months-long PD effect. Plozasiran in",
    "the biophase inhibits zero-order serum APOC3 synthesis (Imax fixed to",
    "100%, IC50 = 2.185 mg of biophase amount); serum APOC3 in turn inhibits",
    "first-order serum triglyceride clearance (Imax fixed to 100%,",
    "half-maximal APOC3 level 17.4 mg/dL). Body mass index acts on kel and on",
    "the APOC3-on-TG potency in opposing directions, and background",
    "triglyceride-lowering therapy potentiates the APOC3-on-TG step.",
    "Two simultaneous endpoints: serum APOC3 and serum triglycerides, both in",
    "mg/dL."
  )
  reference <- paste(
    "Ousey J, Gosselin NH, Ta A, Shi J. Population Pharmacodynamic Modeling",
    "of Plozasiran for Treatment of Familial Chylomicronemia Syndrome.",
    "J Clin Pharmacol. 2026;66(4). doi:10.1002/jcph.70190.",
    "Model equations from the Supporting Information (Supplementary",
    "Equations); parameter estimates from Table 2."
  )
  vignette <- "Ousey_2026_plozasiran"
  units <- list(time = "day", dosing = "mg", concentration = "mg/dL")

  # The dose does not enter `depot` or `central`; declare the K-PD biophase as
  # the dosing target explicitly (Ousey 2026 Results: "the full dose of
  # plozasiran is assumed to be transferred instantaneously into the biophase
  # compartment").
  dosing <- c("depot_kpd")

  # `apoc3` and `tg` are the paper's own symbols for the two serum-biomarker
  # states. They are not yet canonical compartments in
  # inst/references/compartment-names.md (which registers the sibling lipid
  # biomarkers `ldl`, `nefa` and `hc24`); declared here as paper-specific
  # pending a future canonicalisation of the two names.
  paper_specific_compartments <- c("apoc3", "tg")

  compartmentData <- list(
    depot_kpd = list(
      analyte  = "plozasiran",
      units    = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    apoc3 = list(
      analyte  = "apolipoprotein C-III",
      units    = "mg/dL",
      specimen = "serum",
      verified = TRUE
    ),
    tg = list(
      analyte  = "triglycerides",
      units    = "mg/dL",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    APOC3 = list(
      description        = "Baseline (pre-first-dose) fasting serum apolipoprotein C-III concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (the pre-dose measurement), used as an",
        "a-priori continuous covariate on the estimated APOC3 baseline:",
        "(APOC3 / 34)^e_apoc3_rbase_apoc3. Ousey 2026 Methods: 'The",
        "subject-specific baselines of APOC3 and triglycerides seen in the",
        "study were incorporated a priori as continuous covariates on",
        "BASE_APOC3 and BASE_TG, respectively.' The centring value printed in",
        "Table 2 is 34 mg/dL, which is close to but not equal to the Table 1",
        "cohort median of 35.4 mg/dL; the printed 34 is used here (see",
        "vignette Errata). Cohort median 35.4 mg/dL [10.0, 88.0] (Table 1).",
        "Assay LLOQ 0.94 mg/dL, ULOQ 88 mg/dL; BLQ imputed as LLOQ/2 and AQL",
        "imputed as ULOQ in the analysis dataset."
      ),
      source_name        = "Baseline APOC3"
    ),
    TRIG = list(
      description        = "Baseline (pre-first-dose) fasting serum triglyceride concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject, used as an a-priori continuous covariate on",
        "the estimated triglyceride baseline: (TRIG / 2048)^e_trig_rbase_tg.",
        "The centring value printed in Table 2 is 2048 mg/dL versus a Table 1",
        "cohort median of 2044 mg/dL; the printed 2048 is used here (see",
        "vignette Errata). Cohort median 2044 mg/dL [747, 6597] (Table 1).",
        "Units are mg/dL in this model, NOT mmol/L."
      ),
      source_name        = "Baseline TG"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. The only demographic covariate retained in",
        "the final model, acting on TWO parameters in opposing directions:",
        "(BMI / 25)^1.39 on kel (higher BMI implies faster biophase",
        "elimination, hence less APOC3 reduction) and (BMI / 25)^-7.17 on the",
        "APOC3-on-TG half-maximal level (higher BMI implies greater TG",
        "sensitivity for a given APOC3 reduction). Ousey 2026 concludes the",
        "net clinical impact is 'muted and clinically unimportant' and that a",
        "fixed dose is appropriate regardless of BMI or body weight.",
        "Centring value 25 kg/m^2 equals the Table 1 cohort median.",
        "Cohort median 25.0 kg/m^2 [18.5, 35.9] (Table 1)."
      ),
      source_name        = "BMI"
    ),
    CONMED_TGLOWER = list(
      description        = "Indicator of one or more stable background triglyceride-lowering therapies",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no background triglyceride-lowering therapy)",
      notes              = paste(
        "Ousey 2026 Table 1 footnote a: 'Background TG-lowering therapies are",
        "defined as one or any combination of: fibrates, Vascepa or omega-3",
        "fatty acids or fish oil, and statins.' 55/75 patients (73.3%) were on",
        ">= 1 such therapy. Applied multiplicatively on the exponential scale",
        "to the APOC3-on-TG half-maximal level:",
        "ki50_apoc3 = 17.4 * exp(-1.79 * CONMED_TGLOWER). The paper pins the",
        "reference category unambiguously: 'IC50,APOC3 is estimated to be 2.91",
        "and 17.4 mg/dL with and without concomitant TG-lowering therapy,",
        "respectively', and 17.4 * exp(-1.79) = 2.905. Concomitant therapy",
        "therefore potentiates (rather than antagonises) the plozasiran effect.",
        "The paper reports only the composite indicator and explicitly declines",
        "to attribute the effect to any single class: 'there was no definitive",
        "evidence to indicate which, if any, specific TG-lowering therapy was",
        "primarily responsible'. Per-model class composition is documented here",
        "rather than in the canonical register."
      ),
      source_name        = "TG-lowering therapy"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search but not retained; highly correlated with the retained BMI. Cohort median 70.2 kg [43.5, 118] (Table 1)."
    ),
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (including a geriatric vs non-geriatric < 65 / >= 65 split) but not retained. Cohort median 44.0 years [22, 76] (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Cohort 38/75 (50.7%) female (Table 1)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (White vs non-White and Asian vs non-Asian) but not retained. Cohort 16/75 (21.3%) Asian (Table 1)."
    ),
    EGFR = list(
      description = "Baseline estimated glomerular filtration rate",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened, together with a categorical normal/mild/moderate renal-impairment classification, but not retained. Cohort median 113 mL/min [43.0, 640] (Table 1)."
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as part of the liver-function panel but not retained. Cohort median 0.530 mg/dL [0.240, 1.17] (Table 1). All 75 patients had normal hepatic function by Child-Pugh and NCI-ODWG criteria, so hepatic impairment could not be evaluated."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as part of the liver-function panel but not retained. Cohort median 22.0 U/L [12.0, 43.0] (Table 1)."
    ),
    HBA1C = list(
      description = "Baseline hemoglobin A1c",
      units       = "%",
      type        = "continuous",
      notes       = "Screened as a < 6.5% vs >= 6.5% split but not retained. Values not tabulated in Table 1."
    ),
    DIS_FCS_GENETIC = list(
      description = "Genetically confirmed (vs clinically diagnosed) FCS indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained -- a headline negative finding of the paper, supporting plozasiran use in both genetically confirmed and clinically diagnosed FCS. Cohort 41/75 (54.7%) genetically confirmed (Table 1). Not a registered canonical covariate; documented here only as part of the covariate screen."
    ),
    COUNTRY_JAPAN = list(
      description = "Study conducted in Japan (vs ex-Japan) indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "The only extrinsic factor screened; not retained. Only 3 Japanese subjects were in the active-treatment posterior analysis. Not a registered canonical covariate; documented here only as part of the covariate screen."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 75L,
    n_studies      = 1L,
    n_observations = 2106L,
    age_range      = "22-76 years",
    age_median     = "44.0 years",
    weight_range   = "43.5-118 kg",
    weight_median  = "70.2 kg",
    bmi_range      = "18.5-35.9 kg/m^2",
    bmi_median     = "25.0 kg/m^2",
    sex_female_pct = 50.7,
    race_ethnicity = c(White = 73.3, Asian = 21.3, `American Indian or Alaska Native` = 1.3, Other = 4.0),
    disease_state  = paste(
      "Familial chylomicronemia syndrome (FCS) in adults: 41/75 (54.7%)",
      "genetically confirmed and 34/75 (45.3%) clinically diagnosed.",
      "Median baseline fasting serum APOC3 35.4 mg/dL and triglycerides",
      "2044 mg/dL, characteristic of FCS. All 75 patients had normal hepatic",
      "function (Child-Pugh and NCI-ODWG); renal function was normal in 59",
      "(78.7%), mildly impaired in 12 (16.0%) and moderately impaired in 4",
      "(5.3%). 55/75 (73.3%) were on >= 1 stable background",
      "triglyceride-lowering therapy (fibrates, omega-3 fatty acids / fish",
      "oil, statins)."
    ),
    dose_range     = paste(
      "Four Q3M subcutaneous doses of 25 mg (n = 26, 35%) or 50 mg (n = 24,",
      "32%) plozasiran, or volume-matched placebo (n = 25, 33%), over the",
      "12-month double-blinded period. Vial formulation (200 mg/mL) was used",
      "in the trial; the commercial prefilled syringe (50 mg/mL) may deliver",
      "up to 22% higher relative bioavailability, i.e. an effective 30.5 mg",
      "for a nominal 25 mg dose."
    ),
    regions        = "Multinational (PALISADE, NCT05089084), including Japan.",
    renal_function = "Normal 59 (78.7%), mild 12 (16.0%), moderate 4 (5.3%) by eGFR (normal >= 90, mild 60-90, moderate 30-60 mL/min).",
    hepatic_function = "Normal in all 75 patients by Child-Pugh and NCI-ODWG criteria; hepatic impairment therefore not evaluable as a PD covariate.",
    notes          = paste(
      "Demographics from Ousey 2026 Table 1. Only the 12-month randomized,",
      "double-blinded period of PALISADE was used for model development; the",
      "open-label extension was excluded. 1044 serum APOC3 observations",
      "(87.1% within quantification limits, 11.5% BLQ imputed to 0.47 mg/dL,",
      "1.4% above ULOQ imputed to 88 mg/dL) and 1062 serum triglyceride",
      "observations (none BLQ or AQL) were included, giving 2106 total. Fitted",
      "in NONMEM 7.4 with FOCEI; final-model condition number 105."
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters. All from Ousey 2026 Table 2 ("Population
    # Parameter Estimates for the Final PD Model"). Model equations are from
    # the Supporting Information (Supplementary Equations):
    #
    #   dCe/dt    = -Ce * Ke
    #   dAPOC3/dt = Ksyn_APOC3 * DrugEffect_ploz - Kdeg_APOC3 * APOC3
    #   DrugEffect_ploz = 1 - Imax_ploz * Ce^gamma / (IC50_ploz^gamma + Ce^gamma)
    #   dTG/dt    = Ksyn_TG - Kdeg_TG * TG *
    #                 (1 - Imax_APOC3 * APOC3^gamma / (IC50_APOC3^gamma + APOC3^gamma))
    #
    # The Hill coefficient gamma was tested and rejected ("Attempts to
    # introduce a more complex sigmoidal Imax/IC50 relationship via Hill
    # coefficients did not statistically improve the model fit and were
    # rejected"), so gamma == 1 throughout and is not carried as a parameter.
    #
    # UNITS NOTE: Table 2 reports kel in /day but ksyn_APOC3 and kdeg_TG per
    # HOUR. The `ini()` values below are the literal Table 2 numbers in the
    # paper's own units; the 24 h/day conversion is applied explicitly at the
    # top of model() so that every number here is directly auditable against
    # Table 2.
    # --------------------------------------------------------------------

    lkel <- log(0.0144); label("Plozasiran first-order elimination rate constant from the K-PD biophase (1/day)") # Table 2: ke = 0.0144 /day, RSE 8.82%, 95% CI 0.0119-0.0169 (t1/2 = 48.1 days)

    # --- Serum APOC3 arm ---
    lrbase_apoc3 <- log(30.9); label("Typical baseline serum APOC3 concentration (mg/dL)")                       # Table 2: BASE_APOC3 = 30.9 mg/dL, RSE 5.34%, 95% CI 27.7-34.2
    lksyn_apoc3  <- fixed(log(1.13)); label("Zero-order serum APOC3 synthesis rate (mg/dL/h)")                   # Table 2: ksyn,APOC3 = 1.13 mg/dL/h, held constant "to ensure model stability and achieve a successful covariance step ... based on literature and prior estimates from K-PD models established for Phase 2 clinical studies"
    limax_ploz   <- fixed(log(1)); label("Maximum fractional inhibition of APOC3 synthesis by plozasiran (unitless)") # Table 2: Imax,ploz = 100%; "the estimates for Imax,ploz and Imax,APOC3 approached the upper boundary of 100%, and were therefore fixed to 100% for model simplification"
    lic50_ploz   <- log(2.185); label("Plozasiran biophase amount giving half-maximal inhibition of APOC3 synthesis (mg)") # Table 2: IC50,ploz = 2.185 mg, RSE 18.3%, 95% CI 1.401-2.969

    # --- Serum triglyceride arm ---
    lrbase_tg   <- log(1940); label("Typical baseline serum triglyceride concentration (mg/dL)")                 # Table 2: BASE_TG = 1940 mg/dL, RSE 7.39%, 95% CI 1659-2221
    lkdeg_tg    <- fixed(log(0.724)); label("First-order serum triglyceride degradation rate constant, uninhibited (1/h)") # Table 2: kdeg,TG = 0.724 /h, held constant alongside ksyn,APOC3 for the same stated reason
    limax_apoc3 <- fixed(log(1)); label("Maximum fractional inhibition of triglyceride clearance by APOC3 (unitless)") # Table 2: Imax,APOC3 = 100%; fixed for the same reason as Imax,ploz
    lki50_apoc3 <- log(17.4); label("Serum APOC3 concentration giving half-maximal inhibition of triglyceride clearance (mg/dL)") # Table 2: IC50,APOC3 = 17.4 mg/dL, RSE 59.9%, 95% CI -3.02-37.8 (the only major parameter not estimated to <= 18.3% RSE)

    # --------------------------------------------------------------------
    # Covariate effects (Ousey 2026 Table 2, "Covariate effects" block).
    # Continuous covariates use a power function centred on the (median)
    # value printed in Table 2; the categorical covariate is a multiplicative
    # exp(theta) factor. Reference/centring values are the ones PRINTED in
    # Table 2, which for the two baseline biomarkers differ slightly from the
    # Table 1 cohort medians (34 vs 35.4 mg/dL; 2048 vs 2044 mg/dL) -- see
    # vignette Errata.
    # --------------------------------------------------------------------
    e_apoc3_rbase_apoc3      <- 0.888; label("Power exponent on (APOC3 / 34) for baseline serum APOC3 (unitless)")                    # Table 2: Baseline APOC3 on BASE_APOC3, (baseline APOC3/34)^theta, theta = 0.888, RSE 7.25%, 95% CI 0.762-1.01
    e_trig_rbase_tg          <- 1.05;  label("Power exponent on (TRIG / 2048) for baseline serum triglycerides (unitless)")           # Table 2: Baseline TG on BASE_TG, (baseline TG/2048)^theta, theta = 1.05, RSE 9.15%, 95% CI 0.864-1.24
    e_bmi_kel                <- 1.39;  label("Power exponent on (BMI / 25) for the biophase elimination rate constant (unitless)")    # Table 2: BMI on ke, (BMI/25)^theta, theta = 1.39, RSE 29.8%, 95% CI 0.579-2.20
    e_bmi_ki50_apoc3         <- -7.17; label("Power exponent on (BMI / 25) for the APOC3-on-triglyceride half-maximal level (unitless)") # Table 2: BMI on IC50,APOC3, (BMI/25)^theta, theta = -7.17, RSE 33.8%, 95% CI -11.9 - -2.42
    e_conmed_tglower_ki50_apoc3 <- -1.79; label("Log-scale effect of background triglyceride-lowering therapy on the APOC3-on-triglyceride half-maximal level (unitless)") # Table 2: TG-lowering therapy on IC50,APOC3, exp(theta), theta = -1.79, RSE 35.6%, 95% CI -3.05 - -0.542; 17.4 * exp(-1.79) = 2.905, matching the paper's stated 2.91 mg/dL on-therapy value

    # --------------------------------------------------------------------
    # Inter-individual variability (Ousey 2026 Table 2, "Between subject
    # variability (Standard deviation, CV%)" block). The paper reports the
    # log-scale STANDARD DEVIATION with the log-normal CV% in parentheses;
    # nlmixr2 `ini()` takes VARIANCES, so each entry below is SD^2. The
    # reported CV% values are reproduced exactly by
    # CV = sqrt(exp(SD^2) - 1), confirming the SD-on-log-scale reading:
    #   0.316 -> 32.4%, 0.176 -> 17.7%, 0.984 -> 128%, 0.337 -> 34.6%,
    #   0.765 -> 89.2%.
    # --------------------------------------------------------------------
    etalkel        ~ 0.099856                                                                       # Table 2: On ke, SD 0.316 (32.4% CV), RSE 17.9%, 95% CI 0.205-0.427; shrinkage 3.8%
    etalic50_ploz  ~ 0.968256                                                                       # Table 2: On IC50,ploz, SD 0.984 (128% CV), RSE 17.1%, 95% CI 0.654-1.31; shrinkage 4.4%
    etalki50_apoc3 ~ 0.585225                                                                       # Table 2: On IC50,APOC3, SD 0.765 (89.2% CV), RSE 23.3%, 95% CI 0.416-1.11; shrinkage 35.5%

    # Correlated baseline IIV. Table 2 reports "Correlation BASE_APOC3,
    # BASE_TG = 0.150 (15.1%)" with RSE 129% and 95% CI -0.229-0.530, i.e. a
    # CORRELATION COEFFICIENT of 0.150 (the parenthetical 15.1% is the same
    # quantity as a percentage, not a CV). The off-diagonal covariance is
    # therefore r * SD_apoc3 * SD_tg = 0.150 * 0.176 * 0.337 = 0.0088968.
    # The correlation's 95% CI spans zero, so this term is weakly identified.
    etalrbase_apoc3 + etalrbase_tg ~ c(0.030976,
                                       0.0088968, 0.113569)                                         # Table 2: On BASE_APOC3 SD 0.176 (17.7% CV, RSE 28.4%, shrinkage 31.7%); On BASE_TG SD 0.337 (34.6% CV, RSE 20.1%, shrinkage 15.2%); correlation 0.150

    # --------------------------------------------------------------------
    # Residual unexplained variability. "Proportional error models were used
    # to describe the residual unexplained variability for both APOC3 and TG
    # measurements", and the values are "36.2% and 52.7% on the standard
    # deviation scale".
    # --------------------------------------------------------------------
    propSd_apoc3 <- 0.362; label("Proportional residual SD on serum APOC3 (fraction)")              # Table 2: APOC3 proportional error 36.2%, RSE 2.04%, 95% CI 34.8-37.7
    propSd_tg    <- 0.527; label("Proportional residual SD on serum triglycerides (fraction)")      # Table 2: TG proportional error 52.7%, RSE 2.66%, 95% CI 49.9-55.4
  })

  model({
    # ------------------------------------------------------------------
    # 0. Unit conversion. The model time unit is days (dosing is Q3M and the
    #    paper simulates on a day grid: "Q3M dosing was simulated as dosing
    #    every 90 days", "One month is calculated as 30 days"). Table 2
    #    reports ksyn_APOC3 in mg/dL/h and kdeg_TG in /h, so both are scaled
    #    by 24 h/day. kel is already /day and is not scaled.
    # ------------------------------------------------------------------
    hPerDay <- 24  # h/day

    # ------------------------------------------------------------------
    # 1. Individual parameters. Log-normal IIV; covariate effects are power
    #    functions of the centred covariate, except the categorical
    #    triglyceride-lowering-therapy effect which is exp(theta * indicator).
    # ------------------------------------------------------------------
    kel <- exp(lkel + etalkel) * (BMI / 25)^e_bmi_kel

    rbase_apoc3 <- exp(lrbase_apoc3 + etalrbase_apoc3) *
      (APOC3 / 34)^e_apoc3_rbase_apoc3
    ksyn_apoc3 <- exp(lksyn_apoc3) * hPerDay          # mg/dL/day
    imax_ploz  <- exp(limax_ploz)
    ic50_ploz  <- exp(lic50_ploz + etalic50_ploz)

    rbase_tg <- exp(lrbase_tg + etalrbase_tg) *
      (TRIG / 2048)^e_trig_rbase_tg
    kdeg_tg     <- exp(lkdeg_tg) * hPerDay            # 1/day
    imax_apoc3  <- exp(limax_apoc3)
    ki50_apoc3  <- exp(lki50_apoc3 + etalki50_apoc3) *
      (BMI / 25)^e_bmi_ki50_apoc3 *
      exp(e_conmed_tglower_ki50_apoc3 * CONMED_TGLOWER)

    # ------------------------------------------------------------------
    # 2. Turnover reparameterisation. Table 2 estimates the two BASELINES and
    #    fixes ksyn_APOC3 and kdeg_TG, but reports neither kdeg_APOC3 nor
    #    ksyn_TG. Both are therefore forced by the drug-free steady-state
    #    balance of the Supporting-Information ODEs -- this is arithmetic, not
    #    an unreported parameter:
    #
    #      APOC3: at t = 0 (Ce = 0, DrugEffect_ploz = 1) the ODE
    #             0 = ksyn_APOC3 - kdeg_APOC3 * BASE_APOC3
    #             gives kdeg_APOC3 = ksyn_APOC3 / BASE_APOC3.
    #             (1.13 * 24 / 30.9 = 0.878 /day, i.e. t1/2 = 19 h.)
    #
    #      TG:    at t = 0 the APOC3 inhibition term is already active at the
    #             baseline APOC3 level, so
    #             0 = ksyn_TG - kdeg_TG * BASE_TG *
    #                   (1 - Imax_APOC3 * BASE_APOC3 / (IC50_APOC3 + BASE_APOC3))
    #             gives the ksyn_TG below. Omitting the baseline inhibition
    #             factor here would leave tg(0) off its own steady state.
    #
    #    Both are computed per subject from that subject's own baselines and
    #    ki50, so every individual starts exactly at their own baseline.
    # ------------------------------------------------------------------
    kdeg_apoc3 <- ksyn_apoc3 / rbase_apoc3
    ksyn_tg <- kdeg_tg * rbase_tg *
      (1 - imax_apoc3 * rbase_apoc3 / (ki50_apoc3 + rbase_apoc3))

    # ------------------------------------------------------------------
    # 3. Cascading K-PD indirect-response system (Supporting Information,
    #    Supplementary Equations; gamma == 1 as the Hill coefficient was
    #    rejected). `depot_kpd` holds the plozasiran AMOUNT in mg in the
    #    volume-less biophase, so IC50_ploz is an amount in mg, not a
    #    concentration.
    # ------------------------------------------------------------------
    drug_effect_ploz <- 1 - imax_ploz * depot_kpd / (ic50_ploz + depot_kpd)
    apoc3_effect_tg  <- 1 - imax_apoc3 * apoc3 / (ki50_apoc3 + apoc3)

    d/dt(depot_kpd) <- -kel * depot_kpd
    d/dt(apoc3)     <- ksyn_apoc3 * drug_effect_ploz - kdeg_apoc3 * apoc3
    d/dt(tg)        <- ksyn_tg - kdeg_tg * tg * apoc3_effect_tg

    apoc3(0) <- rbase_apoc3
    tg(0)    <- rbase_tg

    # ------------------------------------------------------------------
    # 4. Observations. Both endpoints are the serum concentrations of the
    #    states themselves (mg/dL), each with a proportional residual.
    # ------------------------------------------------------------------
    apoc3 ~ prop(propSd_apoc3)
    tg    ~ prop(propSd_tg)
  })
}
