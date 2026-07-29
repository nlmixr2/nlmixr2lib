Bijleveld_2016_gentamicin <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in term neonates with hypoxic-ischaemic encephalopathy undergoing controlled hypothermia (Bijleveld 2016), with fixed allometric body-weight scaling (exponents 0.75 on CL and Q, 1 on Vc and Vp), an estimated gestational-age power effect on CL, and a categorical post-rewarming (study-day-5, > 96 h PNA) multiplicative increase in CL."
  reference <- paste(
    "Bijleveld YA, de Haan TR, van der Lee HJH, Groenendaal F, Dijk PH,",
    "van Heijst A, de Jonge RCJ, Dijkman KP, van Straaten HLM, Rijken M,",
    "Zonnenberg IA, Cools F, Zecic A, Nuytemans DHGM, van Kaam AH,",
    "Mathot RAA; PharmaCool study group.",
    "Altered gentamicin pharmacokinetics in term neonates undergoing",
    "controlled hypothermia.",
    "Br J Clin Pharmacol. 2016;81(6):1067-1077.",
    "doi:10.1111/bcp.12883.",
    sep = " "
  )
  vignette <- "Bijleveld_2016_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (birth weight in the source cohort)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median 3.4 kg (range 2.09-5.07) in the gentamicin cohort (Table 1).",
        "Allometric scaling uses a 70 kg adult reference: fixed exponent",
        "0.75 on CL and Q, fixed exponent 1 on Vc and Vp."
      ),
      source_name        = "BW"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median 40 weeks (range 36-42) in the cohort (Table 1). Time-fixed",
        "per subject. Enters CL as a power function (GA/40)^e_ga_cl where",
        "40 weeks is the cohort median used as the reference."
      ),
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. The source paper reports PNA in hours and defines",
        "Study Day 5 (SD5) as PNA > 96 h, the period after the 72 h",
        "controlled-hypothermia phase plus 8 h rewarming = > 96 h PNA",
        "(post-rewarming normothermia). Canonical PNA in nlmixr2lib",
        "carries months, so the model converts internally as",
        "PNA_hours = PNA * 24 * 30.4375 and sets the SD5 indicator =",
        "(PNA_hours > 96).  PNA itself is not a retained continuous",
        "covariate on CL; only the binary post-96 h step is."
      ),
      source_name        = "PNA"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (female indicator)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in covariate analysis; not retained (no relationship with CL or V)."
    ),
    BMTEMP = list(
      description = "Core body temperature",
      units       = "degC",
      type        = "continuous",
      notes       = "Screened; not retained on CL (paper Results 'Pharmacokinetic model building')."
    ),
    COOLING = list(
      description = "Cooling on/off indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a categorical covariate; not retained on CL. The",
        "model nonetheless captures the post-rewarming clearance change",
        "via the derived SD5 step (PNA > 96 h)."
      )
    ),
    INOTROPIC = list(
      description = "Inotropic co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained on CL."
    ),
    MOF = list(
      description = "Multi-organ failure indicator (renal or liver dysfunction per Shah et al.)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained on CL."
    ),
    SCR = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained on CL."
    ),
    UREA = list(
      description = "Serum urea",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened; not retained on CL."
    ),
    ASAT = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained on CL."
    ),
    ALAT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; not retained on CL."
    ),
    URINE_OUT = list(
      description = "Daily urine output",
      units       = "mL/day",
      type        = "continuous",
      notes       = "Screened; not retained on CL."
    ),
    PMA = list(
      description = "Postmenstrual age (GA + PNA)",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and produced a larger OFV drop than PNA alone,",
        "but GA was preferred as the most influential variable (paper",
        "Results 'Pharmacokinetic model building'). PMA was therefore",
        "not retained in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 1L,
    n_observations = 612L,
    ga_range       = "36-42 weeks (median 40)",
    pna_range      = "2.3-5.2 days at end of study (median 4.7)",
    weight_range   = "2.09-5.07 kg birth weight (median 3.4)",
    sex_female_pct = 41.3,
    disease_state  = paste(
      "Term newborns (GA > 37 weeks at birth, with three patients of",
      "GA 36 wk admitted under the protocol) meeting criteria for",
      "perinatal asphyxia and hypoxic-ischaemic encephalopathy,",
      "treated with controlled whole-body hypothermia (target core",
      "temperature 33.5 degC) for 72 h beginning within 6 h of birth,",
      "followed by 8 h of rewarming to 36.5 degC and then normothermia."
    ),
    dose_range     = paste(
      "Gentamicin 4 mg/kg once daily per the Dutch Paediatric Formulary,",
      "with TDM-driven dose adjustments. Cohort range 3.5-5.1 mg/kg per",
      "dose; intervals 24-48 h. 70% of patients received 4 mg/kg q24h,",
      "17% 4 mg/kg q36h, 6% 4 mg/kg q48h."
    ),
    co_medications = "Inotropic support in 63.8% of patients; analgesia and antiepileptics per ICU standard of care.",
    regions        = "10 Dutch and 2 Belgian neonatal intensive care units (PharmaCool Study), November 2010 to October 2014.",
    notes          = paste(
      "Baseline demographics from Bijleveld 2016 Table 1. Final two-",
      "compartment model with allometric body-weight scaling (fixed",
      "exponents per West allometry), gestational age as a power",
      "covariate on CL, and a categorical multiplicative increase in CL",
      "after the rewarming phase ends (PNA > 96 h). The bioavailability",
      "of an assumed pre-study 4 mg/kg dose was estimated at 142% with",
      "IIV 75.4% to extrapolate the pre-study time profile and is NOT",
      "carried in this library model (the user supplies observed doses",
      "directly; see vignette Errata). Estimation: NONMEM 7.2, FOCE-I."
    )
  )

  ini({
    # Structural population PK parameters (Bijleveld 2016 Table 2 Final Model
    # column). All four typical values are reported on the linear scale per
    # 70 kg reference; nlmixr2lib log-transforms positive-constrained
    # parameters.
    lcl <- log(1.89);  label("Clearance (CL, L/h per 70 kg)")                     # Bijleveld 2016 Table 2 Final model: 1.89 L/h/70 kg
    lvc <- log(32.5);  label("Central volume of distribution (Vc, L per 70 kg)")   # Bijleveld 2016 Table 2 Final model: 32.5 L/70 kg
    lq  <- log(2.01);  label("Intercompartmental clearance (Q, L/h per 70 kg)")    # Bijleveld 2016 Table 2 Final model: 2.01 L/h/70 kg
    lvp <- log(30.3);  label("Peripheral volume of distribution (Vp, L per 70 kg)") # Bijleveld 2016 Table 2 Final model: 30.3 L/70 kg

    # Allometric body-weight exponents (Bijleveld 2016 Methods 'Structural
    # model' and Results 'Pharmacokinetic model building': 'Parameters
    # were normalized for a body weight of 70 kg using fixed exponents
    # according to the 3/4 rule', citing West et al. The paper also tested
    # estimating the exponents but the fit did not improve).
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Bijleveld 2016 Methods/Results: fixed 0.75
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)") # Bijleveld 2016 Methods/Results: fixed 1

    # Gestational-age power effect on CL (Bijleveld 2016 Table 2 Final
    # Model footer: theta_CLGA = 3.00, CV 16%). Enters as
    # (GA / 40 weeks)^e_ga_cl where 40 weeks is the cohort median.
    e_ga_cl <- 3.00;  label("Power exponent of (GA/40 wk) on CL (unitless)") # Bijleveld 2016 Table 2 Final model: theta_CLGA = 3.00

    # Study-day-5 (post-rewarming, > 96 h PNA) multiplicative effect on CL
    # (Bijleveld 2016 Table 2 Final Model footer: theta_SD5 = 1.29, CV 12%).
    # When the SD5 indicator is 1 (PNA > 96 h), CL is multiplied by 1.29
    # (i.e. a 29% increase). When the indicator is 0, the factor reduces
    # to 1 via e_sd5_cl^0 = 1.
    e_sd5_cl <- 1.29; label("Multiplicative CL factor when PNA > 96 h (unitless)") # Bijleveld 2016 Table 2 Final model: theta_SD5 = 1.29

    # Inter-individual variability (Bijleveld 2016 Table 2 Final Model
    # column). Reported as %CV; for log-normal eta the variance on the
    # log scale is omega^2 = log(1 + CV^2).
    #   IIV CL : 26.6% CV  -> log(1 + 0.266^2) = 0.06837
    #   IIV Vc : 40.8% CV  -> log(1 + 0.408^2) = 0.15399
    #   IIV Vp : 53.3% CV  -> log(1 + 0.533^2) = 0.25022
    # CL and Vc are correlated (r = 0.81, paper Results 'Pharmacokinetic
    # model building': 'Variability between CL and Vc was correlated
    # (r = 0.81)'); the lower-triangular covariance is
    #   cov(CL, Vc) = r * omega_CL * omega_Vc
    #              = 0.81 * sqrt(0.06837) * sqrt(0.15399) = 0.08316.
    etalcl + etalvc ~ c(0.06837,
                        0.08316, 0.15399)  # Bijleveld 2016 Table 2 Final: IIV CL 26.6%, IIV Vc 40.8%, corr 0.81
    etalvp ~ 0.25022                       # Bijleveld 2016 Table 2 Final: IIV Vp 53.3%

    # Residual error. The paper log-transformed the observations and
    # estimated an additive residual on the log scale with sigma = 0.15
    # (Table 2 Final Model). The NONMEM LTBS pattern
    #   log(Cobs) = log(Cpred) + eps,  eps ~ N(0, expSd^2)
    # maps directly to nlmixr2's lnorm() residual with expSd = sigma.
    # The paper additionally reports a 50.2% IIV on the additive residual
    # error (subject-level scaling of the residual magnitude); that
    # second-level structure is not encoded here (see vignette Errata).
    expSd <- 0.15;  label("Lognormal residual SD (additive on log-transformed concentration)") # Bijleveld 2016 Table 2 Final model: additive sigma on log-data = 0.15
  })

  model({
    # Derive Study-Day-5 indicator from canonical PNA (months). The
    # source paper anchors SD5 at PNA > 96 h, i.e. after the 72 h
    # controlled-hypothermia phase plus the 8 h rewarming window plus
    # a 16 h post-rewarming buffer (4 days = 96 h PNA). The threshold
    # 96 h / (24 h/d * 30.4375 d/month) = 0.13143... months.
    pna_hours    <- PNA * 24 * 30.4375
    sd5_ind      <- (pna_hours > 96) * 1.0

    # Individual PK parameters (Bijleveld 2016 Table 2 Final Model footer
    # TVCL/TVVc/TVQ/TVVp equations).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q   * (e_sd5_cl)^sd5_ind * (GA / 40)^e_ga_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK. Doses are administered as short intravenous
    # infusions (typically 30 min) in the source paper; the library model
    # does not hard-code an infusion duration so users supply rate or
    # dur per dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
