Chen_2024_combinedOralContraceptives_btb_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of BREAKTHROUGH BLEEDING (BTB) - ",
    "unscheduled bleeding during active hormone intake - in women taking ",
    "combined oral contraceptives (COCs), fit to 228 aggregate ",
    "observations from 33 treatment arms of 25 published trials of four ",
    "progestins (desogestrel DSG, drospirenone DRSP, gestodene GSD, ",
    "levonorgestrel LNG) combined with ethinyl estradiol (EE) at 15-35 ",
    "ug (Chen 2024 Table 1, Table S2). The response is the FRACTION of ",
    "women in a treatment arm experiencing BTB, described by a ",
    "bi-exponential decline from COC initiation (Chen 2024 Eq. 3): a ",
    "rapid initial phase (intercept A, slope alpha) plus a slow second ",
    "phase (intercept B, slope beta). EE dose enters as a power effect ",
    "on B (normalized to 30 ug) and the molar progestin dose as a power ",
    "effect on alpha (normalized to 0.48 umol); progestin TYPE enters ",
    "only through that molar dose. Variability is INTER-STUDY (ISV) on ",
    "all four macro-constants, not between-subject. Simulation scope is ",
    "the AGGREGATE treatment-arm BTB fraction for a COC regimen, NOT an ",
    "individual woman's bleeding outcome and NOT a concentration; the ",
    "model contains no PK, no ODE state and no dose events. NOTE THE ",
    "TIME ORIGIN: t = 0 is the completion of ONE full COC cycle, so t ",
    "in months equals (months since COC initiation) - 1 (Chen 2024 ",
    "Methods 'Model-based meta-analysis'; see the vignette)."
  )

  reference <- paste(
    "Chen H, Chun D, Lingineni K, Guzy S, Cristofoletti R, Hoechel J,",
    "Jiao T, Cicali B, Vozmediano V, Schmidt S. Development of",
    "breakthrough bleeding model of combined-oral contraceptives",
    "utilizing model-based meta-analysis. CPT Pharmacometrics Syst",
    "Pharmacol. 2024;13(11):2016-2025. doi:10.1002/psp4.13261.",
    sep = " "
  )

  vignette <- "Chen_2024_combinedOralContraceptives_btb"

  # The four random effects are INTER-STUDY variance (ISV) terms on the
  # bi-exponential macro-constants (Chen 2024 Eqs. 4-7, eta_i ~ N(0,
  # omega^2) with theta_i = theta_pop * exp(eta_i)), i.e. how much a
  # study's underlying BTB trajectory differs from the typical study.
  # They are NOT between-SUBJECT effects on a structural PK parameter,
  # so they carry the explicit `eta_isv_` prefix rather than the popPK
  # `eta<transformed-param>` convention. Declared paper-specific per the
  # Goteti_2024_SLE_mbma / Hanan_2026_peginterferon_alfa_eot_mbma
  # precedent.
  paper_specific_etas <- c(
    "eta_isv_a", "eta_isv_alpha", "eta_isv_b", "eta_isv_beta"
  )

  units <- list(
    time          = "month",
    dosing        = paste(
      "(no rxode2 dose events; the COC regimen enters only as the two",
      "covariate columns DOSE_EE_UG, the daily ethinyl estradiol dose",
      "in ug, and DOSE_PROGESTIN_UMOL, the daily progestin dose in",
      "umol. Every arm in the analysis set used a conventional 21/7",
      "monophasic regimen - multiphasic dosing and deviations from the",
      "21/7 schedule were exclusion criteria - so the regimen is fully",
      "described by those two daily dose levels.)"
    ),
    concentration = paste(
      "(dimensionless fraction, 0-1: the proportion of women in a",
      "treatment arm reporting breakthrough bleeding. Aggregate",
      "arm-level summary, not an individual-subject concentration.)"
    )
  )

  covariateData <- list(
    DOSE_EE_UG = list(
      description        = paste(
        "Daily ethinyl estradiol (EE) dose of the combined oral",
        "contraceptive, in ug. Constant within a treatment arm.",
        "Enters as a power effect on the second-phase intercept B,",
        "normalized to a 30 ug reference."
      ),
      units              = "ug",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA treatment-arm-level covariate. Chen 2024 Eq. 6:",
        "log(B_i) = log(theta_B,pop) + log(DOSE_EE_UG / 30) *",
        "theta_EEdose,B + eta_B,i, i.e. B_i = theta_B,pop *",
        "(DOSE_EE_UG / 30)^theta_EEdose,B * exp(eta_B,i). The",
        "coefficient is NEGATIVE (-2.45, Table 2 row 'EE dose on B'),",
        "so LOWER EE doses raise the slow-phase intercept and lengthen",
        "the time for BTB to subside - the paper's central finding.",
        "Relative to the 30 ug reference, B is (20/30)^-2.45 = 2.70x",
        "higher at 20 ug and (15/30)^-2.45 = 5.46x higher at 15 ug.",
        "",
        "The 30 ug reference is the highest EE dose among the",
        "FDA-approved combinations the paper simulates, and the level",
        "at which the paper reports BTB returning to baseline within",
        "3-4 months.",
        "",
        "Observed range in the analysis set is 15, 20, 30 and 35 ug",
        "(Chen 2024 Table 1 and Table S2); the 35 ug level comes from a",
        "single GSD 75/35 arm. Extrapolating the power term outside",
        "15-35 ug is not supported by the data, and it diverges as the",
        "dose approaches 0 (a progestin-only pill is NOT this model's",
        "domain - the paper excluded progestin-only-pill data from",
        "model development).",
        "",
        "This is a dose-as-covariate column rather than an rxode2",
        "`amt` / `EVID = 1` dosing event because the model is purely",
        "algebraic with no PK compartment, following the DOSE_AGT_UG /",
        "DOSE_ISOPROTERENOL_UG precedent."
      ),
      source_name        = "EE dose (mcg)"
    ),
    DOSE_PROGESTIN_UMOL = list(
      description        = paste(
        "Daily progestin dose of the combined oral contraceptive",
        "expressed on a MOLAR basis, in umol. Constant within a",
        "treatment arm. Enters as a power effect on the initial-phase",
        "slope alpha, normalized to a 0.48 umol reference."
      ),
      units              = "umol",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA treatment-arm-level covariate. Chen 2024 Eq. 5:",
        "log(alpha_i) = log(theta_alpha,pop) +",
        "log(DOSE_PROGESTIN_UMOL / 0.48) * theta_ProgestinAMT,alpha +",
        "eta_alpha,i, i.e. alpha_i = theta_alpha,pop *",
        "(DOSE_PROGESTIN_UMOL / 0.48)^theta_ProgestinAMT,alpha *",
        "exp(eta_alpha,i).",
        "",
        "MOLAR, not mass: Chen 2024 Methods states 'Given that",
        "progestins differ in molecular weight, molar doses were used",
        "in the analysis to mitigate collinearity issues.' Populate the",
        "column by dividing the daily progestin dose in ug by the",
        "progestin's molecular weight in g/mol:",
        "desogestrel 310.48, drospirenone 366.50, gestodene 310.43,",
        "levonorgestrel 312.45. The molar doses of the arms the paper",
        "simulates are therefore DSG 150 ug = 0.483 umol,",
        "DRSP 3000 ug = 8.186 umol, GSD 60 ug = 0.193 umol,",
        "GSD 75 ug = 0.242 umol, LNG 100 ug = 0.320 umol and",
        "LNG 150 ug = 0.480 umol.",
        "",
        "The 0.48 umol reference is not named in the text, but it is",
        "exactly the levonorgestrel 150 ug molar dose (150 / 312.45 =",
        "0.4800) - and, to two decimals, the desogestrel 150 ug dose",
        "(0.4831). That coincidence is what pins the molecular-weight",
        "convention above: no other MW set reproduces 0.48 for the",
        "paper's reference regimen.",
        "",
        "Progestin TYPE is not a separate model covariate. Chen 2024",
        "retained the molar dose as the only progestin descriptor, so",
        "the four progestins are distinguished purely by the umol value",
        "this column carries; see covariatesDataExcluded for the",
        "type-level screen. The paper kept this term despite it losing",
        "statistical significance once EE dose entered on B, 'due the",
        "pharmacological interplay between progestin and EE' (Results).",
        "",
        "This is a dose-as-covariate column rather than an rxode2",
        "`amt` / `EVID = 1` dosing event because the model is purely",
        "algebraic with no PK compartment, following the DOSE_AGT_UG /",
        "DOSE_ISOPROTERENOL_UG precedent."
      ),
      source_name        = "Progestin dose by MW"
    )
  )

  # Covariates Chen 2024 SCREENED against the structural parameters but
  # did NOT retain in the final model. None has a published point
  # estimate, so there is nothing to encode; they are recorded here to
  # preserve the provenance of the covariate screen without raising a
  # "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    BMI = list(
      description        = "Treatment-arm mean body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY excluded covariate that reached statistical",
        "significance: Chen 2024 Results states 'Mean BMI showed a",
        "significant effect on A. However, BMI was not included in the",
        "final covariate model due to the limited BMI range, with the",
        "majority of values laying in a relatively narrow range",
        "(21.4-23.5 kg/m2), with the exception of DRSP/EE 3000/30 ug",
        "group, where subjects had a significantly higher BMI of 27.6",
        "kg/m2.' No coefficient is reported. The Discussion names a",
        "broader BMI range as the key requirement for future studies.",
        "",
        "Arm-level means span 21.4-27.6 kg/m2 (Table S2)."
      ),
      source_name        = "Mean BMI (kg/m2)"
    ),
    AGE = list(
      description        = "Treatment-arm mean age.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as one of the study-level average demographics",
        "(Chen 2024 Methods 'Model-based meta-analysis') and not",
        "retained; no point estimate reported. Arm-level means span",
        "16.4-40.2 years (Table 1, Table S2)."
      ),
      source_name        = "Mean age (year)"
    ),
    WT = list(
      description        = "Treatment-arm mean body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained; no point estimate reported.",
        "Arm-level means span 46.0-68.3 kg (Table S2)."
      ),
      source_name        = "Mean BW (kg)"
    ),
    HT = list(
      description        = "Treatment-arm mean height.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained; no point estimate reported.",
        "Arm-level means span 153.1-168.0 cm (Table S2). Table 1",
        "labels this row 'Mean HT (m)' but prints centimetre values",
        "(e.g. 166.20); the unit label is a typo - see the vignette",
        "Errata."
      ),
      source_name        = "Mean HT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112188L,
    n_studies      = 25L,
    n_arms         = 33L,
    n_observations = 228L,
    age_range      = paste(
      "treatment-arm mean age 16.4-40.2 years (Chen 2024 Table 1",
      "total row and Table S2); the 16.4-year arm is an adolescent",
      "cohort. Individual ages are not available - this is aggregate",
      "published data."
    ),
    weight_range   = "treatment-arm mean body weight 46.0-68.3 kg (Table S2)",
    sex_female_pct = 100,
    disease_state  = paste(
      "healthy women using a combined oral contraceptive for pregnancy",
      "prevention. Alternative COC indications were an exclusion",
      "criterion."
    ),
    dose_range     = paste(
      "conventional 21/7 monophasic COC regimens: ethinyl estradiol",
      "15, 20, 30 or 35 ug/day combined with desogestrel 150 ug,",
      "drospirenone 3000 ug, gestodene 60 or 75 ug, or levonorgestrel",
      "100 or 150 ug per day (Chen 2024 Table 1, Table S2)."
    ),
    notes          = paste(
      "Aggregate published data, not individual patient data. The",
      "observation unit is the TREATMENT ARM within a study: 228 BTB",
      "observations from 33 arms of 25 trials (Chen 2024 Results",
      "'Data'). n_subjects is the sum of the per-arm subject counts",
      "printed in Table S2; it is not stated anywhere in the paper and",
      "is dominated by a single very large arm (Brill et al., n =",
      "95,906 of the 112,188 total), so it should be read as an",
      "enrolment total rather than an effective sample size. Arm sizes",
      "span 30 to 95,906 women.",
      "",
      "Analysis-set construction (Methods 'Data'): a systematic search",
      "of PubMed, Cochrane and EMBASE returned 1147 papers, of which 46",
      "carried BTB data; the final set kept only the four progestins",
      "with sufficient longitudinal data (DSG, DRSP, GSD, LNG).",
      "Progestin-only-pill data were EXCLUDED from model development",
      "(too few studies), as were alternative COC indications,",
      "non-observational studies, arms without active hormone,",
      "multiphasic dosing, deviations from the 21/7 regimen and review",
      "articles. Bleeding within the first 7 days of cycle 1 or days",
      "1-4 of later cycles, and bleeding during hormone-free intervals,",
      "were excluded from the BTB definition.",
      "",
      "Missing arm-level covariates were imputed by 'imputation by",
      "individual parameter' (IIP), a log-linear regression of the",
      "covariate on the base-model individual parameter, chosen over",
      "median and multiple imputation (Methods 'Data'; Table S2 carries",
      "the imputed values). Cycles with fewer than 100 subjects in",
      "which the observed BTB probability was < 1% were set missing and",
      "imputed at HALF the 1.68% untreated baseline rate, i.e. 0.84%.",
      "",
      "The 1.68% untreated baseline unscheduled-bleeding rate itself",
      "comes from 44,420 pre-treatment subjects in seven internal Bayer",
      "studies - an external anchor, not part of the 25-study analysis",
      "set."
    )
  )

  ini({
    # ================================================================
    # Chen 2024 Eq. 3, the final bi-exponential ("2-compartment") BTB
    # model:
    #
    #   BTB(x) = A * exp(-alpha * x) + B * exp(-beta * x)
    #
    # where x is TIME IN MONTHS SINCE THE FIRST OBSERVATION, and the
    # first observation is the completion of one full COC cycle
    # (Methods: "the time of the first observation was set to 0,
    # corresponding to the completion of one full cycle of COC
    # treatment"). So x = (months since COC initiation) - 1.
    #
    # A / alpha describe the rapid initial decline, B / beta the slow
    # second phase. A 1-compartment form gave biased estimates at later
    # times; a 3-compartment form did not improve the fit and was
    # unstable (Results).
    #
    # All point estimates below are Chen 2024 Table 2, "Estimate"
    # column. Eqs. 4-7 log-transform every macro-constant, which is why
    # each is stored here on the log scale.
    # ================================================================

    la_btb <- log(0.0383)
    label("Initial-phase BTB intercept A (fraction of women, 0-1)")
    # Table 2, Fixed effect parameters, "A" = 0.0383 (RSE 21.7%;
    # bootstrap median 0.037, 95% CI 0.026-0.050). Eq. 4:
    # log(A_i) = log(theta_A,pop) + eta_A,i -- no covariate retained.

    lalpha_btb <- log(0.922)
    label("Initial-phase BTB decline rate alpha at the 0.48 umol reference progestin dose (1/month)")
    # Table 2, Fixed effect parameters, "alpha" = 0.922 (RSE 35.4%;
    # bootstrap median 1.050, 95% CI 0.032-1.813). Eq. 5 carries the
    # molar-progestin-dose power term; see e_dose_progestin_umol_alpha.

    lb_btb <- log(0.0134)
    label("Second-phase BTB intercept B at the 30 ug reference EE dose (fraction of women, 0-1)")
    # Table 2, Fixed effect parameters, "B" = 0.0134 (RSE 34.2%;
    # bootstrap median 0.014, 95% CI 0.006-0.021). Eq. 6 carries the
    # EE-dose power term; see e_dose_ee_ug_b.

    lbeta_btb <- log(0.0524)
    label("Second-phase BTB decline rate beta (1/month)")
    # Table 2, Fixed effect parameters, "beta" = 0.0524 (RSE 32.5%;
    # bootstrap median 0.056, 95% CI 0.028-0.077). Eq. 7:
    # log(beta_i) = log(theta_beta,pop) + eta_beta,i -- no covariate.

    # ---------------- Covariate effects (Eqs. 5-6) ------------------
    # Both are POWER exponents on a normalized covariate, NOT
    # multiplicative factors: Eqs. 5 and 6 multiply theta by the LOG of
    # the normalized covariate, which exponentiates to
    # (cov / ref)^theta. The alternative reading -- theta inside the
    # log, i.e. B = theta_B,pop * (EEdose/30) * theta_EEdose,B -- is
    # arithmetically impossible here: with theta_EEdose,B = -2.45 it
    # would make B negative at every EE dose, and it would put the
    # typical value of alpha at 0.922 * 0.576 = 0.531 rather than the
    # 0.922 Table 2 reports at the reference dose.

    e_dose_ee_ug_b <- -2.45
    label("Power exponent of (DOSE_EE_UG / 30) on the second-phase intercept B (unitless)")
    # Table 2, Fixed effect parameters, "EE dose on B" = -2.45 (RSE
    # 32.6%; bootstrap median -2.422, 95% CI -3.873 to -1.028).

    e_dose_progestin_umol_alpha <- 0.576
    label("Power exponent of (DOSE_PROGESTIN_UMOL / 0.48) on the initial-phase decline rate alpha (unitless)")
    # Table 2, Fixed effect parameters, "Progestin dose by MW on alpha"
    # = 0.576 (RSE 22.4%; bootstrap median 0.573, 95% CI -0.231 to
    # 1.383 -- the bootstrap CI spans zero, consistent with the
    # Results statement that this term lost significance once EE dose
    # was added on B but was retained for pharmacological reasons).

    # ================================================================
    # Inter-study variance (ISV). Eq. 2: theta_i = theta_pop *
    # exp(eta_i), so these are VARIANCES on the natural-log scale.
    # Table 2 prints each alongside a %CV that confirms the scale:
    # sqrt(exp(0.824) - 1) = 113%, sqrt(exp(0.959) - 1) = 127%,
    # sqrt(exp(0.223) - 1) = 49.9% and sqrt(exp(0.138) - 1) = 38.5%,
    # matching the four printed CVs exactly. They are BETWEEN-STUDY,
    # not between-subject: they describe how much a new STUDY-ARM's
    # BTB trajectory differs from the typical one.
    #
    # Chen 2024 Methods states "Covariances among the inter-study
    # variance (ISV) were also estimated to account for potential
    # correlations between parameters", but Table 2 reports only the
    # four diagonal terms and no covariance or correlation matrix
    # appears anywhere in the paper or its supplement. The off-diagonal
    # elements are therefore encoded as zero; the vignette shows that a
    # diagonal-only omega still reproduces the paper's own Table S3 and
    # Table S4 simulation summaries. See the vignette Errata.
    # ================================================================

    eta_isv_a ~ 0.824
    # Table 2, Random effect parameters, "ISV-A" = 0.824 (113% CV; RSE
    # 40.3%; bootstrap median 0.859, 95% CI 0.255-1.392). The Results
    # flag the high RSE on this term and on ISV-alpha as reflecting
    # sparse BTB data immediately after COC initiation.

    eta_isv_alpha ~ 0.959
    # Table 2, Random effect parameters, "ISV-alpha" = 0.959 (127% CV;
    # RSE 42.3%; bootstrap median 0.987, 95% CI -0.095 to 2.012).

    eta_isv_b ~ 0.223
    # Table 2, Random effect parameters, "ISV-B" = 0.223 (49.9% CV; RSE
    # 67.6%; bootstrap median 0.205, 95% CI -0.124 to 0.569).

    eta_isv_beta ~ 0.138
    # Table 2, Random effect parameters, "ISV-beta" = 0.138 (38.5% CV;
    # RSE 84.3%; bootstrap median 0.126, 95% CI -0.009 to 0.284).

    # ================================================================
    # Residual error. Chen 2024 Eq. 8 is additive on the BTB fraction,
    # Y_ijk = F_ijk + eps_ijk, and Eq. 9 gives the ARM-SIZE-WEIGHTED
    # distribution
    #
    #   eps_ijk ~ N(0, sigma^2 / N_ik)
    #
    # where N_ik is the number of individuals in treatment arm k of
    # study i. Table 2 reports sigma^2 = 0.0410 with the label "(20.2%
    # CV)", and sqrt(0.0410) = 0.2025, so 0.0410 is the VARIANCE and
    # 0.2025 the unit-arm SD.
    #
    # nlmixr2's add() takes a single population parameter and cannot
    # read a per-record arm size, so addSd below is the UNIT-ARM
    # residual SD sigma (i.e. the SD an arm of N = 1 would carry). The
    # residual SD actually applying to an arm of size N is
    # addSd / sqrt(N) -- for the median arm in this analysis (308
    # women, Table 1) that is 0.2025 / sqrt(308) = 0.0115, i.e. about
    # 1.2 percentage points of BTB. A user simulating a specific arm
    # must apply that scaling themselves; the vignette shows how. Same
    # device as Hanan_2026_peginterferon_alfa_eot_mbma, which likewise
    # carries an arm-size variance function nlmixr2 cannot express
    # natively.
    # ================================================================

    addSd <- 0.2025
    label("Unit-arm additive residual SD on the BTB fraction; an arm of N women carries addSd / sqrt(N)")
    # Table 2, Residual unexplained variance, "Unweighted additive
    # error" = 0.0410 (20.2% CV; RSE 35.2%; bootstrap median 0.035, 95%
    # CI 0.024-0.058). sqrt(0.0410) = 0.20248.
  })

  model({
    # ---------------- Untreated reference ---------------------------
    # Rate of unscheduled bleeding in the ABSENCE of any hormonal
    # contraceptive, from 44,420 pre-treatment subjects across seven
    # internal Bayer studies (Chen 2024 Methods 'Data'). It is not a
    # fitted parameter - it is the external anchor the paper uses both
    # to impute low-count cycles (at half this value) and to define
    # "time to return to baseline" in Table S4. Emitted so that the
    # return-to-baseline calculation is self-contained.
    btb_baseline_untreated <- 0.0168

    # ---------------- Study-level macro-constants (Eqs. 4-7) --------
    # Each is theta_pop * (covariate power term) * exp(eta_isv), the
    # exponentiated form of the paper's log-scale equations.

    # Eq. 4 - no covariate on A.
    a_btb <- exp(la_btb + eta_isv_a)

    # Eq. 5 - molar progestin dose as a power effect on alpha,
    # normalized to the 0.48 umol reference.
    alpha_btb <- exp(lalpha_btb + eta_isv_alpha) *
      (DOSE_PROGESTIN_UMOL / 0.48)^e_dose_progestin_umol_alpha

    # Eq. 6 - EE dose as a power effect on B, normalized to the 30 ug
    # reference. The exponent is negative, so lower EE raises B.
    b_btb <- exp(lb_btb + eta_isv_b) *
      (DOSE_EE_UG / 30)^e_dose_ee_ug_b

    # Eq. 7 - no covariate on beta.
    beta_btb <- exp(lbeta_btb + eta_isv_beta)

    # ---------------- Eq. 3, the BTB time course --------------------
    # `time` is months since the FIRST OBSERVATION, which is the
    # completion of one full COC cycle. To read a published
    # "month M of COC treatment" value, evaluate at time = M - 1.
    btb <- a_btb * exp(-alpha_btb * time) +
      b_btb * exp(-beta_btb * time)

    # ---------------- Observation (Eqs. 8-9) ------------------------
    # Cc is the aggregate treatment-arm BTB FRACTION (0-1), not a
    # concentration. addSd is the unit-arm residual SD; an arm of N
    # women carries addSd / sqrt(N) (see the ini() note).
    Cc <- btb
    Cc ~ add(addSd)
  })
}
