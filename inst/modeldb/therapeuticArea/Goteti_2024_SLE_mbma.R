Goteti_2024_SLE_mbma <- function() {
  description <- paste0(
    "MBMA. Bayesian latent-variable disease-trajectory model-based ",
    "meta-analysis (DTM-MBMA) of composite systemic lupus erythematosus ",
    "(SLE) endpoints across 25 studies, 81 study arms, and 16 SLE ",
    "treatments (Goteti 2024 Table 1). A uni-dimensional latent SLE ",
    "disease activity theta(t) is modelled as a mono-exponential ",
    "approach from a baseline mu_ijk to a long-term level mu_ijk + ",
    "delta_ijk with time-constant log(2)/lam (Goteti 2024 Eq 1). The ",
    "paper reparametrizes mu and delta in terms of the more-observable ",
    "theta_wk4 (latent SLE at week 4) and delta_wk24 (increment from ",
    "week 4 to week 24), separately for the SRI (SRI-4/5/6) and BICLA ",
    "endpoints. Mean baseline SLEDAI (centered at 10.5) modifies both ",
    "mu (via beta_sledai_mu) and delta (via beta_sledai_delta). ",
    "Continuous-dose Emax dose-response was fit for six drugs whose ",
    "regimen had three or more doses at the same route and frequency: ",
    "anifrolumab (i.v. q4w), belimumab (i.v. q2wx3+q4w), CC-220 ",
    "(oral qd), epratuzumab (i.v. q2w), lulizumab pegol (s.c. q2w), ",
    "and sifalimumab (i.v. q2w). All 24 other drug-dose-regimen ",
    "combinations from Table 4 (ALX-0061, atacicept, baricitinib, the ",
    "s.c. q.w. belimumab arm, the i.v. q.w. epratuzumab arm, evobrutinib, ",
    "fenebrutinib, IPP-201101, the s.c. q.w. lulizumab pegol arm, lupuzor, ",
    "PF-04236921, tabalumab, ustekinumab) were modelled as discrete ",
    "per-regimen additive effects on the latent-scale delta multiplier ",
    "(1 + xi_ij); the operator passes the reference discrete-effect ",
    "value through the CONMED_XI_DISCRETE covariate (documented lookup ",
    "in the vignette). SRI-4/5/6 responder probabilities are modelled ",
    "as a generalized proportional-odds transform of the latent SLE ",
    "activity (thresholds c_SRI5 = 0.586, c_SRI6 = 0.693), and BICLA ",
    "responder probability is the logistic transform of the latent ",
    "SLE activity directly. Between-study variability is encoded as ",
    "study-level etas on theta_wk4 and delta_wk24 for both endpoints ",
    "(no subject-level IIV); the simulation scope is aggregate ",
    "study-arm-mean responder probability, not individual patient ",
    "trajectories. All Table 3 point estimates are the posterior mean ",
    "under the reference Bayesian fit (model 106-bm-1 in Table S1 of ",
    "the supplement)."
  )

  reference <- paste(
    "Goteti K, Garcia R, Gillespie WR, French J, Klopp-Schulze L, Li Y,",
    "Vazquez Mateo C, Roy S, Guenther O, Benincosa L, Venkatakrishnan K.",
    "Model-based meta-analysis using latent variable modeling to set",
    "benchmarks for new treatments of systemic lupus erythematosus.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(2):281-295.",
    "doi:10.1002/psp4.13083.",
    sep = " "
  )

  vignette <- "Goteti_2024_SLE_mbma"

  # Study-level etas (between-study random effects on the latent-scale
  # theta_wk4_k and delta_wk24_k) do not follow the popPK
  # eta<transformed-param-name> convention because the paper's
  # random-effect names (eta_4wk,i and eta_delta,i) are shared across
  # endpoints and describe study variability rather than subject
  # variability. Declared as paper-specific per Wojciechowski_2015 /
  # Lee_2011 precedent.
  paper_specific_etas <- c(
    "eta_wk4_SRI", "eta_delta_SRI",
    "eta_wk4_BICLA", "eta_delta_BICLA"
  )

  units <- list(
    time          = "week (visit week relative to baseline)",
    dosing        = paste(
      "(no rxode2 dose events; drug-arm dose is a covariate.",
      "Continuous Emax drugs: DOSE_ANIFROLUMAB_MG (mg IV q4w),",
      "DOSE_BELIMUMAB_MGKG (mg/kg IV q2wx3+q4w),",
      "DOSE_CC220_MG (mg oral qd),",
      "DOSE_EPRATUZUMAB_IVQ2W_MG (mg IV q2w),",
      "DOSE_LULIZUMAB_SCQ2W_MG (mg SC q2w),",
      "DOSE_SIFALIMUMAB_MG (mg IV q2w).",
      "Discrete-effect regimens: single CONMED_XI_DISCRETE latent-scale",
      "effect passed in per lookup table)"
    ),
    concentration = paste(
      "(dimensionless responder probability, 0-1; two primary outputs -",
      "p_SRI4 and p_BICLA; two derived outputs - p_SRI5 and p_SRI6 via",
      "the generalized proportional-odds transform. Model is aggregate",
      "study-arm-mean; not for individual-subject simulation.)"
    )
  )

  covariateData <- list(
    SCORE_SLEDAI = list(
      description        = paste(
        "Mean baseline Systemic Lupus Erythematosus Disease Activity",
        "Index (SLEDAI) score at the study-arm level. Enters the DTM as",
        "a centered (SLEDAI - 10.5) covariate on both the latent baseline",
        "mu_ijk (via beta_sledai_mu) and the latent long-term change",
        "delta_ijk (via beta_sledai_delta) per Goteti 2024 Equation 8."
      ),
      units              = "(score, 0-105 possible range; typical trial-arm means 8.4-14)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate: the aggregate mean baseline",
        "SLEDAI across all patients in the study arm, not an individual",
        "score. The covariate register generally holds individual-level",
        "definitions; this per-arm aggregate use is documented here.",
        "Reference (centering) value 10.5 is the across-arm average per",
        "Goteti 2024 Table 2 (mean 10.5, range 8.8-14.0 across 81 arms).",
        "The paper reports beta_sledai_mu = 0.208 (per unit SLEDAI) and",
        "beta_sledai_delta = 0.255 (per unit SLEDAI). Higher baseline",
        "SLEDAI increases both baseline latent disease activity and the",
        "long-term latent change, so the model predicts a higher placebo",
        "response rate for higher-SLEDAI cohorts. Source paper column:",
        "study-arm mean SLEDAI."
      ),
      source_name        = "SLEDAI (mean)"
    ),
    DOSE_ANIFROLUMAB_MG = list(
      description        = paste(
        "Per-arm anifrolumab dose amount (mg IV q4w); 0 if the arm did",
        "not receive anifrolumab."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. The paper fit a continuous",
        "Emax dose-response for the anifrolumab i.v. q4w regimen across",
        "three doses (150, 300, 1000 mg; Goteti 2024 Table 1). Enters",
        "the drug-effect xi_ij via (Dose * Emax) / (Dose + ED50) with",
        "Emax_ANIF = 0.257 and ED50_ANIF = 259 mg (Table 4). Only used",
        "for the anifrolumab i.v. q4w regimen; other routes / frequencies",
        "are not covered by this continuous dose-response."
      ),
      source_name        = "Anifrolumab dose (Goteti 2024 Table 1 / Table 4)"
    ),
    DOSE_BELIMUMAB_MGKG = list(
      description        = paste(
        "Per-arm belimumab dose amount (mg/kg IV q2w x 3 + q4w); 0 if",
        "the arm did not receive belimumab by that regimen. Do NOT use",
        "for the 200 mg s.c. q.w. belimumab arm (which uses a discrete",
        "effect via CONMED_XI_DISCRETE)."
      ),
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. The paper fit a continuous",
        "Emax dose-response for the belimumab i.v. q2wx3+q4w regimen",
        "across three doses (1, 4, 10 mg/kg; Goteti 2024 Table 1). Enters",
        "the drug-effect xi_ij via (Dose * Emax) / (Dose + ED50) with",
        "Emax_BEL = 0.165 and ED50_BEL = 0.787 mg/kg (Table 4). The",
        "200 mg s.c. q.w. belimumab regimen has its own discrete-effect",
        "value 0.16 (Table 4 first belimumab row); pass that as",
        "CONMED_XI_DISCRETE."
      ),
      source_name        = "Belimumab dose (Goteti 2024 Table 1 / Table 4)"
    ),
    DOSE_CC220_MG = list(
      description        = paste(
        "Per-arm CC-220 (iberdomide) dose amount (mg oral qd); 0 if the",
        "arm did not receive CC-220."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. Continuous Emax dose-response",
        "for the CC-220 oral qd regimen across three doses (0.15, 0.30,",
        "0.45 mg; Goteti 2024 Table 1). Enters xi_ij via (Dose * Emax)",
        "/ (Dose + ED50) with Emax_CC220 = 0.415 and ED50_CC220 = 0.415",
        "mg (Table 4). Note the paper's ED50 credible interval crosses",
        "zero on the Emax point estimate, reflecting high posterior",
        "uncertainty for this small dose-response cohort."
      ),
      source_name        = "CC-220 dose (Goteti 2024 Table 1 / Table 4)"
    ),
    DOSE_EPRATUZUMAB_IVQ2W_MG = list(
      description        = paste(
        "Per-arm epratuzumab dose amount for the i.v. q2w regimen (mg);",
        "0 if the arm did not receive epratuzumab i.v. q2w. Do NOT use",
        "for the 600 mg i.v. q.w. epratuzumab arm (which uses a discrete",
        "effect via CONMED_XI_DISCRETE)."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. Continuous Emax dose-response",
        "for the epratuzumab i.v. q2w regimen across three doses (100,",
        "400, 1200, 1800 mg; Goteti 2024 Table 1). Enters xi_ij via",
        "(Dose * Emax) / (Dose + ED50) with Emax_EPRA = 0.397 and",
        "ED50_EPRA = 827 mg (Table 4). The 600 mg i.v. q.w. epratuzumab",
        "regimen has its own discrete-effect value 0.324 (Table 4);",
        "pass that as CONMED_XI_DISCRETE."
      ),
      source_name        = "Epratuzumab dose i.v. q2w (Goteti 2024 Table 1 / Table 4)"
    ),
    DOSE_LULIZUMAB_SCQ2W_MG = list(
      description        = paste(
        "Per-arm lulizumab pegol dose amount for the s.c. q2w regimen",
        "(mg); 0 if the arm did not receive lulizumab pegol s.c. q2w.",
        "Do NOT use for the 12.5 mg s.c. q.w. lulizumab pegol arm (which",
        "uses a discrete effect via CONMED_XI_DISCRETE)."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. Continuous Emax dose-response",
        "for the lulizumab pegol s.c. q2w regimen across three doses",
        "(1.25, 5, 12.5 mg; Goteti 2024 Table 1). Enters xi_ij via",
        "(Dose * Emax) / (Dose + ED50) with Emax_LULI = 0.165 and",
        "ED50_LULI = 9.47 mg (Table 4). The 12.5 mg s.c. q.w. lulizumab",
        "pegol regimen has its own discrete-effect value 0.206 (Table",
        "4); pass that as CONMED_XI_DISCRETE."
      ),
      source_name        = "Lulizumab pegol dose s.c. q2w (Goteti 2024 Table 1 / Table 4)"
    ),
    DOSE_SIFALIMUMAB_MG = list(
      description        = paste(
        "Per-arm sifalimumab dose amount (mg IV q2w x 3 + q4w); 0 if",
        "the arm did not receive sifalimumab."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-arm-level covariate. Continuous Emax dose-response",
        "for the sifalimumab i.v. q2wx3+q4w regimen across three doses",
        "(200, 600, 1200 mg; Goteti 2024 Table 1). Enters xi_ij via",
        "(Dose * Emax) / (Dose + ED50) with Emax_SIFA = 0.0131 and",
        "ED50_SIFA = 942 mg (Table 4). Effect size is small and the",
        "posterior interval crosses zero on Emax."
      ),
      source_name        = "Sifalimumab dose (Goteti 2024 Table 1 / Table 4)"
    ),
    CONMED_XI_DISCRETE = list(
      description        = paste(
        "Per-arm latent-scale discrete-treatment effect xi passed",
        "directly for discrete-effect regimens (Goteti 2024 Table 4);",
        "0 for placebo + SoC arms and for arms whose drug is covered",
        "by a continuous-Emax dose-response above."
      ),
      units              = "(latent-scale relative delta multiplier; dimensionless)",
      type               = "continuous",
      reference_category = "0 (placebo + standard of care)",
      notes              = paste(
        "MBMA study-arm-level covariate documenting the arm's",
        "discrete-effect xi value. The paper fit per-regimen additive",
        "effects on the latent-scale (1 + xi) multiplier for 24 drug-",
        "dose-regimen combinations that did not have three-plus doses",
        "at the same route and frequency (Goteti 2024 Table 4). Users",
        "reproduce a specific arm by passing that arm's tabulated xi",
        "value here.",
        "",
        "Lookup table (Goteti 2024 Table 4; values are posterior means):",
        "  ALX-0061 150 mg s.c. q2w:                -0.148",
        "  ALX-0061 150 mg s.c. q4w:                -0.11",
        "  ALX-0061 225 mg s.c. q2w:                -0.0935",
        "  ALX-0061 75  mg s.c. q4w:                +0.269",
        "  Atacicept 150 mg s.c. q.w.:              +0.00109",
        "  Atacicept 75  mg s.c. q.w.:              -0.0347",
        "  Baricitinib 2 mg oral q.d.:              +0.131",
        "  Baricitinib 4 mg oral q.d.:              +0.24",
        "  Belimumab 200 mg s.c. q.w.:              +0.16",
        "  Epratuzumab 600 mg i.v. q.w.:            +0.324",
        "  Evobrutinib 25 mg oral q.d.:             +0.352",
        "  Evobrutinib 50 mg oral b.i.d.:           +0.213",
        "  Evobrutinib 75 mg oral q.d.:             +0.314",
        "  Fenebrutinib 150 mg oral q.d.:           +0.0978",
        "  Fenebrutinib 200 mg oral b.i.d.:         +0.00941",
        "  IPP-201101 200 mcg s.c. q4w:             +0.134",
        "  Lulizumab pegol 12.5 mg s.c. q.w.:       +0.206",
        "  Lupuzor 200 mg s.c. q2w:                 +0.12",
        "  Lupuzor 200 mg s.c. q4w:                 +0.205",
        "  (Table 4 prints 201 mg; the source Table 1 uses 200 mg -",
        "   Lupuzor is dosed at 200 mg. Treat 201 as a Table 4 typo.)",
        "  PF-04236921 10 mg s.c. q8wx3:            +0.351",
        "  PF-04236921 50 mg s.c. q8wx3:            -0.101",
        "  Tabalumab 120 mg s.c. q2w:               +0.0699",
        "  Tabalumab 120 mg s.c. q4w:               +0.0886",
        "  Ustekinumab 6 mg/kg i.v. once + 90 mg s.c. q8w: -0.182",
        "",
        "Placebo + SoC arms and continuous-Emax arms pass",
        "CONMED_XI_DISCRETE = 0."
      ),
      source_name        = "xi_delta_d (Goteti 2024 Table 4 discrete-effect column)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Mean age of study-arm patients (years).",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Reported per-arm (Goteti 2024 Table 2: mean 40.2, range",
        "31.7-46.4 across arms) but not retained as a covariate in the",
        "final MBMA DTM. Excluded to preserve the screening provenance",
        "without triggering an unreferenced-covariate warning."
      )
    ),
    WT = list(
      description = "Mean body weight of study-arm patients (kg).",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Reported per-arm for 17% of arms (Goteti 2024 Table 2: mean",
        "70.6, range 66.2-75.4) but not retained as a covariate in the",
        "final MBMA DTM."
      )
    ),
    SEXF = list(
      description = "Percentage female per study arm (%).",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Reported per-arm (Goteti 2024 Table 2: mean 94.1%, range",
        "83.3-100%) but not retained as a covariate in the final MBMA",
        "DTM."
      )
    ),
    RACE_WHITE = list(
      description = "Percentage White race per study arm (%).",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Reported per-arm for 85% of arms (Goteti 2024 Table 2: mean",
        "64.4%, range 0-100%) but not retained as a covariate in the",
        "final MBMA DTM. Similarly for Black, Asian, and other race",
        "percentages (Table 2)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12136L,
    n_studies      = 25L,
    n_arms         = 81L,
    n_treatments   = 16L,
    age_range      = paste(
      "adult; per-arm mean age range 31.7-46.4 years (Goteti 2024 Table",
      "2), overall across-arm mean 40.2 years"
    ),
    age_median     = "across-arm median 40.6 years (Goteti 2024 Table 2)",
    weight_range   = paste(
      "per-arm mean weight range 66.2-75.4 kg reported in 17% of arms",
      "(Goteti 2024 Table 2), overall across-arm mean 70.6 kg"
    ),
    sex_female_pct = 94.1,
    race_ethnicity = c(
      White = 64.4,
      Asian = 13.6,
      Black = 11.2,
      Other = 14.9
    ),
    regions        = paste(
      "International; North America 34.2%, Central-South America 19.5%,",
      "Asia 17.5%, other regions 27.0% (Goteti 2024 Table 2 mean across",
      "arms that reported region)"
    ),
    disease_state  = paste(
      "Systemic lupus erythematosus (SLE); across-arm mean baseline",
      "SLEDAI 10.5 (range 8.8-14.0), physician global assessment 1.46,",
      "prednisone-equivalent daily dose 11.7 mg/day, low C3 in 37.1%,",
      "low C4 in 29.6%, anti-dsDNA positive in 55.1%, on concurrent",
      "antimalarials 66.3%, on concurrent immunosuppressants 46.1%",
      "(Goteti 2024 Table 2)"
    ),
    dose_range     = paste(
      "16 investigational and approved SLE therapeutics across a range",
      "of doses, routes, and frequencies (Goteti 2024 Table 1). See",
      "the CONMED_XI_DISCRETE lookup and the DOSE_* covariate notes for",
      "the six continuous-Emax regimens for the per-drug dose ranges."
    ),
    notes          = paste(
      "MBMA at the study-arm level: each modelled data point is the",
      "aggregate mean response rate (SRI-4/5/6 or BICLA) in a group of",
      "patients at a particular visit time in one arm of one trial,",
      "with total sample size N_ij(t) used as the binomial /",
      "multinomial denominator in the joint likelihood (Goteti 2024",
      "supplement Eq 11-13). The model is intended for simulating",
      "study-arm-mean responder probabilities and is NOT suitable for",
      "individual-subject simulation. The paper fit the model in Stan",
      "(Bayesian, NUTS sampler; Stan Dev. Team v2.24) with 1000+",
      "effective posterior samples for all reported parameters. Table 3",
      "reports posterior means of structural parameters and Table 4",
      "reports posterior means of treatment effects; all values encoded",
      "here are posterior means."
    )
  )

  ini({
    # ============================================================
    # Disease-trajectory model (Goteti 2024 Eq 1, Eq 8):
    #
    #   theta_ijk(t) = mu_ijk*(1 + beta_sledai_mu * sledai_dev)
    #                + delta_ijk*(1 + beta_sledai_delta * sledai_dev)
    #                            *(1 + xi_ij)
    #                            *(1 - exp(-t * lam))
    #
    # Reparametrized (Goteti 2024 supplement Eq 3, 5):
    #   mu_ijk    = theta_wk4 - delta_ijk * (1 - exp(-4  * lam))
    #   delta_ijk = delta_wk24 / (exp(-4 * lam) - exp(-24 * lam))
    #
    # theta_wk4_k and (theta_wk4_k + delta_wk24_k) are reported on the
    # probability scale in Goteti 2024 Table 3 (footnote), but the
    # underlying Bayesian estimation is on the latent (logit) scale on
    # which the covariate / drug / study effects operate. Store the
    # latent-scale values here and note the probability-scale
    # equivalent in each label. Study effects (etas below) also live on
    # the latent scale.
    # ============================================================

    theta_wk4_SRI <- log(0.170 / (1 - 0.170))
    label(paste(
      "SRI: latent SLE disease at week 4, placebo + SoC, SLEDAI = 10.5",
      "(logit scale; probability-scale equivalent 0.170 per Table 3)"
    ))
    # Goteti 2024 Table 3 theta_wk4,SRI = 0.170 (90% CI 0.147, 0.197);
    # logit(0.170) = log(0.170/0.830) = -1.5872.

    delta_wk24_SRI <- log(0.469 / (1 - 0.469)) -
                      log(0.170 / (1 - 0.170))
    label(paste(
      "SRI: latent-scale increment from week 4 to week 24, placebo +",
      "SoC, SLEDAI = 10.5 (logit scale; probability-scale wk-24 value",
      "0.469 per Table 3 for theta_wk4,SRI + delta_wk24,SRI)"
    ))
    # Table 3 theta_wk4,SRI + delta_wk24,SRI = 0.469 (90% CI 0.419,
    # 0.520); on latent scale = logit(0.469) - logit(0.170) = -0.1241 -
    # (-1.5872) = 1.4631.

    theta_wk4_BICLA <- log(0.247 / (1 - 0.247))
    label(paste(
      "BICLA: latent SLE disease at week 4, placebo + SoC, SLEDAI =",
      "10.5 (logit scale; probability-scale equivalent 0.247 per Table",
      "3)"
    ))
    # Table 3 theta_wk4,BICLA = 0.247 (90% CI 0.233, 0.262);
    # logit(0.247) = log(0.247/0.753) = -1.1155.

    delta_wk24_BICLA <- log(0.412 / (1 - 0.412)) -
                        log(0.247 / (1 - 0.247))
    label(paste(
      "BICLA: latent-scale increment from week 4 to week 24, placebo +",
      "SoC, SLEDAI = 10.5 (logit scale; probability-scale wk-24 value",
      "0.412 per Table 3 for theta_wk4,BICLA + delta_wk24,BICLA)"
    ))
    # Table 3 theta_wk4,BICLA + delta_wk24,BICLA = 0.412 (90% CI 0.354,
    # 0.468); on latent scale = logit(0.412) - logit(0.247) = -0.3557 -
    # (-1.1155) = 0.7598.

    lthalf <- log(4.51)
    label(paste(
      "Log of the mono-exponential half-life to the long-term latent",
      "disease level (weeks); back-transform thalf = 4.51 wk"
    ))
    # Table 3 t_1/2 = 4.51 wk (90% CI 4.26, 4.78); the rate constant is
    # lam = log(2) / thalf, computed in model().

    # ============================================================
    # SLEDAI covariate (Goteti 2024 Eq 8) - multiplicative on latent
    # scale.
    # ============================================================
    beta_sledai_mu <- 0.208
    label(paste(
      "SLEDAI relative-change coefficient on latent mu_ijk (per unit",
      "SLEDAI above 10.5)"
    ))
    # Table 3 beta_SLEDAI,0 = 0.208 (90% CI 0.165, 0.251).

    beta_sledai_delta <- 0.255
    label(paste(
      "SLEDAI relative-change coefficient on latent delta_ijk (per",
      "unit SLEDAI above 10.5)"
    ))
    # Table 3 beta_SLEDAI,delta = 0.255 (90% CI 0.210, 0.302).

    # ============================================================
    # Generalized proportional-odds thresholds for SRI (Goteti 2024
    # supplement Model for SRI):
    #   p_SRI-l(t) = expit(theta_SRI(t) - c_{l-4})     l = 4, 5, 6
    #   with c_0 = 0, c_1 < c_2 (equivalently, the thresholds shift the
    #   latent activity to yield the stricter SRI-5 / SRI-6 responder
    #   probabilities).
    #
    # Table 3 reports c_2,SRI = 0.586 and c_3,SRI = 0.693; these
    # correspond to the SRI-5 and SRI-6 thresholds respectively in the
    # supplement's indexing (Table 3 uses a 1-based label; the
    # supplement uses a 0-based label - see vignette Assumptions and
    # deviations for the mapping verification).
    # ============================================================
    c_SRI5 <- 0.586
    label("SRI-5 proportional-odds threshold (latent-scale shift; larger = stricter)")
    # Table 3 c_2,SRI = 0.586 (90% CI 0.572, 0.600).

    c_SRI6 <- 0.693
    label("SRI-6 proportional-odds threshold (latent-scale shift; larger = stricter)")
    # Table 3 c_3,SRI = 0.693 (90% CI 0.678, 0.707).

    # ============================================================
    # Continuous-Emax drug parameters (Goteti 2024 Table 4).
    # xi_d_continuous = (Dose * Emax_d) / (Dose + ED50_d)  on the
    # latent-scale delta multiplier (1 + xi).
    #
    # Note the ED50 values live at very different scales across drugs
    # (mg / mg/kg / very-small mg), so each is log-transformed here.
    # ============================================================
    Emax_ANIF <- 0.257
    label("Anifrolumab i.v. q4w Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,Anifrolumab = 0.257 (95% CrI 0.173, 0.370).

    lED50_ANIF <- log(259)
    label("Log anifrolumab i.v. q4w ED50 (mg); back-transform ED50 = 259 mg")
    # Table 4 ED50 Anifrolumab (mg) = 259 (95% CrI 70.4, 544).

    Emax_BEL <- 0.165
    label("Belimumab i.v. q2wx3+q4w Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,Belimumab = 0.165 (95% CrI 0.136, 0.196).

    lED50_BEL <- log(0.787)
    label("Log belimumab i.v. q2wx3+q4w ED50 (mg/kg); back-transform ED50 = 0.787 mg/kg")
    # Table 4 ED50 Belimumab (mg/kg) = 0.787 (95% CrI 0.334, 1.52).

    Emax_CC220 <- 0.415
    label("CC-220 (iberdomide) oral qd Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,CC220 = 0.415 (95% CrI -0.0856, 1.02).

    lED50_CC220 <- log(0.415)
    label("Log CC-220 oral qd ED50 (mg); back-transform ED50 = 0.415 mg")
    # Table 4 ED50 CC220 (mg) = 0.415 (95% CrI 0.0168, 1.14).

    Emax_EPRA <- 0.397
    label("Epratuzumab i.v. q2w Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,Epratuzumab = 0.397 (95% CrI 0.164, 0.749).

    lED50_EPRA <- log(827)
    label("Log epratuzumab i.v. q2w ED50 (mg); back-transform ED50 = 827 mg")
    # Table 4 ED50 Epratuzumab (mg) = 827 (95% CrI 49.3, 2140).

    Emax_LULI <- 0.165
    label("Lulizumab pegol s.c. q2w Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,Lulizumab pegol = 0.165 (95% CrI -0.0873, 0.530).

    lED50_LULI <- log(9.47)
    label("Log lulizumab pegol s.c. q2w ED50 (mg); back-transform ED50 = 9.47 mg")
    # Table 4 ED50 Lulizumab pegol (mg) = 9.47 (95% CrI 0.516, 24.4).

    Emax_SIFA <- 0.0131
    label("Sifalimumab i.v. q2wx3+q4w Emax (latent-scale relative-delta multiplier)")
    # Table 4 Emax,Sifalimumab = 0.0131 (95% CrI -0.0701, 0.112).

    lED50_SIFA <- log(942)
    label("Log sifalimumab i.v. q2wx3+q4w ED50 (mg); back-transform ED50 = 942 mg")
    # Table 4 ED50 Sifalimumab (mg) = 942 (95% CrI 30.7, 2430).

    # ============================================================
    # Between-study variability (Goteti 2024 Table 3). These are
    # STUDY-level random-effect standard deviations on the latent-scale
    # theta_wk4_k and delta_wk24_k, applied additively per the source
    # supplement Eq 6-7. There are NO subject-level etas in the source
    # (this is aggregate summary-level MBMA); the 'eta' names below
    # therefore reflect between-study variability, distinct from the
    # popPK between-subject variability convention.
    #
    # Diagonals (variances):
    #   var(eta_wk4_SRI)    = sigma_wk4,SRI^2   = 0.365^2  = 0.133225
    #   var(eta_delta_SRI)  = sigma_delta,SRI^2 = 0.290^2  = 0.0841
    #   var(eta_wk4_BICLA)  = sigma_wk4,BICLA^2 = 0.102^2  = 0.010404
    #   var(eta_delta_BICLA)= sigma_delta,BICLA^2 = 0.524^2 = 0.274576
    #
    # The four etas are modelled as independent (Goteti 2024 Methods:
    # "no off-diagonal elements of variance matrix were used between
    # the random effects").
    # ============================================================
    eta_wk4_SRI     ~ 0.133225   # Table 3 sigma_4w,SRI = 0.365 (CI 0.252, 0.538); variance = 0.365^2.
    eta_delta_SRI   ~ 0.0841     # Table 3 sigma_delta,SRI = 0.290 (CI 0.161, 0.487); variance = 0.290^2.
    eta_wk4_BICLA   ~ 0.010404   # Table 3 sigma_4w,BICLA = 0.102 (CI 0.0533, 0.185); variance = 0.102^2.
    eta_delta_BICLA ~ 0.274576   # Table 3 sigma_delta,BICLA = 0.524 (CI 0.338, 0.803); variance = 0.524^2.

    # ============================================================
    # Residual error placeholder. The paper's likelihood is a binomial
    # / multinomial at the study-arm level (supplement Eq 11-13), NOT a
    # continuous residual on the response probability. rxode2 / nlmixr2
    # do not natively encode binomial / multinomial residuals for
    # aggregate response probabilities, so we expose the study-arm-mean
    # responder probability as a continuous algebraic output and attach
    # a small additive residual placeholder to satisfy the nlmixr2
    # observation-model requirement. The stochastic scope of the model
    # comes from the study-level etas above, not from this residual.
    # See vignette Assumptions and deviations for the deviation
    # rationale.
    # ============================================================
    addSd <- fixed(0.001)
    label(paste(
      "Additive residual placeholder on the study-arm-mean responder",
      "probability (FIXED; the source paper's likelihood is binomial /",
      "multinomial at the arm level, not encoded here)"
    ))
  })

  model({
    # -------- Rate constant from the half-life --------
    thalf <- exp(lthalf)
    lam   <- log(2) / thalf

    # -------- Reparametrization exponents (Eq 3, 5 of supplement) --------
    a4  <- exp(-4  * lam)
    a24 <- exp(-24 * lam)

    # -------- SLEDAI covariate (centered at 10.5) ------
    # Goteti 2024 Eq 8: the covariate x_ij is the mean baseline SLEDAI
    # centered at 10.5; enters multiplicatively as (1 + beta * x_ij) on
    # both mu and delta.
    sledai_dev <- SCORE_SLEDAI - 10.5
    sledai_mu_fx    <- 1 + beta_sledai_mu    * sledai_dev
    sledai_delta_fx <- 1 + beta_sledai_delta * sledai_dev

    # -------- Continuous-Emax drug effects (Eq 9) ------
    # Each per-drug xi is (Dose * Emax) / (Dose + ED50); zero when the
    # corresponding dose column is 0 (drug not in the arm). Because
    # only ONE drug is administered per study arm in the source data,
    # the sum collapses to the single active drug. If a downstream user
    # supplies non-zero doses for two continuous drugs simultaneously,
    # the model returns the ADDITIVE sum of their xi values, which is
    # outside the paper's calibration range.
    ED50_ANIF  <- exp(lED50_ANIF)
    ED50_BEL   <- exp(lED50_BEL)
    ED50_CC220 <- exp(lED50_CC220)
    ED50_EPRA  <- exp(lED50_EPRA)
    ED50_LULI  <- exp(lED50_LULI)
    ED50_SIFA  <- exp(lED50_SIFA)

    xi_anif  <- (DOSE_ANIFROLUMAB_MG      * Emax_ANIF ) /
                (DOSE_ANIFROLUMAB_MG      + ED50_ANIF )
    xi_bel   <- (DOSE_BELIMUMAB_MGKG      * Emax_BEL  ) /
                (DOSE_BELIMUMAB_MGKG      + ED50_BEL  )
    xi_cc220 <- (DOSE_CC220_MG            * Emax_CC220) /
                (DOSE_CC220_MG            + ED50_CC220)
    xi_epra  <- (DOSE_EPRATUZUMAB_IVQ2W_MG* Emax_EPRA ) /
                (DOSE_EPRATUZUMAB_IVQ2W_MG+ ED50_EPRA )
    xi_luli  <- (DOSE_LULIZUMAB_SCQ2W_MG  * Emax_LULI ) /
                (DOSE_LULIZUMAB_SCQ2W_MG  + ED50_LULI )
    xi_sifa  <- (DOSE_SIFALIMUMAB_MG      * Emax_SIFA ) /
                (DOSE_SIFALIMUMAB_MG      + ED50_SIFA )

    xi_ij <- xi_anif + xi_bel + xi_cc220 + xi_epra + xi_luli + xi_sifa +
             CONMED_XI_DISCRETE

    # -------- Study-arm latent parameters (Eq 6, 7 with study etas) ---
    # SRI endpoint
    theta_wk4_SRI_i    <- theta_wk4_SRI    + eta_wk4_SRI
    delta_wk24_SRI_i   <- delta_wk24_SRI   + eta_delta_SRI
    delta_SRI          <- delta_wk24_SRI_i / (a4 - a24)
    mu_SRI             <- theta_wk4_SRI_i  - delta_SRI * (1 - a4)

    # BICLA endpoint
    theta_wk4_BICLA_i  <- theta_wk4_BICLA  + eta_wk4_BICLA
    delta_wk24_BICLA_i <- delta_wk24_BICLA + eta_delta_BICLA
    delta_BICLA        <- delta_wk24_BICLA_i / (a4 - a24)
    mu_BICLA           <- theta_wk4_BICLA_i  - delta_BICLA * (1 - a4)

    # -------- DTM with SLEDAI covariate + drug effect (Eq 8) ----------
    # theta_k(t) = mu_k*(1 + beta_mu*x) + delta_k*(1 + beta_delta*x)*(1 + xi)*(1 - exp(-t*lam))
    theta_SRI_t   <- mu_SRI   * sledai_mu_fx +
                     delta_SRI   * sledai_delta_fx * (1 + xi_ij) *
                     (1 - exp(-time * lam))

    theta_BICLA_t <- mu_BICLA * sledai_mu_fx +
                     delta_BICLA * sledai_delta_fx * (1 + xi_ij) *
                     (1 - exp(-time * lam))

    # -------- Responder probabilities (supplement Section 1.3) -------
    # SRI: generalized proportional-odds transform with c_0 = 0.
    p_SRI4 <- exp(theta_SRI_t             ) / (1 + exp(theta_SRI_t             ))
    p_SRI5 <- exp(theta_SRI_t - c_SRI5    ) / (1 + exp(theta_SRI_t - c_SRI5    ))
    p_SRI6 <- exp(theta_SRI_t - c_SRI6    ) / (1 + exp(theta_SRI_t - c_SRI6    ))
    # BICLA: logistic transform of the latent SLE activity.
    p_BICLA <- exp(theta_BICLA_t) / (1 + exp(theta_BICLA_t))

    # -------- Observation ---------------------------------------------
    # A single Cc endpoint is declared for nlmixr2 compatibility (the
    # p_SRI-4 responder probability). All four algebraic outputs
    # (p_SRI4, p_SRI5, p_SRI6, p_BICLA) are exposed as derived model
    # variables and are emitted in rxSolve output; downstream users
    # read them all from that data frame. Cc = p_SRI-4 by convention;
    # see vignette Assumptions and deviations for the choice.
    Cc <- p_SRI4
    Cc ~ add(addSd)
  })
}
