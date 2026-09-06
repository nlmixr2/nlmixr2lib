Sethi_2024_acetaminophen_diclofenac_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis of the opioid-sparing effect of combined acetaminophen and diclofenac in acute postoperative pain, fit to arm-level summary data from five randomized controlled trials (353 adults, 16 treatment arms) that permitted opioid patient-controlled analgesia (PCA). The endpoint is the arm-mean cumulative opioid PCA consumption (mg) over each trial's primary observation window; opioid use is described as a typical placebo response plus an additive acetaminophen effect, an additive diclofenac effect, and a multiplicative interaction term: Cc = e0 + f(ace) + f(dic) + gamma * f(ace) * f(dic) (paper Eqs 1-2). The interaction coefficient is positive (gamma = 0.025 /mg), i.e. the combination is SUB-additive: the combined opioid-sparing effect (-33.9 mg) is smaller in magnitude than the sum of the two monotherapy effects (-47.3 mg). IMPORTANT SCOPE LIMITS. (1) The drug effects are DOSE-INDEPENDENT binary study-arm indicators, not dose-response functions: the source trials spanned acetaminophen 1000-2400 mg and diclofenac 75-100 mg but were too few to identify a dose-response, so every acetaminophen arm receives the same -18.93 mg effect regardless of dose. (2) Despite the paper's title, every diclofenac arm in the fitted dataset used SYSTEMIC diclofenac (oral, rectal or intravenous); the extrapolation to TOPICAL diclofenac and to chronic osteoarthritis pain is a qualitative argument in the paper's Discussion and is NOT encoded here. (3) Opioid PCA is pooled in raw mg across morphine (four trials) and oxycodone (one trial) with no potency normalization. (4) The trial-specific placebo response was estimated by an unstructured (non-parametric) per-trial model whose individual estimates are not reported; e0 is fixed to the single typical placebo response (64.7 mg) the paper used for its own simulations, and only one of the five trials carried a placebo arm. (5) There is no time course, no dose event, and no PK: the model is a static algebraic MBMA suitable for simulating study-arm-mean opioid PCA consumption, NOT individual-patient predictions."

  reference <- paste(
    "Sethi V, Qin L, Cox E, Troconiz IF, Della Pasqua O.",
    "Model-Based Meta-Analysis Supporting the Combination of Acetaminophen and",
    "Topical Diclofenac in Acute Pain: A Therapy for Mild-to-Moderate Osteoarthritis Pain?",
    "Pain Ther. 2024 Feb;13(1):145-159.",
    "doi:10.1007/s40122-023-00569-z.",
    sep = " "
  )
  vignette <- "Sethi_2024_acetaminophen_diclofenac_mbma"

  units <- list(
    time          = "hour (placeholder; the endpoint is cumulative opioid PCA consumption over each trial's primary observation window of roughly 6-72 h and the model itself is time-independent)",
    dosing        = "mg (acetaminophen and diclofenac doses are NOT model inputs; both drug effects are dose-independent binary study-arm indicators. The fitted trials used acetaminophen 1000-2400 mg and diclofenac 75-100 mg. The model does not consume rxode2 dose events)",
    concentration = "mg/mg (arm-mean cumulative opioid PCA consumption in mg of morphine or oxycodone; the output Cc is NOT a drug concentration. The slash in the unit string satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    ACETAMINOPHEN = list(
      description        = "Binary study-arm treatment indicator: 1 = the arm received systemic acetaminophen, 0 = it did not. A property of the trial arm in a model-based meta-analysis, not of an individual patient.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (arm did not receive acetaminophen)",
      notes              = "MBMA study-arm-level treatment indicator, a family-conforming member of the bare-INN MBMA arm-indicator family established by NAPROXEN (Boucher_2018_naproxen_mbma) and TRAMADOL / TAPENTADOL (Mercier_2014_tramadol_tapentadol_mbma). Enters the drug-effect term as e_acetaminophen_e0 * ACETAMINOPHEN. NOT dose-dependent: Sethi 2024 could not identify an acetaminophen dose-response from five trials, so arms at 1000, 1200, 1500, 2000 and 2400 mg all receive the identical effect. Note in particular that Beck 2000 contributed both a 1200 mg and a 2400 mg acetaminophen arm and the model predicts the same opioid PCA use for both. Set to 1 on the combination arms as well as the acetaminophen-monotherapy arms.",
      source_name        = "ace (Sethi 2024 Table 1 Treatments column; e.ace in Table 3)"
    ),
    DICLOFENAC = list(
      description        = "Binary study-arm treatment indicator: 1 = the arm received systemic diclofenac, 0 = it did not. A property of the trial arm in a model-based meta-analysis, not of an individual patient.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (arm did not receive diclofenac)",
      notes              = "MBMA study-arm-level treatment indicator; sibling of ACETAMINOPHEN and a family-conforming member of the bare-INN MBMA arm-indicator family (NAPROXEN / TRAMADOL / TAPENTADOL). Enters the drug-effect term as e_diclofenac_e0 * DICLOFENAC. SYSTEMIC route only: every diclofenac arm in the five fitted trials used oral, rectal or intravenous diclofenac at 75-100 mg (Sethi 2024 Table 1 ROA column), and Supplementary Table S5 labels the rows 'Systemic diclofenac'. The paper's title refers to TOPICAL diclofenac because the Discussion extrapolates the systemic-diclofenac result to a topical osteoarthritis setting on mechanistic grounds; that extrapolation is not part of the fitted model. NOT dose-dependent (75 mg and 100 mg arms receive the identical effect). Set to 1 on the combination arms as well as the diclofenac-monotherapy arms.",
      source_name        = "dic (Sethi 2024 Table 1 Treatments column; e.dic in Table 3)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 353L,
    n_studies      = 5L,
    n_arms         = 16L,
    age_range      = "adults; per-trial age distributions are not reported in the meta-analysis",
    weight_range   = "not reported at arm level",
    race_ethnicity = "not reported",
    disease_state  = "acute postoperative pain: elective gynecological surgery (Montgomery 1996, n = 59), hysterectomy (Beck 2000, n = 65), cesarean section (Siddik 2001, n = 80; Munishankar 2008, n = 78) and tonsillectomy (Hiller 2004, n = 71)",
    dose_range     = "acetaminophen 1000-2400 mg and diclofenac 75-100 mg across the five included trials (Sethi 2024 Table 1). Doses are not model inputs -- see the ACETAMINOPHEN / DICLOFENAC covariate notes",
    regions        = "not reported",
    notes          = "Study-arm-level MBMA: each modeled data point is the mean cumulative opioid PCA consumption in one trial arm. The five included trials are the subset of Sethi 2024 Table 2 marked 'Included in the final analysis = Yes'; their sample sizes sum to exactly the 353 patients the Abstract reports. Two further PCA-reporting trials in Table 2 were excluded: Breivik 1999 (n = 72; only a limited percentage of subjects used PCA in each arm, and PCA was reported as a percentage rather than in mg) and Riad 2007 (n = 108; a pediatric population). The 16-arm count is derived by counting the Treatments column of Table 1 for the five included trials (Montgomery 3, Beck 3, Siddik 4, Hiller 3, Munishankar 3); the paper itself does not print an arm count. At least 80 percent of subjects were female -- four of the five trials (Montgomery, Beck, Siddik, Munishankar; 282 of 353 subjects) enrolled women only, and the sex distribution of the Hiller 2004 tonsillectomy cohort (n = 71 adults) is not reported. Opioid PCA was morphine in four trials and oxycodone in one (Hiller 2004) and is pooled in raw mg with no potency normalization. Only Siddik 2001 carried a placebo arm, so a single trial informs the whole placebo response; the paper flags this as a source of estimation bias in its Discussion. Type of surgical intervention was recognised as a driver of opioid PCA use (cesarean pain required a higher PCA dose than tonsillectomy) but could not be included as a covariate because too few trials were available. Residual variance was fixed to the observed per-arm precision sigma_ij^2 / N_ij rather than estimated -- see the addSd note in ini()."
  )

  ini({
    # ==================================================================
    # Sethi 2024 opioid-sparing-effect MBMA (Eqs 1-2, Methods p148):
    #
    #   dY_ij,ose = eo_i,ose + f(Drug_ij, theta) + e_ij,ose            (1)
    #   f(Drug)_ose = f(acet) + f(diclof) + gamma * f(acet) * f(diclof) (2)
    #
    # Fit with the gnls() generalized-least-squares function of the R
    # nlme package (Methods, Model Evaluation).
    # ==================================================================

    e0 <- fixed(64.7)
    label("Typical placebo response: arm-mean cumulative opioid PCA consumption (mg) in the absence of acetaminophen or diclofenac. Not estimated -- the paper modelled the trial-specific placebo response eo_i,ose with an unstructured (non-parametric) per-trial model and does not report the individual per-trial estimates. 64.7 mg is the single typical placebo response the paper itself assumed for the Fig 2 / Table S5 simulations, and it reproduces the Table S5 percent-difference column to two decimal places (see the vignette source-trace table).")  # Sethi 2024 Fig 2 caption: "assuming a typical placebo response (64.7 mg)"

    e_acetaminophen_e0 <- -18.93
    label("Additive acetaminophen effect on arm-mean opioid PCA consumption (mg). Negative = opioid sparing. Applied as e_acetaminophen_e0 * ACETAMINOPHEN; dose-independent (see covariateData). Estimated, RSE 29 percent; Table 3 asymptotic 95 percent CI -31.4 to -6.43 mg, Table S5 resampling 95 percent CI -29.47 to -8.33 mg.")  # Sethi 2024 Table 3 e.ace = -18.9; Supplementary Table S5 prints the same quantity to one more significant figure as -18.93 and that higher-precision value is used here

    e_diclofenac_e0 <- -28.41
    label("Additive systemic-diclofenac effect on arm-mean opioid PCA consumption (mg). Negative = opioid sparing. Applied as e_diclofenac_e0 * DICLOFENAC; dose-independent (see covariateData). Estimated, RSE 19 percent; Table 3 asymptotic 95 percent CI -40.7 to -16.2 mg, Table S5 resampling 95 percent CI -39.15 to -17.68 mg.")  # Sethi 2024 Table 3 e.dic = -28.4; Supplementary Table S5 prints the same quantity to one more significant figure as -28.41 and that higher-precision value is used here

    gamma_ad <- 0.025
    label("Acetaminophen x diclofenac interaction coefficient (1/mg) in Eq 2. POSITIVE means the combination is sub-additive: the combined effect is smaller in magnitude than the sum of the two monotherapy effects. A value of 0 would mean exact additivity and a negative value a synergistic (more-than-additive) combination. Estimated, RSE 18 percent, 95 percent CI 0.0148 to 0.0353 /mg; adding this term to the additive model improved the fit with p = 0.028 (Supplementary Table S4: AIC 115 -> 113).")  # Sethi 2024 Table 3 c = 0.025 [0.0148 to 0.0353]

    # ==================================================================
    # Residual error. Sethi 2024 Methods define the within-trial
    # residual as e_ij,ose ~ N(0, sigma_ij^2 / N_ij), where sigma_ij is
    # the OBSERVED standard deviation of the outcome in arm j of trial i
    # and N_ij is that arm's sample size -- i.e. the residual variance is
    # supplied by the data as the squared standard error of each arm
    # mean, and NO residual magnitude is estimated by the model. Table 3
    # accordingly reports no sigma. addSd is therefore fixed to zero and
    # downstream simulation code must supply the per-arm standard error
    # explicitly if arm-level noise is wanted (same downstream-weighting
    # convention as Mercier_2014_tramadol_tapentadol_mbma and
    # Vargo_2014_statins_ezetimibe_mbma, which do estimate a sigma).
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual SD on arm-mean opioid PCA consumption (mg); ZERO. The source model does not estimate a residual magnitude: its within-trial residual variance is the observed per-arm squared standard error sigma_ij^2 / N_ij (Methods, Model Development). See the ini() comment above and the vignette Assumptions and deviations section.")  # Sethi 2024 Methods: var(e_ij,ose) = sigma_ij^2 / N_ij, data-supplied; no sigma is reported in Table 3
  })

  model({
    # Study-arm covariates supplied per row:
    #   ACETAMINOPHEN -- 1 if the arm received systemic acetaminophen, else 0
    #   DICLOFENAC    -- 1 if the arm received systemic diclofenac, else 0
    # Both are 1 on a combination arm and both are 0 on a placebo arm.

    # 1. Monotherapy drug effects, gated by the arm indicators (Eq 2
    #    f(acet) and f(diclof)). Each is dose-independent.
    f_acetaminophen <- e_acetaminophen_e0 * ACETAMINOPHEN
    f_diclofenac    <- e_diclofenac_e0 * DICLOFENAC

    # 2. Combination drug effect (Eq 2). The interaction product is
    #    non-zero only when BOTH indicators are 1, so this single
    #    expression covers the placebo, both monotherapy and the
    #    combination arms without branching.
    f_drug <- f_acetaminophen + f_diclofenac +
      gamma_ad * f_acetaminophen * f_diclofenac

    # 3. Arm-mean opioid PCA consumption (Eq 1). f_drug is also the
    #    paper's "placebo-adjusted opioid PCA use" plotted in Fig 2 and
    #    tabulated in Supplementary Table S5.
    Cc <- e0 + f_drug

    # 4. Residual is fixed at zero -- see the ini() comment. Per-arm
    #    standard errors are applied downstream.
    Cc ~ add(addSd)
  })
}
