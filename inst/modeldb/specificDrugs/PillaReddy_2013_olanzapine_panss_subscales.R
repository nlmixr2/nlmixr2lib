# Population PK/PD model for olanzapine against the three PANSS subscales
# (positive, negative, general) in adults with schizophrenia (Pilla Reddy 2013
# Part II; Schizophrenia Research 146:153-161; doi:10.1016/j.schres.2013.02.010).

PillaReddy_2013_olanzapine_panss_subscales <- function() {
  description <- paste(
    "Population PK/PD model for olanzapine against the three PANSS",
    "subscales (positive, negative, general) in adults with schizophrenia",
    "from Pilla Reddy 2013 Part II. The PK sub-model is the one-compartment",
    "olanzapine structural model from Part I (PMID 23473810) Table 2:",
    "first-order absorption ka = 0.30 1/h, apparent oral clearance CL/F =",
    "21.8 L/h, apparent central volume of distribution Vc/F = 700 L. The",
    "PD sub-model has three outputs that share the Weibull placebo time-",
    "course form Pplacebo = Pmax * (1 - exp(-(t/TD)^POW)) but each subscale",
    "carries its own placebo Pmax, TD, POW (Part II Table 1) and",
    "olanzapine's own Emax / EC50 / KT triplet per subscale (Part II Table",
    "2). The KT for olanzapine PANSS positive and general (0.048 and 0.035",
    "1/day) is the common-across-atypical-antipsychotic value (Part II",
    "Methods); the KT for the negative subscale (0.028 1/day) was estimated",
    "separately per drug. Olanzapine was numerically superior to the other",
    "SGAs for the negative subscale (Emax = 0.33 vs 0.14-0.17 for the other",
    "SGAs; Part II Results), making this a clinically meaningful subscale-",
    "specific comparison. The exponential time-to-event dropout sub-model",
    "from Part II Table 4 is documented in population$dropout_model but is",
    "not encoded in this model body.",
    sep = " "
  )
  reference <- paste(
    "Pilla Reddy V, Kozielska M, Suleiman AA, Johnson M, Vermeulen A, Liu J,",
    "de Greef R, Groothuis GMM, Danhof M, Proost JH (2013).",
    "Pharmacokinetic-pharmacodynamic modelling of antipsychotic drugs in",
    "patients with schizophrenia: Part II: The use of subscales of the",
    "PANSS score. Schizophrenia Research 146(1-3):153-161.",
    "doi:10.1016/j.schres.2013.02.010. PK structure inherited from Part I",
    "(PMID 23473810; doi:10.1016/j.schres.2013.02.011).",
    sep = " "
  )
  vignette <- "PillaReddy_2013_panss_subscales"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part I Table 2 reports a 25% lower olanzapine CL/F in females vs males ('CL (gender) -0.25 (-0.323 to -0.173)'); not retained in the Part II PD model() body (typical-individual reference is the pooled male+female estimate)."
    ),
    DIS = list(
      description = "Disease state at entry (acute vs chronic schizophrenia)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part II placebo-model covariate; not implemented."
    ),
    USA = list(
      description = "Study geographic origin (USA vs non-USA)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part II placebo-model covariate on Pmax negative and residual error; not implemented."
    ),
    REG = list(
      description = "Dosing regimen (qd vs bid)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part II residual-error covariate on PANSS negative; not implemented."
    ),
    DUR = list(
      description = "Study duration (short vs long)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part II residual-error covariate on PANSS positive; not implemented."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 741L,
    n_studies       = 12L,
    age_range       = "Adults with schizophrenia (specific range not tabulated in Part II).",
    weight_range    = "Adult schizophrenia population.",
    sex_female_pct  = NA_real_,
    disease_state   = paste(
      "Adults with schizophrenia recruited into acute and chronic-stable",
      "Phase II / III double-blind clinical trials. Olanzapine arms in",
      "Part I Table 1: SCH-303 (10 mg qd), SCH-304 (10 mg qd), SCH-305 (10",
      "mg qd), SCH-2002 (15 mg qd). Baseline PANSS positive typical value",
      "22.4 (95% CI 22.4-22.6), PANSS negative 23.8 (23.5-24.2), PANSS",
      "general 44 (43.7-44.2) for the atypical-antipsychotic pool (Part",
      "II Table 2)."
    ),
    dose_range      = "Oral olanzapine 10-15 mg/day (Part I Table 1).",
    regions         = "Pooled across multinational schizophrenia trials 1989-2009.",
    notes           = paste(
      "Olanzapine PANSS-negative Emax (0.33) was numerically the highest",
      "across the five compared antipsychotics (Part II Figure 1a and",
      "Results), associated with olanzapine's broader receptor binding",
      "profile (D2, 5-HT2A, muscarinic, histaminergic, alpha-adrenergic;",
      "Part II Discussion)."
    ),
    dropout_model = paste(
      "Part II Table 4: olanzapine placebo BHAZ values are 0.00096",
      "(positive), 0.00304 (negative), 0.00054 (general) 1/day. BETA",
      "coefficients shared across drugs: -0.111 (positive), -0.0443",
      "(negative), -0.0672 (general). Not encoded in this model body."
    )
  )

  ini({
    # =============== Pharmacokinetic structural parameters ===================
    lka <- fixed(log(0.30));  label("Log first-order absorption rate constant (1/h)")              # Part I Table 2: Ka olanzapine = 0.30 1/h (95% CI not narrow; reported with 11% RSE)
    lcl <- fixed(log(21.8));  label("Log apparent oral clearance CL/F (L/h)")                       # Part I Table 2: CL/F = 21.8 L/h (95% CI 20.5-23.1)
    lvc <- fixed(log(700));   label("Log apparent central volume of distribution Vc/F (L)")         # Part I Table 2: Vc/F = 700 L (95% CI 560-814)

    # =============== Inter-individual variability on PK ======================
    # Part I Table 2: IIV CL not reported for olanzapine (the row is dashed),
    # IIV Vc 39% CV. IIV Ka not reported. Encode IIV Vc only.
    etalvc ~ fixed(0.142)  # Part I Table 2: IIV Vc 39% CV; omega^2 = log(1 + 0.39^2) = 0.142

    # =============== Residual error on plasma concentration ==================
    propSd <- fixed(0.28); label("Proportional residual error on plasma olanzapine concentration (fraction)")  # Part I Table 2: RUV proportional olanzapine = 0.28 (95% CI 0.26-0.31)

    # =============== Placebo Weibull model (per subscale) ====================
    # PANSS positive subscale placebo
    basl_pos <- 22.4; label("Baseline PANSS positive subscale (units)")                                       # Part II Table 2 (olanzapine): BASL PANSS positive = 22.4 (95% CI 22.4-22.6)
    pmax_pos <- 0.094; label("Maximum placebo improvement fraction for PANSS positive (unitless)")            # Part II Table 1: Pmax positive = 0.094
    td_pos   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS positive (days)")               # Part II Table 1: TD positive = 15 days
    pow_pos  <- 1.26; label("Weibull shape parameter for PANSS positive placebo (unitless)")                  # Part II Table 1: POW positive = 1.26

    # PANSS negative subscale placebo
    basl_neg <- 23.8; label("Baseline PANSS negative subscale (units)")                                       # Part II Table 2 (olanzapine): BASL PANSS negative = 23.8 (95% CI 23.5-24.2)
    pmax_neg <- 0.052; label("Maximum placebo improvement fraction for PANSS negative (unitless)")            # Part II Table 1: Pmax negative = 0.052
    td_neg   <- 19;   label("Time to reach 63.2% of max placebo effect, PANSS negative (days)")               # Part II Table 1: TD negative = 19 days
    pow_neg  <- 1.39; label("Weibull shape parameter for PANSS negative placebo (unitless)")                  # Part II Table 1: POW negative = 1.39

    # PANSS general subscale placebo
    basl_gen <- 44;   label("Baseline PANSS general subscale (units)")                                        # Part II Table 2 (olanzapine): BASL PANSS general = 44 (95% CI 43.7-44.2)
    pmax_gen <- 0.048; label("Maximum placebo improvement fraction for PANSS general (unitless)")             # Part II Table 1: Pmax general = 0.048
    td_gen   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS general (days)")                # Part II Table 1: TD general = 15 days
    pow_gen  <- 1.32; label("Weibull shape parameter for PANSS general placebo (unitless)")                   # Part II Table 1: POW general = 1.32

    # =============== Drug effect (olanzapine-specific, per subscale) ========
    emax_pos <- 0.43; label("Maximum olanzapine drug-effect fraction on PANSS positive (unitless)")           # Part II Table 2: Emax PANSS positive = 0.43 (95% CI 0.30-0.78)
    ec50_pos <- 12.4; label("Olanzapine Css EC50 for PANSS positive (ng/mL)")                                 # Part II Table 2: EC50 PANSS positive = 12.4 ng/mL (95% CI 1.2-43.3)
    kt_pos   <- 0.048; label("Onset rate constant of olanzapine effect on PANSS positive (1/day)")            # Part II Table 2: KT PANSS positive = 0.048 1/day (common across SGAs)

    emax_neg <- 0.33; label("Maximum olanzapine drug-effect fraction on PANSS negative (unitless)")           # Part II Table 2: Emax PANSS negative = 0.33 (95% CI 0.22-0.52); highest across SGAs
    ec50_neg <- 10.1; label("Olanzapine Css EC50 for PANSS negative (ng/mL)")                                 # Part II Table 2: EC50 PANSS negative = 10.1 ng/mL (95% CI 2.2-25)
    kt_neg   <- 0.028; label("Onset rate constant of olanzapine effect on PANSS negative (1/day)")            # Part II Table 2: KT PANSS negative = 0.028 1/day (95% CI 0.016-0.05). Estimated separately per drug.

    emax_gen <- 0.34; label("Maximum olanzapine drug-effect fraction on PANSS general (unitless)")            # Part II Table 2: Emax PANSS general = 0.34 (95% CI 0.25-0.55)
    ec50_gen <- 12.2; label("Olanzapine Css EC50 for PANSS general (ng/mL)")                                  # Part II Table 2: EC50 PANSS general = 12.2 ng/mL (95% CI 2.7-37)
    kt_gen   <- 0.035; label("Onset rate constant of olanzapine effect on PANSS general (1/day)")             # Part II Table 2: KT PANSS general = 0.035 1/day (common across SGAs)

    # =============== IIV on PD (per subscale) ================================
    etabasl_pos ~ 0.0515   # Part II Table 1: IIV BASL positive 23% CV; omega^2 = log(1 + 0.23^2) = 0.0515
    etabasl_neg ~ 0.0473   # Part II Table 1: IIV BASL negative 22% CV
    etabasl_gen ~ 0.0319   # Part II Table 1: IIV BASL general  18% CV

    etapmax_pos ~ 0.0576   # Part II Table 1: IIV Pmax positive SD = 0.24
    etapmax_neg ~ 0.0289   # Part II Table 1: IIV Pmax negative SD = 0.17
    etapmax_gen ~ 0.0441   # Part II Table 1: IIV Pmax general  SD = 0.21

    etaemax_pos ~ 0.0625   # Part II Table 2 olanzapine: IIV Emax positive SD = 0.25
    etaemax_neg ~ 0.1600   # Part II Table 2 olanzapine: IIV Emax negative SD = 0.40
    etaemax_gen ~ 0.0484   # Part II Table 2 olanzapine: IIV Emax general  SD = 0.22

    etaec50_pos ~ fixed(0.2231)  # Part II Table 2 footnote: IIV EC50 positive and general fixed at 50% CV
    etaec50_gen ~ fixed(0.2231)
    etaec50_neg ~ 1.0948         # Part II Table 2 olanzapine: IIV EC50 negative = 141% CV; omega^2 = log(1 + 1.41^2) = 1.0948

    # =============== Residual error on PANSS subscales =======================
    # Per-drug RUV from Part II Table 2 (joint PKPD model); Part II Table 1
    # reports a placebo-arm-only RUV that is superseded here.
    addSd_PANSS_pos <- 2.0; label("Additive residual error on PANSS positive (PANSS units)")                  # Part II Table 2 (olanzapine): RUV PANSS positive SD = 2.0 (95% CI 1.9-2.1)
    addSd_PANSS_neg <- 2.2; label("Additive residual error on PANSS negative (PANSS units)")                  # Part II Table 2 (olanzapine): RUV PANSS negative SD = 2.2 (95% CI 2.2-2.3)
    addSd_PANSS_gen <- 3.6; label("Additive residual error on PANSS general (PANSS units)")                   # Part II Table 2 (olanzapine): RUV PANSS general  SD = 3.6 (95% CI 3.5-3.7)
  })

  model({
    t_days <- t / 24

    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- 1000 * central / vc

    basl_pos_i <- basl_pos * exp(etabasl_pos)
    basl_neg_i <- basl_neg * exp(etabasl_neg)
    basl_gen_i <- basl_gen * exp(etabasl_gen)

    pmax_pos_i <- pmax_pos + etapmax_pos
    pmax_neg_i <- pmax_neg + etapmax_neg
    pmax_gen_i <- pmax_gen + etapmax_gen

    emax_pos_i <- emax_pos + etaemax_pos
    emax_neg_i <- emax_neg + etaemax_neg
    emax_gen_i <- emax_gen + etaemax_gen

    ec50_pos_i <- ec50_pos * exp(etaec50_pos)
    ec50_neg_i <- ec50_neg * exp(etaec50_neg)
    ec50_gen_i <- ec50_gen * exp(etaec50_gen)

    Pplacebo_pos <- pmax_pos_i * (1 - exp(-(t_days / td_pos)^pow_pos))
    Pplacebo_neg <- pmax_neg_i * (1 - exp(-(t_days / td_neg)^pow_neg))
    Pplacebo_gen <- pmax_gen_i * (1 - exp(-(t_days / td_gen)^pow_gen))

    Drug_pos <- emax_pos_i * Cc / (ec50_pos_i + Cc) * (1 - exp(-kt_pos * t_days))
    Drug_neg <- emax_neg_i * Cc / (ec50_neg_i + Cc) * (1 - exp(-kt_neg * t_days))
    Drug_gen <- emax_gen_i * Cc / (ec50_gen_i + Cc) * (1 - exp(-kt_gen * t_days))

    PANSS_pos <- basl_pos_i * (1 - Pplacebo_pos) * (1 - Drug_pos)
    PANSS_neg <- basl_neg_i * (1 - Pplacebo_neg) * (1 - Drug_neg)
    PANSS_gen <- basl_gen_i * (1 - Pplacebo_gen) * (1 - Drug_gen)

    Cc        ~ prop(propSd)
    PANSS_pos ~ add(addSd_PANSS_pos)
    PANSS_neg ~ add(addSd_PANSS_neg)
    PANSS_gen ~ add(addSd_PANSS_gen)
  })
}
