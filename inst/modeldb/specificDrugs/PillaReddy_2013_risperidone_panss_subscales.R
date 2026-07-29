# Population PK/PD model for risperidone active moiety (risperidone +
# 9-hydroxy-risperidone) against the three PANSS subscales (positive,
# negative, general) in adults with schizophrenia (Pilla Reddy 2013 Part II;
# Schizophrenia Research 146:153-161; doi:10.1016/j.schres.2013.02.010).

PillaReddy_2013_risperidone_panss_subscales <- function() {
  description <- paste(
    "Population PK/PD model for risperidone against the three PANSS",
    "subscales (positive, negative, general) in adults with schizophrenia",
    "from Pilla Reddy 2013 Part II. The driving exposure variable is the",
    "active moiety (parent risperidone + the equipotent metabolite",
    "9-hydroxy-risperidone), following the Part I (PMID 23473810) and",
    "Vermeulen 2007 (Eur J Clin Pharmacol 63:1063-1077) methodology. The",
    "PK sub-model is a simplified one-compartment representation of the",
    "active moiety: the published Part I model for risperidone is a",
    "two-compartment parent-plus-metabolite system with a lag-time and a",
    "consecutive zero-then-first-order absorption process, and clearance",
    "stratified by CYP2D6 phenotype (poor / medium / fast metabolizers).",
    "This nlmixr2lib representation uses a single CL/F that maps a typical",
    "oral risperidone dose to the active-moiety steady-state",
    "concentration Css consumed by the PD model: CL_AM/F = 6.3 L/h",
    "(derived from Part II Table 3 effective-dose / effective-Css pair",
    "0.8 mg/day / 5.3 ng/mL at 30% PANSS reduction, CL/F = Dose / (Css *",
    "tau)), Vc = 144 L (parent central volume from Part I Table 2;",
    "approximate active-moiety value), and Ka = 2.37 1/h (parent first-",
    "order rate from Part I Table 2). The zero-order absorption duration",
    "DUR = 0.47 h and lag time ALAG1 = 0.16 h are omitted in this Part-II",
    "simplification (Css is approximately invariant to short-window",
    "absorption details at steady state). The PD sub-model has three",
    "outputs that share the Weibull placebo time-course form Pplacebo = Pmax",
    "* (1 - exp(-(t/TD)^POW)) but each subscale carries its own placebo",
    "Pmax, TD, POW (Part II Table 1) and risperidone's own Emax / EC50 /",
    "KT triplet per subscale (Part II Table 2). The KT for risperidone",
    "PANSS positive and general subscales (0.048 / 0.035 1/day) was",
    "estimated as a common value across all atypical antipsychotic drugs",
    "(Part II Methods, 'A common model ... was developed for PANSS positive",
    "and general scales'); the KT for the negative subscale (0.16 1/day)",
    "was estimated separately per drug because the cross-drug pooled fit",
    "did not converge. The exponential time-to-event dropout sub-model",
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
    "(PMID 23473810; doi:10.1016/j.schres.2013.02.011), with active-moiety",
    "scaling rationale per Vermeulen et al. 2007 (Eur J Clin Pharmacol",
    "63:1063-1077; doi:10.1007/s00228-007-0358-5).",
    sep = " "
  )
  vignette <- "PillaReddy_2013_panss_subscales"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    CYP2D6_PHENO = list(
      description = "CYP2D6 metabolizer phenotype (poor / medium / fast)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Part I Table 2 reports separate risperidone parent CL/F values",
        "by CYP2D6 phenotype (poor 0.44 L/h, medium 2.81 L/h, fast 18.4",
        "L/h). Because the active moiety Css depends on the relative",
        "balance between parent and 9-OH-risperidone elimination, the",
        "active-moiety CL_AM/F observed in this Part II simplification",
        "(6.3 L/h, see ini comments) corresponds to the population-typical",
        "phenotype mix and not a single CYP2D6 stratum. CYP2D6 phenotype",
        "would shift Css from a fixed dose but is not modelled in the",
        "Part II model() body."
      )
    ),
    DIS = list(
      description = "Disease state at entry (acute vs chronic schizophrenia)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part II placebo-model covariate; not implemented (typical-individual reference simulation targets acute stratum)."
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
    age_range       = "Adults with schizophrenia (specific range not tabulated in Part II; pooled from 12 industry-sponsored Phase II / III trials 1989-2009).",
    weight_range    = "Adult schizophrenia population; specific demographics not tabulated in Part II (cross-reference Part I Table 1 for per-study summaries).",
    sex_female_pct  = NA_real_,
    disease_state   = paste(
      "Adults with schizophrenia recruited into acute and chronic-stable",
      "Phase II / III double-blind clinical trials. Risperidone arms in",
      "Part I Table 1: INT-2 (0.5, 2, 4, 6, 8 mg) and INT-3 (1, 3, 5, 8",
      "mg). Baseline PANSS positive typical value 22.4 (95% CI 22.4-22.6),",
      "PANSS negative 24.7 (24.3-25), PANSS general 44 (43.7-44.2) for the",
      "atypical-antipsychotic pool (Part II Table 2)."
    ),
    dose_range      = "Oral risperidone 0.5-8 mg/day, qd or bid (Part I Table 1, studies INT-2 and INT-3).",
    regions         = "Pooled across multinational schizophrenia trials 1989-2009 (Part II Methods).",
    notes           = paste(
      "The driving exposure variable in the PD model is the active moiety",
      "(parent risperidone + 9-hydroxy-risperidone); see description.",
      "The active-moiety Css value used for effective-dose calculation is",
      "5.3 ng/mL at 30% PANSS total reduction (Part II Table 3); the",
      "PANSS-subscale-specific EC50 values in ini() use the same active-",
      "moiety scale."
    ),
    dropout_model = paste(
      "Part II Table 4 reports an exponential time-to-event dropout",
      "sub-model with subscale-specific baseline hazards (BHAZ) and BETA",
      "coefficients. For risperidone the placebo BHAZ values are 0.00072",
      "(positive), 0.00219 (negative), 0.00038 (general) 1/day; the BETA",
      "coefficients are -0.111 (positive), -0.0443 (negative), -0.0672",
      "(general). Not encoded in this nlmixr2lib model body."
    )
  )

  ini({
    # =============== Pharmacokinetic structural parameters ===================
    # PK held fixed at the Part I-derived active-moiety values; see
    # description for the derivation of CL_AM/F.
    lka <- fixed(log(2.37));  label("Log first-order absorption rate constant (1/h)")              # Part I Table 2: Ka risperidone parent = 2.37 1/h (95% CI 1.65-3.08); used as approximate active-moiety appearance rate
    lcl <- fixed(log(6.3));   label("Log active-moiety apparent oral clearance CL_AM/F (L/h)")     # Part II Table 3 derivation: CL_AM/F = Dose / (Css * tau) = 0.8 mg/day / (5.3 ng/mL * 24 h) = 6.29 L/h (rounded to 6.3 L/h)
    lvc <- fixed(log(144));   label("Log active-moiety apparent central volume Vc/F (L)")          # Part I Table 2: Vc risperidone parent = 144 L (95% CI 122-165); approximate active-moiety value

    # =============== Inter-individual variability on PK ======================
    # Parent IIV from Part I Table 2 (169% CV CL, 54% CV Vc); applied here
    # as approximate active-moiety IIV. These are large CVs reflecting the
    # CYP2D6-phenotype mixing baked into the parent CL distribution.
    etalcl ~ fixed(1.350)  # Part I Table 2: IIV CL risperidone parent 169% CV; omega^2 = log(1 + 1.69^2) = 1.350
    etalvc ~ fixed(0.256)  # Part I Table 2: IIV Vc parent 54% CV; omega^2 = log(1 + 0.54^2) = 0.256

    # =============== Residual error on plasma concentration ==================
    propSd <- fixed(0.30); label("Proportional residual error on risperidone parent plasma concentration (fraction)")  # Part I Table 2: RUV proportional risperidone parent = 0.30 (95% CI 0.25-0.35)

    # =============== Placebo Weibull model (per subscale) ====================
    # Risperidone baseline PANSS values are the atypical-antipsychotic pool
    # values (Part II Table 2 lists 22.4 / 24.7 / 44.0 for risperidone).
    # Pmax, TD, POW are common across drugs per subscale (Part II Table 1).

    # PANSS positive subscale placebo
    basl_pos <- 22.4; label("Baseline PANSS positive subscale (units)")                                       # Part II Table 2 (risperidone): BASL PANSS positive = 22.4 (95% CI 22.4-22.6)
    pmax_pos <- 0.094; label("Maximum placebo improvement fraction for PANSS positive (unitless)")            # Part II Table 1: Pmax positive = 0.094 (95% CI 0.058-0.125)
    td_pos   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS positive (days)")               # Part II Table 1: TD positive = 15 days (95% CI 13-18)
    pow_pos  <- 1.26; label("Weibull shape parameter for PANSS positive placebo (unitless)")                  # Part II Table 1: POW positive = 1.26 (95% CI 1.06-1.53)

    # PANSS negative subscale placebo
    basl_neg <- 24.7; label("Baseline PANSS negative subscale (units)")                                       # Part II Table 2 (risperidone): BASL PANSS negative = 24.7 (95% CI 24.3-25)
    pmax_neg <- 0.052; label("Maximum placebo improvement fraction for PANSS negative (unitless)")            # Part II Table 1: Pmax negative = 0.052 (95% CI 0.032-0.076)
    td_neg   <- 19;   label("Time to reach 63.2% of max placebo effect, PANSS negative (days)")               # Part II Table 1: TD negative = 19 days (95% CI 14.8-36.8)
    pow_neg  <- 1.39; label("Weibull shape parameter for PANSS negative placebo (unitless)")                  # Part II Table 1: POW negative = 1.39 (95% CI 1.51-1.99)

    # PANSS general subscale placebo
    basl_gen <- 44;   label("Baseline PANSS general subscale (units)")                                        # Part II Table 2 (risperidone): BASL PANSS general = 44 (95% CI 43.7-44.2)
    pmax_gen <- 0.048; label("Maximum placebo improvement fraction for PANSS general (unitless)")             # Part II Table 1: Pmax general = 0.048 (95% CI 0.028-0.066)
    td_gen   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS general (days)")                # Part II Table 1: TD general = 15 days (95% CI 13-18)
    pow_gen  <- 1.32; label("Weibull shape parameter for PANSS general placebo (unitless)")                   # Part II Table 1: POW general = 1.32 (95% CI 1.13-1.54)

    # =============== Drug effect (risperidone-specific, per subscale) ========
    emax_pos <- 0.32; label("Maximum risperidone drug-effect fraction on PANSS positive (unitless)")          # Part II Table 2: Emax PANSS positive = 0.32 (95% CI 0.28-0.38)
    ec50_pos <- 9.4;  label("Risperidone active-moiety Css EC50 for PANSS positive (ng/mL)")                  # Part II Table 2: EC50 PANSS positive = 9.4 ng/mL (95% CI 3.4-21.7)
    kt_pos   <- 0.048; label("Onset rate constant of risperidone effect on PANSS positive (1/day)")           # Part II Table 2: KT PANSS positive = 0.048 1/day (95% CI 0.039-0.057). Common value across SGAs.

    emax_neg <- 0.14; label("Maximum risperidone drug-effect fraction on PANSS negative (unitless)")          # Part II Table 2: Emax PANSS negative = 0.14 (95% CI 0.09-0.15)
    ec50_neg <- 18.5; label("Risperidone active-moiety Css EC50 for PANSS negative (ng/mL)")                  # Part II Table 2: EC50 PANSS negative = 18.5 ng/mL (95% CI 3.3-34)
    kt_neg   <- 0.16; label("Onset rate constant of risperidone effect on PANSS negative (1/day)")            # Part II Table 2: KT PANSS negative = 0.16 1/day (95% CI 0.12-0.21). Estimated separately per drug.

    emax_gen <- 0.19; label("Maximum risperidone drug-effect fraction on PANSS general (unitless)")           # Part II Table 2: Emax PANSS general = 0.19 (95% CI 0.16-0.23)
    ec50_gen <- 3.97; label("Risperidone active-moiety Css EC50 for PANSS general (ng/mL)")                   # Part II Table 2: EC50 PANSS general = 3.97 ng/mL (95% CI 0.49-11.1)
    kt_gen   <- 0.035; label("Onset rate constant of risperidone effect on PANSS general (1/day)")            # Part II Table 2: KT PANSS general = 0.035 1/day (95% CI 0.021-0.045). Common value across SGAs.

    # =============== IIV on PD (per subscale) ================================
    etabasl_pos ~ 0.0515   # Part II Table 1: IIV BASL positive 23% CV; omega^2 = log(1 + 0.23^2) = 0.0515
    etabasl_neg ~ 0.0473   # Part II Table 1: IIV BASL negative 22% CV
    etabasl_gen ~ 0.0319   # Part II Table 1: IIV BASL general  18% CV

    etapmax_pos ~ 0.0576   # Part II Table 1: IIV Pmax positive SD = 0.24 (additive normal); variance 0.24^2 = 0.0576
    etapmax_neg ~ 0.0289   # Part II Table 1: IIV Pmax negative SD = 0.17
    etapmax_gen ~ 0.0441   # Part II Table 1: IIV Pmax general  SD = 0.21

    etaemax_pos ~ 0.0625   # Part II Table 2 risperidone: IIV Emax positive SD = 0.25; variance 0.0625
    etaemax_neg ~ 0.0576   # Part II Table 2 risperidone: IIV Emax negative SD = 0.24
    etaemax_gen ~ 0.0484   # Part II Table 2 risperidone: IIV Emax general  SD = 0.22

    etaec50_pos ~ fixed(0.2231)  # Part II Table 2 footnote: IIV EC50 positive and general fixed at 50% CV
    etaec50_gen ~ fixed(0.2231)
    etaec50_neg ~ 2.367          # Part II Table 2 risperidone: IIV EC50 negative = 311% CV; omega^2 = log(1 + 3.11^2) = 2.367

    # =============== Residual error on PANSS subscales =======================
    # Per-drug RUV from Part II Table 2 (joint PKPD model); Part II Table 1
    # reports a placebo-arm-only RUV that is superseded here.
    addSd_PANSS_pos <- 2.0; label("Additive residual error on PANSS positive (PANSS units)")                  # Part II Table 2 (risperidone): RUV PANSS positive SD = 2.0 (95% CI 1.9-2.1)
    addSd_PANSS_neg <- 2.3; label("Additive residual error on PANSS negative (PANSS units)")                  # Part II Table 2 (risperidone): RUV PANSS negative SD = 2.3 (95% CI 2.2-2.3)
    addSd_PANSS_gen <- 3.6; label("Additive residual error on PANSS general (PANSS units)")                   # Part II Table 2 (risperidone): RUV PANSS general  SD = 3.6 (95% CI 3.5-3.7)
  })

  model({
    t_days <- t / 24

    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma active-moiety risperidone concentration in ng/mL.
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
