# Population PK/PD model for paliperidone ER against the three PANSS subscales
# (positive, negative, general) in adults with schizophrenia (Pilla Reddy 2013
# Part II; Schizophrenia Research 146:153-161; doi:10.1016/j.schres.2013.02.010).

PillaReddy_2013_paliperidone_panss_subscales <- function() {
  description <- paste(
    "Population PK/PD model for paliperidone extended release against the",
    "three PANSS subscales (positive, negative, general) in adults with",
    "schizophrenia from Pilla Reddy 2013 Part II. The PK sub-model is the",
    "one-compartment paliperidone structural model from Part I (PMID",
    "23473810) Table 2 with the sequential zero-order plus first-order",
    "absorption simplified to first-order only (the ER absorption profile",
    "is approximately steady at steady state because the OROS-extended-",
    "release tablet is dosed once daily): first-order absorption ka = 0.57",
    "1/h, apparent oral clearance CL/F = 14.1 L/h, apparent central volume",
    "of distribution Vc/F = 475 L. The PD sub-model has three outputs that",
    "share the Weibull placebo time-course form Pplacebo = Pmax * (1 -",
    "exp(-(t/TD)^POW)) but each subscale carries its own placebo Pmax, TD,",
    "POW (Part II Table 1) and paliperidone's own Emax / EC50 / KT triplet",
    "per subscale (Part II Table 2). The KT for paliperidone PANSS positive",
    "and general (0.048 and 0.035 1/day) is the common-across-atypical-",
    "antipsychotic value; the KT for the negative subscale (0.13 1/day) was",
    "estimated separately per drug. The lag time ALAG1 = 0.67 h and the",
    "zero-order absorption duration DUR = 23.6 h reported in Part I Table 2",
    "for paliperidone ER are not encoded here because the Css used by the",
    "PD model is approximately invariant to within-dose absorption details",
    "at steady state. The exponential time-to-event dropout sub-model from",
    "Part II Table 4 is documented in population$dropout_model but not",
    "encoded in the model body.",
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
    LBM = list(
      description = "Lean body mass",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Part I Table 2 reports a power-law allometric effect of LBM on",
        "paliperidone CL/F: exponent 0.82 (95% CI 0.44-1.25). Not retained",
        "in the Part II PD model() body (the typical-individual reference",
        "simulation uses the reported typical CL/F = 14.1 L/h)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part I Table 2 reports a 14% lower paliperidone CL/F in females vs males ('CL (gender) -0.14 (-0.25 to -0.024)'); not retained in the Part II PD model() body."
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
      notes       = "Part II placebo-model covariate; not implemented."
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
      "Adults with schizophrenia recruited into acute Phase III double-",
      "blind clinical trials. Paliperidone ER arms in Part I Table 1:",
      "SCH-303 (6, 9, 12 mg qd), SCH-304 (6, 12 mg qd), SCH-305 (3, 9, 15",
      "mg qd). Baseline PANSS positive typical value 22.4, PANSS negative",
      "23.8 (95% CI 23.5-24.1), PANSS general 44 for the atypical-",
      "antipsychotic pool (Part II Table 2)."
    ),
    dose_range      = "Oral paliperidone ER 3-15 mg/day, qd (Part I Table 1).",
    regions         = "Pooled across multinational schizophrenia trials 1989-2009.",
    notes           = paste(
      "Paliperidone ER is the OROS extended-release oral formulation; the",
      "paper's intramuscular paliperidone palmitate model was held out for",
      "external validation in Part I and is not addressed in this Part II",
      "model file."
    ),
    dropout_model = paste(
      "Part II Table 4: paliperidone placebo BHAZ values are 0.00119",
      "(positive), 0.00409 (negative), 0.00067 (general) 1/day. Not",
      "encoded in this model body."
    )
  )

  ini({
    # =============== Pharmacokinetic structural parameters ===================
    lka <- fixed(log(0.57));  label("Log first-order absorption rate constant (1/h)")              # Part I Table 2: Ka paliperidone ER = 0.57 1/h (NONMEM final estimate; *flag indicates separately fit)
    lcl <- fixed(log(14.1));  label("Log apparent oral clearance CL/F (L/h)")                       # Part I Table 2: CL/F = 14.1 L/h (95% CI 13.2-14.9)
    lvc <- fixed(log(475));   label("Log apparent central volume of distribution Vc/F (L)")         # Part I Table 2: Vc/F = 475 L (95% CI 325-678)

    # =============== Inter-individual variability on PK ======================
    etalcl ~ fixed(0.231)  # Part I Table 2: IIV CL paliperidone 51% CV; omega^2 = log(1 + 0.51^2) = 0.231

    # =============== Residual error on plasma concentration ==================
    propSd <- fixed(0.38); label("Proportional residual error on plasma paliperidone concentration (fraction)")  # Part I Table 2: RUV proportional paliperidone = 0.38 (95% CI 0.35-0.40)

    # =============== Placebo Weibull model (per subscale) ====================
    # PANSS positive subscale placebo
    basl_pos <- 22.4; label("Baseline PANSS positive subscale (units)")                                       # Part II Table 2 (paliperidone): BASL PANSS positive = 22.4 (95% CI 22.4-22.6)
    pmax_pos <- 0.094; label("Maximum placebo improvement fraction for PANSS positive (unitless)")            # Part II Table 1: Pmax positive = 0.094
    td_pos   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS positive (days)")               # Part II Table 1: TD positive = 15 days
    pow_pos  <- 1.26; label("Weibull shape parameter for PANSS positive placebo (unitless)")                  # Part II Table 1: POW positive = 1.26

    # PANSS negative subscale placebo
    basl_neg <- 23.8; label("Baseline PANSS negative subscale (units)")                                       # Part II Table 2 (paliperidone): BASL PANSS negative = 23.8 (95% CI 23.5-24.1)
    pmax_neg <- 0.052; label("Maximum placebo improvement fraction for PANSS negative (unitless)")            # Part II Table 1: Pmax negative = 0.052
    td_neg   <- 19;   label("Time to reach 63.2% of max placebo effect, PANSS negative (days)")               # Part II Table 1: TD negative = 19 days
    pow_neg  <- 1.39; label("Weibull shape parameter for PANSS negative placebo (unitless)")                  # Part II Table 1: POW negative = 1.39

    # PANSS general subscale placebo
    basl_gen <- 44;   label("Baseline PANSS general subscale (units)")                                        # Part II Table 2 (paliperidone): BASL PANSS general = 44 (95% CI 43.7-44.2)
    pmax_gen <- 0.048; label("Maximum placebo improvement fraction for PANSS general (unitless)")             # Part II Table 1: Pmax general = 0.048
    td_gen   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS general (days)")                # Part II Table 1: TD general = 15 days
    pow_gen  <- 1.32; label("Weibull shape parameter for PANSS general placebo (unitless)")                   # Part II Table 1: POW general = 1.32

    # =============== Drug effect (paliperidone-specific, per subscale) =======
    emax_pos <- 0.33; label("Maximum paliperidone drug-effect fraction on PANSS positive (unitless)")         # Part II Table 2: Emax PANSS positive = 0.33 (95% CI 0.28-0.39)
    ec50_pos <- 5.84; label("Paliperidone Css EC50 for PANSS positive (ng/mL)")                               # Part II Table 2: EC50 PANSS positive = 5.84 ng/mL (95% CI 1.6-12.3)
    kt_pos   <- 0.048; label("Onset rate constant of paliperidone effect on PANSS positive (1/day)")          # Part II Table 2: KT PANSS positive = 0.048 1/day (common across SGAs)

    emax_neg <- 0.15; label("Maximum paliperidone drug-effect fraction on PANSS negative (unitless)")         # Part II Table 2: Emax PANSS negative = 0.15 (95% CI 0.12-0.20)
    ec50_neg <- 17.3; label("Paliperidone Css EC50 for PANSS negative (ng/mL)")                               # Part II Table 2: EC50 PANSS negative = 17.3 ng/mL (95% CI 11.7-57)
    kt_neg   <- 0.13; label("Onset rate constant of paliperidone effect on PANSS negative (1/day)")           # Part II Table 2: KT PANSS negative = 0.13 1/day (95% CI 0.10-0.22)

    emax_gen <- 0.24; label("Maximum paliperidone drug-effect fraction on PANSS general (unitless)")          # Part II Table 2: Emax PANSS general = 0.24 (95% CI 0.20-0.30)
    ec50_gen <- 3.07; label("Paliperidone Css EC50 for PANSS general (ng/mL)")                                # Part II Table 2: EC50 PANSS general = 3.07 ng/mL (95% CI 0.23-10)
    kt_gen   <- 0.035; label("Onset rate constant of paliperidone effect on PANSS general (1/day)")           # Part II Table 2: KT PANSS general = 0.035 1/day (common across SGAs)

    # =============== IIV on PD (per subscale) ================================
    etabasl_pos ~ 0.0515   # Part II Table 1: IIV BASL positive 23% CV; omega^2 = log(1 + 0.23^2) = 0.0515
    etabasl_neg ~ 0.0473   # Part II Table 1: IIV BASL negative 22% CV
    etabasl_gen ~ 0.0319   # Part II Table 1: IIV BASL general  18% CV

    etapmax_pos ~ 0.0576   # Part II Table 1: IIV Pmax positive SD = 0.24
    etapmax_neg ~ 0.0289   # Part II Table 1: IIV Pmax negative SD = 0.17
    etapmax_gen ~ 0.0441   # Part II Table 1: IIV Pmax general  SD = 0.21

    etaemax_pos ~ 0.0625   # Part II Table 2 paliperidone: IIV Emax positive SD = 0.25
    etaemax_neg ~ 0.0729   # Part II Table 2 paliperidone: IIV Emax negative SD = 0.27
    etaemax_gen ~ 0.0484   # Part II Table 2 paliperidone: IIV Emax general  SD = 0.22

    etaec50_pos ~ fixed(0.2231)  # Part II Table 2 footnote: IIV EC50 positive and general fixed at 50% CV
    etaec50_gen ~ fixed(0.2231)
    etaec50_neg ~ 1.8095         # Part II Table 2 paliperidone: IIV EC50 negative = 226% CV; omega^2 = log(1 + 2.26^2) = 1.8095

    # =============== Residual error on PANSS subscales =======================
    # Per-drug RUV from Part II Table 2 (joint PKPD model); Part II Table 1
    # reports a placebo-arm-only RUV that is superseded here.
    addSd_PANSS_pos <- 2.0; label("Additive residual error on PANSS positive (PANSS units)")                  # Part II Table 2 (paliperidone): RUV PANSS positive SD = 2.0
    addSd_PANSS_neg <- 2.0; label("Additive residual error on PANSS negative (PANSS units)")                  # Part II Table 2 (paliperidone): RUV PANSS negative SD = 2.0 (95% CI 1.9-2.2)
    addSd_PANSS_gen <- 3.6; label("Additive residual error on PANSS general (PANSS units)")                   # Part II Table 2 (paliperidone): RUV PANSS general  SD = 3.6
  })

  model({
    t_days <- t / 24

    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
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
