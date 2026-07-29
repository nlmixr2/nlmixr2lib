# Population PK/PD model for haloperidol against the three PANSS subscales
# (positive, negative, general) in adults with schizophrenia (Pilla Reddy 2013
# Part II; Schizophrenia Research 146:153-161; doi:10.1016/j.schres.2013.02.010).

PillaReddy_2013_haloperidol_panss_subscales <- function() {
  description <- paste(
    "Population pharmacokinetic / pharmacodynamic (PK/PD) model for",
    "haloperidol against the three PANSS subscales (positive, negative,",
    "general) in adults with schizophrenia from Pilla Reddy 2013 Part II",
    "(pooled dataset of 1338 placebo-treated patients with subscale data",
    "available for n=741 of them; 12 industry-sponsored Phase II/III",
    "trials between 1989 and 2009 plus one open-label study). The PK",
    "sub-model uses the haloperidol structural model inherited from Part I",
    "(PMID 23473810): a simplified one-compartment representation",
    "parameterised with the typical apparent oral clearance CL/F = 88 L/h,",
    "apparent central volume Vc/F = 669 L, and first-order absorption rate",
    "ka = 0.23 1/h reported in Part I Table 2 (Pilla Reddy 2012a originally",
    "described the haloperidol PK as two-compartment with Q/F = 233 L/h",
    "and Vp/F = 2500 L; the one-compartment simplification preserves the",
    "steady-state average concentration Css = Dose / (CL * tau) used as the",
    "PD driver and is documented as a deliberate Part II simplification in",
    "the vignette Assumptions and deviations). The PD sub-model has three",
    "outputs that share the Weibull placebo time-course form Pplacebo = Pmax",
    "* (1 - exp(-(t/TD)^POW)) but each subscale carries its own placebo",
    "baseline, Pmax, TD, and POW (Pilla Reddy 2013 Part II Table 1) and",
    "haloperidol's own Emax / EC50 / KT triplet per subscale (Part II Table",
    "2). The combined PANSS subscale prediction is BASL * (1 - Pplacebo) *",
    "(1 - Drug), where Drug = Emax * Cc / (EC50 + Cc) * (1 - exp(-KT * t))",
    "with the KT delay capturing the onset time to the maximum drug effect.",
    "Concentration drives the PD with a steady-state assumption: the paper",
    "uses Css from the PK model and feeds it into the Emax equation; in the",
    "rxode2 implementation Cc is the time-varying plasma concentration",
    "derived from the one-compartment ODE, which is approximately constant",
    "at steady state for the once- and twice-daily regimens studied here.",
    "An exponential time-to-event dropout sub-model with subscale- and",
    "drug-specific baseline hazards (Part II Table 4) is reported in the",
    "paper but is not encoded in this nlmixr2lib model body; the dropout",
    "parameters are documented in population$dropout_model and discussed in",
    "the vignette.",
    sep = " "
  )
  reference <- paste(
    "Pilla Reddy V, Kozielska M, Suleiman AA, Johnson M, Vermeulen A, Liu J,",
    "de Greef R, Groothuis GMM, Danhof M, Proost JH (2013).",
    "Pharmacokinetic-pharmacodynamic modelling of antipsychotic drugs in",
    "patients with schizophrenia: Part II: The use of subscales of the PANSS",
    "score. Schizophrenia Research 146(1-3):153-161.",
    "doi:10.1016/j.schres.2013.02.010. PK structure carried over from Part I",
    "(PMID 23473810; doi:10.1016/j.schres.2013.02.011).",
    sep = " "
  )
  vignette <- "PillaReddy_2013_panss_subscales"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    DIS = list(
      description = "Disease state at entry (acute vs chronic schizophrenia)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Reported by Part II as a covariate of baseline PANSS positive",
        "(coefficient -0.15) and of Pmax positive (-1.15) on the placebo",
        "Weibull model; not retained in the PANSS subscale model() because",
        "the typical-individual reference simulation targets the acute",
        "stratum (reference category). See vignette Assumptions and",
        "deviations for the encoding the paper uses."
      )
    ),
    USA = list(
      description = "Study geographic origin (USA vs non-USA)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Reported by Part II as a covariate of Pmax negative (+1.34) and",
        "of residual-error magnitude for all three subscales (+0.29 to",
        "+0.64 on log-RUV); not implemented in model() because the",
        "typical-individual reference simulation targets the USA-study",
        "stratum (reference category)."
      )
    ),
    REG = list(
      description = "Dosing regimen (qd vs bid)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported by Part II as a covariate of residual error on PANSS negative (-0.20 on log-RUV); not implemented."
    ),
    DUR = list(
      description = "Study duration (short vs long)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported by Part II as a covariate of residual error on PANSS positive (-0.28 on log-RUV); not implemented."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 741L,
    n_studies       = 12L,
    age_range       = "Adults with schizophrenia (specific range not tabulated in Part II; pooled from 12 industry-sponsored Phase II / III studies 1989-2009 and one open-label haloperidol study).",
    weight_range    = "Adult schizophrenia population; specific demographics not tabulated in Part II (cross-reference Part I Table 1 for per-study summaries).",
    sex_female_pct  = NA_real_,
    disease_state   = paste(
      "Adults with schizophrenia recruited into acute and chronic-stable",
      "Phase II / III double-blind clinical trials and one open-label",
      "haloperidol study (Part II Methods). Baseline PANSS positive typical",
      "value 23.3 (95% CI 22.9-23.6), PANSS negative 24.1 (23.7-24.4),",
      "PANSS general 45.1 (44.6-45.6) for the haloperidol arm (Part II",
      "Table 2)."
    ),
    dose_range      = paste(
      "Oral haloperidol 2.5-40 mg/day (Part I Table 1 LMU open-label study);",
      "5-10 mg/day Phase III arms (Part I Table 1, studies INT-2, INT-3,",
      "128-115)."
    ),
    regions         = "Pooled across multinational schizophrenia trials 1989-2009 (Part II Methods); USA vs non-USA reported as a placebo-effect covariate.",
    notes           = paste(
      "The n_subjects field (741) is the placebo-arm count of Part II Methods",
      "for whom PANSS subscale data were available; the active-drug cohort",
      "spans haloperidol, risperidone, olanzapine, ziprasidone, and",
      "paliperidone arms. The haloperidol-specific PD parameter set used here",
      "is from Part II Table 2."
    ),
    dropout_model = paste(
      "Part II Table 4 reports an exponential time-to-event dropout",
      "sub-model with subscale-specific baseline hazards (BHAZ) and BETA",
      "coefficients describing the dependence of dropout hazard on each",
      "PANSS subscale: h(t) = BHAZ * exp(BETA * PANSS_subscale). For",
      "haloperidol the placebo BHAZ values are 0.00119 (positive),",
      "0.00438 (negative), 0.00057 (general) 1/day, and the BETA",
      "coefficients are -0.111 (positive), -0.0443 (negative), -0.0672",
      "(general). The Part II 'positive + negative + general' joint",
      "dropout model is the preferred final model (lowest OFV). This",
      "nlmixr2lib model body simulates only the PANSS time course; the",
      "dropout hazards are documented here for reference."
    )
  )

  ini({
    # =============== Pharmacokinetic structural parameters ===================
    # PK parameters inherited from Part I (PMID 23473810) Table 2 for
    # haloperidol; held fixed in this Part II PD-focused model because the
    # Part II analysis did not re-fit PK -- it consumed Part I's predicted
    # Css as a covariate driver of the Emax equation.
    lka <- fixed(log(0.23)); label("Log first-order absorption rate constant (1/h)")              # Part I Table 2 (Pilla Reddy et al., 2012a haloperidol PK): Ka = 0.23 1/h (95% CI 0.056-0.39)
    lcl <- fixed(log(88));   label("Log apparent oral clearance CL/F (L/h)")                       # Part I Table 2: CL/F = 88 L/h (77-101)
    lvc <- fixed(log(669));  label("Log apparent central volume of distribution Vc/F (L)")         # Part I Table 2: Vc/F = 669 L (91-1143). Part I used a 2-cmt model with Q/F = 233 L/h and Vp/F = 2500 L; this Part II implementation collapses to 1-cmt for Css-based PD simulation (vignette Assumptions and deviations).

    # =============== Inter-individual variability on PK ======================
    # Part I Table 2 reports IIV as %CV; convert via omega^2 = log(1 + CV^2).
    # IIV on Ka was not estimated for haloperidol (no value in Part I Table
    # 2 row 'IIV Ka').
    etalcl ~ fixed(0.177)  # Part I Table 2: IIV CL 44% CV; omega^2 = log(1 + 0.44^2) = 0.177
    etalvc ~ fixed(0.853)  # Part I Table 2: IIV Vc 116% CV; omega^2 = log(1 + 1.16^2) = 0.853

    # =============== Residual error on plasma concentration ==================
    propSd <- fixed(0.44); label("Proportional residual error on plasma haloperidol concentration (fraction)")  # Part I Table 2: RUV proportional = 0.44 (95% CI 0.38-0.50)

    # =============== Placebo Weibull model (per subscale) ====================
    # All four parameters (BASL, Pmax, TD, POW) reported in Part II Table 1
    # for the typical reference patient (acute disease, USA, bid, short
    # study). Haloperidol-specific BASL values are taken from Part II Table 2
    # rather than Table 1 -- they slightly differ from the SGA pool (Table 2
    # footnote indicates BASL is estimated per-drug per-subscale).

    # PANSS positive subscale placebo
    basl_pos <- 23.3; label("Baseline PANSS positive subscale (units)")                                       # Part II Table 2 (haloperidol): BASL PANSS positive = 23.3 (95% CI 22.9-23.6)
    pmax_pos <- 0.094; label("Maximum placebo improvement fraction for PANSS positive (unitless)")            # Part II Table 1: Pmax positive = 0.094 (95% CI 0.058-0.125)
    td_pos   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS positive (days)")               # Part II Table 1: TD positive = 15 days (95% CI 13-18)
    pow_pos  <- 1.26; label("Weibull shape parameter for PANSS positive placebo (unitless)")                  # Part II Table 1: POW positive = 1.26 (95% CI 1.06-1.53)

    # PANSS negative subscale placebo
    basl_neg <- 24.1; label("Baseline PANSS negative subscale (units)")                                       # Part II Table 2 (haloperidol): BASL PANSS negative = 24.1 (95% CI 23.7-24.4)
    pmax_neg <- 0.052; label("Maximum placebo improvement fraction for PANSS negative (unitless)")            # Part II Table 1: Pmax negative = 0.052 (95% CI 0.032-0.076)
    td_neg   <- 19;   label("Time to reach 63.2% of max placebo effect, PANSS negative (days)")               # Part II Table 1: TD negative = 19 days (95% CI 14.8-36.8)
    pow_neg  <- 1.39; label("Weibull shape parameter for PANSS negative placebo (unitless)")                  # Part II Table 1: POW negative = 1.39 (95% CI 1.51-1.99)

    # PANSS general subscale placebo
    basl_gen <- 45.1; label("Baseline PANSS general subscale (units)")                                        # Part II Table 2 (haloperidol): BASL PANSS general = 45.1 (95% CI 44.6-45.6)
    pmax_gen <- 0.048; label("Maximum placebo improvement fraction for PANSS general (unitless)")             # Part II Table 1: Pmax general = 0.048 (95% CI 0.028-0.066)
    td_gen   <- 15;   label("Time to reach 63.2% of max placebo effect, PANSS general (days)")                # Part II Table 1: TD general = 15 days (95% CI 13-18)
    pow_gen  <- 1.32; label("Weibull shape parameter for PANSS general placebo (unitless)")                   # Part II Table 1: POW general = 1.32 (95% CI 1.13-1.54)

    # =============== Drug effect (haloperidol-specific, per subscale) ========
    # Emax (maximum drug-effect fraction), EC50 (steady-state concentration
    # giving half Emax), KT (onset rate constant for the (1 - exp(-KT * t))
    # delay). All from Part II Table 2 haloperidol column.
    emax_pos <- 0.41; label("Maximum haloperidol drug-effect fraction on PANSS positive (unitless)")          # Part II Table 2: Emax PANSS positive = 0.41 (95% CI 0.25-0.69)
    ec50_pos <- 1.2;  label("Haloperidol Css EC50 for PANSS positive (ng/mL)")                                # Part II Table 2: EC50 PANSS positive = 1.2 ng/mL (95% CI 0.22-3.12)
    kt_pos   <- 0.11; label("Onset rate constant of haloperidol effect on PANSS positive (1/day)")            # Part II Table 2: KT PANSS positive = 0.11 1/day (95% CI 0.07-0.18)

    emax_neg <- 0.21; label("Maximum haloperidol drug-effect fraction on PANSS negative (unitless)")          # Part II Table 2: Emax PANSS negative = 0.21 (95% CI 0.11-0.38)
    ec50_neg <- 6.4;  label("Haloperidol Css EC50 for PANSS negative (ng/mL)")                                # Part II Table 2: EC50 PANSS negative = 6.4 ng/mL (95% CI 3.7-13.9)
    kt_neg   <- 0.19; label("Onset rate constant of haloperidol effect on PANSS negative (1/day)")            # Part II Table 2: KT PANSS negative = 0.19 1/day (95% CI 0.14-0.27)

    emax_gen <- 0.27; label("Maximum haloperidol drug-effect fraction on PANSS general (unitless)")           # Part II Table 2: Emax PANSS general = 0.27 (95% CI 0.17-0.48)
    ec50_gen <- 2.58; label("Haloperidol Css EC50 for PANSS general (ng/mL)")                                 # Part II Table 2: EC50 PANSS general = 2.58 ng/mL (95% CI 0.73-6.31)
    kt_gen   <- 0.15; label("Onset rate constant of haloperidol effect on PANSS general (1/day)")             # Part II Table 2: KT PANSS general = 0.15 1/day (95% CI 0.09-0.23)

    # =============== IIV on PD (per subscale) ================================
    # Part II Table 1 / Table 2 reports IIV as either CV% (translated below
    # via omega^2 = log(1 + CV^2) for log-normal IIV on baseline) or as SD
    # (for normal-distributed IIV on Pmax and Emax that allow negative
    # values).
    etabasl_pos ~ 0.0515   # Part II Table 1: IIV BASL positive 23% CV; omega^2 = log(1 + 0.23^2) = 0.0515
    etabasl_neg ~ 0.0473   # Part II Table 1: IIV BASL negative 22% CV; omega^2 = log(1 + 0.22^2) = 0.0473
    etabasl_gen ~ 0.0319   # Part II Table 1: IIV BASL general 18% CV; omega^2 = log(1 + 0.18^2) = 0.0319

    etapmax_pos ~ 0.0576   # Part II Table 1: IIV Pmax positive SD = 0.24 (additive normal); variance 0.24^2 = 0.0576
    etapmax_neg ~ 0.0289   # Part II Table 1: IIV Pmax negative SD = 0.17; variance 0.17^2 = 0.0289
    etapmax_gen ~ 0.0441   # Part II Table 1: IIV Pmax general  SD = 0.21; variance 0.21^2 = 0.0441

    etaemax_pos ~ 0.0900   # Part II Table 2 haloperidol: IIV Emax positive SD = 0.30; variance 0.30^2 = 0.09
    etaemax_neg ~ 0.0625   # Part II Table 2 haloperidol: IIV Emax negative SD = 0.25; variance 0.25^2 = 0.0625
    etaemax_gen ~ 0.0625   # Part II Table 2 haloperidol: IIV Emax general  SD = 0.25; variance 0.25^2 = 0.0625

    # IIV on EC50 was estimable only for the PANSS negative subscale (Part
    # II Table 2 footnote). For positive and general it was fixed at 50% CV
    # = log(1 + 0.5^2) = 0.2231 -- encoded as a fixed IIV to preserve the
    # paper's choice.
    etaec50_pos ~ fixed(0.2231)  # Part II Table 2 footnote: IIV EC50 positive and general fixed at 50% CV
    etaec50_gen ~ fixed(0.2231)  # Part II Table 2 footnote: 50% CV nominal value because parameter not estimable
    etaec50_neg ~ 2.1085         # Part II Table 2 haloperidol: IIV EC50 negative = 269% CV; omega^2 = log(1 + 2.69^2) = 2.1085

    # =============== Residual error on PANSS subscales =======================
    # Per-drug RUV from Part II Table 2 (joint PKPD model); Part II Table 1
    # also reports a placebo-arm-only RUV (1.94 / 1.75 / 2.9) but the joint
    # PKPD model values take precedence here because the model body includes
    # the drug-effect term.
    addSd_PANSS_pos <- 2.2; label("Additive residual error on PANSS positive (PANSS units)")                  # Part II Table 2 (haloperidol): RUV PANSS positive SD = 2.2 (95% CI 2.0-2.3)
    addSd_PANSS_neg <- 2.3; label("Additive residual error on PANSS negative (PANSS units)")                  # Part II Table 2 (haloperidol): RUV PANSS negative SD = 2.3 (95% CI 2.2-2.5)
    addSd_PANSS_gen <- 3.8; label("Additive residual error on PANSS general (PANSS units)")                   # Part II Table 2 (haloperidol): RUV PANSS general  SD = 3.8 (95% CI 3.6-4.0)
  })

  model({
    # The paper reports time-dependent PD parameters in DAYS (TD in days, KT
    # in 1/day) while units$time = "hour" for the PK ODE. Convert model time
    # t (hours since first dose) to days for use in the placebo Weibull and
    # the Emax onset.
    t_days <- t / 24

    # ---------------- Individual PK parameters -------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # ---------------- One-compartment ODE ------------------------------------
    # Part I described haloperidol as 2-cmt; the Part II implementation
    # collapses to 1-cmt for the steady-state Css computation that drives
    # the PD model. This is a deliberate scope simplification of the
    # PART-I-inherited PK (vignette Assumptions and deviations).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma haloperidol concentration in ng/mL: amount (mg) / Vc (L) gives
    # mg/L = ug/mL; multiply by 1000 to obtain ng/mL.
    Cc <- 1000 * central / vc

    # ---------------- Individual PD parameters -------------------------------
    basl_pos_i <- basl_pos * exp(etabasl_pos)  # log-normal IIV on baseline
    basl_neg_i <- basl_neg * exp(etabasl_neg)
    basl_gen_i <- basl_gen * exp(etabasl_gen)

    pmax_pos_i <- pmax_pos + etapmax_pos       # additive normal IIV on Pmax (allows worsening / improvement)
    pmax_neg_i <- pmax_neg + etapmax_neg
    pmax_gen_i <- pmax_gen + etapmax_gen

    emax_pos_i <- emax_pos + etaemax_pos       # additive normal IIV on Emax (allows worsening / improvement)
    emax_neg_i <- emax_neg + etaemax_neg
    emax_gen_i <- emax_gen + etaemax_gen

    ec50_pos_i <- ec50_pos * exp(etaec50_pos)  # log-normal IIV on EC50
    ec50_neg_i <- ec50_neg * exp(etaec50_neg)
    ec50_gen_i <- ec50_gen * exp(etaec50_gen)

    # ---------------- Placebo Weibull (per subscale) -------------------------
    Pplacebo_pos <- pmax_pos_i * (1 - exp(-(t_days / td_pos)^pow_pos))
    Pplacebo_neg <- pmax_neg_i * (1 - exp(-(t_days / td_neg)^pow_neg))
    Pplacebo_gen <- pmax_gen_i * (1 - exp(-(t_days / td_gen)^pow_gen))

    # ---------------- Drug effect (per subscale) -----------------------------
    # Emax * Cc / (EC50 + Cc) * (1 - exp(-KT * t)) -- the (1 - exp(-KT t))
    # factor encodes the time delay to reach maximum drug effect.
    Drug_pos <- emax_pos_i * Cc / (ec50_pos_i + Cc) * (1 - exp(-kt_pos * t_days))
    Drug_neg <- emax_neg_i * Cc / (ec50_neg_i + Cc) * (1 - exp(-kt_neg * t_days))
    Drug_gen <- emax_gen_i * Cc / (ec50_gen_i + Cc) * (1 - exp(-kt_gen * t_days))

    # ---------------- Combined PANSS subscale prediction ---------------------
    # PANSS_subscale = BASL * (1 - Pplacebo) * (1 - Drug)
    PANSS_pos <- basl_pos_i * (1 - Pplacebo_pos) * (1 - Drug_pos)
    PANSS_neg <- basl_neg_i * (1 - Pplacebo_neg) * (1 - Drug_neg)
    PANSS_gen <- basl_gen_i * (1 - Pplacebo_gen) * (1 - Drug_gen)

    # ---------------- Multi-output observation model -------------------------
    # PK (haloperidol plasma concentration) and three PANSS subscales as
    # algebraic observables. Multi-output residual error follows the
    # parameter-name-then-output-suffix convention.
    Cc        ~ prop(propSd)
    PANSS_pos ~ add(addSd_PANSS_pos)
    PANSS_neg ~ add(addSd_PANSS_neg)
    PANSS_gen ~ add(addSd_PANSS_gen)
  })
}
