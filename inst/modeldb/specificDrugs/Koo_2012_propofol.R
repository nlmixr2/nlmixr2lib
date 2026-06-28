Koo_2012_propofol <- function() {
  description <- paste(
    "Pharmacodynamic sigmoid Emax model for the probability of recovery of",
    "consciousness (ROC) versus propofol effect-site concentration (Ce, ug/mL)",
    "during emergence from propofol-remifentanil target-controlled-infusion",
    "(TCI) general anesthesia in 94 ASA I-II adult patients undergoing",
    "elective minor eye or ENT surgery (Koo 2012). Age modulates both the",
    "effect-site concentration at 50% probability of ROC (Ce50) and the Hill",
    "exponent lambda via linear-additive age-centred forms",
    "Ce50 = 1.15 - 0.0128 * (AGE - 43) and lambda = 9.69 - 0.141 * (AGE - 43).",
    "Inter-individual variability is log-normal on Ce50 (CV 26.0%); IIV on",
    "lambda was dropped from the final model. Propofol PK is not fit in the",
    "source paper -- the per-record effect-site propofol concentration is",
    "supplied as the time-varying covariate CEFFECT (driven by the Schnider",
    "1998/1999 TCI controller in the source study; Keo = 0.459 /min). NONMEM",
    "Bernoulli LAPLACE likelihood in the source paper; this implementation",
    "exposes the typical-value probability of ROC with a placeholder additive",
    "residual error (see vignette Assumptions and deviations), following the",
    "Shin_2014_sevoflurane.R precedent from the same Yonsei research group.",
    sep = " "
  )
  reference <- paste(
    "Koo BN, Lee JR, Noh GJ, Lee JH, Kang YR, Han DW. (2012).",
    "A pharmacodynamic analysis of factors affecting recovery from anesthesia",
    "with propofol-remifentanil target controlled infusion.",
    "Acta Pharmacologica Sinica 33(8):1080-1084.",
    "doi:10.1038/aps.2012.85. Published online 30 Jul 2012.",
    sep = " "
  )
  vignette <- "Koo_2012_propofol"
  units <- list(
    time          = "minute",
    dosing        = "(not applicable; effect-site propofol concentration is supplied as the time-varying covariate CEFFECT)",
    concentration = "ug/mL (propofol effect-site concentration via CEFFECT)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters Ce50 and lambda via additive-linear",
        "centred forms Ce50 = 1.15 - 0.0128 * (AGE - 43) and",
        "lambda = 9.69 - 0.141 * (AGE - 43) per Koo 2012 Table 2 'Final' row.",
        "Cohort mean +/- SD: 42.8 +/- 16.5 years (Koo 2012 Table 1).",
        "Eligibility was >= 20 years per Materials and methods.",
        "Predicted Ce50 values at 25 / 50 / 75 years are 1.38 / 1.06 / 0.74",
        "ug/mL and predicted lambda values are 12.23 / 8.70 / 5.18,",
        "respectively (Koo 2012 Results). The linear additive form will give",
        "non-physical (zero / negative) Ce50 or lambda at extreme ages",
        "beyond the cohort range; the natural applicability window is the",
        "cohort age range and a modest extrapolation."
      ),
      source_name        = "AGE"
    ),
    CEFFECT = list(
      description        = "Time-varying propofol effect-site concentration (ug/mL) supplied as the PD driver of the sigmoid Emax probability of ROC",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "For this model the canonical effect-site PD driver CEFFECT carries",
        "the propofol effect-site concentration in ug/mL. In the source study,",
        "CEFFECT is the Ce predicted by the Schnider 1998/1999 propofol TCI",
        "controller (Orchestra Base Primea pump) operating against per-patient",
        "anthropometric input, with Keo = 0.459 /min (Koo 2012 Discussion).",
        "Cohort-observed CEFFECT range was approximately 4.4 +/- 1.1 ug/mL at",
        "loss of consciousness, 3.2 +/- 1.0 ug/mL at end of surgery, and",
        "1.1 +/- 0.3 ug/mL at return of consciousness (Koo 2012 Results).",
        "For downstream simulations, users supply CEFFECT either from a",
        "Schnider-style TCI simulation, from a propofol PK + effect-site",
        "model in the registry (e.g. modellib('Knibbe_2005_propofol_human')",
        "chained with an external Keo = 0.459 /min effect-site link), or",
        "from observed plasma concentrations equilibrated with an",
        "effect-site rate constant."
      ),
      source_name        = "(none; effect-site Ce is computed externally by the Schnider TCI controller and not stored as a named NONMEM data column with a standard alias)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 94L,
    n_studies       = 1L,
    age_range       = ">= 20 years (ASA I-II); cohort mean +/- SD 42.8 +/- 16.5 years (Koo 2012 Table 1)",
    weight_range    = "cohort mean +/- SD 71.1 +/- 14.4 kg (Koo 2012 Table 1)",
    height_range    = "cohort mean +/- SD 167 +/- 10.6 cm (Koo 2012 Table 1)",
    bmi_range       = "cohort mean +/- SD 24.8 +/- 4.4 kg/m^2 (Koo 2012 Table 1; BMI > 30 excluded by inclusion criterion)",
    sex_female_pct  = 100 * 41 / 94,
    disease_state   = paste(
      "ASA I-II adults scheduled for elective minor surgery (55 eye and",
      "39 ENT procedures) at the Eye and ENT Severance Hospital. Exclusion",
      "criteria: cardiac / pulmonary / hepatic / renal disease; hearing",
      "loss or other neurological deficit; allergy or adverse reaction to",
      "medication; CNS-affecting medications; BMI > 30."
    ),
    dose_range      = paste(
      "Propofol effect-site TCI (Orchestra Base Primea pump, Fresenius Vial,",
      "France) with the Schnider 1998/1999 PK + effect-site model; initial",
      "induction target Ce = 4 ug/mL with 0.5 ug/mL up-titration if loss of",
      "consciousness was not achieved; intraoperative Ce titrated to keep",
      "BIS between 40 and 60. Concurrent remifentanil effect-site TCI (Minto",
      "1997 model) initial Ce = 2 ng/mL adjusted to intraoperative",
      "hemodynamics. Glycopyrrolate (0.1 mg iv) premedication; rocuronium",
      "(0.6 mg/kg iv) for endotracheal intubation; ramosetron (0.3 mg) and",
      "ketorolac (60 mg) for PONV / pain prophylaxis pre-emergence."
    ),
    regions         = "Korea (Yonsei University College of Medicine, Seoul; Eye and ENT Severance Hospital)",
    notes           = paste(
      "Single-centre prospective study, enrolment January-September 2011.",
      "Duration of surgery 66.9 +/- 53.0 min; duration of anesthesia 97.5",
      "+/- 54.0 min. BIS at baseline 92.8 +/- 4.9, at LOC 67.0 +/- 15.1, at",
      "end of surgery 43.8 +/- 10.6, at ROC 75.7 +/- 6.0. Effect-site",
      "concentrations at LOC: propofol 4.4 +/- 1.1 ug/mL, remifentanil 2.0",
      "+/- 0.3 ng/mL; at end of surgery: 3.2 +/- 1.0 ug/mL and 2.3 +/- 0.4",
      "ng/mL; at ROC: 1.1 +/- 0.3 ug/mL and 0.8 +/- 1.0 ng/mL. Consciousness",
      "was assessed every 10 seconds during emergence (OAA/S = 3 threshold).",
      "Forward inclusion / backward elimination covariate selection in",
      "NONMEM VII (LIKELIHOOD LAPLACE METHOD=conditional); delta-OFV",
      "threshold 3.85 points (P < 0.05, df = 1). Only AGE met the selection",
      "criterion among sex, height, weight, BMI, propofol Ce at LOC,",
      "remifentanil Ce at ROC, duration of propofol infusion, and mean",
      "propofol / remifentanil doses (Koo 2012 Table 1 correlation",
      "analysis). Ethics committee approval Yonsei University Health System",
      "(4-2010-0580); written informed consent."
    )
  )

  ini({
    # Fixed-effect estimates from Koo 2012 Table 2 'Final' row. Source paper
    # reports values without standard errors or RSE; the values are stored
    # on the canonical log scale for the typical-value anchor but enter
    # model() in the printed-in-Table-2 additive-linear age form, NOT a
    # multiplicative power form. Centring age at 43 years (the cohort mean)
    # reproduces the Table 2 'Final' row equations exactly:
    #   Ce50 (ug/mL) = 1.15 - 0.0128 * (AGE - 43)
    #   lambda       = 9.69 - 0.141  * (AGE - 43)
    # %CV(Ce50) = 26.0 (log-normal); lambda has no IIV in the final model
    # (Table 2 %CV column shows '-' on lambda; dropped from the basic
    # model's 32.1% CV via forward-inclusion / backward-elimination
    # covariate selection per Methods).

    lec50 <- log(1.15)
    label("Log of effect-site propofol Ce50 at AGE = 43 years (ug/mL)")
    # Koo 2012 Table 2 'Final' Ce50 intercept at AGE = 43 = 1.15 ug/mL.

    e_age_ec50 <- -0.0128
    label("Additive linear AGE slope on Ce50 (ug/mL per year, centred at AGE = 43)")
    # Koo 2012 Table 2 'Final' Ce50 = 1.15 - 0.0128 * (AGE - 43).

    lhill <- log(9.69)
    label("Log of Hill exponent lambda at AGE = 43 years (unitless)")
    # Koo 2012 Table 2 'Final' lambda intercept at AGE = 43 = 9.69.

    e_age_hill <- -0.141
    label("Additive linear AGE slope on lambda (per year, centred at AGE = 43)")
    # Koo 2012 Table 2 'Final' lambda = 9.69 - 0.141 * (AGE - 43).

    # IIV on Ce50: log-normal, CV(%) = 26.0 (Koo 2012 Table 2 'Final' %CV
    # column on Ce50; Methods: 'The inter-individual random variability of
    # Ce50 and lambda was modeled using a log-normal model'). Internal
    # log-normal variance: omega^2 = log(CV^2 + 1) = log(0.26^2 + 1) =
    # 0.0654.
    etalec50 ~ 0.0654
    # Koo 2012 Table 2 'Final' %CV(Ce50) = 26.0; omega^2 = log(1 + 0.26^2).

    # Residual-error placeholder. Koo 2012 likelihood is Bernoulli ($EST
    # LIKELIHOOD LAPLACE METHOD=conditional in NONMEM VII; Methods,
    # paragraph after the L equation). nlmixr2 / rxode2 do not natively
    # express a Bernoulli observation for a probability output in ini() /
    # model(); the same Yonsei research group's Shin_2014_sevoflurane.R
    # extraction handled this by exposing the typical-value probability
    # with a placeholder additive residual error the user can tune at fit
    # time. Following that precedent here. See vignette 'Assumptions and
    # deviations'.
    addSd_prob_roc <- 0.05
    label("Placeholder additive residual error on the typical-value probability of ROC (unitless; NOT from Koo 2012)")
  })

  model({
    # 1. Per-record effect-site propofol concentration (ug/mL) from the
    #    canonical CEFFECT covariate column.
    conc <- CEFFECT

    # 2. Individual structural parameters. Both Ce50 and lambda use
    #    additive-linear AGE effects exactly as printed in Koo 2012 Table 2
    #    ('Final' row). IIV on Ce50 is log-normal (multiplicative on the
    #    linear-typical value), per Methods ('The inter-individual random
    #    variability of Ce50 and lambda was modeled using a log-normal
    #    model'); IIV on lambda was dropped in the final model.
    #    Precedent for the '(exp(lparam) + e * (cov - ref)) * exp(etalparam)'
    #    additive-on-linear with log-normal IIV encoding is
    #    Chi_2018_propofol.R (additive WT on CL) and other PD-on-Ce
    #    models in the registry.
    ec50 <- (exp(lec50) + e_age_ec50 * (AGE - 43)) * exp(etalec50)
    hill <- exp(lhill) + e_age_hill * (AGE - 43)

    # 3. Sigmoid Emax probability of recovery of consciousness (ROC). The
    #    source paper's formula is rendered as a non-decoded image in the
    #    journal PDF; the consistent textbook interpretation for the
    #    probability of an event PREVENTED by drug (here ROC suppressed by
    #    high propofol effect-site exposure) is the unit-E0 / unit-Emax /
    #    inhibitor sigmoid form
    #         P(ROC) = Ce50^lambda / (Ce50^lambda + Ce^lambda).
    #    P decreases monotonically from 1 at Ce = 0 to 0 as Ce -> infinity,
    #    with P = 0.5 at Ce = Ce50. The Shin_2014_sevoflurane.R model from
    #    the same Yonsei research group uses the identical algebraic form
    #    for the sevoflurane analog (ETSEVO as PD driver,
    #    'PROB = 1 - C^gamma / (CE50^gamma + C^gamma) = C50^gamma /
    #    (C50^gamma + C^gamma)' per its model body comment); the Koo 2012
    #    paper's reported Ce50 / lambda / age effects are consistent with
    #    this form across the 25 / 50 / 75 year predictions in Results.
    prob_roc <- ec50^hill / (ec50^hill + conc^hill)

    # 4. Observation. Placeholder additive residual on the typical-value
    #    probability; see ini() block comment and vignette 'Assumptions and
    #    deviations' for the deviation rationale relative to the source
    #    Bernoulli likelihood.
    prob_roc ~ add(addSd_prob_roc)
  })
}
