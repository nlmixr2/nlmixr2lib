Maseda_2018_micafungin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for micafungin in adult nonobese critically ill,",
    "obese noncritically ill, and morbidly obese critically ill patients (Maseda 2018), fitted",
    "NONPARAMETRICALLY with the nonparametric adaptive grid (NPAG) algorithm in Pmetrics. Clearance",
    "carries a joint power covariate model on total body weight (70 kg reference) and age (60 year",
    "reference), both with an exponent of 0.75. The source parameterises distribution as a central",
    "volume plus the first-order rate constants kcp and kpc; those are carried here as k12 and k21 and",
    "converted to the equivalent intercompartmental clearance and peripheral volume so the model is",
    "usable in either parameterisation. NPAG estimates a discrete joint distribution of individual",
    "parameters that has no closed form; it is approximated here by independent lognormal marginals",
    "matched to the reported coefficients of variation, so the shape of the joint density and any",
    "parameter correlations are not recoverable from this encoding. The source reports no",
    "residual-error model, so the residual term is present but fixed to zero. Sibling micafungin",
    "models: modellib('Martial_2017_micafungin') and modellib('Leroux_2018_micafungin').",
    sep = " "
  )
  reference <- paste(
    "Maseda E, Grau S, Luque S, Castillo-Mafla MP, Suarez-de-la-Rica A, Montero-Feijoo A,",
    "Salgado P, Gimenez MJ, Garcia-Bernedo CA, Gilsanz F, Roberts JA.",
    "Population pharmacokinetics/pharmacodynamics of micafungin against Candida species in obese,",
    "critically ill, and morbidly obese critically ill patients.",
    "Crit Care. 2018;22(1):94. doi:10.1186/s13054-018-2019-8.",
    "Structural model and covariate equation from Results/'Pharmacokinetic model';",
    "parameter estimates from Table 2; demographics from Table 1.",
    sep = " "
  )
  vignette <- "Maseda_2018_micafungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Maseda 2018 Methods/'Drug assay' (total
  # micafungin in plasma, measured by UHPLC-MS/MS over 0.2-30 ug/ml) and
  # Results/'Pharmacokinetic model' (two-compartment model with zero-order input
  # into the central compartment).
  compartmentData <- list(
    central = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "micafungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results/'Pharmacokinetic model': the final covariate equation is",
        "CL = TVCL * (Wt/70)^0.75 * (Age/60)^0.75, where 'Wt is the total body weight (kg)'. Total",
        "body weight, NOT an adjusted or ideal weight -- the cohort was selected for obesity and the",
        "paper's dosing conclusions are stated per total body weight. Table 1: median 95 kg",
        "(range 44-193) overall; 65 kg (44.0-92.5) in the critically ill nonobese stratum and",
        "157.5 kg (142-170) in the morbidly obese stratum receiving 150 mg. The Monte-Carlo dosing",
        "simulations spanned 45, 80, 115, 150 and 185 kg.",
        sep = " "
      ),
      source_name = "Wt"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results/'Pharmacokinetic model': second term of the final covariate equation,",
        "CL = TVCL * (Wt/70)^0.75 * (Age/60)^0.75, where 'Age is the patient's age (years)'.",
        "Methods text describes it as 'age (normalized to 60 years old to an exponential value of",
        "0.75)'. Note the sign: as printed, clearance INCREASES with age, which is the opposite of the",
        "usual renal-maturation-then-decline direction; see the vignette Errata. Table 1: median",
        "58 years (range 27-85). Simulations spanned 30, 50, 70 and 90 years.",
        sep = " "
      ),
      source_name = "Age"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "unitless",
      type = "binary",
      notes = "Methods/'Population pharmacokinetics covariate screening' lists gender among the screened covariates; not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened per Methods; not retained. Table 1 median 34.7 kg/m^2 (range 19.6-60.0). Used only to stratify the cohort and to set the clinical dosing rule (100 mg for BMI <= 45, 150 mg for BMI > 45)."
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened per Methods; not retained. Table 1 median 1.0 mg/dL (range 0.4-3.9)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Both measured and Cockcroft-Gault estimated creatinine clearance were screened per Methods; neither was retained. Table 1 mean 93.4 +/- 51.4 mL/min/1.73m^2. Micafungin is not renally cleared, so this is the expected result."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/dL",
      type = "continuous",
      notes = "Screened per Methods; not retained. Table 1 median 3 g/dL (range 1.2-4.0)."
    ),
    APACHE2 = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units = "unitless",
      type = "continuous",
      notes = "Screened per Methods; not retained. Discussion: 'Unlike the introduction of the severity score as a covariate, introducing the patient's age improved the model.'"
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "unitless",
      type = "continuous",
      notes = "Collected (Table 1, median 6, range 0-12) and discussed as a severity covariate; not retained in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    n_observations = 242L,
    age_range = "27-85 years (median 58)",
    age_median = "58 years",
    weight_range = "44-193 kg (median 95)",
    weight_median = "95 kg",
    bmi_range = "19.6-60.0 kg/m^2 (median 34.7)",
    sex_female_pct = 71.0,
    race_ethnicity = "Not reported (two-centre Spanish cohort).",
    disease_state = paste(
      "Adults receiving micafungin as empirical or directed treatment for invasive candidiasis, in",
      "three strata (Table 1): 11 morbidly obese critically ill patients (Hospital Universitario La",
      "Paz, Madrid), 10 nonobese critically ill patients and 10 obese noncritically ill patients",
      "(Hospital del Mar, Barcelona). Patients admitted to the ICU were those considered critically",
      "ill. Overall severity: SOFA median 6 (range 0-12), SAPS II median 34 (range 9-57), Candida",
      "score median 3 (range 2-4, for the 21 patients without microbiologically documented infection).",
      "Baseline renal function was largely preserved (creatinine clearance 93.4 +/- 51.4",
      "mL/min/1.73m^2) and albumin was low (median 3 g/dL, range 1.2-4.0).",
      sep = " "
    ),
    dose_range = paste(
      "Micafungin 100 mg or 150 mg once daily, diluted in 100 mL isotonic saline and infused",
      "intravenously over 60 min, at the discretion of the treating physician. In practice 100 mg was",
      "given for BMI <= 45 kg/m^2 and 150 mg for BMI > 45 kg/m^2, with three documented exceptions",
      "(two noncritically ill obese patients with BMI ~35 kg/m^2 received 150 mg; one critically ill",
      "patient with BMI > 45 kg/m^2 received 100 mg).",
      sep = " "
    ),
    regions = "Spain (Hospital Universitario La Paz, Madrid; Hospital del Mar, Barcelona).",
    indication = "Empirical or directed treatment of invasive candidiasis.",
    notes = paste(
      "Sampling (Methods): on day 3, at baseline (predose) and 1, 3, 5, 8, 18 and 24 h after the dose,",
      "with additional day-0 and day-7 samples when feasible; 242 total plasma concentrations entered",
      "the model. Total (not unbound) micafungin was assayed by UHPLC-MS/MS over 0.2-30 ug/ml.",
      "Estimation used the NPAG algorithm in Pmetrics for R; model selection was by -2 log-likelihood",
      "(structural -2LL 596.4; weight alone 595.2, p = 0.0586; age alone 587.5, p = 0.597; weight and",
      "age together 415.6, p = 0.0238). Internal validation used a 1000-sample bootstrap and NPDE.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: Table 2, 'Estimated micafungin parameters'. NPAG
    # returns a discrete joint density over support points rather than a point
    # estimate plus a parametric variance; Table 2 summarises that density as
    # mean +/- SD, %CV, variance and median. The MEAN column is taken as the
    # typical value here -- see the note below and the vignette Errata, where
    # the paper's own Table 3 PTA values are used to discriminate the mean and
    # median readings (the mean reading reproduces them, the median does not).
    #
    # Table 2 is internally consistent: for every row the reported %CV equals
    # 100*SD/mean and the reported 'Variance' equals SD^2 (e.g. clearance,
    # 0.49/0.80 = 61.3% vs 61.78 reported, 0.49^2 = 0.240 vs 0.24 reported), so
    # the SD is on the natural scale of each parameter, not a log scale.
    # ------------------------------------------------------------------------
    lcl <- log(0.80)
    label("Typical clearance TVCL at 70 kg and 60 years (L/h)")
    # Table 2, 'Clearance (l/h)' row: mean 0.80, SD 0.49, median 0.73. This is
    # TVCL of the Results covariate equation
    # CL = TVCL * (Wt/70)^0.75 * (Age/60)^0.75, so it is the clearance of a
    # 70 kg, 60-year-old subject, NOT of the median study subject (95 kg,
    # 58 years), whose typical clearance is 0.80*(95/70)^0.75*(58/60)^0.75 =
    # 0.98 L/h.
    lvc <- log(16.34)
    label("Central volume (L)")
    # Table 2, 'Central volume (l)' row: mean 16.34, SD 5.87, median 16.34. No
    # covariate is reported on volume -- Results names weight and age as
    # covariates 'for micafungin clearance' only.
    lk12 <- log(0.38)
    label("Transfer rate constant central -> peripheral1, kcp (1/h)")
    # Table 2, 'k cp (h-1)' row: mean 0.38, SD 0.37, median 0.26. Table 2
    # footnote: 'kcp, rate constant for drug distribution from the central to
    # peripheral compartment'.
    lk21 <- log(0.32)
    label("Transfer rate constant peripheral1 -> central, kpc (1/h)")
    # Table 2, 'k pc (h-1)' row: mean 0.32, SD 0.31, median 0.14. Table 2
    # footnote: 'kpc, rate constant for drug distribution from the peripheral to
    # central compartment'.

    # ------------------------------------------------------------------------
    # Covariate effects on clearance. Results/'Pharmacokinetic model' prints the
    # final covariate model as
    #     Micafungin CL = TVCL * (Wt/70)^0.75 * (Age/60)^0.75
    # (the two 0.75 superscripts are dropped by PDF-to-text conversion of the
    # display equation; they are recovered from the raw PDF text stream, and the
    # Methods prose agrees: 'body weight (normalized to 70 kg) and age
    # (normalized to 60 years old to an exponential value of 0.75) for micafungin
    # clearance'). Both exponents are structural constants written into the
    # equation, not estimated quantities -- Table 2 reports no estimate, SD or
    # %CV for either -- so both are encoded as fixed().
    # ------------------------------------------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Power exponent on (WT/70) for CL (unitless)")
    # Results covariate equation, first term.
    e_age_cl <- fixed(0.75)
    label("Power exponent on (AGE/60) for CL (unitless)")
    # Results covariate equation, second term. Positive, i.e. clearance rises
    # with age as printed; see covariateData$AGE and the vignette Errata.

    # ------------------------------------------------------------------------
    # Between-subject variability. NPAG places every model parameter in the
    # joint density, and Table 2 reports a %CV for all four; each is encoded
    # here as an independent lognormal marginal using the exact conversion
    #     omega^2 = log(CV^2 + 1)
    # which reproduces each reported %CV exactly. Two limitations are inherent
    # to this approximation and are recorded in the vignette Errata: (a) the
    # NPAG joint density is discrete and generally not lognormal, and its
    # parameter correlations are not reported, so the marginals are taken as
    # independent; (b) a lognormal cannot match a reported mean, median and %CV
    # simultaneously when the two centres disagree, and for kpc they disagree
    # strongly (mean 0.32 vs median 0.14, a ratio of 2.3 against the 1.39 a
    # lognormal with this %CV implies).
    # ------------------------------------------------------------------------
    etalcl ~ 0.32330
    # Table 2 clearance %CV = 61.78; log(0.6178^2 + 1) = 0.32330.
    etalvc ~ 0.12155
    # Table 2 central volume %CV = 35.95; log(0.3595^2 + 1) = 0.12155.
    etalk12 ~ 0.67075
    # Table 2 kcp %CV = 97.76; log(0.9776^2 + 1) = 0.67075.
    etalk21 ~ 0.65825
    # Table 2 kpc %CV = 96.51; log(0.9651^2 + 1) = 0.65825.

    # ------------------------------------------------------------------------
    # Residual error. The source reports NO residual-error model: Pmetrics
    # carries assay noise as a fixed error polynomial supplied by the analyst,
    # and neither its coefficients nor any estimated proportional/additive term
    # appear in the paper, its tables or its figures (the paper has no
    # supplement -- EuropePMC reports hasSuppl 'N'). Per the standing policy for
    # an unreported residual term with the structural values present, the term
    # is declared and fixed to zero rather than invented; Cc is therefore an
    # individual prediction with no measurement noise. The only quantitative
    # anchor the paper gives is the assay range, 0.2-30 ug/ml
    # (Methods/'Drug assay').
    # ------------------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction); not reported by the source")
  })

  model({
    # 1. Individual PK parameters. Clearance carries the joint weight/age power
    #    covariate model from Results; the central volume and the two
    #    distribution rate constants carry no covariate.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AGE / 60)^e_age_cl
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. Clearance-parameterisation equivalents of the source's rate-constant
    #    parameterisation, carried as derived quantities so the familiar
    #    intercompartmental clearance and peripheral volume are available
    #    without re-deriving them: Q = k12*Vc and Vp = Q/k21 = Vc*k12/k21 give
    #    the two-compartment system identical to the source's kcp/kpc form.
    kel <- cl / vc
    q <- k12 * vc
    vp <- q / k21

    # 3. ODE system. Two-compartment intravenous disposition. Results: 'A
    #    two-compartment linear model (including zero order input of drug into
    #    the central compartment) best described the time course'; the zero-order
    #    input is the clinical 60-min infusion (Methods), which is supplied
    #    through the event table as an infusion rate rather than as an estimated
    #    duration, because the source estimates no input-duration parameter.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation and error. Total micafungin in plasma; dose in mg over
    #    volume in L gives mg/L, which equals the paper's ug/ml.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
