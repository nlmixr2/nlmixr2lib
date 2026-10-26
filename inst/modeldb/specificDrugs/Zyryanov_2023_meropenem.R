Zyryanov_2023_meropenem <- function() {
  description <- paste(
    "One-compartment IV population PK model for meropenem in 66 preterm",
    "newborns in a Moscow neonatal intensive care unit, fitted to peak-trough",
    "therapeutic drug monitoring data (Zyryanov 2023). Pmetrics NPAG",
    "non-parametric fit parameterised by an elimination rate constant and an",
    "absolute volume of distribution, with zero-order (infusion) input and",
    "no covariates in the model: gestational, postnatal and postconceptual",
    "age, body weight and Schwartz GFR were related to the individual",
    "Bayesian estimates by univariate regression only, and no coefficients",
    "were reported. The NPAG joint density is approximated by independent",
    "log-normal marginals at the published medians and CV%. The assay-error",
    "polynomial used for weighting was not reported, so residual error is",
    "carried as fixed zero.",
    sep = " "
  )
  reference <- paste(
    "Zyryanov S, Bondareva I, Butranova O, Kazanova A.",
    "Population PK/PD modelling of meropenem in preterm newborns based on",
    "therapeutic drug monitoring data.",
    "Front Pharmacol. 2023;14:1079680. doi:10.3389/fphar.2023.1079680",
    sep = " "
  )
  vignette <- "Zyryanov_2023_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Related to the individual Bayesian posterior Kel / V / T1/2 by univariate
  # regression (Statistical analysis; Results). No regression coefficient is
  # printed, and the population model itself carries no covariate.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at TDM sampling",
      units = "kg",
      type = "continuous",
      notes = "Significant univariate predictor of individual T1/2 and V (Results; R^2 8-40% across the significant regressions), coefficients not reported; not in the population model. Table 1: 1.66 +/- 0.464 kg (median 1.58, IQR 1.30-2.04). V and CL in this model are therefore absolute (L, L/h), not per-kg."
    ),
    PAGE = list(
      description = "Postconceptual (postmenstrual) age at TDM sampling",
      units = "months",
      type = "continuous",
      notes = "Source alias 'PCA', reported in WEEKS (Table 1: 33.6 +/- 3.69, median 33.6, IQR 30.7-35.0); the canonical carries months. Significant univariate predictor of individual T1/2 and V, coefficients not reported; not in the population model. Correlated with body weight (r = 0.84, Discussion)."
    ),
    CRCL = list(
      description = "Glomerular filtration rate by the Schwartz equation, k * height / Scr with k = 0.33 for preterm infants",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Source alias 'creatinine clearance (GFR)' (Methods, Schwartz equation after Brion 1986: height in cm, Scr in mg/dL). Significant univariate predictor of individual T1/2, coefficient not reported; not in the population model. Table 1: 27.31 +/- 9.676 (median 25.08, IQR 19.53-32.53)."
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened by univariate regression (Statistical analysis); not named among the significant predictors. Table 1: 28.9 +/- 2.93 weeks (median 28.5, IQR 27-31). Used only to split the PK/PD target-attainment results (GA < 32 vs 32-36 weeks, Figures 3 and 4)."
    ),
    PNA = list(
      description = "Postnatal age at TDM sampling",
      units = "months",
      type = "continuous",
      notes = "Reported in DAYS (Table 1: 33.0 +/- 20.88, median 30.5, IQR 15-47.3); the canonical carries months. Screened by univariate regression; not named among the significant predictors."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_studies = 1L,
    n_concentrations = 132L,
    age_range = "Preterm newborns; postnatal age at sampling 33.0 +/- 20.9 days (median 30.5, IQR 15.0-47.3)",
    ga_range = "Gestational age 28.9 +/- 2.93 weeks (median 28.5, IQR 27-31)",
    pma_range = "Postconceptual age at sampling 33.6 +/- 3.69 weeks (median 33.6, IQR 30.7-35.0)",
    weight_range = "1.66 +/- 0.464 kg (median 1.58, IQR 1.30-2.04) at sampling",
    height_range = "38.3 +/- 4.26 cm (median 39.0, IQR 35.0-40.3)",
    sex_female_pct = 39.4,
    race_ethnicity = "Not reported (single-centre Russian cohort)",
    disease_state = paste(
      "Preterm newborns in intensive care treated with meropenem for bacterial",
      "infection: pneumonia 84.8%, necrotizing enterocolitis 60.6%, congenital",
      "infection 18.2%, bacterial meningitis 9.1%, bronchopulmonary dysplasia",
      "18.2% (more than one per patient allowed). Excluded: renal replacement",
      "therapy, ECMO, and interacting co-medication within 3 days of TDM."
    ),
    renal_function = "Schwartz GFR 27.31 +/- 9.676 mL/min/1.73 m^2 (median 25.08); serum creatinine 46.66 +/- 20.235 umol/L (median 41.15) at TDM sampling (Table 1).",
    dose_range = paste(
      "Empirical meropenem 40-120 mg/kg/day (20-40 mg/kg per dose every 8-12 h)",
      "as 1-3 h IV infusions 2-3 times a day; daily dose at TDM 94.6 +/- 40.4",
      "mg/day (median 94.3). TDM on therapy day 4.2 +/- 1.77 (range 2-10)."
    ),
    regions = "Russia (City Clinical Hospital 24, Moscow); admitted January 2018 - February 2022",
    notes = paste(
      "Demographics from Zyryanov 2023 Table 1. One TDM occasion per patient:",
      "a peak within 30 min after dosing and a trough just before the next",
      "dose (132 samples). HPLC-UV assay, LLOQ 0.2 mg/L. Fitted with the",
      "Pmetrics NPAG algorithm; the final population is a discrete set of",
      "support points (Figure 1)."
    )
  )

  ini({
    # Structural parameters: Results paragraph after Table 1 and Abstract,
    # 'mean +/- SD (median)' of the NPAG population distribution. The MEDIAN
    # is encoded as the typical value (the median of a log-normal marginal is
    # exp(theta)). The medians, not the means, reproduce the paper's own
    # Table 2 derived medians: 0.693 / 0.30 = 2.31 h vs T1/2 median 2.33 h,
    # and 0.30 * 1.2 = 0.36 L/h vs CL median 0.355 L/h. The Results sentence
    # prints the Kel unit as 'L/h'; the Abstract's '1/h' is the correct one.
    lkel <- log(0.30); label("Elimination rate constant Kel (1/h)")    # Results / Abstract: Kel 0.31 +/- 0.13 (median 0.3) 1/h
    lvc  <- log(1.2);  label("Volume of distribution V (L)")           # Results / Abstract: V 1.2 +/- 0.4 (median 1.2) L

    # Between-subject variability: Results / Abstract 'interindividual
    # variability (CV) of 42% and 33%' for Kel and V (= SD / mean of the
    # NPAG distribution, 0.13 / 0.31 and 0.4 / 1.2). Carried as log-normal
    # marginals with omega^2 = log(CV^2 + 1). Correlations of the joint
    # density are not reported, so the etas are independent.
    etalkel ~ 0.162375 # Results: Kel CV 42%; log(0.42^2 + 1)
    etalvc  ~ 0.103347 # Results: V CV 33%; log(0.33^2 + 1)

    # Residual error. Methods 'PK analysis': 'For observation weighting for
    # PK modeling, the inverse of the assay variance was used.' This is the
    # Pmetrics assay-error polynomial SD = C0 + C1 * C (a linear sum,
    # nlmixr2's combined1() form), but neither the coefficients nor a gamma
    # or lambda term are reported anywhere in the paper, which has no
    # supplementary material. The assay validation figures (intra-/interday
    # precision 1.2-5.3%, LLOQ 0.2 mg/L) describe the assay, not a fitted
    # residual, so both terms are carried as fixed(0) rather than invented.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })
  model({
    # Methods 'PK analysis': linear one-compartment model with zero-order
    # input and first-order elimination, parameterised by V and Kel; the
    # individual CL is Kel * V and T1/2 is 0.693 / Kel.
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    cl <- kel * vc

    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L (total meropenem; the paper assumes 98%
    # unbound for the PK/PD targets).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
