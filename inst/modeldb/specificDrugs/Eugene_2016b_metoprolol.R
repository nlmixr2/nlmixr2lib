Eugene_2016b_metoprolol <- function() {
  description <- "One-compartment population PK model for oral S-metoprolol with first-order absorption and an absorption lag time in healthy young adults; every structural parameter is stratified by sex, and the kinetics are flip-flop (absorption is rate-limiting). Estimated by re-fitting the paper's own 100-subject Clinical Trial Simulation, so this is the variant of the Eugene 2016 Med Sci model that carries inter-individual variability (Eugene 2016)."
  reference <- paste(
    "Eugene AR. (2016).",
    "Metoprolol Dose Equivalence in Adult Men and Women Based on Gender",
    "Differences: Pharmacokinetic Modeling and Simulations.",
    "Med Sci (Basel) 4(4):18.",
    "doi:10.3390/medsci4040018.",
    sep = " "
  )
  vignette <- "Eugene_2016b_metoprolol_gender_dose_equivalence"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "S-metoprolol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "S-metoprolol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). Eugene 2016 Table 2 reports men and women as two fully separate parameter columns rather than a reference plus an effect; the male column is taken as the structural reference here and each female value is reproduced exactly by a log-additive offset.",
      notes = "Eugene 2016 Table 2 (Clinical Trial Simulation) reports Tlag, Ka, V and CL separately for men and for women, with no common reference value and no published covariate coefficient. To express this under the canonical SEXF (1 = female, 0 = male), the male column is the reference and each e_sexf_* term is computed as log(female value / male value), which reproduces the published female column to full printed precision. The single omega block is shared across sexes: Table 2 prints identical CV(%) values in the men and women columns (Tlag 0.20%, Ka 42%, V 43%, CL 55%) and a single un-stratified variance column, so the inter-individual variability is not sex-specific.",
      source_name = "gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    age_range = "20-36 years (the age range of the Luzier 1999 cohort whose digitized profiles underlie the model)",
    age_median = "not reported",
    weight_range = "men 83.9 +/- 10.7 kg, women 62.0 +/- 7.3 kg (mean +/- SD, Luzier 1999 as quoted in Eugene 2016 Discussion)",
    weight_median = "not reported (means reported instead)",
    sex_female_pct = 50,
    race_ethnicity = "not reported",
    disease_state = "Healthy volunteers",
    dose_range = "Single 100 mg oral dose of S-metoprolol per simulated subject (Eugene 2016 Methods 2.4). The underlying clinical study (Luzier 1999) administered nine 100 mg doses of racemic metoprolol every 12 h.",
    regions = "United States (Luzier 1999 cohort; Eugene 2016 analysis performed at Mayo Clinic, Rochester MN)",
    notes = "Simulated cohort: 50 men and 50 women, each given a single 100 mg dose of S-metoprolol, sampled at 0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 12, 14, 16, 18 and 24 h, yielding 1700 plasma samples (Eugene 2016 Methods 2.4 and Results 3.3). These parameters are therefore a re-estimate from SIMULATED data, not a fit to observed concentrations: the Clinical Trial Simulation drew subjects using the typical values of Eugene 2016 Table 1 together with per-sex coefficients of variation carried over from Luzier 1999 (S-CL: men 59%, women 49%; S-V: men 44%, women 34%) and an assumed 40% variability on Ka, and MONOLIX 4.3.3 (SAEM-MCMC) was then re-run on the simulated dataset. The recovered variabilities (Ka 42%, V 43%, CL 55%) sit close to those inputs, as expected. The companion model Eugene_2016b_metoprolol_enantiomers carries the Table 1 fit to the digitized Luzier 1999 data, which is the model that generates every published Cmax / Tmax / AUC / half-life value in Eugene 2016 Results 3.2."
  )

  ini({
    # Structural parameters, Eugene 2016 Table 2 ('Gender-stratified population
    # pharmacokinetic model estimates from the Clinical Trial Simulation (CTS)'),
    # Men column. The Men column is used as the structural reference; the Women
    # column is reproduced through the e_sexf_* offsets below.
    ltlag <- log(0.677); label("Absorption lag time in men (h)") # Eugene 2016 Table 2, Men: Tlag = 0.677 h (SE 0.0021, RSE 0%)
    lka <- log(0.233); label("First-order absorption rate constant in men (1/h)") # Eugene 2016 Table 2, Men: Ka = 0.233 1/h (SE 0.0058, RSE 2%)
    lvc <- log(49); label("Apparent central volume of distribution Vc/F in men (L)") # Eugene 2016 Table 2, Men: V = 49 L (SE 1.5, RSE 3%)
    lcl <- log(231); label("Apparent clearance CL/F in men (L/h)") # Eugene 2016 Table 2, Men: CL = 231 L/h (SE 10, RSE 4%)

    # Sex effects. Eugene 2016 Table 2 publishes the women's values directly
    # rather than a covariate coefficient, so each offset below is the log ratio
    # of the published Women value to the published Men value and reproduces the
    # Women column exactly.
    e_sexf_ltlag <- log(0.38 / 0.677); label("Log-additive female-sex effect on absorption lag time") # Eugene 2016 Table 2, Women: Tlag = 0.38 h (SE 0.00013, RSE 0%); exp(ltlag + e_sexf_ltlag) = 0.38
    e_sexf_lka <- log(0.149 / 0.233); label("Log-additive female-sex effect on Ka") # Eugene 2016 Table 2, Women: Ka = 0.149 1/h (SE 0.0037, RSE 2%); exp(lka + e_sexf_lka) = 0.149
    e_sexf_lvc <- log(33.3 / 49); label("Log-additive female-sex effect on Vc/F") # Eugene 2016 Table 2, Women: V = 33.3 L (SE 0.88, RSE 3%); exp(lvc + e_sexf_lvc) = 33.3
    e_sexf_lcl <- log(92.9 / 231); label("Log-additive female-sex effect on CL/F") # Eugene 2016 Table 2, Women: CL = 92.9 L/h (SE 4, RSE 4%); exp(lcl + e_sexf_lcl) = 92.9

    # Inter-individual variability, Eugene 2016 Table 2 lower block. The printed
    # values are variances on the log scale: MONOLIX's CV(%) column equals the
    # log-scale SD, and sqrt() of each variance reproduces the printed CV to the
    # printed precision -- sqrt(0.176) = 42.0% (Ka, printed 42%),
    # sqrt(0.182) = 42.7% (V, printed 43%), sqrt(0.305) = 55.2% (CL, printed 55%).
    # A single un-stratified block is reported and the Men and Women CV(%) columns
    # are identical, so the same variances apply to both sexes.
    etaltlag ~ 0.0003 # Eugene 2016 Table 2: 'omega Tlag, variance for Tlag' = 0.0003 (SE 0.0021, RSE 809%); essentially unidentified -- see the vignette Errata
    etalka ~ 0.176 # Eugene 2016 Table 2: 'omega Ka, variance for Ka' = 0.176 (SE 0.012, RSE 7%); printed CV 42%
    etalvc ~ 0.182 # Eugene 2016 Table 2: 'omega Vd, variance for V' = 0.182 (SE 0.013, RSE 7%); printed CV 43%
    etalcl ~ 0.305 # Eugene 2016 Table 2: 'omega CL, variance for CL' = 0.305 (SE 0.022, RSE 7%); printed CV 55%

    # Residual error. Eugene 2016 Table 2 reports a single row labelled
    # 'Proportional error model' = 0.0281, i.e. the MONOLIX constant-CV parameter
    # b, which is a fraction of the predicted concentration.
    propSd <- 0.0281; label("Proportional residual SD (fraction)") # Eugene 2016 Table 2: 'Proportional error model' = 0.0281 (SE 0.0007, RSE 2%)
  })

  model({
    # Individual parameters. Every structural parameter is sex-stratified, so the
    # SEXF term appears on all four.
    tlag <- exp(ltlag + e_sexf_ltlag * SEXF + etaltlag)
    ka <- exp(lka + e_sexf_lka * SEXF + etalka)
    vc <- exp(lvc + e_sexf_lvc * SEXF + etalvc)
    cl <- exp(lcl + e_sexf_lcl * SEXF + etalcl)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    lag(depot) <- tlag

    # Dose is in mg and vc is in L, so central / vc is mg/L = 1000 ng/mL.
    # Scale to ng/mL to match the paper's reported concentration units.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
