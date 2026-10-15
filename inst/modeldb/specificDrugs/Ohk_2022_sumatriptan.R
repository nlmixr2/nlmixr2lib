Ohk_2022_sumatriptan <- function() {
  description <- "One-compartment population PK model for oral sumatriptan in healthy Korean adults (Ohk 2022): two parallel absorption routes (first-order absorption with lag time, and a transit-compartment chain with the Savic 2007 analytical input form) into a single central compartment with linear elimination, with separate typical apparent clearances for males and females."
  reference <- paste(
    "Ohk B, Seong S, Lee J, Gwon M, Kang W, Lee H, Yoon Y, Yoo H.",
    "Evaluation of sex differences in the pharmacokinetics of oral sumatriptan",
    "in healthy Korean subjects using population pharmacokinetic modeling.",
    "Biopharm Drug Dispos. 2022;43(1):24-31.",
    "doi:10.1002/bdd.2307.",
    sep = " "
  )
  vignette <- "Ohk_2022_sumatriptan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Selects between the separate male and female typical values of CL/F (Ohk 2022 Table 3 model 4 and Table 4). Ohk 2022 Section 2.5 codes SEX as 0 = male, 1 = female, which matches the canonical SEXF polarity directly (no value inversion needed). Time-invariant.",
      source_name = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant on CL/F alone (Table 3 model 3, dOFV = -10.541) but dropped once sex was on CL/F (model 5 vs model 4, dOFV = -0.001); not in the final model."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F and V/F (Table 3 models 2 and 9); not significant."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened on CL/F (with sex) and V/F (Table 3 models 6 and 12); not significant."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on CL/F (with sex) and V/F (Table 3 models 7 and 13); not significant."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on CL/F (with sex) and V/F (Table 3 models 8 and 14); not significant."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sumatriptan", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "sumatriptan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sumatriptan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 38L,
    n_studies = 2L,
    age_range = "21-31 years",
    age_median = "24.8 years (mean)",
    weight_range = "50-84 kg",
    weight_median = "64.0 kg (mean)",
    sex_female_pct = 23.7,
    race_ethnicity = "Korean",
    disease_state = "Healthy adult volunteers",
    dose_range = "Single 50 mg oral dose of sumatriptan (Sumatran 50 mg tablet) under fasting conditions",
    regions = "Republic of Korea (Kyungpook National University Hospital Clinical Trial Center, Daegu)",
    notes = "Pooled analysis of two single-dose studies run with the same protocol: the male cohort previously analysed by Lee 2015 (Transl Clin Pharmacol 23:66) and an additional study (CRIS KCT0001784) run to evaluate sex differences. 29 males and 9 females; 532 plasma concentrations; sampling to 12 h. Height 171.1 +/- 7.8 cm, BMI 21.8 +/- 2.2 kg/m^2, CrCL 119.2 +/- 17.9 mL/min (Table 1 and Section 3.1)."
  )

  # Implementation notes (see the vignette's 'Assumptions and deviations'):
  # * Structure follows Lee 2015 (same group, same protocol), which Ohk 2022
  #   cites as its basic model (Section 2.4 and Figure 1): depot carries the
  #   first-order arm (fraction 1 - f, rate ka1, lag ALAG1); the transit arm
  #   (fraction f) enters depot2 via the Savic 2007 analytical input with
  #   ktr = (n + 1) / MTT (Section 2.4 equations) and depot2 then absorbs
  #   into central at rate ka2.
  # * rxode2's transit(n, mtt, bio) evaluates the same gamma kernel for a
  #   non-integer n (NN = 6.7 is estimated as a continuous quantity). It is
  #   compartment-aware: a user data set must carry two dose records per
  #   administration with the same amount, one to depot and one to depot2.
  #   f(depot) = 1 - fr and f(depot2) = 0 perform the dose split.

  ini({
    # Structural parameters - Ohk 2022 Table 4 ('Parameter estimates from the
    # final population pharmacokinetic model'). BSV enters exponentially
    # (Section 2.4, Table 3 model 4: CL = ((1 - SEX) * theta1,male +
    # SEX * theta1,female) * exp(eta1)), so typical values are log-transformed.
    lcl_male <- log(444); label("Apparent clearance CL/F in males (L/h)") # Table 4 CL/F M = 444 L/h (RSE 4%)
    lcl_female <- log(281); label("Apparent clearance CL/F in females (L/h)") # Table 4 CL/F F = 281 L/h (RSE 7%)
    lvc <- log(68.7); label("Apparent volume of distribution V/F (L)") # Table 4 V3/F = 68.7 L (RSE 32%)
    lka1 <- log(0.568); label("First-order absorption rate constant from depot (1/h)") # Table 4 ka1 = 0.568 1/h (RSE 12%)
    lka2 <- log(0.295); label("Absorption rate constant from the final transit compartment to central (1/h)") # Table 4 ka2 = 0.295 1/h (RSE 5%)
    lmtt <- log(1.52); label("Mean transit time of the transit-compartment chain (h)") # Table 4 MTT = 1.52 h (RSE 11%)
    lntr <- log(6.7); label("Number of transit compartments (continuous)") # Table 4 NN = 6.7 (RSE 48%)
    ltlag <- log(0.239); label("Lag time for the first-order absorption arm (h)") # Table 4 ALAG1 = 0.239 h (RSE 1%)
    lfr <- log(0.558); label("Fraction of the dose absorbed through the transit-compartment arm (unitless)") # Table 4 f = 0.558 (RSE 7%)

    # Between-subject variability: Table 4 reports BSV as CV% for an
    # exponential (log-normal) eta (Section 2.4), converted here with
    # omega^2 = log(1 + CV^2). No BSV on ka1, NN, ALAG1 or f ('-' in Table 4);
    # no correlation reported, so the etas are independent.
    etalcl ~ 0.03016 # log(1 + 0.175^2); Table 4 BSV CL/F = 17.5 CV% (RSE 25%, shrinkage 7%)
    etalvc ~ 0.69215 # log(1 + 0.999^2); Table 4 BSV V3/F = 99.9 CV% (RSE 30%, shrinkage 2%)
    etalka2 ~ 0.05172 # log(1 + 0.2304^2); Table 4 BSV ka2 = 23.04 CV% (RSE 31%, shrinkage 17%)
    etalmtt ~ 0.16554 # log(1 + 0.4243^2); Table 4 BSV MTT = 42.43 CV% (RSE 32%, shrinkage 12%)

    # Residual error: combined proportional and additive (Section 2.4).
    propSd <- 0.249; label("Proportional residual error SD (fraction)") # Table 4 proportional error = 0.249 (RSE 6%)
    addSd <- 0.276; label("Additive residual error SD (ng/mL)") # Table 4 additive error = 0.276 (RSE 21%)
  })

  model({
    # Individual parameters. SEXF = 1 is female and SEXF = 0 is male,
    # matching the source SEX coding (Section 2.5).
    cl <- exp(lcl_female * SEXF + lcl_male * (1 - SEXF) + etalcl)
    vc <- exp(lvc + etalvc)
    ka1 <- exp(lka1)
    ka2 <- exp(lka2 + etalka2)
    mtt <- exp(lmtt + etalmtt)
    ntr <- exp(lntr)
    tlag <- exp(ltlag)
    fr <- exp(lfr)

    kel <- cl / vc

    # Figure 1 scheme. transit(ntr, mtt, fr) is the Savic analytical input
    # f * Dose * ktr * (ktr * t)^n / n! * exp(-ktr * t), ktr = (n + 1) / MTT
    # (Section 2.4 equations).
    d/dt(depot) <- -ka1 * depot
    d/dt(depot2) <- transit(ntr, mtt, fr) - ka2 * depot2
    d/dt(central) <- ka1 * depot + ka2 * depot2 - kel * central

    f(depot) <- 1 - fr
    f(depot2) <- 0
    alag(depot) <- tlag

    # central in mg and vc in L give mg/L; x 1000 gives ng/mL (the assay
    # units, Section 2.2).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
