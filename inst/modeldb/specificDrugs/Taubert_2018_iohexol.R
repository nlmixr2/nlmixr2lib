Taubert_2018_iohexol <- function() {
  description <- "Three-compartment intravenous population PK model for iohexol in elderly adults (aged 70 years or older) of the Berlin Initiative Study, developed to show that a third, rapidly equilibrating compartment is needed to describe early iohexol concentrations and that the customary two-compartment (Schwartz) slope-intercept approach overestimates iohexol clearance, i.e. measured GFR. Linear (first-order) elimination from the central compartment, no covariates, log-normal IIV on all six structural parameters, and a combined additive + proportional residual error. Taubert 2018, 563 of 570 patients analysed, 8-10 samples each after a single 3235 mg IV dose."
  reference <- "Taubert M, Ebert N, Martus P, van der Giet M, Fuhr U, Schaeffner E. Using a three-compartment model improves the estimation of iohexol clearance to assess glomerular filtration rate. Sci Rep. 2018;8(1):17723. doi:10.1038/s41598-018-35989-x. PMC6286374."
  vignette <- "Taubert_2018_iohexol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Iohexol is given as a single intravenous injection of 3235 mg and total
  # plasma iohexol is measured (Methods; Figure 1 y-axis 'Iohexol plasma
  # concentration [mg/L]').
  compartmentData <- list(
    central = list(analyte = "iohexol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iohexol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "iohexol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 563L,
    n_studies = 1L,
    age_range = "70 years or older (Berlin Initiative Study inclusion criterion, Methods); median 76.9 years (Discussion).",
    weight_range = "Not reported in this paper (body weight and height were available but were not covariates in the final model).",
    race_ethnicity = "Not reported; participants were randomly sampled from a large German health insurance database.",
    disease_state = "Community-dwelling elderly adults of the Berlin Initiative Study (BIS) who agreed to an iohexol measured-GFR investigation; renal function ranged from normal to reduced (individual iohexol clearance EBEs 15-136 mL/min). Median serum creatinine 0.91 mg/dL and cystatin C 1.05 mg/L (Discussion).",
    renal_function = "Individual iohexol clearance (empirical Bayes estimates) 15-136 mL/min, median 61 mL/min (Results).",
    dose_range = "Single intravenous injection of 3235 mg iohexol (5 mL Omnipaque 300); syringes were not weighed (Discussion).",
    regions = "Germany (Berlin)",
    sampling = "Plasma samples at 30, 60, 90, 120, 150, 180, 240 and 300 min post-dose in all subjects; additional 10 and 20 min samples in 52 subjects (Methods).",
    notes = "570 patients with complete iohexol measurement were available; 7 with implausible concentration-time curves (suspected para-vascular administration or implausibly fast early decline) were excluded, leaving 563 (Results). Estimation with NONMEM 7.4.2 FOCE; individual clearances were re-estimated with the NONMEM nonparametric option for the BIS-equation step, which is not part of this parametric model."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters: Supp. Table 1 ('Parameter estimates for the
    # iohexol base model as obtained from bootstrap statistics with 1000
    # samples'), median column. The final model adds no covariates to this
    # three-compartment base (Results: 'The final model comprised three
    # compartments, first-order elimination, a combined error model and a
    # full variance-covariance matrix of random effects'), and the Results
    # quote the same values (CL 57.4 mL/min = 3.44 L/h; Q 12.0 L/h with V
    # 3.3 L; Q 3.5 L/h with V 5.4 L). The paper prints only bootstrap
    # medians, no separate NONMEM point estimates.
    #
    # Supp. Table 1 V2/Q1 is the slowly equilibrating peripheral
    # (peripheral1) and V3/Q2 the rapidly equilibrating one (peripheral2);
    # Results: 'Equilibration to the third compartment ... (Q) of 12.0 L/h
    # and peripheral volume (V) of 3.3 L, was fast compared to the second
    # compartment (Q of 3.5 L/h, V of 5.4 L)'.
    # ---------------------------------------------------------------------
    lcl <- log(3.44); label("Iohexol clearance CL (L/h)")                                  # Supp. Table 1: CL 3.44 L/h [3.33-3.55]; Results: 57.4 mL/min
    lvc <- log(6.40); label("Central volume of distribution V1 (L)")                       # Supp. Table 1: V1 6.40 L [5.81-7.19]
    lq <- log(3.45); label("Intercompartmental clearance to peripheral1, Q1 (L/h)")        # Supp. Table 1: Q1 3.45 L/h [2.95-3.78]
    lvp <- log(5.43); label("Peripheral volume V2, slow peripheral (L)")                   # Supp. Table 1: V2 5.43 L [5.23-5.67]
    lq2 <- log(12.0); label("Intercompartmental clearance to peripheral2, Q2 (L/h)")       # Supp. Table 1: Q2 12.0 L/h [9.14-14.4]
    lvp2 <- log(3.26); label("Peripheral volume V3, fast peripheral (L)")                  # Supp. Table 1: V3 3.26 L [2.68-3.83]

    # ---------------------------------------------------------------------
    # IIV: Supp. Table 1 'CV%' column, converted with the paper's own
    # formula CV% = sqrt(exp(omega^2) - 1) (Methods), i.e.
    # omega^2 = log(1 + CV^2). The final model estimated a FULL omega
    # block (Results), but no covariance or correlation is printed, so the
    # etas are encoded as uncorrelated (see the vignette Assumptions).
    # ---------------------------------------------------------------------
    etalcl ~ 0.10336   # Supp. Table 1 CL CV 33%: log(1 + 0.33^2)
    etalvc ~ 0.14842   # Supp. Table 1 V1 CV 40%: log(1 + 0.40^2)
    etalq ~ 0.15537    # Supp. Table 1 Q1 CV 41%: log(1 + 0.41^2)
    etalvp ~ 0.07036   # Supp. Table 1 V2 CV 27%: log(1 + 0.27^2)
    etalq2 ~ 0.23119   # Supp. Table 1 Q2 CV 51%: log(1 + 0.51^2)
    etalvp2 ~ 0.33425  # Supp. Table 1 V3 CV 63%: log(1 + 0.63^2)

    # ---------------------------------------------------------------------
    # Residual error: Supp. Table 1 prints 'prop. error (*10-3) 0.468
    # [0.302-0.660]' and 'add. error [ng/mL] 2.22 [0.79-3.98]'. Both are
    # read as NONMEM SIGMA variances on the mg/L concentration scale:
    # - the proportional entry is only plausible as a variance (SD
    #   sqrt(0.000468) = 2.16%); read as an SD it would be a 0.05% error,
    #   below any HPLC assay's precision;
    # - the additive entry sits in the same column and the analysis
    #   concentrations are in mg/L (Figures 1 and 5, 30-450 mg/L); an
    #   additive error of a few ng/mL would be unidentifiable against
    #   concentrations of 10^4-10^5 ng/mL, yet its bootstrap CI excludes
    #   zero, so the printed '[ng/mL]' is taken as the mg/L scale.
    # sqrt(2.22) = 1.490 mg/L. The vignette documents the competing reading.
    # ---------------------------------------------------------------------
    propSd <- 0.021633; label("Proportional residual error (fraction)")   # Supp. Table 1: prop. error 0.468e-3 (variance); sqrt = 0.021633
    addSd <- 1.490; label("Additive residual error (mg/L)")               # Supp. Table 1: add. error 2.22 (variance); sqrt = 1.490
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
