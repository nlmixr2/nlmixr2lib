Bustinduy_2020_praziquantel <- function() {
  description <- paste(
    "Two-compartment oral population PK model with a first-order absorption (gut) compartment, an",
    "absorption lag and a reversible breast-milk compartment for praziquantel (racemic, total PZQ)",
    "in 45 Filipino women with Schistosoma japonicum infection -- 15 in early pregnancy (12-16",
    "weeks gestation), 15 in late pregnancy (30-36 weeks) and 15 lactating postpartum women (5-7",
    "months postpartum) -- given 60 mg/kg as two 30 mg/kg oral doses 3 h apart. Plasma and breast",
    "milk were co-modelled non-parametrically with NPAG in Pmetrics. Drug moves from central to",
    "milk and back with first-order rate constants and is NOT eliminated through milk (the authors'",
    "identifiability choice), so the milk compartment behaves as a sampled second peripheral",
    "compartment with its own apparent volume. Clearance and volumes are apparent (CL/F, Vc/F,",
    "Vmilk/F); bioavailability was not estimated. No covariate was retained: weight showed no",
    "relationship with CL/F or Vc/F, and the higher CL/F of early-pregnancy women was reported",
    "only as a post hoc comparison of Bayesian posteriors, not built into the model. Residual",
    "variability is fixed(0) because the Pmetrics assay-error model was not published.",
    sep = " "
  )
  reference <- paste(
    "Bustinduy AL, Kolamunnage-Dona R, Mirochnick MH, Capparelli EV, Tallo V, Acosta LP,",
    "Olveda RM, Friedman JF, Hope WW. Population pharmacokinetics of praziquantel in pregnant and",
    "lactating Filipino women infected with Schistosoma japonicum. Antimicrob Agents Chemother.",
    "2020;64(9):e00566-20. doi:10.1128/AAC.00566-20. PMCID: PMC7449211.",
    sep = " "
  )
  vignette <- "Bustinduy_2020_praziquantel"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "praziquantel", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "praziquantel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "praziquantel", units = "mg", specimen = "plasma", verified = TRUE),
    milk = list(analyte = "praziquantel", units = "mg", specimen = "milk", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened but NOT retained. Results: 'There was no relationship between weight and",
        "Bayesian estimates for the apparent clearance ... or between weight and the apparent",
        "volume of the central compartment' (r = 0.0303, P = 0.8435 and r = 0.1617, P = 0.2885;",
        "Fig. 5); 'Hence, covariates were not incorporated into the structural model.' Weight",
        "still sets the administered dose, which is prescribed as 2 x 30 mg/kg.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    age_range = "18-44 years",
    age_median = "24.0 years (mean 25.5, SD 6.39; 47 enrolled)",
    weight_range = "approximately 36-63 kg (read from Fig. 5; not tabulated)",
    weight_median = "47.9 kg (median used for the Monte Carlo simulations); mean 48.5 kg (SD 7.69; 47 enrolled)",
    sex_female_pct = 100,
    race_ethnicity = "Asian (Filipino), 100%",
    disease_state = paste(
      "Stool-positive Schistosoma japonicum infection, otherwise healthy; infection intensity",
      "low (< 100 eggs per gram) in 46 of 47 enrolled and moderate in 1.",
      sep = " "
    ),
    dose_range = paste(
      "60 mg/kg praziquantel given orally as two 30 mg/kg doses approximately 3 h apart, after a",
      "carbohydrate-rich snack. Absolute total doses about 2100-3800 mg (Fig. 6A).",
      sep = " "
    ),
    regions = "Philippines (northeastern Leyte)",
    reproductive_status = paste(
      "Early pregnancy 12-16 weeks gestation (n = 15 analysed; 17 enrolled, 2 vomited and were",
      "not sampled), late pregnancy 30-36 weeks gestation (n = 15), lactating 5-7 months",
      "postpartum (n = 15).",
      sep = " "
    ),
    notes = paste(
      "Baseline demographics from Table 1 (47 enrolled). Plasma sampled pre-dose and at 1, 2, 3",
      "(before the second dose), 4, 5, 6, 7, 8, 9, 12, 15 and 24 h after the first dose; breast",
      "milk hand-expressed at 3, 6, 9, 12, 15 and 24 h in the lactating group only. LC-MS assay;",
      "LLOQ 31.3 ng/mL in plasma and 4.3 ng/mL in milk (Methods). Posterior CL/F was higher in",
      "early pregnancy (median about 425 L/h vs about 245 L/h; Fig. 6C) and AUC0-24 lower (median",
      "about 7 vs 11-12 mg*h/L; Fig. 6D), but this was not encoded as a covariate effect.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # All structural values are the MEAN of the NPAG non-parametric parameter
    # distribution in Table 2 (45 women, plasma + breast milk co-modelled);
    # the Table 2 median is recorded in each trailing comment. The mean is used
    # because it is the parameter vector the paper's population Monte Carlo
    # simulation (Fig. 7) reproduces -- the median vector under-predicts the
    # Fig. 7A centiles about two-fold because its Ka is 5-fold lower -- and for
    # consistency with the sibling extraction Bustinduy_2016_praziquantel. See
    # the vignette for the mean-vs-median evidence.
    #
    # Bioavailability was not estimated, so CL, Vc and Vmilk are apparent.
    # ------------------------------------------------------------------------
    lka <- log(2.012); label("Absorption rate constant from the gut to the central compartment (1/h)")
    # Table 2, 'Ka (h-1)' mean = 2.012 (median 0.395, SD 4.301, CV 213.750%)
    lcl <- log(324.075); label("Apparent clearance from the central compartment, SCL/F (L/h)")
    # Table 2, 'SCL/F (liter/h)' mean = 324.075 (median 277.447, SD 175.373, CV 54.115%)
    lvc <- log(183.006); label("Apparent central volume of distribution, Vc/F (L)")
    # Table 2, 'Vc/F (liter)' mean = 183.006 (median 142.618, SD 93.211, CV 50.933%)

    lk12 <- log(19.313); label("Transfer rate constant central -> peripheral1, Kcp (1/h)")
    # Table 2, 'Kcp (h-1)' mean = 19.313 (median 18.941, SD 10.167, CV 52.644%)
    lk21 <- log(15.816); label("Transfer rate constant peripheral1 -> central, Kpc (1/h)")
    # Table 2, 'Kpc (h-1)' mean = 15.816 (median 13.996, SD 9.447, CV 59.733%)

    lk_central_milk <- log(18.750); label("Transfer rate constant central -> breast milk, Kcb (1/h)")
    # Table 2, 'Kcb (h-1)' mean = 18.750 (median 19.301, SD 9.387, CV 50.067%)
    lk_milk_central <- log(17.816); label("Transfer rate constant breast milk -> central, Kbc (1/h)")
    # Table 2, 'Kbc (h-1)' mean = 17.816 (median 17.077, SD 7.845, CV 44.031%)
    lvmilk <- log(612.130); label("Apparent volume of the breast-milk compartment, Vb/F (L)")
    # Table 2, 'Vb/F (liter)' mean = 612.130 (median 563.802, SD 395.661, CV 64.637%)

    ltlag <- log(0.772); label("Absorption lag time (h)")
    # Table 2, 'Lag (h)' mean = 0.772 (median 0.868, SD 0.233, CV 30.202%)

    lfdepot <- fixed(log(1)); label("Oral bioavailability of the depot (unitless)")
    # Not estimated: Table 2 reports only apparent (/F) clearance and volumes, and the Results
    # attribute the early-pregnancy exposure difference to CL and/or F without separating them.

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution, not a parametric omega. Table 2 CV% equals SD / mean on the
    # linear scale (e.g. 175.373 / 324.075 = 54.115%), so each CV% is carried
    # as a LOG-NORMAL approximation, omega^2 = log(CV^2 + 1). Covariances are
    # not reported (the Results give only the CL/F-V/F posterior correlation,
    # r = 0.636) and are therefore absent.
    # ------------------------------------------------------------------------
    etalka ~ 1.717199 # Table 2 Ka CV 213.750% -> log(2.1375^2 + 1)
    etalcl ~ 0.256844 # Table 2 SCL/F CV 54.115% -> log(0.54115^2 + 1)
    etalvc ~ 0.230649 # Table 2 Vc/F CV 50.933% -> log(0.50933^2 + 1)
    etalk12 ~ 0.244622 # Table 2 Kcp CV 52.644% -> log(0.52644^2 + 1)
    etalk21 ~ 0.305131 # Table 2 Kpc CV 59.733% -> log(0.59733^2 + 1)
    etalk_central_milk ~ 0.223680 # Table 2 Kcb CV 50.067% -> log(0.50067^2 + 1)
    etalk_milk_central ~ 0.177203 # Table 2 Kbc CV 44.031% -> log(0.44031^2 + 1)
    etalvmilk ~ 0.349102 # Table 2 Vb/F CV 64.637% -> log(0.64637^2 + 1)
    etaltlag ~ 0.087293 # Table 2 Lag CV 30.202% -> log(0.30202^2 + 1)

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported: neither the Pmetrics
    # assay-error polynomial nor a lambda/gamma term appears in the paper, which
    # reports only weighted-residual diagnostics (Fig. 4). Carried as fixed(0)
    # for both outputs rather than invented -- see vignette Errata.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD, plasma (fraction; 0 -- not reported in the source)")
    addSd <- fixed(0); label("Additive residual SD, plasma (mg/L; 0 -- not reported in the source)")
    propSd_Cmilk <- fixed(0); label("Proportional residual SD, breast milk (fraction; 0 -- not reported in the source)")
    addSd_Cmilk <- fixed(0); label("Additive residual SD, breast milk (mg/L; 0 -- not reported in the source)")
  })

  model({
    # 1. Individual parameters. No covariate was retained (see
    #    covariatesDataExcluded): the base model is the final model.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    k_central_milk <- exp(lk_central_milk + etalk_central_milk)
    k_milk_central <- exp(lk_milk_central + etalk_milk_central)
    vmilk <- exp(lvmilk + etalvmilk)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    # 2. Methods equations (1)-(4), with the three typesetting errors of the
    #    printed equations corrected so that mass is conserved (vignette Errata):
    #    eq. (3) prints '- Kcp x X(3)' for the peripheral return (Kpc intended),
    #    and eq. (4) prints '- Kcb * X(2)' for the milk inflow (+ intended, the
    #    mirror of the '- Kcb * X(2)' loss in eq. (2)). Breast milk has no
    #    terminal elimination: 'PZQ was allowed to redistribute back to the
    #    maternal plasma without terminal elimination via expression of breast
    #    milk.'
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k_central_milk * central + k_milk_central * milk
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(milk) <- k_central_milk * central - k_milk_central * milk

    # 3. Absorption lag ('A lag function ... was applied between the oral
    #    administration of PZQ and the appearance of drug in the central
    #    compartment') and the unit bioavailability anchor.
    alag(depot) <- tlag
    f(depot) <- exp(lfdepot)

    # 4. Outputs. Printed Y(1) = X(1) / Vc divides the GUT amount; the plasma
    #    output is the central amount X(2) / Vc (vignette Errata).
    #    Y(2) = X(4) / Vb is the breast-milk concentration.
    Cc <- central / vc
    Cmilk <- milk / vmilk

    Cc ~ add(addSd) + prop(propSd)
    Cmilk ~ add(addSd_Cmilk) + prop(propSd_Cmilk)
  })
}
