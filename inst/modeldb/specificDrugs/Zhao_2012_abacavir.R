Zhao_2012_abacavir <- function() {
  description <- "Two-compartment population PK model for oral abacavir in HIV-infected infants and toddlers (Zhao 2012) developed on the PENTA 15 crossover trial of 8 mg/kg twice-daily vs 16 mg/kg once-daily dosing; CL/F scales with body weight via an estimated power exponent (1.14) referenced to the population median weight of 12 kg, and inter-occasion variability on CL/F is multiplexed by the binary OCC indicator across the BID (occasion 1) and QD (occasion 2) study phases."
  reference <- "Zhao W, Cella M, Della Pasqua O, Burger D, Jacqz-Aigrain E, on behalf of Pediatric European Network for Treatment of AIDS (PENTA) 15 study group. Population pharmacokinetics and maximum a posteriori probability Bayesian estimator of abacavir: application of individualized therapy in HIV-infected infants and toddlers. Br J Clin Pharmacol. 2012;73(4):641-648. doi:10.1111/j.1365-2125.2011.04121.x"
  vignette <- "Zhao_2012_abacavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight on the day of pharmacokinetic sampling (Table 1 footnote). The covariate effect is an estimated power model on CL/F with reference weight 12 kg (population median): CL/F = 13.4 * (WT/12)^1.14. Cohort range 7.4-15.9 kg (mean 11.6 kg).",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing across the PENTA 15 crossover phases.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1 and 2 identify the two PENTA 15 crossover phases within subject: 1 = abacavir 8 mg/kg twice daily (BID, weeks 0-4), 2 = abacavir 16 mg/kg once daily (QD, weeks 4-8). Decomposed inside `model()` into binary indicators `oc1` and `oc2` that multiplex the two IOV etas on log-CL (NONMEM `$OMEGA BLOCK(1)` + `SAME` translation; both etas share the same variance per Zhao 2012 Table 3).",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_studies      = 1L,
    age_range      = "0.43-2.89 years (mean 1.8)",
    age_median     = "1.8 years (mean reported)",
    weight_range   = "7.4-15.9 kg (mean 11.6)",
    weight_median  = "12 kg",
    sex_female_pct = 47.8,
    disease_state  = "HIV type-1-infected infants and toddlers aged 3-36 months on antiretroviral therapy, enrolled in the PENTA 15 open-label crossover study comparing once- and twice-daily abacavir + lamivudine dosing.",
    dose_range     = "Oral abacavir 8 mg/kg BID (weeks 0-4) then 16 mg/kg QD (weeks 4-8) per the PENTA 15 crossover schedule; intensive PK sampling at steady state in each phase.",
    regions        = "France, Germany, Italy, Spain, United Kingdom",
    n_observations = 347L,
    notes          = "Baseline demographics from Zhao 2012 Table 1. 12 male and 11 female (47.8% female). Height (cm) mean 81 (range 62-98, n=22). Body mass index (kg/m^2) mean 17.8 (range 15.0-24.2, n=22). Serum creatinine (mg/L) mean 34.7 (range 22.1-53.9, n=21). 347 plasma abacavir concentrations were available for population modelling; 13.5% were below the lower limit of quantification (LLQ 0.015 mg/L) and were imputed as LLQ/2. 18 of 23 patients had a full PK profile; 5 had incomplete profiles. The covariate screen tested age, gender, weight, height, body mass index, serum creatinine, and drug administration frequency; only WT on CL/F was retained in the final model (Tables 2-3)."
  )

  ini({
    # Structural PK parameters (Zhao 2012 Table 3 "Final estimate" column).
    # Reference weight for the CL/F covariate term is the population median (12 kg).
    lka <- log(0.758); label("Absorption rate constant (ka, 1/h)")                                              # Table 3 Ka = 0.758 1/h
    lcl <- log(13.4);  label("Apparent oral clearance at 12 kg reference (CL/F, L/h)")                          # Table 3 theta_4 = 13.4 L/h
    lvc <- log(4.94);  label("Apparent central volume of distribution (V1/F, L)")                               # Table 3 V1/F = 4.94 L
    lvp <- log(8.12);  label("Apparent peripheral volume of distribution (V2/F, L)")                            # Table 3 V2/F = 8.12 L
    lq  <- log(1.25);  label("Apparent inter-compartmental clearance (Q/F, L/h)")                               # Table 3 Q/F = 1.25 L/h

    # Estimated power-form covariate effect on CL/F.
    # CL/F_j = 13.4 * (WT_j / 12)^e_wt_cl with reference WT = 12 kg (cohort median).
    e_wt_cl <- 1.14; label("Power exponent for body weight on CL/F (unitless)")                                 # Table 3 theta_5 = 1.14

    # Inter-individual variability (IIV). Paper reports an exponential model:
    # P_i = TV(P) * exp(eta_i), CV% summarized in Table 3. omega^2 = log(1 + CV^2)
    # for a log-normal random effect. No off-diagonal covariances are reported in
    # Table 3, so the IIV omegas are diagonal here.
    etalcl ~ 0.03884  # Table 3 IIV CL/F = 19.9% CV; omega^2 = log(1 + 0.199^2) = 0.038835
    etalvp ~ 0.14977  # Table 3 IIV V2/F = 40.2% CV; omega^2 = log(1 + 0.402^2) = 0.149770
    etalq  ~ 0.09120  # Table 3 IIV Q/F  = 30.9% CV; omega^2 = log(1 + 0.309^2) = 0.091197

    # Inter-occasion variability (IOV) on log-CL across the BID and QD crossover
    # phases (Zhao 2012 Methods: "Interoccasion variability on CL/F was coupled to
    # interindividual variability by an additive model"). Both occasions share the
    # same variance per Table 3 (NONMEM `$OMEGA BLOCK(1)` + `SAME` translation);
    # nlmixr2 has no `SAME` shortcut so each occasion gets its own eta, with the
    # second fixed equal to the first (matching the Jonsson 2011 ethambutol pattern).
    etaiov_cl_1 ~ 0.04559         # Table 3 IOV CL/F = 21.6% CV; omega^2 = log(1 + 0.216^2) = 0.045591
    etaiov_cl_2 ~ fix(0.04559)    # fixed equal to occasion-1 variance per Table 3 (`SAME` translation)

    # Residual error: paper Methods state "Residual variability was best described
    # by a proportional model"; Table 3 reports 14.1% (the residual proportional SD
    # on the linear concentration scale).
    propSd <- 0.141; label("Proportional residual error (fraction)")                                            # Table 3 residual proportional = 14.1%
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for IOV
    # multiplexing on log-CL across the PENTA 15 BID and QD crossover phases.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Individual PK parameters; estimated power-form weight effect on CL/F.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 12)^e_wt_cl
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma abacavir concentration: dose in mg, volumes in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
