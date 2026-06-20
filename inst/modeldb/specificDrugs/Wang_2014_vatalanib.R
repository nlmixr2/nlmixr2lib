Wang_2014_vatalanib <- function() {
  description <- "One-compartment population PK model for oral vatalanib in adults with myelodysplastic syndrome (CALGB 10105); apparent oral clearance carries a first-order auto-induction term that rises from a pre-induction value toward a steady-state post-induction value over the first 7 days of therapy, lagged first-order absorption, log-normal residual error on the natural-log-transformed concentration; no covariates retained in the final model (Wang 2014)."
  reference <- "Wang X, Owzar K, Gupta P, Larson RA, Mulkey F, Miller AA, Lewis LD, Hurd D, Vij R, Ratain MJ, Murry DJ; for the Alliance for Clinical Trials in Oncology. Vatalanib population pharmacokinetics in patients with myelodysplastic syndrome: CALGB 10105 (Alliance). Br J Clin Pharmacol. 2014;78(5):1005-1013. doi:10.1111/bcp.12427"
  vignette <- "Wang_2014_vatalanib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Actual body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL_induced/F and on Vd/F; not retained (Wang 2014 Results 'Model development': 'none of the covariate effects was significant after the stepwise covariate model-building process'). Cohort median 80 kg (range 48-128)."
    ),
    IBW = list(
      description = "Ideal body weight (calculated using a standard equation)",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL_induced/F; not retained. Wang 2014 Methods 'Covariate model development': 'Ideal bodyweight and dosing weight were calculated using standard equations.'"
    ),
    DOSEWT = list(
      description = "Dosing weight (calculated using a standard equation)",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL_induced/F; not retained. Wang 2014 Methods 'Covariate model development'."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Patient demographic characteristic assessed during covariate screening (Wang 2014 Methods 'Covariate model development'); cohort median 170 cm (range 149-193). Not retained."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Patient demographic characteristic assessed during covariate screening. Highly correlated with body weight; selected against based on physiological plausibility per Wang 2014 Methods. Cohort median 1.91 m^2 (range 1.46-2.46). Not retained."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL_induced/F; not retained. 86 males / 51 females (62.8% / 37.2%) in the analysis cohort."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL_induced/F; not retained. Cohort median 70 years (range 20-91). Wang 2014 Discussion notes that the relatively narrow age range likely limited the power to detect a covariate effect."
    ),
    RACE = list(
      description = "Race (Caucasian vs. Other)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Patient demographic characteristic assessed during covariate screening; not retained. 128 Caucasian / 9 Other (93.4% / 6.6%) in the analysis cohort."
    ),
    TBIL = list(
      description = "Total bilirubin (liver function test)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on CL_induced/F; not retained. Cohort median 0.7 mg/dL (range 0.2-2.0)."
    ),
    AST = list(
      description = "Aspartate aminotransferase (liver function test)",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened on CL_induced/F; not retained. Cohort median 23 IU/L (range 7-92)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 137L,
    n_studies      = 1L,
    age_range      = "20-91 years",
    age_median     = "70 years",
    weight_range   = "48-128 kg",
    weight_median  = "80 kg",
    sex_female_pct = 37.2,
    race_ethnicity = c(Caucasian = 93.4, Other = 6.6),
    disease_state  = "Adult patients with primary or therapy-related (secondary) myelodysplastic syndrome (MDS) enrolled in the Cancer and Leukemia Group B (CALGB) 10105 study (Alliance), a multicentre open-label phase II study.",
    dose_range     = "Oral vatalanib 750 or 1250 mg once daily in 28-day cycles; the protocol was amended after enrolment began to allow a 750 mg starting dose.",
    regions        = "United States (multicentre Alliance / CALGB sites).",
    notes          = "564 vatalanib plasma concentration measurements across 137 patients (Wang 2014 Results 'Database description', Table 1). Sampling clusters: 66.5% within 24 h of the first dose (predose, 15-45 min, 1-3 h, 4-6 h), 26.6% between days 7 and 14 (two samples at least 1 h apart), 6.9% on day 28 (pre-cycle-2 trough). Patients dosed on an empty stomach or at least 30 min after breakfast; grapefruit and grapefruit juice avoided. NONMEM 7.2 with FO estimation (FOCE failed). Final model: 1-compartment with lagged first-order absorption and first-order auto-induction of oral clearance."
  )

  ini({
    # Structural PK parameters (Wang 2014 Table 2, p. 1009).
    # Final-model parameterisation per Wang 2014 Results 'Model development':
    #   the base model with no retained covariates IS the final model.
    # The autoinduction structure:
    #   CL/F(t) = CL_induced/F - delta_CL/F * exp(-K_induct * t)
    # so that CL/F(t=0) = CL_initial/F = CL_induced/F - delta_CL/F
    # and    CL/F(t -> Inf) = CL_induced/F (the post-induction steady state).
    # Implemented with the canonical CL-component names lcl_ss
    # (post-induction asymptote = CL_induced/F) and lcl_time
    # (auto-induction offset = delta_CL/F) following the
    # Hussein_1997_lamotrigine precedent.

    # Apparent post-induction (steady-state) oral clearance CL_induced/F.
    lcl_ss   <- log(54.9);  label("Apparent oral clearance after full induction CL_induced/F (L/h)")  # Wang 2014 Table 2: CL_induced/F = 54.9 L/h (95% CI 42.5-67.2; 11.7% RSE)

    # Auto-induction offset: difference between post-induction steady-state
    # CL and pre-induction (t=0) CL, in the same units as CL/F (L/h).
    lcl_time <- log(30.1);  label("Auto-induction offset on CL/F (CL_induced/F - CL_initial/F) (L/h)")  # Wang 2014 Table 2: delta_CL/F = 30.1 L/h (95% CI 14.8-45.4; 27.6% RSE)

    # First-order auto-induction rate constant. Wang 2014 fixed this
    # parameter at 0.023 /h because samples were collected only on day 1
    # and after day 7, leaving insufficient information to estimate it
    # accurately; the value was chosen so that maximal induction is reached
    # on day 7 (Wang 2014 Table 2 footnote * and Discussion p. 1011).
    lkdes    <- fixed(log(0.023)); label("Auto-induction first-order rate constant K_induct (1/h)")    # Wang 2014 Table 2: K_induct = 0.023 /h (fixed; footnote * 'maximal clearance induction was reached on day 7')

    # Apparent volume of distribution V/F. One-compartment model.
    lvc      <- log(53.8);  label("Apparent volume of distribution Vd/F (L)")  # Wang 2014 Table 2: Vd/F = 53.8 L (95% CI 38.4-69.1; 14.9% RSE)

    # Absorption rate constant for the lagged first-order absorption.
    lka      <- log(0.172); label("First-order absorption rate constant Ka (1/h)")  # Wang 2014 Table 2: Ka = 0.172 /h (95% CI 0.141-0.203; 9.2% RSE)

    # Absorption lag time A_lag (h). Wang 2014 Methods 'Base model
    # building' notes that the lagged first-order absorption model
    # reduced the objective function value by 88.7 compared with plain
    # first-order absorption.
    ltlag    <- log(0.178); label("Absorption lag time A_lag (h)")  # Wang 2014 Table 2: A_lag = 0.178 h (95% CI 0.136-0.220; 11.7% RSE)

    # Inter-individual variability. Wang 2014 Methods 'Base model building'
    # uses log-normal IIV: P_i = P_pop * exp(eta_i). Table 2 reports
    # variability as %CV (back-transformed from the omega^2 estimates).
    # Conversion: omega^2 = log(1 + CV^2).
    # The IIV on delta_CL/F (cl_time) was fixed to zero (Wang 2014 Results
    # 'Model development': "the IIV value associated with delta CL/F was
    # fixed to zero due to extremely small estimate and failure of
    # convergence of the covariance step"); cl_time therefore enters
    # the model without an eta term.
    # No IIV on K_induct: Wang 2014 fixed K_induct itself and does not
    # report a between-subject random effect on it.
    etalcl_ss ~ 0.0507  # Wang 2014 Table 2: %CV(CL_induced/F) = 22.8% (36.7% RSE); omega^2 = log(1 + 0.228^2) = 0.0507
    etalvc    ~ 0.5380  # Wang 2014 Table 2: %CV(Vd/F)        = 84.4% (20.3% RSE); omega^2 = log(1 + 0.844^2) = 0.5380
    etalka    ~ 0.1181  # Wang 2014 Table 2: %CV(Ka)          = 35.4% (24.1% RSE); omega^2 = log(1 + 0.354^2) = 0.1181
    etaltlag  ~ 1.0577  # Wang 2014 Table 2: %CV(A_lag)       = 137.1% (76.1% RSE); omega^2 = log(1 + 1.371^2) = 1.0577

    # Residual error. Wang 2014 Methods 'Base model building': "Residual
    # variability (RV) was modelled using an additive error model for
    # natural logarithm-transformed data as follows:
    # ln C_ij = ln C_pred,ij + epsilon_ij." This is the log-normal
    # residual structure, encoded in nlmixr2 as `Cc ~ lnorm(expSd)` with
    # expSd = log-domain SD = sqrt(sigma^2_additive). Wang 2014 Table 2
    # reports sigma^2_additive = 0.596 (17.8% RSE), so expSd = sqrt(0.596).
    expSd <- 0.772;  label("Log-normal residual error (log-scale SD)")  # Wang 2014 Table 2: sigma^2 additive = 0.596 (17.8% RSE); sqrt(0.596) = 0.7720
  })

  model({
    # Individual structural parameters. Wang 2014 puts a single eta on
    # CL_induced/F; the auto-induction offset delta_CL/F (cl_time)
    # has its IIV fixed to zero per Results 'Model development'.
    cl_ss   <- exp(lcl_ss + etalcl_ss)
    cl_time <- exp(lcl_time)
    kdes    <- exp(lkdes)
    vc      <- exp(lvc + etalvc)
    ka      <- exp(lka + etalka)
    tlag    <- exp(ltlag + etaltlag)

    # Composite time-dependent apparent oral clearance:
    #   CL/F(t) = CL_induced/F - delta_CL/F * exp(-K_induct * t)
    # At t = 0:        CL/F = cl_ss - cl_time = CL_initial/F (pre-induction).
    # At t -> infinity: CL/F = cl_ss (post-induction steady state).
    # `time` is the rxode2 simulation-time variable in hours since the start
    # of therapy (consistent with the dose schedule used in the cohort:
    # first dose on study day 1 hour 0).
    cl <- cl_ss - cl_time * exp(-kdes * time)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag time applied to the depot compartment per Wang 2014's
    # lagged first-order absorption model.
    alag(depot) <- tlag

    # Observation: vatalanib plasma concentration in ng/mL. Doses are in
    # mg and vc is in L, so central/vc has units mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
  })
}
