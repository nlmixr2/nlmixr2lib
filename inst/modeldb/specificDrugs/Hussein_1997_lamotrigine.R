Hussein_1997_lamotrigine <- function() {
  description <- "One-compartment population PK model for oral lamotrigine monotherapy in adults and adolescents newly diagnosed with epilepsy; apparent oral clearance carries a first-order auto-induction term that decays toward a steady-state value over treatment duration and a multiplicative race effect for Asians vs Caucasians; apparent volume of distribution and absorption rate constant are time-invariant with no covariate effects retained in the final model (Hussein 1997)."
  reference <- "Hussein Z, Posner J. Population pharmacokinetics of lamotrigine monotherapy in patients with epilepsy: retrospective analysis of routine monitoring data. Br J Clin Pharmacol. 1997;43(5):457-465. doi:10.1046/j.1365-2125.1997.00594.x"
  vignette <- "Hussein_1997_lamotrigine"
  paper_specific_etas <- c("etalcl")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian)",
      notes              = "Time-fixed per subject. Multiplicative effect on apparent oral clearance: CL_o (Asian) = CL_o (Caucasian) x (1 + theta_ASIANS_CL). Hussein 1997 Table 4 reports theta_ASIANS_CL = -0.287 (5 Asian / 158 Caucasian patients).",
      source_name        = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened univariately on CL_o (Table 1; continuous and power forms, DOBJF 8.9 / 9.0) and V/F (Table 1; DOBJF 4.7 / 4.2 / 3.7); only the power form on CL_o reached statistical significance, but the effect did not survive multivariate elimination once duration of therapy and Asians were in the model (Table 3 power-of-WT on CL_o: 0.231, DOBJF 6.55, not retained). Not in the final model."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened univariately on CL_o (Table 2 DOBJF 4.08) and V/F (Table 2 DOBJF 0.10); neither reached statistical significance and AGE is not in the final model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened univariately on CL_o (Table 2 DOBJF 8.48, P<0.005 univariate) and V/F (Table 2 DOBJF 5.34); the female effect on CL_o did not survive multivariate elimination once duration of therapy and Asians were in the model (Table 3 females-on-CL_o: -0.088, DOBJF 7.67, not retained). Not in the final model."
    ),
    OC = list(
      description = "Concomitant oral-contraceptive use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened univariately on CL_o (Table 2 DOBJF 2.95); did not reach statistical significance and is not in the final model."
    ),
    DOSE = list(
      description = "Total daily lamotrigine dose",
      units       = "mg",
      type        = "continuous",
      notes       = "Screened univariately as an exponential decay form on CL_o (Table 2 DOBJF 7.35, theta_DOSE1 = 5.39, theta_DOSE2 = 0.0014); did not reach the 2-DF significance threshold (10.6) and is not in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 163L,
    n_studies      = 3L,
    age_range      = "14-76 years",
    weight_range   = "40.5-106.5 kg",
    weight_typical = "approximately 68 kg (population average per Discussion)",
    sex_female_pct = 50.3,
    race_ethnicity = c(Caucasian = 96.9, Asian = 3.1),
    disease_state  = "Newly diagnosed epilepsy (partial or generalised tonic-clonic seizures); not previously treated with antiepileptic drugs; normal hepatic and renal function.",
    dose_range     = "50-200 mg/day oral lamotrigine monotherapy; dosing frequency once, twice, or three times daily (24, 12, or 8 h interval) over 48 weeks; common dose-escalation 50 mg QD week 1 -> 50 mg BID week 2 -> 50 mg AM + 100 mg PM weeks 3-48, with investigator-led titration on efficacy/safety.",
    regions        = "Three multicentre Phase II/III trials sponsored by GlaxoWellcome (UK); regions not stated in the paper.",
    notes          = "Population per Hussein 1997 Results 'Demographic characteristics' (page 459) and Figure 1: 158 Caucasians + 5 Asians; 81 males + 82 females; 20 Caucasian females received concomitant oral contraceptives. Concentration data are retrospective routine monitoring samples from three Phase II/III trials. Dosing times were not captured in the case report forms except for the last dose pre-visit, so dosing intervals of 24/12/8 h were assumed for once/twice/thrice daily regimens (Methods 'Database construction')."
  )

  ini({
    # Structural PK parameters (Hussein 1997 Table 4, p. 464).
    # The published random-effects model uses an additive-on-fraction
    # form (CL_o = TVCL x (1 + eta_CL); V/F = TV(V/F) x (1 + eta_V)).
    # We translate to the canonical nlmixr2 log-normal form
    # (CL_o = TVCL x exp(eta_CL); V/F = TV(V/F) x exp(eta_V)) via
    # omega^2 = log(CV^2 + 1); for the published CV magnitudes
    # (~32 % / 34 %) the two parameterisations are within ~5 %.
    # See vignette Assumptions and deviations.

    # Asymptotic (post-induction) typical CL_o for a Caucasian patient.
    lcl_ss <- log(2.28);    label("Apparent oral clearance after full induction CL_o,ss (L/h)")              # Table 4: theta_CL = 2.28 L/h (95% CI 2.14-2.41)

    # Magnitude of the time-decaying auto-induction offset, in the same
    # units as CL_o (L/h): cl_o (t=0) = cl_ss - cl_time.
    lcl_time <- log(0.338); label("Auto-induction offset on CL_o (CL_o,ss - CL_o,t=0) (L/h)")               # Table 4: theta_CL^TIME1 = 0.338 L/h (95% CI 0.191-0.485)

    # First-order rate constant for the auto-induction decay; defines how
    # fast cl_o approaches cl_ss after starting therapy. Reported in
    # 1/week; converted to 1/h inside model() via division by 168.
    lkdes <- log(0.119);    label("Auto-induction first-order rate constant kdes (1/week)")                  # Table 4: theta_CL^TIME2 = 0.119 /wk (95% CI 0.035-0.203)

    # Apparent volume of distribution and absorption rate constant.
    lvc <- log(77.4);       label("Apparent volume of distribution V/F (L)")                                  # Table 4: theta_V = 77.4 L (95% CI 67.4-87.4)
    lka <- log(3.18);       label("First-order absorption rate constant Ka (1/h)")                            # Table 4: theta_Ka = 3.18 /h (95% CI 2.31-4.05)

    # Multiplicative Asian-race effect on the composite CL_o.
    # Negative value indicates lower CL_o in Asians: CL_o (Asian) =
    # CL_o (Caucasian) x (1 + (-0.287)) = 0.713 x CL_o (Caucasian).
    e_race_asian_cl <- -0.287; label("Race effect on CL_o: Asian vs Caucasian (fraction; applied as (1 + e * RACE_ASIAN))")  # Table 4: theta_CL^ASIANS = -0.287 (95% CI -0.410 to -0.164)

    # Inter-individual variability. Hussein 1997 reports CV(CL_o) =
    # 32.1 % and CV(V/F) = 33.6 % as Table 4 omega^2 entries. The
    # paper applies a single eta per parameter, with eta_CL multiplying
    # the composite TVCL (i.e. the same eta scales cl_ss and cl_time
    # in lockstep); we encode this with the paper-specific eta etalcl,
    # which is applied to BOTH lcl_ss and lcl_time inside model().
    # No IIV is reported for Ka or for kdes; Hussein 1997 Discussion
    # notes that the magnitude of variability in absorption was
    # "inestimable" because few plasma samples were collected during
    # the absorption phase.
    etalcl ~ 0.0981  # Table 4: CV(CL_o) = 32.1% (95% CI 26.6-36.8); omega^2 = log(1 + 0.321^2) = 0.0981
    etalvc ~ 0.1070  # Table 4: CV(V/F)  = 33.6% (95% CI 13.0-45.7); omega^2 = log(1 + 0.336^2) = 0.1070

    # Residual error. Hussein 1997 Table 4 reports a proportional
    # residual model (Cobs = Cpred x (1 + epsilon)) with CV = 20.8%.
    propSd <- 0.208;        label("Proportional residual error (fraction)")                                   # Table 4: sigma^2 row CV = 20.8% (95% CI 18.8-22.6)
  })

  model({
    # Convert the rxode2 simulation time (h) to weeks since start of
    # therapy; the auto-induction term is parameterised in 1/week.
    time_wk <- time / 168

    # Individual structural parameters. The single composite eta etalcl
    # multiplies both the asymptotic and the time-varying CL arms in
    # lockstep, reproducing the paper's CL_o = TVCL(t) x exp(eta).
    cl_ss   <- exp(lcl_ss   + etalcl)
    cl_time <- exp(lcl_time + etalcl)
    kdes    <- exp(lkdes)
    vc      <- exp(lvc + etalvc)
    ka      <- exp(lka)

    # Composite apparent oral clearance: asymptotic CL minus the
    # time-decaying auto-induction offset, scaled multiplicatively by
    # the Asian race factor. At time = 0: CL_o = (cl_ss - cl_time) x
    # (1 + e * RACE_ASIAN). As time -> infinity: CL_o -> cl_ss x
    # (1 + e * RACE_ASIAN).
    cl <- (cl_ss - cl_time * exp(-kdes * time_wk)) * (1 + e_race_asian_cl * RACE_ASIAN)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
