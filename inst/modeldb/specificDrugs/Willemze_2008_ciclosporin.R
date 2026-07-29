Willemze_2008_ciclosporin <- function() {
  description <- "Two-compartment population PK model for ciclosporin in children (aged 1.8-16.1 years) after allogeneic haematopoietic stem cell transplantation (Willemze 2008). First-order absorption with lag time and partial bioavailability for oral Neoral microemulsion; intravenous Sandimmune is given as a 2-hour infusion to the central compartment. The 'alternative parameterization' (CL, Q, Vp, plus Ka, Vc, Tlag, F) reported in Table 2 is used directly because it is the more physiologically interpretable set. IIV on Vc was fixed to zero; IIVs on Ka, CL, Q, Vp, Tlag, and F are estimated. Residual error is proportional. No covariate (body weight, length, age, or estimated GFR) was retained in the final model; those covariates are documented in covariatesDataExcluded."
  reference <- "Willemze AJ, Cremers SC, Schoemaker RC, Lankester AC, den Hartigh J, Burggraaf J, Vossen JM. Ciclosporin kinetics in children after stem cell transplantation. Br J Clin Pharmacol. 2008;66(4):539-545. doi:10.1111/j.1365-2125.2008.03217.x"
  vignette <- "Willemze_2008_ciclosporin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the empirical Bayes regression vs CL and V (Willemze 2008 Figure 2; Results 'Ciclosporin pharmacokinetics'). No correlation between clearance (or distribution volume) and body weight was observed (Pearson r = -0.10, Spearman r = -0.00, n = 17). The authors explicitly conclude that 'dosing ciclosporin per kg bodyweight is not supported by the results of this study' (Conclusion). Cohort range 10 to ~60 kg, mean ~35 kg."
    ),
    HT = list(
      description = "Body length (height)",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the empirical Bayes regression vs CL and V (Willemze 2008 Results 'Ciclosporin pharmacokinetics'). No correlation found; not retained."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the empirical Bayes regression vs CL and V (Willemze 2008 Results 'Ciclosporin pharmacokinetics'). Cohort range 1.8-16.1 years. No correlation found; not retained."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (Schwartz and Cockcroft-Gault formulas)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Screened in the empirical Bayes regression vs CL and V (Willemze 2008 Methods 'Pharmacokinetics and statistical analysis'). No correlation between renal function and CL/V was observed; not retained. Patients with severe renal dysfunction (serum creatinine > 2x ULN) were excluded by protocol."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 17L,
    n_studies      = 1L,
    age_range      = "1.8-16.1 years",
    weight_range   = ">10 kg (cohort mean approximately 35 kg per Discussion paragraph 4)",
    sex_female_pct = NA_real_,
    disease_state  = "Paediatric allogeneic haematopoietic stem cell transplantation (SCT); ciclosporin given as GVHD prophylaxis. Exclusions: severe renal dysfunction (serum creatinine > 2x ULN), severe liver dysfunction (bilirubin > 50 umol/L) or coagulopathy, veno-occlusive disease, respiratory insufficiency, haemodynamic instability.",
    dose_range     = "Intravenous Sandimmune 2 mg/kg/day in two short 2-h infusions starting the day before graft infusion, then oral Neoral microemulsion twice daily (every 12 h) with the daily dose tripled to compensate for oral bioavailability; doses adjusted to keep trough between 50 and 200 ug/L.",
    regions        = "Netherlands (single centre, Leiden University Medical Centre)",
    study_period   = "Prospective enrolment January 2002 to October 2005",
    formulations   = "Sandimmune (intravenous) and Neoral microemulsion (oral); Novartis Pharma, Basel, Switzerland",
    notes          = "Patients enrolled in the paediatric SCT unit at Leiden University Medical Centre. Inclusion required body weight > 10 kg and at least 2 days of ciclosporin treatment before sampling. Whole-blood ciclosporin was measured by fluorescence polarization immunoassay (Abbott AxSYM). Trough was sampled before each dose; additional samples were taken at multiple post-dose time points (see Figure 1 of the source). Baseline demographics in Table 1 of the source."
  )

  ini({
    # Structural parameters at the population reference subject. Values are
    # the 'alternative parameterization' point estimates from Willemze 2008
    # Table 2 (CL, Q, Vp) plus the primary-parameterization values for Ka,
    # Vc, t_lag, and F which are identical across the two parameterizations.
    lka     <- log(0.831); label("First-order oral absorption rate constant Ka (1/h)")      # Table 2: Ka = 0.831 1/h (SEM 0.156); IIV 38 %CV
    lcl     <- log(11.3);  label("Apparent clearance CL (L/h)")                              # Table 2 alternative parameterization: CL = 11.3 L/h (SEM 1.74); IIV 36 %CV
    lvc     <- log(16.5);  label("Central volume of distribution Vc (L)")                    # Table 2: Vc = 16.5 L (SEM 4.72); IIV fixed at 0 %CV
    lq      <- log(12.9);  label("Inter-compartmental clearance Q (L/h)")                    # Table 2 alternative parameterization: Q = 12.9 L/h (SEM 2.81); IIV 52 %CV
    lvp     <- log(59.9);  label("Peripheral volume of distribution Vp (L)")                 # Table 2 alternative parameterization: Vp = 59.9 L (SEM 9.00); IIV 19 %CV
    ltlag   <- log(0.638); label("Oral absorption lag time t_lag (h)")                       # Table 2: t_lag = 0.638 h (SEM 0.0912); IIV 37 %CV
    lfdepot <- log(0.386); label("Oral bioavailability F (fraction)")                        # Table 2: F = 0.386 (SEM 0.0787); IIV 28 %CV

    # IIV: omega^2 = log(1 + CV^2). Variance values shown to 5 decimal places.
    etalka      ~ 0.13491  # Ka IIV 38 %CV; log(1 + 0.38^2) = 0.13491
    etalcl      ~ 0.12182  # CL IIV 36 %CV (alternative parameterization); log(1 + 0.36^2) = 0.12182
    etalq       ~ 0.23931  # Q  IIV 52 %CV (alternative parameterization); log(1 + 0.52^2) = 0.23931
    etalvp      ~ 0.03551  # Vp IIV 19 %CV (alternative parameterization); log(1 + 0.19^2) = 0.03551
    etaltlag    ~ 0.12833  # t_lag IIV 37 %CV; log(1 + 0.37^2) = 0.12833
    etalfdepot  ~ 0.07549  # F     IIV 28 %CV; log(1 + 0.28^2) = 0.07549
    # IIV on Vc was fixed to 0 (Table 2 footnote 'IIV (%) 0% (fixed)'); no etalvc.

    # Proportional residual error (Methods 'Pharmacokinetics and statistical analysis').
    propSd <- 0.193; label("Proportional residual error (fraction)")                          # Table 2: residual variability 19.3 %CV
  })
  model({
    # Individual PK parameters. No covariate effects were retained in the
    # final model (Results 'Ciclosporin pharmacokinetics'; Figure 2; Conclusion).
    ka         <- exp(lka + etalka)
    cl         <- exp(lcl + etalcl)
    vc         <- exp(lvc)
    q          <- exp(lq  + etalq)
    vp         <- exp(lvp + etalvp)
    tlag_depot <- exp(ltlag + etaltlag)
    fdepot     <- exp(lfdepot + etalfdepot)

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Oral absorption with lag and first-order Ka (depot); IV infusion routes
    # the dose directly to central via the data column cmt = central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag_depot

    # Whole-blood ciclosporin (FPIA). Dose in mg, Vc in L gives mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
