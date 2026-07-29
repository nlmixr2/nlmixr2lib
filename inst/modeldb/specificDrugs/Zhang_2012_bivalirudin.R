Zhang_2012_bivalirudin <- function() {
  description <- "Population PK and PK-PD model for bivalirudin, a synthetic bivalent direct thrombin inhibitor, in young healthy Chinese volunteers (Zhang 2012). PK: two-compartment intravenous disposition with body-weight-normalised structural parameters (CL = 0.323 L/h/kg, V1 = 0.086 L/kg, Q = 0.0957 L/h/kg, V2 = 0.0554 L/kg); no covariates retained after a 30-covariate screen; log-normal IIV on CL, V1, and V2 with IIV on Q fixed to zero. PD: direct-response sigmoid Emax (Hill coefficient fixed at 1) on activated clotting time (ACT) using the central-compartment concentration as the effect site (E0 = 134 s, Emax = 318 s, EC50 = 2.44 mg/L); one covariate retained -- red blood cell count (RBC, 10^12 cells/L) on EC50 via the linear-deviation form EC50_i = theta_EC50 * exp(eta_EC50) * (1 + 1.70 * (RBC - 4.40)) centred at the cohort median 4.40."
  reference <- paste(
    "Zhang DM, Wang K, Zhao X, Li YF, Zheng QS, Wang ZN, Cui YM.",
    "Population pharmacokinetics and pharmacodynamics of bivalirudin in young healthy Chinese volunteers.",
    "Acta Pharmacologica Sinica (2012) 33: 1387-1394.",
    "doi:10.1038/aps.2012.37.",
    sep = " "
  )
  vignette <- "Zhang_2012_bivalirudin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; time-fixed within the single study).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-kg multiplier on CL, V1, Q, and V2. Zhang 2012 reports all PK structural parameters in per-kg units (Table 2: 0.323 L/h/kg, 0.086 L/kg, 0.0957 L/h/kg, 0.0554 L/kg). The model recovers the individual values as cl = exp(lcl + etalcl) * WT (and analogously for vc, q, vp). Bivalirudin dosing in the four study arms is also weight-adjusted (0.5 / 0.75 / 1.05 mg/kg IV bolus and 1.75 mg/kg/h IV infusion). Cohort total-body-weight range 50-78 kg, group means 55.4-61.6 kg (Zhang 2012 Table 1).",
      source_name        = "Weight"
    ),
    RBC = list(
      description        = "Red blood cell (erythrocyte) count.",
      units              = "10^12 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline value (Zhang 2012 Table 1 demographic table; the asterisk in the paper's column heading 'RBC*' marks the per-subject demographic value). Reference value 4.40 x 10^12 cells/L is the across-cohort median used by Zhang 2012 in the centred linear-deviation effect on EC50: EC50_i = theta_EC50 * exp(eta_EC50) * (1 + 1.70 * (RBC - 4.40)) per the equation printed on page 1391. Cohort range 3.79-5.17 x 10^12 cells/L across the four dose groups (Table 1). The biological rationale offered in the Discussion (page 1391-1392) is dilution of clotting-factor concentration during the in-vitro ACT measurement when RBC is low, not a direct pharmacological effect of erythrocytes on bivalirudin potency.",
      source_name        = "RBC*"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "29-37 years",
    age_median     = "33 years (cohort mean across the four dose groups)",
    weight_range   = "50-78 kg",
    weight_median  = "approximately 60 kg (group means 55.4-61.6 kg, Zhang 2012 Table 1)",
    sex_female_pct = 52.8,
    race_ethnicity = "100% Chinese Han",
    disease_state  = "Healthy adult volunteers in a phase I study; no anticoagulation indication; placebo-controlled randomisation.",
    dose_range     = "Bivalirudin 0.5 mg/kg IV bolus, 0.75 mg/kg IV bolus, 1.05 mg/kg IV bolus, and 0.75 mg/kg IV bolus followed by 1.75 mg/kg/h IV infusion for 4 h (n = 9 per arm).",
    regions        = "China (Peking University First Hospital, Beijing).",
    rbc_distribution = "Group mean +/- SD (range), x 10^12 cells/L: 4.4 +/- 0.1 (4.22-4.61); 4.5 +/- 0.4 (4.12-5.17); 4.5 +/- 0.4 (3.85-5.02); 4.3 +/- 0.2 (3.79-4.58); covariate reference value used in the EC50 equation is 4.40 (Zhang 2012 Table 1).",
    sampling_design = "423 plasma bivalirudin concentrations and a paired ACT time course per subject (LC-MS/MS for bivalirudin concentration; Hemochron Response device for ACT). Placebo-treated n = 12 subjects were not included in modelling but informed the placebo-corrected ACT baseline.",
    notes          = "48 subjects randomised (36 bivalirudin, 12 placebo). The placebo arm is excluded from the population PK-PD modelling but contributes to baseline / placebo correction (Zhang 2012 Methods). All subjects have normal renal function (GFR > 90 mL/min); Discussion notes that renal-impairment effects on bivalirudin CL reported in references 19-20 could not be assessed in this cohort. Baseline demographics from Zhang 2012 Table 1."
  )

  ini({
    # ====================================================================
    # Population PK -- two-compartment IV with no retained covariates.
    # All structural parameters are reported per kilogram body weight in
    # Zhang 2012 Table 2; the model() block multiplies each by WT to obtain
    # the individual L/h or L value.
    # ====================================================================
    lcl <- log(0.323);  label("Clearance per kg body weight (L/h/kg)")                          # Zhang 2012 Table 2 (CL = 0.323 L/h/kg, SE 2.6%; bootstrap median 0.323, 95% CI 0.307-0.339)
    lvc <- log(0.086);  label("Central volume of distribution per kg body weight (L/kg)")       # Zhang 2012 Table 2 (V1 = 0.086 L/kg, SE 4.2%; bootstrap median 0.0855, 95% CI 0.079-0.094)
    lq  <- log(0.0957); label("Inter-compartmental clearance per kg body weight (L/h/kg)")      # Zhang 2012 Table 2 (Q = 0.0957 L/h/kg, SE 3.0%; bootstrap median 0.0959, 95% CI 0.0908-0.102)
    lvp <- log(0.0554); label("Peripheral volume of distribution per kg body weight (L/kg)")    # Zhang 2012 Table 2 (V2 = 0.0554 L/kg, SE 3.3%; bootstrap median 0.0555, 95% CI 0.052-0.059)

    # PK inter-individual variability -- log-normal exponential model
    # (Zhang 2012 Methods 'Random effect model': P_i = theta * exp(eta_i),
    # eta ~ N(0, omega)). Reported CV% values from Table 2 are translated
    # to NONMEM variance via omega^2 = log(1 + CV^2). Q has no IIV
    # (Table 2 footnote a: 'IIV of Q was fixed at 0 during model
    # estimation, i.e., Q has no IIV'); the model has no etalq term.
    etalcl ~ log(0.148^2 + 1)   # Zhang 2012 Table 2 (CV 14.8%, SE 19.4%; bootstrap median 14.2%, 95% CI 8.63-20.0%)
    etalvc ~ log(0.242^2 + 1)   # Zhang 2012 Table 2 (CV 24.2%, SE 29.1%; bootstrap median 23.3%, 95% CI 8.18-36.2%)
    etalvp ~ log(0.156^2 + 1)   # Zhang 2012 Table 2 (CV 15.6%, SE 14.7%; bootstrap median 15.1%, 95% CI 9.75-19.4%)

    # PK residual error -- proportional only. Zhang 2012 Table 2 reports
    # 'Proportional error (%) 0.08' with SE 10.2% (bootstrap median 0.08,
    # 95% CI 0.07-0.09). The numeric value 0.08 is read as a fraction
    # (8% proportional residual), consistent with the bootstrap CI scale,
    # typical LC-MS/MS assay precision for bivalirudin, and the PD table's
    # standalone-percent reporting (PD residual 4.67%, encoded as 0.0467
    # below).
    propSd <- 0.08;  label("Proportional residual error on plasma Cc (fraction)")               # Zhang 2012 Table 2 (proportional error = 0.08 fraction, SE 10.2%)

    # ====================================================================
    # Population PK-PD -- direct-response sigmoid Emax (Hill coefficient
    # fixed at 1) using the central-compartment plasma concentration as
    # the effect site (Zhang 2012 Methods: 'Ce is the bivalirudin
    # concentration of the central compartment, which was also the effect
    # compartment'). The paper writes the equation as
    #   E = E0 + [(Emax - E0) * Ce^gamma] / (Ce^gamma + EC50^gamma)
    # with gamma = 1 in the final model (the abstract and Conclusion
    # explicitly note 'sigmoid Emax model without the Hill coefficient').
    # ====================================================================
    lemax <- log(318);  label("Maximum bivalirudin effect on activated clotting time (s)")     # Zhang 2012 Table 3 (Emax = 318 s, SE 2.39%; bootstrap median 320, 95% CI 304-340)
    lec50 <- log(2.44); label("EC50 of bivalirudin on ACT at reference RBC = 4.40 x 10^12/L (mg/L)") # Zhang 2012 Table 3 (EC50 = 2.44 mg/L, SE 11.8%; bootstrap median 2.51, 95% CI 1.88-2.99)
    le0   <- log(134);  label("Baseline ACT in absence of bivalirudin (s)")                    # Zhang 2012 Table 3 (E0 = 134 s, SE 0.98%; bootstrap median 134, 95% CI 131-136)

    # Covariate effect of RBC on EC50 -- linear-deviation form printed on
    # page 1391: EC50_i = theta_EC50 * exp(eta_EC50) * (1 + theta_RBC *
    # (RBC_i - 4.40)). theta_RBC has units (L per 10^12 cells), is applied
    # as a multiplicative factor in model(), and is centred on the cohort
    # median 4.40 x 10^12 cells/L. See vignette Assumptions and deviations
    # for the negative-EC50 behaviour at RBC below ~3.81.
    e_rbc_ec50 <- 1.70; label("Linear-deviation coefficient of RBC on EC50 (per 10^12 cells/L)") # Zhang 2012 Table 3 (theta_RBC = 1.70, SE 3.54%; bootstrap median 1.69, 95% CI 0.23-1.73)

    # PD inter-individual variability -- log-normal exponential model on
    # Emax, EC50, and E0 (Zhang 2012 Table 3). Same CV% -> omega^2
    # conversion as for PK IIV: omega^2 = log(1 + CV^2).
    etalemax ~ log(0.068^2 + 1)   # Zhang 2012 Table 3 (CV 6.80%, SE 28.5%; bootstrap median 6.33%, 95% CI 0.50-9.33%)
    etalec50 ~ log(0.464^2 + 1)   # Zhang 2012 Table 3 (CV 46.4%, SE 28.6%; bootstrap median 27.1%, 95% CI 0.45-62.6%)
    etale0   ~ log(0.041^2 + 1)   # Zhang 2012 Table 3 (CV 4.10%, SE 19.2%; bootstrap median 4.17%, 95% CI 2.27-5.58%)

    # PD residual error -- proportional. Zhang 2012 Table 3 reports
    # 'Proportional error, % 4.67' with SE 17.0% (bootstrap 2.45-6.06%).
    # Value encoded as fraction = 0.0467.
    propSd_ACT <- 0.0467;  label("Proportional residual error on activated clotting time (fraction)") # Zhang 2012 Table 3 (proportional error = 4.67% = 0.0467 fraction, SE 17.0%)
  })

  model({
    # Individual PK parameters. Zhang 2012 reports per-kg disposition
    # parameters in Table 2; the model multiplies each by total body
    # weight (WT, kg) so that downstream simulations consume per-subject
    # rather than per-kg dosing tables. Q has no IIV (Table 2 footnote).
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq)           * WT
    vp <- exp(lvp + etalvp) * WT

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition. Dose enters the central compartment
    # directly via the user data set (cmt = central / amt rows for IV
    # bolus; rate-coded rows for the IV infusion arm). Bivalirudin is
    # administered intravenously only (Zhang 2012 Methods).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Central-compartment plasma concentration in mg/L (dose units mg, vc
    # units L). Zhang 2012 reports bivalirudin concentrations in the same
    # mg/L scale that the EC50 = 2.44 mg/L is given in.
    Cc <- central / vc

    # Individual PD parameters.
    emax <- exp(lemax + etalemax)
    e0   <- exp(le0   + etale0)
    # Linear-deviation RBC effect on EC50, page 1391:
    #   EC50_i = theta_EC50 * exp(eta_EC50) * (1 + theta_RBC * (RBC - 4.40))
    # The covariate factor can become non-positive for RBC < ~3.81; the
    # cohort minimum (RBC = 3.79, Table 1 group 4 range) falls just inside
    # that region. This is a property of the published linear-deviation
    # form; see the vignette Assumptions and deviations section for the
    # implication and the recommended simulated-cohort RBC bounds.
    ec50 <- exp(lec50 + etalec50) * (1 + e_rbc_ec50 * (RBC - 4.40))

    # Direct-response sigmoid Emax with Hill = 1 (Zhang 2012 Conclusion:
    # 'sigmoid Emax model without the Hill coefficient'):
    #   ACT = E0 + (Emax - E0) * Cc / (Cc + EC50)
    # The central-compartment concentration drives the effect with no
    # delay (Zhang 2012 Methods: central compartment is the effect site).
    ACT <- e0 + (emax - e0) * Cc / (Cc + ec50)

    Cc  ~ prop(propSd)
    ACT ~ prop(propSd_ACT)
  })
}
