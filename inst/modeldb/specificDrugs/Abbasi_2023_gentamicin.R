Abbasi_2023_gentamicin <- function() {
  description <- "One-compartment intravenous pooled meta-analytic population PK model for gentamicin in critically ill adult ICU patients (n=1215 pooled from 21 published studies; Abbasi 2023). Vd (0.33 +/- 0.20 L/kg), CL (4.70 +/- 2.89 L/h) and total body weight (70.8 +/- 19.9 kg) were pooled as mean +/- SD across the 21 studies and used as normal-distribution inputs to a Monte Carlo Simulation (10,000 virtual patients, Crystal Ball, Oracle) for probability-of-target-attainment (PTA) analysis of gentamicin once-daily doses 5-10 mg/kg (0.5-h IV infusion, 24-h dosing interval) against AUC24h/MIC and Cmax/MIC efficacy targets. The reported variability reflects combined between-study and between-patient effects rather than a NONMEM-fitted OMEGA."
  reference <- "Abbasi MY, Chaijamorn W, Wiwattanawongsa K, Charoensareerat T, Doungngern T. Recommendations of Gentamicin Dose Based on Different Pharmacokinetic/Pharmacodynamic Targets for Intensive Care Adult Patients: A Redefining Approach. Clin Pharmacol Adv Appl. 2023;15:67-76. doi:10.2147/CPAA.S417298"
  vignette <- "Abbasi_2023_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (dosing weight)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Abbasi 2023 Results: pooled dosing weight 70.8 +/- 19.9 kg (mean +/- SD) across the 21 published studies (n=1215 critically ill adult patients). WT enters the model as a linear (exponent = 1) structural scaler on Vc because the pooled typical Vd was reported per kg (0.33 L/kg). WT does NOT scale CL: the paper pooled CL as total L/h (4.70 +/- 2.89 L/h) rather than L/h/kg, so CL is encoded as WT-independent in this extraction. In the Monte Carlo Simulation WT was truncated at >= 40 kg (Methods, Monte Carlo Simulations: 'total body weight (kg) was set at the range 40 to infinity in MCS analysis'). Dose was computed as mg/kg * WT.",
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1215L,
    n_studies      = 21L,
    age_range      = "Not reported at the pooled level; underlying studies enrolled adult intensive care unit patients (Methods, Pharmacokinetic Model Development: 'critically ill adult patients'; pregnant women were excluded).",
    weight_range   = "Pooled dosing weight 70.8 +/- 19.9 kg (mean +/- SD across 21 studies); MCS truncated total body weight at >= 40 kg.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported at the pooled level.",
    disease_state  = "Critically ill adult patients admitted to medical, surgical, or traumatic intensive care units. Among the 1215 pooled patients, 90.5% were confirmed to have severe infections (Results). Exclusions: pregnant women, patients receiving extracorporeal membrane oxygenation, and patients receiving renal replacement therapy (Methods, Pharmacokinetic Model Development).",
    dose_range     = "5-10 mg/kg once daily as a 0.5-h intravenous infusion, 24-h dosing interval, evaluated over 72 h (3 daily doses).",
    regions        = "Global (systematic pooling across 21 previously published critically-ill-adult ICU gentamicin PK studies indexed in PubMed / EMBASE / SCOPUS / CINAHL / EBSCO up to September 2022; individual-study geography was not aggregated).",
    renal_function = "Heterogeneous ICU renal function; patients requiring renal replacement therapy were excluded. Seven studies contributed acute-phase (first-dose, within 48-72 h) PK; two studies contributed steady-state PK.",
    notes          = "Pooled meta-analytic PK parameters (Table S1 lists the 21 contributing studies; Table S2 lists 'additional parameters used in the model' [dosing / target / MIC inputs] and is not on disk in this extraction). No individual-level baseline demographic table appears in the paper. Search cut-off: September 2022. See references 14-34 of the paper for the 21 contributing studies. Additional (non-critically-ill) Vd range reported in Discussion: 0.27-0.83 L/kg."
  )

  ini({
    # Structural PK parameters (Abbasi 2023 Results, pooled means +/- SD across
    # the 21 contributing studies; n=1215 critically ill adult ICU patients).
    #
    # Reported pooled acute-phase values:
    #   Vd = 0.33 +/- 0.20 L/kg   (per kg of total body weight)
    #   CL = 4.70 +/- 2.89 L/h    (total; not per kg)
    #   Ke = 0.18 +/- 0.10 h^-1   (redundant with CL / (Vd * TBW))
    #
    # This extraction parameterises the model in CL and Vc (canonical popPK
    # form) and derives Ke = cl / vc inside model(). At the typical WT of
    # 70.8 kg the implied Ke is 4.70 / (0.33 * 70.8) = 0.201 h^-1, slightly
    # higher than the paper's independently pooled Ke of 0.18 h^-1; the
    # ~12% inconsistency is a known artefact of the paper's independent
    # pooling of three redundant parameters (see vignette Assumptions and
    # deviations).
    lcl <- log(4.70); label("Clearance CL (L/h; pooled total, not per kg)")  # Abbasi 2023 Results, pooled acute-phase mean
    lvc <- log(0.33); label("Volume of distribution Vd (L/kg; per kg of total body weight)")  # Abbasi 2023 Results, pooled acute-phase mean

    # Pooled variability. The paper reports mean +/- SD on the natural
    # (linear) scale and used those directly as normal-distribution inputs
    # to Crystal Ball (Methods, Monte Carlo Simulations). Because the
    # nlmixr2 canonical is log-normal random effects on structural PK
    # parameters, the reported %CV is transformed to variance on the log
    # scale via omega^2 = log(CV^2 + 1). The reported variability reflects
    # combined between-study and between-patient effects (not a
    # NONMEM-fitted OMEGA); see vignette Assumptions and deviations.
    #   CL : SD/mean = 2.89/4.70 = 0.6149 -> omega^2 = log(0.6149^2 + 1) = 0.31989
    #   Vd : SD/mean = 0.20/0.33 = 0.6061 -> omega^2 = log(0.6061^2 + 1) = 0.31278
    etalcl ~ 0.31989  # Abbasi 2023 Results, pooled CL SD 2.89 L/h (CV 61.5%)
    etalvc ~ 0.31278  # Abbasi 2023 Results, pooled Vd SD 0.20 L/kg (CV 60.6%)

    # Residual error: not reported. Abbasi 2023 does not fit an OMEGA/SIGMA
    # to observed concentrations; the Monte Carlo Simulation is
    # deterministic given each virtual patient's sampled PK parameters and
    # total body weight (Methods, Monte Carlo Simulations). propSd is
    # therefore fixed at 0 to reproduce the paper's noise-free simulation.
    propSd <- fixed(0); label("Proportional residual error (fraction; 0 - MCS is deterministic given sampled PK parameters)")  # Abbasi 2023 Methods (no RUV reported)
  })

  model({
    # Individual PK parameters. The published pooled Vd is reported per kg
    # of total body weight (Results); multiplying by actual WT returns the
    # absolute volume of distribution (L). CL is reported as total L/h
    # rather than L/h/kg and is encoded WT-independent to reproduce the
    # paper's pooled parameterisation exactly.
    cl <- exp(lcl + etalcl)       # L/h (total; not scaled by WT)
    vc <- exp(lvc + etalvc) * WT  # L   (per-kg * WT)

    # Derived first-order elimination rate.
    kel <- cl / vc

    # One-compartment IV disposition. Dosing enters 'central' via the input
    # event (rate / duration set on the dose row for the 0.5-h infusion).
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central / vc has units mg/L (matches the
    # paper's stated concentration units).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
