Liu_2020_cefepime <- function() {
  description <- "Two-compartment population PK model for intravenous cefepime fitted jointly to adults (febrile neutropenia with haematological malignancy, and critical illness) and hospitalised children aged 2 months to 16 years. Central volume scales linearly with total body weight standardised to 70 kg. The elimination rate constant is an additive non-renal intercept plus a renal component proportional to raw Cockcroft-Gault creatinine clearance standardised to 120 mL/min, the whole sum multiplied by an allometric (WT/70)^-0.25 scaler. Inter-compartmental distribution is parameterised directly as the micro-constants KCP and KPC. Estimated with the Pmetrics non-parametric adaptive grid (NPAG); typical values are the medians of the non-parametric joint density and the diagonal omegas are log-normal reconstructions from the published per-parameter CV% (see the vignette Assumptions and deviations). Liu 2020, n = 70 subjects (34 adults, 36 children), 683 plasma concentrations."
  reference <- "Liu J, Neely M, Lipman J, Sime F, Roberts JA, Kiel PJ, Avedissian SN, Rhodes NJ, Scheetz MH. Development of Population and Bayesian Models for Applied Use in Patients Receiving Cefepime. Clin Pharmacokinet. 2020;59(8):1027-1036. doi:10.1007/s40262-020-00873-3. PMID 32140940; PMC7222999. Preliminary results presented as IDWeek 2019 abstract 1573, Open Forum Infect Dis. 2019;6(Suppl 2):S575 (PMC6810068)."
  vignette <- "Liu_2020_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefepime is given as an intravenous infusion and
  # total plasma cefepime was measured (Liu 2020 Methods, 'Study
  # Populations'; the 20% protein binding is applied only afterwards, in
  # 'Simulations and Probability of Target Attainment').
  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Standardised to 70 kg (Liu 2020 Methods, 'Pharmacokinetic Models': 'Weight and CrCL were standardized to 70 kg and 120 mL/min, respectively'). Enters central volume LINEARLY as (WT / 70) and the elimination rate constant through the allometric scaler (WT / 70)^-0.25 (ESM 'Supplementary Materials' equations 3-4; Table 1 final row). Cohort means [SD]: adults 82.7 [21.5] kg, children 16.0 [16.1] kg (Results, 'Demographics').",
      source_name = "weight"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation, applied to all subjects including children",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RAW Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalised to mL/min/1.73 m^2 (Liu 2020 Methods, 'Model Building Populations': 'The Cockcroft-Gault formula (applied to all subjects) served as a standardized descriptor for the elimination rate constant', citing Cockcroft and Gault 1976). The adult formula was used for the children too, which is why the paediatric mean (195.6 mL/min) exceeds the adult mean. Enters the renal arm of the elimination rate constant linearly as (CRCL / 120). Cohort means [SD]: adults 106.7 [58.4] mL/min, children 195.6 [40.5] mL/min (Results, 'Demographics'). Matches the convention of the same group's Valadez_2025_cefepime.R.",
      source_name = "CrCL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 4L,
    n_samples = 683L,
    age_range = "Adults 22-82 years (mean 55.4); children approximately 2 months to 16 years (mean 3.9) (Results, 'Demographics')",
    weight_mean_sd = "Adults 82.7 [21.5] kg; children 16.0 [16.1] kg (mean [SD]; Results, 'Demographics')",
    crcl_mean_sd = "Adults 106.7 [58.4] mL/min; children 195.6 [40.5] mL/min, raw Cockcroft-Gault (mean [SD]; Results, 'Demographics')",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = "Adults: febrile neutropenia with haematological malignancy after chemotherapy and/or stem cell transplant (Sime, n = 12; Whited, n = 9) and critical illness (Roberts DALI, n = 13). Children: hospitalised with presumed or documented bacterial infection, cefepime monotherapy (Reed 1997, n = 36).",
    dose_range = "Maximum labelled dosing in the adult febrile-neutropenia studies; the paediatric study used intravenous cefepime only. The simulations evaluated 2 g every 8 h (adults) and 50 mg/kg (800 mg at 16 kg) every 8 h (children), each infused over 0.5 h.",
    regions = "United States and Australia (pooled studies; Sime and Roberts DALI are Australian-led, Whited and Reed are US).",
    protein_binding = "Not fitted. The model outputs TOTAL plasma cefepime as Cc; the paper applies a literature protein binding of 20% (unbound fraction 0.80) when computing fT>MIC (Methods, 'Simulations and Probability of Target Attainment').",
    notes = "Model-building set: 9 adults (Whited) and 24 randomly selected children (Reed), 428 observations, 0.5-249.7 ug/mL. Evaluation set: 12 adults (Sime), 13 adults (Roberts DALI) and the 12 remaining children. The four studies contributed 683 concentrations in the Results text; the Abstract states 604. Estimated with the NPAG algorithm in Pmetrics 1.5.2 for R. Final model: population R^2 0.87, bias 0.53, imprecision 7.75; Bayesian posterior R^2 0.965, bias -0.15, imprecision 1.07 (Table 1)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters: Liu 2020 Table 2, 'Median (95% CI)' column of
    # the Pmetrics NPAG non-parametric joint density (CI is a credibility
    # interval). The median is the natural typical value for the
    # log-normal random effects reconstructed below. Equations are ESM
    # 'Supplementary Materials' equations 1-4, repeated in the Table 1
    # final row.
    # ---------------------------------------------------------------------

    lvc <- log(11.172)
    label("Central volume of distribution V0 at WT = 70 kg (L)")
    # Liu 2020 Table 2: V0 = 11.172 L (95% CI 9.40, 12.50), CV 22.66%,
    # shrinkage 14.7%. ESM equation 3: Vc = V0 * (weight/70).

    lkel_nonren <- log(0.506)
    label("Non-renal elimination rate constant KeIntcpt at WT = 70 kg (1/h)")
    # Liu 2020 Table 2: KeIntcpt = 0.506 1/h (95% CI 0.33, 0.71), CV
    # 59.24%, shrinkage 15.1%. Methods, 'Pharmacokinetic Models':
    # 'non-renal elimination (KeIntcpt)'.

    lkel_renal <- log(0.236)
    label("Renal elimination rate constant Ke0 at CRCL = 120 mL/min and WT = 70 kg (1/h)")
    # Liu 2020 Table 2: Ke0 = 0.236 1/h (95% CI 0.11, 0.40), CV 133.95%,
    # shrinkage 6.8%. Methods: 'renal elimination descriptor (Ke0
    # vectorized as a function of glomerular filtration estimates)'. The
    # paper's symbol Ke0 is NOT the effect-compartment rate constant
    # (registered as lke0); it is the renal arm of the elimination rate
    # constant, hence lkel_renal.

    lk12 <- log(1.716)
    label("Central-to-peripheral rate constant KCP (1/h)")
    # Liu 2020 Table 2: KCP = 1.716 1/h (95% CI 1.02, 3.11), CV 68.82%,
    # shrinkage 7.7%. ESM: 'KCP is the rate constant from central to
    # peripheral compartment'.

    lk21 <- log(1.502)
    label("Peripheral-to-central rate constant KPC (1/h)")
    # Liu 2020 Table 2: KPC = 1.502 1/h (95% CI 1.35, 1.90), CV 60.63%,
    # shrinkage 8.8%. ESM: 'KPC is the rate constant from peripheral to
    # central compartment'.

    e_wt_vc <- fixed(1)
    label("Exponent of (WT / 70) on central volume (unitless)")
    # ESM equation 3, Vc = V0 * (weight/70); Table 1 final row 'V0 *
    # (WT/70)'. A structural constant of the model, not estimated.

    e_wt_kel <- fixed(-0.25)
    label("Allometric exponent of (WT / 70) on the elimination rate constant (unitless)")
    # ESM equation 4, Ke = (KeIntcpt + (CrCL/120) * Ke0) * (weight/70)^-0.25;
    # Methods, 'Pharmacokinetic Models': 'an allometric scaler was applied
    # to standardized weight ... raised to the negative 0.25th power'.
    # Multiplies the WHOLE bracket (ESM equation 4 and Table 1 final row),
    # although the Table 2 footnote mentions only V0 and KeIntcpt for the
    # weight standardisation.

    # ---------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete
    # non-parametric joint density; Table 2 summarises each marginal by
    # its CV%, carried here as a log-normal random effect with
    # omega^2 = log(CV^2 + 1). This is a parametric approximation of a
    # non-parametric distribution (precedent: Tsai_2023_ceftriaxone.R).
    # The CV% is consistent with the Table 3 covariance diagonal
    # (SD / CV recovers means of 11.1, 0.49, 0.53, 2.23, 1.96). Table 3's
    # off-diagonal covariances are NOT carried: transformed to the log
    # scale with sigma_ij = log(1 + Cov_ij / (m_i * m_j)) they give a
    # matrix that is not positive definite. See the vignette.
    # ---------------------------------------------------------------------
    etalvc ~ 0.050073 # Table 2 V0 CV 22.66% -> log(0.2266^2 + 1)
    etalkel_nonren ~ 0.300799 # Table 2 KeIntcpt CV 59.24% -> log(0.5924^2 + 1)
    etalkel_renal ~ 1.027567 # Table 2 Ke0 CV 133.95% -> log(1.3395^2 + 1)
    etalk12 ~ 0.387721 # Table 2 KCP CV 68.82% -> log(0.6882^2 + 1)
    etalk21 ~ 0.313057 # Table 2 KPC CV 60.63% -> log(0.6063^2 + 1)

    # ---------------------------------------------------------------------
    # Residual error. Pmetrics assay-error polynomial SD = C0 + C1 * Y with
    # C0 = 0.5 and C1 = 0.15 (Methods, 'Pharmacokinetic Models'), scaled by
    # a multiplicative factor gamma that was 'initially set at 4'; the
    # final gamma is not reported, so the polynomial is carried unscaled
    # (gamma = 1). The SDs add linearly in Pmetrics, hence combined1().
    # ---------------------------------------------------------------------
    addSd <- fixed(0.5)
    label("Additive residual error C0 (mg/L)")
    # Liu 2020 Methods: 'C0 and C1 equal to 0.5 and 0.15, respectively'.

    propSd <- fixed(0.15)
    label("Proportional residual error C1 (fraction)")
    # Liu 2020 Methods: 'C0 and C1 equal to 0.5 and 0.15, respectively'.
  })

  model({
    # Covariate reference values (Methods, 'Pharmacokinetic Models').
    wt_ref <- 70 # kg
    crcl_ref <- 120 # mL/min, raw Cockcroft-Gault

    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc
    kel_nonren <- exp(lkel_nonren + etalkel_nonren)
    kel_renal <- exp(lkel_renal + etalkel_renal)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # ESM equation 4: Ke = (KeIntcpt + (CrCL/120) * Ke0) * (weight/70)^-0.25
    kel <- (kel_nonren + kel_renal * CRCL / crcl_ref) * (WT / wt_ref)^e_wt_kel

    # ESM equations 1-2.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma cefepime (mg / L = ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
