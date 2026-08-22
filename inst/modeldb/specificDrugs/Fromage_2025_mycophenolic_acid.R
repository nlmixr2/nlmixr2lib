Fromage_2025_mycophenolic_acid <- function() {
  description <- "One-compartment population PK model for mycophenolic acid (MPA) after oral enteric-coated mycophenolate sodium (EC-MPS, Myfortic) in a multi-indication cohort of solid-organ transplant, haematopoietic cell transplant and autoimmune-disease patients (Fromage 2025), with double-gamma absorption describing the characteristic double concentration peak, first-order elimination, a steady-state trough offset, and an indication effect on the second gamma rate constant. Evaluated in the closed form published by the authors (Monolix fitted the analytic solution), so the model is a steady-state single-dosing-interval model."
  reference <- "Fromage Y, Sayadi H, Koloskoff K, Labriffe M, Monchaud C, Marquet P, Woillard JB. Killing several birds with one stone: A multi-indication population pharmacokinetic model and Bayesian estimator for enteric-coated mycophenolate sodium. Br J Clin Pharmacol. 2025;91(5):1396-1408. doi:10.1111/bcp.16374"
  vignette <- "Fromage_2025_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte  = "enteric-coated mycophenolate sodium",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes    = paste(
        "Dose-registration state only, and it is deliberately always empty.",
        "Fromage 2025 published C(t) in closed form (Section 2.2.1) and fitted",
        "that analytic solution in Monolix, so this file evaluates the same",
        "closed form algebraically rather than integrating an equivalent ODE",
        "system. rxode2 nevertheless requires a dosed compartment with a",
        "defined d/dt() for podo() and tad() to resolve, which is what this",
        "state provides. f(depot) <- 0 keeps it at zero, so depot is NOT a",
        "meaningful amount; the plasma concentration is Cc."
      )
    )
  )

  covariateData <- list(
    TX_HEART = list(
      description = "Heart (cardiac) transplant recipient indicator",
      units       = "(binary)",
      type        = "binary",
      source_name = "Indication (cardiac transplantation)",
      notes       = paste(
        "Fromage 2025 Section 2.2.2 pools the indication into three groups and",
        "estimates ONE effect for group (ii): cardiac + pulmonary + bone marrow",
        "transplantation. This file keeps the three transplant types as separate",
        "canonical columns and forms the paper's group-(ii) indicator inside",
        "model(), so the paper-specific pooling stays out of the column names",
        "and a future model that separates the three effects can reuse the same",
        "columns. Reference group (i) = renal or hepatic transplantation."
      )
    ),
    TX_LUNG = list(
      description = "Lung (pulmonary) transplant recipient indicator",
      units       = "(binary)",
      type        = "binary",
      source_name = "Indication (pulmonary transplantation)",
      notes       = "Member of the pooled group (ii); see TX_HEART notes."
    ),
    TX_HCT = list(
      description = "Haematopoietic cell (bone marrow / stem cell) transplant recipient indicator",
      units       = "(binary)",
      type        = "binary",
      source_name = "Indication (bone marrow transplantation)",
      notes       = paste(
        "Fromage 2025 writes 'bone marrow transplantation' in the covariate",
        "definition (Section 2.2.2) and 'haematopoietic stem cell (HSC)",
        "recipients' in the Abstract; both denote the same group (ii) member.",
        "Registered under the canonical TX_HCT, which the TX_ANY register entry",
        "had already reserved for this concept. Member of the pooled group (ii);",
        "see TX_HEART notes."
      )
    ),
    DIS_AUTOIMMUNE = list(
      description = "Autoimmune-disease indication indicator (pooled)",
      units       = "(binary)",
      type        = "binary",
      source_name = "Indication (autoimmune diseases)",
      notes       = paste(
        "Fromage 2025 group (iii): systemic lupus erythematosus, lupus nephritis",
        "and nephrotic syndrome, pooled into a single indicator (Section 2.1 and",
        "Section 2.2.2). Pooled rather than DIS_SLE because the group also",
        "contains nephrotic-syndrome patients. Reference group (i) = renal or",
        "hepatic transplantation. NOTE: Fromage 2025 Section 3.2.2 reports this",
        "effect as non-significant (188% RSE) and states it 'can be disregarded',",
        "but the authors retained it in the model used for the pcVPC simulations,",
        "so it is retained here; see the vignette Errata."
      )
    )
  )

  # Covariates that Fromage 2025 screened but did NOT retain in the final
  # model (Section 2.2.2 lists the tested set; Section 3.2.2 reports that only
  # the indication effect on b2 met the retention criterion). Documented here
  # so the covariate screen's provenance is preserved without declaring
  # covariates that model() never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a qualitative covariate (Section 2.2.2); not retained. Sex was unknown for 34-37% of occasions (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested as a continuous covariate via the median-normalised power model (Section 2.2.2); not retained. Cohort median 45 years (range 6-80), Table 1."
    ),
    TAT = list(
      description = "Time since transplantation (time since EC-MPS initiation)",
      units       = "days",
      type        = "continuous",
      notes       = "Tested as a continuous covariate (Section 2.2.2); not retained. Cohort median 771 days (range 6-10488), Table 1."
    ),
    CONMED_CICLOSPORIN = list(
      description = "Ciclosporin co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as one level of a three-level immunosuppressive co-treatment covariate (none / ciclosporin / other) because ciclosporin inhibits the enterohepatic recirculation of MPA (Section 3.1); not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 129L,
    n_studies      = 1L,
    n_profiles     = 153L,
    n_observations = 863L,
    age_range      = "Median 45 years (range 6-80); 20 of 127 patients were paediatric (< 18 years)",
    weight_range   = "Not reported (body weight was not available in the ISBA extraction; Section 4)",
    sex_female_pct = 27,
    race_ethnicity = "Not reported (single-centre French cohort)",
    disease_state  = paste(
      "Multi-indication: renal, hepatic, cardiac and pulmonary transplantation,",
      "bone marrow / haematopoietic cell transplantation, and autoimmune disease",
      "(systemic lupus erythematosus, lupus nephritis, nephrotic syndrome).",
      "116 patients with solid-organ or haematopoietic cell transplantation and",
      "13 with autoimmune disease."
    ),
    dose_range     = "Oral EC-MPS (Myfortic) 180-1440 mg per intake; median 360 mg",
    administration = "Oral, enteric-coated delayed-release tablet",
    regions        = "France (single centre: CHU de Limoges, ISBA / ABIS 3.0 platform)",
    notes          = paste(
      "Retrospective real-world therapeutic-drug-monitoring data extracted from",
      "the ISBA Bayesian dosing-adaptation platform. Development set 75% /",
      "validation set 25%, split on PK profiles, plus a separate later",
      "extraction as an external validation set. 3-11 samples per profile",
      "(median 5). Quantification limit 0.5 mg/L with below-limit data retained",
      "and handled by the M4 method. Estimation by SAEM in Monolix.",
      "Concomitant immunosuppressants: none (12 occasions), ciclosporin (16),",
      "other (131). INTER-OCCASION VARIABILITY: the source also estimated IOV",
      "on b1 (SD 1.08, Table 2), which this file does NOT encode, per the",
      "nlmixr2lib convention that library model files omit IOV (no occasion",
      "column is assumed); see the vignette Errata."
    )
  )

  # Implementation notes (the vignette 'Assumptions and deviations' section
  # carries the full justification for each item):
  #
  # * STRUCTURE. Fromage 2025 Section 2.2.1 gives the absorption input rate as
  #   a two-component gamma mixture,
  #     Va(t) = F*D*[r*f1(t) + (1-r)*f2(t)],
  #     fi(t) = bi^ai / Gamma(ai) * t^(ai-1) * exp(-bi*t),
  #   and gives the resulting plasma concentration in closed form,
  #     C(t) = C0 + FAIV*D*( r    *(b1/(b1-lam))^a1 * P[a1,(b1-lam)t]
  #                        + (1-r)*(b2/(b2-lam))^a2 * P[a2,(b2-lam)t] )*exp(-lam*t)
  #   where P is the regularised lower incomplete gamma function. That closed
  #   form is what Monolix fitted, and it is what this file evaluates, using
  #   rxode2's gammap(a, z) for P[a, z]. Deriving it from the ODE
  #   d/dt(central) = Va(t) - lam*central with central(0) = 0 is an exact
  #   convolution identity, so the two are the same model; the closed form is
  #   used because it is the published one AND because it is numerically
  #   robust (see NUMERICS below).
  #
  # * NUMERICS. a1 = 1280.97 makes the first gamma component a near-delta
  #   pulse of width ~sqrt(a1)/b1 = 0.086 h. Integrating the equivalent ODE
  #   drives the central amount through denormal magnitudes (~1e-181) before
  #   the pulse arrives, so lsoda's relative-tolerance test collapses the step
  #   size: with rxode2 5.1.7 the ODE form returns all-NaN (liblsoda, which
  #   then reports Cc == the trough offset with NO error) on output grids finer
  #   than ~0.02 h unless atol is loosened to >= 1e-4. The closed form has no
  #   such failure mode and matches an independent R pgamma() evaluation of the
  #   published equation to 7e-15 mg/L at every grid density tested
  #   (0.05 / 0.01 / 0.005 / 0.001 h).
  #
  # * PARAMETERISATION. The gamma components are stored in the canonical
  #   transit parameterisation (ntr, mtt) rather than as the paper's (a, b),
  #   following Debord_2001_cyclosporin.R (the same Limoges group's gamma
  #   absorption model): ai = ntri + 1 and bi = ai/mtti, which is a bijection,
  #   so mtti = ai/bi is exactly the paper's mean absorption time MATi.
  #   Because ai carries no IIV, the paper's IIV on log(bi) transfers to
  #   log(mtti) with a sign flip and an unchanged variance, and likewise for
  #   the indication effect on log(b2) (see the e_* parameters below).
  #
  # * FAIV -> vc. FAIV is reported in mg/L, i.e. normalised to a 100 mg dose,
  #   and the paper prints V/F = 100/FAIV, so the amplitude FAIV*D of the
  #   closed form equals D/vc with vc = 100/FAIV = 12.71 L. That reproduces
  #   the Abstract's Vd/F = 12.60 L to within 0.9%. The paper's IIV on
  #   log(FAIV) transfers to log(vc) with a sign flip and unchanged variance.
  #
  # * C0 -> lrbase. C0 is defined as "the estimated trough concentration at
  #   steady state for a reference dose" and enters the closed form as a
  #   CONSTANT additive offset, NOT as a decaying initial condition: with
  #   lam = 1.81 1/h (t1/2 = 0.38 h) a decaying C0 would be ~0 by 12 h, which
  #   contradicts both the paper's own definition and the observed troughs of
  #   ~2 mg/L in Figure 2. It is stored under the canonical baseline prefix
  #   lrbase and added to the observed concentration, following
  #   Larsen_2018_factorviii_monkey.R (Cc <- central/vc + rbase).
  #
  # * SINGLE-INTERVAL SCOPE. tad() and podo() refer to the most recent dose,
  #   so each dose restarts the absorption input and there is no
  #   superposition. This matches the authors' framing (C0 IS the steady-state
  #   trough carried in from previous doses), but it means the model describes
  #   ONE steady-state dosing interval; it must not be used to build up
  #   accumulation across doses.
  #
  # * r. The source constrains r to (0, 1) with a logit-normal distribution
  #   (Section 2.2.2), but Table 2 reports no IIV on r, so only the typical
  #   value is needed and the log transform used here is exact for it.
  ini({
    # Double-gamma absorption, component 1 (the fast phase). Stored as
    # (ntr, mtt); a1 = ntr1 + 1 and b1 = a1/mtt1 recover Table 2 exactly.
    lntr1 <- log(1280.97 - 1); label("Gamma component 1 shape minus 1 (ntr1, dimensionless)")  # Fromage 2025 Table 2: a1 = 1280.97 (RSE 0.25%); ntr1 = a1 - 1
    lmtt1 <- log(1280.97 / 418.42); label("Gamma component 1 mean absorption time MAT1 (h)")   # Fromage 2025 Table 2: a1 = 1280.97, b1 = 418.42 1/h; MAT1 = a1/b1 = 3.061 h

    # Double-gamma absorption, component 2 (the slow phase; carries the
    # indication covariate).
    lntr2 <- log(395.34 - 1); label("Gamma component 2 shape minus 1 (ntr2, dimensionless)")   # Fromage 2025 Table 2: a2 = 395.34 (RSE 2.30%); ntr2 = a2 - 1
    lmtt2 <- log(395.34 / 95.69); label("Gamma component 2 mean absorption time MAT2 (h)")     # Fromage 2025 Table 2: a2 = 395.34, b2 = 95.69 1/h; MAT2 = a2/b2 = 4.131 h (the paper's Discussion quotes 4.09 h for this group; see vignette Errata)

    # Fraction of the dose absorbed through the fast (component 1) phase.
    lfdepot <- log(0.79); label("Fraction absorbed via the fast gamma phase (r, dimensionless)")  # Fromage 2025 Table 2: r = 0.79 (RSE 2.14%)

    # Disposition.
    lvc  <- log(100 / 7.87); label("Apparent central volume Vd/F (L)")                        # Fromage 2025 Table 2: FAIV = 7.87 mg/L (RSE 3.35%); V/F = 100/FAIV = 12.71 L (Abstract: 12.60 L)
    lkel <- log(1.81); label("First-order elimination rate constant lambda (1/h)")            # Fromage 2025 Table 2: lambda = 1.81 1/h (RSE 7.48%)

    # Steady-state trough offset (the paper's C0).
    lrbase <- log(2.21); label("Steady-state trough concentration offset C0 (mg/L)")          # Fromage 2025 Table 2: C0 = 2.21 mg/L (RSE 6.50%)

    # Indication effect on the second gamma component. Fromage 2025 estimates
    # these on log(b2); because mtt2 = a2/b2 with a2 covariate-free, the
    # effect on log(mtt2) is the NEGATIVE of the reported value.
    e_tx_heart_lung_hct_mtt2 <- -0.59; label("Cardiac / pulmonary / haematopoietic-cell transplant effect on log(MAT2)")  # Fromage 2025 Table 2: beta_b2_cardiac,pulmonary,bone marrow = +0.59 (RSE 47.1%) on log(b2); sign inverts on log(MAT2). Recovers b2 = 95.69*exp(0.59) = 172.62 1/h and MAT2 = 2.29 h (Discussion)
    e_dis_autoimmune_mtt2 <- 0.17; label("Autoimmune-disease effect on log(MAT2)")                                        # Fromage 2025 Table 2: beta_b2_autoimmune = -0.17 (RSE 188%) on log(b2); sign inverts on log(MAT2). Recovers b2 = 95.69*exp(-0.17) = 80.73 1/h and MAT2 = 4.90 h (Discussion)

    # Inter-individual variability. Fromage 2025 Table 2 reports IIV as
    # STANDARD DEVIATIONS of the exponential (log-scale) random effects
    # (table footnote: "IIV and IOV are expressed as SD"), so the variances
    # below are the reported SDs squared. IIV on log(b1) / log(b2) /
    # log(FAIV) transfers to log(mtt1) / log(mtt2) / log(vc) with a sign flip
    # and an unchanged variance (see implementation notes above).
    etalmtt1  ~ 0.8464  # Fromage 2025 Table 2 IIV b1 = 0.92 SD (RSE 44.9%); 0.92^2 = 0.8464
    etalmtt2  ~ 0.2704  # Fromage 2025 Table 2 IIV b2 = 0.52 SD (RSE 16.1%); 0.52^2 = 0.2704
    etalvc    ~ 0.0144  # Fromage 2025 Table 2 IIV FAIV = 0.12 SD (RSE 11.0%); 0.12^2 = 0.0144. NOTE: shrinkage of the conditional distribution was 83.6% for FAIV
    etalkel   ~ 0.1936  # Fromage 2025 Table 2 IIV lambda = 0.44 SD (RSE 14.2%); 0.44^2 = 0.1936
    etalrbase ~ 0.2401  # Fromage 2025 Table 2 IIV C0 = 0.49 SD (RSE 10.8%); 0.49^2 = 0.2401

    # Residual error: additive was preferred (Section 3.2.1).
    addSd <- 1.38; label("Additive residual SD (mg/L)")  # Fromage 2025 Table 2: residual additive error sigma = 1.38 mg/L (RSE 3.47%)
  })
  model({
    # Individual parameters
    ntr1   <- exp(lntr1)
    mtt1   <- exp(lmtt1 + etalmtt1)
    ntr2   <- exp(lntr2)
    vc     <- exp(lvc + etalvc)
    kel    <- exp(lkel + etalkel)
    fdepot <- exp(lfdepot)
    rbase  <- exp(lrbase + etalrbase)

    # Fromage 2025 group (ii) = cardiac OR pulmonary OR haematopoietic-cell
    # transplantation, which the source treats as a single level with one
    # estimated effect. Formed here as an inclusive-OR of the three canonical
    # transplant indicators, written as a product of complements so it stays
    # exact 0/1 arithmetic.
    txgrp2 <- 1 - (1 - TX_HEART) * (1 - TX_LUNG) * (1 - TX_HCT)

    mtt2 <- exp(lmtt2 +
                  e_tx_heart_lung_hct_mtt2 * txgrp2 +
                  e_dis_autoimmune_mtt2 * DIS_AUTOIMMUNE +
                  etalmtt2)

    # Gamma shape / rate parameters of the two absorption components
    # (Fromage 2025 Section 2.2.1). ai = ntri + 1, bi = ai/mtti.
    ga1 <- ntr1 + 1
    gb1 <- ga1 / mtt1
    ga2 <- ntr2 + 1
    gb2 <- ga2 / mtt2

    # Dose-registration state; always empty (see compartmentData notes).
    d/dt(depot) <- -kel * depot
    f(depot) <- 0

    # Published closed form (Fromage 2025 Section 2.2.1). gammap(a, z) is the
    # regularised lower incomplete gamma function P[a, z]. The FAIV*D
    # amplitude is D/vc because FAIV is normalised per 100 mg and
    # vc = 100/FAIV.
    tad1 <- tad()
    arm1 <- fdepot * (gb1 / (gb1 - kel))^ga1 * gammap(ga1, (gb1 - kel) * tad1)
    arm2 <- (1 - fdepot) * (gb2 / (gb2 - kel))^ga2 * gammap(ga2, (gb2 - kel) * tad1)

    Cc <- rbase + (podo() / vc) * (arm1 + arm2) * exp(-kel * tad1)
    Cc ~ add(addSd)
  })
}
