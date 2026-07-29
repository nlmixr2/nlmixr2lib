Lacy_2018_cabozantinib_tumor <- function() {
  description <- "Longitudinal sum-of-tumor-diameter (SOD) growth-inhibition PD model for oral cabozantinib in adults with advanced renal cell carcinoma (RCC) enrolled in the phase III METEOR study (Lacy 2018 exposure-response analysis, n=319 patients with 1637 evaluable tumor-diameter measurements). The tumor diameter Y follows first-order exponential growth at rate k_grow, with a saturable cabozantinib drug-effect of the form Cavg/(EC50 + Cavg) modulating a time-dependent decay rate decay(t) = k_dmax + k_dmax_tot * exp(-k_tol * t). The k_dmax term is the non-attenuating asymptotic decay rate, k_dmax_tot is the magnitude of the resistance-driven loss of decay rate, and k_tol governs the attenuation kinetics (paper-reported attenuation half-life 25.6 days). The drug input Cavg is the individual predicted daily average plasma cabozantinib concentration (ng/mL) carried as a time-varying CAV data column; the upstream popPK model is `Lacy_2018_cabozantinib` (Lacy 2018 popPK companion paper). Residual error is additive on Y (mm); IIV is exponential on Y(0), k_grow, k_dmax, and k_dmax_tot, with IIV on EC50 and k_tol fixed at a near-zero variance (paper Supplemental Table 3 footnote b)."
  reference <- paste(
    "Lacy S, Nielsen J, Yang B, Miles D, Nguyen L, Hutmacher M.",
    "Population exposure-response analysis of cabozantinib efficacy",
    "and safety endpoints in patients with renal cell carcinoma.",
    "Cancer Chemother Pharmacol. 2018;81(6):1061-1070.",
    "doi:10.1007/s00280-018-3579-7.",
    "Upstream popPK (driver of CAV): Lacy S et al.,",
    "Cancer Chemother Pharmacol. 2018;81(6):1071-1082;",
    "doi:10.1007/s00280-018-3581-0; see modellib('Lacy_2018_cabozantinib').",
    sep = " "
  )
  vignette <- "Lacy_2018_cabozantinib_exposure_response"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the time-varying drug input is supplied as the CAV data covariate, in ng/mL, derived from the upstream Lacy 2018 popPK)",
    concentration = "mm (tumor diameter; not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Time-varying individual predicted daily average plasma cabozantinib concentration (ng/mL).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Carried as a per-record column; updated when dose changes (dose holds set CAV = 0; reductions and escalations recompute CAV from the current daily dose and the subject-specific popPK clearance). Derived in the source paper from individual empirical-Bayes post-hoc parameters of the upstream Lacy 2018 popPK model fit to nine pooled clinical studies (Lacy S et al., Cancer Chemother Pharmacol 2018;81(6):1071-1082; doi:10.1007/s00280-018-3581-0). Source paper's predicted steady-state Cavg values: 375 ng/mL at 20 mg/day, 750 ng/mL at 40 mg/day, 1125 ng/mL at 60 mg/day in RCC patients. The validation vignette shows the recommended cohort construction using the upstream popPK model.",
      source_name        = "Cavg"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 319L,
    n_observations   = 1637L,
    n_studies        = 1L,
    age_range        = "ages reported only at cohort level for METEOR (paper Methods cites enrolment criteria of >=18 years; demographic breakdown for the 319-patient tumor-model subset is not separately reported)",
    weight_range     = "not separately reported for the tumor-model subset",
    sex_female_pct   = NA_real_,
    race_ethnicity   = "not separately reported for the tumor-model subset",
    disease_state    = "Advanced or metastatic renal cell carcinoma (RCC) with clear-cell histology and measurable disease per RECIST. Subjects had received at least one prior VEGFR-TKI therapy. METEOR randomised 658 patients 1:1 to cabozantinib 60-mg tablet QD vs everolimus 10-mg QD; the tumor-model fit uses 319 patients on the cabozantinib arm with at least one evaluable post-baseline tumor diameter.",
    dose_range       = "Cabozantinib 60-mg tablet (Cabometyx) QD starting dose; dose reductions allowed to 40 mg and 20 mg per protocol AE-management rules; dose interruptions also permitted.",
    regions          = "Multinational phase III trial (NCT01865747); regional breakdown not reported in the Lacy 2018 ER paper.",
    notes            = "Baseline cohort characteristics for the dose-modification analysis (n=317 with at least one dose record) are in Lacy 2018 ER Supplemental Table 2. The tumor-model cohort (n=319) is the largest analysis subset (one fewer than the dose-modification subset because two patients had no evaluable post-baseline tumor diameter)."
  )

  ini({
    # ---- Structural parameters (Lacy 2018 ER Supplemental Table 3, "Longitudinal Tumor Growth" rows) ----
    # The "Transformed Estimate" column reports back-transformed typical values
    # in their natural units. All structural parameters are log-transformed
    # here in ini() so the model() back-transformation `exp(lX + etaX)`
    # reproduces the table's transformed estimate at eta = 0.
    lrbase_tumor <- log(63.1)    ; label("Baseline tumor diameter Y(0) (mm)")                                            # ER Supplemental Table 3: Baseline Tumor Size (mm) = 63.1 (90% CI 58.9, 67.5)
    lkgrow       <- log(0.00155) ; label("First-order tumor growth rate constant k_grow (1/day)")                        # ER Supplemental Table 3: k_grow (1/day) = 0.00155 (90% CI 0.00133, 0.0018)
    lkdmax       <- log(0.00125) ; label("Maximum non-attenuating drug-induced tumor decay rate k_dmax (1/day)")         # ER Supplemental Table 3: k_dmax (1/day) = 0.00125 (90% CI 0.000984, 0.00158)
    lkdmaxtot    <- log(0.00835) ; label("Magnitude of attenuating drug-induced decay-rate loss k_dmax_tot (1/day)")     # ER Supplemental Table 3: k_dmax_tot (1/day) = 0.00835 (90% CI 0.00689, 0.0101)
    lktol        <- log(0.0271)  ; label("Attenuation (resistance) rate constant k_tol (1/day)")                         # ER Supplemental Table 3: k_tol (1/day) = 0.0271 (90% CI 0.0238, 0.0308); attenuation half-life = log(2)/k_tol = 25.6 days
    lec50        <- log(251)     ; label("Half-maximal effective cabozantinib concentration EC50 (ng/mL)")               # ER Supplemental Table 3: EC50 (ng/mL) = 251 (90% CI 169, 375)

    # ---- Inter-individual variability ----
    # Supplemental Table 3 reports IIV under "(omega)" column header in
    # variance units (omega^2 on the log-normal scale); see paper text
    # "Baseline tumor size variability between patients was large (CV = 72%)"
    # which matches sqrt(omega^2) = sqrt(0.522) = 0.722 as the linearised CV
    # approximation. The companion popPK paper (Lacy_2018_cabozantinib) uses
    # the same omega^2-in-place convention (footnote d, "omega^2_Ka = 2.063").
    etalrbase_tumor ~ 0.522                                  # ER Supplemental Table 3: IIV Base (omega) = 0.522 (90% CI 0.45, 0.594)
    etalkgrow       ~ 0.313                                  # ER Supplemental Table 3: IIV k_grow (omega) = 0.313 (90% CI 0.218, 0.408)
    etalkdmax       ~ 0.353                                  # ER Supplemental Table 3: IIV k_dmax (omega) = 0.353 (90% CI 0.224, 0.482)
    etalkdmaxtot    ~ 0.641                                  # ER Supplemental Table 3: IIV k_dmax_tol (omega) = 0.641 (90% CI 0.469, 0.814); supplement spells the IIV row k_dmax_tol but the structural row is k_dmax_tot -- treated here as a supplement typo for the same parameter
    etalec50        ~ fixed(0.02)                            # ER Supplemental Table 3 footnote b: IIV EC50 (omega) = 0.02 (fixed value)
    etalktol        ~ fixed(0.02)                            # ER Supplemental Table 3 footnote b: IIV k_tol (omega) = 0.02 (fixed value)

    # ---- Residual error ----
    # Paper Methods: "an additive error model for residual variability"
    # (Results section "Longitudinal sum of tumor diameter model"). The
    # supplement column labels this row "Residual Variability (SD) (mm)" so
    # 5.75 is already the SD in mm.
    addSd <- 5.75 ; label("Additive residual SD on tumor diameter (mm)")  # ER Supplemental Table 3: Residual Variability (SD) (mm) = 5.75 (90% CI 5.52, 6)
  })

  model({
    # ---- Individual parameters (log-normal IIV) ----
    rbase_tumor <- exp(lrbase_tumor + etalrbase_tumor)
    kgrow       <- exp(lkgrow       + etalkgrow)
    kdmax       <- exp(lkdmax       + etalkdmax)
    kdmaxtot    <- exp(lkdmaxtot    + etalkdmaxtot)
    ktol        <- exp(lktol        + etalktol)
    ec50        <- exp(lec50        + etalec50)

    # ---- Time-dependent drug-effect decay rate ----
    # Eq. 5 of Lacy 2018 ER decomposes the cabozantinib-induced tumor decay
    # rate into a non-attenuating asymptote k_dmax and an exponentially
    # attenuating component that starts at +k_dmax_tot and decays with
    # half-life log(2)/k_tol = 25.6 days, modelling acquired resistance:
    #   decay(t) = k_dmax + k_dmax_tot * exp(-k_tol * t)
    # At t = 0 the total decay rate is k_dmax + k_dmax_tot (maximum effect,
    # no resistance yet); as t -> infinity the rate asymptotes to k_dmax.
    # The simulation-time variable `t` is rxode2's running solver time and
    # should be aligned with t = 0 at the start of cabozantinib therapy in
    # downstream simulations.
    decay_rate <- kdmax + kdmaxtot * exp(-ktol * t)

    # ---- Saturable cabozantinib drug effect ----
    # Cavg enters as a Hill-1 / Michaelis-Menten function with EC50:
    # at Cavg = EC50 the effect is half-maximal, and at Cavg = 4 * EC50
    # (~ EC80) the effect is approximately 0.8 * decay(t). For the paper's
    # 60-mg starting dose Cavg = 1125 ng/mL > 4 * EC50 = 1004 ng/mL, so the
    # drug effect is near plateau.
    drug_effect <- decay_rate * CAV / (ec50 + CAV)

    # ---- Tumor growth ODE (Eq. 5) ----
    # dY/dt = k_grow * Y - drug_effect(t, Cavg) * Y
    d/dt(tumor_size) <- kgrow * tumor_size - drug_effect * tumor_size
    tumor_size(0)    <- rbase_tumor

    # ---- Observation ----
    tumor_size ~ add(addSd)
  })
}
