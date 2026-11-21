Wang_2025_ritlecitinib_eyelash <- function() {
  description <- paste0(
    "Longitudinal exposure-response model of the eyelash assessment (ELA) ",
    "score during ritlecitinib treatment in adolescents and adults with ",
    "alopecia areata and abnormal eyelashes at baseline (ELA < 3; study ",
    "B7981015, phase 2b/3). ELA is a four-level ordinal score (0 = no ",
    "eyelashes, 3 = normal eyelashes) described by a cumulative-probit ",
    "ordinal regression: P(ELA <= j) = Phi(alpha_j + placebo effect - drug ",
    "effect), with alpha_2 = B1, alpha_1 = B1 - B2 and alpha_0 = B1 - B2 - B3 ",
    "(B2, B3 > 0). The drug effect is a stimulated-kin indirect-response ",
    "latent variable (kin = kout) driven by an Emax/EC50 function of the ",
    "average ritlecitinib concentration between adjacent ELA records; Asian ",
    "race raises Emax by 64.7% and female sex lowers B1 by 1.19 probit units ",
    "(less severe baseline eyelash loss). The typical placebo effect is fixed ",
    "to zero, with additive between-subject variability on its maximum and ",
    "the same onset rate as the drug effect. There is no PK component -- ",
    "exposure enters as the CAV column, the empirical-Bayes average ",
    "concentration from the companion ritlecitinib population PK model ",
    "(modellib('Wojciechowski_2023_ritlecitinib_final')). Companion eyebrow ",
    "model from the same paper: modellib('Wang_2025_ritlecitinib_eyebrow')."
  )

  reference <- paste(
    "Wang Y, Huh Y, Lejeune A.",
    "Population exposure-response analysis of the effect of ritlecitinib",
    "on eyebrow assessment and eyelash assessment in patients with",
    "alopecia areata.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):389-402.",
    "doi:10.1002/psp4.13283.",
    sep = " "
  )

  vignette <- "Wang_2025_ritlecitinib_eyebrow_eyelash"

  # There is no PK component: ritlecitinib exposure must be supplied from
  # outside the model as the per-interval average concentration CAV. Wang 2025
  # took it from the empirical Bayes estimates of the companion population PK
  # model, packaged here as `Wojciechowski_2023_ritlecitinib_final`.
  depends <- "CAV"

  # `etab1` and `etapmax` are ADDITIVE on linear-scale parameters (Section
  # 2.3.2: "The IIV on baseline B1 and Pmax was evaluated using an additive
  # error model since these parameters do not need to be constrained to
  # positivity"), so the canonical `etal<name>` + `l<name>` log-normal pairing
  # does not apply. Same pattern as `Huh_2024_ritlecitinib`'s `etapmax`.
  paper_specific_etas <- c("etab1", "etapmax")

  # The single state is the dimensionless drug-effect latent variable A_drug
  # of Wang 2025 Eq 6-8, starting at its untreated steady state kin/kout = 1.
  compartmentData <- list(
    effect = list(
      analyte = "drug-effect latent variable A_drug of Wang 2025 Eq 6",
      units = "(dimensionless)",
      specimen = "not applicable",
      verified = FALSE
    )
  )

  units <- list(
    time = "week (time since the first study treatment record)",
    dosing = "(none; no PK component -- ritlecitinib exposure enters as the CAV covariate column, not as dosing events)",
    concentration = "(ordinal eyelash assessment score 0-3; the observation eyelash_assessment is a category and the outputs p0-p3 are the probabilities of each score, not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description = "Average ritlecitinib plasma concentration over the interval between the previous ELA record and the current ELA record.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying per record. Wang 2025 Section 2.2: derived from the empirical Bayes estimates of the final ritlecitinib population PK model (Wojciechowski 2023; this package's `Wojciechowski_2023_ritlecitinib_final`) combined with the dosing diary, then averaged over the interval between two adjacent ELA records. Set to 0 for placebo periods; the Emax term then collapses to 0. The same exposure metric was used for the SALT-score model of `Huh_2024_ritlecitinib`. For orientation, Wang 2025 Section 3.2 quotes the steady-state Cavg of 50 mg QD as 52 ng/mL and of 30 mg QD as 27 ng/mL.",
      source_name = "CAVE"
    ),
    RACE_ASIAN = list(
      description = "Asian-race indicator (1 = Asian, 0 = any other race).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or Other race).",
      notes = "Multiplicative fractional effect on Emax: emax * (1 + 0.647 * RACE_ASIAN). Derived from the supplementary control stream's `IF(RACE.EQ.3) COVASIAN = 1 + THETA(7)`, where RACE = 3 is Asian (the Table 2 race order White, Black, Asian, Other). Table 2: 142/535 (26.5%) Asian.",
      source_name = "RACE"
    ),
    SEXF = list(
      description = "Female-sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male).",
      notes = "Additive shift on the probit-scale baseline B1: b1 + (-1.19) * SEXF. Derived from the supplementary control stream's `IF(SEX.EQ.2) SEXB1 = THETA(8)` (source SEX 1 = male, 2 = female, so SEXF = as.integer(SEX == 2)); Table 4 labels the row 'Effect of Female on B1'. A lower B1 lowers every cumulative probability P(ELA <= j), so females are less likely to have severe baseline eyelash loss (Section 3.2; Table S3). Table 2: 313/535 (58.5%) female.",
      source_name = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened on B1, Emax and Pmax (Wang 2025 Section 2.3.3) but not retained. Table 2: median 33 years (range 12-72)."
    ),
    REGION = list(
      description = "Geographic location of the enrolling site.",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened on B1, Emax and Pmax (Wang 2025 Section 2.3.3) but not retained. Table 2: North America 221 (41.3%), Europe 121 (22.6%), Asia 110 (20.6%), Rest of World 83 (15.5%)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 535L,
    n_studies = 1L,
    n_observations = 4929L,
    age_range = "12-72 years (Table 2 median 33, mean 34.4, SD 14.5)",
    weight_range = "35.0-200.0 kg (Table 2 median 68.7, mean 70.9, SD 17.8)",
    sex_female_pct = 58.5,
    race_ethnicity = "White 363 (67.9%), Asian 142 (26.5%), Black 17 (3.2%), Other 13 (2.4%) (Table 2)",
    disease_state = "Alopecia areata with >= 50% scalp hair loss and abnormal eyelashes at baseline (ELA < 3). Baseline ELA score 0 in 283 (52.9%), 1 in 166 (31.0%) and 2 in 86 (16.1%) participants (Table 2).",
    dose_range = "Study B7981015 (phase 2b/3): ritlecitinib 200 mg QD x 4 weeks then 50 or 30 mg QD, 50, 30 or 10 mg QD for 48 weeks, or placebo for 24 weeks followed by 200/50 mg or 50 mg QD for 24 weeks.",
    regions = "North America 41.3%, Europe 22.6%, Asia 20.6%, Rest of World 15.5% (Table 2).",
    notes = "Participants with a normal ELA score (= 3) at baseline were excluded from the analysis. ELA was assessed at pre-dose and weeks 4, 8, 12, 18, 24, 28, 34, 40 and 48 (Table S1)."
  )

  ini({
    # ================================================================
    # Cumulative-probit ordinal regression -- Wang 2025 Eq 1-5 and the
    # supplementary ELA control stream (run56.mod, Data S2).
    #
    #   P(ELA <= 2) = Phi(B1            + PBOEFF - DEEFF)
    #   P(ELA <= 1) = Phi(B1 - B2       + PBOEFF - DEEFF)
    #   P(ELA <= 0) = Phi(B1 - B2 - B3  + PBOEFF - DEEFF)
    #
    # with B2 = EXP(THETA(2)) and B3 = EXP(THETA(3)) (control stream $PK),
    # so the Table 4 'B2' and 'B3' rows are the log-scale THETA estimates
    # (symmetric asymptotic CIs; B2 = exp(1.12) = 3.06 and B3 = exp(1.18) =
    # 3.25 on the probit scale). The drug effect is SUBTRACTED, as in the
    # control stream's $ERROR block; see the vignette for the sign of Eq 4.
    #
    # Time is in weeks here; the control stream works in hours and converts
    # the half-life with log(2)/(168 * THETA(6)).
    # ================================================================

    b1 <- 7.31
    label("Probit-scale cut-off of P(ELA <= 2) in males, the baseline parameter B1 (probit units)")
    # Table 4 'Baseline (B1)' = 7.31 (95% CI 6.74, 7.89; SIR median 7.29).

    lb2 <- 1.12
    label("Log of the probit-scale gap between the P(ELA <= 2) and P(ELA <= 1) cut-offs (B2 = exp(lb2) = 3.06)")
    # Table 4 'Difference in baselines for score 3 and score 2 (B2)' = 1.12 (95% CI 1.06, 1.19; SIR median 1.12); control stream B2 = EXP(THETA(2)).

    lb3 <- 1.18
    label("Log of the probit-scale gap between the P(ELA <= 1) and P(ELA <= 0) cut-offs (B3 = exp(lb3) = 3.25)")
    # Table 4 'Difference in baselines for score 2 and score 1 (B3)' = 1.18 (95% CI 1.12, 1.23; SIR median 1.18); control stream B3 = EXP(THETA(3)).

    e_sexf_b1 <- -1.19
    label("Additive shift in B1 for female participants (probit units)")
    # Table 4 'Effect of Female on B1' = -1.19 (95% CI -1.77, -0.610; SIR median -1.18); control stream IF(SEX.EQ.2) SEXB1 = THETA(8).

    pmax <- fixed(0)
    label("Typical maximum placebo effect on the probit scale (probit units)")
    # Section 3.2: the estimate of theta_Pmax 'was not significantly different from zero, and therefore theta_Pmax was fixed to zero'; control stream POPPMAX = 0.

    lemax <- log(7.75)
    label("Log of the maximum stimulation of the drug-effect latent variable (Emax = 7.75, probit units)")
    # Table 4 'Maximum Drug effect (Emax)' = 7.75 (95% CI 5.45, 10.1; SIR median 7.78).

    lec50 <- log(exp(10.2) / 1000)
    label("Log of the average ritlecitinib concentration producing half of Emax (ng/mL; EC50 = 26.9 ng/mL)")
    # Table 4 'Drug concentration producing half of the maximum stimulation (EC50)' = 10.2 (95% CI 9.69, 10.7); Table 4 note: 'Transformed EC50 = exp(10.2)/1000 = 26.90 (ng/mL)', matching the control stream EC50 = EXP(POPLNEC50)/1000.

    lthalfrec_drug <- log(27.3)
    label("Log of the drug-effect onset half-life (weeks; kin = kout = log(2)/thalfrec_drug, shared by the placebo effect)")
    # Table 4 'Half-Life of Drug Onset (weeks, kin,drug)' = 27.3 (95% CI 20.9, 33.7; SIR median 27.6); Table 4 note: kin,drug = log(2)/(168 * half-time) per hour.

    e_race_asian_emax <- 0.647
    label("Fractional change in Emax for Asian participants (unitless)")
    # Table 4 'Effect of Asian on Emax' = 0.647 (95% CI 0.323, 0.970; SIR median 0.644); control stream COVASIAN = 1 + THETA(7).

    # ================================================================
    # Inter-individual variability (Table 4; variances). B1 and Pmax are
    # additive on the probit scale, Emax is exponential.
    # ================================================================
    etab1 ~ 9.05 # Table 4 'omega2 B1' = 9.05 (95% CI 7.36, 10.7); shrinkage 16.1%
    etapmax ~ 12.9 # Table 4 'omega2 Pmax' = 12.9 (95% CI 7.63, 18.3); shrinkage 42.5%
    etalemax ~ 0.408 # Table 4 'omega2 Emax' = 0.408 (95% CI 0.285, 0.531); shrinkage 40.0%
  })

  model({
    # --- Individual parameters ----------------------------------------
    b1_i <- b1 + e_sexf_b1 * SEXF + etab1
    b2 <- exp(lb2)
    b3 <- exp(lb3)
    pmax_i <- pmax + etapmax
    emax_i <- exp(lemax + etalemax) * (1 + e_race_asian_emax * RACE_ASIAN)
    ec50 <- exp(lec50)

    # kin,drug = kout,drug (Section 3.2), and the placebo onset rate equals
    # the drug onset rate (Section 3.2; control stream POPPBOKOFF).
    thalfrec_drug <- exp(lthalfrec_drug)
    kout <- log(2) / thalfrec_drug

    # --- Drug effect: stimulated-kin indirect response (Eq 6-8) -------
    # CAV = 0 collapses the stimulation, which is how placebo periods are
    # encoded.
    d/dt(effect) <- kout * (1 + emax_i * CAV / (ec50 + CAV)) - kout * effect
    effect(0) <- 1
    drug_effect <- effect - 1

    # --- Placebo effect (Eq 9-11) --------------------------------------
    # The ELA control stream integrates A_pbo as a second compartment,
    # dA_pbo/dt = -kout * A_pbo with A_pbo(0) = 1; its exact solution
    # exp(-kout * t) is used here, as in the EBA control stream.
    placebo_effect <- pmax_i * (1 - exp(-kout * t))

    # --- Cumulative probit and per-score probabilities (Eq 4-5) -------
    lp <- b1_i + placebo_effect - drug_effect
    pcum2 <- phi(lp)
    pcum1 <- phi(lp - b2)
    pcum0 <- phi(lp - b2 - b3)
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- pcum2 - pcum1
    p3 <- 1 - pcum2

    # rxode2 ordinal endpoint: categories 1-4 carry probabilities p0, p1, p2
    # and 1 - p0 - p1 - p2, i.e. simulated category k is ELA score k - 1.
    eyelash_assessment ~ c(p0, p1, p2)
  })
}
