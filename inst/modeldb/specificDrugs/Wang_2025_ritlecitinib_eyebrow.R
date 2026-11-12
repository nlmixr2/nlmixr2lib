Wang_2025_ritlecitinib_eyebrow <- function() {
  description <- paste0(
    "Longitudinal exposure-response model of the eyebrow assessment (EBA) ",
    "score during ritlecitinib treatment in adolescents and adults with ",
    "alopecia areata and abnormal eyebrows at baseline (EBA < 3; study ",
    "B7981015, phase 2b/3). EBA is a four-level ordinal score (0 = no ",
    "eyebrows, 3 = normal eyebrows) described by a cumulative-probit ",
    "ordinal regression: P(EBA <= j) = Phi(alpha_j + placebo effect - drug ",
    "effect), with alpha_2 = B1, alpha_1 = B1 - B2 and alpha_0 = B1 - B2 - B3 ",
    "(B2, B3 > 0). The drug effect is a stimulated-kin indirect-response ",
    "latent variable (kin = kout) driven by an Emax/EC50 function of the ",
    "average ritlecitinib concentration between adjacent EBA records; Asian ",
    "race raises Emax by 48.7%. The typical placebo effect is fixed to zero, ",
    "with additive between-subject variability on its maximum and the same ",
    "onset rate as the drug effect. There is no PK component -- exposure ",
    "enters as the CAV column, the empirical-Bayes average concentration ",
    "from the companion ritlecitinib population PK model ",
    "(modellib('Wojciechowski_2023_ritlecitinib_final')). Companion eyelash ",
    "model from the same paper: modellib('Wang_2025_ritlecitinib_eyelash')."
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
    concentration = "(ordinal eyebrow assessment score 0-3; the observation eyebrow_assessment is a category and the outputs p0-p3 are the probabilities of each score, not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description = "Average ritlecitinib plasma concentration over the interval between the previous EBA record and the current EBA record.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying per record. Wang 2025 Section 2.2: derived from the empirical Bayes estimates of the final ritlecitinib population PK model (Wojciechowski 2023; this package's `Wojciechowski_2023_ritlecitinib_final`) combined with the dosing diary, then averaged over the interval between two adjacent EBA records. Set to 0 for placebo periods; the Emax term then collapses to 0. The same exposure metric was used for the SALT-score model of `Huh_2024_ritlecitinib`. For orientation, Wang 2025 Section 3.2 quotes the steady-state Cavg of 50 mg QD as 52 ng/mL and of 30 mg QD as 27 ng/mL.",
      source_name = "CAVE"
    ),
    RACE_ASIAN = list(
      description = "Asian-race indicator (1 = Asian, 0 = any other race).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or Other race).",
      notes = "Multiplicative fractional effect on Emax: emax * (1 + 0.487 * RACE_ASIAN). Derived from the supplementary control stream's `IF(RACE.EQ.3) COVASIAN = 1 + THETA(7)`, where RACE = 3 is Asian (the Table 1 race order White, Black, Asian, Other). Table 1: 160/594 (26.9%) Asian.",
      source_name = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened on B1, Emax and Pmax (Wang 2025 Section 2.3.3) but not retained. Table 1: median 33 years (range 12-72)."
    ),
    SEXF = list(
      description = "Female-sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on B1, Emax and Pmax (Wang 2025 Section 2.3.3) and retained for the eyelash model only, not for this eyebrow model. Table 1: 355/594 (59.8%) female."
    ),
    REGION = list(
      description = "Geographic location of the enrolling site.",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened on B1, Emax and Pmax (Wang 2025 Section 2.3.3) but not retained. Table 1: North America 239 (40.2%), Europe 130 (21.9%), Asia 128 (21.5%), Rest of World 97 (16.3%)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 594L,
    n_studies = 1L,
    n_observations = 5456L,
    age_range = "12-72 years (Table 1 median 33, mean 34.3, SD 14.4)",
    weight_range = "35.0-200.0 kg (Table 1 median 68.3, mean 70.5, SD 17.8)",
    sex_female_pct = 59.8,
    race_ethnicity = "White 398 (67.0%), Asian 160 (26.9%), Black 21 (3.5%), Other 15 (2.5%) (Table 1)",
    disease_state = "Alopecia areata with >= 50% scalp hair loss and abnormal eyebrows at baseline (EBA < 3). Baseline EBA score 0 in 328 (55.2%), 1 in 182 (30.6%) and 2 in 84 (14.1%) participants (Table 1).",
    dose_range = "Study B7981015 (phase 2b/3): ritlecitinib 200 mg QD x 4 weeks then 50 or 30 mg QD, 50, 30 or 10 mg QD for 48 weeks, or placebo for 24 weeks followed by 200/50 mg or 50 mg QD for 24 weeks.",
    regions = "North America 40.2%, Europe 21.9%, Asia 21.5%, Rest of World 16.3% (Table 1).",
    notes = "Participants with a normal EBA score (= 3) at baseline were excluded from the analysis. EBA was assessed at pre-dose and weeks 4, 8, 12, 18, 24, 28, 34, 40 and 48 (Table S1)."
  )

  ini({
    # ================================================================
    # Cumulative-probit ordinal regression -- Wang 2025 Eq 1-5 and the
    # supplementary EBA control stream (run56.mod, Data S1).
    #
    #   P(EBA <= 2) = Phi(B1            + PBOEFF - DEEFF)
    #   P(EBA <= 1) = Phi(B1 - B2       + PBOEFF - DEEFF)
    #   P(EBA <= 0) = Phi(B1 - B2 - B3  + PBOEFF - DEEFF)
    #
    # with B2 = EXP(THETA(2)) and B3 = EXP(THETA(3)) (control stream $PK),
    # so the Table 3 'B2' and 'B3' rows are the log-scale THETA estimates
    # (symmetric asymptotic CIs; B2 = exp(1.01) = 2.75 and B3 = exp(1.11) =
    # 3.03 on the probit scale). The drug effect is SUBTRACTED, as in the
    # control stream's $ERROR block; see the vignette for the sign of Eq 4.
    #
    # Time is in weeks here; the control stream works in hours and converts
    # the half-life with log(2)/(168 * THETA(6)).
    # ================================================================

    b1 <- 6.05
    label("Probit-scale cut-off of P(EBA <= 2), the baseline parameter B1 (probit units)")
    # Table 3 'Baseline (B1)' = 6.05 (95% CI 5.70, 6.40; SIR median 6.04).

    lb2 <- 1.01
    label("Log of the probit-scale gap between the P(EBA <= 2) and P(EBA <= 1) cut-offs (B2 = exp(lb2) = 2.75)")
    # Table 3 'Difference in baselines for score 3 and score 2 (B2)' = 1.01 (95% CI 0.948, 1.07; SIR median 1.01); control stream B2 = EXP(THETA(2)).

    lb3 <- 1.11
    label("Log of the probit-scale gap between the P(EBA <= 1) and P(EBA <= 0) cut-offs (B3 = exp(lb3) = 3.03)")
    # Table 3 'Difference in baselines for score 2 and score 1 (B3)' = 1.11 (95% CI 1.06, 1.16; SIR median 1.11); control stream B3 = EXP(THETA(3)).

    pmax <- fixed(0)
    label("Typical maximum placebo effect on the probit scale (probit units)")
    # Section 3.2: the estimate of theta_Pmax 'was not significantly different from zero, and therefore theta_Pmax was fixed to zero'; control stream POPPMAX = 0.

    lemax <- log(11.0)
    label("Log of the maximum stimulation of the drug-effect latent variable (Emax = 11.0, probit units)")
    # Table 3 'Maximum Drug effect (Emax)' = 11.0 (95% CI 6.98, 15.0; SIR median 10.7).

    lec50 <- log(exp(10.6) / 1000)
    label("Log of the average ritlecitinib concentration producing half of Emax (ng/mL; EC50 = 40.1 ng/mL)")
    # Table 3 'Drug concentration producing half of the maximum stimulation (EC50)' = 10.6 (95% CI 10.1, 11.1); Table 3 note: 'Transformed EC50 = exp(10.6)/1000 = 40.13 (ng/mL)', matching the control stream EC50 = EXP(POPLNEC50)/1000.

    lthalfrec_drug <- log(34.1)
    label("Log of the drug-effect onset half-life (weeks; kin = kout = log(2)/thalfrec_drug, shared by the placebo effect)")
    # Table 3 'Half-Life of Drug Onset (weeks, kin,drug)' = 34.1 (95% CI 25.1, 43.2; SIR median 33.9); Table 3 note: kin,drug = log(2)/(168 * half-time) per hour.

    e_race_asian_emax <- 0.487
    label("Fractional change in Emax for Asian participants (unitless)")
    # Table 3 'Effect of Asian on Emax' = 0.487 (95% CI 0.195, 0.778; SIR median 0.490); control stream COVASIAN = 1 + THETA(7).

    # ================================================================
    # Inter-individual variability (Table 3; variances). B1 and Pmax are
    # additive on the probit scale, Emax is exponential.
    # ================================================================
    etab1 ~ 6.53 # Table 3 'omega2 B1' = 6.53 (95% CI 5.36, 7.69); shrinkage 15.8%
    etalemax ~ 0.409 # Table 3 'omega2 Emax' = 0.409 (95% CI 0.285, 0.533); shrinkage 40.9%
    etapmax ~ 21.8 # Table 3 'omega2 Pmax' = 21.8 (95% CI 11.6, 31.9); shrinkage 37.5%
  })

  model({
    # --- Individual parameters ----------------------------------------
    b1_i <- b1 + etab1
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

    # --- Placebo effect (Eq 9-11, solved in closed form as in the
    # control stream's $ERROR block) -----------------------------------
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
    # and 1 - p0 - p1 - p2, i.e. simulated category k is EBA score k - 1.
    eyebrow_assessment ~ c(p0, p1, p2)
  })
}
