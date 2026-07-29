Bisaso_2014_albumin <- function() {
  description <- "Semi-mechanistic disease-progression model for plasma albumin in Ugandan adults co-infected with HIV (on efavirenz-based ART) with or without active TB (on rifampicin-based anti-TB co-treatment). Indirect-response form with zero-order hepatic albumin secretion that follows a Verhulst logistic transition from baseline rate Q0 to steady-state rate Qss, and first-order plasma elimination at a fixed literature rate."
  reference <- "Bisaso KR, Owen JS, Ojara FW, Namuwenge PM, Mugisha A, Mbuagbaw L, Luboobi LS, Mukonzo JK. Characterizing plasma albumin concentration changes in TB/HIV patients on anti retroviral and anti-tuberculosis therapy. In Silico Pharmacology. 2014;2:3. doi:10.1186/s40203-014-0003-9"
  vignette <- "Bisaso_2014_albumin"
  units <- list(
    time = "day",
    dosing = "g",
    concentration = "g/dL"
  )
  # Note on dosing: the Bisaso 2014 model does not consume rxode2 dose events
  # (no drug dosing into 'central'); the population were on efavirenz-based
  # HAART and rifampicin-based anti-TB therapy and the model's t = 0 is the
  # time of ART initiation. The dosing-mass placeholder 'g' is recorded for
  # convention compatibility (matches the concentration numerator g/dL) and to
  # document the units a user would supply if they ever introduced exogenous
  # albumin (e.g., IV human albumin therapy) via cmt = 'central'.

  covariateData <- list(
    TB_POS = list(
      description = "Active tuberculosis co-infection at study entry indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = active TB co-infection at study entry (treated with rifampicin-based anti-TB therapy -- ethambutol/isoniazid/rifampicin/pyrazinamide induction then rifampicin/isoniazid maintenance -- initiated at least 2 weeks before ART); 0 = HIV monoinfection (no active TB). Time-fixed per subject. Reduces baseline albumin secretion rate Q0 by 30.8% relative to the HIV-only reference (Bisaso 2014 Table 3 derived: 0.0864 / 0.1248 - 1; paper text frames the same effect as a 44.2% reduction relative to the TB-HIV cohort).",
      source_name = "TB"
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 c.3435C>T (rs1045642, exon 26, synonymous Ile1145Ile) mutant allele carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = subject carries at least one T allele at ABCB1 rs1045642 (heterozygous CT or homozygous TT; pooled per Bisaso 2014 Table 3 -- the 3435TT cohort was only n=1); 0 = homozygous wild type (CC). Time-fixed per subject (germline genotype). Increases baseline albumin secretion rate Q0 by 16.7% relative to the CC wild-type reference (Bisaso 2014 Table 3 derived: 0.1008 / 0.0864 - 1; paper text reports as a 16% increase).",
      source_name = "ABCB13435"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 174,
    n_studies = 1,
    age_range = "adults; whole-cohort median 33 years (IQR 29-39); TB-HIV median 31 (28-37); HIV-only median 37 (31-42)",
    weight_range = "whole-cohort median 51 kg (IQR 47-58); TB-HIV median 50 (45-53); HIV-only median 55 (50-60)",
    sex_female_pct = 52.9,
    race_ethnicity = c(Black = 100),
    disease_state = "HIV monoinfection (n=104 in full cohort) or HIV-TB co-infection (n=158 in full cohort), ART-naive, initiated on combination ART (efavirenz / lamivudine / zidovudine); TB co-infected received ethambutol / isoniazid / rifampicin / pyrazinamide for 2 months followed by isoniazid + rifampicin for 4 months",
    regions = "Uganda (Mulago National Referral Hospital, Kampala)",
    dose_range = "(no drug PK modelled; the population was on efavirenz-based HAART and, for TB co-infected subjects, rifampicin-based anti-TB co-treatment)",
    notes = "Demographics from Bisaso 2014 Table 1 (full cohort n=262). Per Methods 'Model evaluation' the dataset was randomly split: a model-building subset of ~2/3 (n=174) used here for the population parameter estimates, and a validation subset (n=88, 269 observations) used to compute mean / RMSE prediction error. Average of 3 albumin observations per subject; serum albumin sampled on days 1, 3, 7, 14, 21, 42, 56, and 84 after ART initiation. Baseline median serum albumin was 3.02 g/dL overall (2.57 in TB-HIV, 3.91 in HIV only). Genotype counts (full cohort, Bisaso 2014 Table 1): ABCB1 3435CC=205, CT=56, TT=1. CD4 count and viral load were tested as time-varying covariates on Qss and R but were not retained in the final model."
  )

  ini({
    # ----- Structural parameters (Bisaso 2014 Table 3, original-dataset column) -----
    # Reference for the typical-value Q0 is the HIV-only / ABCB1-CC subgroup
    # (TB_POS = 0, SNP_ABCB1_RS1045642 = 0). The paper reports Q0 = 0.0864 with
    # TB-HIV as the implicit reference; we re-anchor to the canonical-register
    # reference categories so that "TB_POS = 0" maps to the no-comorbidity
    # group (matching HIV_POS / SNP_* register conventions). The model
    # behaviour is identical to the paper -- only the labelling of the
    # reference subgroup changes.
    lq0  <- log(0.1248)
    label("Baseline albumin secretion rate at HIV-only / ABCB1-CC reference Q0 (g/dL/day)")  # Bisaso 2014 Table 3 'Q 0 (HIV only)' row
    lqss <- log(0.1464)
    label("Steady-state albumin secretion rate following treatment Qss (g/dL/day)")           # Bisaso 2014 Table 3 'Q ss (g/dl/day)' row
    lR   <- log(0.0072)
    label("Rate of change from Q0 to Qss R (1/day)")                                          # Bisaso 2014 Table 3 'R (1/day)' row
    lK   <- fixed(log(0.0336))
    label("First-order plasma albumin elimination rate constant K (1/day; half-life 20.6 d)") # Bisaso 2014 Table 2 / Table 3 'K (1/day) ... FIX' and Methods 'Data analysis' (fixed to literature value, half-life 20.6 d)

    # ----- Covariate effects on Q0 (multiplicative, linear on natural scale) -----
    # Encoded as Q0 = exp(lq0) * (1 + e_tb_pos_q0 * TB_POS) * (1 + e_snp_abcb1_rs1045642_q0 * SNP_ABCB1_RS1045642)
    # so that TB_POS = 0 and SNP_ABCB1_RS1045642 = 0 reproduce exp(lq0) = 0.1248.
    e_tb_pos_q0            <- -0.308
    label("Multiplicative shift on Q0 for TB co-infection (fraction)")
    # Bisaso 2014 derived from Table 3: 0.0864 / 0.1248 - 1 = -0.308 (paper text frames the same effect as a 44.2% reduction relative to the TB-HIV cohort)
    e_snp_abcb1_rs1045642_q0 <- 0.167
    label("Multiplicative shift on Q0 for ABCB1 c.3435C>T (rs1045642) carriage (fraction)")
    # Bisaso 2014 derived from Table 3: 0.1008 / 0.0864 - 1 = +0.167 (paper text reports as a 16% increase)

    # ----- IIV -- only Q0 retained in the final model -----
    # Per Results: "The random effects on parameters Qss and R had high shrinkage
    # (>40%) and had a very low variability of less than 10^-6 and were
    # therefore dropped from the model." Only IIV on Q0 (15.0% CV per
    # Table 3 original-dataset column) is included.
    # Variance on the log scale: log(CV^2 + 1) = log(0.15^2 + 1) = 0.02225.
    etalq0 ~ 0.02225  # Bisaso 2014 Table 3 IIV_Q0 = 15.0% CV; log-normal IIV per Methods 'Data analysis'

    # ----- Residual error -- proportional (Bisaso 2014 Methods 'Data analysis': "modeled as proportional but additive and additive plus proportional error models were also tested") -----
    propSd <- 0.182
    label("Proportional residual error on plasma albumin (fraction)")  # Bisaso 2014 Table 3 'Residual error (proportional) (%CV)' = 18.2
  })

  model({
    # ----- Typical-value parameters (paper symbols on natural scale) -----
    # Q0 with multiplicative TB and ABCB1 effects.
    q0_typ  <- exp(lq0) * (1 + e_tb_pos_q0 * TB_POS) *
      (1 + e_snp_abcb1_rs1045642_q0 * SNP_ABCB1_RS1045642)
    qss_typ <- exp(lqss)
    R_typ   <- exp(lR)
    K       <- exp(lK)

    # ----- Individual Q0 (log-normal IIV) -----
    q0 <- q0_typ * exp(etalq0)

    # ----- Time-varying albumin secretion rate Q(t) -----
    # Verhulst logistic transition from q0 (at t = 0) to qss (as t -> inf)
    # at rate R, per Bisaso 2014 Eq. 5 (analytical solution of
    # dQ/dt = R * Q * (1 - Q/Qss)):
    #   Q(t) = q0 * qss / (q0 + (qss - q0) * exp(-R * t))
    # Validations: Q(0) = q0; Q(inf) = qss; reduces to qss if q0 == qss.
    Qt <- q0 * qss_typ / (q0 + (qss_typ - q0) * exp(-R_typ * t))

    # ----- Plasma albumin ODE (Bisaso 2014 Eq. 6) -----
    # dX/dt = Q(t) - K * X
    # Units: [g/dL/day] = [g/dL/day] - [1/day] * [g/dL]   (balanced)
    d/dt(central) <- Qt - K * central

    # Initial condition X(0) = Q0 / K (Bisaso 2014, derivation following Eq. 7).
    central(0) <- q0 / K

    # ----- Observation: plasma albumin concentration (g/dL) -----
    # Cc here is the plasma-compartment endogenous-biomarker concentration
    # (albumin), not a drug concentration. Single-output naming Cc per
    # nlmixr2lib convention; the canonical name is uniform across drug-PK
    # and endogenous-biomarker single-output models.
    Cc <- central
    Cc ~ prop(propSd)
  })
}
