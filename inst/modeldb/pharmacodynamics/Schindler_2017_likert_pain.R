Schindler_2017_likert_pain <- function() {
  description <- "Population PD model of 11-point Likert pain score (0-10) in adults with painful distal diabetic neuropathy (placebo arm of three phase III studies), using the minimal continuous-time Markov model (mCTMM) form of Schindler & Karlsson 2017. Eleven ordered-categorical states (scores 0-10) evolve via ODEs (Eq. 6) with transfer rate constants (Eqs. 8-9) parameterised through a mean equilibration time (MET) and the steady-state score probabilities Pss,k. Steady-state cumulative probabilities Pss(Y >= k) follow a proportional-odds model driven by a linear effect of time (on logit scale, accounting for placebo effect) and a shift for concomitant acetaminophen use (Eq. 13). MET is itself time-varying with a linear slope on t (in weeks), consistent with the time-dependent increase in score stability reported in the paper. The model has no PK and no drug-exposure input -- the placebo effect is modelled directly on the score-probability logit."
  reference <- paste(
    "Schindler E, Karlsson MO.",
    "A Minimal Continuous-Time Markov Pharmacometric Model.",
    "AAPS J. 2017;19(5):1424-1435.",
    "doi:10.1208/s12248-017-0109-1.",
    "Companion count-model form for the same endpoint / cohort:",
    "Plan EL, Elshoff JP, Stockis A, Sargentini-Maier ML, Karlsson MO.",
    "Likert pain score modeling: a Markov integer model and an",
    "autoregressive continuous model.",
    "Clin Pharmacol Ther. 2012;91(4):820-828.",
    "doi:10.1038/clpt.2011.301. DDMODEL00000194.",
    "See modellib('Plan_2012_pain') for the count-model comparator.",
    "Sister mCTMM model files from the same paper:",
    "modellib('Schindler_2017_sunitinib_fatigue'),",
    "modellib('Schindler_2017_sunitinib_hfs').",
    sep = " "
  )
  vignette <- "Schindler_2017_mCTMM"
  paper_specific_compartments <- c(
    "pscore0", "pscore1", "pscore2", "pscore3", "pscore4",
    "pscore5", "pscore6", "pscore7", "pscore8", "pscore9", "pscore10"
  )

  units <- list(
    time = "h",
    dosing = "(none; placebo arm)",
    concentration = "(11-point Likert pain score 0-10, ordinal)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    pscore0  = list(analyte = "11-point Likert pain score 0", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore1  = list(analyte = "11-point Likert pain score 1", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore2  = list(analyte = "11-point Likert pain score 2", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore3  = list(analyte = "11-point Likert pain score 3", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore4  = list(analyte = "11-point Likert pain score 4", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore5  = list(analyte = "11-point Likert pain score 5", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore6  = list(analyte = "11-point Likert pain score 6", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore7  = list(analyte = "11-point Likert pain score 7", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore8  = list(analyte = "11-point Likert pain score 8", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore9  = list(analyte = "11-point Likert pain score 9", units = NA_character_, specimen = "administration site", verified = FALSE),
    pscore10 = list(analyte = "11-point Likert pain score 10", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list(
    CONMED_PARA = list(
      description        = "Concomitant paracetamol (acetaminophen) rescue-medication indicator on the current day: 1 = subject took acetaminophen, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant paracetamol)",
      notes              = "Time-varying per-observation flag. Enters the steady-state cumulative logit as an additive shift theta_ace * CONMED_PARA; a positive theta_ace shifts the logit upward, meaning acetaminophen use is associated with higher pain scores (rescue medication is taken because of pain, not vice versa; see paper Discussion). Matches the CONMED_PARA convention already used in modellib('Plan_2012_pain').",
      source_name        = "ACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 231L,
    n_studies      = 3L,
    n_observations = 22492L,
    age_range      = "adults with painful distal diabetic neuropathy",
    weight_range   = "not reported in Schindler 2017 (data reference is Plan 2012)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Painful distal diabetic neuropathy, placebo arm of three phase III trials.",
    dose_range     = "Placebo arm only (no active-drug exposure). Rescue acetaminophen recorded via the CONMED_PARA indicator.",
    regions        = "multinational (three pooled phase III trials).",
    biomarkers     = "Daily Likert pain score (0-10, 11-point scale) recorded in the evening. Median follow-up 124.5 days (range 2.5-125.5).",
    notes          = "Same 231-patient placebo-arm cohort as modellib('Plan_2012_pain') (Plan 2012 count model / DDMODEL00000194). Detailed baseline demographics are in Plan 2012 CPT."
  )

  ini({
    # ---------------------------------------------------------------------
    # mCTMM structural parameters (Schindler 2017 Table III).
    # 11 categories (0-10) -> alpha_1 estimated plus 9 negative b_k
    # increments (b_2 through b_10). Cumulative alphas per Eq. 11 are
    # alpha_k = alpha_1 + sum_{j=2..k} b_j.
    # ---------------------------------------------------------------------
    alpha1 <- 9.56;   label("Baseline cumulative logit for score >= 1 (unitless, logit scale)") # Table III alpha_1 = 9.56 (RSE 2.5%)
    b2     <- -1.90;  label("Increment b_2 (constrained negative)")                             # Table III b_2 = -1.90 (RSE 6.8%)
    b3     <- -2.48;  label("Increment b_3 (constrained negative)")                             # Table III b_3 = -2.48 (RSE 3.5%)
    b4     <- -1.90;  label("Increment b_4 (constrained negative)")                             # Table III b_4 = -1.90 (RSE 2.7%)
    b5     <- -1.68;  label("Increment b_5 (constrained negative)")                             # Table III b_5 = -1.68 (RSE 2.5%)
    b6     <- -1.66;  label("Increment b_6 (constrained negative)")                             # Table III b_6 = -1.66 (RSE 2.4%)
    b7     <- -1.64;  label("Increment b_7 (constrained negative)")                             # Table III b_7 = -1.64 (RSE 2.3%)
    b8     <- -1.76;  label("Increment b_8 (constrained negative)")                             # Table III b_8 = -1.76 (RSE 2.4%)
    b9     <- -2.24;  label("Increment b_9 (constrained negative)")                             # Table III b_9 = -2.24 (RSE 3.1%)
    b10    <- -2.07;  label("Increment b_10 (constrained negative)")                            # Table III b_10 = -2.07 (RSE 5.4%)

    # Mean equilibration time and its time-varying slope. Paper Table III
    # reports MET_0 in hours and SlopeMET in 1/week; the time-varying MET
    # is MET(t) = MET_0 * (1 + SlopeMET * t_weeks) (linear-fractional
    # form consistent with the paper's "linear increase in MET"
    # description and the 1/week unit on SlopeMET; the exact functional
    # form is documented in vignette Errata since the paper's Eq. 13
    # gives Pss but does not print MET(t) explicitly).
    lmet0    <- log(10.8);  label("Log baseline MET at t = 0 (h); Table III MET_0 = 10.8 h")   # Table III MET_0 = 10.8 h (RSE 7.6%)
    slopemet <- 0.740;      label("SlopeMET (1/week); fractional slope of the linear MET(t)")  # Table III SlopeMET = 0.740 week^-1 (RSE 22%)

    # Predictor slopes (Eq. 13). theta_time is per week; theta_ace is
    # unitless additive shift on the cumulative logit.
    theta_time <- -0.141; label("Slope of time (weeks) on the cumulative logits (Eq. 13)")           # Table III theta_time = -0.141 week^-1 (RSE 13%)
    theta_ace  <-  1.84;  label("Additive effect of concomitant acetaminophen use on the logits") # Table III theta_ace = 1.84 (RSE 4.7%)

    # ---------------------------------------------------------------------
    # Inter-individual variability (Table III footnote a: additive;
    # b: exponential). Omega values are standard deviations (Discussion
    # footnote); variances = SD^2. Subsequent pop parameters in the
    # cumulative logits (b_2 ... b_10) are wrapped in parentheses inside
    # model() to bypass rxode2's mu-reference syntax parser.
    # ---------------------------------------------------------------------
    etaalpha1     ~ 10.1761   # SD 3.19  -> var 10.176; additive BSV on alpha_1        (Table III omega_alpha_1 = 3.19, RSE 4.6%)
    etalmet0      ~  0.7762   # SD 0.881 -> var 0.7762; exponential BSV on MET_0       (Table III omega_MET0 = 0.881, RSE 6.8%)
    etatheta_time ~  0.0581   # SD 0.241 -> var 0.0581; additive BSV on theta_time     (Table III omega_theta_time = 0.241, RSE 4.9%)
    etaslopemet   ~  4.6225   # SD 2.15  -> var 4.6225; exponential BSV on SlopeMET    (Table III omega_SlopeMET = 2.15, RSE 7.7%)

    # ---------------------------------------------------------------------
    # Placeholder additive residual error on the expected-score output.
    # The Schindler 2017 fit uses NONMEM $ESTIMATION LAPLACE LIKE (the
    # likelihood IS the categorical mCTMM density). A placeholder
    # additive residual is attached for nlmixr2 compatibility.
    # ---------------------------------------------------------------------
    addSd_score <- 0.5
    label("Placeholder additive residual error on the expected-Likert-score output (score units); not from the source paper -- see vignette Assumptions")
  })

  model({
    # 1. Time conversion (rxode2 t is in hours; paper predictors use weeks).
    t_weeks <- t / 168

    # 2. Individual parameters. Each definition uses one mu-reference
    #    (pop + eta or pop * exp(eta)). Subsequent pop params (b_2 ...
    #    b_10) are wrapped in parentheses in step 4 to bypass rxode2's
    #    mu-reference syntax parser.
    alpha1_ind      <- alpha1     + etaalpha1
    theta_time_ind  <- theta_time + etatheta_time
    met0_ind        <- exp(lmet0  + etalmet0)
    slopemet_ind    <- slopemet   * exp(etaslopemet)

    # 4. Steady-state cumulative probabilities (Eq. 13); cumulative
    #    alphas from Eq. 11 constructed inline with subsequent pop
    #    params in parentheses (see step 2 comment).
    linpred  <- theta_time_ind * t_weeks + theta_ace * CONMED_PARA
    pge1_ss  <- expit(alpha1_ind                                                                        + linpred)
    pge2_ss  <- expit(alpha1_ind + (b2)                                                                 + linpred)
    pge3_ss  <- expit(alpha1_ind + (b2) + (b3)                                                          + linpred)
    pge4_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4)                                                   + linpred)
    pge5_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5)                                            + linpred)
    pge6_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5) + (b6)                                     + linpred)
    pge7_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5) + (b6) + (b7)                              + linpred)
    pge8_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5) + (b6) + (b7) + (b8)                       + linpred)
    pge9_ss  <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5) + (b6) + (b7) + (b8) + (b9)                + linpred)
    pge10_ss <- expit(alpha1_ind + (b2) + (b3) + (b4) + (b5) + (b6) + (b7) + (b8) + (b9) + (b10)        + linpred)

    # 5. Steady-state score probabilities Pss,k (Eq. 11).
    pss0  <- 1        - pge1_ss
    pss1  <- pge1_ss  - pge2_ss
    pss2  <- pge2_ss  - pge3_ss
    pss3  <- pge3_ss  - pge4_ss
    pss4  <- pge4_ss  - pge5_ss
    pss5  <- pge5_ss  - pge6_ss
    pss6  <- pge6_ss  - pge7_ss
    pss7  <- pge7_ss  - pge8_ss
    pss8  <- pge8_ss  - pge9_ss
    pss9  <- pge9_ss  - pge10_ss
    pss10 <- pge10_ss

    # 6. Time-varying MET (linear-fractional; see ini() comment on the form).
    met <- met0_ind * (1 + slopemet_ind * t_weeks)

    # 7. Transfer rate constants (Eqs. 8-9). Small eps guards against
    #    division by zero at extreme covariate values.
    eps    <- 1e-12
    lam01  <- pss1  / (met * (pss0 + pss1 + eps))
    lam10  <- pss0  / (met * (pss0 + pss1 + eps))
    lam12  <- pss2  / (met * (pss1 + pss2 + eps))
    lam21  <- pss1  / (met * (pss1 + pss2 + eps))
    lam23  <- pss3  / (met * (pss2 + pss3 + eps))
    lam32  <- pss2  / (met * (pss2 + pss3 + eps))
    lam34  <- pss4  / (met * (pss3 + pss4 + eps))
    lam43  <- pss3  / (met * (pss3 + pss4 + eps))
    lam45  <- pss5  / (met * (pss4 + pss5 + eps))
    lam54  <- pss4  / (met * (pss4 + pss5 + eps))
    lam56  <- pss6  / (met * (pss5 + pss6 + eps))
    lam65  <- pss5  / (met * (pss5 + pss6 + eps))
    lam67  <- pss7  / (met * (pss6 + pss7 + eps))
    lam76  <- pss6  / (met * (pss6 + pss7 + eps))
    lam78  <- pss8  / (met * (pss7 + pss8 + eps))
    lam87  <- pss7  / (met * (pss7 + pss8 + eps))
    lam89  <- pss9  / (met * (pss8 + pss9 + eps))
    lam98  <- pss8  / (met * (pss8 + pss9 + eps))
    lam9_10 <- pss10 / (met * (pss9 + pss10 + eps))
    lam10_9 <- pss9  / (met * (pss9 + pss10 + eps))

    # 8. Probability-compartment ODEs (Eq. 6). Boundary states 0 and 10
    #    have only one neighbour; interior states have two.
    pscore0(0)  <- pss0
    pscore1(0)  <- pss1
    pscore2(0)  <- pss2
    pscore3(0)  <- pss3
    pscore4(0)  <- pss4
    pscore5(0)  <- pss5
    pscore6(0)  <- pss6
    pscore7(0)  <- pss7
    pscore8(0)  <- pss8
    pscore9(0)  <- pss9
    pscore10(0) <- pss10

    d/dt(pscore0)  <- -lam01 * pscore0 + lam10 * pscore1
    d/dt(pscore1)  <-  lam01 * pscore0 - (lam10 + lam12) * pscore1 + lam21 * pscore2
    d/dt(pscore2)  <-  lam12 * pscore1 - (lam21 + lam23) * pscore2 + lam32 * pscore3
    d/dt(pscore3)  <-  lam23 * pscore2 - (lam32 + lam34) * pscore3 + lam43 * pscore4
    d/dt(pscore4)  <-  lam34 * pscore3 - (lam43 + lam45) * pscore4 + lam54 * pscore5
    d/dt(pscore5)  <-  lam45 * pscore4 - (lam54 + lam56) * pscore5 + lam65 * pscore6
    d/dt(pscore6)  <-  lam56 * pscore5 - (lam65 + lam67) * pscore6 + lam76 * pscore7
    d/dt(pscore7)  <-  lam67 * pscore6 - (lam76 + lam78) * pscore7 + lam87 * pscore8
    d/dt(pscore8)  <-  lam78 * pscore7 - (lam87 + lam89) * pscore8 + lam98 * pscore9
    d/dt(pscore9)  <-  lam89 * pscore8 - (lam98 + lam9_10) * pscore9 + lam10_9 * pscore10
    d/dt(pscore10) <-  lam9_10 * pscore9 - lam10_9 * pscore10

    # 9. Expected score (weighted sum). The observation name `score` uses
    #    the canonical generic-pain-score PD output (see
    #    inst/references/compartment-names.md); an alias `expected_score`
    #    is exposed for vignette convenience.
    score <- 0 * pscore0 + 1 * pscore1 + 2 * pscore2 + 3 * pscore3 +
             4 * pscore4 + 5 * pscore5 + 6 * pscore6 + 7 * pscore7 +
             8 * pscore8 + 9 * pscore9 + 10 * pscore10
    expected_score <- score

    # 10. Observation model (placeholder -- see vignette Assumptions).
    score ~ add(addSd_score)
  })
}
