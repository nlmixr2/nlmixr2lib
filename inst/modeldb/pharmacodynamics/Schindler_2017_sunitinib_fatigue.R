Schindler_2017_sunitinib_fatigue <- function() {
  description <- "Population PD model of NCI-CTCAE fatigue grade (0 / 1 / 2 / 3+) in adults with imatinib-resistant gastrointestinal stromal tumours (GIST) on sunitinib, using the minimal continuous-time Markov model (mCTMM) form of Schindler & Karlsson 2017. Four ordered-categorical states (grades 0, 1, 2, 3+) evolve via ODEs (Eq. 6) with transfer rate constants (Eqs. 8-9) parameterised through a single mean equilibration time (MET) and the steady-state score probabilities Pss,k. Steady-state cumulative probabilities Pss(Y >= k) follow a proportional-odds model (Eq. 10) driven by the relative change from baseline over time in soluble vascular endothelial growth factor receptor 3 (sVEGFR-3), a sunitinib pharmacodynamic biomarker. sVEGFR-3 is generated inline as a one-compartment indirect-response turnover with simple-Imax inhibition of Kin, driven by a per-cycle drug-exposure summary AUC = DOSE / CLI (following the same convention used in the sister Hansson_2013_sunitinib_hfs and Hansson_2013a_sunitinib model files); the sVEGFR-3 upstream parameters BAS_SVEGFR3, MRT_SVEGFR3, and EC50_SVEGFR3 are consumed as data covariates. The model has no PK ODE. Outputs are the four marginal state probabilities pscore0-pscore3, the cumulative probabilities pge1-pge3, and the expected grade."
  reference <- paste(
    "Schindler E, Karlsson MO.",
    "A Minimal Continuous-Time Markov Pharmacometric Model.",
    "AAPS J. 2017;19(5):1424-1435.",
    "doi:10.1208/s12248-017-0109-1.",
    "Companion DTMM form for the same endpoint / cohort:",
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO. PKPD modeling of predictors for adverse effects and",
    "overall survival in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "Sister mCTMM model file from the same paper:",
    "modellib('Schindler_2017_sunitinib_hfs').",
    "Upstream sVEGFR-3 biomarker dynamics adapted from",
    "Hansson EK et al. CPT Pharmacometrics Syst Pharmacol 2013;2:e84,",
    "doi:10.1038/psp.2013.61, DDMODEL00000197; see",
    "modellib('Hansson_2013a_sunitinib').",
    sep = " "
  )
  vignette <- "Schindler_2017_mCTMM"
  paper_specific_compartments <- c("svegfr3", "pscore0", "pscore1", "pscore2", "pscore3")

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "(NCI-CTCAE fatigue grade 0-3+, ordinal)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    svegfr3 = list(analyte = "soluble vascular endothelial growth factor receptor 3 (sVEGFR-3)", units = "mg", specimen = "plasma", verified = FALSE),
    pscore0 = list(analyte = "fatigue grade 0", units = "mg", specimen = "not applicable", verified = FALSE),
    pscore1 = list(analyte = "fatigue grade 1", units = "mg", specimen = "not applicable", verified = FALSE),
    pscore2 = list(analyte = "fatigue grade 2", units = "mg", specimen = "not applicable", verified = FALSE),
    pscore3 = list(analyte = "fatigue grade 3+", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles of a 4-weeks-on / 2-weeks-off schedule or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "For typical-cohort vignette simulations the value is held at 50 mg during on-cycles of a 4/2 schedule, matching the standard GIST regimen (Demetri 2006). Time-varying by design.",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject, time-fixed input. Required. Vignette virtual cohort uses the Hansson 2013a typical value of 32.819 L/h.",
      source_name        = "CL"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197).",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject, time-fixed input. Required. Used as the initial condition for the in-model svegfr3 state and as the denominator in the driver sVEGFR3rel(t) = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3. Hansson 2013a typical value 63,900 pg/mL.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject, time-fixed input. Used as kout3 = 1 / MRT_SVEGFR3 inside model(). Hansson 2013a typical value 401 h.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject, time-fixed input. Appears in the drug-effect term eff3 = auc / (EC50_SVEGFR3 + auc). Hansson 2013a typical IC50 is 1.0 mg*h/L.",
      source_name        = "EC53"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    n_observations = 55027L,
    age_range      = "adults with GIST",
    weight_range   = "not reported in Schindler 2017; upstream Hansson 2013 pooled four sunitinib studies",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Imatinib-resistant gastrointestinal stromal tumours (GIST) treated with placebo or sunitinib.",
    dose_range     = "sunitinib 25-75 mg PO QD on 4-weeks-on / 2-weeks-off, 2/2, 2/1 or continuous schedules; placebo run-in and active treatment.",
    regions        = "multinational (four pooled clinical trials).",
    biomarkers     = "NCI-CTCAE fatigue grade (v3.0), collected daily; grades 3 and 4 lumped into a single grade 3+ category because grade 4 was rare (~1%). Median follow-up 157 days (range 7-687).",
    notes          = "Same 303-patient GIST cohort as the sister Schindler_2017_sunitinib_hfs (HFS) mCTMM and Hansson_2013_sunitinib_hfs / _fatigue DTMM. The upstream sVEGFR-3 biomarker dynamics are consumed as data covariates (BAS_SVEGFR3 + MRT_SVEGFR3 + EC50_SVEGFR3) rather than re-fitted; see Hansson_2013a_sunitinib for the source popPD."
  )

  ini({
    # ---------------------------------------------------------------------
    # mCTMM structural parameters (Schindler 2017 Table II, Fatigue column).
    # Final estimates with %RSE parenthesised; source-trace per parameter.
    # ---------------------------------------------------------------------
    alpha1     <- -1.85;   label("Baseline cumulative logit for score >= 1 (unitless, logit scale)")                            # Table II Fatigue: alpha_1 = -1.85 (RSE 7.5%)
    b2         <- -1.08;   label("Increment b_2 such that alpha_2 = alpha_1 + b_2 (constrained negative)")                     # Table II Fatigue: b_2 = -1.08 (RSE 7.9%)
    b3         <- -1.92;   label("Increment b_3 such that alpha_3 = alpha_2 + b_3 (constrained negative)")                     # Table II Fatigue: b_3 = -1.92 (RSE 10%)
    lmet       <- log(909);  label("Log mean equilibration time (h) between two consecutive states; Eq. 7")                    # Table II Fatigue: MET = 909 h (= 37.9 days; RSE 11%)
    theta_svegfr3 <- -1.87;  label("Slope of the linear effect of sVEGFR-3 relative change on the cumulative logits (Eq. 12)") # Table II Fatigue: theta_sVEGFR3 = -1.87 (RSE 15%)

    # ---------------------------------------------------------------------
    # Inter-individual variability (Table II footnote a: additive;
    # b: exponential). The paper reports omega values as standard
    # deviations (Discussion footnote: "omega_P standard deviation of
    # the between-subject variability associated with parameter P");
    # variances = SD^2. Subsequent pop parameters in the cumulative
    # logits (b_2, b_3) are wrapped in parentheses inside model() to
    # keep rxode2 5.1.3's mu-reference syntax parser happy (it flags
    # bare `alpha_1 + b_2` in an mu-referenced context, but accepts
    # `alpha_1 + (b_2)`).
    # ---------------------------------------------------------------------
    etaalpha1        ~ 0.593   # SD 0.770 -> var 0.593; additive BSV on alpha_1  (Table II Fatigue omega_alpha1 = 0.770, RSE 19%)
    etatheta_svegfr3 ~ 1.1449  # SD 1.07  -> var 1.145; exponential BSV on theta_sVEGFR3 (Table II Fatigue omega_theta_sVEGFR3 = 1.07, RSE 7.7%)

    # ---------------------------------------------------------------------
    # Residual-error placeholder on the derived expected-grade output.
    # The Schindler 2017 fit uses NONMEM $ESTIMATION LAPLACE LIKE (the
    # likelihood IS the categorical mCTMM density); there is no
    # observation residual-error component in the source model. This
    # model file exposes a continuous expected_grade output (weighted
    # sum of the four state probabilities) for typical-value
    # visualisation; a placeholder additive residual error is attached
    # so nlmixr2 accepts the model. Not from the source -- see vignette.
    # ---------------------------------------------------------------------
    addSd_fatigue_grade <- 0.5
    label("Placeholder additive residual error on the expected-grade output (grade units); not from the source paper -- see vignette Assumptions")
  })

  model({
    # ------------------------------------------------------------
    # 1. Per-cycle drug-exposure summary (mg*h/L AUC).
    #    Matches the Hansson 2013a / 2013 (e84 / e85) idiom for
    #    driving sVEGFR-3 turnover from a daily-dose input.
    # ------------------------------------------------------------
    auc <- DOSE / CLI

    # ------------------------------------------------------------
    # 2. Simple-Imax drug effect on sVEGFR-3 (Hansson 2013a Eq. 3).
    # ------------------------------------------------------------
    eff3 <- auc / (EC50_SVEGFR3 + auc)

    # ------------------------------------------------------------
    # 3. sVEGFR-3 indirect-response turnover (Kin inhibition).
    #    kin3 * (1 - eff3) - kout3 * svegfr3, with baseline
    #    kin3 = BAS_SVEGFR3 * kout3 so that svegfr3(0) = BAS_SVEGFR3.
    # ------------------------------------------------------------
    kout3         <- 1 / MRT_SVEGFR3
    kin3          <- BAS_SVEGFR3 * kout3
    svegfr3(0)    <- BAS_SVEGFR3
    d/dt(svegfr3) <- kin3 * (1 - eff3) - kout3 * svegfr3

    # ------------------------------------------------------------
    # 4. sVEGFR-3 relative-change driver (drug effect on fatigue is
    #    driven by depletion of sVEGFR-3 from its baseline).
    # ------------------------------------------------------------
    sVEGFR3rel <- (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3

    # ------------------------------------------------------------
    # 5. Individual proportional-odds intercept and drug-effect slope.
    #    Additive BSV on alpha_1 is applied via a single mu-referenced
    #    variable `alpha1_ind`; higher-category cumulative intercepts
    #    inherit the same eta by construction (paper: "BSV terms are
    #    incorporated additively on alpha_1 so that the PO assumption
    #    is respected"). Exponential BSV on theta_sVEGFR3 keeps the
    #    sign of the negative pop estimate.
    # ------------------------------------------------------------
    alpha1_ind <- alpha1 + etaalpha1
    theta_sv   <- theta_svegfr3 * exp(etatheta_svegfr3)

    # ------------------------------------------------------------
    # 6. Steady-state cumulative probabilities (Eq. 12).
    #    Pss(Y >= k) = expit(alpha_k + theta_sVEGFR3 * sVEGFR3rel)
    #    with alpha_k = alpha_1 + sum_{j=2..k} b_j (Eq. 11).
    #    Subsequent pop params b_2 and b_3 are wrapped in parentheses
    #    to bypass rxode2's mu-reference syntax check (see ini()
    #    comment above).
    # ------------------------------------------------------------
    linpred_svr <- theta_sv * sVEGFR3rel
    pge1_ss <- expit(alpha1_ind                     + linpred_svr)
    pge2_ss <- expit(alpha1_ind + (b2)              + linpred_svr)
    pge3_ss <- expit(alpha1_ind + (b2) + (b3)       + linpred_svr)

    # ------------------------------------------------------------
    # 7. Steady-state score probabilities Pss,k (Eq. 11).
    # ------------------------------------------------------------
    pss0 <- 1       - pge1_ss
    pss1 <- pge1_ss - pge2_ss
    pss2 <- pge2_ss - pge3_ss
    pss3 <- pge3_ss

    # ------------------------------------------------------------
    # 8. Transfer rate constants (Eqs. 8 and 9). The constant-MET
    #    assumption combined with the steady-state balance
    #    lambda_{k,k+1} * Pss,k = lambda_{k+1,k} * Pss,{k+1} gives
    #      lambda_{k,k+1} = Pss,{k+1} / (MET * (Pss,k + Pss,{k+1}))
    #      lambda_{k+1,k} = Pss,k    / (MET * (Pss,k + Pss,{k+1}))
    #    A small floor eps prevents division by zero at extreme
    #    covariate values that could collapse a Pss,k pair to 0.
    # ------------------------------------------------------------
    met  <- exp(lmet)
    eps  <- 1e-12
    lam01 <- pss1 / (met * (pss0 + pss1 + eps))
    lam10 <- pss0 / (met * (pss0 + pss1 + eps))
    lam12 <- pss2 / (met * (pss1 + pss2 + eps))
    lam21 <- pss1 / (met * (pss1 + pss2 + eps))
    lam23 <- pss3 / (met * (pss2 + pss3 + eps))
    lam32 <- pss2 / (met * (pss2 + pss3 + eps))

    # ------------------------------------------------------------
    # 9. Probability-compartment ODEs (Eq. 6). Boundary cases at
    #    k = 0 and k = K only transfer to / from the neighbouring
    #    state; interior states have both neighbours. Compartment
    #    amounts sum to 1 at all times (invariant).
    # ------------------------------------------------------------
    pscore0(0) <- pss0
    pscore1(0) <- pss1
    pscore2(0) <- pss2
    pscore3(0) <- pss3

    d/dt(pscore0) <- -lam01 * pscore0 + lam10 * pscore1
    d/dt(pscore1) <-  lam01 * pscore0 - (lam10 + lam12) * pscore1 + lam21 * pscore2
    d/dt(pscore2) <-  lam12 * pscore1 - (lam21 + lam23) * pscore2 + lam32 * pscore3
    d/dt(pscore3) <-  lam23 * pscore2 - lam32 * pscore3

    # ------------------------------------------------------------
    # 10. Derived cumulative-probability outputs and expected grade.
    #     pgeK <- sum of state probabilities >= K.
    # ------------------------------------------------------------
    pge1 <- pscore1 + pscore2 + pscore3
    pge2 <- pscore2 + pscore3
    pge3 <- pscore3

    expected_grade <- 0 * pscore0 + 1 * pscore1 + 2 * pscore2 + 3 * pscore3
    fatigue_grade  <- expected_grade

    # ------------------------------------------------------------
    # 11. Observation model (placeholder -- see vignette Assumptions).
    # ------------------------------------------------------------
    fatigue_grade ~ add(addSd_fatigue_grade)
  })
}
