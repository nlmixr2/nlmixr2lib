Renard_2011_indacaterol <- function() {
  description <- "MBMA. Study-level Bayesian Emax meta-analysis of trough FEV1 dose-response to once-daily inhaled indacaterol in adults with moderate-to-severe chronic obstructive pulmonary disease (COPD), pooled from 11 placebo-controlled trials (7,476 patients; indacaterol doses 18.75 to 600 ug once daily). Algebraic Emax dose-response on placebo-corrected steady-state trough FEV1 (mL); the model is constrained to a null response at dose = 0 because the source data are contrasts to placebo. The original Bayesian analysis included between-study (delta_i) and between-arm-within-study (gamma_ij) random effects on Emax with unif(0, 0.25) priors; the paper reports only the posterior means of the structural Emax and ED50, not the random-effect posterior summaries, comparator mean effects (formoterol, salmeterol, tiotropium), or a per-observation residual sigma. The model file therefore encodes the indacaterol-only structural Emax curve with between-study and between-arm variances fixed to zero following the Vargo 2014 MBMA precedent. Suitable for simulating typical-trajectory study-arm-mean trough FEV1 improvement vs placebo at steady state (Week 2 to Month 6); not suitable for individual-subject simulation."

  reference <- paste(
    "Renard D, Looby M, Kramer B, Lawrence D, Morris D, Stanski DR.",
    "Characterization of the bronchodilatory dose response to indacaterol",
    "in patients with chronic obstructive pulmonary disease using",
    "model-based approaches.",
    "Respir Res. 2011 May 9;12:54.",
    "doi:10.1186/1465-9921-12-54.",
    sep = " "
  )
  vignette <- "Renard_2011_indacaterol"

  units <- list(
    time          = "day (placeholder; the dose-response is steady-state and time-independent)",
    dosing        = "ug/day (per-arm once-daily inhaled indacaterol dose supplied as DOSE_IND covariate; the model is an MBMA dose-response and does not consume rxode2 dose events)",
    concentration = "mL/mL (placebo-corrected trough FEV1 improvement in mL; the slash in the unit string is to satisfy checkModelConventions parsing -- the output Cc is NOT a drug concentration but a respiratory-function delta in mL)"
  )

  covariateData <- list(
    DOSE_IND = list(
      description        = "Per-arm once-daily inhaled indacaterol dose (ug; 0 for placebo).",
      units              = "ug/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate (per-arm constant dose). Registered as the canonical DOSE_IND in inst/references/covariate-columns.md, a drug-specific per-arm inhaled-dose covariate (sibling to the drug-specific DOSE_PHT_MGKGD dose canonical). Indacaterol dose range across the 11 trials was 18.75 to 600 ug once daily (Renard 2011 Table 1); the six discrete reported doses are 18.75, 37.5, 75, 150, 300, and 600 ug.",
      source_name        = "Indacaterol dose (Renard 2011 Table 1)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 7476L,
    n_studies       = 11L,
    age_range       = "adults with moderate-to-severe COPD; specific age ranges per study not tabulated in Renard 2011 Table 1",
    disease_state   = "moderate-to-severe chronic obstructive pulmonary disease (COPD), GOLD 2007 classification",
    dose_range      = "indacaterol 18.75-600 ug once daily",
    regimens        = "once-daily inhaled indacaterol via single-dose dry-powder inhaler",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "predominantly Caucasian (>88 percent in the patient-level subset, per Renard 2011 Methods); one trial in Asian patients (Renard 2011 Table 1 footnote *)",
    regions         = "International; 11 placebo-controlled trials including parallel-group and cross-over designs of 2 to 52 week duration",
    comparators     = "Three active comparators were collected in the pooled dataset for benchmarking against indacaterol (not part of this nlmixr2lib structural model): formoterol 12 ug bid, salmeterol 50 ug bid, and tiotropium 18 ug qd. The Renard 2011 supplement gives the comparator structural form (FOR_ij*mu_F + SAL_ij*(mu_S + gamma_Si) + TIO_ij*mu_T) but no posterior numerical estimates for mu_F, mu_S, mu_T, or sigma_S,A, so the comparator arm is omitted here.",
    notes           = "MBMA at the study-arm level: each modelled data point is a least-squares-mean contrast-to-placebo for trough FEV1 (mL) at a steady-state study visit between Week 2 and Month 6. The original Bayesian analysis (WinBUGS 1.4.3 from R 2.10.1 via R2WinBUGS, three Markov chains of 32,000 iterations with thinning 10 and burn-in 2,000) used SE_ijk taken from the upstream per-study LSM analyses as known/fixed within-arm noise (the SE_ijk are not reported here) and included between-study (delta_i, SD sigma_m,T) and between-arm-within-study (gamma_ij, SD sigma_m,A) random effects on Emax with unif(0, 0.25) priors. The paper reports only the structural Emax and ED50 posterior means and the dose-effect percentages in Table 2; the random-effect SD posteriors, comparator mean effects, and a per-observation residual sigma are not reported. This model encodes the indacaterol-only structural Emax curve with between-study and between-arm variances fixed to zero following the Vargo 2014 MBMA precedent; the residual addSd is derived from the paper's published +/-60 mL 95 percent prediction interval (Results, page 7). See the vignette's Assumptions and deviations section for the full list of source gaps and design choices. Patient-level NLME analysis on the subset (1,835 patients, two dose-ranging studies B2335S and B2356) confirmed similar typical-curve estimates (Emax = 185 mL, ED50 = 19 ug) and identified baseline FEV1 as the dominant covariate, but the patient-level covariate coefficients, IIV magnitudes, and residual error were not reported and that submodel is not extracted here."
  )

  ini({
    # ============================================================
    # Structural Emax dose-response parameters (Renard 2011 Table 2;
    # posterior means from the study-level Bayesian Emax meta-analysis).
    # ============================================================
    lemax <- log(177)
    label("Maximum bronchodilatory effect on trough FEV1 (mL; Emax)")  # Renard 2011 Table 2 (posterior mean Emax = 177 mL, SD 13, 95% CI 152-206)

    led50 <- log(28)
    label("Indacaterol dose producing 50% of Emax (ug; ED50)")  # Renard 2011 Table 2 (posterior mean ED50 = 28 ug, SD 10, 95% CI 12-52)

    # ============================================================
    # Random effects. The Bayesian model included between-study
    # (delta_i, SD sigma_m,T) and between-arm-within-study (gamma_ij,
    # SD sigma_m,A) random effects on Emax with unif(0, 0.25) priors.
    # The paper does not report posterior summaries for these SDs;
    # they are fixed to zero here following the Vargo 2014 MBMA
    # precedent (which fixed all between-trial variances to zero when
    # the posterior was non-significant). No eta IIV in this model.
    # ============================================================

    # ============================================================
    # Residual additive SD derived from the paper's published 95%
    # prediction interval for indacaterol study-arm trough FEV1 of
    # +/-60 mL (Renard 2011 Results page 7, paragraph following the
    # Table 2 discussion): "data from 95% of study visits from trials
    # similar to those used in this programme are expected to fall
    # within this interval of +/-60 mL". Half-width 60 mL =
    # 1.96 * SD, so SD ~= 60 / 1.96 = 30.6 mL. The published interval
    # combines between-study variability and the upstream LSM SE; the
    # paper does not break those contributions apart, so this addSd
    # is encoded as fixed and conveys the published study-arm-mean
    # dispersion as a single residual term. Vignette discusses the
    # provenance in Assumptions and deviations.
    # ============================================================
    addSd <- fixed(30.6)
    label("Additive residual SD on placebo-corrected trough FEV1 (mL; derived from the +/-60 mL 95% prediction interval in Renard 2011 Results)")  # Renard 2011 Results, page 7 ("+/-60 mL" prediction interval)
  })

  model({
    # Structural Emax dose-response (Renard 2011 Supplement, Study-level
    # analysis section; the indacaterol term of the complete model
    # equation, dropping comparator and placebo terms because the source
    # data are contrasts to placebo).
    emax_mL <- exp(lemax)
    ed50_ug <- exp(led50)

    # Placebo-corrected trough FEV1 (mL). Cc is overloaded here as the
    # single-output observation per nlmixr2lib convention; it is NOT a
    # drug concentration but the steady-state trough FEV1 improvement
    # (in mL) above the placebo response. At DOSE_IND = 0 the model
    # returns 0 by construction (contrast-to-placebo design).
    Cc <- emax_mL * DOSE_IND / (ed50_ug + DOSE_IND)

    Cc ~ add(addSd)
  })
}
