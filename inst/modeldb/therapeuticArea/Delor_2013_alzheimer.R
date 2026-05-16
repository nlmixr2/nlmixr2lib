Delor_2013_alzheimer <- function() {
  description <- "Disease-progression model (no drug input) for Alzheimer's disease (AD) progression on the Clinical Dementia Rating scale - Sum of Boxes (CDR-SOB, 0-18 score) over time, fit by Delor et al. (2013) to 2,700 CDR-SOB observations from 380 mild cognitive impairment (MCI) plus 180 AD patients in the Alzheimer's Disease Neuroimaging Initiative (ADNI) database with up to 4 years of follow-up. The model embeds an individual disease-onset time (DOT) and a logit-domain disease trajectory A(1) with a smooth-step activation (T^30 / (DOT^30 + T^30)) that turns the disease progression on at the subject's DOT; a two-component mixture model assigns each subject to either a fast-progression branch (rate plus accelerating term alpha*A(1)) or a slow-progression branch (rate alone, alpha = 0). A study-entry 'placebo' term PL*(1 - exp(-KPL*(t - T_ENTRY))) captures an early drop or delay observed after enrollment. Baseline CDR-SOB and ADAS-cog explain most of the DOT variability; baseline MMSE modifies alpha; baseline CDR-SOB, FAQ, and normalized hippocampal volume (RHPNM) modify the mixture probability. The published mixture probability and the residual-error scale-domain choice are clarified in the validation vignette's Assumptions and deviations section. No drug input is dosed in this model."
  reference <- paste(
    "Delor I, Charoin J-E, Gieschke R, Retout S, Jacqmin P;",
    "for the Alzheimer's Disease Neuroimaging Initiative. (2013).",
    "Modeling Alzheimer's Disease Progression Using Disease Onset Time",
    "and Disease Trajectory Concepts Applied to CDR-SOB Scores From ADNI.",
    "CPT Pharmacometrics Syst Pharmacol 2(10):e78.",
    "doi:10.1038/psp.2013.54.",
    sep = " "
  )
  vignette <- "Delor_2013_alzheimer"
  units <- list(
    time          = "year",
    dosing        = "(none; disease-progression model, no drug input)",
    concentration = "(CDR-SOB score, 0-18, unitless)"
  )

  covariateData <- list(
    CDR_SOB = list(
      description        = "Clinical Dementia Rating scale - Sum of Boxes score (0-18; integer-valued with half-unit increments in ADNI). Used here as a time-fixed baseline covariate (the source paper's CDR_bsl column). Drives both the per-subject disease-onset time (multiplicative power form, centred on CDR_SOB = 2) and the per-subject slow-progression-subpopulation probability (logit-additive shift, centred on CDR_SOB = 1).",
      units              = "(CDR-SOB units, 0-18 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the baseline value per subject in this model (source column CDR_bsl). Centring values come from the source paper (CDR_SOB/2 on DOT per Table 1 Part I; CDR_SOB - 1 on the mixture-logit per Table 1 Part II). Higher baseline score implies an earlier disease onset and a higher probability of being in the fast-progressing subpopulation.",
      source_name        = "CDR_bsl"
    ),
    ADAS_COG = list(
      description        = "Alzheimer's Disease Assessment Scale - cognitive subscale (ADAS-cog total 11) score. Used here as a time-fixed baseline covariate (the source paper's ADAS_bsl column). Drives the per-subject disease-onset time via a multiplicative power form centred on ADAS_COG = 12.67 (the dataset median per the source paper Table 1 Part I).",
      units              = "(ADAS-cog total 11 units, 0-70 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the baseline value per subject in this model (source column ADAS_bsl). The ADAS-cog total-11 form (not the modernized ADAS-cog 13) is used here per the source paper; this matches the ADNI total-11 score available in 2011. Centring value 12.67 is the dataset median.",
      source_name        = "ADAS_bsl"
    ),
    MMSE = list(
      description        = "Mini Mental State Examination score (0-30; higher = better). Used here as a time-fixed baseline covariate (the source paper's MMSE_bsl column). Modifies the per-subject disease-progression acceleration parameter alpha via a multiplicative power form centred on MMSE = 26 (Table 1 Part III).",
      units              = "(MMSE units, 0-30 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the baseline value per subject in this model (source column MMSE_bsl). Centring value 26 is the dataset median per the source.",
      source_name        = "MMSE_bsl"
    ),
    FAQ = list(
      description        = "Functional Assessment Questionnaire score (0-30 range commonly; higher = more functional impairment). Used here as a time-fixed baseline covariate (the source paper's FAQ_bsl column). Modifies the per-subject slow-progression-subpopulation probability via a logit-additive shift centred on FAQ = 1 (Table 1 Part II).",
      units              = "(FAQ units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the baseline value per subject in this model (source column FAQ_bsl). Centring value 1 is the dataset median per the source Table 1 Part II.",
      source_name        = "FAQ_bsl"
    ),
    RHPNM = list(
      description        = "Normalized hippocampal volume: the subject's average left+right hippocampal volume divided by the value expected for a healthy subject with the same age and estimated intracranial volume (head size). 1.0 corresponds to a healthy reference; values below 1.0 indicate atrophy. Used here as a time-fixed baseline covariate (the source paper's RHPNM column). Drives a logit-additive shift on the slow-progression-subpopulation probability centred on RHPNM = 1 (Table 1 Part II).",
      units              = "(unitless ratio; 1.0 = healthy reference)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at the baseline value per subject in this model. Derived in the source paper from the regression HPNMbsl_i = Age_i * (-26.6268 + EICVbsl_i * 0.0016 + 3340.4395) computed on healthy subjects, with RHPNMbsl_i = HIPVbsl_i / HPNMbsl_i. Source paper section 'ADNI data extraction and assembly'. Centring value 1 = healthy reference.",
      source_name        = "RHPNM"
    ),
    T_ENTRY = list(
      description        = "Per-subject time of study entry on the same global disease-time axis the model is integrated on (years). The dataset's TIME column represents global disease time (with TIME = 0 the integration origin and the activation function turning on around the subject's individual DOT). Observations occur at TIME = T_ENTRY to T_ENTRY + study_duration; the placebo term uses (TIME - T_ENTRY) as its clock.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Necessary because the source paper's disease-progression equation (1) uses a single time variable T anchored to a global disease-onset reference, while the placebo equation (2) uses time since study entry; for a typical-value reproduction the user must supply T_ENTRY per subject. A reasonable virtual-cohort choice is T_ENTRY = DOT_individual + a few years of established disease, so the patient is being observed during disease progression. See the validation vignette for the construction used to reproduce Figures 2-4.",
      source_name        = "(derived; not a directly observed column in ADNI -- in the source paper's NONMEM dataset, study time and disease-onset shift are handled internally)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 560L,
    n_studies      = 1L,
    age_range      = "ADNI cohort 55-91 years at enrollment (paper Methods: 'ADNI database consisted of 819 adults, aged 55-91 years'; the analysis dataset retained 380 MCI plus 180 AD with multiple CDR-SOB records).",
    age_median     = "(not reported in the source paper)",
    weight_range   = "(not a covariate in this disease-progression model and not tabulated in the source)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Alzheimer's disease (n = 180 at enrollment) and mild cognitive impairment (n = 380 at enrollment), with 166 MCI-to-AD conversions during the observation period in the wider ADNI download (6 AD-to-MCI and 19 NL-to-MCI, of which 3 subsequently AD, also recorded). Control subjects (normal cognition) had CDR-SOB scores at or near zero and were excluded from the modelled dataset because they did not carry information. The retained analysis cohort therefore mixes MCI and AD diagnoses; the mixture model partitions them into slow-progressing (about 44 percent of MCI and 3 percent of AD) and fast-progressing subpopulations.",
    dose_range     = "(not applicable; disease-progression model with no drug input)",
    regions        = "United States and Canada (ADNI sites; paper Methods: 'Subjects have been enrolled from over 50 sites across the USA and Canada').",
    notes          = "Subject counts and disease-state breakdown come directly from the source paper's Results section ('The final analysis dataset consisted of 2,700 CDR-SOB values collected from 380 MCI and 180 AD patients for up to 4 years'). Age range and regions come from the Methods 'ADNI data extraction and assembly' section. Sex distribution, ethnicity / race breakdown, and median ages are not reported in the publication tables and were not extracted. Biomarker covariates investigated in the source paper that are NOT included in the final model -- and therefore not in this file's covariateData list -- are: APOE-epsilon-4 carrier status, age (as a covariate on the model parameters), sex (as a covariate on the model parameters), serum cholesterol, Trail B test, MCI-to-AD conversion indicator, anticholinesterase co-medication, CSF tau / p-tau-181P / A-beta-1-42, plasma A-beta-1-42 / A-beta-1-40, ventricular volume, and unnormalised hippocampal volume; all were tested but did not enter the final model in Table 2. ADNI download date: 1 August 2011. Source paper acknowledges that CSF biomarkers (p-tau-181P, A-beta-1-42, tau) were only available in about half of the patients so their impact could not be fully explored."
  )

  ini({
    # Final-model parameter estimates from Delor 2013 Table 2 (right column,
    # 'Final model with covariates on $MIX'). Where the same parameter has a
    # value in both the base and final columns, the final column is used.

    # Structural disease-progression rates (in the logit-domain disease state).
    lrate_fast <- log(0.373); label("Log fast disease progression rate (1/year, on logit-domain CDR-SOB scale)")
    # Delor 2013 Table 2 final 'Fast disease progression rate, 1/year' = 0.373.

    lrate_slow <- log(0.260); label("Log slow disease progression rate (1/year, on logit-domain CDR-SOB scale)")
    # Delor 2013 Table 2 final 'Slow disease progression rate, 1/year' = 0.26.

    # Disease-onset time (DOT). Typical value reflects how many years ago the
    # disease started for a median patient (CDR_SOB_BL = 2, ADAS_COG_BL = 12.67).
    ldot <- log(16.1); label("Log typical disease onset time (years, global disease-time axis)")
    # Delor 2013 Table 2 final 'DOT, year' = 16.1.

    e_cdr_sob_dot <- -0.072; label("Power exponent for baseline CDR-SOB on DOT (DOT scales as (CDR_SOB/2)^e_cdr_sob_dot)")
    # Delor 2013 Table 2 final 'COV CDR-SOB on DOT' = -0.072. Centring 2 per
    # Table 1 Part I Start-model row.

    e_adas_cog_dot <- -0.0439; label("Power exponent for baseline ADAS-cog on DOT (DOT scales as (ADAS_COG/12.67)^e_adas_cog_dot)")
    # Delor 2013 Table 2 final 'COV ADAS on DOT' = -0.0439. Centring 12.67
    # per Table 1 Part I 'Base' row.

    # Multiplicative disease-progression acceleration term alpha (1/year),
    # active in the fast-progression branch only.
    lalpha <- log(0.0499); label("Log of typical multiplicative disease-progression acceleration alpha (1/year, fast branch only)")
    # Delor 2013 Table 2 final 'Multiplicative term (alpha)' = 0.0499.

    e_mmse_alpha <- -2.01; label("Power exponent for baseline MMSE on alpha (alpha scales as (MMSE/26)^e_mmse_alpha)")
    # Delor 2013 Table 2 final 'COV MMSE on alpha' = -2.01. Centring 26 per
    # Table 1 Part III 'Final' row. (The estimated effect was tested last by
    # the source authors via GAM-selected covariates; numerically unstable
    # alpha-covariate runs were excluded in Table 1 Part III.)

    # Baseline correction in the logit-to-CDR-SOB back-transform: predicted
    # CDR-SOB = 18 / (1 + exp(blc - A)). At A = 0 (pre-disease initial state)
    # this gives 18*exp(-4.48)/(1+exp(-4.48)) ~= 0.20, i.e. effectively zero,
    # which matches the source paper's structural assumption that CDR-SOB
    # starts near zero before disease onset.
    blc <- 4.48; label("Baseline correction shift in the logit-to-CDR-SOB back-transform (unitless)")
    # Delor 2013 Table 2 final 'Baseline correction (BLC)' = 4.48.

    # Placebo effect: an early drop or delay in CDR-SOB after study entry,
    # parameterised as CDR_observed = CDR_predicted - PL * (1 - exp(-KPL*(t-T_ENTRY))).
    lpl <- log(0.434); label("Log placebo-effect magnitude (CDR-SOB units; subtracted from predicted score after study entry)")
    # Delor 2013 Table 2 final 'Placebo magnitude (PL)' = 0.434.

    lkpl <- fixed(log(1)); label("Log placebo-effect equilibration rate constant (1/year; FIXED 1 in source -- IIV not supported by the data)")
    # Delor 2013 Table 2 final 'Placebo rate (KPL), 1/year' = 1 with no
    # estimated CV%, indicating FIXED. The source text 'estimation of IIV on
    # KPL was not supported by the data' confirms that the typical value
    # itself is also held constant in the final model.

    # Mixture-model intercept on the logit-scale slow-progression probability.
    # The source paper does NOT report the typical theta_5 in Table 2 for
    # either the base or final model; only the COV1/COV2/COV3 effects and the
    # post-fit MCI/AD classification fractions (44/3) are tabulated. The value
    # used here is the MCI-cohort slow fraction (0.44), which corresponds to
    # a typical-MCI patient (CDR_SOB_BL = 1, FAQ_BL = 1, RHPNM_BL = 1) -- the
    # MCI cohort is roughly the 'typical' patient at the centring values for
    # the three mixture covariates. See vignette Errata for the rationale and
    # for the alternative defaults (overall-cohort 0.31, AD-cohort 0.03)
    # implied by the source paper's reported MCI/AD slow breakdown.
    logit_p_slow <- log(0.44 / (1 - 0.44)); label("Logit of typical slow-progression-subpopulation probability at median covariate values (CDR_SOB = 1, FAQ = 1, RHPNM = 1)")
    # Delor 2013: typical theta_5 not directly tabulated; imputed from the
    # MCI slow fraction (44/3 per Table 2 final 'MCI/AD slow (alpha = 0)' row).
    # See vignette Assumptions and deviations / Errata.

    e_cdr_sob_slow <- -1.27; label("Logit-additive effect of (CDR_SOB - 1) on slow-subpopulation probability")
    # Delor 2013 Table 2 final 'COV CDR-SOB on $MIX' = -1.27.

    e_faq_slow <- -0.341; label("Logit-additive effect of (FAQ - 1) on slow-subpopulation probability")
    # Delor 2013 Table 2 final 'COV FAQ on $MIX' = -0.341.

    e_rhpnm_slow <- 7.5; label("Logit-additive effect of (RHPNM - 1) on slow-subpopulation probability")
    # Delor 2013 Table 2 final 'COV RHPNM on $MIX' = 7.5.

    # Inter-individual variability. The source paper reports IIV as a percent
    # CV (paper Table 2 column 'IIV (%)' with footnote b 'Interindividual
    # variability'); the variance on the log-normal scale is computed as
    # omega^2 = log(1 + (CV/100)^2). For DOT the reported CV is 1.5%, giving
    # an effectively negligible residual IIV (the source paper notes that
    # 'log normally distributed interindividual variability could only be
    # implemented on one of the two parameters, RATE or DOT, forcing either
    # the same RATE or the same DOT for every subject; the model was stable
    # when implementing IIV on parameter DOT'). For alpha and PL the CVs are
    # 64% and 86% respectively.
    etaldot   ~ 0.000225  # omega^2 = log(1 + 0.015^2) ~ 0.000225; Table 2 final IIV(DOT) = 1.5%.
    etalalpha ~ 0.3434    # omega^2 = log(1 + 0.64^2) ~ 0.3434;   Table 2 final IIV(alpha) = 64%.
    etalpl    ~ 0.5535    # omega^2 = log(1 + 0.86^2) ~ 0.5535;   Table 2 final IIV(PL) = 86%.

    # Residual error. The source paper performs the fit in the logit domain
    # with additive error on A(1) (paper Structural base model development:
    # 'modeling was performed in the logit domain and additive residual error
    # was used'); the reported SD is 0.377 on the logit-domain disease state.
    # In nlmixr2 / rxode2 the observation is most natural on the back-
    # transformed CDR-SOB scale, so the additive error is applied there.
    # This is a deviation from the source-paper residual model; the typical-
    # value trajectory is unaffected. See the vignette Assumptions and
    # deviations section for the deviation rationale.
    #
    # Only ONE observation endpoint is declared in the model() block (the
    # mixture-weighted prediction cdr_mix); the per-branch trajectories
    # cdr_fast and cdr_slow are exposed as derived model variables so a
    # downstream user can extract them from rxSolve output and compose
    # custom mixture-aware likelihoods. This single-endpoint encoding is
    # what allows rxSolve to attach observation records to the model
    # without requiring a per-row CMT / DVID specifier from the caller.
    addSd <- 0.377; label("Additive residual error standard deviation on the mixture-weighted CDR-SOB prediction (paper applies it on the logit-domain disease state; this model applies it on the back-transformed CDR-SOB scale)")
    # Delor 2013 Table 2 final 'RES Error (SD)' = 0.377.
  })

  model({
    # 1. Per-subject parameters. The covariate-effect forms reproduce the
    # source paper's Table 1 covariate relationships verbatim, with the
    # canonical covariate column names from inst/references/covariate-columns.md.
    dot_indiv <- exp(ldot + etaldot) *
                 (CDR_SOB / 2)^e_cdr_sob_dot *
                 (ADAS_COG / 12.67)^e_adas_cog_dot
    alpha_indiv <- exp(lalpha + etalalpha) *
                   (MMSE / 26)^e_mmse_alpha
    pl_indiv  <- exp(lpl + etalpl)
    kpl_indiv <- exp(lkpl)
    rate_fast <- exp(lrate_fast)
    rate_slow <- exp(lrate_slow)

    # 2. Mixture-probability of being in the slow-progression subpopulation
    # for this subject. The mixture intercept logit_p_slow plus three logit-
    # additive covariate shifts (Table 1 Part II structure) determine the
    # per-subject probability. Exposed as model variables so a downstream
    # user can simulate from both branches and combine the predictions with
    # the source-paper weights (see vignette).
    logit_p_slow_indiv <- logit_p_slow +
                          e_cdr_sob_slow * (CDR_SOB - 1) +
                          e_faq_slow     * (FAQ     - 1) +
                          e_rhpnm_slow   * (RHPNM   - 1)
    p_slow_indiv <- 1 / (1 + exp(-logit_p_slow_indiv))
    p_fast_indiv <- 1 - p_slow_indiv

    # 3. Smooth-step disease-activation function. The exponent 30 follows the
    # source paper exactly (T^30 / (DOT^30 + T^30) approximates a step from 0
    # to 1 around T = DOT); the function is ~0 for time << DOT, 0.5 at
    # time = DOT, and ~1 for time >> DOT.
    activation <- time^30 / (dot_indiv^30 + time^30)

    # 4. Disease-trajectory ODE in the logit domain. Two state variables
    # carry the fast and slow branches; both initialised at A = 0 (pre-
    # disease, mapping to CDR-SOB ~= exp(-blc)/(1+exp(-blc))*18 ~= 0.2 after
    # back-transformation, which approximates the floor of the score).
    d/dt(a_fast) <- (rate_fast + alpha_indiv * a_fast) * activation
    d/dt(a_slow) <-  rate_slow                          * activation

    # 5. Back-transformation from the logit-domain disease state to the
    # bounded 0-18 CDR-SOB scale.
    cdr_fast_predisease <- 18 / (1 + exp(blc - a_fast))
    cdr_slow_predisease <- 18 / (1 + exp(blc - a_slow))

    # 6. Study-entry 'placebo' subtraction. The clock is (time - T_ENTRY),
    # zeroed before study entry so the term is exactly zero in the pre-
    # observation simulation window. After study entry the term grows
    # asymptotically toward pl_indiv as t -> infinity at rate kpl_indiv.
    t_pl_raw    <- time - T_ENTRY
    post_entry  <- t_pl_raw > 0
    t_pl        <- t_pl_raw * post_entry
    placebo_term <- pl_indiv * (1 - exp(-kpl_indiv * t_pl)) * post_entry

    cdr_fast <- cdr_fast_predisease - placebo_term
    cdr_slow <- cdr_slow_predisease - placebo_term

    # 7. Mixture-weighted prediction. The user obtains both branches and the
    # mixture probability from rxSolve output (cdr_fast, cdr_slow, and
    # p_slow_indiv are all derived model variables and are emitted in the
    # rxSolve data frame); cdr_mix below is the typical-value expectation
    # across the mixture and is the single declared observation endpoint
    # for the model so that rxSolve can attach observation records without
    # requiring a per-row CMT / DVID specifier.
    cdr_mix <- p_slow_indiv * cdr_slow + p_fast_indiv * cdr_fast

    cdr_mix ~ add(addSd)
  })
}
