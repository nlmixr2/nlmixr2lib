Mercier_2014_tramadol_tapentadol_mbma <- function() {
  description <- "MBMA. Longitudinal model-based meta-analysis pain-intensity time-course model comparing tramadol and tapentadol in adults with chronic non-malignant pain (osteoarthritis, back pain, neuropathic pain, and other chronic non-malignant pain), fit to arm-level summary data from 45 double-blind Phase II/III randomized clinical trials representing 81 treatment arms and approximately 12,985 patients. Pain intensity on a normalized 0-10 scale is described by a logistic (expit) latent-scale model: PI(t) = 10 * expit(Base + R * (1 - exp(-k*t))), where the extent-of-reduction term R has an Emax-in-dose form for tramadol: R = R_Pla * (1 + theta_Base * logit(PI0/10) + theta_Trm * Dose_trm / (Dose_trm + ED50) * TRAMADOL + theta_Tap * TAPENTADOL). Between-study-arm variability enters additively on Base (SD 0.313 on the latent scale) and exponentially on R (SD 0.065). The residual variance is inversely proportional to per-arm sample size N; the paper reports two scale-dependent residual SDs (VAS/continuous SD 0.260 and categorical SD 0.205) and this file exposes the dominant VAS SD as the primary residual, with the categorical SD documented in ini() and in the vignette Assumptions/Errata. Suitable simulation scope is study-arm-mean pain-intensity time-course over 0-15 weeks (one 52-week trial contributed but the paper's simulations use 12 weeks); NOT individual-patient predictions. Tramadol was studied over a wide dose range (adequate to characterize a dose-response); tapentadol was only studied over 100-250 mg bid (no dose-response estimable), so the tapentadol effect is a single per-arm indicator effect. The paper also fit separate logistic MBMA sub-models for adverse events (constipation, nausea, vomiting, dizziness, somnolence) and drop-outs (due to adverse event, lack of efficacy); those sub-models are reported only as graphical odds ratios in Figs 3-4 with no tabulated logistic intercepts or slopes, so they are NOT extracted here (see vignette Errata for the omission)."

  reference <- paste(
    "Mercier F, Claret L, Prins K, Bruno R.",
    "A Model-Based Meta-analysis to Compare Efficacy and Tolerability of",
    "Tramadol and Tapentadol for the Treatment of Chronic Non-Malignant Pain.",
    "Pain Ther. 2014 Jun;3(1):31-44.",
    "doi:10.1007/s40122-014-0023-5.",
    sep = " "
  )
  vignette <- "Mercier_2014_tramadol_tapentadol"

  units <- list(
    time          = "week (time since randomization; the paper reports k in 1/week)",
    dosing        = "mg/day (per-arm daily tramadol dose supplied as CONMED_TRAMADOL_DOSE covariate; this MBMA does not consume rxode2 dose events)",
    concentration = "score/score (arm-mean pain intensity on the 0-10 normalized scale; the observation `score` is NOT a drug concentration; the slash in the unit string satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    PAIN = list(
      description        = "Per-arm mean baseline pain intensity on the 0-10 normalized scale used across the pooled Mercier 2014 trials (converted per Table 1 from 13+ raw scales including VAS 0-100 mm, VAS 0-10 cm, various Likert scales, BS11 11-point box score, NAS 0-100 mm, and 0-3/0-4/1-5 categorical scales; all rescaled linearly to 0-10 per the paper's conversion rules).",
      units              = "score (0-10)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-arm baseline (pre-treatment). Enters the R equation as logit(PAIN/10) (see paper Table 3 covariate effect theta_Base = 0.158). Note the SCALE: the canonical PAIN register entry (inst/references/covariate-columns.md) documents PAIN in mm on a 0-100 VAS scale; this MBMA uses arm-level values already converted to the 0-10 normalized scale, so downstream users assembling covariate data for this model must supply PAIN on the 0-10 scale (multiply mm VAS values by 0.1). Mercier 2014 Results reports the mean baseline pain intensity across studies was 6.9 (SD 0.72). PAIN can equal 0 or 10 in edge cases; the logit(PAIN/10) form maps 0 -> -Inf and 10 -> +Inf, so simulation code should clamp PAIN to (0 + eps, 10 - eps) if extreme values are supplied.",
      source_name        = "PI0 (Mercier 2014 Eq. R and Table 3 theta_Base row)"
    ),
    TRAMADOL = list(
      description        = "Binary study-arm treatment indicator: 1 = the arm received tramadol, 0 otherwise (placebo, tapentadol, or oxycodone active-control arm).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (arm did not receive tramadol)",
      notes              = "MBMA study-arm-level treatment indicator (a property of the trial arm, not of an individual patient). Documented inline per the Vargo_2014_statins_ezetimibe_mbma / Sadouki_2025 in-file-documentation precedent for MBMA drug-arm indicators. Enters the extent-of-reduction R via the tramadol Emax-in-dose term theta_Trm * Dose_trm / (Dose_trm + ED50) * TRAMADOL. The tramadol dose (mg/day; qd and bid arms have been rescaled to the total daily dose in the CONMED_TRAMADOL_DOSE column). Of the 81 study arms in the Mercier 2014 database, 43 were tramadol arms (30 osteoarthritis + 7 back pain + 2 neuropathic pain + 4 miscellaneous per Table 2).",
      source_name        = "Tramadol arm indicator (Mercier 2014 Eq. R)"
    ),
    TAPENTADOL = list(
      description        = "Binary study-arm treatment indicator: 1 = the arm received tapentadol, 0 otherwise (placebo, tramadol, or oxycodone active-control arm).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (arm did not receive tapentadol)",
      notes              = "MBMA study-arm-level treatment indicator. Documented inline per the Vargo_2014_statins_ezetimibe_mbma / Sadouki_2025 precedent. Enters R additively via theta_Tap * TAPENTADOL. Mercier 2014 could not fit a tapentadol dose-response because tapentadol was studied only over the narrow 100-250 mg bid range; the effect is therefore encoded as a fixed per-arm additive term rather than an Emax-in-dose form. Of the 81 study arms, 8 were tapentadol arms (5 osteoarthritis + 2 back pain + 1 neuropathic pain per Table 2). Simulations with TAPENTADOL = 1 should restrict the tapentadol daily dose to the 100-250 mg bid range that constitutes the domain of validity per the paper Discussion.",
      source_name        = "Tapentadol arm indicator (Mercier 2014 Eq. R)"
    ),
    CONMED_TRAMADOL_DOSE = list(
      description        = "Per-arm daily tramadol dose (mg/day; 0 for arms that did not receive tramadol).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Documented inline per the Vargo_2014_statins_ezetimibe_mbma CONMED_<INN>_DOSE precedent. Enters the tramadol dose-response as theta_Trm * Dose_trm / (Dose_trm + ED50) * TRAMADOL. The Mercier 2014 database included tramadol arms at doses from ~100 mg/day up to ~400 mg/day (the paper's Discussion cites tramadol 300 mg qd as a typical dose for the tramadol-vs-tapentadol simulation comparison); the estimated ED50 was 184 mg (SE 66). Set to 0 for non-tramadol arms; the TRAMADOL indicator handles the arm-selection.",
      source_name        = "Tramadol daily dose (Mercier 2014 Table 3 ED50 row and Discussion)"
    )
  )

  covariatesDataExcluded <- list(
    SCALE_CATEGORICAL = list(
      description        = "Study-arm scale-type indicator: 1 = arm-mean pain intensity was reported on a categorical scale (Likert / box-score / 0-3 / 0-4 / 1-5 / 16-point), 0 = arm-mean was reported on a continuous VAS scale.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (continuous VAS scale)",
      notes              = "Mercier 2014 Methods reports two scale-dependent residual variances (VAS/continuous sigma_1 = 0.260 and categorical sigma_2 = 0.205); the scale-type indicator selects between them. This model uses the dominant VAS residual as the encoded addSd; SCALE_CATEGORICAL is documented in covariatesDataExcluded so that users assembling arm-level data can preserve the scale-type provenance without triggering an 'unused covariate' checkModelConventions warning. See the vignette Assumptions and Deviations section for the dominant-residual rationale."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12985L,
    n_studies      = 45L,
    n_arms         = 81L,
    n_observations = 534L,
    age_range      = "median 58 years (range across trial medians 47-72 years)",
    weight_range   = "not reported at arm level",
    sex_female_pct = 64,
    race_ethnicity = "not reported",
    disease_state  = "adults with chronic non-malignant pain: osteoarthritis (63.2% of arms), back pain (22.6%), neuropathic pain (8.5%), and other chronic non-malignant pain (5.7% -- rheumatoid arthritis, non-cancer chronic pain, fibromyalgia)",
    dose_range     = "tramadol arms ranged from ~100 to ~400 mg/day (paper's typical-dose simulation uses 300 mg qd); tapentadol arms 100-250 mg bid; placebo arms carried no drug",
    regimens       = "tramadol qd (once-daily controlled-release) and bid/tid (immediate-release) rescaled to total daily dose; tapentadol 100-250 mg bid",
    treatments     = "placebo, tramadol, tapentadol (six of the tapentadol trials were active-controlled against oxycodone; only the tramadol/tapentadol/placebo arms are modeled)",
    timepoints     = "arm-mean pain intensity reported at multiple timepoints per trial; median follow-up 9.0 weeks (SD 6.8), range up to 52 weeks (Wild et al. 2010)",
    regions        = "not reported at arm level",
    notes          = "MBMA at the study-arm level: each modeled data point is the arm-mean pain intensity in one trial arm at one timepoint, weighted by arm sample size N (residual variance scales as sigma^2 / N per Mercier 2014 Methods). Total 534 arm-timepoint observations pooled across 81 arms and 45 trials. The model is intended for simulating arm-mean pain-intensity time-courses and is NOT suitable for individual-subject simulation. Adverse-event and drop-out logistic sub-models are described in the paper but their coefficient values are reported only as graphical odds ratios in Figures 3-4 with no tabulated intercepts or slopes; those sub-models are NOT extracted -- see vignette Errata for the omission and the digitized Fig 3/4 odds-ratio ranges."
  )

  ini({
    # ============================================================
    # Latent-scale logistic pain-intensity model (Mercier 2014 Eq. 1):
    #   PI_ijk = 10 * expit(Base_ij + R_ij * (1 - exp(-k*t))) + e_ijk
    # with expit(x) = exp(x)/(1 + exp(x)).
    # All values are from Mercier 2014 Table 3 (Parameter estimates of
    # the final pain intensity model; page 38 / Estimate (SE) column).
    # ============================================================

    e0 <- 0.812
    label("Typical latent-scale baseline pain intensity (paper: Base). The observed 0-10 baseline is 10*expit(e0); at e0 = 0.812 that is 6.92, matching the paper's reported mean baseline 6.9 (Results, first Pain Intensity paragraph).")  # Mercier 2014 Table 3 Base = 0.812 (SE 0.053)

    emax_pla <- -0.819
    label("Placebo extent-of-reduction on the latent scale (paper: R_Pla). Negative because pain intensity decreases from baseline under placebo. Enters R additively before the drug and covariate modifiers.")  # Mercier 2014 Table 3 RPla = -0.819 (SE 0.058)

    e_pain_emax <- 0.158
    label("Baseline-pain-intensity covariate coefficient on R (paper: theta_Base). Enters as theta_Base * logit(PAIN/10); positive value means subjects with higher baseline pain intensity have a greater reduction in pain intensity (Mercier 2014 Discussion; expected relationship also noted in Katz 2008 ref [27]).")  # Mercier 2014 Table 3 thetaBase = 0.158 (SE 0.111)

    e_tramadol_emax <- 0.980
    label("Tramadol Emax-in-dose amplitude on R (paper: theta_Trm). Enters R as theta_Trm * CONMED_TRAMADOL_DOSE / (CONMED_TRAMADOL_DOSE + ED50) * TRAMADOL. Together with theta_Base and TRAMADOL = 1 at Dose = 300 mg/day this gives the paper's 46% CFB reduction for tramadol 300 mg qd at week 12 assuming PAIN = 6.9.")  # Mercier 2014 Table 3 thetaTrm = 0.980 (SE 0.144)

    e_tapentadol_emax <- 0.259
    label("Tapentadol per-arm additive amplitude on R (paper: theta_Tap). Enters R as theta_Tap * TAPENTADOL; no dose-response was estimable for tapentadol because it was studied only over 100-250 mg bid (Mercier 2014 Results). Domain of validity: 100-250 mg bid.")  # Mercier 2014 Table 3 thetaTap = 0.259 (SE 0.018)

    led50_tramadol <- log(184)
    label("Log tramadol ED50 (mg/day). Back-transformed value 184 mg/day is the dose at which the tramadol Emax-in-dose term reaches half its maximum amplitude.")  # Mercier 2014 Table 3 ED50 = 184 mg (SE 66)

    lkel <- log(0.571)
    label("Log first-order approach-to-plateau rate constant (paper: k, 1/week). The onset of the pain-intensity time course; 50% of the maximum effect is reached at ln(2)/k = 1.21 weeks (paper Results: t1/2 = 1.2 weeks). The paper's Methods paragraph describes drug-specific rates kdrug = kpbo - kDdrug but reports a single common k in Table 3 because 'the onset of effect was found to be as fast in the active groups (tapentadol and tramadol) as in placebo' (Results); this file therefore uses a single k for all treatment groups.")  # Mercier 2014 Table 3 k = 0.571 /week (SE 0.015)

    # ============================================================
    # Between-study-arm random effects (Mercier 2014 Methods, second
    # Pain Intensity paragraph): additive on Base (eta_Base ~ N(0,
    # omega_Base^2)) and exponential on R (eta_R ~ N(0, omega_R^2)).
    # Table 3 reports omega values (standard deviations, per nlme
    # reporting convention); the ini() values below are the variances
    # (omega^2). Encoded as MBMA study-level etas (NOT individual
    # between-subject variability) per SKILL Phase-1 Step-3a MBMA
    # guidance.
    # ============================================================
    eta_study_e0   ~ 0.098     # Mercier 2014 Table 3 omega_Base = 0.313 (SD); variance = 0.313^2 = 0.098
    eta_study_emax ~ 0.004225  # Mercier 2014 Table 3 omega_RPla = 0.065 (SD); variance = 0.065^2 = 0.004225

    # ============================================================
    # Residual error. Mercier 2014 Methods reports the residual as
    # e_ijk ~ N(0, sigma_res^2 / N_ijk) where N_ijk is the per-arm
    # sample size (residual variance inversely proportional to arm
    # size). Two scale-dependent variances are reported:
    #   sigma_1 = 0.260 for VAS / continuous scales
    #   sigma_2 = 0.205 for categorical scales
    # This file exposes the dominant VAS residual (the majority of
    # arms) as addSd; sigma_2 is documented alongside for source
    # trace. Per-arm sample-size weighting (dividing by sqrt(N)) is
    # left to downstream simulation code (same pattern as
    # Boucher_2018_naproxen_mbma and Vargo_2014_statins_ezetimibe_mbma).
    # ============================================================
    addSd <- 0.260
    label("Residual SD on the 0-10 pain-intensity scale (Mercier 2014 Table 3 sigma_1 for VAS/continuous scales). Per-observation SD is addSd / sqrt(N_arm) per Methods; N weighting is applied downstream. Categorical-scale sigma_2 = 0.205 is reported by the paper for arms measured on non-VAS categorical scales; this file uses sigma_1 as the primary residual and documents sigma_2 in vignette Errata.")  # Mercier 2014 Table 3 sigma_1 = 0.260
  })

  model({
    # Study-arm covariates supplied per row:
    #   PAIN                    -- per-arm baseline pain intensity on 0-10 scale
    #   TRAMADOL                -- 1 if the arm received tramadol, 0 otherwise
    #   TAPENTADOL              -- 1 if the arm received tapentadol, 0 otherwise
    #   CONMED_TRAMADOL_DOSE    -- per-arm daily tramadol dose (mg/day; 0 for non-tramadol arms)

    # 1. Back-transform log-scale structural parameters.
    kel  <- exp(lkel)                    # 1/week; paper's k
    ed50 <- exp(led50_tramadol)          # mg/day; paper's ED50

    # 2. Baseline pain-intensity covariate term inside R.
    # Paper Eq. R uses logit(PI0/10); PAIN is on 0-10 scale in this MBMA
    # so PAIN/10 is in (0, 1). PAIN = 6.9 -> logit(0.69) = 0.7985.
    logit_pain <- log((PAIN / 10) / (1 - PAIN / 10))

    # 3. Individual (study-arm-level) baseline on the latent scale
    # (paper: Base_ij = Base + eta_Base_ij, additive random effect).
    base_arm <- e0 + eta_study_e0

    # 4. Extent-of-reduction R on the latent scale (paper's Eq. R).
    # Drug-selection is via the TRAMADOL and TAPENTADOL binary
    # indicators (mutually exclusive in the fitted dataset); the
    # eta_study_emax random effect enters exponentially per paper Methods.
    r_arm_raw <- emax_pla * (
      1 +
        e_pain_emax * logit_pain +
        e_tramadol_emax * (CONMED_TRAMADOL_DOSE / (CONMED_TRAMADOL_DOSE + ed50)) * TRAMADOL +
        e_tapentadol_emax * TAPENTADOL
    )
    r_arm <- r_arm_raw * exp(eta_study_emax)

    # 5. Observed pain intensity on the 0-10 scale (paper Eq. 1).
    latent  <- base_arm + r_arm * (1 - exp(-kel * time))
    score   <- 10 * exp(latent) / (1 + exp(latent))

    # 6. Additive residual on the observed 0-10 scale. Per-arm
    # sample-size weighting (SD -> addSd / sqrt(N)) is applied
    # downstream in the simulation event-table code; see the
    # vignette Assumptions and Deviations section.
    score ~ add(addSd)
  })
}
