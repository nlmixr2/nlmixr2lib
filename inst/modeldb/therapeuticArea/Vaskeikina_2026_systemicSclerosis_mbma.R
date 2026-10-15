Vaskeikina_2026_systemicSclerosis_mbma <- function() {
  description <- "MBMA. Longitudinal model-based meta-analysis of targeted therapies in systemic sclerosis (SSc), fit in Monolix 2024R1 to study-arm-level aggregate data from 32 randomized controlled trials (73 arms, 2036 patients, 23 targeted agents). Two endpoints are modelled jointly as change from baseline: the modified Rodnan skin score (mrsscfb, points; 30 trials) and percent-predicted forced vital capacity (fvcppcfb, % predicted; 20 trials). Each follows an Emax-in-time function CFB(t) = Emax * t / (ET50 + t) with ONE shared ET50 of 27.49 weeks. The arm Emax is the comparator (standard-of-care) typical value plus an ADDITIVE per-drug shift for the six drugs with a distinguishable mRSS effect (abatacept, baricitinib, guselkumab, inebilizumab, rituximab, tofacitinib) and the four with a distinguishable FVC effect (belimumab, nintedanib, rituximab, tocilizumab), selected by the TRT_<INN> arm indicators. Variability is two-level and additive on Emax: between-study (eta_study_*) and between-treatment-arm (eta_arm_*), the latter and the residual both divided by sqrt(N_ARM), the arm sample size, which is therefore a required input column. No baseline covariate was retained. Suitable for simulating study-arm-mean change-from-baseline trajectories, NOT individual-patient scores."

  reference <- paste(
    "Vaskeikina M, Ugolkov Y, Kireev B, Peskov K, Volkova A.",
    "Systematic Review and Model-Based Meta-Analysis of Targeted Drugs for",
    "Systemic Sclerosis. Pharmaceutics. 2026;18(2):250.",
    "doi:10.3390/pharmaceutics18020250.",
    sep = " "
  )
  vignette <- "Vaskeikina_2026_systemicSclerosis"

  units <- list(
    time = "week",
    dosing = "none",
    concentration = "change from baseline in an arm-mean clinical score: mrsscfb in mRSS points (0-51 scale) and fvcppcfb in percentage points of predicted FVC; NOT a drug concentration"
  )

  trtNotes <- paste(
    "MBMA study-arm-level treatment indicator (a property of the published trial arm, not of an individual patient); 0 on a comparator arm.",
    "Vaskeikina 2026 Equation 11 adds a per-drug effect beta_TRT to the comparator Emax, with the control arm as the reference category.",
    "There is no dose or exposure term: the MBMA pools every dose arm of a drug into one effect (Discussion), so this indicator, not a concentration, is what switches the effect on.",
    "Only drugs whose effect was estimated with RSE <= 50% were kept in the final model (Section 2.3.2); arms of the other agents were removed from the fit while their control arms were retained."
  )

  covariateData <- list(
    TRT_ABATACEPT = list(
      description = "1 = study arm received abatacept, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the mRSS Emax only (Table 1)."),
      source_name = "beta_mRSS_Abatacept (Table 1)"
    ),
    TRT_BARICITINIB = list(
      description = "1 = study arm received baricitinib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the mRSS Emax only (Table 1)."),
      source_name = "beta_mRSS_Baricitinib (Table 1)"
    ),
    TRT_GUSELKUMAB = list(
      description = "1 = study arm received guselkumab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the mRSS Emax only (Table 1)."),
      source_name = "beta_mRSS_Guselkumab (Table 1)"
    ),
    TRT_INEBILIZUMAB = list(
      description = "1 = study arm received inebilizumab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the mRSS Emax only (Table 1)."),
      source_name = "beta_mRSS_Inebilizumab (Table 1)"
    ),
    TRT_RITUXIMAB = list(
      description = "1 = study arm received rituximab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(
        trtNotes,
        "The only drug with an effect on BOTH endpoints (Table 1), so the one indicator enters both Emax expressions."
      ),
      source_name = "beta_mRSS_Rituximab and beta_FVC_Rituximab (Table 1)"
    ),
    TRT_TOFACITINIB = list(
      description = "1 = study arm received tofacitinib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the mRSS Emax only (Table 1)."),
      source_name = "beta_mRSS_Tofacitinib (Table 1)"
    ),
    TRT_BELIMUMAB = list(
      description = "1 = study arm received belimumab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the FVC Emax only (Table 1)."),
      source_name = "beta_FVC_Belimumab (Table 1)"
    ),
    TRT_NINTEDANIB = list(
      description = "1 = study arm received nintedanib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the FVC Emax only (Table 1)."),
      source_name = "beta_FVC_Nintedanib (Table 1)"
    ),
    TRT_TOCILIZUMAB = list(
      description = "1 = study arm received tocilizumab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (comparator / standard of care)",
      notes = paste(trtNotes, "Acts on the FVC Emax only (Table 1)."),
      source_name = "beta_FVC_Tocilizumab (Table 1)"
    ),
    N_ARM = list(
      description = "Number of participants contributing to the study-arm-level mean at that observation; the meta-analytic weight.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Used twice, which is why it is a data column rather than a scaling applied after the solve. (1) Equation 10: the between-treatment-arm random effect on Emax has variance gamma^2 / N, so its SD is gamma / sqrt(N_ARM). (2) Equation 9: the residual has variance sigma^2 / N, so its SD is sigma / sqrt(N_ARM). The between-study random effect is NOT scaled (Section 2.3.1). Treatment-arm sizes in the fitted database have median 15, mean 45-46 and maximum 288 participants (Sections 2.4 and 3.3). Must be strictly positive.",
      source_name = "N_ijk^a (Equations 9, 10, 12 and 13)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Study-arm mean patient age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on Emax as a median-normalised power function (Section 2.3.2) and not retained: 'None of the tested covariates improved model performance or reduced BSV' (Section 3.2). No estimate is published."
    ),
    SEXF_PCT = list(
      description = "Study-arm percentage of enrolled participants who are female.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper's 'proportion of females'. Screened on Emax and not retained (Section 3.2)."
    ),
    DIS_DURATION = list(
      description = "Study-arm mean SSc disease duration.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on Emax and not retained (Section 3.2). Documentation only, so not registered in inst/references/covariate-columns.md; the paper does not state the unit."
    ),
    MRSS_BL = list(
      description = "Study-arm mean baseline modified Rodnan skin score.",
      units = "points",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on Emax and not retained (Section 3.2); the across-study mean baseline was 22.2 points (Section 3.1). Documentation only, so not registered."
    ),
    DIS_DCSSC = list(
      description = "Type of SSc enrolled in the arm (diffuse cutaneous SSc versus other or unspecified).",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = "Screened on Emax as an indicator variable and not retained (Section 3.2). Section 3.1: 16 studies enrolled dcSSc and 11 did not specify the subtype. The exact coding is not published. Documentation only, so not registered."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2036L,
    n_studies = 32L,
    n_arms = 73L,
    age_range = "study-arm means; overall mean age 50.5 years (Section 3.1)",
    weight_range = "not reported",
    sex_female_pct = 75.5,
    race_ethnicity = "collected per arm but not reported in the publication",
    disease_state = "adults with systemic sclerosis enrolled in randomized controlled trials of targeted therapies: diffuse cutaneous SSc in 16 studies, unspecified SSc in 11, SSc-associated pulmonary arterial hypertension in 1 (excluded from modelling because it reported neither mRSS nor FVC); only 4 studies explicitly reported interstitial lung disease",
    dose_range = "not modelled; every dose arm of a drug is pooled into a single per-drug effect",
    regimens = "23 targeted agents (cytokine, surface-receptor, intracellular-kinase and antifibrotic targets) versus a standard-of-care comparator arm",
    baseline_scores = "across-study mean baseline mRSS 22.2 points and percent-predicted FVC 83.4% (Section 3.1)",
    timepoints = "follow-up up to about two years; 16 of 30 mRSS trials and 14 of 20 FVC trials contributed a single post-baseline assessment (Section 3.1)",
    regions = "not reported",
    notes = "MBMA at the study-arm level. mRSS was modelled from 30 trials and percent-predicted FVC from 20 (21 reported FVC; STRATUS was excluded because its mL values could not be converted to percent predicted). Arm means and SDs were harmonised from medians, IQRs, SEs and CIs with published conversion formulas, and data were digitised from figures with PlotDigitizer. The model simulates arm-mean change-from-baseline trajectories and is NOT suitable for individual-subject simulation."
  )

  ini({
    # ================================================================
    # Vaskeikina 2026 Equations 8, 10 and 11 (display equations; the
    # values below are Table 1, 'Value' column, with RSE in comments):
    #
    #   CFB_ik(t) = Emax_ik * t / (ET50 + t)
    #   Emax_ik   = Emaxpop + sum(beta_TRT) + eta_BSV_i + eta_BTAV_ik
    #
    # Emax and the drug effects are signed and additive on the LINEAR
    # scale, so none of them is log-transformed. mRSS falls as skin
    # fibrosis improves (negative = benefit); percent-predicted FVC
    # rises as lung function improves (positive = benefit).
    # ================================================================

    emax_mrss <- -5.29
    label("Maximum change from baseline in mRSS on the comparator (standard-of-care) arm, Emaxpop_mRSS (points)") # Table 1 'Maximum change in mRSS of the comparator' = -5.29 (RSE 22.2%)

    e_trt_abatacept_emax_mrss <- -5.41
    label("Additive effect of abatacept on the mRSS Emax, beta_mRSS_Abatacept (points)") # Table 1 'Effect of abatacept on mRSS reduction' = -5.41 (RSE 44.6%)

    e_trt_baricitinib_emax_mrss <- -13.87
    label("Additive effect of baricitinib on the mRSS Emax, beta_mRSS_Baricitinib (points)") # Table 1 'Effect of baricitinib on mRSS reduction' = -13.87 (RSE 35.4%)

    e_trt_guselkumab_emax_mrss <- -23.55
    label("Additive effect of guselkumab on the mRSS Emax, beta_mRSS_Guselkumab (points)") # Table 1 'Effect of guselkumab on mRSS reduction' = -23.55 (RSE 16.7%)

    e_trt_inebilizumab_emax_mrss <- -13.45
    label("Additive effect of inebilizumab on the mRSS Emax, beta_mRSS_Inebilizumab (points)") # Table 1 'Effect of inebilizumab on mRSS reduction' = -13.45 (RSE 38.0%)

    e_trt_rituximab_emax_mrss <- -8.03
    label("Additive effect of rituximab on the mRSS Emax, beta_mRSS_Rituximab (points)") # Table 1 'Effect of rituximab on mRSS reduction' = -8.03 (RSE 26.3%)

    e_trt_tofacitinib_emax_mrss <- -14.35
    label("Additive effect of tofacitinib on the mRSS Emax, beta_mRSS_Tofacitinib (points)") # Table 1 'Effect of tofacitinib on mRSS reduction' = -14.35 (RSE 18.4%)

    emax_fvc <- -2.29
    label("Maximum change from baseline in percent-predicted FVC on the comparator arm, Emaxpop_FVC (% predicted)") # Table 1 'Maximum change in FVC of the comparator in %' = -2.29 (RSE 35.3%)

    e_trt_nintedanib_emax_fvc <- 2.06
    label("Additive effect of nintedanib on the FVC Emax, beta_FVC_Nintedanib (% predicted)") # Table 1 'Effect of nintedanib on FVC improvement' = 2.06 (RSE 45.5%)

    e_trt_rituximab_emax_fvc <- 9.81
    label("Additive effect of rituximab on the FVC Emax, beta_FVC_Rituximab (% predicted)") # Table 1 'Effect of rituximab on FVC improvement' = 9.81 (RSE 20.1%)

    e_trt_tocilizumab_emax_fvc <- 4.49
    label("Additive effect of tocilizumab on the FVC Emax, beta_FVC_Tocilizumab (% predicted)") # Table 1 'Effect of tocilizumab on FVC improvement' = 4.49 (RSE 31.7%)

    e_trt_belimumab_emax_fvc <- 10.38
    label("Additive effect of belimumab on the FVC Emax, beta_FVC_Belimumab (% predicted)") # Table 1 'Effect of belimumab on FVC improvement' = 10.38 (RSE 43.0%)

    # One ET50 shared by both endpoints; the Discussion reports that an
    # FVC-specific ET50 drifted to an implausible 77.61 weeks, so the
    # authors kept the single shared value. No random effect on ET50.
    let50 <- log(27.49)
    label("Time to 50% of the maximum change from baseline, shared by mRSS and FVC, ET50 (log week)") # Table 1 'Time to 50% maximum change in mRSS and FVC in weeks' = 27.49 (RSE 16.6%)

    # ================================================================
    # Two-level random effects on Emax (Equation 10). Table 1 labels
    # all four as 'Standard deviation of ...', so each ini() value is
    # the SD squared. Both levels are ADDITIVE on the linear Emax.
    #   eta_study_* : between-study variability (BSV, omega), unscaled
    #   eta_arm_*   : between-treatment-arm variability (BTAV, gamma),
    #                 for an arm of ONE participant; model() divides it
    #                 by sqrt(N_ARM) per the gamma^2 / N variance.
    # ================================================================
    eta_study_emax_mrss ~ 20.5209 # Table 1 'Standard deviation of BSV for change in mRSS' omega = 4.53 (RSE 19.6%); variance 4.53^2
    eta_arm_emax_mrss ~ 104.6529 # Table 1 'Standard deviation of BTAV for change in mRSS' gamma = 10.23 (RSE 28.7%); variance 10.23^2, scaled by 1/N_ARM in model()
    eta_study_emax_fvc ~ 3.8809 # Table 1 'Standard deviation of BSV for change in FVC' omega = 1.97 (RSE 38.0%); variance 1.97^2
    eta_arm_emax_fvc ~ 98.8036 # Table 1 'Standard deviation of BTAV for change in FVC' gamma = 9.94 (RSE 30.1%); variance 9.94^2, scaled by 1/N_ARM in model()

    # ================================================================
    # Residual error (Equation 9): constant, with variance sigma^2 / N,
    # so the per-observation SD is addSd / sqrt(N_ARM).
    # ================================================================
    addSd_mrsscfb <- 4.12
    label("Additive residual SD of the mRSS change from baseline for an arm of one participant; the per-observation SD is addSd_mrsscfb / sqrt(N_ARM) (points)") # Table 1 'Standard deviation of residual error for change in mRSS' sigma = 4.12 (RSE 10.0%)

    addSd_fvcppcfb <- 4.42
    label("Additive residual SD of the percent-predicted FVC change from baseline for an arm of one participant; the per-observation SD is addSd_fvcppcfb / sqrt(N_ARM) (% predicted)") # Table 1 'Standard deviation of residual error for change in FVC' sigma = 4.42 (RSE 12.0%)
  })

  model({
    # Algebraic model with no ODE and no dosing: `time` is weeks since
    # the start of treatment. Study-arm inputs supplied per row are the
    # nine TRT_<INN> arm indicators (all 0 on a comparator arm) and
    # N_ARM, the arm sample size.

    # Shared time course (Equation 8).
    et50 <- exp(let50)
    tcourse <- time / (et50 + time)

    # Arm Emax for mRSS (Equations 10 and 11). The between-arm term is
    # divided by sqrt(N_ARM); the between-study term is not.
    btav_mrss <- eta_arm_emax_mrss / sqrt(N_ARM)
    emax_mrss_arm <- emax_mrss +
      e_trt_abatacept_emax_mrss * TRT_ABATACEPT +
      e_trt_baricitinib_emax_mrss * TRT_BARICITINIB +
      e_trt_guselkumab_emax_mrss * TRT_GUSELKUMAB +
      e_trt_inebilizumab_emax_mrss * TRT_INEBILIZUMAB +
      e_trt_rituximab_emax_mrss * TRT_RITUXIMAB +
      e_trt_tofacitinib_emax_mrss * TRT_TOFACITINIB +
      eta_study_emax_mrss + btav_mrss

    # Arm Emax for percent-predicted FVC (Equations 10 and 11).
    btav_fvc <- eta_arm_emax_fvc / sqrt(N_ARM)
    emax_fvc_arm <- emax_fvc +
      e_trt_belimumab_emax_fvc * TRT_BELIMUMAB +
      e_trt_nintedanib_emax_fvc * TRT_NINTEDANIB +
      e_trt_rituximab_emax_fvc * TRT_RITUXIMAB +
      e_trt_tocilizumab_emax_fvc * TRT_TOCILIZUMAB +
      eta_study_emax_fvc + btav_fvc

    # Arm-mean change from baseline in each endpoint.
    mrsscfb <- emax_mrss_arm * tcourse
    fvcppcfb <- emax_fvc_arm * tcourse

    # Arm-size-weighted residual (Equation 9).
    sd_mrsscfb <- addSd_mrsscfb / sqrt(N_ARM)
    sd_fvcppcfb <- addSd_fvcppcfb / sqrt(N_ARM)
    mrsscfb ~ add(sd_mrsscfb)
    fvcppcfb ~ add(sd_fvcppcfb)
  })
}
