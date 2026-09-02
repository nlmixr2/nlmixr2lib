Rayner_2013_oseltamivir_symptomalleviation <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for the time to alleviation",
    "of the composite influenza symptom score in healthy adults experimentally",
    "inoculated with influenza A/Texas/36/91 (H1N1) or B/Yamagata/16/88 and",
    "treated with oral oseltamivir or placebo (Rayner 2013, pooled phase 2",
    "inoculation studies PV15616 and NP15717).",
    "The exposure driver is the steady-state (day-5) AUC from 0 to 24 h of the",
    "active metabolite oseltamivir carboxylate (OC), taken as a per-subject input",
    "covariate from post hoc estimates of the companion population PK model",
    "(Kamal 2013; available as modellib('Kamal_2013_oseltamivir')).",
    "OC AUC0-24 enters as a 3-group categorical variable with cutoffs 1,568 and",
    "13,638 ng*h/mL; the neuraminidase-inhibition IC50 of the infecting strain",
    "(0.18 vs 16.76 nM, a surrogate for study and virus type) is retained as a",
    "second covariate although it is not statistically significant.",
    "A higher hazard means EARLIER symptom alleviation, so the hazard ratios",
    "above 1 for the middle (1.85) and high (5.34) exposure groups mean faster",
    "recovery with increasing OC exposure.",
    "Because the fit is a semiparametric Cox regression, the baseline hazard is",
    "unspecified by construction and is not reported; the model therefore returns",
    "the RELATIVE hazard hr against the lowest-exposure influenza A reference",
    "cohort, not an absolute survivor function. Multiply hr by any user-supplied",
    "baseline hazard h0(t) to obtain a subject hazard.",
    "The model is algebraic and deterministic (no ODE state, no drug input, no",
    "IIV, no residual error).",
    "Companion models from the same paper are",
    "modellib('Rayner_2013_oseltamivir_symptomscore') and",
    "modellib('Rayner_2013_oseltamivir_shedding')."
  )
  reference <- paste(
    "Rayner CR, Bulik CC, Kamal MA, Reynolds DK, Toovey S, Hammel JP, Smith PF,",
    "Bhavnani SM, Van Wart SA, Ambrose PG, Forrest A.",
    "Pharmacokinetic-pharmacodynamic determinants of oseltamivir efficacy using",
    "data from phase 2 inoculation studies.",
    "Antimicrob Agents Chemother. 2013;57(8):3478-3487. doi:10.1128/AAC.02440-12.",
    "Hazard ratios are from Table 6 (final multivariable model); the univariable",
    "cross-check hazard ratios are from Table 3 and the subgroup Kaplan-Meier",
    "quartiles are from Table 4.",
    "The upstream population PK model that generates the OC AUC0-24 covariate is",
    "Kamal MA et al. Antimicrob Agents Chemother. 2013;57(8):3470-3477,",
    "doi:10.1128/AAC.02438-12; see modellib('Kamal_2013_oseltamivir').",
    sep = " "
  )
  vignette <- "Rayner_2013_oseltamivir_exposure_response"
  units <- list(
    time = "day",
    dosing = "n/a (exposure-response model; OC exposure enters as the AUC_OSELCARB covariate, not as a dose record)",
    concentration = "hr (hazard of symptom alleviation relative to the reference cohort, unitless); not a drug concentration"
  )

  covariateData <- list(
    AUC_OSELCARB = list(
      description        = paste(
        "Steady-state area under the plasma concentration-time curve from 0 to 24 h",
        "of oseltamivir carboxylate (OC), the active metabolite of the prodrug",
        "oseltamivir. Per-subject, time-fixed."
      ),
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Rayner 2013 Methods, 'Determination of plasma drug exposures': individual",
        "post hoc PK parameter estimates from the companion Kamal 2013 population PK",
        "model were used to generate individual predicted steady-state concentrations",
        "every 0.1 h during day 5 of therapy; AUC0-24 was then computed by the linear",
        "trapezoidal rule. Placebo subjects were assigned AUC0-24 = 0.",
        "The model consumes this column ONLY through the two printed cutoffs",
        "auc_cut1 = 1,568 and auc_cut2 = 13,638 ng*h/mL, so it is effectively a",
        "3-level categorical predictor. Note that these cutoffs differ from the",
        "1,495 / 14,497 ng*h/mL pair used for the composite symptom score AUC",
        "endpoint: each pair was optimised separately for its own endpoint",
        "(Rayner 2013 Results). For orientation, the average AUC0-24 at the approved",
        "75 mg BID regimen is about 6,000 ng*h/mL (Discussion).",
        "Generate this column with modellib('Kamal_2013_oseltamivir')."
      ),
      source_name        = "AUC0-24"
    ),
    IC50_NEURAMINIDASE = list(
      description        = paste(
        "Concentration of the neuraminidase inhibitor that reduces neuraminidase",
        "activity of the infecting influenza virus by 50%, measured in vitro.",
        "Per-subject, time-fixed (a property of the inoculating strain)."
      ),
      units              = "nM",
      type               = "continuous",
      reference_category = "0.18 nM (influenza A/Texas/36/91 (H1N1), study 1 / PV15616)",
      notes              = paste(
        "Takes exactly two values in this dataset: 0.18 +/- 0.11 nM for influenza",
        "A/Texas/36/91 (H1N1) inoculated in study 1, and 16.76 +/- 4.10 nM for",
        "influenza B/Yamagata/16/88 inoculated in study 2 (Rayner 2013 Methods,",
        "'NA inhibition assay'). Perfectly collinear with study identity; the paper",
        "describes it as 'a surrogate for study name or virus type' (Discussion).",
        "Fitted as a two-level CATEGORICAL variable (Table 2 footnote c), so the",
        "model applies the full effect to any value above the printed reference of",
        "0.18 nM rather than interpolating a slope. Only the two studied levels are",
        "supported.",
        "The hazard ratio of 0.82 for the influenza B stratum is not significant",
        "(P = 0.494) and its confidence interval spans 1; it is retained because",
        "IC50 is collinear with OC AUC0-24 across the two studies."
      ),
      source_name        = "IC50"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 2L,
    age_range      = "18 years and older by inclusion criterion; mean 23.4 years (CV 25.4%) pooled, 22.3 years (CV 20.3%) in study 1 and 25.0 years (CV 29.3%) in study 2 (Rayner 2013 Table 1)",
    weight_range   = "mean 70.6 kg (CV 18.0%) pooled, 68.6 kg (CV 20.5%) in study 1 and 73.7 kg (CV 13.4%) in study 2 (Rayner 2013 Table 1)",
    height_mean    = "172 cm (CV 5.45%) pooled (Rayner 2013 Table 1)",
    sex_female_pct = 51.0,
    race_ethnicity = c(White = 84.4, Black = 6.96, Other = 8.70),
    renal_function = "Creatinine clearance mean 114 mL/min/1.73 m^2 (CV 22.2%) pooled (Rayner 2013 Table 1). Healthy volunteers with normal renal function.",
    disease_state  = paste(
      "Healthy adult volunteers experimentally inoculated with influenza virus by",
      "nasal drops (TCID50 10^6) and housed in an isolation unit. Study 1 (PV15616)",
      "used influenza A/Texas/36/91 (H1N1) with treatment started 28 h after",
      "inoculation; study 2 (NP15717) used influenza B/Yamagata/16/88. Entry",
      "required a pre-existing influenza antibody titre of 1:8 or less (study 1) or",
      "below 1:10 (study 2)."
    ),
    dose_range     = paste(
      "Study 1: placebo, oseltamivir 20 mg BID, 100 mg BID, 200 mg QD or 200 mg BID",
      "orally for 5 days. Study 2: placebo, oseltamivir 75 mg BID or 150 mg BID",
      "orally for 5 days (Rayner 2013 Table 1)."
    ),
    notes          = paste(
      "115 subjects were evaluable overall (69 in study 1, 46 in study 2).",
      "n_subjects = 64 here is the number contributing to THIS endpoint, from the",
      "Table 4 subgroup counts for time to alleviation of composite symptom score",
      "(20 + 30 + 14); 20 of those 64 fell in the lowest OC AUC0-24 group and 18 of",
      "the 20 were placebo recipients. The endpoint is only defined for subjects who",
      "developed symptoms, which is why this subset is smaller than the 112 subjects",
      "evaluable for composite symptom score AUC and the 92 for cessation of viral",
      "shedding.",
      "The event is defined in Methods, 'Efficacy endpoints': time from the point at",
      "which any of the seven individual symptom scores was rated above 1 until the",
      "time at which all applicable individual scores were 1 or lower. Symptoms were",
      "self-rated twice daily, so observed event times are quantised to roughly",
      "half-day steps.",
      "This is a retrospective pharmacometric analysis of pooled phase 2",
      "experimental-inoculation data; the authors caution that it may not represent",
      "efficacy in natural influenza infection (Discussion)."
    )
  )

  ini({
    # ==================================================================
    # Rayner 2013 Table 6, final multivariable Cox proportional-hazards
    # model for the time to alleviation of the composite symptom score,
    # with OC AUC0-24 as a 3-group variable and IC50 as a two-level
    # categorical variable:
    #
    #   h(t) = h0(t) * exp(
    #       e_auc_oselcarb_mid_haz  * I(auc_cut1 < AUC <= auc_cut2)
    #     + e_auc_oselcarb_high_haz * I(AUC > auc_cut2)
    #     + e_ic50_neuraminidase_haz * I(IC50 > ic50_ref) )
    #
    # The coefficients below are the natural logs of the printed hazard
    # ratios, which is the scale the Cox model estimates on.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is left completely unspecified by the method, so it is not
    # an unreported parameter but a quantity the fit never produced.
    # This model therefore returns the RELATIVE hazard hr only, and
    # deliberately does not offer a survivor function. Table 4 gives the
    # paper's per-subgroup Kaplan-Meier 25th / 50th / 75th percentiles
    # (reference group 2, 3.5 and 4.5 days) as an empirical description
    # of the absolute time course; a Weibull baseline calibrated to those
    # quartiles was evaluated during extraction and rejected, because
    # combining it with the published hazard ratios overpredicts the
    # middle group's own published 25th and 50th percentiles by 46% and
    # 53%. Shipping such a baseline would embed a construction that
    # contradicts the paper's own data. See the vignette's "Assumptions
    # and deviations" section.
    #
    # The regression was fitted in R 2.11.1, not NONMEM; the paper
    # reports hazard ratios with 95% confidence intervals and no
    # variance components, so there is no IIV and no residual error to
    # encode and no observation endpoint is declared. This matches the
    # deterministic algebraic pattern of
    # Nagy_2017_obiltoxaximab_survival.R and
    # Lin_2020_glasdegib_decitabine.R.
    # ==================================================================

    # ----- OC AUC0-24 3-group cutoffs -----
    # Optimally determined by the authors as the pair of cutoffs
    # minimising the log rank P value from Cox proportional-hazards
    # regression (Methods, 'Univariable analyses'), not estimated
    # jointly with the hazard ratios, hence fixed().
    auc_cut1 <- fixed(1568);  label("Lower OC AUC0-24 cutoff separating the low and middle exposure groups (ng*h/mL)")   # Rayner 2013 Table 6, time-to-alleviation rows: reference <= 1,568
    auc_cut2 <- fixed(13638); label("Upper OC AUC0-24 cutoff separating the middle and high exposure groups (ng*h/mL)")  # Rayner 2013 Table 6, time-to-alleviation rows: > 13,638

    # ----- Reference IC50 -----
    ic50_ref <- fixed(0.18); label("Reference neuraminidase-inhibition IC50, influenza A/Texas/36/91 (H1N1) (nM)")  # Rayner 2013 Table 6, IC50 (nM) reference group 0.18

    # ----- Estimated log hazard ratios (Table 6) -----
    # Both AUC0-24 coefficients are positive: higher OC exposure raises
    # the hazard of alleviation, i.e. shortens the time to recovery.
    # The high-vs-middle hazard ratio implied by these two estimates is
    # 5.34 / 1.85 = 2.887, matching the 2.89 printed for that pairwise
    # comparison in the same table.
    e_auc_oselcarb_mid_haz  <- log(1.85); label("Log hazard ratio for symptom alleviation, middle OC AUC0-24 group vs low group (log scale; HR 1.85)")  # Rayner 2013 Table 6: HR 1.85 (95% CI 1.01, 3.39), P = 0.045
    e_auc_oselcarb_high_haz <- log(5.34); label("Log hazard ratio for symptom alleviation, high OC AUC0-24 group vs low group (log scale; HR 5.34)")    # Rayner 2013 Table 6: HR 5.34 (95% CI 2.37, 12.07), P < 0.001

    # Negative: subjects inoculated with influenza B/Yamagata recovered
    # marginally later, but the effect is not significant (P = 0.494)
    # and its CI spans 1.
    e_ic50_neuraminidase_haz <- log(0.82); label("Log hazard ratio for symptom alleviation, influenza B/Yamagata (IC50 16.76 nM) vs influenza A/Texas (IC50 0.18 nM) (log scale; HR 0.82)")  # Rayner 2013 Table 6: HR 0.82 (95% CI 0.46, 1.45), P = 0.494

    # No IIV: a Cox regression on one event time per subject has no
    # subject-level random effect.

    # No residual error: a partial-likelihood Cox fit has no residual
    # variance component, and none is reported.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Recover the paper's 3-group OC AUC0-24 categorical variable
    #    from the continuous exposure column. The two indicators are
    #    mutually exclusive; when both are 0 the subject is in the
    #    reference (lowest-exposure) group, OC AUC0-24 <= 1,568 ng*h/mL,
    #    which is mostly placebo. Boundary convention follows Table 6
    #    exactly: the middle group is "> auc_cut1 to <= auc_cut2".
    # ------------------------------------------------------------------
    aucMid  <- (AUC_OSELCARB > auc_cut1) * (AUC_OSELCARB <= auc_cut2)
    aucHigh <- (AUC_OSELCARB > auc_cut2)

    # ------------------------------------------------------------------
    # 2. Recover the two-level IC50 categorical variable.
    # ------------------------------------------------------------------
    highIc50 <- (IC50_NEURAMINIDASE > ic50_ref)

    # ------------------------------------------------------------------
    # 3. Linear predictor and relative hazard against the reference
    #    cohort (lowest OC AUC0-24 group, influenza A/Texas). Under the
    #    proportional-hazards assumption hr is constant in time, so the
    #    subject hazard is h(t) = h0(t) * hr for any baseline h0(t) the
    #    user supplies. hr is the quantity to gate against Table 6.
    # ------------------------------------------------------------------
    lhr <- e_auc_oselcarb_mid_haz * aucMid +
      e_auc_oselcarb_high_haz * aucHigh +
      e_ic50_neuraminidase_haz * highIc50

    hr <- exp(lhr)
  })
}
