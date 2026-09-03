Rayner_2013_oseltamivir_shedding <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for the time to cessation",
    "of influenza virus shedding in healthy adults experimentally inoculated with",
    "influenza A/Texas/36/91 (H1N1) or B/Yamagata/16/88 and treated with oral",
    "oseltamivir or placebo (Rayner 2013, pooled phase 2 inoculation studies",
    "PV15616 and NP15717).",
    "The exposure driver is the steady-state (day-5) AUC from 0 to 24 h of the",
    "active metabolite oseltamivir carboxylate (OC), taken as a per-subject input",
    "covariate from post hoc estimates of the companion population PK model",
    "(Kamal 2013; available as modellib('Kamal_2013_oseltamivir')).",
    "OC AUC0-24 enters as a 3-group categorical variable whose reference level is",
    "zero exposure (placebo) and whose upper cutoff is 14,180 ng*h/mL; the",
    "neuraminidase-inhibition IC50 of the infecting strain (0.18 vs 16.76 nM, a",
    "surrogate for study and virus type) is retained as a second covariate",
    "although it is not statistically significant.",
    "A higher hazard means EARLIER cessation of shedding, so the hazard ratios",
    "above 1 for the drug-treated (1.72) and high-exposure (2.42) groups mean",
    "faster viral clearance with increasing OC exposure.",
    "Because the fit is a semiparametric Cox regression, the baseline hazard is",
    "unspecified by construction and is not reported; the model therefore returns",
    "the RELATIVE hazard hr against the placebo reference cohort, not an absolute",
    "survivor function. Multiply hr by any user-supplied baseline hazard h0(t) to",
    "obtain a subject hazard.",
    "The model is algebraic and deterministic (no ODE state, no drug input, no",
    "IIV, no residual error).",
    "Companion models from the same paper are",
    "modellib('Rayner_2013_oseltamivir_symptomscore') and",
    "modellib('Rayner_2013_oseltamivir_symptomalleviation')."
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
    concentration = "hr (hazard of cessation of viral shedding relative to the reference cohort, unitless); not a drug concentration"
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
      reference_category = "0 ng*h/mL (placebo)",
      notes              = paste(
        "Rayner 2013 Methods, 'Determination of plasma drug exposures': individual",
        "post hoc PK parameter estimates from the companion Kamal 2013 population PK",
        "model were used to generate individual predicted steady-state concentrations",
        "every 0.1 h during day 5 of therapy; AUC0-24 was then computed by the linear",
        "trapezoidal rule. Placebo subjects were assigned AUC0-24 = 0.",
        "For this endpoint the lower cutoff the authors selected is exactly zero, so",
        "the reference group is the 22 placebo recipients and the middle group is",
        "every drug-treated subject up to 14,180 ng*h/mL. That makes the",
        "reference-to-middle contrast a placebo-versus-drug comparison rather than a",
        "low-versus-moderate exposure comparison, which the paper notes explicitly",
        "(Results, 'Univariable analyses'). The upper cutoff of 14,180 ng*h/mL is",
        "close to the 14,497 and 13,638 ng*h/mL upper cutoffs found independently for",
        "the two symptom endpoints, which is the paper's headline finding of a common",
        "roughly 14,000 ng*h/mL threshold. For orientation, the average AUC0-24 at",
        "the approved 75 mg BID regimen is about 6,000 ng*h/mL (Discussion).",
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
        "The hazard ratio of 0.66 for the influenza B stratum is not significant",
        "(P = 0.136) and its confidence interval spans 1; it is retained because",
        "IC50 is collinear with OC AUC0-24 across the two studies. This is the",
        "largest of the three endpoints' IC50 effects, and the univariable",
        "relationship for this endpoint was nominally significant (P = 0.026,",
        "Table 2), so the multivariable analysis is what establishes that strain",
        "identity does not drive the exposure-response relationship."
      ),
      source_name        = "IC50"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 92L,
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
      "n_subjects = 92 here is the number contributing to THIS endpoint, from the",
      "Table 4 subgroup counts for time to cessation of viral shedding",
      "(22 + 49 + 21); all 22 subjects in the reference group were placebo",
      "recipients, by construction of the zero-exposure cutoff. The endpoint is only",
      "defined for subjects with a positive viral culture, which is why this subset",
      "differs from the 112 subjects evaluable for composite symptom score AUC and",
      "the 64 for symptom alleviation.",
      "The event is defined in Methods, 'Efficacy endpoints': time from the first",
      "positive nasal-lavage viral culture until the first negative culture. Nasal",
      "washes were taken twice daily early in the study and once daily thereafter,",
      "so observed event times are coarsely quantised.",
      "This is a retrospective pharmacometric analysis of pooled phase 2",
      "experimental-inoculation data; the authors caution that it may not represent",
      "efficacy in natural influenza infection (Discussion)."
    )
  )

  ini({
    # ==================================================================
    # Rayner 2013 Table 6, final multivariable Cox proportional-hazards
    # model for the time to cessation of viral shedding, with OC AUC0-24
    # as a 3-group variable and IC50 as a two-level categorical
    # variable:
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
    # (reference group 2.5, 4.75 and 7 days) as an empirical description
    # of the absolute time course. See the vignette's "Assumptions and
    # deviations" section for why a calibrated parametric baseline was
    # rejected.
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
    # jointly with the hazard ratios, hence fixed(). The lower cutoff is
    # exactly zero for this endpoint, so the reference group is placebo.
    auc_cut1 <- fixed(0);     label("Lower OC AUC0-24 cutoff separating the placebo reference and drug-treated groups (ng*h/mL)")  # Rayner 2013 Table 6, viral-shedding rows: reference group 0
    auc_cut2 <- fixed(14180); label("Upper OC AUC0-24 cutoff separating the middle and high exposure groups (ng*h/mL)")            # Rayner 2013 Table 6, viral-shedding rows: > 14,180

    # ----- Reference IC50 -----
    ic50_ref <- fixed(0.18); label("Reference neuraminidase-inhibition IC50, influenza A/Texas/36/91 (H1N1) (nM)")  # Rayner 2013 Table 6, IC50 (nM) reference group 0.18

    # ----- Estimated log hazard ratios (Table 6) -----
    # Both AUC0-24 coefficients are positive: higher OC exposure raises
    # the hazard of cessation of shedding, i.e. shortens the shedding
    # duration. The high-vs-middle hazard ratio implied by these two
    # estimates is 2.42 / 1.72 = 1.407, matching the 1.40 printed for
    # that pairwise comparison in the same table.
    e_auc_oselcarb_mid_haz  <- log(1.72); label("Log hazard ratio for cessation of viral shedding, drug-treated OC AUC0-24 group vs placebo (log scale; HR 1.72)")  # Rayner 2013 Table 6: HR 1.72 (95% CI 0.98, 3.03), P = 0.059
    e_auc_oselcarb_high_haz <- log(2.42); label("Log hazard ratio for cessation of viral shedding, high OC AUC0-24 group vs placebo (log scale; HR 2.42)")           # Rayner 2013 Table 6: HR 2.42 (95% CI 1.20, 4.84), P = 0.013

    # Negative: subjects inoculated with influenza B/Yamagata cleared
    # virus later, but the effect is not significant (P = 0.136) and its
    # CI spans 1.
    e_ic50_neuraminidase_haz <- log(0.66); label("Log hazard ratio for cessation of viral shedding, influenza B/Yamagata (IC50 16.76 nM) vs influenza A/Texas (IC50 0.18 nM) (log scale; HR 0.66)")  # Rayner 2013 Table 6: HR 0.66 (95% CI 0.38, 1.14), P = 0.136

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
    #    reference group, which for this endpoint is OC AUC0-24 = 0,
    #    i.e. placebo. Boundary convention follows Table 6 exactly: the
    #    middle group is "> auc_cut1 to <= auc_cut2".
    # ------------------------------------------------------------------
    aucMid  <- (AUC_OSELCARB > auc_cut1) * (AUC_OSELCARB <= auc_cut2)
    aucHigh <- (AUC_OSELCARB > auc_cut2)

    # ------------------------------------------------------------------
    # 2. Recover the two-level IC50 categorical variable.
    # ------------------------------------------------------------------
    highIc50 <- (IC50_NEURAMINIDASE > ic50_ref)

    # ------------------------------------------------------------------
    # 3. Linear predictor and relative hazard against the reference
    #    cohort (placebo, influenza A/Texas). Under the proportional-
    #    hazards assumption hr is constant in time, so the subject
    #    hazard is h(t) = h0(t) * hr for any baseline h0(t) the user
    #    supplies. hr is the quantity to gate against Table 6.
    # ------------------------------------------------------------------
    lhr <- e_auc_oselcarb_mid_haz * aucMid +
      e_auc_oselcarb_high_haz * aucHigh +
      e_ic50_neuraminidase_haz * highIc50

    hr <- exp(lhr)
  })
}
