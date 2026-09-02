Rayner_2013_oseltamivir_symptomscore <- function() {
  description <- paste(
    "Exposure-response linear-regression model for the area under the composite",
    "influenza symptom score curve (AUCSC) over 9 days in healthy adults",
    "experimentally inoculated with influenza A/Texas/36/91 (H1N1) or",
    "B/Yamagata/16/88 and treated with oral oseltamivir or placebo",
    "(Rayner 2013, pooled phase 2 inoculation studies PV15616 and NP15717).",
    "The exposure driver is the steady-state (day-5) AUC from 0 to 24 h of the",
    "active metabolite oseltamivir carboxylate (OC), taken as a per-subject",
    "input covariate from post hoc estimates of the companion population PK",
    "model (Kamal 2013; available as modellib('Kamal_2013_oseltamivir')).",
    "OC AUC0-24 enters as a 3-group categorical variable with cutoffs 1,495 and",
    "14,497 ng*h/mL, which the authors determined by minimising the likelihood-",
    "ratio P value; the neuraminidase-inhibition IC50 of the infecting strain",
    "(0.18 vs 16.76 nM, a surrogate for study and virus type) is retained as a",
    "second covariate although it is not statistically significant.",
    "The model is algebraic and deterministic (no ODE state, no drug input, no",
    "IIV, no residual error): it returns the predicted AUCSC and the predicted",
    "contrast against the lowest-exposure reference cohort.",
    "Companion time-to-event models from the same paper are",
    "modellib('Rayner_2013_oseltamivir_symptomalleviation') and",
    "modellib('Rayner_2013_oseltamivir_shedding')."
  )
  reference <- paste(
    "Rayner CR, Bulik CC, Kamal MA, Reynolds DK, Toovey S, Hammel JP, Smith PF,",
    "Bhavnani SM, Van Wart SA, Ambrose PG, Forrest A.",
    "Pharmacokinetic-pharmacodynamic determinants of oseltamivir efficacy using",
    "data from phase 2 inoculation studies.",
    "Antimicrob Agents Chemother. 2013;57(8):3478-3487. doi:10.1128/AAC.02440-12.",
    "Parameter estimates are from Table 6 (final multivariable model);",
    "the reference-cohort level is from Table 4 and the univariable",
    "cross-check estimates are from Table 3.",
    "The upstream population PK model that generates the OC AUC0-24 covariate is",
    "Kamal MA et al. Antimicrob Agents Chemother. 2013;57(8):3470-3477,",
    "doi:10.1128/AAC.02438-12; see modellib('Kamal_2013_oseltamivir').",
    sep = " "
  )
  vignette <- "Rayner_2013_oseltamivir_exposure_response"
  units <- list(
    time = "day",
    dosing = "n/a (exposure-response model; OC exposure enters as the AUC_OSELCARB covariate, not as a dose record)",
    concentration = "aucsc (area under the composite symptom score curve over 9 days, score*day); not a drug concentration"
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
        "model (13 pooled clinical studies, including the two studies analysed here)",
        "were used to generate individual predicted steady-state concentrations every",
        "0.1 h during day 5 of therapy; AUC0-24 was then computed by the linear",
        "trapezoidal rule. Placebo subjects were assigned AUC0-24 = 0 (Methods,",
        "'Pharmacokinetic-pharmacodynamic analyses').",
        "The model consumes this column ONLY through the two printed cutoffs",
        "auc_cut1 = 1,495 and auc_cut2 = 14,497 ng*h/mL, so it is effectively a",
        "3-level categorical predictor; the paper reports no continuous-slope",
        "parameter for the final model. Cmax and Cmin were also derived and were",
        "correlated with AUC0-24 at Spearman rho > 0.96, so AUC0-24 was the only",
        "exposure measure carried into the analyses (Results, 'Summary of exposure",
        "measures'). For orientation, the average AUC0-24 at the approved 75 mg BID",
        "regimen is about 6,000 ng*h/mL (Rayner 2013 Discussion).",
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
        "influenza B/Yamagata/16/88 inoculated in study 2 (mean +/- SD of triplicate",
        "fluorometric NA-inhibition assays; Rayner 2013 Methods, 'NA inhibition",
        "assay'). Because each study inoculated one strain, this column is perfectly",
        "collinear with study identity, and the paper describes it as 'a surrogate",
        "for study name or virus type' (Discussion).",
        "The paper fitted it as a two-level CATEGORICAL variable (Table 2 footnote c:",
        "'IC50 was evaluated as a categorical variable, since it only had values of",
        "0.18 and 16.76 nM'), so the model applies the full effect to any value above",
        "the printed reference of 0.18 nM rather than interpolating a slope. Only the",
        "two studied levels are supported; do not supply an intermediate value and",
        "expect a graded effect.",
        "The effect is retained in the final model despite being non-significant",
        "(P = 0.257) because of its collinearity with OC AUC0-24 (Results,",
        "'Multivariable analyses')."
      ),
      source_name        = "IC50"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112L,
    n_studies      = 2L,
    age_range      = "18 years and older by inclusion criterion; mean 23.4 years (CV 25.4%) pooled, 22.3 years (CV 20.3%) in study 1 and 25.0 years (CV 29.3%) in study 2 (Rayner 2013 Table 1)",
    weight_range   = "mean 70.6 kg (CV 18.0%) pooled, 68.6 kg (CV 20.5%) in study 1 and 73.7 kg (CV 13.4%) in study 2 (Rayner 2013 Table 1)",
    height_mean    = "172 cm (CV 5.45%) pooled (Rayner 2013 Table 1)",
    sex_female_pct = 51.0,
    race_ethnicity = c(White = 84.4, Black = 6.96, Other = 8.70),
    renal_function = "Creatinine clearance mean 114 mL/min/1.73 m^2 (CV 22.2%) pooled; 108 (CV 23.8%) in study 1 and 124 (CV 17.8%) in study 2 (Rayner 2013 Table 1). Healthy volunteers with normal renal function.",
    disease_state  = paste(
      "Healthy adult volunteers experimentally inoculated with influenza virus by",
      "nasal drops (TCID50 10^6) and housed in an isolation unit. Study 1 (PV15616)",
      "used influenza A/Texas/36/91 (H1N1) with treatment started 28 h after",
      "inoculation; study 2 (NP15717) used influenza B/Yamagata/16/88. Entry",
      "required a pre-existing influenza antibody titre of 1:8 or less (study 1) or",
      "below 1:10 (study 2). Both were single-centre, multiple-dose, double-blind,",
      "randomised, placebo-controlled studies."
    ),
    dose_range     = paste(
      "Study 1: placebo, oseltamivir 20 mg BID, 100 mg BID, 200 mg QD or 200 mg BID",
      "orally for 5 days. Study 2: placebo, oseltamivir 75 mg BID or 150 mg BID",
      "orally for 5 days. Total daily doses 40 to 400 mg (Rayner 2013 Table 1)."
    ),
    notes          = paste(
      "115 subjects were evaluable overall (69 in study 1, 46 in study 2), of whom",
      "86 received oseltamivir and had PK data (Rayner 2013 Results, 'Subject",
      "population'). n_subjects = 112 here is the number contributing to THIS",
      "endpoint, from the Table 4 subgroup counts for composite symptom score AUC",
      "(30 + 64 + 18); 30 of those 112 fell in the lowest OC AUC0-24 group and 29 of",
      "the 30 were placebo recipients. The companion time-to-event endpoints used",
      "different evaluable subsets (64 for symptom alleviation, 92 for cessation of",
      "viral shedding), which is why the three model files carry different",
      "n_subjects.",
      "This is a retrospective pharmacometric analysis of pooled phase 2",
      "experimental-inoculation data; the authors caution that it may not represent",
      "efficacy in natural influenza infection (Discussion)."
    )
  )

  ini({
    # ==================================================================
    # Rayner 2013 Table 6, final multivariable linear-regression model
    # for composite symptom score AUC (AUCSC), with OC AUC0-24 as a
    # 3-group variable and IC50 as a two-level categorical variable:
    #
    #   AUCSC = aucsc_int
    #           + e_auc_oselcarb_mid_aucsc  * I(auc_cut1 < AUC <= auc_cut2)
    #           + e_auc_oselcarb_high_aucsc * I(AUC > auc_cut2)
    #           + e_ic50_neuraminidase_aucsc * I(IC50 > ic50_ref)
    #
    # AUCSC is the 9-day area under the composite symptom score curve
    # (linear trapezoidal rule; the composite score is the sum of seven
    # individual 0-3 symptom scores self-rated twice daily), so its units
    # are score*day.
    #
    # The regression was fitted in R 2.11.1, not NONMEM; the paper
    # reports point estimates with 95% confidence intervals and no
    # variance components. There is therefore no IIV and no residual
    # error to encode, and no observation endpoint is declared -- the
    # predicted AUCSC falls out as a derived variable. This matches the
    # deterministic algebraic pattern of
    # Nagy_2017_obiltoxaximab_survival.R and
    # Lin_2020_glasdegib_decitabine.R.
    # ==================================================================

    # ----- Reference level -----
    # NOT PAPER-DERIVED AS A MODEL INTERCEPT. Table 6 reports only the
    # pairwise CONTRASTS of the final multivariable model; the intercept
    # is not printed and cannot be recovered exactly, because that would
    # need the IC50-by-AUC-group cross-tabulation the paper does not
    # give. The value used here, 14.6 score*day, is the Table 4 MEAN
    # AUCSC of the lowest OC AUC0-24 subgroup (n = 30, 29 of them
    # placebo), which is the fitted value of that cell in the
    # UNIVARIABLE 3-group model: Table 4's 14.6 minus Table 3's -5.50
    # and -7.88 reproduces Table 4's other two subgroup means (9.1 and
    # 6.7) exactly, which confirms the identification.
    #
    # Consequence: every CONTRAST the model predicts is exact, while
    # absolute AUCSC predictions carry a bias of at most the size of the
    # IC50 effect (1.77 score*day) because 14.6 is marginal over the two
    # strata whereas the model intercept is conditional on the
    # influenza A stratum. Use d_aucsc, not aucsc, when an exact
    # published quantity is needed. See the vignette's "Assumptions and
    # deviations" section.
    aucsc_int <- fixed(14.6); label("Composite symptom score AUC in the reference cohort: lowest OC AUC0-24 group, influenza A/Texas (score*day over 9 days)")  # Rayner 2013 Table 4, composite symptom score AUC, subgroup <= 1,495 ng*h/mL; NOT the Table 6 intercept, which is unreported

    # ----- OC AUC0-24 3-group cutoffs -----
    # Optimally determined by the authors as the pair of cutoffs
    # minimising the likelihood-ratio P value from linear regression
    # (Methods, 'Univariable analyses'), not estimated jointly with the
    # regression coefficients, hence fixed().
    auc_cut1 <- fixed(1495);  label("Lower OC AUC0-24 cutoff separating the low and middle exposure groups (ng*h/mL)")   # Rayner 2013 Table 6, composite symptom score AUC rows: reference <= 1,495
    auc_cut2 <- fixed(14497); label("Upper OC AUC0-24 cutoff separating the middle and high exposure groups (ng*h/mL)")  # Rayner 2013 Table 6, composite symptom score AUC rows: > 14,497

    # ----- Reference IC50 -----
    ic50_ref <- fixed(0.18); label("Reference neuraminidase-inhibition IC50, influenza A/Texas/36/91 (H1N1) (nM)")  # Rayner 2013 Table 6, IC50 (nM) reference group 0.18; assay value 0.18 +/- 0.11 nM per Methods

    # ----- Estimated regression coefficients (Table 6) -----
    # Both AUC0-24 coefficients are negative: higher OC exposure lowers
    # the symptom burden. The high-vs-middle contrast implied by these
    # two estimates is -7.03 - (-5.36) = -1.67, matching the -1.68
    # printed for that pairwise comparison in the same table.
    e_auc_oselcarb_mid_aucsc  <- -5.36; label("Change in composite symptom score AUC for the middle OC AUC0-24 group vs the low group (score*day)")  # Rayner 2013 Table 6: -5.36 (95% CI -8.73, -1.99), P = 0.0021
    e_auc_oselcarb_high_aucsc <- -7.03; label("Change in composite symptom score AUC for the high OC AUC0-24 group vs the low group (score*day)")   # Rayner 2013 Table 6: -7.03 (95% CI -11.8, -2.27), P = 0.0042

    # Positive: subjects inoculated with influenza B/Yamagata (the
    # 100-fold less NA-inhibitor-susceptible strain) had a slightly
    # higher symptom burden, but the effect is not significant
    # (P = 0.257) and its CI spans zero. It is retained because IC50 is
    # collinear with OC AUC0-24 across the two studies.
    e_ic50_neuraminidase_aucsc <- 1.77; label("Change in composite symptom score AUC for influenza B/Yamagata (IC50 16.76 nM) vs influenza A/Texas (IC50 0.18 nM) (score*day)")  # Rayner 2013 Table 6: 1.77 (95% CI -1.31, 4.85), P = 0.257

    # No IIV: the analysis is a subject-level linear regression with one
    # AUCSC observation per subject, so there is no random effect to
    # estimate.

    # No residual error: the paper reports no residual standard
    # deviation, residual sum of squares or R^2 for this regression
    # (Tables 3, 5 and 6 give coefficients, confidence intervals,
    # P values and AICc only). Per the house convention for an
    # exposure-response model with no published variance component, no
    # observation endpoint is declared and no placeholder residual is
    # invented.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Recover the paper's 3-group OC AUC0-24 categorical variable
    #    from the continuous exposure column. The two indicators are
    #    mutually exclusive; when both are 0 the subject is in the
    #    reference (lowest-exposure) group, which for this endpoint is
    #    OC AUC0-24 <= 1,495 ng*h/mL and is almost entirely placebo.
    #    Boundary convention follows Table 6 exactly: the middle group
    #    is "> auc_cut1 to <= auc_cut2".
    # ------------------------------------------------------------------
    aucMid  <- (AUC_OSELCARB > auc_cut1) * (AUC_OSELCARB <= auc_cut2)
    aucHigh <- (AUC_OSELCARB > auc_cut2)

    # ------------------------------------------------------------------
    # 2. Recover the two-level IC50 categorical variable. The paper
    #    fitted IC50 as a categorical contrast between its only two
    #    observed values, so any value above the printed 0.18 nM
    #    reference selects the influenza B/Yamagata level.
    # ------------------------------------------------------------------
    highIc50 <- (IC50_NEURAMINIDASE > ic50_ref)

    # ------------------------------------------------------------------
    # 3. Exposure contrast against the reference cohort. This is the
    #    quantity the paper actually estimates and the one to gate
    #    against Table 6; it is free of the unreported intercept.
    # ------------------------------------------------------------------
    d_aucsc <- e_auc_oselcarb_mid_aucsc * aucMid +
      e_auc_oselcarb_high_aucsc * aucHigh

    # ------------------------------------------------------------------
    # 4. Predicted composite symptom score AUC. Carries the reference-
    #    level caveat documented in ini(): the contrasts are exact, the
    #    absolute level is anchored on the Table 4 subgroup mean.
    # ------------------------------------------------------------------
    aucsc <- aucsc_int + d_aucsc + e_ic50_neuraminidase_aucsc * highIc50
  })
}
