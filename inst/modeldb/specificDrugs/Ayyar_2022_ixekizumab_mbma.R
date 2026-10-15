Ayyar_2022_ixekizumab_mbma <- function() {
  description <- "MBMA. Dose-based model-based meta-analysis of the week-12 placebo-adjusted PASI75 and PASI90 responder rates for ixekizumab in moderate-to-severe plaque psoriasis, as a sigmoid Emax function of the average weekly dose over the first 12 weeks (Ayyar 2022 Eq. 1, Figure 2). Trial-arm-level trend line with no between-study variability; no PK layer and no rxode2 dose events. The companion TE-based MBMA, which uses predicted skin IL-17A suppression instead of dose, is part of Ayyar_2022_ixekizumab_mpbpk."
  reference <- "Ayyar VS, Lee JB, Wang W, Pryor M, Zhuang Y, Wilde T, Vermeulen A. Minimal Physiologically-Based Pharmacokinetic (mPBPK) Metamodeling of Target Engagement in Skin Informs Anti-IL17A Drug Development in Psoriasis. Front Pharmacol. 2022;13:862291. doi:10.3389/fphar.2022.862291"
  vignette <- "Ayyar_2022_il17a_target_engagement"
  units <- list(
    time = "week (single week-12 landmark read-out; the model has no time dependence)",
    dosing = "mg/week (average weekly ixekizumab dose over the first 12 weeks, supplied in the DOSE_IXEKIZUMAB_MGWK covariate column. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi75_pbo_adj and prob_pasi90_pbo_adj are the STUDY-ARM responder fractions minus the placebo-arm fraction, on a 0-1 scale; neither is a drug concentration)"
  )

  covariateData <- list(
    DOSE_IXEKIZUMAB_MGWK = list(
      description = "Average weekly ixekizumab dose in the study arm over the first 12 weeks of treatment",
      units = "mg/week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TRIAL-ARM-LEVEL. Methods 'Model-Based Meta-Analysis': 'the total dose administered",
        "during the 12 weeks divided by the total duration (12 weeks) to obtain dose per unit",
        "time (mg/week)'. The arm x-positions plotted in Figure 2 sit about 20% below that",
        "definition (e.g. the Phase 3 160 mg then 80 mg q2w arms near 38 mg/week rather",
        "than 560 mg / 12 weeks = 47 mg/week); see the vignette."
      ),
      source_name = "x, dose (mg/week) (Ayyar 2022 Eq. 1, Figure 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 5L,
    age_range = "adults (not reported at arm level)",
    weight_range = "not reported at arm level",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "Moderate-to-severe plaque psoriasis",
    dose_range = "15 mg IV and 5-150 mg SC at weeks 0, 2, 4 (Phase 1); 10-150 mg SC at weeks 0, 2, 4, 8 (Phase 2); 160 mg SC then 80 mg SC q2w or q4w (Phase 3)",
    timepoints = "single landmark read-out 12 weeks after the first dose",
    regions = "international",
    notes = paste(
      "MBMA at the STUDY-ARM level fitted by non-linear least squares (R nls) to the",
      "placebo-adjusted week-12 PASI75 and PASI90 rates of the randomised placebo-controlled",
      "ixekizumab trials: Phase 2 dose-ranging (Leonardi 2012, N = 141, Table 1) and the",
      "Phase 3 UNCOVER-1/2/3 studies (Papp 2018) for the 80 mg q2w and q4w maintenance arms.",
      "n_subjects is NA because the paper does not tabulate the Phase 3 arm sizes. Not for",
      "simulating individual patients."
    )
  )

  ini({
    # Methods Eq. 1: Response = E0 + x^hill * (Emax - E0) / (x^hill + E50^hill), with
    # x = average dose (mg/week). The paper prints no parameter values: each set was
    # digitised by the maintainers from the ixekizumab trend line of Figure 2 (PDF
    # image, about 190 points per panel over 0.5-55 mg/week) and Eq. 1 refitted by
    # nls; the refit reproduces the digitised line with residual SD 0.5 (PASI75) and
    # 0.4 (PASI90) percentage points. Below 3 mg/week the plotted line is a polyline
    # through a 0.5 mg/week simulation grid, so only its vertices (0.5, 1, 1.5, 2,
    # 2.5 mg/week) were used there. E0 is not identified by the plotted range (the
    # line starts at 0.5 mg/week) and absorbs the steep low-dose rise, hence its
    # negative value; use the model only over the plotted dose range.
    e0_pasi75 <- -7.944; label("Dose-MBMA PASI75 response at zero dose E0 (% placebo-adjusted)") # digitised Figure 2 left panel (ixekizumab line), Eq. 1 refit
    emax_pasi75 <- 82.12; label("Dose-MBMA PASI75 maximum response Emax (% placebo-adjusted)") # digitised Figure 2 left panel (ixekizumab line), Eq. 1 refit
    led50_pasi75 <- log(2.605); label("Dose-MBMA PASI75 average weekly dose for half-maximal response E50 (mg/week)") # digitised Figure 2 left panel (ixekizumab line), Eq. 1 refit
    hill_pasi75 <- 1.653; label("Dose-MBMA PASI75 Hill coefficient (unitless)") # digitised Figure 2 left panel (ixekizumab line), Eq. 1 refit
    e0_pasi90 <- -6.973; label("Dose-MBMA PASI90 response at zero dose E0 (% placebo-adjusted)") # digitised Figure 2 right panel (ixekizumab line), Eq. 1 refit
    emax_pasi90 <- 68.56; label("Dose-MBMA PASI90 maximum response Emax (% placebo-adjusted)") # digitised Figure 2 right panel (ixekizumab line), Eq. 1 refit
    led50_pasi90 <- log(2.888); label("Dose-MBMA PASI90 average weekly dose for half-maximal response E50 (mg/week)") # digitised Figure 2 right panel (ixekizumab line), Eq. 1 refit
    hill_pasi90 <- 1.474; label("Dose-MBMA PASI90 Hill coefficient (unitless)") # digitised Figure 2 right panel (ixekizumab line), Eq. 1 refit

    # Residual error: unweighted nls on the response scale; residual SE not reported
    addSd_prob_pasi75_pbo_adj <- fixed(0); label("Additive residual error, placebo-adjusted PASI75 fraction (not reported)") # Methods: nls fit; residual SE not reported
    addSd_prob_pasi90_pbo_adj <- fixed(0); label("Additive residual error, placebo-adjusted PASI90 fraction (not reported)") # Methods: nls fit; residual SE not reported
  })

  model({
    ed50_pasi75 <- exp(led50_pasi75)
    ed50_pasi90 <- exp(led50_pasi90)
    xdose <- DOSE_IXEKIZUMAB_MGWK
    prob_pasi75_pbo_adj <- (e0_pasi75 + xdose^hill_pasi75 * (emax_pasi75 - e0_pasi75) / (xdose^hill_pasi75 + ed50_pasi75^hill_pasi75)) / 100
    prob_pasi90_pbo_adj <- (e0_pasi90 + xdose^hill_pasi90 * (emax_pasi90 - e0_pasi90) / (xdose^hill_pasi90 + ed50_pasi90^hill_pasi90)) / 100

    prob_pasi75_pbo_adj ~ add(addSd_prob_pasi75_pbo_adj)
    prob_pasi90_pbo_adj ~ add(addSd_prob_pasi90_pbo_adj)
  })
}
