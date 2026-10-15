Ayyar_2022_secukinumab_mbma <- function() {
  description <- "MBMA. Dose-based model-based meta-analysis of the week-12 placebo-adjusted PASI75 and PASI90 responder rates for secukinumab in moderate-to-severe plaque psoriasis, as a sigmoid Emax function of the average weekly dose over the first 12 weeks (Ayyar 2022 Eq. 1, Figure 2). Trial-arm-level trend line with no between-study variability; no PK layer and no rxode2 dose events. The companion TE-based MBMA, which uses predicted skin IL-17A suppression instead of dose, is part of Ayyar_2022_secukinumab_mpbpk."
  reference <- "Ayyar VS, Lee JB, Wang W, Pryor M, Zhuang Y, Wilde T, Vermeulen A. Minimal Physiologically-Based Pharmacokinetic (mPBPK) Metamodeling of Target Engagement in Skin Informs Anti-IL17A Drug Development in Psoriasis. Front Pharmacol. 2022;13:862291. doi:10.3389/fphar.2022.862291"
  vignette <- "Ayyar_2022_il17a_target_engagement"
  units <- list(
    time = "week (single week-12 landmark read-out; the model has no time dependence)",
    dosing = "mg/week (average weekly secukinumab dose over the first 12 weeks, supplied in the DOSE_SECUKINUMAB_MGWK covariate column. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi75_pbo_adj and prob_pasi90_pbo_adj are the STUDY-ARM responder fractions minus the placebo-arm fraction, on a 0-1 scale; neither is a drug concentration)"
  )

  covariateData <- list(
    DOSE_SECUKINUMAB_MGWK = list(
      description = "Average weekly secukinumab dose in the study arm over the first 12 weeks of treatment",
      units = "mg/week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TRIAL-ARM-LEVEL. Methods 'Model-Based Meta-Analysis': 'the total dose administered",
        "during the 12 weeks divided by the total duration (12 weeks) to obtain dose per unit",
        "time (mg/week)'. The arm x-positions plotted in Figure 2 sit about 25% below that",
        "definition for the q4w regimens (e.g. the 300 mg Phase 3 arms near 110 mg/week",
        "rather than 1800 mg / 12 weeks = 150 mg/week); see the vignette."
      ),
      source_name = "x, dose (mg/week) (Ayyar 2022 Eq. 1, Figure 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2549L,
    n_studies = 7L,
    age_range = "adults (not reported at arm level)",
    weight_range = "not reported at arm level",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "Moderate-to-severe plaque psoriasis",
    dose_range = "3 and 10 mg/kg IV; 25-300 mg SC single, q4w, weeks 0, 1, 2, 4, or weeks 0, 1, 2, 3, 4 then q4w (Table 1)",
    timepoints = "single landmark read-out 12 weeks after the first dose",
    regions = "international",
    notes = paste(
      "MBMA at the STUDY-ARM level fitted by non-linear least squares (R nls) to the",
      "placebo-adjusted week-12 PASI75 and PASI90 rates of the randomised placebo-controlled",
      "secukinumab trials in Table 1 (Phase 2 proof-of-concept, low and high dose-ranging,",
      "regimen-finding; Phase 3 ERASURE, FIXTURE, FEATURE). n_subjects sums the Table 1 N",
      "column, which excludes active-comparator arms. Not for simulating individual patients."
    )
  )

  ini({
    # Methods Eq. 1: Response = E0 + x^hill * (Emax - E0) / (x^hill + E50^hill), with
    # x = average dose (mg/week). The paper prints no parameter values: each set was
    # digitised by the maintainers from the secukinumab trend line of Figure 2 (PDF
    # image, about 400 points per panel over 0.6-300 mg/week) and Eq. 1 refitted by
    # nls; the refit reproduces the digitised line with residual SD 0.9 (PASI75) and
    # 0.7 (PASI90) percentage points. E0 is not identified by the plotted range and
    # absorbs curvature at low dose, hence its small negative value.
    e0_pasi75 <- -1.245; label("Dose-MBMA PASI75 response at zero dose E0 (% placebo-adjusted)") # digitised Figure 2 left panel (secukinumab line), Eq. 1 refit
    emax_pasi75 <- 75.82; label("Dose-MBMA PASI75 maximum response Emax (% placebo-adjusted)") # digitised Figure 2 left panel (secukinumab line), Eq. 1 refit
    led50_pasi75 <- log(16.94); label("Dose-MBMA PASI75 average weekly dose for half-maximal response E50 (mg/week)") # digitised Figure 2 left panel (secukinumab line), Eq. 1 refit
    hill_pasi75 <- 1.898; label("Dose-MBMA PASI75 Hill coefficient (unitless)") # digitised Figure 2 left panel (secukinumab line), Eq. 1 refit
    e0_pasi90 <- -0.821; label("Dose-MBMA PASI90 response at zero dose E0 (% placebo-adjusted)") # digitised Figure 2 right panel (secukinumab line), Eq. 1 refit
    emax_pasi90 <- 66.61; label("Dose-MBMA PASI90 maximum response Emax (% placebo-adjusted)") # digitised Figure 2 right panel (secukinumab line), Eq. 1 refit
    led50_pasi90 <- log(36.91); label("Dose-MBMA PASI90 average weekly dose for half-maximal response E50 (mg/week)") # digitised Figure 2 right panel (secukinumab line), Eq. 1 refit
    hill_pasi90 <- 1.868; label("Dose-MBMA PASI90 Hill coefficient (unitless)") # digitised Figure 2 right panel (secukinumab line), Eq. 1 refit

    # Residual error: unweighted nls on the response scale; residual SE not reported
    addSd_prob_pasi75_pbo_adj <- fixed(0); label("Additive residual error, placebo-adjusted PASI75 fraction (not reported)") # Methods: nls fit; residual SE not reported
    addSd_prob_pasi90_pbo_adj <- fixed(0); label("Additive residual error, placebo-adjusted PASI90 fraction (not reported)") # Methods: nls fit; residual SE not reported
  })

  model({
    ed50_pasi75 <- exp(led50_pasi75)
    ed50_pasi90 <- exp(led50_pasi90)
    xdose <- DOSE_SECUKINUMAB_MGWK
    prob_pasi75_pbo_adj <- (e0_pasi75 + xdose^hill_pasi75 * (emax_pasi75 - e0_pasi75) / (xdose^hill_pasi75 + ed50_pasi75^hill_pasi75)) / 100
    prob_pasi90_pbo_adj <- (e0_pasi90 + xdose^hill_pasi90 * (emax_pasi90 - e0_pasi90) / (xdose^hill_pasi90 + ed50_pasi90^hill_pasi90)) / 100

    prob_pasi75_pbo_adj ~ add(addSd_prob_pasi75_pbo_adj)
    prob_pasi90_pbo_adj ~ add(addSd_prob_pasi90_pbo_adj)
  })
}
