Loprete_2016_safinamide_ontime <- function() {
  description <- paste0(
    "Linear disease-progression and treatment-effect model for daily ON-time (ON-time ",
    "plus ON-time with minor dyskinesia, recorded in the Hauser patient diary over an ",
    "18 h window) in patients with Parkinson's disease and motor fluctuations on stable ",
    "levodopa, from the safinamide phase 3 Study 016 (Loprete 2016 model 915). ON-time ",
    "is the individual observed baseline plus an intercept offset plus a slope times ",
    "study time in 4-week months: ontime = ONTIME_BL + (intplac + inttreat * ",
    "ON_TREATMENT) + slopplac * time. The placebo intercept offset was small, poorly ",
    "estimated and set to zero; the safinamide effect on the slope was also set to zero, ",
    "so safinamide acts purely as an instantaneous 0.728 h upward shift attained by the ",
    "first post-baseline visit at week 4, after which ON-time rises with the same slope ",
    "as placebo. Interindividual variability is additive on the intercept and on the ",
    "slope; the baseline carries no random effect because individual observed baselines ",
    "are supplied as data. Residual error is additive. This model has NO drug PK input: ",
    "safinamide average 24 h concentration, safinamide dose, age and levodopa exposure ",
    "were all screened as covariates and none improved the fit, so the drug effect is ",
    "carried entirely by a treatment-arm indicator. It is therefore a companion to, but ",
    "not coupled with, the population PK model in the same paper ",
    "(modellib('Loprete_2016_safinamide'))."
  )
  reference <- paste(
    "Loprete L, Leuratti C, Cattaneo C, Thapar MM, Farrell C, Sardina M.",
    "Population pharmacokinetic and pharmacodynamic analyses of safinamide",
    "in subjects with Parkinson's disease.",
    "Pharmacol Res Perspect 2016;4(5):e00251.",
    "doi:10.1002/prp2.251.",
    sep = " "
  )
  vignette <- "Loprete_2016_safinamide"

  # Paper-specific etas: Loprete 2016 equation 8 places the random effects
  # additively on the linear-scale parameters (P_i = TVP + eta_i), and the
  # Results state that one additive term describes the variability on the
  # sum of INTTREAT and INTPLAC and another the variability on the sum of
  # SLOPTREAT and SLOPPLAC. Because intercept and slope are signed
  # linear-scale parameters rather than log-transformed positive ones, the
  # canonical eta + l<param> pairing does not apply; declare the eta names
  # as paper-specific so checkModelConventions() accepts them. Same pattern
  # as Lee_2011_parkinson_progression.R.
  paper_specific_etas <- c("etaintercept", "etaslope")

  units <- list(
    time = paste0(
      "month (4 weeks = 672 h; the source dataset carries TIME in hours from the ",
      "baseline visit while SLOPPLAC is per month, and the paper's own arithmetic ",
      "fixes the month at 4 weeks -- see description and the validation vignette)"
    ),
    dosing = "(none; the drug effect is encoded via the ON_TREATMENT covariate, not via drug input)",
    concentration = "(daily ON-time in hours; observation variable ontime)"
  )

  covariateData <- list(
    ONTIME_BL = list(
      description = paste0(
        "Per-subject baseline daily ON-time, i.e. the sum of ON-time and ON-time with ",
        "minor dyskinesia averaged over 2 to 5 diary recording days during the 18 h ",
        "(0600-2400) recording window, at the week 1 baseline visit."
      ),
      units = "hr",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed per subject and entered UNCENTRED as an additive offset, not as a ",
        "normalised power term. Loprete 2016 Results, PKPD analysis and final model: 'the ",
        "typical value for BL was not estimated. Instead, BL utilized the individual ",
        "observed ON-time baseline values'; Table 4 accordingly reports the ",
        "interindividual variance on BL as 0 FIX. The model therefore predicts a change ",
        "from each subject's own observed baseline, and ONTIME_BL must be supplied in the ",
        "event data. Loprete 2016 does not tabulate the cohort baseline ON-time ",
        "distribution, so no reference value is available from this source; the validation ",
        "vignette documents the value used for simulation and where it came from."
      ),
      source_name = "BL"
    ),
    ON_TREATMENT = list(
      description = paste0(
        "Binary treatment-arm indicator. 1 = randomized to safinamide (50 or 100 mg/day in ",
        "Study 016); 0 = randomized to placebo."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm)",
      notes = paste0(
        "Time-fixed per subject; Study 016 was a randomized parallel-group trial. The ",
        "source column TRT is categorical and distinguishes the 50 mg and 100 mg ",
        "safinamide arms from placebo, but the final model does not separate the two ",
        "active dose levels: Loprete 2016 screened the categorical safinamide dose (and ",
        "the individual average 24 h safinamide concentration SAAV) on both the intercept ",
        "and the slope and retained neither, reporting that 'there was very little ",
        "difference between the safinamide treatment arms'. The pooled active-versus-",
        "placebo switch is therefore the whole drug-effect model here, which is exactly ",
        "the ON_TREATMENT use case. Applied at all times including TIME = 0; see the ",
        "validation vignette's Assumptions and deviations section for what that implies ",
        "at the baseline visit."
      ),
      source_name = "TRT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste0(
        "Screened univariately on SLOPPLAC, INTTREAT and INTPLAC with a 60-year reference ",
        "and not retained (Loprete 2016 Results, PKPD analysis and final model). No point ",
        "estimate is published."
      ),
      source_name = "AGE"
    ),
    DOSE_LEVODOPA_MGD = list(
      description = "Levodopa dose rate",
      units = "mg/24 h",
      type = "continuous",
      notes = paste0(
        "Screened as LEVO (levodopa dose rate per 24 h at each visit, reference ",
        "500 mg/24 h) and LEVR (rate of change from baseline) on SLOPPLAC, INTTREAT and ",
        "INTPLAC, and not retained. Model 931 put LEVO on SLOPPLAC with an effect of ",
        "-0.779 and an OFV drop of 7.6 points, but the authors traced that signal to 8 ",
        "patients with LEVO < 200 mg/24 h and found the effect disappeared without them, ",
        "so it was rejected. Not registered in inst/references/covariate-columns.md ",
        "because it is screened-only and never referenced in model()."
      ),
      source_name = "LEVO / LEVR"
    ),
    AUC_SAFINAMIDE = list(
      description = "Individual average safinamide plasma concentration over 24 h",
      units = "ng/mL",
      type = "continuous",
      notes = paste0(
        "The paper's exposure metric SAAV, computed as the safinamide dose at each visit ",
        "divided by the individual clearance estimate from the companion population PK ",
        "analysis. Screened on SLOPPLAC and INTTREAT with a stated reference of 13 ng/mL ",
        "(the reported median) and not retained: Loprete 2016 Results record that the 95% ",
        "confidence intervals of the SAAV effect included the null and that Figure 5 shows ",
        "no exposure-response trend. This is the parameter that severs the PK and PD ",
        "layers of the paper: the final PD model has no exposure term. The printed 13 ",
        "ng/mL reference is not reproducible from the paper's own PK estimates -- a ",
        "100 mg/day dose divided by a CL/F of 4.96 L/h gives an average concentration near ",
        "840 ng/mL -- so the reference value appears to carry a unit or scale error; this ",
        "has no effect on the final model, which contains no SAAV term. Not registered in ",
        "inst/references/covariate-columns.md because it is screened-only and never ",
        "referenced in model()."
      ),
      source_name = "SAAV"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 668L,
    n_studies = 1L,
    age_range = "34.0-80.0 years",
    age_median = "60.0 years",
    weight_range = "33.5-120 kg",
    weight_median = "62.0 kg",
    sex_female_pct = 28.1,
    race_ethnicity = c(Asian = 80.5, White = 19.3, Other = 0.1),
    disease_state = paste0(
      "Parkinson's disease with motor fluctuations, on a stable dose of levodopa for at ",
      "least 4 weeks before baseline (phase 3 Study 016)."
    ),
    dose_range = paste0(
      "Placebo, safinamide 50 mg/day, or safinamide 100 mg/day, for 24 weeks; the 100 ",
      "mg/day dose could be reduced to 50 mg/day for intolerance. The final model does not ",
      "distinguish the two active dose levels -- see the ON_TREATMENT covariate notes."
    ),
    regions = "Not reported in Loprete 2016; the cohort was 80.5% Asian.",
    notes = paste0(
      "Demographics from Loprete 2016 Table 2 (Study 016, n = 668). The analysis dataset ",
      "held 669 patients contributing 3607 observation records; 668 patients contributing ",
      "3603 records entered the population PKPD analysis. Observations were taken at ",
      "weeks 1 (baseline), 4, 8, 12, 18 and 24, each the average over 2 to 5 diary ",
      "recording days. Levodopa dose was largely stable: 11/668 (2%) of patients increased ",
      "it during the trial and 51/668 (8%) decreased it. Table 2 prints the 'Other' race ",
      "cell as 1 (< 0.01), which appears to be a fraction rather than a percentage; 1/668 ",
      "is 0.1%, the value recorded here."
    )
  )

  ini({
    # Population-mean structural parameters from Loprete 2016 Table 4 (final
    # model 915). All three are signed linear-scale parameters on the ON-time
    # axis, so none is log-transformed: the source places the random effects
    # additively (equation 8, P_i = TVP + eta_i) rather than exponentially.
    #
    # Time is in 4-week months. The source dataset carries TIME in hours from
    # the baseline visit, but SLOPPLAC is reported per month, and the paper's
    # own arithmetic pins the month at 4 weeks: the Discussion states that a
    # slope of 0.117 h/month 'would lead to an average increase of 0.70 h over
    # the duration of the study', and 0.70 / 0.117 = 5.98, i.e. the 24-week
    # study is 6 months, so 1 month = 4 weeks = 672 h. The same arithmetic
    # reproduces the paper's quoted 1.43 h total safinamide increase at the
    # end of the trial: 0.728 + 0.117 * 6 = 1.430 h.

    intplac <- fixed(0)
    label("Placebo offset on the ON-time intercept (h)")
    # Table 4 row 'Intercept (INT PLAC) (h)': 0 FIX. Model 903 estimated it at 0.231 h (%RSE 62.3, 95% CI -0.0512 to 0.513) and model 907 at 0.238 h (95% CI -0.0207 to 0.497); both intervals included the null, so model 915 set it to zero, which raised the OFV by 3 points.

    inttreat <- 0.728
    label("Safinamide instantaneous effect on the ON-time intercept (h)")
    # Table 4 row 'INT TREAT (h)': 0.728, %RSE 15.0, 95% CI 0.514-0.942. Bootstrap 0.727 (95% CI 0.505-0.943). Applied only when ON_TREATMENT = 1.

    slopplac <- 0.117
    label("Rate of ON-time change, shared by the placebo and safinamide arms (h/month)")
    # Table 4 row 'Slope (SLOP PLAC) (h/month)': 0.117, %RSE 17.0, 95% CI 0.0780-0.156. Bootstrap 0.116 (95% CI 0.0769-0.154). Table 4 caption: 'SLOPPLAC, slope of the effect, the same for placebo and Safinamide treatments'; SLOPTREAT was set to zero in model 908 with no change in OFV.

    # Additive interindividual variability, Loprete 2016 equation 8. The
    # Results state that one additive term carries the variability on the sum
    # of INTTREAT and INTPLAC and another the variability on the sum of
    # SLOPTREAT and SLOPPLAC, so the etas sit on the combined intercept and
    # on the combined slope rather than on the placebo components alone. The
    # Table 4 'x2' column holds variances: its 'SD' column is their square
    # root (sqrt(3.29) = 1.81 and sqrt(0.130) = 0.361, both matching the
    # printed values). There is no random effect on the baseline: Table 4
    # reports 'Interindividual variability x2 BL' as 0 FIX because the
    # individual observed baselines are supplied as data.
    etaintercept ~ 3.29
    # Table 4 row 'x2 INT PLAC (h)': 3.29, %RSE 7.93, 95% CI 2.78-3.80, SD 1.81. Bootstrap 3.29 (95% CI 2.79-3.81). Carries the variability on INTPLAC + INTTREAT per the Results text.
    etaslope ~ 0.130
    # Table 4 row 'x2 SLOP PLAC (h/month)': 0.130, %RSE 10.5, 95% CI 0.103-0.157, SD 0.361. Bootstrap 0.129 (95% CI 0.104-0.156). Shrinkage 25% on SLOPPLAC and 12% on INTPLAC.

    # Additive residual error on the ON-time scale, Loprete 2016 equation 9.
    # Table 4 reports the variance and its square root as the SD column
    # (sqrt(1.19) = 1.0909, matching the printed 1.09 h).
    addSd <- sqrt(1.19)
    label("Additive residual error standard deviation on daily ON-time (h)")
    # Table 4 row 'Residual variability r2 add': 1.19, %RSE 5.59, 95% CI 1.06-1.32, SD 1.09 h. Bootstrap 1.18 (95% CI 1.05-1.32). addSd = sqrt(1.19) = 1.0909.
  })

  model({
    # Individual intercept and slope. Both random effects are additive on the
    # linear ON-time scale (Loprete 2016 equation 8). The safinamide effect
    # enters the intercept only; the slope is shared with placebo because
    # SLOPTREAT was set to zero in the final model.
    intercept <- intplac + inttreat * ON_TREATMENT + etaintercept
    slope <- slopplac + etaslope

    # Loprete 2016 equation 7, reduced to the final model 915 by INTPLAC = 0
    # and SLOPTREAT = 0, and restated verbatim in the Table 4 caption as
    # 'PD = BL + INTPLAC + INTTREAT + SLOPPLAC * TIMEmonths'. The baseline is
    # the subject's own observed ON-time, supplied as the ONTIME_BL covariate
    # rather than estimated.
    ontime <- ONTIME_BL + intercept + slope * time

    ontime ~ add(addSd)
  })
}
