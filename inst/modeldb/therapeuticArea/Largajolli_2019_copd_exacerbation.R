Largajolli_2019_copd_exacerbation <- function() {
  description <- paste(
    "Two-state semi-Markov model of COPD exacerbations (Largajolli 2019):",
    "each subject alternates between an asymptomatic state 1 and a",
    "symptomatic (exacerbation) state 2, and the hazard of leaving a state",
    "is log-logistic in the time since entry into that state (the clock",
    "restarts at every transition). Log-logistic hazard",
    "h(t) = (shape / scale) * (t / scale)^(shape - 1) / (1 + (t / scale)^shape),",
    "survival S(t) = 1 / (1 + (t / scale)^shape), scale in days. Disease",
    "stage (ATS 1995 moderate or severe, FEV1 below 50 percent of",
    "predicted, versus mild) shortens both sojourns and the time to",
    "dropout; an active-treatment indicator lengthens the asymptomatic",
    "sojourn (default values are salmeterol 50 ug twice daily versus",
    "placebo, Data 1). A separate log-logistic time-to-dropout model in",
    "time since randomisation describes missing-at-random dropout. Fitted",
    "by Bayesian MCMC in OpenBUGS to a 52-week phase III trial (n = 619).",
    "Closed-form outputs only, no ODE: haz12 / sur1 (exacerbation onset),",
    "haz21 / sur2 (exacerbation resolution), haz_drop / sur_drop. The",
    "transition scales were not printed and were digitised from the",
    "supplementary figures, and the between-subject variances were not",
    "reported (carried at zero)."
  )
  reference <- paste(
    "Largajolli A, Beerahee M, Yang S.",
    "Bayesian approach to investigate a two-state mixed model of COPD",
    "exacerbations.",
    "J Pharmacokinet Pharmacodyn. 2019;46(4):371-384.",
    "doi:10.1007/s10928-019-09643-6. PMCID: PMC6848253.",
    "Correction (licence change to open access only, no values changed):",
    "J Pharmacokinet Pharmacodyn. 2019;46(6):627.",
    "doi:10.1007/s10928-019-09664-1.",
    "Shapes from Table 2, covariate effects from Table 3 and the Results",
    "text, drug effects from Table 4; transition scales digitised from",
    "Supplementary Figure S4 and the dropout shape and scale from the",
    "Supplementary Figure S3 trace plots."
  )
  vignette <- "Largajolli_2019_copd_exacerbation"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; treatment enters as the ON_TREATMENT indicator)",
    concentration = "n/a (outputs are hazards in 1/day and survival probabilities)"
  )

  covariateData <- list(
    FEV1_PCTPRED = list(
      description = "Baseline FEV1 as percent of the predicted value, used to assign the ATS 1995 COPD disease stage",
      units = "% predicted",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the binary ATS 1995 stage indicator",
        "copd_nonmild = (FEV1_PCTPRED < 50): ATS 1995 grades FEV1 >= 50",
        "percent predicted as mild, 35-50 as moderate and below 35 as",
        "severe (Largajolli 2019 Table 1 footnote), and the covariate",
        "analysis tested disease stage as mild versus the rest (Methods,",
        "Two-state semi Markov model). The indicator is centred at the",
        "Data 1 non-mild fraction 0.5428 = (224 + 112) / 619 (Table 1)",
        "because the transition scales come from the covariate-free base",
        "model; see the vignette. Data 1 used the ATS staging; Data 2 and",
        "Data 3 used GOLD, which the model does not reproduce."
      ),
      source_name = "disease stage (ATS)"
    ),
    ON_TREATMENT = list(
      description = "Active-treatment arm indicator (1 = active, 0 = reference arm)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Default ini() values are the Data 1 contrast, salmeterol 50 ug",
        "inhaled twice daily (1) versus placebo (0), parallel-group and",
        "time-fixed (Largajolli 2019 Table 4). Multiplies the log-logistic",
        "scales by exp(e_trt_scale12 * (ON_TREATMENT - 0.5089)) and",
        "exp(e_trt_scale21 * (ON_TREATMENT - 0.5089)), centred at the Data",
        "1 salmeterol fraction 315 / 619 (Table 1). Table 4 also gives the",
        "contrasts fluticasone furoate / vilanterol versus vilanterol (Data",
        "2) and losmapimod 15 mg versus placebo and lower doses (Data 3);",
        "those were estimated against each data set's own unreported",
        "baseline scales, so substituting them here transfers a relative",
        "effect onto the Data 1 baseline (see the vignette)."
      ),
      source_name = "drug effect (placebo vs Salmeterol)"
    ),
    T_STATE_ENTRY = list(
      description = "Study time at which the subject entered its current state",
      units = "day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the sojourn clock of the semi-Markov process:",
        "tsoj = time - T_STATE_ENTRY. Time-varying: set it to the study",
        "day of the latest observed transition (0 before the first one)",
        "and carry it forward. The default of 0 makes time itself the",
        "sojourn clock, which gives the hazard of a single sojourn as a",
        "function of time since entry (Supplementary Figure S4). The",
        "dropout hazard ignores it and runs on time since randomisation."
      ),
      source_name = "time since entry into the state (Methods, Two-state semi Markov model)"
    )
  )

  covariatesDataExcluded <- list(
    SMOKE = list(
      description = "Smoking status (ex-smoker versus the rest)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on both scales in Data 1 and not retained: DIC 7025 versus 7026 for the base model; beta 0.07 (95% CI -0.30 to 0.42) on the state 1 scale and 0.26 (0.08 to 0.44) on the state 2 scale (Largajolli 2019 Table 3)."
    ),
    SMOKE_PACKYEARS = list(
      description = "Cumulative cigarette pack-years (current smokers)",
      units = "pack-years",
      type = "continuous",
      notes = "Screened and not retained: DIC 7025; beta -0.002 (-0.01 to 0.01) on the state 1 scale and 0.005 (0.001 to 0.01) on the state 2 scale (Largajolli 2019 Table 3)."
    ),
    SEASON_AUTUMN_WINTER = list(
      description = "Baseline season, spring-summer versus autumn-winter",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained: DIC 7028; beta -0.08 (-0.44 to 0.29) on the state 1 scale and -0.06 (-0.25 to 0.12) on the state 2 scale (Largajolli 2019 Table 3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 619L,
    n_studies = 1L,
    age_range = "Adults; not reported in Largajolli 2019.",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "COPD with at least 10 pack-years of smoking and at least one acute",
      "exacerbation per year in the previous 3 years; ATS 1995 stage mild",
      "283, moderate 224, severe 112 (Table 1)."
    ),
    dose_range = "Salmeterol 50 ug inhaled twice daily (n = 315) or placebo (n = 304) for 52 weeks.",
    regions = "Not reported.",
    notes = paste(
      "Data 1 of Largajolli 2019: the salmeterol-alone and placebo arms of",
      "the 52-week phase III trial SFCB3024, which also had salmeterol /",
      "fluticasone propionate and fluticasone-alone arms that were not",
      "analysed. Current smokers 298, ex-smokers 314, unknown 7; mean",
      "pack-years in current smokers 45.9 (salmeterol) and 41.8 (placebo);",
      "enrolled in autumn-winter 275, spring-summer 344 (Table 1). Base",
      "model, covariate and dropout analyses used Data 1 only. Drug effects",
      "were also estimated in Data 2 (3238 subjects, two replicate 52-week",
      "fluticasone furoate / vilanterol trials HZC102970 and HZC102871)",
      "and Data 3 (602 subjects, 24-week losmapimod phase 2b MKI113006)."
    )
  )

  ini({
    # Transition scales. Not printed: digitised by the maintainers from
    # Supplementary Figure S4 (base log-logistic model, typical transition
    # rate against time since entry into the state, box medians at 6 and 10
    # times) and fitted with the shape held at Table 2. A free fit of both
    # scale and shape to the same points returns shape 1.345 and 2.962,
    # against Table 2's 1.34 and 2.97, which confirms the reading.
    lscale12 <- log(222); label("Log-logistic scale of the asymptomatic sojourn, state 1 to 2 (day)") # Digitised: Supplementary Fig S4 left panel, least-squares fit with shape 1.34
    lshape12 <- log(1.34); label("Log-logistic shape of the asymptomatic sojourn, state 1 to 2 (unitless)") # Table 2 Loglogistic Alpha1 median 1.34 (95% CI 1.22-1.48)
    lscale21 <- log(11.6); label("Log-logistic scale of the exacerbation sojourn, state 2 to 1 (day)") # Digitised: Supplementary Fig S4 right panel, least-squares fit with shape 2.97
    lshape21 <- log(2.97); label("Log-logistic shape of the exacerbation sojourn, state 2 to 1 (unitless)") # Table 2 Loglogistic Alpha2 median 2.97 (95% CI 2.67-3.27)

    # Dropout. Not printed: read by the maintainers off the final-model
    # trace plots (Supplementary Figure S3, alpha[3] and logtheta0[3]) as the
    # median of the plotted chain. The same reading reproduces the printed
    # values of every chain that has one (alpha1 1.345, log alpha2 1.084,
    # beta1[2] -0.110, beta1[3] -0.682).
    lscale_drop <- 7.02; label("Log-logistic scale of the time to dropout for a mild-stage subject (day)") # Digitised: Supplementary Fig S3 logtheta0[3] chain median 7.02
    lshape_drop <- log(0.81); label("Log-logistic shape of the time to dropout (unitless)") # Digitised: Supplementary Fig S3 alpha[3] chain median 0.81

    # Disease stage (ATS 1995 moderate or severe versus mild), log-scale
    # coefficients on the scale (Eq. 6). Negative values shorten the sojourn.
    e_copdnonmild_scale12 <- -0.40; label("Disease-stage coefficient on the log state 1 scale (unitless)") # Table 3 Loglogistic + disease stage (ATS) Beta11 median -0.40 (95% CI -0.78 to -0.05)
    e_copdnonmild_scale21 <- -0.11; label("Disease-stage coefficient on the log state 2 scale (unitless)") # Table 3 Loglogistic + disease stage (ATS) Beta12 median -0.11 (95% CI -0.29 to 0.09)
    e_copdnonmild_scale_drop <- -0.68; label("Disease-stage coefficient on the log dropout scale (unitless)") # Results, Drop out model: median -0.68 (95% CI -1.28 to -0.1)

    # Drug effect, Data 1 (salmeterol versus placebo), log-scale coefficients
    # on the scale (Eq. 6). Positive values lengthen the sojourn.
    e_trt_scale12 <- 0.310; label("Treatment coefficient on the log state 1 scale (unitless)") # Table 4 Data 1 Beta11 median 0.310 (95% CI 0.046-0.559)
    e_trt_scale21 <- -0.031; label("Treatment coefficient on the log state 2 scale (unitless)") # Table 4 Data 1 Beta12 median -0.031 (95% CI -0.162 to 0.104)

    # Random effects on the two transition scales were retained (Results,
    # Model selection) but their variances were not reported.
    etalscale12 ~ fixed(0) # random effect on the state 1 scale retained; variance not reported
    etalscale21 ~ fixed(0) # random effect on the state 2 scale retained; variance not reported
  })
  model({
    # ATS 1995 stage: mild is FEV1 >= 50 percent predicted (Table 1 footnote)
    copd_nonmild <- 0
    if (FEV1_PCTPRED < 50) copd_nonmild <- 1

    # Eq. 6 scale model. The covariate-free scales come from the base model
    # (Supplementary Fig S4), so the covariates are centred at their Data 1
    # cohort means: non-mild 336 / 619 = 0.5428, salmeterol 315 / 619 =
    # 0.5089 (Table 1). The study-length term beta_k0 * log(T_j) of Eq. 6
    # was not reported and is absorbed into the digitised scales.
    scale12 <- exp(lscale12 + etalscale12 +
      e_copdnonmild_scale12 * (copd_nonmild - 0.5428) +
      e_trt_scale12 * (ON_TREATMENT - 0.5089))
    shape12 <- exp(lshape12)
    scale21 <- exp(lscale21 + etalscale21 +
      e_copdnonmild_scale21 * (copd_nonmild - 0.5428) +
      e_trt_scale21 * (ON_TREATMENT - 0.5089))
    shape21 <- exp(lshape21)
    scale_drop <- exp(lscale_drop + e_copdnonmild_scale_drop * copd_nonmild)
    shape_drop <- exp(lshape_drop)

    # Semi-Markov clock: time since entry into the current state
    tsoj <- t - T_STATE_ENTRY
    if (tsoj < 0) tsoj <- 0

    # Eq. 5 log-logistic hazard and survival. The 1e-6 day offset keeps the
    # dropout hazard finite at time 0, where a shape below 1 diverges.
    x12 <- (tsoj + 1e-6) / scale12
    haz12 <- (shape12 / scale12) * x12^(shape12 - 1) / (1 + x12^shape12)
    sur1 <- 1 / (1 + (tsoj / scale12)^shape12)
    x21 <- (tsoj + 1e-6) / scale21
    haz21 <- (shape21 / scale21) * x21^(shape21 - 1) / (1 + x21^shape21)
    sur2 <- 1 / (1 + (tsoj / scale21)^shape21)
    xdrop <- (t + 1e-6) / scale_drop
    haz_drop <- (shape_drop / scale_drop) * xdrop^(shape_drop - 1) / (1 + xdrop^shape_drop)
    sur_drop <- 1 / (1 + (t / scale_drop)^shape_drop)
  })
}
