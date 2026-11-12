Du_2025_repotrectinib_orr_ros1 <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-efficacy model for the probability",
    "of a confirmed objective response (complete or partial response) to",
    "repotrectinib, a ROS1 / TRK / ALK tyrosine kinase inhibitor, in adults",
    "with ROS1-positive non-small cell lung cancer from the TRIDENT-1 trial",
    "(Du 2025, n = 127, 124 evaluable for response). The logit of the",
    "response probability is linear in the natural log of Cavg56, the",
    "average repotrectinib plasma concentration over the first 56 days of",
    "treatment computed from the actual dosing record, plus age, body",
    "weight, baseline tumor size, sex, Asian race, ECOG performance status",
    "and prior ROS1 TKI treatment. Odds ratio 4.71 per e-fold increase in",
    "Cavg56; TKI-pretreated patients have 85% lower odds of response than",
    "TKI-naive patients. There is no PK layer and no ODE: Cavg56 is supplied",
    "as a data column, derived in the source analysis from the individual",
    "post hoc parameters of the companion population PK model, packaged as",
    "Du_2025_repotrectinib. No random effect and no residual error are",
    "estimated (Bernoulli likelihood). The intercept is not printed in the",
    "paper; it was back-solved from the observed response counts (see",
    "logit_ref and the vignette). One of nine exposure-response models in",
    "the Du_2025_repotrectinib_* family."
  )
  reference <- paste(
    "Du S, Hu Z, Shen J, Zhu L, Roy A, Lam J, Lu M, Kondic A, Hamuro L.",
    "Exposure-Response Analysis of Repotrectinib to Support the Dose",
    "Recommendation for Patients With ROS1-Positive NSCLC or NTRK-Positive",
    "Solid Tumors. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):1993-2005.",
    "doi:10.1002/psp4.70102. Individual repotrectinib exposures derive from",
    "the companion population PK model (doi:10.1002/psp4.70036); see",
    "modellib('Du_2025_repotrectinib')."
  )
  vignette <- "Du_2025_repotrectinib_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-efficacy regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr (probability of a confirmed objective response, 0-1; also logit_orr)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cavg56: average repotrectinib plasma concentration over the first 56 days of treatment, per subject",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Du 2025 Methods, Analysis Data and Table S3: 'Simulated average",
        "concentration over the first 56 days using actual dosing records',",
        "computed from the empirical Bayes estimates of the companion popPK",
        "model (Du_2025_repotrectinib), so it reflects the day-15 switch from",
        "160 mg once daily to twice daily and any dose interruption or",
        "reduction in the first two cycles. Enters on the NATURAL log scale,",
        "uncentred (Table S5 formula ORR ~ LCAVG56 + ...). Figure 1A pins the",
        "log base: (523/273)^1.551 = 2.74 and (136/273)^1.551 = 0.339 against",
        "the printed P95 and P5 odds ratios 2.74 and 0.338. Analysis-set",
        "median 273 ng/mL, P5-P95 136-523 ng/mL (Figure 1A). Units are",
        "load-bearing: supplying ug/mL would shift the logit by",
        "-1.551 * log(1000) = -10.7."
      ),
      source_name = "LCAVG56 (log Cavg56)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear term on the logit (Table S5 column AGE2). ROS1 analysis set median 57 years, range 28-80 (Table S1); Figure 1A reference 57 (P5-P95 35-75).",
      source_name = "AGE2"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear term on the logit (Table S5 column WTB). ROS1 analysis set median 63.5 kg, range 39.5-119.3 (Table S1); Figure 1A reference 64 (P5-P95 45.1-89).",
      source_name = "WTB"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Du 2025 reports baseline tumor size in cm (Table S5 column TSIZEB,",
        "Table S1 'Baseline Tumor Burden [cm]'); the register's canonical unit",
        "is mm, so model() divides by 10 and the coefficient keeps its",
        "printed per-cm value. The paper does not define the measurement",
        "further; a RECIST sum of target-lesion diameters is the usual",
        "construct. ROS1 analysis set median 5.0 cm, range 1.0-24.7",
        "(Table S1); Figure 1A reference 4.9 cm (P5-P95 1.44-14.7)."
      ),
      source_name = "TSIZEB"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table S6 'Sex [Female: Male]': the female level carries the coefficient, same orientation as SEXF. 81 female of 127 (63.8%; Table S1, Figure 1A).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White or other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or other race)",
      notes = "Table S6 'Race [Asian:Caucasian/Others]'. 75 Asian of 127 (59.1%; Table S1, Figure 1A).",
      source_name = "RACE2"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes = "Table S6 'Baseline Performance Status [1:0]'; Table S1 lists only PS 0 and PS 1, so the indicator equals PS = 1. 85 of 127 with PS 1 (66.9%; Table S1, Figure 1A).",
      source_name = "PS"
    ),
    PRIOR_TKI = list(
      description = "Prior treatment with a tyrosine kinase inhibitor targeting the tumor's driver fusion (here a ROS1 TKI); 1 = TKI-pretreated, 0 = TKI-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (TKI-naive)",
      notes = paste(
        "Table S6 'TKI pretreatment [TKI-pretreated:TKI-naive]'. 56 of 127",
        "pretreated (44.1%; Figure 1A, Table S1 PFS column). In the ROS1",
        "cohort the pretreated group had one prior ROS1 TKI and no",
        "chemotherapy or immunotherapy (Table S1 population-type rows)."
      ),
      source_name = "PRVTRTC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 127L,
    n_studies = 1L,
    n_observations = "124 binary confirmed-response records (one per evaluable patient; landmark analysis)",
    age_range = "28-80 years (median 57)",
    weight_range = "39.5-119.3 kg (median 63.5)",
    sex_female_pct = 63.8,
    race_ethnicity = c(Asian = 59.1, `White or other` = 40.9),
    disease_state = "Locally advanced or metastatic ROS1-positive non-small cell lung cancer; 71 ROS1 TKI-naive and 56 pretreated with one prior ROS1 TKI and no chemotherapy or immunotherapy",
    dose_range = "Repotrectinib orally; 114 of 124 response-evaluable patients on 160 mg once daily for 14 days then 160 mg twice daily (the recommended phase 2 dose), the rest on phase 1 regimens of 40-240 mg once daily or 160 mg twice daily",
    regions = "Multinational (TRIDENT-1, NCT03093116)",
    notes = paste(
      "Du 2025 Table S1 (ER-efficacy analysis dataset, data cutoff March",
      "2023) and Figure 1A. Table S1's ORR 'Responders' and 'Non-responders'",
      "column headers are transposed: as printed they give TKI-naive ORR",
      "13/69 = 19% and pretreated 34/55 = 62%, contradicting the fitted odds",
      "ratio of 0.147 for pretreatment; with the headers swapped the ROS1 set",
      "has 77 responders and 47 non-responders (naive 56/69 = 81%, pretreated",
      "21/55 = 38%). See the vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Du 2025 Table S6, 'Parameter Estimates of the E-R of OR (Full
    # Models)', ROS1-positive NSCLC block. The full model was fitted with
    # rms::lrm (Methods; Supplement S1 code: lrm(ORR ~ LCAVG56 + AGE2 +
    # WTB + TSIZEB + SEX + RACE2 + PS + PRVTRTC)), so every covariate is
    # UNCENTRED and the intercept is the logit at all-zero covariates.
    #
    # THE INTERCEPT IS NOT PRINTED anywhere in the paper or supplement
    # (Table S6 lists only the slopes). It is back-solved from the
    # logistic maximum-likelihood score equation for the intercept and
    # for the TKI-pretreatment dummy, which force the fitted response
    # probabilities to sum to the observed responder count overall and
    # within each TKI stratum:
    #   - Observed responders (Table S1, with its transposed
    #     Responders / Non-responders headers corrected): 56/69 TKI-naive,
    #     21/55 pretreated, 77/124 overall.
    #   - Covariate distribution: Table S1 totals (age 55.6 +/- 12.4 y,
    #     WT 65.4 +/- 14.1 kg, tumor size 6.1 +/- 4.3 cm, 64.5% female,
    #     58.1% Asian, 66.9% PS 1) and Figure 1A Cavg56 (median 273,
    #     P5-P95 136-523 ng/mL, log-normal), drawn independently
    #     (200,000 draws).
    #   - Solving mean(expit(logit)) = observed proportion gives -5.42
    #     overall, -5.30 within TKI-naive and -5.53 within pretreated.
    #     The overall value is used.
    #   - INDEPENDENT CHECK: driving this model with Cavg56 simulated from
    #     Du_2025_repotrectinib reproduces the paper's model-predicted
    #     median Pr(OR) at every simulated regimen (Table S7), to a median
    #     logit difference of about 0.01 across the 12 regimen-by-stratum
    #     cells; see the vignette.
    # ==================================================================
    logit_ref <- -5.42 ; label("Logit of the probability of a confirmed objective response at all-zero covariates (CAV = 1 ng/mL; unitless logit)")  # back-solved, not printed: logistic score equation on Du 2025 Table S1 response counts (77/124) over the Table S1 / Figure 1A covariate distribution; see block comment

    e_cav_logit        <- 1.551    ; label("Log-odds of response per natural-log unit of Cavg56 (unitless)")                    # Du 2025 Table S6 ROS1 'Repo Log Cavg56 [ng/mL]' 1.551, SE 0.5871; OR 4.714 (1.492, 14.9)
    e_age_logit        <- -0.0113  ; label("Log-odds of response per year of age (1/year)")                                     # Du 2025 Table S6 ROS1 'Age [yr]' -0.0113, SE 0.01733; OR 0.9888
    e_wt_logit         <- -0.0118  ; label("Log-odds of response per kg of baseline body weight (1/kg)")                        # Du 2025 Table S6 ROS1 'Body Weight [kg]' -0.0118, SE 0.01852; OR 0.9883
    e_tumsz_logit      <- -0.01916 ; label("Log-odds of response per cm of baseline tumor size (1/cm)")                         # Du 2025 Table S6 ROS1 'Baseline Tumor Size (cm)' -0.01916, SE 0.05233; OR 0.981
    e_sexf_logit       <- 0.1005   ; label("Log-odds of response, female vs male (unitless)")                                   # Du 2025 Table S6 ROS1 'Sex [Female: Male]' 0.1005, SE 0.5009; OR 1.106
    e_race_asian_logit <- -0.4202  ; label("Log-odds of response, Asian vs Caucasian or other (unitless)")                      # Du 2025 Table S6 ROS1 'Race [Asian:Caucasian/Others]' -0.4202, SE 0.5204; OR 0.6569
    e_ecog_ge1_logit   <- -0.138   ; label("Log-odds of response, ECOG performance status 1 vs 0 (unitless)")                   # Du 2025 Table S6 ROS1 'Baseline Performance Status [1:0]' -0.138, SE 0.4723; OR 0.8711
    e_prior_tki_logit  <- -1.92    ; label("Log-odds of response, ROS1-TKI-pretreated vs TKI-naive (unitless)")                 # Du 2025 Table S6 ROS1 'TKI pretreatment [TKI-pretreated:TKI-naive]' -1.92, SE 0.4494; OR 0.1467 (0.06079, 0.3539)

    # The source likelihood is Bernoulli: no omega or sigma is estimated.
    # rxode2 needs an observation declaration, so the deterministic
    # probability carries a tiny placeholder additive residual (the
    # Liu_2024_saf189s_orr.R pattern).
    addSd_prob_orr <- fixed(0.001) ; label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Uncentred linear predictor, as fitted. Tumor size enters in cm, so
    # the canonical mm column is divided by 10.
    logit_orr <- logit_ref +
      e_cav_logit * log(CAV) +
      e_age_logit * AGE +
      e_wt_logit * WT +
      e_tumsz_logit * (TUMSZ / 10) +
      e_sexf_logit * SEXF +
      e_race_asian_logit * RACE_ASIAN +
      e_ecog_ge1_logit * ECOG_GE1 +
      e_prior_tki_logit * PRIOR_TKI

    prob_orr <- expit(logit_orr)

    # Deterministic response probability; sample outcomes downstream with
    # rbinom(n, 1, prob_orr).
    prob_orr ~ add(addSd_prob_orr)
  })
}
