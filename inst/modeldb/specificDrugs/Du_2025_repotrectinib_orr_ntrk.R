Du_2025_repotrectinib_orr_ntrk <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-efficacy model for the probability",
    "of a confirmed objective response (complete or partial response) to",
    "repotrectinib, a ROS1 / TRK / ALK tyrosine kinase inhibitor, in adults",
    "with NTRK-fusion-positive advanced solid tumors from the TRIDENT-1",
    "trial (Du 2025, n = 88, 83 evaluable for response). The logit of the",
    "response probability is linear in the natural log of Cavg56, the",
    "average repotrectinib plasma concentration over the first 56 days of",
    "treatment, plus age, body weight, baseline tumor size, sex, Asian race,",
    "ECOG performance status, prior TRK TKI treatment, NTRK3 (vs NTRK1 or",
    "NTRK2) fusion and non-NSCLC tumor type. Exposure is NOT a significant",
    "predictor here (odds ratio 1.33 per e-fold, 95% CI 0.53-3.37); the",
    "dominant effect is NTRK3 fusion (odds ratio 9.57). There is no PK layer",
    "and no ODE: Cavg56 is supplied as a data column, derived in the source",
    "analysis from the companion population PK model, packaged as",
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
        "Du 2025 Methods, Analysis Data and Table S3: simulated average",
        "concentration over the first 56 days using actual dosing records,",
        "from the empirical Bayes estimates of the companion popPK model",
        "(Du_2025_repotrectinib). Enters on the NATURAL log scale, uncentred",
        "(Table S5 formula ORR ~ LCAVG56 + ...). Figure 1B pins the log base:",
        "(478/257)^0.2871 = 1.195 and (79/257)^0.2871 = 0.713 against the",
        "printed P95 and P5 odds ratios 1.19 and 0.713. Analysis-set median",
        "257 ng/mL, P5-P95 79-478 ng/mL (Figure 1B)."
      ),
      source_name = "LCAVG56 (log Cavg56)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear term on the logit (Table S5 column AGE2). NTRK analysis set median 60 years, range 20-84 (Table S1); Figure 1B reference 60 (P5-P95 27-79).",
      source_name = "AGE2"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear term on the logit (Table S5 column WTB). NTRK analysis set median 70 kg, range 42-116.3 (Table S1); Figure 1B reference 70 (P5-P95 48.3-96.1).",
      source_name = "WTB"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Du 2025 reports baseline tumor size in cm (TSIZEB); model() divides the canonical mm column by 10 so the coefficient keeps its printed per-cm value. NTRK analysis set median 5.7 cm, range 1.2-26.6 (Table S1); Figure 1B reference 5.7 cm (P5-P95 1.51-17).",
      source_name = "TSIZEB"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table S6 'Sex [Female: Male]'. 47 female of 88 (53.4%; Table S1, Figure 1B).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White or other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or other race)",
      notes = "Table S6 'Race [Asian:Caucasian/Others]'. 33 Asian of 88 (37.5%; Table S1, Figure 1B).",
      source_name = "RACE2"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes = "Table S6 'Baseline Performance Status [1:0]'. 51 of 88 with PS 1 (58.0%; Table S1, Figure 1B).",
      source_name = "PS"
    ),
    PRIOR_TKI = list(
      description = "Prior treatment with a tyrosine kinase inhibitor targeting the tumor's driver fusion (here a TRK TKI); 1 = TKI-pretreated, 0 = TKI-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (TKI-naive)",
      notes = "Table S6 'TKI pretreatment [TKI-pretreated:TKI-naive]'. 48 of 88 pretreated (54.5%; Table S1, Figure 1B).",
      source_name = "PRVTRTC"
    ),
    TUM_NTRK3_MUT = list(
      description = "NTRK3 gene fusion indicator within NTRK-fusion-positive tumors; 1 = NTRK3 fusion, 0 = NTRK1 or NTRK2 fusion",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NTRK1 or NTRK2 fusion)",
      notes = "Table S6 'Genetic Mutation (NTRK3: NTRK1-2)' (source column GENMTC). 51 NTRK3 of 88 (58%; Figure 1B); by stratum 18 of 40 TKI-naive and 33 of 48 pretreated (Table S11 column counts).",
      source_name = "GENMTC"
    ),
    TUMTP_OTHER = list(
      description = "Non-NSCLC tumor type indicator; 1 = any solid tumor other than non-small cell lung cancer, 0 = NSCLC",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-small cell lung cancer, the only named tumor type in this model)",
      notes = "Table S6 'Tumor Type (Others:NSCLC)' (source column TTYPEF). 53 non-NSCLC of 88 (60%; Figure 1B).",
      source_name = "TTYPEF"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 1L,
    n_observations = "83 binary confirmed-response records (one per evaluable patient; landmark analysis)",
    age_range = "20-84 years (median 60)",
    weight_range = "42-116.3 kg (median 70)",
    sex_female_pct = 53.4,
    race_ethnicity = c(Asian = 37.5, `White or other` = 62.5),
    disease_state = "Locally advanced or metastatic NTRK1-3-fusion-positive solid tumors (53 of 88 non-NSCLC); 40 TRK TKI-naive and 48 TRK TKI-pretreated; 51 NTRK3 and 37 NTRK1 or NTRK2 fusions",
    dose_range = "Repotrectinib orally; 76 of 83 response-evaluable patients on 160 mg once daily for 14 days then 160 mg twice daily, the rest on 40 or 160 mg once daily or 160 mg twice daily",
    regions = "Multinational (TRIDENT-1, NCT03093116)",
    notes = paste(
      "Du 2025 Table S1 (ER-efficacy analysis dataset, data cutoff March",
      "2023) and Figure 1B. As for the ROS1 cohort, Table S1's ORR",
      "'Responders' and 'Non-responders' headers are transposed, and the",
      "Results sentence quoting '36 responders and 47 non-responders' repeats",
      "the transposition: the corrected counts are 47 responders and 36",
      "non-responders (naive 23/37 = 62%, pretreated 24/46 = 52%), the only",
      "reading consistent with the fitted odds ratio of 0.349 for",
      "pretreatment and with the paper's own Table S7 predictions. See the",
      "vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Du 2025 Table S6, 'Parameter Estimates of the E-R of OR (Full
    # Models)', NTRK-positive solid tumors block. rms::lrm full model
    # (Table S5: ORR ~ LCAVG56 + AGE2 + WTB + TSIZEB + SEX + RACE2 + PS +
    # PRVTRTC + GENMTC + TTYPEF), every covariate UNCENTRED.
    #
    # THE INTERCEPT IS NOT PRINTED. It is back-solved from the logistic
    # score equations, which force the fitted probabilities to sum to the
    # observed responder count overall and within each TKI stratum:
    #   - Observed responders (Table S1 with its transposed headers
    #     corrected): 23/37 TKI-naive, 24/46 pretreated, 47/83 overall.
    #   - Covariate distribution: Table S1 totals (age 55.6 +/- 15.6 y,
    #     WT 71.2 +/- 16.3 kg, tumor size 7.2 +/- 5.5 cm, 55.4% female,
    #     36.1% Asian, 56.6% PS 1), Figure 1B (Cavg56 median 257,
    #     P5-P95 79-478 ng/mL, log-normal; 53/88 non-NSCLC) and NTRK3
    #     prevalence by TKI stratum from the Table S11 column counts
    #     (18/40 naive, 33/48 pretreated); 200,000 independent draws.
    #   - Solving mean(expit(logit)) = observed proportion gives -0.56
    #     overall, -0.49 within TKI-naive and -0.62 within pretreated.
    #     The overall value is used.
    #   - INDEPENDENT CHECK: driving this model with Cavg56 simulated from
    #     Du_2025_repotrectinib over a Table S1 virtual cohort reproduces
    #     the paper's model-predicted medians (Table S7) to a median logit
    #     difference of about 0.05 across the 24 regimen-by-subgroup cells;
    #     see the vignette.
    # ==================================================================
    logit_ref <- -0.56 ; label("Logit of the probability of a confirmed objective response at all-zero covariates (CAV = 1 ng/mL; unitless logit)")  # back-solved, not printed: logistic score equation on Du 2025 Table S1 response counts (47/83) over the Table S1 / Figure 1B covariate distribution; see block comment

    e_cav_logit         <- 0.2871   ; label("Log-odds of response per natural-log unit of Cavg56 (unitless)")                   # Du 2025 Table S6 NTRK 'Repo Log Cavg56 [ng/mL]' 0.2871, SE 0.4737; OR 1.333 (0.5266, 3.372)
    e_age_logit         <- -0.03726 ; label("Log-odds of response per year of age (1/year)")                                    # Du 2025 Table S6 NTRK 'Age [yr]' -0.03726, SE 0.0197; OR 0.9634
    e_wt_logit          <- 0.01421  ; label("Log-odds of response per kg of baseline body weight (1/kg)")                       # Du 2025 Table S6 NTRK 'Body Weight [kg]' 0.01421, SE 0.01964; OR 1.014
    e_tumsz_logit       <- -0.07697 ; label("Log-odds of response per cm of baseline tumor size (1/cm)")                        # Du 2025 Table S6 NTRK 'Baseline Tumor Size (cm)' -0.07697, SE 0.04836; OR 0.9259
    e_sexf_logit        <- 0.3882   ; label("Log-odds of response, female vs male (unitless)")                                  # Du 2025 Table S6 NTRK 'Sex [Female: Male]' 0.3882, SE 0.5972; OR 1.474
    e_race_asian_logit  <- 0.7497   ; label("Log-odds of response, Asian vs Caucasian or other (unitless)")                     # Du 2025 Table S6 NTRK 'Race [Asian:Caucasian/Others]' 0.7497, SE 0.6621; OR 2.116
    e_ecog_ge1_logit    <- -0.1154  ; label("Log-odds of response, ECOG performance status 1 vs 0 (unitless)")                  # Du 2025 Table S6 NTRK 'Baseline Performance Status [1:0]' -0.1154, SE 0.5662; OR 0.891
    e_prior_tki_logit   <- -1.054   ; label("Log-odds of response, TRK-TKI-pretreated vs TKI-naive (unitless)")                 # Du 2025 Table S6 NTRK 'TKI pretreatment [TKI-pretreated:TKI-naive]' -1.054, SE 0.6384; OR 0.3486
    e_tum_ntrk3_mut_logit <- 2.259  ; label("Log-odds of response, NTRK3 vs NTRK1 or NTRK2 fusion (unitless)")                  # Du 2025 Table S6 NTRK 'Genetic Mutation (NTRK3: NTRK1-2)' 2.259, SE 0.6189; OR 9.575 (2.847, 32.2)
    e_tumtp_other_logit <- -0.3631  ; label("Log-odds of response, non-NSCLC vs NSCLC tumor type (unitless)")                   # Du 2025 Table S6 NTRK 'Tumor Type (Others:NSCLC)' -0.3631, SE 0.6001; OR 0.6955

    # Bernoulli likelihood in the source: placeholder additive residual
    # only so rxode2 accepts the observation declaration.
    addSd_prob_orr <- fixed(0.001) ; label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- logit_ref +
      e_cav_logit * log(CAV) +
      e_age_logit * AGE +
      e_wt_logit * WT +
      e_tumsz_logit * (TUMSZ / 10) +
      e_sexf_logit * SEXF +
      e_race_asian_logit * RACE_ASIAN +
      e_ecog_ge1_logit * ECOG_GE1 +
      e_prior_tki_logit * PRIOR_TKI +
      e_tum_ntrk3_mut_logit * TUM_NTRK3_MUT +
      e_tumtp_other_logit * TUMTP_OTHER

    prob_orr <- expit(logit_orr)

    prob_orr ~ add(addSd_prob_orr)
  })
}
