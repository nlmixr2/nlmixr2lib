Suzuki_2020_tofacitinib_acr20 <- function() {
  description <- paste(
    "Landmark logit-domain Emax dose-response model for the probability of",
    "an ACR20 response at Week 12 in adults with active rheumatoid arthritis",
    "and an inadequate response to DMARDs given placebo or tofacitinib",
    "monotherapy 1, 3, 5, 10 or 15 mg twice daily, pooled from a mostly",
    "Western (NCT00550446) and a Japanese (NCT00687193) phase 2 dose-ranging",
    "study (Suzuki 2020 final model). logit(p) = E0 + Emax * D / (D + ED50)",
    "with D the tofacitinib dose per administration (mg, DOSE_TOFACITINIB_MG;",
    "0 for placebo) and a multiplicative study effect on Emax for the",
    "Japanese study (STUDY_NCT00687193). Fitted separately from the ACR50 and",
    "ACR70 models. No PK layer, no ODE and no random effects."
  )
  reference <- paste(
    "Suzuki M, Shoji S, Miyoshi S, Krishnaswami S. Model-Based Comparison of",
    "Dose-Response Profiles of Tofacitinib in Japanese Versus Western",
    "Rheumatoid Arthritis Patients. J Clin Pharmacol. 2020;60(2):198-208.",
    "doi:10.1002/jcph.1514. PMID: 31512746. PMCID: PMC6972659."
  )
  vignette <- "Suzuki_2020_tofacitinib"
  units <- list(
    time = "n/a (static landmark model for the Week 12 assessment; no time dimension)",
    dosing = "mg (tofacitinib dose per twice-daily administration, supplied in the DOSE_TOFACITINIB_MG covariate; no dose events)",
    concentration = "prob_acr20 (probability of an ACR20 response at Week 12, 0-1)"
  )

  covariateData <- list(
    DOSE_TOFACITINIB_MG = list(
      description = "Randomized tofacitinib dose per twice-daily administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Suzuki 2020 Methods 'Base Model Description': Di is the randomized",
        "dose (in milligrams) for each administration given twice daily.",
        "Studied levels 0 (placebo), 1, 3, 5, 10 and 15 mg twice daily.",
        "Time-fixed per patient."
      ),
      source_name = "Di"
    ),
    STUDY_NCT00687193 = list(
      description = "Japanese phase 2 dose-ranging study (NCT00687193) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = the mostly Western phase 2 study NCT00550446",
      notes = paste(
        "Suzuki 2020 Methods 'Inclusion of Covariates': x = 1 for the effect",
        "of the Japanese study compared with the Western study. All patients",
        "in NCT00687193 were Japanese and none in NCT00550446 were, so the",
        "study indicator is also the paper's Japanese/Western designation;",
        "the paper notes study and population differences are confounded."
      ),
      source_name = "study"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = non-Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-Asian",
      notes = paste(
        "Tested on Emax, ED50 and E0 (Supplementary Table S1; ACR20 dMOF",
        "-9.6, -8.8 and -3.7) as a surrogate for the study effect; the final",
        "model retains the study effect instead. No estimate is printed."
      ),
      source_name = "race (Asian or non-Asian)"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on Emax, ED50 and E0 (Supplementary Table S2; ACR20 dMOF",
        "-4.4, -5.2 and -3.5), not significant at alpha = 0.01. No estimate",
        "is printed."
      ),
      source_name = "body weight"
    ),
    SCORE_DAS28ESR = list(
      description = "Baseline 28-joint Disease Activity Score with erythrocyte sedimentation rate (DAS28-4(ESR))",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on Emax, ED50 and E0 (Supplementary Table S3; ACR20 dMOF",
        "-1.0, -0.48 and -2.5), not retained. No estimate is printed."
      ),
      source_name = "DAS28-4(ESR)"
    ),
    BLHAQ = list(
      description = "Baseline HAQ-DI score",
      units = "unitless (0-3 composite score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on Emax, ED50 and E0 (Supplementary Table S4; ACR20 dMOF",
        "-0.57, -0.64 and -0.26), not retained. No estimate is printed."
      ),
      source_name = "HAQ-DI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 648L,
    n_studies = 2L,
    age_range = "18-83 years (mean 53.4)",
    weight_range = "31.4-125 kg (mean 62.5; Japanese study mean 54.4, Western study mean 70.3)",
    sex_female_pct = 85.2,
    race_ethnicity = c(White = 38.7, Black = 1.2, Asian = 53.7, Other = 6.3),
    disease_state = paste(
      "Active rheumatoid arthritis with an inadequate response to at least",
      "one DMARD (lack of efficacy or toxicity); tofacitinib monotherapy"
    ),
    dose_range = "Placebo or tofacitinib 1, 3, 5, 10 or 15 mg orally twice daily; ACR response at Week 12 (adalimumab arm of NCT00550446 excluded)",
    regions = "Japan (NCT00687193, n = 317) and mostly Western countries (NCT00550446, n = 331)",
    notes = paste(
      "Suzuki 2020 Table 1 and Study Design. Missing ACR components imputed",
      "by last observation carried forward. Fitted by maximum likelihood in",
      "NONMEM version 6; the Bernoulli likelihood estimates no random effects."
    )
  )

  ini({
    # logit P(ACR20 = 1) = E0 + Emax * D / (D + ED50) (Methods 'Base Model
    # Description'); Emax = Emax0 * theta_study^STUDY (Methods 'Inclusion of
    # Covariates'). Final model = Table 2 ACR20 run 2.
    e0 <- -1.3; label("Logit of the ACR20 response probability on placebo, E0 (logit)") # Table 4, ACR20 'E0' = -1.3 (RSE -16.7)
    emax <- 3.2; label("Maximum tofacitinib effect on the ACR20 logit in the Western study, Emax (logit)") # Table 4, ACR20 'Emax' = 3.2 (RSE 13.7)
    led50 <- log(4.2); label("Log dose per administration giving half the maximum effect, ED50 (log mg)") # Table 4, ACR20 'ED50' = 4.2 mg (RSE 39.2); log ED50 was estimated (Methods)
    e_study_nct00687193_emax <- 1.4; label("Ratio of Emax in the Japanese study to Emax in the Western study (ratio)") # Table 4, ACR20 'Study effect on Emax' = 1.4 (RSE 10.5)

    # Bernoulli likelihood in the source; no residual error was estimated.
    addSd_prob_acr20 <- fixed(0.001); label("Placeholder additive residual SD on the response probability (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    ed50 <- exp(led50)
    emax_i <- emax * e_study_nct00687193_emax^STUDY_NCT00687193

    logit_acr20 <- e0 + emax_i * DOSE_TOFACITINIB_MG / (DOSE_TOFACITINIB_MG + ed50)
    prob_acr20 <- expit(logit_acr20)

    prob_acr20 ~ add(addSd_prob_acr20)
  })
}
