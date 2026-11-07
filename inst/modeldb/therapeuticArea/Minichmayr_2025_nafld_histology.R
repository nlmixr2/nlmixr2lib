Minichmayr_2025_nafld_histology <- function() {
  description <- paste(
    "Multidimensional item response theory (IRT) model of the 13 NASH CRN",
    "liver-biopsy histology scores in adults with nonalcoholic fatty liver",
    "disease (Minichmayr 2025). Five correlated standard-normal latent",
    "variables represent the disease facets: LV1 steatosis, LV2 lobular",
    "inflammation, LV3 hepatocellular ballooning and LV4 fibrosis each drive",
    "their own cardinal ordinal score alone, and LV5 drives only the nine",
    "residual features. The five ordinal scores (steatosis 0-3, lobular",
    "inflammation 0-3, ballooning 0-2, fibrosis 0-4, portal chronic",
    "inflammation 0-2) are graded-response items and the eight present/absent",
    "features are two-parameter logistic items; the nine residual items load on",
    "all five latent variables. Cross-sectional (one biopsy per patient): no",
    "time dependence, no drug and no covariates. The latent-variable values",
    "are the random effects etalv1-etalv5."
  )
  reference <- paste(
    "Minichmayr IK, Plan EL, Weber B, Ueckert S.",
    "A Model-Based Evaluation of Noninvasive Biomarkers to Reflect Histological",
    "Nonalcoholic Fatty Liver Disease Scores.",
    "Pharm Res. 2025;42(1):123.",
    "doi:10.1007/s11095-024-03791-2."
  )
  vignette <- "Minichmayr_2025_nafld_histology"

  # The latent disease facets are random effects with no parameter to scale:
  # each eta IS the subject's latent-variable value (mean 0, variance 1).
  paper_specific_etas <- c("etalv1", "etalv2", "etalv3", "etalv4", "etalv5")

  units <- list(
    time = "(none; cross-sectional model of a single liver biopsy per patient, with no time dependence)",
    dosing = "(none; disease model with no drug input)",
    concentration = "(ordinal and binary NASH CRN histology scores; the observations are biopsy grades, and the outputs p<k>_<item> are the probabilities of each grade, not a drug concentration)"
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 914L,
    n_studies = 1L,
    age_range = "median 51.3 years; 2.5th-97.5th percentile 23.5-69.0 years (Minichmayr 2025 Table I)",
    weight_range = "median 95.0 kg; 2.5th-97.5th percentile 63-139 kg (Minichmayr 2025 Table I)",
    sex_female_pct = 64.2,
    race_ethnicity = "10.1% Hispanic/Latino; race not reported (Minichmayr 2025 Table I)",
    disease_state = paste(
      "Adults with histologically proven or suspected nonalcoholic fatty liver",
      "disease (NAFLD); the cohort spans fibrosis stages 0-4 and NAFLD activity",
      "scores 0-8. 52.3% had a definite and 26.8% a definite negative diagnosis",
      "of NASH; 38.1% had diabetes; median BMI 33.7 kg/m^2."
    ),
    dose_range = "(none; untreated natural-history cohort)",
    regions = "United States (nine NASH CRN medical centres)",
    notes = paste(
      "Observational NAFLD Database study of the NASH Clinical Research Network",
      "(NIDDK NAFLD Adult Database), enrolled 2004-2009; one biopsy per patient",
      "(the one closest to the liver-panel laboratory tests), read centrally by",
      "the NASH CRN Pathology Committee. Fibrosis sub-stages 1a/1b/1c were pooled",
      "as stage 1. Fitted in R 3.6.1 with the mirt package (MHRM estimation).",
      "The paper also relates 69 noninvasive biomarkers to the latent variables",
      "with a full random effects model (NONMEM 7.4); that biomarker-latent",
      "correlation matrix is not reported numerically and is not part of this",
      "model (see the vignette)."
    )
  )

  ini({
    # Minichmayr 2025 Online Resource Table SI. mirt slope-intercept form:
    #   graded items  logit P(score >= k) = a1*LV1 + ... + a5*LV5 + d_k
    #   binary items  logit P(present)    = a1*LV1 + ... + a5*LV5 + d
    # Table SI prints a = 0 for the latent variables a cardinal score does not
    # load on (the model structure, not an estimate); those terms are omitted.

    # Steatosis grade 0-3 (loads on LV1 only)
    a1_steatosis <- 1.287
    label("Steatosis: discrimination on LV1 (logit per latent SD)") # Table SI, Steatosis a1 = 1.287
    d1_steatosis <- 3.026
    label("Steatosis: intercept of logit P(grade >= 1) (logit)") # Table SI, Steatosis d1 = 3.026
    d2_steatosis <- -0.051
    label("Steatosis: intercept of logit P(grade >= 2) (logit)") # Table SI, Steatosis d2 = -0.051
    d3_steatosis <- -1.887
    label("Steatosis: intercept of logit P(grade >= 3) (logit)") # Table SI, Steatosis d3 = -1.887

    # Lobular inflammation grade 0-3 (loads on LV2 only)
    a2_inflammation <- 2.707
    label("Lobular inflammation: discrimination on LV2 (logit per latent SD)") # Table SI, Inflammation a2 = 2.707
    d1_inflammation <- 8.413
    label("Lobular inflammation: intercept of logit P(grade >= 1) (logit)") # Table SI, Inflammation d1 = 8.413
    d2_inflammation <- -0.863
    label("Lobular inflammation: intercept of logit P(grade >= 2) (logit)") # Table SI, Inflammation d2 = -0.863
    d3_inflammation <- -4.212
    label("Lobular inflammation: intercept of logit P(grade >= 3) (logit)") # Table SI, Inflammation d3 = -4.212

    # Hepatocellular ballooning grade 0-2 (loads on LV3 only)
    a3_ballooning <- 6.031
    label("Ballooning: discrimination on LV3 (logit per latent SD)") # Table SI, Ballooning a3 = 6.031
    d1_ballooning <- 1.721
    label("Ballooning: intercept of logit P(grade >= 1) (logit)") # Table SI, Ballooning d1 = 1.721
    d2_ballooning <- -2.057
    label("Ballooning: intercept of logit P(grade >= 2) (logit)") # Table SI, Ballooning d2 = -2.057

    # Fibrosis stage 0-4 (loads on LV4 only)
    a4_fibrosis <- 3.161
    label("Fibrosis: discrimination on LV4 (logit per latent SD)") # Table SI, Fibrosis a4 = 3.161
    d1_fibrosis <- 2.108
    label("Fibrosis: intercept of logit P(stage >= 1) (logit)") # Table SI, Fibrosis d1 = 2.108
    d2_fibrosis <- -0.090
    label("Fibrosis: intercept of logit P(stage >= 2) (logit)") # Table SI, Fibrosis d2 = -0.090
    d3_fibrosis <- -1.686
    label("Fibrosis: intercept of logit P(stage >= 3) (logit)") # Table SI, Fibrosis d3 = -1.686
    d4_fibrosis <- -4.046
    label("Fibrosis: intercept of logit P(stage >= 4) (logit)") # Table SI, Fibrosis d4 = -4.046

    # Portal chronic inflammation grade 0-2 (loads on all five latent variables)
    a1_portinfl <- -0.537
    label("Portal chronic inflammation: discrimination on LV1 (logit per latent SD)") # Table SI, Portal chronic inflammation a1 = -0.537
    a2_portinfl <- 0.668
    label("Portal chronic inflammation: discrimination on LV2 (logit per latent SD)") # Table SI, Portal chronic inflammation a2 = 0.668
    a3_portinfl <- -0.870
    label("Portal chronic inflammation: discrimination on LV3 (logit per latent SD)") # Table SI, Portal chronic inflammation a3 = -0.870
    a4_portinfl <- 2.104
    label("Portal chronic inflammation: discrimination on LV4 (logit per latent SD)") # Table SI, Portal chronic inflammation a4 = 2.104
    a5_portinfl <- 0.526
    label("Portal chronic inflammation: discrimination on LV5 (logit per latent SD)") # Table SI, Portal chronic inflammation a5 = 0.526
    d1_portinfl <- 3.039
    label("Portal chronic inflammation: intercept of logit P(grade >= 1) (logit)") # Table SI, Portal chronic inflammation d1 = 3.039
    d2_portinfl <- -1.817
    label("Portal chronic inflammation: intercept of logit P(grade >= 2) (logit)") # Table SI, Portal chronic inflammation d2 = -1.817

    # Microvesicular steatosis, absent/present
    a1_microsteat <- 3.839
    label("Microvesicular steatosis: discrimination on LV1 (logit per latent SD)") # Table SI, Microvesicular steatosis a1 = 3.839
    a2_microsteat <- -2.268
    label("Microvesicular steatosis: discrimination on LV2 (logit per latent SD)") # Table SI, Microvesicular steatosis a2 = -2.268
    a3_microsteat <- -0.411
    label("Microvesicular steatosis: discrimination on LV3 (logit per latent SD)") # Table SI, Microvesicular steatosis a3 = -0.411
    a4_microsteat <- 3.316
    label("Microvesicular steatosis: discrimination on LV4 (logit per latent SD)") # Table SI, Microvesicular steatosis a4 = 3.316
    a5_microsteat <- 0.435
    label("Microvesicular steatosis: discrimination on LV5 (logit per latent SD)") # Table SI, Microvesicular steatosis a5 = 0.435
    d_microsteat <- -5.031
    label("Microvesicular steatosis: intercept of logit P(present) (logit)") # Table SI, Microvesicular steatosis d = -5.031

    # Microgranulomas, absent/present
    a1_microgran <- 1.201
    label("Microgranulomas: discrimination on LV1 (logit per latent SD)") # Table SI, Microgranulomas a1 = 1.201
    a2_microgran <- 7.752
    label("Microgranulomas: discrimination on LV2 (logit per latent SD)") # Table SI, Microgranulomas a2 = 7.752
    a3_microgran <- -2.456
    label("Microgranulomas: discrimination on LV3 (logit per latent SD)") # Table SI, Microgranulomas a3 = -2.456
    a4_microgran <- -2.463
    label("Microgranulomas: discrimination on LV4 (logit per latent SD)") # Table SI, Microgranulomas a4 = -2.463
    a5_microgran <- 3.009
    label("Microgranulomas: discrimination on LV5 (logit per latent SD)") # Table SI, Microgranulomas a5 = 3.009
    d_microgran <- 5.922
    label("Microgranulomas: intercept of logit P(present) (logit)") # Table SI, Microgranulomas d = 5.922

    # Large lipogranulomas, absent/present
    a1_lipogran <- -0.004
    label("Large lipogranulomas: discrimination on LV1 (logit per latent SD)") # Table SI, Large lipogranulomas a1 = -0.004
    a2_lipogran <- 0.113
    label("Large lipogranulomas: discrimination on LV2 (logit per latent SD)") # Table SI, Large lipogranulomas a2 = 0.113
    a3_lipogran <- -0.346
    label("Large lipogranulomas: discrimination on LV3 (logit per latent SD)") # Table SI, Large lipogranulomas a3 = -0.346
    a4_lipogran <- 0.826
    label("Large lipogranulomas: discrimination on LV4 (logit per latent SD)") # Table SI, Large lipogranulomas a4 = 0.826
    a5_lipogran <- 0.698
    label("Large lipogranulomas: discrimination on LV5 (logit per latent SD)") # Table SI, Large lipogranulomas a5 = 0.698
    d_lipogran <- -0.411
    label("Large lipogranulomas: intercept of logit P(present) (logit)") # Table SI, Large lipogranulomas d = -0.411

    # Acidophil bodies, absent/present
    a1_acidophil <- -0.180
    label("Acidophil bodies: discrimination on LV1 (logit per latent SD)") # Table SI, Acidophil bodies a1 = -0.180
    a2_acidophil <- 1.156
    label("Acidophil bodies: discrimination on LV2 (logit per latent SD)") # Table SI, Acidophil bodies a2 = 1.156
    a3_acidophil <- 0.037
    label("Acidophil bodies: discrimination on LV3 (logit per latent SD)") # Table SI, Acidophil bodies a3 = 0.037
    a4_acidophil <- 0.051
    label("Acidophil bodies: discrimination on LV4 (logit per latent SD)") # Table SI, Acidophil bodies a4 = 0.051
    a5_acidophil <- 0.544
    label("Acidophil bodies: discrimination on LV5 (logit per latent SD)") # Table SI, Acidophil bodies a5 = 0.544
    d_acidophil <- -1.264
    label("Acidophil bodies: intercept of logit P(present) (logit)") # Table SI, Acidophil bodies d = -1.264

    # Pigmented macrophages, absent/present
    a1_pigmac <- 0.142
    label("Pigmented macrophages: discrimination on LV1 (logit per latent SD)") # Table SI, Pigmented macrophages a1 = 0.142
    a2_pigmac <- 3.076
    label("Pigmented macrophages: discrimination on LV2 (logit per latent SD)") # Table SI, Pigmented macrophages a2 = 3.076
    a3_pigmac <- -1.238
    label("Pigmented macrophages: discrimination on LV3 (logit per latent SD)") # Table SI, Pigmented macrophages a3 = -1.238
    a4_pigmac <- -1.098
    label("Pigmented macrophages: discrimination on LV4 (logit per latent SD)") # Table SI, Pigmented macrophages a4 = -1.098
    a5_pigmac <- 2.892
    label("Pigmented macrophages: discrimination on LV5 (logit per latent SD)") # Table SI, Pigmented macrophages a5 = 2.892
    d_pigmac <- 4.402
    label("Pigmented macrophages: intercept of logit P(present) (logit)") # Table SI, Pigmented macrophages d = 4.402

    # Megamitochondria, absent/present
    a1_megamito <- 1.258
    label("Megamitochondria: discrimination on LV1 (logit per latent SD)") # Table SI, Megamitochondria a1 = 1.258
    a2_megamito <- -1.063
    label("Megamitochondria: discrimination on LV2 (logit per latent SD)") # Table SI, Megamitochondria a2 = -1.063
    a3_megamito <- 0.305
    label("Megamitochondria: discrimination on LV3 (logit per latent SD)") # Table SI, Megamitochondria a3 = 0.305
    a4_megamito <- 1.274
    label("Megamitochondria: discrimination on LV4 (logit per latent SD)") # Table SI, Megamitochondria a4 = 1.274
    a5_megamito <- 0.297
    label("Megamitochondria: discrimination on LV5 (logit per latent SD)") # Table SI, Megamitochondria a5 = 0.297
    d_megamito <- -2.208
    label("Megamitochondria: intercept of logit P(present) (logit)") # Table SI, Megamitochondria d = -2.208

    # Mallory's hyaline, absent/present
    a1_mallory <- -0.881
    label("Mallory's hyaline: discrimination on LV1 (logit per latent SD)") # Table SI, Mallory's hyaline a1 = -0.881
    a2_mallory <- 1.379
    label("Mallory's hyaline: discrimination on LV2 (logit per latent SD)") # Table SI, Mallory's hyaline a2 = 1.379
    a3_mallory <- 5.754
    label("Mallory's hyaline: discrimination on LV3 (logit per latent SD)") # Table SI, Mallory's hyaline a3 = 5.754
    a4_mallory <- 0.621
    label("Mallory's hyaline: discrimination on LV4 (logit per latent SD)") # Table SI, Mallory's hyaline a4 = 0.621
    a5_mallory <- -0.403
    label("Mallory's hyaline: discrimination on LV5 (logit per latent SD)") # Table SI, Mallory's hyaline a5 = -0.403
    d_mallory <- -4.300
    label("Mallory's hyaline: intercept of logit P(present) (logit)") # Table SI, Mallory's hyaline d = -4.300

    # Glycogen nuclei, absent/present
    a1_glycogen <- 0.098
    label("Glycogen nuclei: discrimination on LV1 (logit per latent SD)") # Table SI, Glyogen nuclei a1 = 0.098
    a2_glycogen <- -0.201
    label("Glycogen nuclei: discrimination on LV2 (logit per latent SD)") # Table SI, Glyogen nuclei a2 = -0.201
    a3_glycogen <- 0.350
    label("Glycogen nuclei: discrimination on LV3 (logit per latent SD)") # Table SI, Glyogen nuclei a3 = 0.350
    a4_glycogen <- -0.321
    label("Glycogen nuclei: discrimination on LV4 (logit per latent SD)") # Table SI, Glyogen nuclei a4 = -0.321
    a5_glycogen <- 0.029
    label("Glycogen nuclei: discrimination on LV5 (logit per latent SD)") # Table SI, Glyogen nuclei a5 = 0.029
    d_glycogen <- 0.167
    label("Glycogen nuclei: intercept of logit P(present) (logit)") # Table SI, Glyogen nuclei d = 0.167

    # Latent variables LV1 steatosis, LV2 inflammation, LV3 ballooning,
    # LV4 fibrosis, LV5 residual items only. Variances are 1 by the IRT
    # identification constraint (Table SI footnote a: mean 0, variance 1),
    # so the off-diagonals are the correlations printed in Figure 3. The
    # whole block is left unfixed because rxode2 cannot fix the diagonal of
    # a block alone and the correlations were estimated.
    etalv1 + etalv2 + etalv3 + etalv4 + etalv5 ~ c(
      1,
      0.60, 1,
      0.22, 0.62, 1,
      -0.17, 0.36, 0.69, 1,
      0.09, 0.15, 0.23, 0.02, 1
    )
  })

  model({
    lv1 <- etalv1
    lv2 <- etalv2
    lv3 <- etalv3
    lv4 <- etalv4
    lv5 <- etalv5

    # Graded-response items: cumulative P(grade >= k), then per-grade
    # probabilities. The intercepts decrease with k, so the cumulative
    # probabilities are ordered without a constraint.
    lp_steatosis <- a1_steatosis * lv1
    pge1_steatosis <- expit(lp_steatosis + d1_steatosis)
    pge2_steatosis <- expit(lp_steatosis + d2_steatosis)
    pge3_steatosis <- expit(lp_steatosis + d3_steatosis)
    p0_steatosis <- 1 - pge1_steatosis
    p1_steatosis <- pge1_steatosis - pge2_steatosis
    p2_steatosis <- pge2_steatosis - pge3_steatosis
    p3_steatosis <- pge3_steatosis
    escore_steatosis <- pge1_steatosis + pge2_steatosis + pge3_steatosis

    lp_inflammation <- a2_inflammation * lv2
    pge1_inflammation <- expit(lp_inflammation + d1_inflammation)
    pge2_inflammation <- expit(lp_inflammation + d2_inflammation)
    pge3_inflammation <- expit(lp_inflammation + d3_inflammation)
    p0_inflammation <- 1 - pge1_inflammation
    p1_inflammation <- pge1_inflammation - pge2_inflammation
    p2_inflammation <- pge2_inflammation - pge3_inflammation
    p3_inflammation <- pge3_inflammation
    escore_inflammation <- pge1_inflammation + pge2_inflammation + pge3_inflammation

    lp_ballooning <- a3_ballooning * lv3
    pge1_ballooning <- expit(lp_ballooning + d1_ballooning)
    pge2_ballooning <- expit(lp_ballooning + d2_ballooning)
    p0_ballooning <- 1 - pge1_ballooning
    p1_ballooning <- pge1_ballooning - pge2_ballooning
    p2_ballooning <- pge2_ballooning
    escore_ballooning <- pge1_ballooning + pge2_ballooning

    lp_fibrosis <- a4_fibrosis * lv4
    pge1_fibrosis <- expit(lp_fibrosis + d1_fibrosis)
    pge2_fibrosis <- expit(lp_fibrosis + d2_fibrosis)
    pge3_fibrosis <- expit(lp_fibrosis + d3_fibrosis)
    pge4_fibrosis <- expit(lp_fibrosis + d4_fibrosis)
    p0_fibrosis <- 1 - pge1_fibrosis
    p1_fibrosis <- pge1_fibrosis - pge2_fibrosis
    p2_fibrosis <- pge2_fibrosis - pge3_fibrosis
    p3_fibrosis <- pge3_fibrosis - pge4_fibrosis
    p4_fibrosis <- pge4_fibrosis
    escore_fibrosis <- pge1_fibrosis + pge2_fibrosis + pge3_fibrosis + pge4_fibrosis

    lp_portinfl <- a1_portinfl * lv1 + a2_portinfl * lv2 + a3_portinfl * lv3 +
      a4_portinfl * lv4 + a5_portinfl * lv5
    pge1_portinfl <- expit(lp_portinfl + d1_portinfl)
    pge2_portinfl <- expit(lp_portinfl + d2_portinfl)
    p0_portinfl <- 1 - pge1_portinfl
    p1_portinfl <- pge1_portinfl - pge2_portinfl
    p2_portinfl <- pge2_portinfl
    escore_portinfl <- pge1_portinfl + pge2_portinfl

    # Expected NAFLD activity score (NAS = steatosis + lobular inflammation +
    # ballooning, 0-8).
    escore_nas <- escore_steatosis + escore_inflammation + escore_ballooning

    # Two-parameter logistic items: probability the feature is present.
    p1_microsteat <- expit(a1_microsteat * lv1 + a2_microsteat * lv2 + a3_microsteat * lv3 +
      a4_microsteat * lv4 + a5_microsteat * lv5 + d_microsteat)
    p1_microgran <- expit(a1_microgran * lv1 + a2_microgran * lv2 + a3_microgran * lv3 +
      a4_microgran * lv4 + a5_microgran * lv5 + d_microgran)
    p1_lipogran <- expit(a1_lipogran * lv1 + a2_lipogran * lv2 + a3_lipogran * lv3 +
      a4_lipogran * lv4 + a5_lipogran * lv5 + d_lipogran)
    p1_acidophil <- expit(a1_acidophil * lv1 + a2_acidophil * lv2 + a3_acidophil * lv3 +
      a4_acidophil * lv4 + a5_acidophil * lv5 + d_acidophil)
    p1_pigmac <- expit(a1_pigmac * lv1 + a2_pigmac * lv2 + a3_pigmac * lv3 +
      a4_pigmac * lv4 + a5_pigmac * lv5 + d_pigmac)
    p1_megamito <- expit(a1_megamito * lv1 + a2_megamito * lv2 + a3_megamito * lv3 +
      a4_megamito * lv4 + a5_megamito * lv5 + d_megamito)
    p1_mallory <- expit(a1_mallory * lv1 + a2_mallory * lv2 + a3_mallory * lv3 +
      a4_mallory * lv4 + a5_mallory * lv5 + d_mallory)
    p1_glycogen <- expit(a1_glycogen * lv1 + a2_glycogen * lv2 + a3_glycogen * lv3 +
      a4_glycogen * lv4 + a5_glycogen * lv5 + d_glycogen)

    p0_microsteat <- 1 - p1_microsteat
    p0_microgran <- 1 - p1_microgran
    p0_lipogran <- 1 - p1_lipogran
    p0_acidophil <- 1 - p1_acidophil
    p0_pigmac <- 1 - p1_pigmac
    p0_megamito <- 1 - p1_megamito
    p0_mallory <- 1 - p1_mallory
    p0_glycogen <- 1 - p1_glycogen

    # rxode2 ordinal endpoints: the listed probabilities belong to categories
    # 1, 2, ... and the last category takes the remainder, so a simulated
    # category k is grade k - 1 (for a binary feature, 1 = absent and
    # 2 = present).
    steatosis ~ c(p0_steatosis, p1_steatosis, p2_steatosis)
    inflammation ~ c(p0_inflammation, p1_inflammation, p2_inflammation)
    ballooning ~ c(p0_ballooning, p1_ballooning)
    fibrosis ~ c(p0_fibrosis, p1_fibrosis, p2_fibrosis, p3_fibrosis)
    portal_inflammation ~ c(p0_portinfl, p1_portinfl)
    microvesicular_steatosis ~ c(p0_microsteat)
    microgranulomas ~ c(p0_microgran)
    lipogranulomas ~ c(p0_lipogran)
    acidophil_bodies ~ c(p0_acidophil)
    pigmented_macrophages ~ c(p0_pigmac)
    megamitochondria ~ c(p0_megamito)
    mallory_hyaline ~ c(p0_mallory)
    glycogen_nuclei ~ c(p0_glycogen)
  })
}
