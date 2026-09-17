Ambery_2015_batefenterol_fev1 <- function() {
  description <- paste(
    "Empirical Emax dose-response model for day-29 trough forced expiratory",
    "volume in 1 s (FEV1) after 4 weeks of inhaled GSK961081 (batefenterol)",
    "DISKUS in patients with moderate-to-severe chronic obstructive pulmonary",
    "disease (COPD) (Ambery 2015, substudy of GSK MAB115032 / NCT01319019).",
    "A landmark (static, time-independent) model: the response is driven by",
    "the TOTAL DAILY DOSE rather than by a plasma concentration, so once- and",
    "twice-daily arms with the same daily dose are pooled, and the model",
    "carries no ODE state. The zero-dose FEV1 is an intercept scaled by the",
    "patient's own day-1 baseline trough FEV1 normalised to the population",
    "median, and the drug effect is an Emax term in total daily dose. The",
    "companion plasma PK model from the same paper is",
    "modellib('Ambery_2015_batefenterol'); the two were fitted separately and",
    "are not linked through exposure."
  )

  reference <- paste(
    "Ambery CL, Wielders P, Ludwig-Sengpiel A, Chan R, Riley JH.",
    "Population Pharmacokinetics and Pharmacodynamics of GSK961081",
    "(Batefenterol), a Muscarinic Antagonist and beta2-Agonist, in",
    "Moderate-to-Severe COPD Patients: Substudy of a Randomized Trial.",
    "Drugs R D. 2015;15(3):281-291. doi:10.1007/s40268-015-0104-x.",
    "PMID 26286203; PMCID PMC4561049.",
    sep = " "
  )

  vignette <- "Ambery_2015_batefenterol"

  units <- list(
    time = "day (placeholder; this is a day-29 landmark dose-response and carries no time dimension)",
    dosing = "ug/day (total daily batefenterol dose supplied as the DOSE_BATEFENTEROL_UGD covariate; the model consumes no rxode2 dose events)",
    concentration = "L (day-29 trough FEV1, a lung volume rather than a drug concentration)"
  )

  covariateData <- list(
    DOSE_BATEFENTEROL_UGD = list(
      description = "Total daily inhaled batefenterol (GSK961081) dose",
      units = "ug/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The dose regressor of the Emax term in Ambery 2015 Eq. 1 ('Total",
        "Dose'). Set to 0 for placebo. Ambery 2015 Fig. 3d and Fig. 4",
        "plot the model against total daily dose on the six observed",
        "levels 0, 100, 200, 400 and 800 ug/day: 100, 400 and 800 ug once",
        "daily and 100, 200 and 400 ug twice daily, so 200 ug/day arises",
        "only from 100 ug twice daily and 800 ug/day only from 800 ug once",
        "daily, while 400 ug/day pools 400 ug once daily with 200 ug twice",
        "daily. Ambery 2015 Sect. 4 states there was 'no apparent influence",
        "of dosing regimen on the PD model', which is why the once-daily",
        "and twice-daily arms were combined onto a single total-daily-dose",
        "axis."
      ),
      source_name = "Total Dose (Ambery 2015 Eq. 1; Fig. 3d and Fig. 4 x-axis)"
    ),
    FEV1_BL = list(
      description = "Patient's own day-1 baseline trough FEV1",
      units = "L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The sole retained covariate of the PD model (Ambery 2015 Sect. 3.2,",
        "Eq. 2), entering as CON * (FEV1_BL / median FEV1_BL) multiplying the",
        "zero-dose intercept. Ambery 2015 Table 1 reports 1.31 +/- 0.46 L",
        "(mean +/- SD) for the n = 347 PD analysis set. The MEDIAN used as",
        "the normalising constant in Eq. 2 is not reported anywhere in the",
        "paper; the reported mean of 1.31 L is used in its place and the",
        "substitution is recorded in the vignette's Assumptions and",
        "deviations section. Because the normalisation cancels at the median,",
        "the substitution moves predictions only for patients away from the",
        "cohort centre. Named FEV1_BL rather than the bare FEV1 covariate",
        "canonical because this model also OUTPUTS FEV1, and a covariate",
        "column cannot share a name with the output state it drives."
      ),
      source_name = "Baseline FEV1 on Day 1 (Ambery 2015 Eq. 2; Table 1)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at screening",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened by stepwise forward addition (p < 0.01) and backward",
        "elimination (p < 0.001) on the PD model (Ambery 2015 Sect. 2.4) and",
        "not retained; only baseline FEV1 on day 1 survived. Reported in",
        "Table 1 as 63 +/- 8.2 years for the n = 347 PD analysis set."
      ),
      source_name = "Age (years)"
    ),
    WT = list(
      description = "Body weight at screening",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on the PD model and not retained (Ambery 2015 Sect. 2.4).",
        "Reported in Table 1 as 76 +/- 14 kg for the n = 347 PD analysis set."
      ),
      source_name = "Weight (kg)"
    ),
    HT = list(
      description = "Height at screening",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on the PD model and not retained (Ambery 2015 Sect. 2.4).",
        "Reported in Table 1 as 171 +/- 8.5 cm for the n = 347 PD analysis",
        "set."
      ),
      source_name = "Height (cm)"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Screened on the PD model and not retained (Ambery 2015 Sect. 2.4).",
        "Table 1 reports 65% male, i.e. 35% female, for the n = 347 PD",
        "analysis set."
      ),
      source_name = "Male (%)"
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = former smoker",
      notes = paste(
        "Screened on the PD model and not retained (Ambery 2015 Sect. 2.4).",
        "Table 1 reports 49% current smokers for the n = 347 PD analysis set."
      ),
      source_name = "Current smoker (%)"
    ),
    CONMED_ICS = list(
      description = "Concurrent inhaled-corticosteroid use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concurrent ICS",
      notes = paste(
        "Screened on the PD model and not retained (Ambery 2015 Sect. 2.4).",
        "Table 1 reports 58% concurrent ICS use for the n = 347 PD analysis",
        "set; patients on a stable inhaled-corticosteroid dose were eligible",
        "for enrollment (Sect. 2.2)."
      ),
      source_name = "Concurrent ICS use (%)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 347L,
    n_studies = 1L,
    n_observations = 347L,
    age_mean_sd = "63 +/- 8.2 years",
    weight_mean_sd = "76 +/- 14 kg",
    height_mean_sd = "171 +/- 8.5 cm",
    sex_female_pct = 35,
    disease_state = paste(
      "Moderate-to-severe stable COPD: post-bronchodilator FEV1/FVC below",
      "70% and FEV1 30-70% of predicted by NHANES III normative values;",
      "current or former smokers aged 40 years or above with at least a",
      "10 pack-year history; diagnosed asthma excluded. Baseline trough",
      "FEV1 on day 1 was 1.31 +/- 0.46 L"
    ),
    dose_range = paste(
      "GSK961081 DISKUS 100, 400 and 800 ug once daily and 100, 200 and",
      "400 ug twice daily for 4 weeks, plus placebo; total daily doses of",
      "0, 100, 200, 400 and 800 ug/day"
    ),
    regimens = "Inhaled dry powder (DISKUS), once or twice daily for 28 days",
    regions = "Multicenter international (GSK MAB115032 / NCT01319019)",
    notes = paste(
      "One trough-FEV1 observation per patient: 347 day-29 trough FEV1",
      "values from 347 patients (Ambery 2015 Sect. 3.2). Day-29 trough FEV1",
      "was defined as the mean of the 11 h and 12 h measurements after the",
      "evening dose on day 28 (Sect. 2.3), measured on a Vitalograph",
      "spirometer. The parent trial was a 4-week multicenter randomized",
      "double-blind double-dummy placebo- and salmeterol-controlled",
      "parallel-group study; only the GSK961081 and placebo arms enter this",
      "model. The PD model was fitted in NONMEM 7 with FOCE-I using all",
      "available data. Table 1 demographics above are for the PD / PK-PD",
      "analysis set."
    )
  )

  ini({
    # ================================================================
    # Ambery 2015 Eq. 1 and Eq. 2 (Sect. 3.2):
    #
    #   Effect = (E0 * COV) + Emax * TotalDose / (TotalDose + ED50)
    #   COV    = CON * (Baseline FEV1 on Day 1 / Median Baseline FEV1
    #                   on Day 1)
    #
    # 'where Emax is the trough FEV1 at the maximum effect, E0 is the
    # FEV1 at zero dose, ED50 is the dose producing 50% of the maximum
    # effect, and CON is the baseline FEV1 covariate effect.'
    #
    # E0 and CON enter as a PRODUCT, which is how the paper writes it:
    # Table 3's 'Intercept for E0' (0.0650 L) and 'Baseline' (19.0)
    # multiply to 1.235 L, the zero-dose day-29 trough FEV1 of a
    # median-baseline patient. That product is reproduced by Fig. 4a,
    # whose model line sits at about 1.28 L at a total daily dose of 0,
    # and the redundancy of the parameterisation is visible in the two
    # rows' similar and relatively poor precision (RSE 32.8% and 36.9%).
    # ================================================================

    le0 <- log(0.0650)
    label("Intercept of the zero-dose day-29 trough FEV1 (L); multiplied by the baseline-FEV1 covariate term COV")  # Table 3, row 'Intercept for E0 (L)' = 0.0650 (RSE 32.8%, 95% CI 0.0233 to 0.107)

    e_fev1_bl_e0 <- 19.0
    label("Baseline-FEV1 covariate coefficient CON on the zero-dose intercept (unitless)")  # Table 3, row 'Baseline (L)' = 19.0 (RSE 36.9%, 95% CI 5.24 to 32.8); CON of Eq. 2

    lemax <- log(0.293)
    label("Maximum increment in day-29 trough FEV1 above the zero-dose value (L)")  # Table 3, row 'Emax (L)' = 0.293 (RSE 14.9%, 95% CI 0.207 to 0.379)

    led50 <- log(152)
    label("Total daily dose producing 50% of the maximum effect ED50 (ug/day)")  # Table 3, row 'ED50 (ug)' = 152 (RSE 50.2%, 95% CI 2.45 to 302)

    # ================================================================
    # Residual error. Ambery 2015 Table 3 reports a proportional
    # residual error of 0.204 with a properly computed RSE column
    # (0.0809-0.327 is 0.204 +/- 1.96 * 0.0628 and 0.0628 / 0.204 =
    # 30.8%), and the value is a standard deviation -- a 20.4% CV --
    # rather than a variance. A variance reading would imply a 45.2%
    # residual CV, which is refuted by two of the paper's own figures:
    # in Fig. 3a the observed day-29 trough FEV1 scatters only about
    # +/-0.4 L about a population prediction of 1.0 L (20% CV or less,
    # and that spread also contains any inter-individual term), and in
    # Fig. 4a the 95% prediction interval at zero dose spans roughly
    # 0.65-2.35 L about a median of 1.28 L, which the variance reading
    # would overshoot on both sides once the +/-0.46 L spread of the
    # baseline-FEV1 covariate is included.
    #
    # NOTE that the companion PK model in the same paper reports its
    # variability terms as VARIANCES (see
    # modellib('Ambery_2015_batefenterol')); Table 2 and Table 3 do not
    # share a convention. Each was adjudicated separately against the
    # matching figure.
    # ================================================================
    propSd_FEV1 <- 0.204
    label("Proportional residual SD on day-29 trough FEV1 (fraction)")  # Table 3, row 'Proportional residual error' = 0.204 (RSE 30.8%, 95% CI 0.0809 to 0.327)

    # ================================================================
    # No inter-individual variability is encoded. Ambery 2015 Table 3
    # lists no IIV term, even though Fig. 3b (observed vs individual
    # predictions) is materially tighter than Fig. 3a (observed vs
    # population predictions), which can only happen if the fitted
    # model carried an eta. With one observation per patient that eta
    # is not separately identifiable from the residual in any case, and
    # no variance is reported for it, so none is invented here. The gap
    # is recorded in the vignette's Assumptions and deviations section.
    # ================================================================
  })

  model({
    emax <- exp(lemax)
    ed50 <- exp(led50)
    e0   <- exp(le0)

    # ----------------------------------------------------------------
    # Ambery 2015 Eq. 2. The median day-1 baseline trough FEV1 is not
    # reported; the Table 1 MEAN for the n = 347 PD analysis set,
    # 1.31 L, is substituted (see the FEV1_BL covariateData notes).
    # ----------------------------------------------------------------
    cov_fev1 <- e_fev1_bl_e0 * (FEV1_BL / 1.31)

    # ----------------------------------------------------------------
    # Ambery 2015 Eq. 1. At DOSE_BATEFENTEROL_UGD = 0 the Emax term is
    # exactly 0, so a placebo patient's predicted day-29 trough FEV1 is
    # e0 * cov_fev1.
    # ----------------------------------------------------------------
    FEV1 <-
      e0 * cov_fev1 +
      emax * DOSE_BATEFENTEROL_UGD / (DOSE_BATEFENTEROL_UGD + ed50)

    FEV1 ~ prop(propSd_FEV1)
  })
}
