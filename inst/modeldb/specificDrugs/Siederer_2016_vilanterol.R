Siederer_2016_vilanterol <- function() {
  description <- "Three-compartment population PK model with zero-order absorption for inhaled vilanterol in subjects with COPD and healthy volunteers, with age, bodyweight, sex, smoking and study effects on apparent clearance and central volume"
  reference <- "Siederer S, Allen A, Yang S. Population Pharmacokinetics of Inhaled Fluticasone Furoate and Vilanterol in Subjects with Chronic Obstructive Pulmonary Disease. Eur J Drug Metab Pharmacokinet. 2016;41(6):743-758. doi:10.1007/s13318-015-0303-4"
  vignette <- "Siederer_2016_fluticasoneFuroate_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Siederer 2016 reports vilanterol concentrations and
  # exposures in pg/mL and pg*h/mL (assay LLQ 10 pg/mL, raised to 20 pg/mL in
  # study HZC110946; Sect. 2.2); multiply `Cc` by 1000 to compare against the
  # published values. No scale factor is applied inside the model so that
  # dose / volume / clearance stay mutually consistent.

  covariateData <- list(
    DIS_COPD = list(
      description = "COPD patient indicator; selects the COPD stratum of CL/F, V1/F and V2/F and gates every covariate effect",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = healthy volunteer (study HZA102936)",
      notes = "1 = subject with COPD (94% of the vilanterol dataset), 0 = healthy volunteer (6%). Siederer 2016 Sect. 3.2.2 separated CL/F, V1/F and V2/F by population in the structural base model, and every covariate effect in Table 2 is labelled 'COPD' and applies to COPD subjects only.",
      source_name = "population (healthy subjects or subjects with COPD)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F and V1/F in COPD subjects, normalised to 60 years. Observed COPD range 41-84 years (Sect. 3.2.2); dataset median 61.0 years (Online Resource Table S2).",
      source_name = "age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-model effect on CL/F in COPD subjects, normalised to 70 kg. Observed COPD range 35-160 kg (Sect. 3.2.2); dataset mean 75.3 kg (Online Resource Table S2).",
      source_name = "weight"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Siederer 2016 codes sex as males = 1, females = 2 and applies the categorical covariate equation theta_COV * (covariate - 1), so the effect multiplier applies exactly when SEXF = 1. SEXF = source sex code - 1.",
      source_name = "sex (males = 1, females = 2)"
    ),
    SMOKE = list(
      description = "Current-smoker indicator at screening",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = former smoker",
      notes = "1 = current smoker at screening, 0 = former smoker. All COPD subjects had a current or prior history of at least 10 pack-years (Sect. 2.1), so the reference group is former rather than never smokers. 53% of the vilanterol dataset were current smokers (Online Resource Table S2).",
      source_name = "smoking status at screening (former or current)"
    ),
    STUDY_HZC110946 = list(
      description = "Study HZC110946 (Study 3) cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = studies HZC112206 and HZC112207 (Studies 1 and 2), the pooled Phase III reference",
      notes = "Phase III 3-way incomplete cross-over in 54 COPD subjects (NCT01072149). Retained as a covariate on V1/F because raw concentration-time data suggested higher vilanterol exposure in this study (Sect. 2.3.1).",
      source_name = "study (Study 3)"
    ),
    STUDY_HZC111348 = list(
      description = "Study HZC111348 (Study 4) cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = studies HZC112206 and HZC112207 (Studies 1 and 2), the pooled Phase III reference",
      notes = "Phase II parallel-group study in 60 COPD subjects receiving fluticasone furoate/vilanterol 400/25 ug (NCT00731822). Retained as a covariate on both CL/F and V1/F; Sect. 4 concludes the marked study difference 'may just reflect between-study variability'.",
      source_name = "study (Study 4)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "vilanterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1167,
    n_studies = 5,
    age_range = "18-84 years",
    age_median = "61.0 years",
    weight_range = "34.6-160.0 kg",
    weight_mean = "75.3 kg",
    sex_female_pct = 68,
    race_ethnicity = c(
      `White/Caucasian/European` = 82,
      `African American/African` = 3,
      `Asian - East Asian` = 5,
      `Asian - Japanese` = 3,
      `Asian - South East Asian` = 5,
      `Asian - Central/South Asian` = 1,
      `American Indian/Native Alaskan` = 1,
      `White - Arabic/North African` = 1,
      Other = 1
    ),
    disease_state = "chronic obstructive pulmonary disease (94% of subjects) pooled with healthy volunteers (6%)",
    dose_range = "vilanterol 25 or 100 ug once daily by oral inhalation, alone or as the fluticasone furoate/vilanterol combination",
    regions = "global (Argentina, Chile, Czech Republic, Estonia, Germany, Japan, Korea, Mexico, Norway, Philippines, Poland, Russia, Sweden, USA and others)",
    notes = "Three Phase III studies (HZC112206, HZC112207, HZC110946) and one Phase II study (HZC111348) in COPD, plus one Phase I study (HZA102936) in healthy volunteers. Demographics from Online Resource Table S2; mean BMI 26.0 kg/m2, mean height 170 cm, mean percent-predicted FEV1 48.8% in the COPD population. 10,807 observations, 30% below the assay LLQ, handled by the NONMEM M3 likelihood method."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Siederer 2016 Table 2 reports each THETA twice:
    # as the estimated log-scale value ('Ln estimate') and as its exponential
    # ('Estimate'). The untransformed column is used here because it carries
    # more significant figures; log() of it reproduces the Ln column
    # (e.g. log(94.6) = 4.550 vs the printed 4.55).
    #
    # CL/F, V1/F and V2/F were each estimated separately in the COPD and the
    # healthy-volunteer stratum of the same joint fit (Sect. 3.2.2: 'CL/F,
    # V1/F, volume of the peripheral compartment (V2/F), and residual error
    # were separated by population'), so each carries an explicit stratum
    # suffix. Q2/F, Q3/F, V3/F and D1 are shared across strata and keep the
    # bare canonical name.
    #
    # Compartment mapping (Table 2 footnote): V1/F -> vc, V2/F -> vp with its
    # partner Q2/F -> q, V3/F -> vp2 with its partner Q3/F -> q2.
    # -----------------------------------------------------------------------
    lcl_copd <- log(94.6); label("Apparent inhaled clearance CL/F in COPD subjects, age 60 y and 70 kg (L/h)") # Table 2 'CL/F, COPD (L/h)' = 94.6 (95% CI 90.9, 98.5; RSE 0.41%); Ln estimate 4.55; Sect. 3.2.2 'The typical value of CL/F was 94.6 L/h for a subject with COPD (aged 60 years and weighing 70 kg)'
    lcl_hvt <- log(135.6); label("Apparent inhaled clearance CL/F in healthy volunteers (L/h)") # Table 2 'CL/F, HVT (L/h)' = 135.6 (95% CI 122.7, 149.9; RSE 1.06%); Ln estimate 4.91
    lvc_copd <- log(639.0); label("Apparent central volume V1/F in COPD subjects, age 60 y, non-smoking male (L)") # Table 2 'V1/F, COPD (L)' = 639.0 (95% CI 584.1, 699.2; RSE 0.74%); Ln estimate 6.46; Sect. 3.2.2 'The typical value for V1/F (non-smoking male, aged 60 years) is predicted to be 639 L (Studies 1 and 2)'
    lvc_hvt <- log(159.2); label("Apparent central volume V1/F in healthy volunteers (L)") # Table 2 'V1/F, HVT (L)' = 159.2 (95% CI 144.0, 175.9; RSE 0.99%); Ln estimate 5.07
    lvp_copd <- log(177.7); label("Apparent first peripheral volume V2/F in COPD subjects (L)") # Table 2 'V2/F, COPD (L)' = 177.7 (95% CI 152.9, 206.4; RSE 1.52%); Ln estimate 5.18
    lvp_hvt <- log(507.8); label("Apparent first peripheral volume V2/F in healthy volunteers (L)") # Table 2 'V2/F, HVT (L)' = 507.8 (95% CI 415.7, 620.2; RSE 1.65%); Ln estimate 6.23
    lq <- log(242.3); label("Apparent intercompartmental clearance Q2/F to the first peripheral compartment (L/h)") # Table 2 'Q2/F (L/h)' = 242.3 (95% CI 219.2, 267.7; RSE 0.94%); Ln estimate 5.49
    lq2 <- log(141.2); label("Apparent intercompartmental clearance Q3/F to the second peripheral compartment (L/h)") # Table 2 'Q3/F (L/h)' = 141.2 (95% CI 125.2, 159.2; RSE 1.26%); Ln estimate 4.95
    lvp2 <- log(2100.6); label("Apparent second peripheral volume V3/F (L)") # Table 2 'V3/F' = 2100.6 (95% CI 1958.6, 2253.0; RSE 0.43%); Ln estimate 7.65. The Table 2 row header misprints the unit as 'L/h'; the footnote 'V2/F, V3/F volumes of peripheral compartment' confirms it is a volume
    ld1 <- log(0.098); label("Zero-order input duration D1 after oral inhalation (h)") # Table 2 'D1 (h)' = 0.098 (95% CI 0.092, 0.105; RSE 1.63%); Ln estimate -2.32

    # -----------------------------------------------------------------------
    # Continuous covariate effects. Siederer 2016 Sect. 3.2.2 Eq. (c) is
    # PRINTED as
    #   Ln CL = theta_1 + theta_COV x (covariate / median)
    # but that form is falsified by the paper's own downstream numbers, which
    # all require a power model, i.e. an inner logarithm:
    #   Ln P = theta_1 + theta_COV x Ln(covariate / reference)
    # 1. Sect. 3.2.2 quotes 'a reduction (47 %) in inhaled clearance ... with
    #    decreasing bodyweight (range of 160-35 kg)'. The power model gives
    #    exp(0.421 * log(35/160)) / 1 = 0.527, a 47.3% reduction; the printed
    #    linear-ratio form gives exp(0.421 * (35 - 160)/70) = 0.471, 52.9%.
    # 2. exp(4.55) = 94.6 L/h is quoted as the typical CL/F at exactly
    #    60 years and 70 kg, and exp(6.46) = 639 L as the typical V1/F at
    #    exactly 60 years. Both require the covariate term to vanish at the
    #    reference, which only the power form does (the printed form would
    #    leave theta_COV x 1 behind).
    # 3. Sect. 3.2.2 quotes a 30% decrease in V1/F over 41-84 years:
    #    exp(-0.499 * log(84/41)) = 0.699, i.e. 30.1%.
    # The reference values 60 years and 70 kg come from the Sect. 3.2.2
    # typical-value statements above, not from Online Resource Table S2
    # (median age 61.0 years, mean weight 75.3 kg).
    # -----------------------------------------------------------------------
    e_age_cl <- -0.433; label("Power exponent on (AGE/60) for CL/F in COPD subjects (unitless)") # Table 2 'Age on CL/F, COPD' Ln estimate -0.433 (95% CI -0.660, -0.206), multiplier 0.649 (95% CI 0.517, 0.814; RSE 26.8%)
    e_wt_cl <- 0.421; label("Power exponent on (WT/70) for CL/F in COPD subjects (unitless)") # Table 2 'Wt on CL/F, COPD' Ln estimate 0.421 (95% CI 0.286, 0.556), multiplier 1.52 (95% CI 1.33, 1.74; RSE 16.4%)
    e_age_vc <- -0.499; label("Power exponent on (AGE/60) for V1/F in COPD subjects (unitless)") # Table 2 'Age on V1/F, COPD' Ln estimate -0.499 (95% CI -0.911, -0.087), multiplier 0.607 (95% CI 0.402, 0.917; RSE 42.1%)

    # -----------------------------------------------------------------------
    # Categorical covariate effects. Siederer 2016 Sect. 3.2.2 Eq. (b) is
    #   Ln CL = theta_1 + theta_COV x (covariate - 1)
    # with the example coding 'sex (males = 1, females = 2)'. Every categorical
    # covariate here is two-level, so (covariate - 1) is a 0/1 indicator and
    # the coefficient below is added once when the indicator is 1. Checks:
    #   exp(-0.128) = 0.880, Sect. 3.2.2 'to be lower (12 %) in females'
    #   exp( 0.295) = 1.34,  Sect. 3.2.2 'to be increased with smoking (34 %)'
    #   639.0 * exp(-0.358) = 447 L, Sect. 3.2.2 '447 L (Study 3)'
    #   639.0 * exp(-1.24)  = 185 L, Sect. 3.2.2 '185 L (Study 4)'
    #    94.6 * exp(-0.465) = 59.4 L/h, Sect. 3.2.2 'In Study 4, the typical
    #                                    value of CL/F (59.4 L/h)'
    # -----------------------------------------------------------------------
    e_sexf_vc <- -0.128; label("Log-scale effect of female sex on V1/F in COPD subjects (unitless)") # Table 2 'Sex on V1/F, COPD' Ln estimate -0.128 (95% CI -0.25, -0.006), multiplier 0.880 (95% CI 0.779, 0.994; RSE 48.4%)
    e_smoke_vc <- 0.295; label("Log-scale effect of current smoking on V1/F in COPD subjects (unitless)") # Table 2 'Smoking on V1/F, COPD' Ln estimate 0.295 (95% CI 0.179, 0.411), multiplier 1.34 (95% CI 1.20, 1.51; RSE 20.1%)
    e_study_hzc111348_cl <- -0.465; label("Log-scale effect of study HZC111348 on CL/F in COPD subjects (unitless)") # Table 2 'Study 4 on CL/F, COPD' Ln estimate -0.465 (95% CI -0.633, -0.297), multiplier 0.628 (95% CI 0.531, 0.743; RSE 18.5%)
    e_study_hzc110946_vc <- -0.358; label("Log-scale effect of study HZC110946 on V1/F in COPD subjects (unitless)") # Table 2 'Study 3 on V1/F, COPD' Ln estimate -0.358 (95% CI -0.601, -0.115), multiplier 0.699 (95% CI 0.548, 0.891; RSE 34.6%)
    e_study_hzc111348_vc <- -1.24; label("Log-scale effect of study HZC111348 on V1/F in COPD subjects (unitless)") # Table 2 'Study 4 on V1/F, COPD' Ln estimate -1.24 (95% CI -1.51, -0.968), multiplier 0.289 (95% CI 0.221, 0.380; RSE 11.2%)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Sect. 3.2.2 states only that
    # 'Inter-individual variances (exponential model) were estimated with
    # reasonable precision (%RSE <= 25 %), with exception of ETA on volume of
    # the peripheral compartment (V3/F) where %RSE was 64 %' -- Table 2 lists
    # no OMEGA row, and neither the paper nor the Online Resource reports any
    # variance magnitude. Two etas are positively evidenced: on CL/F, because
    # Sect. 2.5 derives 'individual post hoc estimates of CL/F' per subject,
    # and on V3/F, which is named explicitly above. Both are declared here at
    # fixed(0) so the structure is preserved without inventing a variance; see
    # the vignette 'Assumptions and deviations' section. Sect. 3.2.2 also notes
    # that inter-subject variability on D1, Q2/F and V2/F 'was fixed' in the
    # structural base model.
    # -----------------------------------------------------------------------
    etalcl ~ fixed(0) # Sect. 2.5 'Individual post hoc estimates of CL/F from the final population pharmacokinetic models' - magnitude not reported
    etalvp2 ~ fixed(0) # Sect. 3.2.2 'ETA on volume of the peripheral compartment (V3/F) where %RSE was 64 %' - magnitude not reported

    # -----------------------------------------------------------------------
    # Residual error. Sect. 3.2.2: 'An additive error model described the
    # residual variability.' Fig. 3 plots 'Observed log-Value (pg/mL)' against
    # 'Population log-Prediction (pg/mL)' and the Fig. 4 VPC axis is
    # 'LN concentration (pg/mL)', so vilanterol was fitted on the natural-log
    # concentration scale and this additive error is additive on the LOG
    # scale, i.e. log-normal in linear space -- hence lnorm() rather than
    # add(). (Contrast the fluticasone furoate model, whose Fig. 1 axes are
    # untransformed pg/mL and whose additive error therefore stays additive.)
    # Sect. 3.2.2 also separated residual error by population; neither
    # stratum's SIGMA magnitude is reported, so a single log-normal SD is
    # declared at fixed(0).
    # -----------------------------------------------------------------------
    expSd <- fixed(0); label("Log-normal residual SD (fraction)") # Sect. 3.2.2 'An additive error model described the residual variability' applied on the log concentration scale (Figs. 3 and 4) - magnitude not reported
  })

  model({
    # Covariate terms, log-additive and applied only to COPD subjects; every
    # Table 2 covariate row is labelled 'COPD'.
    copdcl <-
      e_age_cl * log(AGE / 60) +
      e_wt_cl * log(WT / 70) +
      e_study_hzc111348_cl * STUDY_HZC111348
    copdvc <-
      e_age_vc * log(AGE / 60) +
      e_sexf_vc * SEXF +
      e_smoke_vc * SMOKE +
      e_study_hzc110946_vc * STUDY_HZC110946 +
      e_study_hzc111348_vc * STUDY_HZC111348

    cl <- exp(DIS_COPD * (lcl_copd + copdcl) + (1 - DIS_COPD) * lcl_hvt + etalcl)
    vc <- exp(DIS_COPD * (lvc_copd + copdvc) + (1 - DIS_COPD) * lvc_hvt)
    vp <- exp(DIS_COPD * lvp_copd + (1 - DIS_COPD) * lvp_hvt)
    vp2 <- exp(lvp2 + etalvp2)
    q <- exp(lq)
    q2 <- exp(lq2)
    d1 <- exp(ld1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Zero-order input into the central compartment. Dose records must carry
    # rate = -2 so that rxode2 uses this modelled duration.
    dur(central) <- d1

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
