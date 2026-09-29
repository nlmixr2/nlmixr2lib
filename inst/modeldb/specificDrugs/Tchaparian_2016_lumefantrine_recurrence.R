# Adjusted Cox proportional-hazards model for 28-day recurrent Plasmodium
# falciparum parasitaemia after artemether-lumefantrine in Ugandan children aged
# 6 months to 2 years, driven by a landmark lumefantrine concentration read off
# the companion population PK fit. Source: Tchaparian E, Sambol NC, Arinaitwe E,
# et al. J Infect Dis. 2016;214(8):1243-1251. doi:10.1093/infdis/jiw338
# (PMC5034953).

Tchaparian_2016_lumefantrine_recurrence <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for the hazard of",
    "recurrent Plasmodium falciparum parasitaemia within 28 days of treatment",
    "with artemether-lumefantrine in 100 Ugandan children aged 6 months to 2",
    "years (222 malaria episodes), adjusted for age group, residence,",
    "underweight status, parasite density and haemoglobin (Tchaparian 2016",
    "Table 3). The exposure driver is the day 7 (168 h) capillary whole-blood",
    "lumefantrine concentration dichotomised at 200 ng/mL, taken as a",
    "per-episode covariate from the companion population PK model",
    "(modellib('Tchaparian_2016_lumefantrine')). The paper's headline finding",
    "is an interaction between that exposure and concurrent",
    "trimethoprim-sulfamethoxazole prophylaxis: among children NOT receiving",
    "prophylaxis a day 7 concentration below 200 ng/mL carries a 2.97-fold",
    "higher hazard of recurrence, whereas among children receiving it the same",
    "threshold carries a hazard ratio of 0.50 that does not reach significance,",
    "so the exposure-response relationship is present only in the absence of",
    "prophylaxis. Two further outputs carry the paper's comparative analysis of",
    "the day 3 (72 h) against the day 7 landmark on a continuous natural-log",
    "concentration scale in the no-prophylaxis subgroup, where a one-unit rise",
    "in log concentration reduces the hazard by 49 percent on day 3 but only 20",
    "percent on day 7 - the basis for the paper's suggestion that the earlier,",
    "logistically easier sampling day is the better predictor.",
    "NO BASELINE HAZARD IS ENCODED: a Cox regression is semiparametric, so",
    "h0(t) is left completely unspecified by the method rather than being an",
    "unreported parameter. This model returns RELATIVE hazards only and",
    "deliberately offers no survivor function; multiply by any user-supplied",
    "h0(t) to obtain a subject hazard. The model is algebraic and deterministic",
    "(no ODE state, no drug input, no IIV, no residual error).",
    "Companion pharmacokinetic model from the same paper:",
    "modellib('Tchaparian_2016_lumefantrine').",
    sep = " "
  )
  reference <- paste(
    "Tchaparian E, Sambol NC, Arinaitwe E, McCormack SA, Bigira V, Wanzira H,",
    "Muhindo M, Creek DJ, Sukumar N, Blessborn D, Tappero JW, Kakuru A,",
    "Bergqvist Y, Aweeka FT, Parikh S.",
    "Population pharmacokinetics and pharmacodynamics of lumefantrine in young",
    "Ugandan children treated with artemether-lumefantrine for uncomplicated",
    "malaria. J Infect Dis. 2016;214(8):1243-1251. doi:10.1093/infdis/jiw338.",
    "The adjusted hazard ratios are from Table 3, 'Adjusted HR (95% CI)'",
    "column; the continuous day 3 and day 7 log-concentration hazard ratios are",
    "from Results, 'Day 3 or 7 Lumefantrine Concentration and Clinical",
    "Outcomes', corroborated by the Discussion which prints them as 0.51 and",
    "0.80. The landmark sampling times in hours after the first dose are from",
    "Methods, 'Sample Collection and Analysis'. The Cox regression used a robust",
    "sandwich estimator to account for repeated episodes within a child",
    "(Methods, 'Association Analysis Between Day 7 Lumefantrine Concentration",
    "and Recurrent Malaria').",
    sep = " "
  )
  vignette <- "Tchaparian_2016_lumefantrine"
  units <- list(
    time = "n/a (semiparametric Cox relative hazard; the baseline hazard and therefore the time scale are unspecified by the fit)",
    dosing = "n/a (exposure-response model; lumefantrine exposure enters through the CONC_LUMEFANTRINE_168H and CONC_LUMEFANTRINE_72H covariates, not as a dose record)",
    concentration = "hr (hazard of 28-day recurrent parasitaemia relative to the reference cohort, unitless); not a drug concentration"
  )

  covariateData <- list(
    CONC_LUMEFANTRINE_168H = list(
      description = paste(
        "Individual capillary whole-blood lumefantrine concentration at the day 7",
        "landmark, 168 h after the first of the six artemether-lumefantrine doses.",
        "Per-episode, time-fixed.",
        sep = " "
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = "greater than or equal to 200 ng/mL (the adequate-exposure group, which is the printed reference row of Table 3)",
      notes = paste(
        "Enters the hazard only through the dichotomy at 200 ng/mL, so the model",
        "is insensitive to the value except through which side of the threshold",
        "it falls on. The threshold is carried as the separate fixed() parameter",
        "conc_cut rather than being baked into this column, following the",
        "Rayner_2013_oseltamivir_shedding precedent. Tchaparian 2016 selected",
        "200 ng/mL by ROC analysis (area under the ROC curve 0.684, best for",
        "children not receiving prophylaxis) and note it reproduces the cutoff",
        "the WWARN individual-patient meta-analysis had already identified.",
        "Methods pin the landmark to '108 hours after the last dose'; the last of",
        "the six twice-daily doses is at 60 h, so the nominal landmark is 168 h",
        "after the first dose. The Supplement reports that the ACTUAL median",
        "recorded collection time was 6.89 days (IQR 6.82-6.94), i.e. about",
        "165 h, marginally earlier than nominal. Observed distribution: median",
        "216 ng/mL (IQR 136-345, n = 216 measured samples); the companion",
        "population PK model predicts a median of 202.4 ng/mL (IQR 143.1-321.6)",
        "at this landmark. Generate this column with",
        "modellib('Tchaparian_2016_lumefantrine').",
        sep = " "
      ),
      source_name = "day 7 capillary whole-blood lumefantrine concentration"
    ),
    CONC_LUMEFANTRINE_72H = list(
      description = paste(
        "Individual capillary whole-blood lumefantrine concentration at the day 3",
        "landmark, 72 h after the first of the six artemether-lumefantrine doses.",
        "Per-episode, time-fixed.",
        sep = " "
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = "n/a -- enters the comparative day 3 output on a continuous natural-log scale, referenced to the observed median of 2777 ng/mL carried as conc_day3_ref",
      notes = paste(
        "Used only by the hr_c3 output, which carries the paper's comparative",
        "day-3-versus-day-7 analysis in the subgroup NOT receiving",
        "trimethoprim-sulfamethoxazole prophylaxis (182 episodes with both",
        "landmarks available). Methods pin the landmark to '12 hours after the",
        "last dose', i.e. 72 h after the first dose; the Supplement reports an",
        "actual median recorded collection time of 2.86 days (IQR 2.80-2.93),",
        "about 69 h. Observed distribution: median 2777 ng/mL (IQR 1672-4760,",
        "n = 187 measured samples), so this landmark sits near the peak of the",
        "profile and is roughly an order of magnitude above the day 7 value.",
        "Unlike the day 7 column this one is NOT dichotomised - the paper reports",
        "it per one-unit change in natural log concentration. Generate this",
        "column with modellib('Tchaparian_2016_lumefantrine').",
        sep = " "
      ),
      source_name = "day 3 capillary whole-blood lumefantrine concentration"
    ),
    CONMED_TMPSMX = list(
      description = "Concurrent daily trimethoprim-sulfamethoxazole prophylaxis indicator; 1 = receiving prophylaxis, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prophylaxis), the printed reference row of Table 3",
      notes = paste(
        "The effect modifier that carries the paper's headline finding. Daily",
        "prophylaxis was given to all HIV-infected participants and to",
        "HIV-exposed participants until breast-feeding ended, after which",
        "HIV-uninfected children were randomised to continue to 2 years of age or",
        "to stop (Methods, 'Study Area and Patient Enrollment'). 43 children over",
        "85 of the 222 analysed episodes received it (Table 1, PK outcomes",
        "column). Tchaparian 2016 report a significant interaction with",
        "lumefantrine exposure (P = .0005 adjusted), which is why the",
        "low-exposure effect is carried as two stratum-specific coefficients",
        "rather than one main effect plus an interaction term. Note that",
        "prophylaxis itself protects against malaria independently, which the",
        "authors give as the likely reason the exposure-response signal",
        "disappears in the treated stratum. This covariate was screened but NOT",
        "retained in the companion population PK model, a point the Discussion",
        "flags explicitly even though measured day 7 concentrations were higher",
        "in children receiving prophylaxis (median 243 vs 206 ng/mL, P = .018).",
        sep = " "
      ),
      source_name = "TMP-SMZ use"
    ),
    AGE = list(
      description = "Subject age at the time of malaria diagnosis, entering as the paper's three-level age group rather than as a continuous term.",
      units = "years",
      type = "continuous",
      reference_category = "6 to less than 12 months (0.5 to less than 1.0 years), the printed reference row of Table 3",
      notes = paste(
        "Table 3 reports age as a three-level categorical covariate with",
        "boundaries at 12 and 18 months, so AGE is supplied in the canonical unit",
        "of years and the model derives the two indicator variables using the",
        "cut points age_cut1 = 1.0 and age_cut2 = 1.5 years. The reference",
        "stratum is the youngest. Neither non-reference level is significant",
        "(P = .13 and P = .34; group P = .29) and the point estimates are",
        "non-monotonic - the middle 12-18 month group carries the larger hazard",
        "ratio (1.87) than the oldest group (1.53) - so this covariate is carried",
        "as a published adjustment term, not as an interpretable age trend.",
        "Cohort median 14.9 months, range 6.6-24.2 months (Table 1, PK outcomes",
        "column), with 49, 126 and 47 episodes in the three strata. Contrast the",
        "companion population PK model, where age is the single retained",
        "covariate and acts as a continuous power function on bioavailability.",
        sep = " "
      ),
      source_name = "Age group, mo"
    ),
    RESIDENCE_RURAL = list(
      description = "Settlement type of the subject's residence; 1 = rural, 0 = urban.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (urban), the printed reference row of Table 3",
      notes = paste(
        "Listed among the covariates evaluated for the recurrence analysis",
        "(Methods, 'Association Analysis Between Day 7 Lumefantrine",
        "Concentration and Recurrent Malaria') and retained in the adjusted",
        "model. Only 16 of 222 episodes (7.2 percent) involved urban-resident",
        "children (Table 1, PK outcomes column), so the reference stratum is much",
        "the smaller of the two and the estimate is correspondingly imprecise:",
        "adjusted hazard ratio 1.21 with a 95 percent confidence interval of",
        "0.43-3.43 and P = .71. Carried because Figure 4 and Supplemental Figure",
        "S3 both state that the cumulative-risk curves were adjusted for",
        "residence. Not recorded at all in the population PK data set, where",
        "Table 1 prints 'Not included' for this row.",
        sep = " "
      ),
      source_name = "Residence"
    ),
    HGB = list(
      description = "Blood haemoglobin concentration measured on the day of malaria diagnosis.",
      units = "g/dL",
      type = "continuous",
      reference_category = "10.0 g/dL, the cohort median (Table 1, PK outcomes column), used as the centring value hgb_ref",
      notes = paste(
        "Table 3 footnote a marks this as a continuous per-unit covariate;",
        "the adjusted hazard ratio of 0.93 is therefore per 1 g/dL increase.",
        "Not significant (95 percent confidence interval 0.83-1.05, P = .25).",
        "The paper does not state a centring value, so the model centres on the",
        "cohort median of 10.0 g/dL (range 5.6-15.9) so that the relative hazard",
        "is 1 at a median child; this is a presentational choice that leaves",
        "every hazard RATIO between two haemoglobin values unchanged. Screened",
        "but not retained in the companion population PK model.",
        sep = " "
      ),
      source_name = "Hemoglobin level (g/dL)"
    ),
    PARA = list(
      description = "Asexual Plasmodium parasite density measured on the day of malaria diagnosis, entering the hazard on a natural-log scale.",
      units = "parasites/uL",
      type = "continuous",
      reference_category = "17603 parasites/uL, the cohort geometric mean (Table 1, PK outcomes column), used as the centring value para_ref",
      notes = paste(
        "Table 3 footnote a marks this as a continuous per-unit covariate and",
        "Methods describe it as 'parasite density (log transformed)'. The paper",
        "does not name the logarithm base for this term. The model uses the",
        "NATURAL log, because that is the base the same paper names explicitly",
        "for the day 3 and day 7 concentration analyses in the immediately",
        "following subsection; see the vignette's 'Assumptions and deviations'",
        "section, which also states what a log10 reading would change. The",
        "adjusted hazard ratio of 1.04 per log unit is not significant (95",
        "percent confidence interval 0.95-1.14, P = .35). Centred on the cohort",
        "geometric mean of 17,603 parasites/uL (95 percent confidence interval",
        "13,762-22,516) so that the relative hazard is 1 at a typical child.",
        "The model guards the logarithm with max(PARA, 1) so that densities below",
        "1 parasite/uL cannot produce a negative-infinite linear predictor, the",
        "same guard Kloprogge_2014_quinine.R applies to this column.",
        sep = " "
      ),
      source_name = "Parasite density (parasites/uL)"
    ),
    WAZ = list(
      description = "Weight-for-age z-score computed against the WHO Child Growth Standards, entering as the paper's dichotomised underweight indicator.",
      units = "unitless (z-score; standard-deviation units)",
      type = "continuous",
      reference_category = "greater than or equal to -2 (not underweight), the printed reference row of Table 3",
      notes = paste(
        "Table 3 footnote b marks underweight as dichotomised; Methods give the",
        "cut point as 'a weight-for-age z score with a -2 cutoff', computed using",
        "World Health Organization standards. WAZ is therefore supplied on the",
        "canonical continuous z-score scale and the model derives the indicator",
        "using the fixed() cut point waz_cut = -2. 29 of 222 episodes (13.1",
        "percent) were in underweight children (Table 1, PK outcomes column).",
        "The adjusted hazard ratio of 1.04 is essentially null (95 percent",
        "confidence interval 0.47-2.31, P = .93), which is notable given that",
        "underweight children received a HIGHER median milligram-per-kilogram",
        "dose (90.0 vs 77.4 mg/kg, P < .001) and that the WWARN meta-analysis the",
        "paper cites had flagged low weight-for-age as the main risk factor for",
        "underexposure. Screened but not retained in the companion population PK",
        "model.",
        sep = " "
      ),
      source_name = "Underweight (dichotomized)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    n_episodes = 222L,
    age_range = "6.6-24.2 months (median 14.9); enrolment window 6 months to 2 years",
    weight_range = "6.1-13.3 kg (median 9.0)",
    sex_female_pct = 52.2,
    race_ethnicity = "Ugandan children resident in Tororo district; race / ethnicity not otherwise reported",
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria diagnosed on a positive thick",
      "blood smear plus documented or recent fever. Unlike the companion",
      "population PK analysis, only P. falciparum episodes were eligible for this",
      "pharmacodynamic analysis. 8 children were HIV-infected and 43 received",
      "daily trimethoprim-sulfamethoxazole prophylaxis because they were",
      "HIV-infected or HIV-exposed.",
      sep = " "
    ),
    dose_range = paste(
      "Artemether-lumefantrine (Coartem, Novartis) 20 mg artemether plus 120 mg",
      "lumefantrine (1 tablet) twice daily for 3 days - 6 doses, 720 mg total",
      "lumefantrine - for every child, because all body weights were under 14 kg.",
      "Median total body-weight-adjusted lumefantrine dose 80.0 mg/kg (range",
      "54.1-118.0).",
      sep = " "
    ),
    regions = "Tororo district, eastern Uganda (perennial high-transmission area, entomological inoculation rate up to 562 infective bites per person-year)",
    notes = paste(
      "The endpoint is the 28-day risk of recurrent parasitaemia, assessed at",
      "clinic visits on days 0, 1, 2, 3, 7, 14, 21 and 28, so observed event",
      "times are coarsely quantised. Of the 249 episodes enrolled in the PK/PD",
      "study, 222 in 100 children were analysed; the Supplement itemises the 27",
      "exclusions as non-falciparum species (14), elevated pre-first-dose",
      "concentrations (9), no outcome classification (3) and medication",
      "non-compliance (1). Outcome classification among the 137 episodes without",
      "prophylaxis was 59.1 percent adequate clinical and parasitological",
      "response, 19.7 percent late clinical failure and 21.2 percent late",
      "parasitological failure; among the 85 with prophylaxis it was 65.9, 28.2",
      "and 5.9 percent, so the prophylaxis stratum fared better overall, which is",
      "the direction of the prophylaxis main effect encoded here. At 63 days",
      "(secondary outcome) 65.8 percent of episodes were followed by recurrent",
      "malaria, of which 95.9 percent genotyped as new infection rather than",
      "recrudescence - so this hazard is dominated by reinfection in a",
      "high-transmission setting and should be read as a post-treatment",
      "prophylaxis endpoint rather than as a cure rate. For the PD analysis",
      "concentrations below the limit of detection were treated as 0 and those",
      "below the limit of quantification were retained as measured, which differs",
      "from the left-censored treatment used in the companion PK fit. Companion",
      "pharmacokinetic model:",
      "modellib('Tchaparian_2016_lumefantrine').",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Tchaparian 2016 Table 3, 'Adjusted HR (95% CI)' column: a Cox
    # proportional-hazards regression on recurrent malaria by day 28,
    # fitted with a robust sandwich estimator to account for repeated
    # episodes within a child.
    #
    #   h(t) = h0(t) * exp(
    #       e_conc_lumefantrine_168h_low_haz_nots * low * (1 - CONMED_TMPSMX)
    #     + e_conc_lumefantrine_168h_low_haz_ts   * low * CONMED_TMPSMX
    #     + e_conmed_tmpsmx_haz                   * CONMED_TMPSMX
    #     + e_age_haz_12to18   * age12to18
    #     + e_age_haz_ge18     * agege18
    #     + e_residence_rural_haz * RESIDENCE_RURAL
    #     + e_hgb_haz  * (HGB - hgb_ref)
    #     + e_para_haz * (log(max(PARA, 1)) - log(para_ref))
    #     + e_waz_underweight_haz * underweight )
    #
    # The coefficients below are natural logs of the printed hazard
    # ratios, which is the scale a Cox model estimates on.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is left completely unspecified by the method, so it is not
    # an unreported parameter but a quantity the fit never produced.
    # This model therefore returns RELATIVE hazards only and
    # deliberately does not offer a survivor function. Figure 4 gives
    # the paper's four-way Kaplan-Meier cumulative-risk curves as an
    # empirical description of the absolute time course. This matches
    # the deterministic algebraic pattern of
    # Rayner_2013_oseltamivir_shedding.R and Liu_2024_saf189s_pfs.R.
    #
    # The regression was fitted in Stata SE12.1 and SAS 9.4, not in
    # NONMEM or Monolix. The paper reports hazard ratios with 95
    # percent confidence intervals and no variance components, so there
    # is no IIV and no residual error to encode, and no observation
    # endpoint is declared.
    # ==================================================================

    # ----- Exposure dichotomy -----
    # Selected by the authors from an ROC analysis rather than estimated
    # jointly with the hazard ratios, hence fixed().
    conc_cut <- fixed(200)
    label("Day 7 lumefantrine concentration cutoff separating the low- and adequate-exposure groups (ng/mL)")
    # Table 3 row heading 'Day 7 lumefantrine concentration <200 ng/mL';
    # Results: 'a day 7 capillary whole-blood lumefantrine concentration
    # of approximately 200 ng/mL was the optimal cutoff for predicting
    # the risk of 28-day recurrence' (area under the ROC curve 0.684).

    # ----- Stratum-specific low-exposure log hazard ratios -----
    # Table 3 reports the low-exposure effect SEPARATELY within each
    # prophylaxis stratum because the interaction is significant
    # (P = .0005 adjusted, P = .001 unadjusted). Carrying the two
    # stratum-specific coefficients directly, rather than a main effect
    # plus an interaction offset, reproduces both printed numbers with
    # no algebra and is the encoding the stratum-suffix grammar of
    # references/parameter-names.md prescribes.
    e_conc_lumefantrine_168h_low_haz_nots <- log(2.97)
    label("Log hazard ratio for 28-day recurrent parasitaemia, day 7 lumefantrine below 200 ng/mL vs at or above it, among children NOT receiving trimethoprim-sulfamethoxazole prophylaxis (log scale; HR 2.97)")
    # Table 3, 'Without TMP-SMZ use' block: adjusted HR 2.97
    # (95 percent CI 1.59-5.55), P = .0007. This is the paper's headline
    # result and is quoted in the Abstract as a '3-fold higher hazard'.

    e_conc_lumefantrine_168h_low_haz_ts <- log(0.50)
    label("Log hazard ratio for 28-day recurrent parasitaemia, day 7 lumefantrine below 200 ng/mL vs at or above it, among children receiving trimethoprim-sulfamethoxazole prophylaxis (log scale; HR 0.50)")
    # Table 3, 'With TMP-SMZ use' block: adjusted HR 0.50
    # (95 percent CI 0.22-1.13), P = .10. Not significant; the
    # Discussion describes it as 'a trend towards an opposite
    # interaction' and cautions against over-reading it.

    # ----- Prophylaxis main effect -----
    e_conmed_tmpsmx_haz <- log(0.75)
    label("Log hazard ratio for 28-day recurrent parasitaemia, receiving trimethoprim-sulfamethoxazole prophylaxis vs not, at a day 7 lumefantrine concentration at or above 200 ng/mL (log scale; HR 0.75)")
    # Table 3, 'TMP-SMZ use / Yes' row: adjusted HR 0.75 (95 percent CI
    # 0.43-1.29), P = .30. NOTE that this row sits in the table's
    # main-effect block while the low-exposure effect above is reported
    # within strata, so the two are not printed on a single consistent
    # parameterisation; see the vignette's 'Assumptions and deviations'
    # section. The direction is corroborated twice independently: the
    # prophylaxis stratum had a higher adequate-response rate (65.9 vs
    # 59.1 percent, Results) and the Discussion attributes it to 'the
    # independent protective effects of TMP-SMZ against malaria'.

    # ----- Age group, reference 6 to <12 months -----
    age_cut1 <- fixed(1)
    label("Lower age-group boundary separating the 6-to-under-12-month reference stratum from the 12-18 month stratum (years)")
    # Table 3 'Age group, mo' rows: strata are 6-12, 12-18 and >=18
    # months, so the boundary is 12 months = 1 year.

    age_cut2 <- fixed(1.5)
    label("Upper age-group boundary separating the 12-18 month stratum from the at-or-over-18-month stratum (years)")
    # Table 3 'Age group, mo' rows: 18 months = 1.5 years.

    e_age_haz_12to18 <- log(1.87)
    label("Log hazard ratio for 28-day recurrent parasitaemia, age 12-18 months vs 6 to under 12 months (log scale; HR 1.87)")
    # Table 3, 'Age group, mo / 12-18' row: adjusted HR 1.87
    # (95 percent CI 0.83-4.22), P = .13.

    e_age_haz_ge18 <- log(1.53)
    label("Log hazard ratio for 28-day recurrent parasitaemia, age at or over 18 months vs 6 to under 12 months (log scale; HR 1.53)")
    # Table 3, 'Age group, mo / >=18' row: adjusted HR 1.53
    # (95 percent CI 0.64-3.63), P = .34. Note the non-monotonic
    # ordering against the 12-18 month estimate above; the group P value
    # is .29.

    # ----- Residence, reference urban -----
    e_residence_rural_haz <- log(1.21)
    label("Log hazard ratio for 28-day recurrent parasitaemia, rural vs urban residence (log scale; HR 1.21)")
    # Table 3, 'Residence / Rural' row: adjusted HR 1.21
    # (95 percent CI 0.43-3.43), P = .71. Coincidence warning for future
    # readers: the table's day-7-concentration MAIN-effect row carries
    # the same printed value of 1.21 (95 percent CI 0.71-2.07, P = .48).
    # The two are different rows with different confidence intervals and
    # must not be conflated; the value used here is the Residence row.

    # ----- Haemoglobin, continuous per g/dL -----
    hgb_ref <- fixed(10)
    label("Centring haemoglobin concentration for the relative-hazard reference child (g/dL)")
    # Table 1, PK outcomes column: 'Hemoglobin level at diagnosis, g/dL,
    # median (range) 10.0 (5.6-15.9)'. A centring value is not given by
    # the paper; the cohort median is used so hr is 1 at a median child.

    e_hgb_haz <- log(0.93)
    label("Log hazard ratio for 28-day recurrent parasitaemia per 1 g/dL increase in haemoglobin (log scale; HR 0.93)")
    # Table 3, 'Hemoglobin level (g/dL)' row with footnote a marking it
    # continuous: adjusted HR 0.93 (95 percent CI 0.83-1.05), P = .25.

    # ----- Parasite density, continuous per natural-log unit -----
    para_ref <- fixed(17603)
    label("Centring parasite density for the relative-hazard reference child (parasites/uL)")
    # Table 1, PK outcomes column: 'Parasite density, parasites,
    # geometric mean no./uL (95% CI) 17 603 (13 762-22 516)'. Used as
    # the centring value so hr is 1 at a typical child.

    e_para_haz <- log(1.04)
    label("Log hazard ratio for 28-day recurrent parasitaemia per one natural-log-unit increase in parasite density (log scale; HR 1.04)")
    # Table 3, 'Parasite density (parasites/uL)' row with footnote a
    # marking it continuous: adjusted HR 1.04 (95 percent CI 0.95-1.14),
    # P = .35. Methods say only 'log transformed' without naming the
    # base; see the covariateData notes and the vignette.

    # ----- Underweight, reference WAZ >= -2 -----
    waz_cut <- fixed(-2)
    label("Weight-for-age z-score cut point defining the underweight stratum (z-score units)")
    # Methods, 'Association Analysis...': 'being underweight, based on a
    # weight-for-age z score with a -2 cutoff', WHO standards.

    e_waz_underweight_haz <- log(1.04)
    label("Log hazard ratio for 28-day recurrent parasitaemia, underweight (weight-for-age z below -2) vs not (log scale; HR 1.04)")
    # Table 3, 'Underweight (dichotomized) / Yes' row with footnote b:
    # adjusted HR 1.04 (95 percent CI 0.47-2.31), P = .93.

    # ==================================================================
    # Comparative continuous-exposure analysis, Results, 'Day 3 or 7
    # Lumefantrine Concentration and Clinical Outcomes'. These are TWO
    # SEPARATE multivariate Cox fits, each restricted to children NOT
    # receiving prophylaxis and each carrying one landmark on a
    # continuous natural-log concentration scale; the paper reports them
    # side by side to argue that the earlier landmark is the better
    # predictor. Their adjustment coefficients are not printed, so each
    # is encoded as an EXPOSURE-ONLY relative hazard between two
    # concentrations - a ratio in which any common adjustment terms
    # cancel exactly. That is the full identifiable content of what the
    # paper reports for these two fits.
    # ==================================================================

    conc_day3_ref <- fixed(2777)
    label("Reference day 3 lumefantrine concentration for the continuous comparative output (ng/mL)")
    # Results: 'the median concentration of lumefantrine on day 3
    # (n = 187 samples) was 2777 ng/mL (interquartile range
    # 1672-4760 ng/mL)'.

    conc_day7_ref <- fixed(216)
    label("Reference day 7 lumefantrine concentration for the continuous comparative output (ng/mL)")
    # Results: 'and on day 7 (n = 216 samples) was 216 ng/mL
    # (IQR 136-345 ng/mL)'.

    e_conc_lumefantrine_72h_haz <- log(0.51)
    label("Log hazard ratio for 28-day recurrent parasitaemia per one natural-log-unit increase in day 3 lumefantrine concentration, children not receiving prophylaxis (log scale; HR 0.51)")
    # Results: 'for children not receiving TMP-SMZ, each 1-unit increase
    # in natural log-transformed day 3 lumefantrine concentration was
    # associated with a 49% reduced hazard of 28-day recurrent malaria
    # (P = .002)'; 1 - 0.49 = 0.51, printed as 0.51 in the Discussion.

    e_conc_lumefantrine_168h_haz <- log(0.80)
    label("Log hazard ratio for 28-day recurrent parasitaemia per one natural-log-unit increase in day 7 lumefantrine concentration, children not receiving prophylaxis (log scale; HR 0.80)")
    # Results: 'each 1-unit increase in the natural log-transformed day
    # 7 lumefantrine concentration was associated with a 20% reduced
    # hazard of 28-day recurrent malaria (P = .002)'; 1 - 0.20 = 0.80,
    # printed as 0.80 in the Discussion. The smaller displacement from 1
    # than the day 3 coefficient above is the paper's evidence that
    # 'day 3 concentrations were stronger predictors of 28-day
    # recurrence than day 7 concentrations'.

    # No IIV: the repeated-episode correlation was handled by a robust
    # sandwich variance estimator, not by a subject-level random effect,
    # so there is no variance component to encode.

    # No residual error: a partial-likelihood Cox fit has no residual
    # variance component and none is reported.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Recover the paper's dichotomised and categorical covariates
    #    from the canonical continuous columns.
    # ------------------------------------------------------------------
    # Low day 7 exposure. Table 3's row heading is '<200 ng/mL', so the
    # comparison is strict and the reference group is at or above the cut.
    low <- (CONC_LUMEFANTRINE_168H < conc_cut)

    # Age group. The two indicators are mutually exclusive; when both are
    # 0 the episode is in the 6-to-under-12-month reference stratum.
    # Table 3's strata are '6-12', '12-18' and '>=18' months, so the
    # boundaries are taken as left-closed.
    age12to18 <- (AGE >= age_cut1) * (AGE < age_cut2)
    agege18 <- (AGE >= age_cut2)

    # Underweight. Methods give a '-2 cutoff' on weight-for-age z, read
    # as strictly below -2 so that the reference stratum is WAZ >= -2.
    underweight <- (WAZ < waz_cut)

    # ------------------------------------------------------------------
    # 2. Linear predictor of the adjusted Table 3 model and the relative
    #    hazard against the reference cohort: a 6-to-under-12-month,
    #    urban-resident, not-underweight child who is not receiving
    #    trimethoprim-sulfamethoxazole prophylaxis, whose day 7
    #    lumefantrine concentration is at or above 200 ng/mL and whose
    #    haemoglobin and parasite density sit at the cohort centring
    #    values. Under the proportional-hazards assumption hr is constant
    #    in time, so the episode hazard is h(t) = h0(t) * hr for any
    #    baseline h0(t) the user supplies.
    #
    #    The two low-exposure coefficients are gated on prophylaxis
    #    status so that exactly one of them is ever active, which is what
    #    reproduces Table 3's two within-stratum hazard ratios directly.
    #    max(PARA, 1) guards the logarithm against densities below one
    #    parasite/uL.
    # ------------------------------------------------------------------
    lhr <- e_conc_lumefantrine_168h_low_haz_nots * low * (1 - CONMED_TMPSMX) +
      e_conc_lumefantrine_168h_low_haz_ts * low * CONMED_TMPSMX +
      e_conmed_tmpsmx_haz * CONMED_TMPSMX +
      e_age_haz_12to18 * age12to18 +
      e_age_haz_ge18 * agege18 +
      e_residence_rural_haz * RESIDENCE_RURAL +
      e_hgb_haz * (HGB - hgb_ref) +
      e_para_haz * (log(max(PARA, 1)) - log(para_ref)) +
      e_waz_underweight_haz * underweight

    hr <- exp(lhr)

    # ------------------------------------------------------------------
    # 3. The two comparative continuous-exposure outputs, each from its
    #    own Cox fit in the no-prophylaxis subgroup. Each is a hazard
    #    RATIO between the episode's landmark concentration and the
    #    reference concentration for that landmark, so the unprinted
    #    adjustment terms of those fits cancel. They are alternative
    #    exposure parameterisations of the SAME endpoint as hr above, not
    #    additional multiplicative factors on it - do not multiply them
    #    together or into hr.
    # ------------------------------------------------------------------
    hr_c3 <- exp(e_conc_lumefantrine_72h_haz *
      (log(max(CONC_LUMEFANTRINE_72H, 1)) - log(conc_day3_ref)))

    hr_c7 <- exp(e_conc_lumefantrine_168h_haz *
      (log(max(CONC_LUMEFANTRINE_168H, 1)) - log(conc_day7_ref)))
  })
}
