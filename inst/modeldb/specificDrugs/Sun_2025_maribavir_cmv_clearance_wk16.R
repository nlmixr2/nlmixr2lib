Sun_2025_maribavir_cmv_clearance_wk16 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-EFFICACY model for the KEY ",
    "SECONDARY endpoint of the phase 3 AURORA study -- confirmed ",
    "clearance of plasma cytomegalovirus (CMV) DNA with no clinical ",
    "findings of tissue-invasive disease at week 8, maintained through ",
    "week 16 -- in allogeneic hematopoietic cell transplant (HCT) ",
    "recipients with a first asymptomatic CMV infection treated with ",
    "maribavir 400 mg orally twice daily (Sun 2025, n = 238, the AURORA ",
    "maribavir arm). THE EXPOSURE TERM IS NOT STATISTICALLY SIGNIFICANT, ",
    "and its point estimate is NEGATIVE (slope -0.0131 per 10 ug*h/mL, SE ",
    "0.0113, p = 0.247) where the primary endpoint's was positive -- the ",
    "two flat relationships do not even agree in sign, which is the ",
    "clearest available statement that neither carries signal. There is ",
    "no PK layer and no ODE: the exposure metric is supplied as a data ",
    "column, derived in the source analysis from the companion population ",
    "PK model packaged as Sun_2025_maribavir. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Fifteen companion exposure-response models in the ",
    "Sun_2025_maribavir_* family."
  )
  reference <- paste(
    "Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y, Bhattacharya I.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Maribavir",
    "in Transplant Recipients With First Episode or Refractory Cytomegalovirus.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1346-1356.",
    "doi:10.1002/psp4.70054.",
    "Logistic-regression coefficients in Supporting Information file s001, Table S3.",
    sep = " "
  )
  vignette <- "Sun_2025_maribavir"
  units <- list(
    time          = "n/a (static landmark exposure-response regression; no time dimension. The transplant-to-treatment interval T_HCT is carried in canonical hours and divided by 24 inside model() to recover the paper's days)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_MBV_SS covariate column)",
    concentration = "prob_cmv_clearance_wk16 (probability of confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8, maintained to week 16, 0-1; also logit_cmv_clearance_wk16)"
  )

  covariateData <- list(
    AUC_MBV_SS = list(
      description        = "Individual maribavir area under the plasma concentration-time curve at steady state, over the dosing interval on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual predictions from the companion maribavir population PK model.",
      units              = "ug*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sun 2025 estimates the coefficient PER INCREMENT OF 10 ug*h/mL, so",
        "model() divides this column by 10; a user supplying the raw AUC needs no",
        "rescaling of their own. Units are load-bearing -- the intercept absorbs",
        "the unit choice, so supplying mg*h/L is correct (it is the same number)",
        "but supplying ng*h/mL would shift the logit enormously. Derive the column",
        "from the companion popPK model, modellib('Sun_2025_maribavir'), which is",
        "the model Sun 2025 itself used. For calibration, the AURORA maribavir arm",
        "median steady-state AUC is about 170 ug*h/mL at 400 mg twice daily",
        "(Table 3), and Figure 3 draws its fitted curve across roughly",
        "50-600 ug*h/mL.",
        "In this model the exposure enters with slope -0.0131 per 10 ug*h/mL (odds",
        "ratio 0.987 (0.965-1.01), p = 0.247)."
      ),
      source_name        = "AUCss of maribavir (increment of 10 h.ug/mL)"
    ),
    TE_RESIST_MBV = list(
      description        = "Treatment-emergent cytomegalovirus mutation conferring resistance to maribavir.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no treatment-emergent maribavir-resistance mutation)",
      notes              = paste(
        "23 of 238 subjects (9.7%) positive; every subject was maribavir-",
        "susceptible at baseline (Table S2, 'CMV resistant at baseline: No",
        "238 (100)'), so this records an ON-TREATMENT event rather than a baseline",
        "stratifier. That makes the coefficient a descriptive association, not a",
        "quantity available for prognosis before treatment starts.",
        "Coefficient -3.51, odds ratio 0.0299 -- an even larger effect than on",
        "the primary endpoint, as expected for an endpoint that additionally",
        "requires the clearance to be MAINTAINED to week 16."
      ),
      source_name        = "TE CMV mutation conferring resistance to maribavir"
    ),
    CD8_PP65_MID = list(
      description        = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 0.5% to < 2% stratum.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (baseline CD8+CD69+ pp65 stimulation < 0.5%, the lowest-immunity stratum)",
      notes              = paste(
        "24 of 238 subjects (10.1%). One of three indicators decomposing the",
        "four-level baseline CMV-specific cell-mediated immunity assay; the other",
        "two are CD8_PP65_HI and CD8_PP65_NR, and all three are 0 for the < 0.5%",
        "reference. A higher percentage means more surviving CMV-directed",
        "cytotoxic T-cell immunity.",
        "Coefficient +0.623, not significant."
      ),
      source_name        = "Baseline CD8+ CD69+ pp65 stimulation group, >= 0.5% to < 2%"
    ),
    CD8_PP65_HI = list(
      description        = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 2% stratum.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (baseline CD8+CD69+ pp65 stimulation < 0.5%, the lowest-immunity stratum)",
      notes              = paste(
        "71 of 238 subjects (29.8%). The highest-immunity stratum of the assay.",
        "Coefficient +1.01, significant. The set is monotone in this endpoint",
        "(+0.623 then +1.01 with rising immunity), unlike in the primary",
        "endpoint model, which is the more biologically interpretable ordering."
      ),
      source_name        = "Baseline CD8+ CD69+ pp65 stimulation group, >= 2%"
    ),
    CD8_PP65_NR = list(
      description        = "Baseline CD8+CD69+ pp65-stimulated T-cell assay not reported.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (baseline CD8+CD69+ pp65 stimulation < 0.5%, the lowest-immunity stratum)",
      notes              = paste(
        "46 of 238 subjects (19.3%). This is a FITTED level with its own estimated",
        "coefficient, not a missing-data code: a subject whose assay was not",
        "reported must be given CD8_PP65_NR = 1 rather than being imputed into the",
        "< 0.5% reference, which would silently move a fifth of the cohort.",
        "Coefficient +1.11, significant."
      ),
      source_name        = "Baseline CD8+ CD69+ pp65 stimulation group, Not reported"
    ),
    T_HCT = list(
      description        = "Time from the onset of the current haematopoietic cell transplant to the start of maribavir treatment.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Canonical units are hours, so model() divides by 24 to recover the DAYS",
        "in which Sun 2025 estimates the coefficient -- the same convention the",
        "T_FIRSTDOSE register entry prescribes. A per-subject scalar fixed before",
        "the analysis window opens, not a clock that runs during it. It is a proxy",
        "for how far into post-transplant immune reconstitution the subject was",
        "when maribavir began.",
        "Coefficient +0.0489 PER DAY, so model() divides the canonical hours by",
        "24. Subjects treated later after transplant were more likely to achieve",
        "sustained clearance, consistent with greater elapsed immune",
        "reconstitution. Retained only in this endpoint, not in the primary."
      ),
      source_name        = "Days from onset of current HCT to treatment"
    ),
    REGION_EUROPE = list(
      description        = "Europe enrolling-region indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (North America, the reference region, when REGION_ASIAPACIFIC is also 0)",
      notes              = paste(
        "138 of 238 subjects (58.0%). Paired with REGION_ASIAPACIFIC against a",
        "NORTH AMERICA reference -- both indicators 0 -- rather than the",
        "REGION_* family's more usual US reference. Enrolling region here is a",
        "proxy for regional differences in transplant and CMV-management practice,",
        "not for ancestry: Table S2 records 138 European enrollments but only 31",
        "subjects of Asian race against 41 Asia Pacific enrollments, so region and",
        "race disagree within this cohort.",
        "Coefficient -0.641, not significant here (it was significant on the",
        "primary endpoint)."
      ),
      source_name        = "Enrolling region, Europe"
    ),
    REGION_ASIAPACIFIC = list(
      description        = "Asia Pacific enrolling-region indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (North America, the reference region, when REGION_EUROPE is also 0)",
      notes              = paste(
        "41 of 238 subjects (17.2%). Paired with REGION_EUROPE against a NORTH",
        "AMERICA reference. Broader than REGION_JAPAN or REGION_EASTASIA: Asia",
        "Pacific as a trial-operations region conventionally spans East Asia,",
        "South-East Asia and Oceania.",
        "Coefficient -1.49 against the North America reference, significant."
      ),
      source_name        = "Enrolling region, Asia Pacific"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 238L,
    n_studies      = 1L,
    n_observations = "238 binary the key secondary endpoint records, one per patient",
    age_range      = "12 to <18 years: 1 (0.4%); 18 to <45: 47 (19.7%); 45 to <65: 142 (59.7%); >=65: 48 (20.2%) (Table S2)",
    weight_range   = "Not reported for the AURORA arm; Table 1 reports 'NR' for the AURORA weights specifically. The pooled PK analysis population spans 36.1-141 kg.",
    sex_female_pct = 45.4,
    race_ethnicity = c(Caucasian = 79.8, Asian = 13.0, Black = 3.4, Other = 2.9, Missing = 0.8),
    disease_state  = "Allogeneic HCT recipients (238/238, 100%) with FIRST asymptomatic cytomegalovirus infection after transplant; all were maribavir-susceptible at baseline. Baseline CMV DNA: very low 71, low 127, high 39, missing 1. CMV serostatus D+/R+ 123 (51.7%), D-/R+ 83 (34.9%), D+/R- 18 (7.6%), D-/R- 8 (3.4%). Reason for transplant: acute myeloid leukaemia 88 (37.0%), myelodysplastic syndrome 37 (15.5%), acute lymphocytic leukaemia 21 (8.8%), non-Hodgkin lymphoma 19 (8.0%), other 73 (30.7%). Conditioning: reduced-intensity 116 (48.7%), myeloablative 86 (36.1%), non-myeloablative 32 (13.4%).",
    dose_range     = "Maribavir 400 mg orally twice daily, the AURORA randomized dose",
    regions        = "North America 59 (24.8%), Europe 138 (58.0%), Asia Pacific 41 (17.2%) (Table S2)",
    notes          = paste0(
      "This is the exposure-response analysis population: the maribavir arm of ",
      "the phase 3 AURORA study, a subset of the 930-subject population PK ",
      "analysis population. Individual exposures were derived from the ",
      "companion population PK model, packaged as Sun_2025_maribavir. ",
      "Sun 2025 selected safety endpoints for formal exposure-response analysis by incidence (TEAEs occurring in more than 10% of participants, plus the most severe TEAEs)."
    )
  )

  ini({
    # ================================================================
    # Sun 2025 Table S3 (Supporting Information file s001),
    # column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16'.
    #
    #   logit(p) = intercept + slope * (AUC_MBV_SS / 10)
    #              + risk-factor terms
    #
    # Every value below is the source's 'Estimate (SE)' column. Each is
    # confirmed a second time by its printed odds ratio, which the source
    # defines as the exponentiated estimate; because both columns are
    # rounded, the check is the INTERVAL one (does a real b exist with
    # round(b, 3) == Estimate and round(exp(b), 3) == OR) rather than
    # exp(printed) == printed. The checks are re-run mechanically in the
    # validation vignette. Covariates are not centred except where noted,
    # so the intercept is an anchor at zero exposure, not a reference-
    # patient probability.
    # ================================================================

    # ----- Logit intercept -----
    logit_ref <- 0.794; label("Logit of the the key secondary endpoint probability at zero maribavir exposure with all risk factors at their reference level (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Intercept 0.794 (SE 0.493), p = 0.10742

    # ----- Exposure effect -----
    e_auc_logit <- -0.0131; label("Log-odds of the key secondary endpoint per 10 ug*h/mL increase in maribavir steady-state AUC (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Estimate -0.0131 (SE 0.0113), odds ratio 0.987 (0.965-1.01), p = 0.247

    # ----- Risk-factor effects on the logit -----
    e_te_resist_mbv_logit <- -3.51; label("Log-odds of the key secondary endpoint for a treatment-emergent maribavir-resistance mutation (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', TE CMV mutation conferring resistance to maribavir Estimate -3.51 (SE 1.05), odds ratio 0.0299 (0.00381-0.235), p < 0.001
    e_cd8_pp65_mid_logit <- 0.623; label("Log-odds of the key secondary endpoint for the middle baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Baseline CD8+ CD69+ pp65 stimulation group, >= 0.5% to < 2% Estimate 0.623 (SE 0.521), odds ratio 1.86 (0.671-5.18), p = 0.232
    e_cd8_pp65_hi_logit <- 1.01; label("Log-odds of the key secondary endpoint for the highest baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Baseline CD8+ CD69+ pp65 stimulation group, >= 2% Estimate 1.01 (SE 0.367), odds ratio 2.74 (1.33-5.62), p = 0.00192
    e_cd8_pp65_nr_logit <- 1.11; label("Log-odds of the key secondary endpoint for an unreported baseline CD8+CD69+ pp65 assay versus the lowest stratum (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Baseline CD8+ CD69+ pp65 stimulation group, Not reported Estimate 1.11 (SE 0.463), odds ratio 3.04 (1.23-7.53), p = 0.0163
    e_t_hct_logit <- 0.0489; label("Log-odds of the key secondary endpoint for each additional DAY between the current HCT and the start of treatment (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Days from onset of current HCT to treatment Estimate 0.0489 (SE 0.0199), odds ratio 1.05 (1.01-1.09), p = 0.0140
    e_region_europe_logit <- -0.641; label("Log-odds of the key secondary endpoint for enrolment in Europe versus North America (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Enrolling region, Europe Estimate -0.641 (SE 0.412), odds ratio 0.527 (0.235-1.18), p = 0.120
    e_region_asiapacific_logit <- -1.49; label("Log-odds of the key secondary endpoint for enrolment in Asia Pacific versus North America (unitless logit)")  # Sun 2025 Table S3 (Supporting Information file s001), column 'Confirmed CMV clearance with no clinical findings of tissue-invasive disease at week 8 followed by maintenance to week 16', Enrolling region, Asia Pacific Estimate -1.49 (SE 0.508), odds ratio 0.225 (0.0832-0.609), p = 0.00334

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli likelihood
    # has no sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_cmv_clearance_wk16 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Sun 2025 estimates every coefficient per 10 ug*h/mL increment, so the
    # exposure column is divided by 10 here rather than in the data.
    auc10 <- AUC_MBV_SS / 10

    # T_HCT is carried in canonical hours; Sun 2025 estimates its
    # coefficient per DAY ('Days from onset of current HCT to treatment').

    # ----- Linear predictor -----
    logit_cmv_clearance_wk16 <- logit_ref +
      e_auc_logit * auc10 +
      e_te_resist_mbv_logit * TE_RESIST_MBV +
      e_cd8_pp65_mid_logit * CD8_PP65_MID +
      e_cd8_pp65_hi_logit * CD8_PP65_HI +
      e_cd8_pp65_nr_logit * CD8_PP65_NR +
      e_t_hct_logit * (T_HCT / 24) +
      e_region_europe_logit * REGION_EUROPE +
      e_region_asiapacific_logit * REGION_ASIAPACIFIC

    prob_cmv_clearance_wk16 <- expit(logit_cmv_clearance_wk16)

    # ----- Observation -----
    prob_cmv_clearance_wk16 ~ add(addSd_prob_cmv_clearance_wk16)
  })
}
