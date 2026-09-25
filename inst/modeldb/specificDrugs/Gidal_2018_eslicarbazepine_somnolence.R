Gidal_2018_eslicarbazepine_somnolence <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for the probability of ",
    "treatment-emergent SOMNOLENCE in adults with focal-onset seizures ",
    "taking adjunctive eslicarbazepine acetate (ESL) (Gidal 2018, ",
    "n = 1,152 patients pooled from the phase 3 trials 2093-301, ",
    "2093-302 and 2093-304). The linear predictor is ",
    "-2.62 + 2.47*I(start 400 mg) + 4.22*I(start 800 mg) ",
    "- 0.000178*Cmax - 0.0247*(WT - 70) + 0.661*LatinAmerica ",
    "+ 0.464*female - 0.459*baselineCBZ (Gidal 2018 Appendix S1 ",
    "Eq. E-5, Table S-5). This is the only one of the paper's three TEAE ",
    "models driven by PEAK concentration rather than by AUC0-24: both ",
    "metrics were screened for every endpoint and Cmax was the ",
    "statistically significant predictor here. Two further contrasts ",
    "with the companion dizziness model: only Latin America separates ",
    "from the other regions (North America and Rest of World were not ",
    "retained), and baseline carbamazepine use is PROTECTIVE here ",
    "(-0.459) where it was a risk factor for dizziness (+0.626) -- Gidal ",
    "2018 suggests baseline carbamazepine users may have developed ",
    "tolerance to the sedative effects of the drug class. As in the ",
    "other two TEAE models the exposure coefficient is NEGATIVE because ",
    "only the FIRST occurrence of the event was modelled and first ",
    "occurrences cluster in the low-exposure titration period. There is ",
    "no PK layer and no ODE: exposure enters as the static per-patient ",
    "column CMAX, an empirical-Bayes prediction from ",
    "modellib('Gidal_2018_eslicarbazepine'). No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Hosmer-Lemeshow chi-squared 11.06 on 10 df (p = 0.3531), C ",
    "statistic 0.7548, minimum objective function 726.433."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equation are in Appendix S1 (supporting",
    "information), Table S-5 and Equation E-5.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-safety regression over the 14-week double-blind period; no time dimension)",
    dosing = "n/a (no dose events; the first-week starting dose enters as the covariate DOSE_ESL_MGD)",
    concentration = "prob_somnolence (probability of treatment-emergent somnolence, 0-1; also logit_somnolence)"
  )

  covariateData <- list(
    CMAX = list(
      description = paste(
        "Individual predicted maximum eslicarbazepine plasma concentration",
        "over the 24-hour dosing interval."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; simulate it with",
        "modellib('Gidal_2018_eslicarbazepine') and take the",
        "steady-state peak of the Cc profile. Unlike AUC0-24 there is no",
        "closed form, because Cmax depends on ka and V/F as well as on",
        "CL/F. The ng/mL scale is the one Table S-2 uses throughout",
        "(its residual-error footnotes quote predicted concentrations of",
        "740-39,200 ng/mL). Set to 0 for placebo patients, which makes",
        "the exposure term vanish and leaves the intercept as the placebo",
        "logit. NOT centred and NOT scaled. The coefficient is NEGATIVE;",
        "see the description. Cmax rather than AUC0-24 is the retained",
        "metric for this endpoint only."
      ),
      source_name = "C_max (Eq. E-5), eslicarbazepine Cmax"
    ),
    DOSE_ESL_MGD = list(
      description = paste(
        "Eslicarbazepine acetate STARTING daily dose during the first week",
        "of the titration period: 0 (placebo), 400 or 800 mg/day."
      ),
      units = "mg/day",
      type = "continuous",
      reference_category = "0 mg/day (placebo)",
      notes = paste(
        "Carried as a dose level rather than as two separate indicator",
        "columns; model() derives the paper's two indicator variables as",
        "(DOSE_ESL_MGD == 400) and (DOSE_ESL_MGD == 800), so only those",
        "three values are meaningful. This is the FIRST-WEEK dose, not the",
        "randomised maintenance dose. Somnolence carries the largest",
        "starting-dose shifts of the paper's three TEAE models (2.47 and",
        "4.22 on the logit)."
      ),
      source_name = "400 mg_i and 800 mg_i (Eq. E-5)"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters linearly, centred at the population median of 70 kg, so",
        "the coefficient is a log-odds change per kg. The slope is",
        "numerically identical to the dizziness model's (-0.0247)."
      ),
      source_name = "WTKG (Eq. E-5)"
    ),
    REGION_LATINAMERICA = list(
      description = "Latin American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Europe, North America or Rest of World pooled)",
      notes = paste(
        "The ONLY region indicator retained for somnolence: Gidal 2018",
        "Results reports 'the risk of somnolence was predicted to be",
        "higher in patients from Latin America than those from Europe,",
        "North America, and Rest of World (p < 0.05)'. The reference group",
        "is therefore the pooled other three regions, which is a different",
        "reference from the companion dizziness model (where Europe alone",
        "is the reference and the other three regions each carry their own",
        "shift). Do not transfer a region encoding between the two models."
      ),
      source_name = "LA_i (Eq. E-5)"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Same orientation as Gidal 2018's SEXF_i indicator, so no value",
        "transformation is needed. Women are predicted to be at higher",
        "risk of somnolence and of dizziness (p < 0.05); the coefficient",
        "here (0.464) is close to the dizziness one (0.486)."
      ),
      source_name = "SEXF_i (Eq. E-5)"
    ),
    CONMED_CBZ = list(
      description = "Carbamazepine use during the baseline period; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (took other antiepileptic drugs at baseline)",
      notes = paste(
        "BASELINE use, i.e. during the 8-week pre-randomisation period,",
        "not concomitant use during treatment. 48% of the safety analysis",
        "population. The coefficient is NEGATIVE here (-0.459), the",
        "opposite sign from the same covariate in the companion dizziness",
        "model (+0.626). Gidal 2018 Discussion states the reason for the",
        "protective direction is unknown but could reflect tolerance to",
        "the sedative effects of the voltage-gated sodium channel",
        "modulator class developed during the baseline period. The sign",
        "flip is a published result, not a transcription error."
      ),
      source_name = "BCARB_i (Eq. E-5)"
    )
  )

  covariatesDataExcluded <- list(
    AUC_ESL = list(
      description = "Individual predicted eslicarbazepine AUC over the 24-hour dosing interval.",
      units = "ng*h/mL",
      type = "continuous",
      notes = paste(
        "Screened as the alternative exposure metric and not retained:",
        "Gidal 2018 Results states that 'eslicarbazepine AUC0-24 was found",
        "to be a statistically significant predictor of the probability of",
        "dizziness and headache, while Cmax was a statistically",
        "significant predictor of the probability of somnolence'. No",
        "AUC0-24 coefficient is printed for somnolence, so the correct",
        "encoding is the omission of the term rather than a fixed(0)",
        "coefficient. Retained in the companion dizziness and headache",
        "models."
      )
    ),
    REGION_NORTHAMERICA = list(
      description = "North American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and not retained for somnolence -- only Latin America",
        "separated from the other regions. Retained in the companion",
        "dizziness model. Pooled into the reference group here."
      )
    ),
    REGION_ROW = list(
      description = "Rest-of-World study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained for somnolence; pooled into the reference group. See REGION_NORTHAMERICA."
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Explicitly screened and not retained: 'Concomitant use of",
        "lamotrigine was predicted to increase the risk of dizziness",
        "(p < 0.01) and headache (p < 0.05), but not that of somnolence'",
        "(Gidal 2018 Results). A published null result rather than a",
        "reporting gap."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1152L,
    n_studies = 3L,
    n_observations = "1,152 binary first-occurrence somnolence records (one per patient; landmark analysis, no repeated measures)",
    age_median = "37 years",
    weight_median = "70 kg (the centring value used by Eq. E-5)",
    race_ethnicity = c(Caucasian = 80.0),
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo, after a 2-week titration from a 400 or 800 mg starting dose",
    regions = "Latin America versus the pooled Europe / North America / Rest of World reference",
    co_medication = "48% took carbamazepine during the baseline period",
    notes = paste(
      "Safety analysis set: 306 patients from study 2093-301, 307 from",
      "2093-302 and 539 from 2093-304. Both linear and power exposure",
      "models and both AUC0-24 and Cmax were screened for every endpoint;",
      "Cmax was retained only here. The Hosmer-Lemeshow statistic is",
      "quoted against 10 degrees of freedom in Appendix S1 while the",
      "dizziness and headache statistics use 8; the difference is not",
      "explained in the source."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-5 and Equation E-5:
    #
    #   logit(p) = -2.62 + 2.47*I(400 mg) + 4.22*I(800 mg)
    #              - 0.000178*Cmax - 0.0247*(WTKG - 70)
    #              + 0.661*LA + 0.464*SEXF - 0.459*BCARB
    #
    # Eight rows, reported once in Table S-5 and once in Eq. E-5; both
    # agree exactly. Only WT is centred (at 70 kg); Cmax is neither
    # centred nor scaled, so the intercept is the logit for a 70 kg man
    # on placebo outside Latin America with no baseline carbamazepine.
    # ==================================================================

    logit_ref <- -2.62 ; label("Logit of the probability of somnolence for a 70 kg man on placebo outside Latin America, with no baseline carbamazepine (unitless logit)")  # Table S-5, 'Logit intercept (placebo effect)' -2.62, 9.1% SEM; Eq. E-5. Corresponds to a 6.8% placebo somnolence probability
    e_dose_esl_mgd400_logit <- 2.47 ; label("Log-odds shift for a 400 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-5, 'Additive shift in logit for 400 mg starting dose' 2.47, 14.6% SEM; Eq. E-5
    e_dose_esl_mgd800_logit <- 4.22 ; label("Log-odds shift for an 800 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-5, 'Additive shift in logit for 800 mg starting dose' 4.22, 12.1% SEM; Eq. E-5
    e_cmax_logit <- -0.000178 ; label("Log-odds of somnolence per 1 ng/mL increase in eslicarbazepine Cmax (unitless logit)")  # Table S-5, 'Slope for eslicarbazepine Cmax' -0.000178, 14.4% SEM; Eq. E-5. Negative by design; see the description and vignette Errata item 3
    e_wt_logit <- -0.0247 ; label("Log-odds of somnolence per kg of body weight above 70 kg (unitless logit)")  # Table S-5, 'Slope for weight' -0.0247, 28.2% SEM; Eq. E-5
    e_region_latinamerica_logit <- 0.661 ; label("Log-odds shift for a Latin American study site versus the pooled Europe / North America / Rest of World reference (unitless logit)")  # Table S-5, 'Additive shift for Latin America' 0.661, 32.2% SEM; Eq. E-5
    e_sexf_logit <- 0.464 ; label("Log-odds shift for female versus male sex (unitless logit)")  # Table S-5, 'Additive shift for females' 0.464, 45.9% SEM; Eq. E-5
    e_conmed_cbz_logit <- -0.459 ; label("Log-odds shift for baseline carbamazepine use (unitless logit)")  # Table S-5, 'Additive shift for baseline carbamazepine use' -0.459, 43.8% SEM; Eq. E-5. Protective direction, opposite to the same covariate in the dizziness model

    # ----- No between-subject variability, no residual error -----
    # Bernoulli likelihood: the source estimates no sigma and no random
    # effects. The tiny fixed additive residual exists only so rxode2 has
    # an error model to attach to the typical-value probability.
    addSd_prob_somnolence <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Starting-dose indicators derived from the dose-level column ---
    dose400 <- (DOSE_ESL_MGD == 400)
    dose800 <- (DOSE_ESL_MGD == 800)

    # ----- Linear predictor (Gidal 2018 Eq. E-5) -----
    logit_somnolence <- logit_ref +
      e_dose_esl_mgd400_logit * dose400 +
      e_dose_esl_mgd800_logit * dose800 +
      e_cmax_logit * CMAX +
      e_wt_logit * (WT - 70) +
      e_region_latinamerica_logit * REGION_LATINAMERICA +
      e_sexf_logit * SEXF +
      e_conmed_cbz_logit * CONMED_CBZ

    prob_somnolence <- expit(logit_somnolence)

    # ----- Observation -----
    prob_somnolence ~ add(addSd_prob_somnolence)
  })
}
