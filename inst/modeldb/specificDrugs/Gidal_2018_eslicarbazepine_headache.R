Gidal_2018_eslicarbazepine_headache <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for the probability of ",
    "treatment-emergent HEADACHE in adults with focal-onset seizures ",
    "taking adjunctive eslicarbazepine acetate (ESL) (Gidal 2018, ",
    "n = 1,152 patients pooled from the phase 3 trials 2093-301, ",
    "2093-302 and 2093-304). The linear predictor is ",
    "-2.43 + 1.42*I(start 400 mg) + 2.64*I(start 800 mg) ",
    "- 0.00000632*AUC0-24 - 0.0231*(WT - 70) (Gidal 2018 Appendix S1 ",
    "Eq. E-4, Table S-4). This is the most parsimonious of the paper's ",
    "three TEAE models: region, baseline carbamazepine use and sex were ",
    "all screened and NOT retained, and only the first-week starting ",
    "dose, eslicarbazepine AUC0-24 and body weight remain. Concomitant ",
    "lamotrigine was reported as a significant predictor of headache in ",
    "the main text (p < 0.05) but does not appear in the final Table S-4 ",
    "coefficient set, so it is documented as screened-not-retained ",
    "rather than encoded; see the vignette Errata. As in the companion ",
    "dizziness model the AUC0-24 coefficient is NEGATIVE, reflecting ",
    "that only the FIRST occurrence of the event was modelled and that ",
    "first occurrences cluster in the low-exposure titration period. ",
    "There is no PK layer and no ODE: exposure enters as the static ",
    "per-patient column AUC_ESL, an empirical-Bayes prediction from ",
    "modellib('Gidal_2018_eslicarbazepine'). No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Hosmer-Lemeshow chi-squared 14.30 on 8 df (p = 0.0743), C statistic ",
    "0.6630, minimum objective function 737.594 -- the paper describes ",
    "the fit as adequate but only somewhat predictive."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equation are in Appendix S1 (supporting",
    "information), Table S-4 and Equation E-4.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-safety regression over the 14-week double-blind period; no time dimension)",
    dosing = "n/a (no dose events; the first-week starting dose enters as the covariate DOSE_ESL_MGD)",
    concentration = "prob_headache (probability of treatment-emergent headache, 0-1; also logit_headache)"
  )

  covariateData <- list(
    AUC_ESL = list(
      description = paste(
        "Individual predicted eslicarbazepine area under the plasma",
        "concentration-time curve over the 24-hour dosing interval."
      ),
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; compute it as dose / CL/F with",
        "modellib('Gidal_2018_eslicarbazepine'). At steady state under",
        "once-daily dosing AUC0-24 = dose / (CL/F), so the reference",
        "subject (CL/F 2.43 L/h) has 164,600 ng*h/mL on 400 mg and",
        "329,200 ng*h/mL on 800 mg. Set to 0 for placebo patients, which",
        "makes the exposure term vanish and leaves the intercept as the",
        "placebo logit. NOT centred and NOT scaled. The coefficient is",
        "NEGATIVE; see the description."
      ),
      source_name = "AUC (Eq. E-4), eslicarbazepine AUC0-24"
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
        "randomised maintenance dose. It is the strongest predictor of",
        "headache, as it is for dizziness and somnolence."
      ),
      source_name = "400 mg_i and 800 mg_i (Eq. E-4)"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters linearly, centred at the population median of 70 kg, so",
        "the coefficient is a log-odds change per kg. The slope (-0.0231)",
        "is within 7% of the dizziness and somnolence slopes (-0.0247",
        "each), so the weight effect is essentially common to all three",
        "TEAE models."
      ),
      source_name = "WTKG (Eq. E-4)"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Gidal 2018 Results states that 'concomitant use of lamotrigine",
        "was predicted to increase the risk of dizziness (p < 0.01) and",
        "headache (p < 0.05)', but Table S-4 prints only five coefficient",
        "rows -- intercept, the two starting-dose shifts, the AUC0-24",
        "slope and the weight slope -- and Eq. E-4 contains no lamotrigine",
        "term. The main text's description of the final model likewise",
        "lists only 'a linear function of the eslicarbazepine AUC0-24, and",
        "weight'. No point estimate exists in any available source, so the term",
        "cannot be encoded; the printed final model is taken as",
        "authoritative over the narrative sentence. See the vignette",
        "Errata. 14.0% of the pooled population took lamotrigine."
      )
    ),
    REGION_NORTHAMERICA = list(
      description = "North American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and explicitly not retained: 'Region and concomitant use",
        "of CBZ were not statistically significant predictors of the",
        "probability of headache' (Gidal 2018 Results). This is a",
        "published null result, not a reporting gap, so the correct",
        "encoding is the omission of the term rather than a fixed(0)",
        "coefficient. Retained in the companion dizziness model."
      )
    ),
    REGION_LATINAMERICA = list(
      description = "Latin American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained for headache; see REGION_NORTHAMERICA."
    ),
    REGION_ROW = list(
      description = "Rest-of-World study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained for headache; see REGION_NORTHAMERICA."
    ),
    CONMED_CBZ = list(
      description = "Carbamazepine use during the baseline period; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and explicitly not retained for headache (Gidal 2018",
        "Results). Retained with opposite signs in the companion dizziness",
        "(+0.626) and somnolence (-0.459) models."
      )
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and not retained for headache; the paper reports a",
        "female excess risk only for dizziness and somnolence (p < 0.05)."
      )
    ),
    CONMED_LEV = list(
      description = "Concomitant levetiracetam; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "'Concomitant use of levetiracetam and valproate were not",
        "statistically significant predictors of the probability of the",
        "TEAEs evaluated' (Gidal 2018 Results) -- screened and not",
        "retained for any of the three TEAE endpoints."
      )
    ),
    CONMED_VPA = list(
      description = "Concomitant valproate; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained for any of the three TEAE endpoints; see CONMED_LEV."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1152L,
    n_studies = 3L,
    n_observations = "1,152 binary first-occurrence headache records (one per patient; landmark analysis, no repeated measures)",
    age_median = "37 years",
    weight_median = "70 kg (the centring value used by Eq. E-4)",
    race_ethnicity = c(Caucasian = 80.0),
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo, after a 2-week titration from a 400 or 800 mg starting dose",
    regions = "Europe, North America, Latin America and Rest of World (region was screened and not retained for this endpoint)",
    co_medication = "48% took carbamazepine during the baseline period",
    notes = paste(
      "Safety analysis set: 306 patients from study 2093-301, 307 from",
      "2093-302 and 539 from 2093-304. Both linear and power exposure",
      "models were screened; AUC0-24 was the significant metric for",
      "headache, as for dizziness. The C statistic of 0.66 is the lowest",
      "of the paper's three TEAE models."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-4 and Equation E-4:
    #
    #   logit(p) = -2.43 + 1.42*I(400 mg) + 2.64*I(800 mg)
    #              - 0.00000632*AUC - 0.0231*(WTKG - 70)
    #
    # Five rows, reported once in Table S-4 and once in Eq. E-4; both
    # agree exactly. Only WT is centred (at 70 kg), so the intercept is
    # the logit for a 70 kg patient on placebo.
    # ==================================================================

    logit_ref <- -2.43 ; label("Logit of the probability of headache for a 70 kg patient on placebo (unitless logit)")  # Table S-4, 'Logit intercept (placebo effect)' -2.43, 7.4% SEM; Eq. E-4. Corresponds to an 8.1% placebo headache probability
    e_dose_esl_mgd400_logit <- 1.42 ; label("Log-odds shift for a 400 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-4, 'Additive shift in logit for 400 mg starting dose' 1.42, 23.3% SEM; Eq. E-4
    e_dose_esl_mgd800_logit <- 2.64 ; label("Log-odds shift for an 800 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-4, 'Additive shift in logit for 800 mg starting dose' 2.64, 15.9% SEM; Eq. E-4
    e_auc_esl_logit <- -0.00000632 ; label("Log-odds of headache per 1 ng*h/mL increase in eslicarbazepine AUC0-24 (unitless logit)")  # Table S-4, 'Slope for eslicarbazepine AUC0-24' -0.00000632, 20.9% SEM; Eq. E-4. Negative by design; see the description and vignette Errata item 3
    e_wt_logit <- -0.0231 ; label("Log-odds of headache per kg of body weight above 70 kg (unitless logit)")  # Table S-4, 'Slope for weight' -0.0231, 27.9% SEM; Eq. E-4

    # ----- No between-subject variability, no residual error -----
    # Bernoulli likelihood: the source estimates no sigma and no random
    # effects. The tiny fixed additive residual exists only so rxode2 has
    # an error model to attach to the typical-value probability.
    addSd_prob_headache <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Starting-dose indicators derived from the dose-level column ---
    dose400 <- (DOSE_ESL_MGD == 400)
    dose800 <- (DOSE_ESL_MGD == 800)

    # ----- Linear predictor (Gidal 2018 Eq. E-4) -----
    logit_headache <- logit_ref +
      e_dose_esl_mgd400_logit * dose400 +
      e_dose_esl_mgd800_logit * dose800 +
      e_auc_esl_logit * AUC_ESL +
      e_wt_logit * (WT - 70)

    prob_headache <- expit(logit_headache)

    # ----- Observation -----
    prob_headache ~ add(addSd_prob_headache)
  })
}
