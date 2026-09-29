Gidal_2018_eslicarbazepine_dizziness <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for the probability of ",
    "treatment-emergent DIZZINESS in adults with focal-onset seizures ",
    "taking adjunctive eslicarbazepine acetate (ESL) (Gidal 2018, ",
    "n = 1,152 patients pooled from the phase 3 trials 2093-301, ",
    "2093-302 and 2093-304). The linear predictor is ",
    "-3.67 + 1.97*I(start 400 mg) + 3.62*I(start 800 mg) ",
    "- 0.00000835*AUC0-24 - 0.0247*(WT - 70) + 1.04*NorthAmerica ",
    "+ 1.13*LatinAmerica + 0.668*RestOfWorld + 0.626*baselineCBZ ",
    "+ 0.589*lamotrigine + 0.486*female (Gidal 2018 Appendix S1 Eq. E-3, ",
    "Table S-3). The dominant predictor is the FIRST-WEEK STARTING DOSE, ",
    "not the maintenance dose. The eslicarbazepine AUC0-24 coefficient is ",
    "NEGATIVE, so higher exposure predicts LESS dizziness once starting ",
    "dose is accounted for; Gidal 2018 calls this finding unexpected and ",
    "attributes it to the models counting only the FIRST occurrence of an ",
    "event, which typically falls in the low-exposure 2-week titration ",
    "period. Europe is the region reference. There is no PK layer and no ",
    "ODE: exposure enters as the static per-patient column AUC_ESL, an ",
    "empirical-Bayes prediction from the companion model ",
    "modellib('Gidal_2018_eslicarbazepine'). No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Hosmer-Lemeshow chi-squared 14.86 on 8 df (p = 0.062), C statistic ",
    "0.7947, minimum objective function 794.541."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equation are in Appendix S1 (supporting",
    "information), Table S-3 and Equation E-3.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-safety regression over the 14-week double-blind period; no time dimension)",
    dosing = "n/a (no dose events; the first-week starting dose enters as the covariate DOSE_ESL_MGD)",
    concentration = "prob_dizziness (probability of treatment-emergent dizziness, 0-1; also logit_dizziness)"
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
        "329,200 ng*h/mL on 800 mg. The ng*h/mL scale is fixed by the",
        "companion serum-sodium model, whose Table S-6 prints its slope",
        "units explicitly as (mmol/L)/(ng x h/mL) and whose main-text",
        "worked example (a 400 mg dose increase raises AUC0-24 by",
        "165 ug*h/mL and lowers sodium by 0.68 mEq/L) reproduces only on",
        "that scale. Set to 0 for placebo patients, which makes the",
        "exposure term vanish and leaves the intercept as the placebo",
        "logit. NOT centred and NOT scaled.",
        "The coefficient is NEGATIVE; see the description."
      ),
      source_name = "AUC (Eq. E-3), eslicarbazepine AUC0-24"
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
        "randomised maintenance dose -- patients randomised to 800 or",
        "1,200 mg titrated up from one of these two starting doses over a",
        "2-week titration period. Gidal 2018 identifies it as the single",
        "strongest predictor of each of the three analysed TEAEs, and the",
        "practical conclusion of the paper (Table 1) is that starting at",
        "400 mg rather than 800 mg gives a better risk-benefit profile."
      ),
      source_name = "400 mg_i and 800 mg_i (Eq. E-3)"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters linearly, centred at the population median of 70 kg, so",
        "the coefficient is a log-odds change per kg. The negative sign",
        "means heavier patients are predicted to be at lower risk",
        "(p < 0.001 in the paper's Results). Gidal 2018 Discussion links",
        "this to the rise of CL/F and V/F with weight, though the",
        "exposure term is already in the model."
      ),
      source_name = "WTKG (Eq. E-3)"
    ),
    REGION_NORTHAMERICA = list(
      description = "North American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with REGION_LATINAMERICA and REGION_ROW also 0 this selects the EUROPE reference group",
      notes = paste(
        "One of three region indicators against a Europe reference. The",
        "three are mutually exclusive: exactly one of",
        "REGION_NORTHAMERICA / REGION_LATINAMERICA / REGION_ROW is 1, or",
        "all three are 0 for a European patient. Study 2093-304 was the",
        "North American trial; studies 301 and 302 were European. Gidal",
        "2018 reports the risk of dizziness as higher in North America,",
        "Latin America and Rest of World than in Europe (p < 0.05) and",
        "attributes it to demographic and clinical differences between",
        "regions rather than to a pharmacological mechanism."
      ),
      source_name = "NA_i (Eq. E-3)"
    ),
    REGION_LATINAMERICA = list(
      description = "Latin American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_NORTHAMERICA for the shared Europe reference",
      notes = paste(
        "Carries the largest of the three regional shifts on the dizziness",
        "logit (1.13, odds ratio 3.10 versus Europe). Mutually exclusive",
        "with the other two region indicators."
      ),
      source_name = "LA_i (Eq. E-3)"
    ),
    REGION_ROW = list(
      description = "Rest-of-World study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_NORTHAMERICA for the shared Europe reference",
      notes = paste(
        "'Rest of World' for Gidal 2018 is any study site outside Western",
        "Europe, North America and Latin America. Note that the four-way",
        "region split used by the efficacy models of the same paper",
        "(Western Europe / Latin America / North America, with Rest of",
        "World as the reference) differs from the split used here",
        "(Europe as the reference); the two are not interchangeable and",
        "each model's covariateData records its own reference group."
      ),
      source_name = "ROW_i (Eq. E-3)"
    ),
    CONMED_CBZ = list(
      description = "Carbamazepine use during the baseline period; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (took other antiepileptic drugs at baseline)",
      notes = paste(
        "BASELINE use, i.e. during the 8-week pre-randomisation period,",
        "not concomitant use during treatment. 48% of the safety analysis",
        "population was on carbamazepine during baseline. The positive",
        "coefficient makes carbamazepine users higher-risk for dizziness",
        "(p < 0.05); Gidal 2018 Discussion reads this as a",
        "pharmacodynamic interaction, both drugs being voltage-gated",
        "sodium channel modulators with dizziness as a class effect. Note",
        "the SIGN REVERSES for somnolence in the companion model, where",
        "baseline carbamazepine use is protective."
      ),
      source_name = "BCARB_i (Eq. E-3)"
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no lamotrigine)",
      notes = paste(
        "Increases the risk of dizziness (p < 0.01) and of headache",
        "(p < 0.05) but not of somnolence. Lamotrigine is a third",
        "voltage-gated sodium channel modulator, which Gidal 2018 offers",
        "as the explanation. Lamotrigine has no effect on eslicarbazepine",
        "PK (it was screened and not retained on CL/F or V/F, consistent",
        "with a dedicated phase 1 interaction study), so this is a",
        "pharmacodynamic term only. 14.0% of the pooled population took",
        "lamotrigine."
      ),
      source_name = "LAM_i (Eq. E-3)"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Same orientation as Gidal 2018's SEXF_i indicator, so no value",
        "transformation is needed. Women are predicted to be at higher",
        "risk of dizziness and of somnolence (p < 0.05); Gidal 2018",
        "Results suggests this may be related to their lower body weight,",
        "although weight is already in the model as a separate term."
      ),
      source_name = "SEXF_i (Eq. E-3)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1152L,
    n_studies = 3L,
    n_observations = "1,152 binary first-occurrence dizziness records (one per patient; landmark analysis, no repeated measures)",
    age_median = "37 years",
    weight_median = "70 kg (the centring value used by Eq. E-3)",
    race_ethnicity = c(Caucasian = 80.0),
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo, after a 2-week titration from a 400 or 800 mg starting dose",
    regions = "Europe (reference), North America, Latin America and Rest of World",
    co_medication = "48% took carbamazepine during the baseline period",
    notes = paste(
      "Safety analysis set: 306 patients from study 2093-301, 307 from",
      "2093-302 and 539 from 2093-304. Logistic models were fitted for the",
      "three TEAEs reported in more than 10% of patients: dizziness,",
      "headache and somnolence. Both linear and power exposure models were",
      "screened; AUC0-24 was the significant metric for dizziness and",
      "headache, Cmax for somnolence. Levetiracetam and valproate were",
      "screened and not retained for any of the three endpoints."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-3 and Equation E-3:
    #
    #   logit(p) = -3.67 + 1.97*I(400 mg) + 3.62*I(800 mg)
    #              - 0.00000835*AUC - 0.0247*(WTKG - 70)
    #              + 1.04*NA + 1.13*LA + 0.668*ROW
    #              + 0.626*BCARB + 0.589*LAM + 0.486*SEXF
    #
    # Table S-3 prints a population mean and a %SEM per row and no
    # confidence intervals or odds ratios, so each value is reported once
    # in the table and once in Eq. E-3; both agree exactly.
    # Only WT is centred (at 70 kg); AUC is neither centred nor scaled,
    # so the intercept is the logit for a 70 kg European man on placebo
    # with no baseline carbamazepine and no lamotrigine.
    # ==================================================================

    # ----- Intercept -----
    logit_ref <- -3.67 ; label("Logit of the probability of dizziness for a 70 kg European man on placebo, with no baseline carbamazepine and no concomitant lamotrigine (unitless logit)")  # Table S-3, 'Logit intercept (placebo effect)' -3.67, 8.5% SEM; Eq. E-3. Corresponds to a 2.5% placebo dizziness probability

    # ----- Starting-dose shifts -----
    e_dose_esl_mgd400_logit <- 1.97 ; label("Log-odds shift for a 400 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-3, 'Additive shift in logit for 400 mg starting dose' 1.97, 16.5% SEM; Eq. E-3
    e_dose_esl_mgd800_logit <- 3.62 ; label("Log-odds shift for an 800 mg/day eslicarbazepine acetate starting dose versus placebo (unitless logit)")  # Table S-3, 'Additive shift in logit for 800 mg starting dose' 3.62, 11.9% SEM; Eq. E-3

    # ----- Exposure -----
    e_auc_esl_logit <- -0.00000835 ; label("Log-odds of dizziness per 1 ng*h/mL increase in eslicarbazepine AUC0-24 (unitless logit)")  # Table S-3, 'Slope for eslicarbazepine AUC0-24' -0.00000835, 16.8% SEM; Eq. E-3. Negative by design, not a transcription error -- see the description and vignette Errata item 3

    # ----- Baseline covariates -----
    e_wt_logit <- -0.0247 ; label("Log-odds of dizziness per kg of body weight above 70 kg (unitless logit)")  # Table S-3, 'Slope for weight' -0.0247, 27.2% SEM; Eq. E-3
    e_region_northamerica_logit <- 1.04 ; label("Log-odds shift for a North American versus a European study site (unitless logit)")  # Table S-3, 'Additive shift for North America' 1.04, 29.1% SEM; Eq. E-3
    e_region_latinamerica_logit <- 1.13 ; label("Log-odds shift for a Latin American versus a European study site (unitless logit)")  # Table S-3, 'Additive shift for Latin America' 1.13, 20.2% SEM; Eq. E-3
    e_region_row_logit <- 0.668 ; label("Log-odds shift for a Rest-of-World versus a European study site (unitless logit)")  # Table S-3, 'Additive shift for Rest of World' 0.668, 43.6% SEM; Eq. E-3
    e_conmed_cbz_logit <- 0.626 ; label("Log-odds shift for baseline carbamazepine use (unitless logit)")  # Table S-3, 'Additive shift for baseline carbamazepine use' 0.626, 33.7% SEM; Eq. E-3
    e_conmed_lamotrigine_logit <- 0.589 ; label("Log-odds shift for concomitant lamotrigine use (unitless logit)")  # Table S-3, 'Additive shift for concomitant lamotrigine use' 0.589, 37.9% SEM; Eq. E-3
    e_sexf_logit <- 0.486 ; label("Log-odds shift for female versus male sex (unitless logit)")  # Table S-3, 'Additive shift for females' 0.486, 42.2% SEM; Eq. E-3

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression fitted in NONMEM with
    # a Bernoulli likelihood: no sigma and no random effects are
    # estimated. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity. See the vignette's
    # Assumptions and deviations.
    addSd_prob_dizziness <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Starting-dose indicators derived from the dose-level column ---
    dose400 <- (DOSE_ESL_MGD == 400)
    dose800 <- (DOSE_ESL_MGD == 800)

    # ----- Linear predictor (Gidal 2018 Eq. E-3) -----
    logit_dizziness <- logit_ref +
      e_dose_esl_mgd400_logit * dose400 +
      e_dose_esl_mgd800_logit * dose800 +
      e_auc_esl_logit * AUC_ESL +
      e_wt_logit * (WT - 70) +
      e_region_northamerica_logit * REGION_NORTHAMERICA +
      e_region_latinamerica_logit * REGION_LATINAMERICA +
      e_region_row_logit * REGION_ROW +
      e_conmed_cbz_logit * CONMED_CBZ +
      e_conmed_lamotrigine_logit * CONMED_LAMOTRIGINE +
      e_sexf_logit * SEXF

    prob_dizziness <- expit(logit_dizziness)

    # ----- Observation -----
    prob_dizziness ~ add(addSd_prob_dizziness)
  })
}
