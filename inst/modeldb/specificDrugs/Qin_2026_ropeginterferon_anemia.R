Qin_2026_ropeginterferon_anemia <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for the probability of ",
    "an anemia adverse event during subcutaneous ropeginterferon ",
    "alfa-2b (ropeg) treatment in 179 patients with essential ",
    "thrombocythemia (SURPASS-ET and EXCEED-ET). logit(P) = 5.20 + ",
    "0.024 * CAV - 0.522 * HGB, where CAV is the individual average ",
    "total serum ropeg concentration over the whole treatment period ",
    "(ng/mL) and HGB is baseline hemoglobin (g/dL). Higher exposure and ",
    "lower baseline hemoglobin raise the risk. No PK layer, no ODE and ",
    "no between-subject random effect. Companion models in the ",
    "Qin_2026_ropeginterferon_* family."
  )
  reference <- paste(
    "Qin A, Masarova L, Chen H, Xiao Z, Jin J, Gill H, Reeves BN,",
    "El Chaer F, Kirito K, Komatsu N, Abu-Zeinah G, Shimoda K,",
    "Zagrijtschuk O, Sato T, Wu D, Wu L, Su X, Gao Y, Cui J, Zhang L,",
    "Mesa RA.",
    "Exposure-Response Analyses of Ropeginterferon Alfa-2b in",
    "Essential Thrombocythemia.",
    "J Clin Pharmacol. 2026;66(10):e70258. doi:10.1002/jcph.70258.",
    sep = " "
  )
  vignette <- "Qin_2026_ropeginterferon"
  units <- list(
    time = "n/a (static whole-treatment-period landmark regression; no time dimension)",
    dosing = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "prob_anemia (probability, 0-1, of an anemia adverse event during treatment)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average total serum ropeginterferon alfa-2b concentration throughout treatment, from the first dose to 2 weeks after the last dose.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL serum ropeg (the PK model's Cc), simulated from the actual",
        "dosing records and individual post hoc PK parameters of",
        "modellib('Qin_2026_ropeginterferon') (Qin 2026 Methods,",
        "Exposure-Safety Analysis). Not centred or scaled. Observed range",
        "about 10-92 ng/mL (Qin 2026 Figure 3e)."
      ),
      source_name = "Cavg, throughout the treatment"
    ),
    HGB = list(
      description = "Baseline hemoglobin concentration.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Uncentred linear covariate on the logit, beta2 = -0.522 per",
        "g/dL (Qin 2026 Table 4). Table 1 labels the row 'Hemoglobin",
        "(g /L)' but prints values of 9.30-17.5 with a median of 14.2,",
        "which are g/dL; Figure 3f labels its 12.75, 13.8 and 14.6",
        "quartile lines in g/dL. With g/dL the model reproduces the",
        "observed 29.1% anemia incidence (Table S2) at a typical exposure,",
        "and with g/L it would predict a probability of zero."
      ),
      source_name = "baseline hemoglobin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 179L,
    n_studies = 2L,
    n_observations = "179 binary anemia records, one per patient; 52 events (29.1%), 27/88 in EXCEED-ET and 25/91 in SURPASS-ET (Qin 2026 Table S2)",
    age_range = "median 57.0 (EXCEED-ET) and 61.0 (SURPASS-ET) years, range 21.0-84.0 (Qin 2026 Table 1)",
    weight_range = "median 76.9 (EXCEED-ET) and 60.0 (SURPASS-ET) kg, range 41.5-153 (Qin 2026 Table 1)",
    sex_female_pct = 57.5,
    race_ethnicity = "SURPASS-ET 94.5% Asian; EXCEED-ET 77.3% Caucasian (Qin 2026 Table 1)",
    disease_state = "Essential thrombocythemia; baseline hemoglobin median 13.9 (EXCEED-ET) and 13.6 (SURPASS-ET) g/dL, range 9.30-17.2",
    dose_range = "HIDAT regimen: 250 ug at week 0, 350 ug at week 2, 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    endpoint_definition = "Treatment-emergent adverse event of anemia of any grade (NCI-CTCAE v5.0), one of the TEAEs with incidence above 20% screened in the exposure-safety analysis; only five were grade 3 and none higher (Qin 2026 Methods and Discussion). First occurrence median 18.1 weeks (Table S2)"
  )

  ini({
    # Qin 2026 Table 4, rows 'Logistic regression of occurrence of
    # anemia'; model form Equation (5):
    #   logit(P) = beta0 + beta1*Exposure_i + beta2*HGB_baseline,i
    # Both regressors uncentred. Check against Figure 3f: at CAV = 0 and
    # HGB 12.75 g/dL, expit(5.20 - 0.522*12.75) = 0.19 (figure about
    # 0.18); at CAV = 150, expit(-1.456 + 3.6) = 0.90 (figure about 0.9).
    logit_ref <- 5.20; label("Logit of the probability of anemia at CAV = 0 ng/mL and baseline hemoglobin 0 g/dL (unitless logit; an extrapolated anchor)") # Qin 2026 Table 4: Intercept beta0 5.20 (SE 1.77), P = .00336
    e_cav_logit <- 0.024; label("Log-odds of anemia per 1 ng/mL of CAV (unitless logit)") # Qin 2026 Table 4: Exposure effect beta1 0.024 (SE 0.0106), P = .02155
    e_hgb_logit <- -0.522; label("Log-odds of anemia per 1 g/dL of baseline hemoglobin (unitless logit)") # Qin 2026 Table 4: Baseline hemoglobin effect beta2 -0.522 (SE 0.1306), P < .001

    # A Bernoulli likelihood has no residual; this fixed placeholder only
    # gives rxode2 an error model for the typical-value probability.
    addSd_prob_anemia <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_anemia <- logit_ref + e_cav_logit * CAV + e_hgb_logit * HGB
    prob_anemia <- expit(logit_anemia)
    prob_anemia ~ add(addSd_prob_anemia)
  })
}
