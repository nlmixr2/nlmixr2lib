Menon_2025_tofacitinib_acr <- function() {
  description <- paste(
    "Landmark ordered-categorical exposure-response model for the joint",
    "probability of ACR20, ACR50 and ACR70 response at Month 3 in adults",
    "with active psoriatic arthritis given placebo or tofacitinib 5 or 10 mg",
    "twice daily in OPAL Broaden and OPAL Beyond (Menon 2025). The three",
    "cumulative logits share one intercept at zero exposure, are separated",
    "by two exponentiated threshold decrements, and share one Emax function",
    "of the individual steady-state average tofacitinib concentration",
    "(CAV, supplied from a population PK model such as",
    "Xie_2019_tofacitinib). No covariates and no random effects."
  )
  reference <- paste(
    "Menon S, Shoji S, Tsuchiwata S, Fallon L, Kanik K. Exposure-Response",
    "Analysis of Tofacitinib in Active Psoriatic Arthritis: Results from Two",
    "Phase 3 Studies. J Clin Pharmacol. 2025;65(3):369-377.",
    "doi:10.1002/jcph.6147. PMID: 39453735. PMCID: PMC11867917."
  )
  vignette <- "Menon_2025_tofacitinib"
  units <- list(
    time = "day (landmark model; the probabilities do not depend on time and apply at Month 3)",
    dosing = "mg",
    concentration = "prob_acr20 / prob_acr50 / prob_acr70 (probability of at least that ACR response at Month 3, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state tofacitinib plasma concentration over the twice-daily dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Patient-specific Cavg predicted from the tofacitinib psoriatic",
        "arthritis population PK model of Xie 2019 (Menon 2025 Methods",
        "'Pharmacokinetics'; modellib('Xie_2019_tofacitinib')), i.e. the",
        "dose divided by the individual CL/F and the 12 h dosing interval.",
        "Set to 0 for placebo, which makes the drug term vanish so the",
        "intercept is the placebo logit. Median (range) 17.6 (6.67-47.2) and",
        "36.1 (15.1-88.1) ng/mL at 5 and 10 mg twice daily (Results,",
        "hemoglobin section). Time-fixed per patient."
      ),
      source_name = "Cavg"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested as a power function (WT/median WT)^theta8 on EC50 (Table 1",
        "runs 5 and 8) and on Emax (run 6); dOFV -0.022 and -0.457, 90% CI",
        "including 0, not retained (Results). No estimate is printed."
      ),
      source_name = "WT"
    ),
    PRIOR_TNF = list(
      description = "Prior tumor necrosis factor inhibitor therapy (OPAL Beyond population)",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = paste(
        "Tested as a shift on Int and Emax (Table 1 runs 7 and 8); dOFV",
        "-1.174 for two parameters, not retained (Results). No estimate is",
        "printed."
      ),
      source_name = "study (OPAL Beyond)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 655L,
    n_studies = 2L,
    age_range = "adults >= 18 years (median 50 years in the hemoglobin analysis set)",
    weight_range = "median 83 kg in the hemoglobin analysis set",
    sex_female_pct = 54.9,
    disease_state = paste(
      "Active psoriatic arthritis (CASPAR criteria, >= 3 swollen and >= 3",
      "tender joints, active plaque psoriasis) with inadequate response to",
      "conventional synthetic DMARDs (OPAL Broaden, TNF-inhibitor naive) or",
      "to TNF inhibitors (OPAL Beyond), on a stable background csDMARD"
    ),
    dose_range = "Placebo or tofacitinib 5 or 10 mg orally twice daily; ACR response at Month 3 (adalimumab arm not analysed)",
    regions = "Multinational (OPAL Broaden NCT01877668, OPAL Beyond NCT01882439)",
    n_observations = "655 Month-3 ordered-categorical ACR observations, one per patient (OPAL Broaden 298, OPAL Beyond 357)",
    notes = paste(
      "Menon 2025 Results; 670 patients were randomised to placebo or",
      "tofacitinib but 15 had no Cavg (Table 1 footnote b). Sex and",
      "demographic medians are those printed for the 685-patient hemoglobin",
      "analysis set; the paper does not print demographics for the ACR set.",
      "Fitted by maximum likelihood in SAS 9.4 PROC NLMIXED."
    )
  )

  ini({
    # logit Pr(sACR >= m) = Int - I(m >= 2) exp(theta2) - I(m >= 3) exp(theta3)
    #                       + Emax * CAV / (CAV + exp(theta5))
    # (Menon 2025 Methods equation). sACR = ACR20 + ACR50 + ACR70 responder
    # flags, so m = 1, 2, 3 are ACR20, ACR50 and ACR70.
    logit_bl_acr20 <- -0.852; label("Logit of ACR20 response at zero exposure (placebo), Month 3 (logit)")  # Table 2, theta1 'Int' = -0.852 (ASE 0.145)
    lgap_acr50 <- -0.0636; label("Log of the logit decrement from the ACR20 to the ACR50 threshold (log logit)")  # Table 2, theta2 'Log ACR50 threshold' = -0.0636 (ASE 0.0807); exp = 0.938
    lgap_acr70 <- -0.201; label("Log of the logit decrement from the ACR50 to the ACR70 threshold (log logit)")  # Table 2, theta3 'Log ACR70 threshold' = -0.201 (ASE 0.106); exp = 0.818
    emax_acr <- 1.76; label("Maximum tofacitinib effect on the cumulative logits (logit)")  # Table 2, theta4 'Emax' = 1.76 (ASE 0.52)
    lec50_acr <- 2.84; label("Log Cavg giving half the maximum effect (log ng/mL)")  # Table 2, theta5 'Log EC50' = 2.84 (ASE 0.79); exp = 17.1 ng/mL (90% CI 4.66 to 62.9)

    # The source likelihood is the multinomial ordered-categorical one and
    # estimates no residual error; these keep the nlmixr2 observation
    # machinery working and are not published quantities.
    addSd_prob_acr20 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr20 (not from source)")  # not from source; see vignette Assumptions and deviations
    addSd_prob_acr50 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr50 (not from source)")  # not from source; see vignette Assumptions and deviations
    addSd_prob_acr70 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr70 (not from source)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50_acr)
    drug <- emax_acr * CAV / (CAV + ec50)

    logit_acr20 <- logit_bl_acr20 + drug
    logit_acr50 <- logit_acr20 - exp(lgap_acr50)
    logit_acr70 <- logit_acr50 - exp(lgap_acr70)

    prob_acr20 <- expit(logit_acr20)
    prob_acr50 <- expit(logit_acr50)
    prob_acr70 <- expit(logit_acr70)

    prob_acr20 ~ add(addSd_prob_acr20)
    prob_acr50 ~ add(addSd_prob_acr50)
    prob_acr70 ~ add(addSd_prob_acr70)
  })
}
