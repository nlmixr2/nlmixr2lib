Li_2025_patritumabDeruxtecan_orr <- function() {
  description <- paste(
    "Bayesian logistic-regression exposure-efficacy model for objective",
    "response (complete or partial response by blinded independent central",
    "review) in adults with previously treated EGFR-mutated non-small-cell",
    "lung cancer given patritumab deruxtecan (HER3-DXd) in the U102 and U201",
    "studies (Li 2025, n = 446; 3.2-6.4 mg/kg IV Q3W or up-titration).",
    "logit(P) is log-linear in the average over Cycles 1-3 of the",
    "anti-HER3-ac-DXd (conjugated payload) trough concentration, with",
    "additive shifts for the up-titration regimen, no prior chemotherapy",
    "and baseline bone metastasis. Li 2025 prints no coefficient table; the",
    "exposure slope, the reference exposure and the regimen and prior",
    "chemotherapy shifts are back-solved from the typical-value forest plot",
    "(Figure 2), and the bone-metastasis shift is the odds ratio printed in",
    "the Results. Covariate-by-exposure interaction terms of the source",
    "model are not recoverable and are omitted. There is no PK layer:",
    "exposure is a data column from the companion population PK model",
    "(see Lee_2023_patritumab). One of nine Li_2025_patritumabDeruxtecan",
    "exposure-response models."
  )
  reference <- paste(
    "Li L, Lee M, Joshi R, Xu Y, Garcia R, Baron K, Anderson E, Waterhouse T,",
    "Sternberg D, Shrestha P, Vaddady P, Abutarif M, Garimella T.",
    "Exposure-Response Relationships in Patients with Non-Small-Cell Lung",
    "Cancer and Other Solid Tumors Treated with Patritumab Deruxtecan",
    "(HER3-DXd). Clin Pharmacol Ther. 2025;118(2):408-417.",
    "doi:10.1002/cpt.3674."
  )
  vignette <- "Li_2025_patritumabDeruxtecan_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_orr_central (probability of objective response, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Anti-HER3-ac-DXd (antibody-conjugated DXd) plasma trough concentration at the end of the dosing interval, AVERAGED over Cycles 1-3, per subject.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Total conjugated-payload analyte, not unconjugated DXd. Averaged over the troughs of Cycles 1-3 so the up-titration regimen (3.2, 4.8, then 6.4 mg/kg) is distinguished from 3.2 mg/kg Q3W (Li 2025 Methods, 'PK exposure metrics'). Empirical-Bayes predictions from the population PK model 'previously developed' (uncited in Li 2025; the same group's model is Lee_2023_patritumab). Centred at 90.5 ng/mL inside model(), the reference exposure implied by Figure 2; Table S2 reports a median of 90.9 ng/mL at 5.6 mg/kg Q3W and 58.4 ng/mL for up-titration. Figure 2 5th/25th/75th/95th percentiles: 19 / 50 / 132 / 206 ng/mL.",
      source_name = "Anti-HER3-ac-DXd Ctrough (Average C1-3 ng/mL)"
    ),
    REGIMEN_UPTITRATION = list(
      description = "Up-titration dosing regimen indicator: 1 = 3.2 mg/kg in Cycle 1, 4.8 mg/kg in Cycle 2, then 6.4 mg/kg Q3W; 0 = a fixed-dose Q3W regimen.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fixed-dose regimen).",
      notes = "Li 2025 added this regimen effect to the ORR model to account for a different exposure-response relationship in patients on up-titration than on fixed doses at the same average Cycle 1-3 trough (Methods). Figure 2 'Titration flag: Up-titration' 22% vs reference 30%.",
      source_name = "Titration flag"
    ),
    PRIOR_CHEMO = list(
      description = "Prior chemotherapy indicator: 1 = received chemotherapy before study entry, 0 = no prior chemotherapy.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 in THIS model: the Figure 2 reference patient had prior chemotherapy, and the forest-plot row is 'No'. The shift is applied as e_noprior_chemo_logit * (1 - PRIOR_CHEMO).",
      notes = "Influential (68% CrI excludes the null) but not significant; only 16 of 446 patients (3.6%) had no prior chemotherapy, so the effect is very uncertain (Li 2025 Results and Discussion). Every patient had received a prior EGFR TKI, so this is chemotherapy specifically, not any systemic therapy.",
      source_name = "Prior chemotherapy"
    ),
    BONE_MET = list(
      description = "Baseline bone metastasis by blinded independent central review: 1 = present, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no bone metastasis).",
      notes = "The only significant ORR covariate: odds ratio 0.48 ('52% (95% CrI, 1.0%, 73%) decrease in the odds of being a responder', Li 2025 Results). 130 of 446 efficacy patients (29.1%) had bone metastasis at baseline.",
      source_name = "Baseline bone metastasis"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 446L,
    n_studies = 2L,
    disease_state = "Locally advanced or metastatic EGFR-mutated non-small-cell lung cancer previously treated with an EGFR tyrosine kinase inhibitor (efficacy-evaluable set).",
    dose_range = "Patritumab deruxtecan 3.2, 4.8, 5.6 or 6.4 mg/kg IV Q3W, or up-titration (3.2 mg/kg Cycle 1, 4.8 mg/kg Cycle 2, then 6.4 mg/kg Q3W).",
    regions = "Global (U102, NCT03260491; U201 HERTHENA-Lung01, NCT04619004).",
    biomarkers = "Baseline bone metastasis by BICR in 29.1%; no prior chemotherapy in 3.6%.",
    notes = "U102 n = 172 and U201 n = 274. Observed ORR 30.5% (136 of 446 responders). The 715-patient exposure-safety population summary (Li 2025 Table S1) is in the safety model files; no separate demographic table is printed for the efficacy set."
  )

  ini({
    # Li 2025 prints no coefficient table. Values below are reconstructed
    # from the typical-value forest plot (Figure 2: printed probabilities at
    # the 5th/25th/75th/95th exposure percentiles and for each covariate
    # level) under the logit-log-linear exposure form of Supplementary
    # methods. The log-linear form reproduces the four exposure rows to
    # within 0.010 logit; the linear form leaves residuals of 0.12 logit.
    logit_ref <- log(0.30 / (1 - 0.30)); label("Logit of objective response probability for the reference patient (unitless logit)") # Figure 2 'Reference patient' 30 [18, 45]
    e_ctrough_logit <- 0.3424; label("Change in the response logit per unit natural log of the average Cycle 1-3 anti-HER3-ac-DXd trough (unitless logit)") # back-solved from Figure 2 exposure rows 19/50/132/206 ng/mL -> 20/26/33/36 (reference 90.5 ng/mL implied)
    e_uptitration_logit <- log(0.22 / 0.78) - log(0.30 / 0.70); label("Shift in the response logit for the up-titration regimen (unitless logit)") # Figure 2 'Titration flag: Up-titration' 22 [11, 38] vs reference 30
    e_noprior_chemo_logit <- log(0.54 / 0.46) - log(0.30 / 0.70); label("Shift in the response logit for no prior chemotherapy (unitless logit)") # Figure 2 'Prior chemotherapy: No' 54 [23, 83] vs reference 30
    e_bone_met_logit <- log(0.48); label("Shift in the response logit for baseline bone metastasis (unitless logit)") # Results: OR 0.48, '52% (95% CrI 1.0%, 73%) decrease in the odds'; Figure 2 row 18 [9, 33]

    # The source likelihood is Bernoulli with no residual error; the
    # placeholder below only lets the nlmixr2 observation machinery accept
    # the model. Not from the source.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Reference exposure: the Figure 2 reference-patient trough implied by
    # the exposure rows (back-solved, 90.5 ng/mL).
    logit_orr_central <- logit_ref +
      e_ctrough_logit * log(CTROUGH / 90.5) +
      e_uptitration_logit * REGIMEN_UPTITRATION +
      e_noprior_chemo_logit * (1 - PRIOR_CHEMO) +
      e_bone_met_logit * BONE_MET
    prob_orr_central <- expit(logit_orr_central)
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
