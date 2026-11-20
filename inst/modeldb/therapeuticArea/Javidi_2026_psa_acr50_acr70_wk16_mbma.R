Javidi_2026_psa_acr50_acr70_wk16_mbma <- function() {
  description <- paste0(
    "MBMA. Week-16 cross-endpoint bridging model ('Model 2', 16-week fit) of ",
    "Javidi 2026 for psoriatic arthritis: the logit of a study arm's ACR50 ",
    "(prob_acr50) and ACR70 (prob_acr70) responder rate at week 16 is a linear ",
    "function of the same arm's week-16 ACR20 responder proportion, ",
    "logit(ACR50 or ACR70) = slope * ACR20 + intercept, with an ",
    "endpoint-specific intercept and a slope that depends on drug class ",
    "(placebo, PDE4 inhibitor, anti-TNF, anti-IL-17A, anti-IL-17A/F, ",
    "anti-IL-23, anti-IL-12/23, TYK2 inhibitor), endpoint, study phase and ",
    "the arm's percentage of prior-TNF-experienced patients. Fitted with nlme ",
    "gnls (fixed effects only, no random effects) to 27 psoriatic arthritis ",
    "trials (78 arms, 11,326 participants) reporting ACR20 together with ACR50 ",
    "and/or ACR70. The ACR20 input is a covariate (RESP_ACR20_PCT): an observed ",
    "arm value, or the prob_acr20 output of ",
    "modellib('Javidi_2026_psoriaticDisease_bridging_mbma') times 100, which is ",
    "how the source chains the two models to project ACR50 / ACR70 in an ",
    "unstudied indication. Valid at week 16 only; the separately fitted ",
    "week-24 model is modellib('Javidi_2026_psa_acr50_acr70_wk24_mbma'). ",
    "The residual is additive on the logit scale with variance ",
    "1 / (ACR20 (1 - ACR20) N), sigma fixed to 1. Simulation scope is ",
    "study-arm responder rates, NOT individual patients."
  )

  reference <- paste(
    "Javidi M, Alexander GE, Kosoglou T, Vermeulen A, Beutler A, Shagroni T,",
    "Valiathan C.",
    "MBMA Bridging Models as a Tool for Exploration of Clinical Endpoints in",
    "Unstudied Indications.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5):726-740.",
    "doi:10.1002/psp4.70207. PMC13142780.",
    "Structural model: main-text Equation 5 and Supporting Information",
    "Equation S.6, together with the authors' gnls code for Model 2",
    "('Model for data at 16 week') in the Supporting Information (Data S1).",
    "Parameter estimates: Supporting Information Table S3 (16-week rows).",
    sep = " "
  )

  vignette <- "Javidi_2026_psoriaticDisease_mbma"

  units <- list(
    time = "week (the model has no time dependence; it applies to week-16 responder rates only)",
    dosing = "n/a (no dose events; drug class enters through TRT_<class> arm indicators)",
    concentration = "probability/arm (prob_acr50 and prob_acr70 are STUDY-ARM responder probabilities on a 0-1 scale, NOT drug concentrations. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    RESP_ACR20_PCT = list(
      description = "Week-16 ACR20 responder rate of the same study arm, in percent (0-100).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "The model's predictor: divided by 100 inside model() because the Data S1 gnls code regresses on the proportion. Javidi 2026 fits the model on observed arm ACR20 values; in the paper's applications the input is the Model 1 projection of ACR20 (prob_acr20 of Javidi_2026_psoriaticDisease_bridging_mbma times 100). It also sets the residual weight 1 / (p (1 - p) N) with p = ACR20 proportion, per the Data S1 definition of mdl.weights. Must lie strictly between 0 and 100.",
      source_name = "ACR20_ij (Javidi 2026 Equation 5); p.x in the Data S1 gnls code (whose comment block labels p.x and p.y the other way round -- see the vignette Assumptions and deviations)"
    ),
    TRT_PDE4_INHIBITOR = list(
      description = "Study-arm drug-class indicator: 1 = the arm received a phosphodiesterase-4 inhibitor (apremilast in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "The eight class levels (placebo plus seven drug classes) are mutually exclusive; placebo is the arm with every class indicator 0. Arm-level.",
      source_name = "PDE4inhibitor (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_ANTITNF = list(
      description = "Study-arm drug-class indicator: 1 = the arm received an anti-TNF agent (adalimumab or etanercept in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; see TRT_PDE4_INHIBITOR.",
      source_name = "antiTNF (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_ANTIIL17 = list(
      description = "Study-arm drug-class indicator: 1 = the arm received an anti-IL-17A agent (secukinumab or ixekizumab in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; the dual IL-17A/F inhibitor bimekizumab has its own class, TRT_ANTIIL17AF. See TRT_PDE4_INHIBITOR.",
      source_name = "antiIL17 (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_ANTIIL17AF = list(
      description = "Study-arm drug-class indicator: 1 = the arm received a dual anti-IL-17A/IL-17F agent (bimekizumab in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; see TRT_PDE4_INHIBITOR.",
      source_name = "antiIL17AF (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_ANTIIL23 = list(
      description = "Study-arm drug-class indicator: 1 = the arm received an anti-IL-23 p19 agent (guselkumab or risankizumab in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; see TRT_PDE4_INHIBITOR.",
      source_name = "antiIL23 (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_ANTIIL12_23 = list(
      description = "Study-arm drug-class indicator: 1 = the arm received an anti-IL-12/23 p40 agent (ustekinumab in the source database), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; see TRT_PDE4_INHIBITOR.",
      source_name = "antiIL1223 (Javidi 2026 Data S1 Model 2 code)"
    ),
    TRT_TYK2_INHIBITOR = list(
      description = "Study-arm drug-class indicator: 1 = the arm received a TYK2 inhibitor (brepocitinib in the source database; deucravacitinib in the paper's external validation), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Arm-level drug-class indicator; see TRT_PDE4_INHIBITOR. Week 16 only: the 24-week data set had no TYK2-inhibitor arm, so the week-24 model has no such class.",
      source_name = "TYK2inhibitor (Javidi 2026 Data S1 Model 2 code)"
    ),
    STUDY_PHASE3 = list(
      description = "Study-level indicator: 1 = phase 3 trial, 0 = phase 2 trial.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Shifts the slope by slp.ph3 on every arm of a phase 3 trial (code column ph3). How the phase 2a/2b, 3b and 4 trials of Table S5 were binned is not stated; see the Javidi_2026_psoriaticDisease_bridging_mbma notes.",
      source_name = "ph3 (Javidi 2026 Data S1 Model 2 code); alpha1_i (Equation S.6)"
    ),
    PRIOR_TNF_PCT = list(
      description = "Percentage (0-100) of the study arm's patients who had previously received an anti-TNF agent.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "The source regresses on a 'normalized' value of this percentage without printing the normalisation. Inside model() it is converted as (PRIOR_TNF_PCT - 57) / 19: the centre 57% is the 'average of % prior TNF-experienced' the paper states as the reference of Figures S15 and S20; the scale of 19 percentage points per unit was recovered by the maintainers by fitting the Figure S19 prediction lines (the logit change across 0-100% implies about 19 percentage points per unit; the centre is not identifiable from that figure). The slope coefficient is 0.04 (SE 0.04), so the covariate moves predictions by well under one percentage point across its range; the normalisation is immaterial to the paper's conclusions but is an assumption of this implementation. Javidi 2026 Table 1 simulates a 50% (mixed) population.",
      source_name = "PTnf (Javidi 2026 Data S1 Model 2 code: 'normalized observed values of % population with prior-TNF experienced'); prTNF_j (Equation S.6)"
    ),
    N_ARM = list(
      description = "Number of participants in the study arm.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Sets the residual weight mdl.weights = 1 / (p (1 - p) N) on the logit scale (Data S1 gnls weights = varFixed(~mdl.weights), sigma fixed to 1). Typical-value simulations do not depend on it.",
      source_name = "n (Javidi 2026 Data S1 Model 2 code, definition of mdl.weights)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11326L,
    n_studies = 27L,
    n_arms = 78L,
    age_range = "not reported at arm level",
    weight_range = "not reported at arm level",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state = "adults with active psoriatic arthritis in randomised, parallel-group phase 2 / phase 3 trials",
    dose_range = "approved and dose-ranging regimens of apremilast, adalimumab, etanercept, secukinumab, ixekizumab, bimekizumab, guselkumab, risankizumab, ustekinumab and brepocitinib, plus placebo (Supporting Information Table S5, rows marked MDL2)",
    timepoints = "week 16 arm responder rates (the 16- and 24-week models were fitted separately; this file is the 16-week fit)",
    regions = "international; trials published 2000-2023 (same curated database as the cross-indication model)",
    notes = "MBMA at the study-arm level; fixed effects only (nlme gnls), so there is no between-trial variability. The 27 studies / 78 arms / 11,326 participants are the Model 2 totals across both time points (Javidi 2026 Results 'Data'); the paper does not split them by week. The deucravacitinib psoriatic arthritis trial IM011-084 was excluded from the fit and used for external validation."
  )

  ini({
    # ============================================================
    # Javidi 2026 Equation 5 and Equation S.6, week-16 fit:
    #   logit(Y_s) = slp_s * ACR20 + intc_s
    #   slp = alpha1*slp.ph3 + alpha2*slp.acr + slp.drg.class
    #         + slp.prtnf * prTNF
    # with alpha2 = 1 for ACR50 (Table S3: 'Shift in slope for ACR50
    # relative to ACR70'). All values Javidi 2026 Supporting
    # Information Table S3, 16-week rows.
    # ============================================================
    slp_placebo <- 3.05
    label("ACR20-to-logit slope for placebo arms, week 16 (paper: slp.drg.class placebo,16wk; unitless)") # Table S3 slp.drg.class placebo,16wk = 3.05 (SE 0.67)
    slp_pde4_inhibitor <- 4.17
    label("ACR20-to-logit slope for PDE4-inhibitor arms, week 16 (paper: slp.drg.class PDE4 inhibitor,16wk; unitless)") # Table S3 slp.drg.class PDE4 inhibitor,16wk = 4.17 (SE 0.59)
    slp_antitnf <- 4.86
    label("ACR20-to-logit slope for anti-TNF arms, week 16 (paper: slp.drg.class anti-TNF,16wk; unitless)") # Table S3 slp.drg.class anti-TNF,16wk = 4.86 (SE 0.44)
    slp_antiil17 <- 4.98
    label("ACR20-to-logit slope for anti-IL-17A arms, week 16 (paper: slp.drg.class anti-IL-17,16wk; unitless)") # Table S3 slp.drg.class anti-IL-17,16wk = 4.98 (SE 0.42)
    slp_tyk2_inhibitor <- 4.49
    label("ACR20-to-logit slope for TYK2-inhibitor arms, week 16 (paper: slp.drg.class TYK2 inhibitor,16wk; unitless)") # Table S3 slp.drg.class TYK2 inhibitor,16wk = 4.49 (SE 0.39)
    slp_antiil17af <- 4.45
    label("ACR20-to-logit slope for anti-IL-17A/F arms, week 16 (paper: slp.drg.class anti-IL-17A/F,16wk; unitless)") # Table S3 slp.drg.class anti-IL-17A/F,16wk = 4.45 (SE 0.37)
    slp_antiil23 <- 4.29
    label("ACR20-to-logit slope for anti-IL-23 arms, week 16 (paper: slp.drg.class anti-IL-23,16wk; unitless)") # Table S3 slp.drg.class anti-IL-23,16wk = 4.29 (SE 0.43)
    slp_antiil12_23 <- 4.92
    label("ACR20-to-logit slope for anti-IL-12/23 arms, week 16 (paper: slp.drg.class anti-IL-12/23,16wk; unitless)") # Table S3 slp.drg.class anti-IL-12/23,16wk = 4.92 (SE 0.54)
    slp_acr50 <- -0.74
    label("Shift in slope for the ACR50 endpoint relative to ACR70, week 16 (paper: slp.acr 16wk; unitless)") # Table S3 slp.acr 16wk = -0.74 (SE 0.36)
    slp_ph3 <- -0.29
    label("Shift in slope for phase 3 trials, week 16 (paper: slp.ph3 16wk; unitless)") # Table S3 slp.ph3 16wk = -0.29 (SE 0.17)
    slp_prior_tnf <- 0.04
    label("Change in slope per unit of normalised percentage prior-TNF-experienced, week 16 (paper: slp.prtnf 16wk; unitless)") # Table S3 slp.prtnf 16wk = 0.04 (SE 0.04)
    intc_acr50 <- -2.95
    label("Logit intercept for ACR50, week 16 (paper: intc acr50,16wk; unitless log-odds)") # Table S3 intc acr50,16wk = -2.95 (SE 0.15)
    intc_acr70 <- -4.24
    label("Logit intercept for ACR70, week 16 (paper: intc acr70,16wk; unitless log-odds)") # Table S3 intc acr70,16wk = -4.24 (SE 0.21)

    # Residual on the logit scale: gnls weights = varFixed(~mdl.weights)
    # with gnlsControl(sigma = 1), mdl.weights = 1/(p (1 - p) n), p =
    # observed ACR20. Carried per output as a fixed(1) multiplier.
    addSd_prob_acr50 <- fixed(1)
    label("Logit-scale residual SD multiplier for ACR50; the record SD is sqrt(1 / (p (1 - p) N_ARM)) times this, p = ACR20 proportion (unitless; the source sets sigma = 1)") # Data S1 Model 2 code gnlsControl(sigma = 1)
    addSd_prob_acr70 <- fixed(1)
    label("Logit-scale residual SD multiplier for ACR70; the record SD is sqrt(1 / (p (1 - p) N_ARM)) times this, p = ACR20 proportion (unitless; the source sets sigma = 1)") # Data S1 Model 2 code gnlsControl(sigma = 1)
  })

  model({
    # No ODEs, no time dependence: an arm-level regression at week 16.
    acr20 <- RESP_ACR20_PCT / 100

    # Drug-class slope; placebo is the arm with every class indicator 0.
    trtClass <- TRT_PDE4_INHIBITOR + TRT_ANTITNF + TRT_ANTIIL17 +
      TRT_TYK2_INHIBITOR + TRT_ANTIIL17AF + TRT_ANTIIL23 + TRT_ANTIIL12_23
    slpClass <- slp_placebo * (1 - trtClass) +
      slp_pde4_inhibitor * TRT_PDE4_INHIBITOR +
      slp_antitnf * TRT_ANTITNF +
      slp_antiil17 * TRT_ANTIIL17 +
      slp_tyk2_inhibitor * TRT_TYK2_INHIBITOR +
      slp_antiil17af * TRT_ANTIIL17AF +
      slp_antiil23 * TRT_ANTIIL23 +
      slp_antiil12_23 * TRT_ANTIIL12_23

    # Normalised prior-TNF percentage: centre 57% (paper-stated average),
    # scale 19 percentage points per unit (recovered from Figure S19).
    priorTnfNorm <- (PRIOR_TNF_PCT - 57) / 19

    slpAcr70 <- slpClass + slp_ph3 * STUDY_PHASE3 + slp_prior_tnf * priorTnfNorm
    slpAcr50 <- slpAcr70 + slp_acr50

    prob_acr50 <- expit(slpAcr50 * acr20 + intc_acr50)
    prob_acr70 <- expit(slpAcr70 * acr20 + intc_acr70)

    # Residual on the logit scale (Data S1 mdl.weights).
    sdLogit <- sqrt(1 / (acr20 * (1 - acr20) * N_ARM))
    sdAcr50 <- sdLogit * addSd_prob_acr50
    sdAcr70 <- sdLogit * addSd_prob_acr70
    prob_acr50 ~ logitNorm(sdAcr50, 0, 1)
    prob_acr70 ~ logitNorm(sdAcr70, 0, 1)
  })
}
