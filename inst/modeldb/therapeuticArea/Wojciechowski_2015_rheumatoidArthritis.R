Wojciechowski_2015_rheumatoidArthritis <- function() {
  description <- paste0(
    "Disease-progression model of 28-joint disease activity score ",
    "(DAS28) over time since initiation of triple disease-modifying ",
    "anti-rheumatic drug (DMARD) therapy (methotrexate + sulfasalazine ",
    "+ hydroxychloroquine) in early rheumatoid arthritis. Wojciechowski ",
    "2015 fit DAS28 from initiation until 60 weeks in 263 patients at ",
    "the Royal Adelaide Hospital Early Arthritis Clinic. The structural ",
    "model is an exponential decay of logit-transformed DAS28 toward a ",
    "new treated steady state, additive to the (logit) baseline: ",
    "DAS28_logit(t) = BASE * (1 + e_conmed_steroid_base * CONMED_STEROID) ",
    "* (AGE/57)^e_age_base ",
    "+ EX1 * (1 + e_conmed_steroid_fu_ex1 * CONMED_STEROID_FU) ",
    "* (1 - exp(-EX2 * (1 + e_smoke_ex2 * SMOKE) * t)). ",
    "The natural-scale DAS28 (0-9.2 plausible range) is recovered by ",
    "DAS28 = 9.2 / (1 + exp(-DAS28_logit)) -- a logit back-transform ",
    "with upper bound 9.2 corresponding to 28 swollen and tender joint ",
    "counts, 100 mm visual analogue scale for patient global assessment, ",
    "and 120 mm/h ESR. Random effects on BASE and EX1 are additive on ",
    "the linear (logit-domain) scale; EX2 carries log-normal IIV; all ",
    "three etas are correlated via a full 3x3 OMEGA BLOCK (Table 2 ",
    "off-diagonal entries reported as correlations). Combined ",
    "proportional + additive residual error is applied to the predicted ",
    "logit-DAS28 (paper Equation 3). The model has no PK component -- ",
    "DMARD doses were titrated to disease severity in a treat-to-target ",
    "protocol and were not retained as covariates. CONMED_STEROID is ",
    "time-varying per record (1 at clinic visits where intramuscular ",
    "corticosteroids were administered) and CONMED_STEROID_FU is time-",
    "fixed per subject (1 if the subject received any systemic ",
    "corticosteroid -- i.m. or oral -- at any point during the 60-week ",
    "follow-up window)."
  )

  reference <- paste(
    "Wojciechowski J, Wiese MD, Proudman SM, Foster DJR, Upton RN.",
    "A population model of early rheumatoid arthritis disease activity",
    "during treatment with methotrexate, sulfasalazine and",
    "hydroxychloroquine.",
    "Br J Clin Pharmacol 2015;79(5):777-88.",
    "doi:10.1111/bcp.12553.",
    sep = " "
  )

  vignette <- "Wojciechowski_2015_rheumatoidArthritis"

  # Paper-specific structural parameters: BASE, EX1, and EX2 are
  # signed linear-scale (BASE, EX1) and log-normal-IIV (EX2)
  # disease-progression parameters on the logit-DAS28 domain. They are
  # not log-transformed PK-style positive parameters, so the canonical
  # `eta + l<base>` pairing does not apply for BASE and EX1. Declare
  # `etabase` and `etaex1` as paper-specific etas (the same pattern as
  # `Lee_2011_parkinson_progression`'s `etaslope` / `etasymeff`); the
  # log-normal `etalex2` follows the standard `eta + l<base>`
  # convention against `lex2`.
  paper_specific_etas <- c("etabase", "etaex1")

  # Paper-mechanistic single-output observation `das28_logit` (the
  # logit-transformed DAS28 score on which residual error is applied;
  # paper Equation 3) and its natural-scale back-transform companion
  # `das28` (0-9.2 dimensionless; paper Methods "Base model
  # development"). The canonical PK-output name `Cc` does not apply --
  # there is no drug concentration in this disease-progression model.
  # The DAS28 score is paper-mechanistic and does not generalise to
  # other models, so it is declared paper-specific here following the
  # same pattern as `Wu_2014_FEV1_asthma`'s `fev1` and
  # `Lee_2011_parkinson_progression`'s `deltaUPDRS`.
  paper_specific_compartments <- c("das28_logit", "das28")

  units <- list(
    time          = "week (time since initiation of triple DMARD therapy)",
    dosing        = "(none; no PK component -- DMARD doses were titrated to disease severity in a treat-to-target protocol and were not retained as covariates)",
    concentration = "(28-joint disease activity score DAS28, dimensionless 0-9.2; observation das28_logit on the logit-transformed scale per paper Equation 3, with the back-transformed das28 also emitted for visualisation)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age at initiation of triple DMARD therapy.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject time-fixed (the source paper takes age at the initiation of triple therapy, not time-varying age). Centered at the median age of the cohort (57 years; Table 1) via a power function: BASE_j = TVBASE * (AGE_j / 57)^e_age_base. The continuous-covariate power-around-median form is paper Equation 5; the AGE specialisation is paper Equation 8 and the worked example in the text (BASE = 5.4 at age 40 and 5.8 at age 70 after the back-transform from the logit scale to the 0-9.2 natural-scale DAS28).",
      source_name        = "AGE"
    ),
    CONMED_STEROID = list(
      description        = "Time-varying per-record indicator for intramuscular corticosteroid administration at the current clinic visit. 1 at visits where an i.m. methylprednisolone acetate dose (40-80 mg) was administered, 0 otherwise. In the source paper this is the column CSIM (case-sensitive).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no i.m. corticosteroid administration at the current visit).",
      notes              = "Time-varying per record. In the source paper this is the i.m.-route subset of the systemic-corticosteroid concept; oral or intra-articular corticosteroid administration does NOT set CONMED_STEROID = 1 in the Wojciechowski 2015 dataset (paper Methods 'Patient data' and Results paragraph following Eq 8: 'where CSIM = 1 ... at times where the patient received i.m. corticosteroids and CSIM = 0 ... if not'). The existing CONMED_STEROID canonical's umbrella description -- 'systemic corticosteroid administration indicator' covering 'time-varying per record, capturing acute corticosteroid pulses for relapse / flare treatment' -- already permits route-restricted i.m.-pulse encodings as a subset of systemic. Enters BASE as a multiplicative shift on the logit scale: BASE_j = TVBASE * (1 + e_conmed_steroid_base * CONMED_STEROID) * (AGE / 57)^e_age_base. Worked example from the paper text: a typical 57-year-old patient (no smoking) has BASE_logit = 0.472 (natural-scale DAS28 = 5.7) when CONMED_STEROID = 0; the same patient at a visit with i.m. corticosteroid administration has BASE_logit = 0.472 * (1 + 0.737) = 0.820, corresponding to a natural-scale DAS28 of 6.4. The route restriction is recorded here in source_name (CSIM) and in this notes block so a downstream user knows that an i.m.-only signal is the source-paper convention.",
      source_name        = "CSIM"
    ),
    CONMED_STEROID_FU = list(
      description        = "Per-subject time-fixed indicator equal to 1 if the subject received any systemic (i.m. or oral) corticosteroid at any time during the 60-week analysis follow-up window, 0 otherwise. In the source paper this is the column CSSYS (case-sensitive).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no systemic corticosteroid administration anywhere in the 60-week follow-up window).",
      notes              = "Per-subject time-fixed (the source paper aggregates i.m. or oral corticosteroid administration across the entire 60-week follow-up into a single per-subject 0/1 indicator). New canonical entry registered alongside this model in inst/references/covariate-columns.md (Operator-approved sidecar 2026-06-19, request-001 / response-001). Distinct from the existing CONMED_STEROID canonical (which is per-record acute administration or per-subject baseline / chronic concurrent use) and from PRICORT (strictly pre-study). Enters EX1 as a multiplicative shift on the logit-domain extent-of-response: EX1_j = TVEX1 * (1 + e_conmed_steroid_fu_ex1 * CONMED_STEROID_FU). Worked example from the paper text: TVEX1 = -1.28, so for CONMED_STEROID_FU = 0 EX1 = -1.28 and for CONMED_STEROID_FU = 1 EX1 = -1.28 * (1 + (-0.237)) = -0.977 (paper Eq 8 text: 'EX1 = -0.977 if the patient received systemic (i.m./oral) corticosteroids at any time point throughout the 60 week period and EX1 = -1.28 if they did not'). The source authors interpret this as reflecting that the subset of patients who needed steroid rescue had a less favourable disease-modifying-response profile rather than as a causal pharmacological effect of corticosteroids (paper Discussion: 'within our cohort, single dose i.m. corticosteroids and low dose oral corticosteroids were usually given to individuals with high disease activity, so the improved ability of the model to describe disease activity may simply be reflecting this practice').",
      source_name        = "CSSYS"
    ),
    SMOKE = list(
      description        = "Current-smoker binary indicator at baseline. 1 = current smoker, 0 = never or past smoker.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (never or past smoker; pooled per the paper's covariate analysis 'EX2 as never/past versus current smoker').",
      notes              = "Per-subject time-fixed (baseline smoking status at the initiation of triple therapy; Table 1 reports 47% never, 18% current, 30% past, 5% missing). The source paper's smoking covariate on EX2 collapses never and past smokers into the reference category and contrasts them against current smokers (Results paragraph 'smoking status ... as never/past versus current smoker'); this matches the SMOKE canonical (2-level current vs non-smoker) rather than the 3-level SMOKE_CURRENT + SMOKE_NEVER pair (which would split former from never). Enters EX2 as a multiplicative shift on the rate constant: EX2_j = TVEX2 * (1 + e_smoke_ex2 * SMOKE). Worked example from the paper text: TVEX2 = 0.111 /week (half-life 6.2 weeks) for non-smokers; current smokers have EX2 = 0.111 * (1 + (-0.398)) = 0.0668 /week (half-life 10.4 weeks; paper Description of the final model: 'a typical population half-life of 6.2 weeks (10.4 weeks for current smokers)').",
      source_name        = "Smoking (current = 1)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as a covariate on the EX1 / EX2 parameters during the multivariate covariate analysis but not retained in the final model (paper Methods 'Covariate analyses'). Excluded here to preserve the screening provenance without triggering an unreferenced-covariate warning."
    ),
    HT = list(
      description = "Standing height at baseline.",
      units       = "cm",
      type        = "continuous",
      notes       = "Tested but not retained in the final model (paper Methods 'Covariate analyses')."
    ),
    SEXF = list(
      description = "Female-sex indicator at baseline.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a covariate on BASE (the paper's example of binary-covariate encoding is the gender effect on BASE; paper Methods 'Covariate analyses' and Eq 4 example) but not retained in the final model. The source authors attribute the lack of a retained sex effect to the treat-to-target protocol's dose-titration absorbing the variability (paper Discussion 'Titrating drug doses based on disease severity could explain why gender ... did not significantly affect disease activity, as their effects were accounted for by higher doses')."
    ),
    RHEUMATOID_FACTOR = list(
      description = "Rheumatoid factor positivity at baseline.",
      units       = "(binary; 1 = positive)",
      type        = "binary",
      notes       = "Tested but not retained in the final model. Paper Discussion: 'The presence of the SE, anti-CCP antibodies, RF and female gender are documented as predictors of poorer response ... However, when added on their own, none was considered to be a significant contributor to improving the fit of the model.'"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 263L,
    n_studies      = 1L,
    age_range      = "18-86 years (Table 1; inclusion criteria required age > 18)",
    age_median     = "57 years (Table 1; this is also the centering value used by the AGE covariate effect on BASE)",
    weight_range   = "40.4-143.5 kg (Table 1; 22% missing)",
    weight_median  = "72.5 kg (Table 1)",
    sex_female_pct = 71,
    race_ethnicity = "Predominantly Caucasian (>=91% per paper Methods 'Covariate analyses'; ethnicity not tested as a covariate due to low power)",
    disease_state  = "Early rheumatoid arthritis, diagnosed per the 1987 Revised ACR Criteria, with no prior use of DMARDs. All subjects enrolled at the Royal Adelaide Hospital Early Arthritis Clinic between September 1998 and March 2012. Baseline DAS28 median 5.6 (range 1.8-8.5), indicating high disease activity at the initiation of triple therapy.",
    dose_range     = "Treat-to-target triple DMARD regimen per Proudman et al. [3]: methotrexate 10 mg/week (with folic acid 0.5 mg/day) titrated up to 25 mg/week; sulfasalazine 500 mg/day titrated up to 3000 mg/day (1500 mg twice daily); hydroxychloroquine 200 mg twice daily. Intra-articular (40-80 mg methylprednisolone acetate) and intramuscular (80-120 mg methylprednisolone acetate) corticosteroid injections were administered at the physician's discretion; oral corticosteroids and NSAIDs were actively discouraged but used if deemed necessary. Doses are not retained as covariates in the final model.",
    regions        = "Australia (single-centre, Royal Adelaide Hospital).",
    smoking_status = c(Never = 47, Current = 18, Past = 30, Missing = 5),
    corticosteroid_use_during_followup_pct = 48,
    notes          = "Total of 2080 DAS28 observations across 263 patients (median 8 observations per patient, range 1-15). Visits typically occurred every 3-6 weeks until disease was stable, then every 3 months. Follow-up median 52 weeks (range 0-60); the model was fit over the 0-60-week window. Missing continuous covariate values were imputed with the median and missing categorical values with the mode prior to population modelling (paper Methods 'Patient data')."
  )

  ini({
    # ============================================================
    # Structural model (paper Equations 7 and 8):
    #
    #   DAS28_logit(t) = BASE * (1 + e_conmed_steroid_base * CONMED_STEROID)
    #                  * (AGE / 57)^e_age_base
    #                  + EX1 * (1 + e_conmed_steroid_fu_ex1 * CONMED_STEROID_FU)
    #                  * (1 - exp(-EX2 * (1 + e_smoke_ex2 * SMOKE) * t))
    #
    # Natural-scale DAS28 is the back-transform from the logit domain:
    #
    #   das28 = 9.2 / (1 + exp(-DAS28_logit))
    #
    # All parameter point estimates and 95% CIs are from Wojciechowski
    # 2015 Table 2 (final model, NONMEM FOCE-I). The OMEGA off-diagonal
    # entries in Table 2 are correlations (not covariances) -- the
    # convention has been verified by Cauchy-Schwarz: the reported
    # -0.376 between sd_BASE = 0.492 and sd_EX1 = 0.740 exceeds 0.492 *
    # 0.740 = 0.364 in absolute value, which is mathematically
    # impossible for a covariance. Treating the off-diagonals as
    # correlations resolves this. See the validation vignette's
    # Assumptions and deviations section for the worked transcription.
    # ============================================================

    base <- 0.472   ; label("Typical baseline disease activity in the logit-DAS28 domain (dimensionless; back-transformed natural-scale DAS28 = 9.2 / (1 + exp(-base)) ~= 5.7 for a 57-year-old non-smoker with no concurrent i.m. corticosteroids)")
    # Table 2 theta_1 = 0.472 (95% CI 0.392, 0.552), the typical BASE in the logit-DAS28 domain.

    ex1  <- -1.28   ; label("Typical asymptotic change in the logit-DAS28 domain from initiation to a new treated steady state (dimensionless; corresponds to a typical natural-scale change in DAS28 of about -2.8 from baseline 5.7 to a treated steady state of about 2.9)")
    # Table 2 theta_2 = -1.28 (95% CI -1.39, -1.17), the typical EX1 in the logit-DAS28 domain.

    lex2 <- log(0.111) ; label("Log of the typical first-order rate constant for the approach to the new treated steady state on the logit-DAS28 domain (1/week; bare-name ex2 = 0.111 corresponds to a typical population half-life of log(2) / 0.111 ~= 6.2 weeks)")
    # Table 2 theta_3 = 0.111 /week (95% CI 0.093, 0.129), the typical EX2 (back-transformed from lex2 here for log-normal IIV).

    e_age_base <- 0.672 ; label("Exponent of the AGE / 57 power function applied to BASE on the logit-DAS28 domain (dimensionless)")
    # Table 2 theta_4 = 0.672 (95% CI 0.406, 0.938). Paper Eq 8 worked example: natural-scale DAS28 = 5.4 at age 40 and 5.8 at age 70 (paper Description of the final model paragraph and Results 'Influence of covariates').

    e_conmed_steroid_base <- 0.737 ; label("Multiplicative shift on the logit-domain BASE for visits with intramuscular corticosteroid administration (CONMED_STEROID = 1) relative to no administration (CONMED_STEROID = 0); enters as base * (1 + e_conmed_steroid_base * CONMED_STEROID)")
    # Table 2 theta_5 = 0.737 (95% CI 0.379, 1.095). Paper Eq 8 worked example: BASE = 0.820 (natural-scale DAS28 ~= 6.4) when CONMED_STEROID = 1 vs BASE = 0.472 (natural-scale DAS28 ~= 5.7) when CONMED_STEROID = 0, for an age-57 non-smoker.

    e_conmed_steroid_fu_ex1 <- -0.237 ; label("Multiplicative shift on EX1 for subjects who received any systemic corticosteroid during the 60-week follow-up (CONMED_STEROID_FU = 1) relative to those who did not (CONMED_STEROID_FU = 0); enters as ex1 * (1 + e_conmed_steroid_fu_ex1 * CONMED_STEROID_FU). Negative sign means corticosteroid-rescued patients had a less favourable extent of response.")
    # Table 2 theta_6 = -0.237 (95% CI -0.394, -0.079). Paper Eq 8 worked example: EX1 = -0.977 if CONMED_STEROID_FU = 1, vs EX1 = -1.28 if CONMED_STEROID_FU = 0.

    e_smoke_ex2 <- -0.398 ; label("Multiplicative shift on EX2 for current smokers (SMOKE = 1) relative to never / past smokers (SMOKE = 0); enters as ex2 * (1 + e_smoke_ex2 * SMOKE). Negative sign means current smokers achieved a slower approach to the new treated steady state (population half-life 10.4 weeks vs 6.2 weeks).")
    # Table 2 theta_7 = -0.398 (95% CI -0.651, -0.145).

    # ============================================================
    # Inter-individual variability -- full 3x3 OMEGA BLOCK.
    #
    # Diagonals (variances on the parameter scale):
    #   var(etabase)  = sd_BASE^2  = 0.492^2  = 0.2421
    #   var(etaex1)   = sd_EX1^2   = 0.740^2  = 0.5476
    #   var(etalex2)  = log(CV^2 + 1) = log(0.861^2 + 1) = 0.5547
    #
    # Off-diagonals (covariances; Table 2 entries are correlations, see
    # Cauchy-Schwarz remark in the structural-model header):
    #   cov(etabase, etaex1)   = rho_12 * sd_BASE  * sd_EX1   = -0.376 * 0.492 * 0.740 = -0.1369
    #   cov(etabase, etalex2)  = rho_13 * sd_BASE  * sd_lex2  = -0.137 * 0.492 * sqrt(0.5547) = -0.0502
    #   cov(etaex1,  etalex2)  = rho_23 * sd_EX1   * sd_lex2  =  0.651 * 0.740 * sqrt(0.5547) =  0.3587
    #
    # The block lower-triangular fill order in nlmixr2's `c(...)`
    # convention is column-major: (1,1) (2,1) (2,2) (3,1) (3,2) (3,3).
    # ============================================================

    etabase + etaex1 + etalex2 ~ c(
      0.2421,
      -0.1369,  0.5476,
      -0.0502,  0.3587,  0.5547
    )
    # Diagonals from Table 2 (sd_BASE = 0.492 (95% CI 0.348, 0.636); sd_EX1 = 0.740 (95% CI 0.536, 0.945); CV_EX2 = 86.1% (95% CI 42.1%, 130%)); off-diagonals derived from Table 2 reported correlations rho_12 = -0.376, rho_13 = -0.137, rho_23 = 0.651.

    # ============================================================
    # Residual error -- combined proportional + additive on the
    # predicted logit-DAS28 (paper Equation 3: LGT3 = DAS28_t,j +
    # DAS28_t,j * eps_t,j,1 + eps_t,j,2). Magnitudes from Table 2.
    # ============================================================

    propSd <- 0.184  ; label("Proportional residual SD applied to the predicted logit-DAS28 (dimensionless CV; paper Equation 3 epsilon_1 ~ N(0, sigma_1^2) with sigma_1 = 18.4% CV)")
    # Table 2 proportional sigma_1 = 18.4% CV (95% CI 0%, 37.5%).

    addSd  <- 0.327  ; label("Additive residual SD applied to the predicted logit-DAS28 (logit-DAS28 units; paper Equation 3 epsilon_2 ~ N(0, sigma_2^2) with sigma_2 = 0.327)")
    # Table 2 additive sigma_2 = 0.327 (95% CI 0.249, 0.405) in logit-DAS28 units.
  })

  model({
    # The model has no PK / no ODE / no dosing events -- the only
    # independent variable is time since initiation of triple DMARD
    # therapy (weeks). The output is the predicted DAS28 (both in the
    # logit-transformed domain used for fitting and the natural-scale
    # 0-9.2 domain used for clinical interpretation and figures).

    # Subject-specific structural parameters. BASE and EX1 carry
    # additive normal IIV on the logit-domain linear scale; EX2 carries
    # log-normal IIV on its log scale. Covariate effects multiply the
    # individual parameters per the worked examples below paper
    # Equation 8.
    base_i <- (base + etabase) *
      (1 + e_conmed_steroid_base * CONMED_STEROID) *
      (AGE / 57)^e_age_base

    ex1_i  <- (ex1  + etaex1) *
      (1 + e_conmed_steroid_fu_ex1 * CONMED_STEROID_FU)

    ex2_i  <- exp(lex2 + etalex2) *
      (1 + e_smoke_ex2 * SMOKE)

    # Predicted DAS28 in the logit-transformed domain (paper Eq 7 / Eq
    # 8). DAS28 in the natural-scale 0-9.2 domain is the logit
    # back-transform with upper bound 9.2 (paper Methods: 'A logit
    # transform was employed to constrain the predicted DAS28 values
    # within a plausible range of 0 and 9.2, corresponding to 28
    # swollen and tender joint counts, 100 mm on the visual analogue
    # scale for patient global assessment and 120 mm h-1 ESR').
    das28_logit <- base_i + ex1_i * (1 - exp(-ex2_i * t))
    das28       <- 9.2 / (1 + exp(-das28_logit))

    # Combined proportional + additive residual error on the predicted
    # logit-DAS28 (paper Eq 3). When fitting against observed DAS28
    # data, the user must transform their natural-scale DAS28
    # observations to the logit-DAS28 scale via
    #   das28_logit_obs = log(DAS28_obs / (9.2 - DAS28_obs))
    # before passing them to nlmixr2 (the observation column in the
    # event dataset is das28_logit).
    das28_logit ~ prop(propSd) + add(addSd)
  })
}
