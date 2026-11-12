Li_2025_mosunetuzumab_cr_cov <- function() {
  description <- paste(
    "Exploratory covariate version of the Emax logistic-regression model",
    "for INVESTIGATOR-ASSESSED COMPLETE RESPONSE versus the two-cycle",
    "mosunetuzumab AUC0-42 (AUC_MOSUNETUZUMAB, day*ug/mL) in adults with",
    "relapsed/refractory follicular lymphoma after >= 2 prior therapies (Li",
    "2025, study GO29781, n = 159; Figure S3). Residual baseline anti-CD20",
    "antibody (CP_ACD20_UGML) lowers the maximal response and baseline",
    "tumour burden (TUMSZ, through its square root) raises E50. Li 2025",
    "prints NO coefficients or equations for this model: every value was",
    "back-solved from the Max CR and E90 annotations and the fitted curves",
    "of Figure S3, with the functional forms ASSUMED (see the ini() block",
    "and the vignette). The base model Li_2025_mosunetuzumab_cr is the one",
    "the paper uses for its dose-regimen predictions."
  )
  reference <- paste(
    "Li C-C, Bender B, Wilkins J, Li F, Turner DC, Wang B, Deng R,",
    "Vadhavkar S, Li Z, Kwan A, Huang H, Peng K, Penuel E, Huw L-Y, Chanu P,",
    "Li C, Yin S, Wei MC. A Novel Step-Up Dosage Regimen for Enhancing the",
    "Benefit-to-Risk Ratio of Mosunetuzumab in Relapsed or Refractory",
    "Follicular Lymphoma. Clin Pharmacol Ther. 2025;117(2):465-474.",
    "doi:10.1002/cpt.3445. Coefficients back-solved from Figure S3."
  )
  vignette <- "Li_2025_mosunetuzumab"
  units <- list(
    time = "n/a (static two-cycle landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_MOSUNETUZUMAB data column)",
    concentration = "prob_cr_investigator (probability of complete response, investigator-assessed, 0-1)"
  )

  covariateData <- list(
    AUC_MOSUNETUZUMAB = list(
      description = "Model-predicted mosunetuzumab serum AUC from day 0 to day 42 (the first two 21-day cycles)",
      units = "day*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Same exposure metric as Li_2025_mosunetuzumab_cr: the empirical-",
        "Bayes AUC0-42 from the Bender 2024 popPK model",
        "(Bender_2024_mosunetuzumab, `auc` state at day 42). Observed range",
        "0.4-251 day*ug/mL (Table S10)."
      ),
      source_name = "AUC0-42"
    ),
    CP_ACD20_UGML = list(
      description = "Composite baseline anti-CD20 drug concentration: the maximum of the predose rituximab and obinutuzumab concentrations",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Li 2025 Table S8 'Baseline anti-CD20 concentration (LOG_aCD20,",
        "ng/mL, natural log transformed)', tested on Emax and E50. The same",
        "composite as Bender 2024 (max of rituximab and obinutuzumab).",
        "Below-LOQ values were set to the 0.5 ug/mL (500 ng/mL) rituximab",
        "LOQ in the analysis data (Table S10 median 500 ng/mL), which is the",
        "reference here; supply 0.5, not 0, for a patient with no",
        "quantifiable residual antibody. Figure S3 shows only 0.5 and 12.1",
        "ug/mL (the median of the ~35 percent with measurable levels), so",
        "the shape between and beyond those two values is assumed."
      ),
      source_name = "LOG_aCD20"
    ),
    TUMSZ = list(
      description = "Baseline tumour burden, sum of the products of perpendicular diameters (SPD, investigator-assessed)",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Li 2025 Table S8 'Square root of baseline SPD (SPDS)', in mm,",
        "tested on E50. model() takes sqrt(TUMSZ) and references it to",
        "55.4 mm, the efficacy-population median (Table S10); Figure S3",
        "panels sit at 24.5, 55.4 and 104 mm (5th, 50th and 95th",
        "percentiles). Same mm^2 SPD convention as Bender_2024_mosunetuzumab."
      ),
      source_name = "SPDS"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 159,
    n_studies = 1,
    age_range = "27-90 years (median 60)",
    weight_range = "42-140 kg (median 80.4)",
    sex_female_pct = 38.4,
    race_ethnicity = c(
      White = 76.1,
      Asian = 15.1,
      Black = 5.03,
      `American Indian or Alaska Native` = 0.629,
      Unknown = 3.14
    ),
    disease_state = "Relapsed/refractory follicular lymphoma after >= 2 prior systemic therapies (median 3)",
    dose_range = "0.2-2.8 mg IV q3w fixed dosing (Group A, n = 7) and Cycle 1 step-up dosing from 0.4/1.0/2.8 to 1/2/60/30 mg IV q3w (Group B, n = 152; cohort B11 1/2/60/30 mg n = 90)",
    regions = "Multinational",
    notes = paste(
      "Efficacy exposure-response population of GO29781 (NCT02500407);",
      "baseline characteristics from Li 2025 Table S10: median sqrt(SPD)",
      "55.4 mm (range 15-149); baseline rituximab median 500 ng/mL (range",
      "500-102000), obinutuzumab median below LOQ (maximum 305000 ng/mL).",
      "Covariates were selected by stepwise forward addition (dOFV >= 3.84)",
      "and backward deletion (dOFV <= 6.63) from the predefined list of",
      "Table S8."
    )
  )

  ini({
    # ======================================================================
    # Li 2025 prints no estimates for this model. Figure S3 shows six panels
    # (aCD20 = 500 or 12100 ng/mL x SPDS = 24.5, 55.4 or 104 mm), each
    # annotated with the maximal CR probability (76.2 percent at 500 ng/mL,
    # 52.3 percent at 12100 ng/mL, the same down every column) and with E90,
    # the AUC0-42 that gives 90 percent of that maximal PROBABILITY
    # (7.07 / 42.9 / 173 and 11.7 / 68.2 / 283 day*ug/mL). The E90
    # definition was confirmed on the base models: Table S7 reproduces the
    # printed base E90 values 92.7 (CR) and 38.4 (OR) day*ug/mL.
    #
    # Back-solve:
    # 1. Max CR = expit(beta0 + Emax) fixes Emax at each aCD20 level once
    #    beta0 is known; each E90 then fixes that panel's E50.
    # 2. beta0 was chosen by least squares of the six digitised fitted
    #    curves (pixel-traced, 1900 points) under those constraints:
    #    beta0 = -6.23 (RMSE 0.021 in probability; the fit is shallow
    #    between -5.5 and -7).
    # 3. The six implied E50 values were regressed on log(SPDS / 55.4) and
    #    log(aCD20 / 500): exponent 2.207 on SPDS (the per-column pairwise
    #    exponents agree to within 2.16-2.26), exponent 0.0061 on aCD20 (E50
    #    is practically unchanged by aCD20), residual SD 0.017 on log E50.
    # 4. The rounded parameters reproduce all twelve printed annotations:
    #    Max CR exactly, E90 within 2.4 percent (worst 69.8 vs 68.2).
    #
    # ASSUMED forms: aCD20 enters linearly in LOG_aCD20 (Table S8 names the
    # natural-log-transformed covariate) -- additively on Emax and
    # log-linearly (a power) on E50; SPDS enters as a power on E50. Only the
    # two plotted aCD20 levels are constrained by the source.
    # ======================================================================
    logite0 <- -6.23; label("Logit of the probability of complete response at AUC0-42 = 0, beta0 (unitless logit)") # back-solved from Figure S3 (least squares of the six digitised curves under the printed Max CR and E90 constraints)
    emax <- 7.394; label("Maximal logit shift produced by AUC0-42 at the reference aCD20 of 0.5 ug/mL, Emax (unitless logit)") # back-solved: logit(0.762) - beta0, Figure S3 'Max CR = 76.2%'
    e_acd20_emax <- -0.3363; label("Change in Emax per unit increase in the natural log of the baseline anti-CD20 concentration (unitless logit)") # back-solved: (logit(0.523) - logit(0.762)) / log(12100/500), Figure S3 'Max CR = 52.3%'
    lec50 <- log(2.346); label("AUC0-42 producing half of Emax at SPDS = 55.4 mm and aCD20 = 0.5 ug/mL, E50 (day*ug/mL)") # back-solved from the Figure S3 E90 annotations
    e_tumsz_e50 <- 2.207; label("Power exponent of sqrt(SPD)/55.4 mm on E50 (unitless)") # back-solved from the Figure S3 E90 annotations (SPDS 24.5 / 55.4 / 104 mm)
    e_acd20_e50 <- 0.0061; label("Power exponent of aCD20/0.5 ug/mL on E50 (unitless)") # back-solved from the Figure S3 E90 annotations (aCD20 500 vs 12100 ng/mL)

    addSd_prob_cr_investigator <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Covariates on the Emax logistic regression (Li 2025 Methods equation;
    # covariate forms assumed, see ini()). LOG_aCD20 is the natural log in
    # ng/mL, so log(aCD20) - log(500) = log(CP_ACD20_UGML / 0.5) in any unit.
    lacd20_rel <- log(CP_ACD20_UGML / 0.5)
    spds <- sqrt(TUMSZ)

    emax_i <- emax + e_acd20_emax * lacd20_rel
    ec50 <- exp(lec50) * (spds / 55.4)^e_tumsz_e50 * exp(e_acd20_e50 * lacd20_rel)

    logit_cr <- logite0 + emax_i * AUC_MOSUNETUZUMAB / (ec50 + AUC_MOSUNETUZUMAB)

    prob_cr_investigator <- expit(logit_cr)

    prob_cr_investigator ~ add(addSd_prob_cr_investigator)
  })
}
