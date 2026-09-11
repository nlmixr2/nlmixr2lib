Sun_2025_maribavir_thrombocytopenia <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-SAFETY model for ",
    "treatment-emergent thrombocytopenia as a function of the maribavir ",
    "AUC on the DAY OF THE EVENT, in allogeneic hematopoietic cell ",
    "transplant (HCT) recipients with a first asymptomatic ",
    "cytomegalovirus infection treated with maribavir 400 mg orally twice ",
    "daily (Sun 2025, n = 238, the AURORA maribavir arm). Note that this ",
    "is the DAY-OF-EVENT exposure metric (AUCday), not the steady-state ",
    "metric: Sun 2025's abstract reports that STEADY-STATE exposures were ",
    "not significantly associated with any adverse event other than ",
    "nausea and vomiting, whereas all fourteen day-of-event relationships ",
    "in Figure S3 are significant. A day-of-event AUC is partly a ",
    "consequence of when in the treatment course the event happened, so ",
    "it is the weaker causal claim of the two. This is the only one of ",
    "the fourteen models carrying a CONTINUOUS risk factor, the ",
    "transplanted nucleated cell number, which enters centred at the ",
    "source's 400-cell reference and scaled per 1000-cell increment. The ",
    "exposure slope is +0.0480 per 10 ug*h/mL. There is no PK layer and ",
    "no ODE: the exposure metric is supplied as a data column, derived in ",
    "the source analysis from the companion population PK model packaged ",
    "as Sun_2025_maribavir. No between-subject random effect and no ",
    "residual error are estimated (Bernoulli likelihood). Fifteen ",
    "companion exposure-response models in the Sun_2025_maribavir_* ",
    "family."
  )
  reference <- paste(
    "Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y, Bhattacharya I.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Maribavir",
    "in Transplant Recipients With First Episode or Refractory Cytomegalovirus.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1346-1356.",
    "doi:10.1002/psp4.70054.",
    "Logistic-regression coefficients in the Figure S3 parameter tables, Supporting Information file s002; the corresponding odds ratios are also tabulated in Table S4, file s001.",
    sep = " "
  )
  vignette <- "Sun_2025_maribavir"
  units <- list(
    time          = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_MBV_DAY covariate column)",
    concentration = "prob_thrombocytopenia (probability of treatment-emergent thrombocytopenia, 0-1; also logit_thrombocytopenia)"
  )

  covariateData <- list(
    AUC_MBV_DAY = list(
      description        = "Individual maribavir area under the plasma concentration-time curve over the calendar day on which the adverse event occurred. Supplied as data: this model has no PK layer, and the source analysis used individual predictions from the companion maribavir population PK model.",
      units              = "ug*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sun 2025 estimates the coefficient PER INCREMENT OF 10 ug*h/mL, so",
        "model() divides this column by 10. This is the DAY-OF-EVENT exposure, not",
        "the steady-state exposure: where a subject had more than one event of the",
        "same type, the source retained the exposure on the day of the MOST SEVERE",
        "event (Methods 2.4), and subjects without the event contribute the",
        "exposure on the corresponding on-treatment day. Do not substitute",
        "AUC_MBV_SS -- that column drives the two EFFICACY models and is a",
        "different quantity. The Figure S3 quartile boundaries run from about",
        "[0, 116) to [313, 951] ug*h/mL, so the observed range spans roughly an",
        "order of magnitude.",
        "In this model the exposure enters with slope 0.0480 per 10 ug*h/mL (odds",
        "ratio 1.05 (1.02-1.08), p = 0.00367)."
      ),
      source_name        = "AUCday of maribavir (increment of 10 h.ug/mL)"
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "108 of 238 subjects (45.4%) female, 130 (54.6%) male. Sun 2025 fits the",
        "effect as 'Male vs. female' with FEMALE as the reference, which is the",
        "opposite polarity from the SEXF canonical. model() therefore multiplies",
        "the published male coefficient by (1 - SEXF) rather than by SEXF; the",
        "published estimate and odds ratio are carried unchanged and the inversion",
        "happens at the call site, where it is visible.",
        "Coefficient +1.11 for MALE sex, with FEMALE as the source's reference;",
        "model() multiplies by (1 - SEXF) to match the SEXF canonical's polarity.",
        "Almost identical in magnitude to the acute-GvHD model's male effect."
      ),
      source_name        = "Sex, Male (reference Female)"
    ),
    CONMED_ANTILYMPHOCYTE = list(
      description        = "Antilymphocyte agent use (class-level composite).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no antilymphocyte agent use)",
      notes              = paste(
        "108 of 238 subjects (45.4%) positive. Overlaps with but is NOT the same",
        "column as the HCT_TCD_* set: those record the graft-preparation protocol",
        "for the current transplant, this records antilymphocyte drug exposure as",
        "a concomitant-medication fact, and their denominators differ (108 here",
        "against 83 across the in-vivo HCT_TCD_* levels). The source reports only",
        "the class, so the composite cannot be decomposed to named agents.",
        "Coefficient -1.30, significant."
      ),
      source_name        = "Antilymphocyte use"
    ),
    HCT_NUCCELL = list(
      description        = "Nucleated cell number transplanted for the current haematopoietic cell transplant.",
      units              = "cells (source scale; see notes -- the absolute scale is unstated by Sun 2025)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sun 2025 Figure S3 panel m prints '400 cells  Reference' and estimates",
        "the coefficient per 'Increment of 1000 cells', so model() forms",
        "(HCT_NUCCELL - 400) / 1000. The ABSOLUTE SCALE of the column is unstated",
        "in the source and read literally is not a plausible graft cell dose (a",
        "real allogeneic graft holds on the order of 10^8-10^10 nucleated cells),",
        "so the source's own numbers are carried verbatim rather than a multiplier",
        "being guessed. The model is self-consistent and reproduces the published",
        "odds ratio whatever the underlying unit turns out to be, but a user",
        "supplying a real cell dose must first rescale it onto the source's scale.",
        "See the vignette Assumptions and deviations section.",
        "Coefficient +0.00689 per 1000-cell increment above the 400-cell",
        "reference, so model() forms (HCT_NUCCELL - 400)/1000. The published odds",
        "ratio 1.01 per increment is exp(0.00689) = 1.0069 rounded, which confirms",
        "the increment basis; the centring at 400 follows the figure's own",
        "'400 cells  Reference' row and shifts only the intercept."
      ),
      source_name        = "Nucleated cell number for current HCT (per increment of 1000 cells), reference 400 cells"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 238L,
    n_studies      = 1L,
    n_observations = "238 binary thrombocytopenia records, one per patient",
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
    # Sun 2025 Figure S3 panel m (Supporting Information file s002).
    #
    #   logit(p) = intercept + slope * (AUC_MBV_DAY / 10)
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
    logit_ref <- -3.56; label("Logit of the thrombocytopenia probability at zero maribavir exposure with all risk factors at their reference level (unitless logit)")  # Sun 2025 Figure S3 panel m (Supporting Information file s002), Intercept -3.56 (SE 0.654), p < 0.001

    # ----- Exposure effect -----
    e_auc_logit <- 0.0480; label("Log-odds of thrombocytopenia per 10 ug*h/mL increase in maribavir AUC on the day of the event (unitless logit)")  # Sun 2025 Figure S3 panel m (Supporting Information file s002), Estimate 0.0480 (SE 0.0165), odds ratio 1.05 (1.02-1.08), p = 0.00367

    # ----- Risk-factor effects on the logit -----
    e_sexm_logit <- 1.11; label("Log-odds of thrombocytopenia for male sex versus female (so the term is multiplied by 1 - SEXF) (unitless logit)")  # Sun 2025 Figure S3 panel m (Supporting Information file s002), Sex, Male (reference Female) Estimate 1.11 (SE 0.498), odds ratio 3.03 (1.14-8.04), p = 0.02605
    e_conmed_antilymphocyte_logit <- -1.30; label("Log-odds of thrombocytopenia for antilymphocyte agent use versus none (unitless logit)")  # Sun 2025 Figure S3 panel m (Supporting Information file s002), Antilymphocyte use Estimate -1.30 (SE 0.505), odds ratio 0.271 (0.101-0.730), p = 0.00976
    e_hct_nuccell_logit <- 0.00689; label("Log-odds of thrombocytopenia for each increment of 1000 transplanted nucleated cells above the 400-cell reference (unitless logit)")  # Sun 2025 Figure S3 panel m (Supporting Information file s002), Nucleated cell number for current HCT (per increment of 1000 cells), reference 400 cells Estimate 0.00689 (SE 0.00301), odds ratio 1.01 (1.00-1.01), p = 0.02195

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli likelihood
    # has no sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_thrombocytopenia <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Sun 2025 estimates every coefficient per 10 ug*h/mL increment of the
    # DAY-OF-EVENT AUC, so the exposure column is divided by 10 here.
    auc10 <- AUC_MBV_DAY / 10

    # ----- Linear predictor -----
    logit_thrombocytopenia <- logit_ref +
      e_auc_logit * auc10 +
      e_sexm_logit * (1 - SEXF) +
      e_conmed_antilymphocyte_logit * CONMED_ANTILYMPHOCYTE +
      e_hct_nuccell_logit * ((HCT_NUCCELL - 400) / 1000)

    prob_thrombocytopenia <- expit(logit_thrombocytopenia)

    # ----- Observation -----
    prob_thrombocytopenia ~ add(addSd_prob_thrombocytopenia)
  })
}
