Wu_2025_idecabtageneVicleucel_dose_crs_corticosteroid <- function() {
  description <- paste0(
    "Exploratory linear logistic DOSE-response model for the probability ",
    "of cytokine release syndrome requiring corticosteroid treatment (sCRS) ",
    "after a single infusion of the BCMA-directed CAR T-cell therapy ",
    "idecabtagene vicleucel (ide-cel) in adults with triple-class-exposed ",
    "relapsed/refractory multiple myeloma (Wu 2025; ide-cel arm of the ",
    "phase 3 KarMMa-3 trial, all 225 infused patients). The actual infused ",
    "dose enters linearly (0.0052 per 10^6 CAR+ T cells, 95% CI -0.0011 to ",
    "0.0122, not significant); body weight lowers the logit. The dose is ",
    "uncentred; body weight and log soluble BCMA are centred on 75 kg and ",
    "144.5 ng/mL (an assumption: the paper does not print the centring of ",
    "this model; these are the reference values it states for its ",
    "dose-exposure model). Companion to the exposure-response model ",
    "Wu_2025_idecabtageneVicleucel_crs_corticosteroid. No ODE; no random ",
    "effect and no residual error (Bernoulli likelihood)."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S15.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark dose-response regression; no time dimension)",
    dosing = "n/a (no dose events; the actual infused dose enters as the DOSE_IDECEL_MCELLS covariate column)",
    concentration = "prob_crs_corticosteroid (probability of CRS requiring corticosteroids, 0-1; also logit_crs_corticosteroid)"
  )

  covariateData <- list(
    DOSE_IDECEL_MCELLS = list(
      description = "Actual infused dose of idecabtagene vicleucel, in millions (10^6) of CAR+ T cells.",
      units = "10^6 CAR+ T cells",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear and uncentred (vignette).",
      source_name = "Actual dose (x10^6 cells)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "Male vs female with female as reference; applied to (1 - SEXF).",
      source_name = "Sex (Male vs Female)"
    ),
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear, centred on 75 kg (assumed; see description).",
      source_name = "Baseline Body Weight (kg)"
    ),
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-specified covariate.",
      source_name = "Product Attributes (Formulation B vs A)"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen, measured on the day of the ide-cel infusion.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Natural log, centred on 144.5 ng/mL (assumed; see description).",
      source_name = "Log(Baseline Soluble BCMA) (ng/mL)"
    ),
    LINE_4L = list(
      description = "1 = exactly three prior multiple-myeloma regimens, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference level = 2 prior regimens.",
      source_name = "No. of prior MM Regimens (3 vs 2)"
    ),
    LINE_5L = list(
      description = "1 = exactly four prior multiple-myeloma regimens, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference level = 2 prior regimens.",
      source_name = "No. of prior MM Regimens (4 vs 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 225L,
    n_studies = 1L,
    n_observations = "225 binary outcomes (all infused patients of the ide-cel arm)",
    age_range = "30-81 years (median 63; Table S3, N = 220 exposure-evaluable subset)",
    weight_range = "41.1-144 kg (median 83; Table S3)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens with progression within 60 days of the last regimen",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted)",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, data cutoff 1 March 2022; exploratory dose-response analysis on all 225 infused patients (Supporting Information D)."
  )

  ini({
    # Wu 2025 Table S15, 'Parameter Estimates of sCRS Dose-Response Model'.
    logit_ref <- -3.5294; label("sCRS logit for the reference subject at zero dose (unitless logit)") # Table S15 Intercept -3.5294 (SE 1.5646)
    e_dose_idecel_logit <- 0.0052; label("Slope of the sCRS logit per 10^6 CAR+ T cells of actual dose (1/(10^6 cells))") # Table S15 'Actual dose (x10^6 cells)' 0.0052 (SE 0.0034)
    e_sexm_logit <- 0.3207; label("Shift in the sCRS logit for male vs female sex (unitless logit)") # Table S15 'Sex (Male vs Female)' 0.3207 (SE 0.3735)
    e_wt_logit <- -0.0248; label("Slope of the sCRS logit per kg of body weight above 75 kg (1/kg)") # Table S15 'Baseline Body Weight (kg)' -0.0248 (SE 0.011)
    e_form_idecel_b_logit <- 0.218; label("Shift in the sCRS logit for Formulation B vs A (unitless logit)") # Table S15 'Product Attributes (Formulation B vs A)' 0.218 (SE 0.3617)
    e_sbcma_logit <- 0.2334; label("Slope of the sCRS logit per unit log(soluble BCMA / 144.5 ng/mL) (unitless logit)") # Table S15 'Log(Baseline Soluble BCMA)' 0.2334 (SE 0.1534)
    e_line_4l_logit <- 0.3502; label("Shift in the sCRS logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S15 'No. of prior MM Regimens (3 vs 2)' 0.3502 (SE 0.383)
    e_line_5l_logit <- 0.0424; label("Shift in the sCRS logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S15 'No. of prior MM Regimens (4 vs 2)' 0.0424 (SE 0.4082)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_corticosteroid <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_crs_corticosteroid <- logit_ref +
      e_dose_idecel_logit * DOSE_IDECEL_MCELLS +
      e_sexm_logit * (1 - SEXF) +
      e_wt_logit * (WT - 75) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_sbcma_logit * log(SBCMA / 144.5) +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L

    prob_crs_corticosteroid <- expit(logit_crs_corticosteroid)

    prob_crs_corticosteroid ~ add(addSd_prob_crs_corticosteroid)
  })
}
