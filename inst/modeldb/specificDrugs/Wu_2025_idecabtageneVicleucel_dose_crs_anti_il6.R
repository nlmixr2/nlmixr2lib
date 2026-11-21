Wu_2025_idecabtageneVicleucel_dose_crs_anti_il6 <- function() {
  description <- paste0(
    "Exploratory linear logistic DOSE-response model for the probability ",
    "of cytokine release syndrome requiring anti-IL-6 treatment ",
    "(tocilizumab or siltuximab; tCRS) after a single infusion of the ",
    "BCMA-directed CAR T-cell therapy idecabtagene vicleucel (ide-cel) in ",
    "adults with triple-class-exposed relapsed/refractory multiple myeloma ",
    "(Wu 2025; ide-cel arm of the phase 3 KarMMa-3 trial, all 225 infused ",
    "patients). The actual infused dose enters linearly (0.0038 per 10^6 ",
    "CAR+ T cells, 95% CI -0.0018 to 0.0095, not significant). The dose is ",
    "uncentred; body weight and log soluble BCMA are centred on 75 kg and ",
    "144.5 ng/mL (an assumption: the paper does not print the centring of ",
    "this model; these are the reference values it states for its ",
    "dose-exposure model). Companion to the exposure-response model ",
    "Wu_2025_idecabtageneVicleucel_crs_anti_il6. No ODE; no random effect ",
    "and no residual error (Bernoulli likelihood)."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S14.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark dose-response regression; no time dimension)",
    dosing = "n/a (no dose events; the actual infused dose enters as the DOSE_IDECEL_MCELLS covariate column)",
    concentration = "prob_crs_anti_il6 (probability of CRS requiring tocilizumab or siltuximab, 0-1; also logit_crs_anti_il6)"
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
    ),
    ADA_POS = list(
      description = "1 = anti-CAR antidrug antibody positive at at least one post-infusion visit, 0 = negative.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "POST-BASELINE covariate, carried from the final tCRS exposure-response model.",
      source_name = "Post-infusion ADA Positive (Yes vs No)"
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
    # Wu 2025 Table S14, 'Parameter Estimates of tCRS Dose-Response Model'.
    logit_ref <- -0.0746; label("tCRS logit for the reference subject at zero dose (unitless logit)") # Table S14 Intercept -0.0746 (SE 1.2906)
    e_dose_idecel_logit <- 0.0038; label("Slope of the tCRS logit per 10^6 CAR+ T cells of actual dose (1/(10^6 cells))") # Table S14 'Actual dose (x10^6 cells)' 0.0038 (SE 0.0028)
    e_sexm_logit <- -0.3518; label("Shift in the tCRS logit for male vs female sex (unitless logit)") # Table S14 'Sex (Male vs Female)' -0.3518 (SE 0.3664)
    e_wt_logit <- -0.0186; label("Slope of the tCRS logit per kg of body weight above 75 kg (1/kg)") # Table S14 'Baseline Body Weight (kg)' -0.0186 (SE 0.0097)
    e_form_idecel_b_logit <- -0.0259; label("Shift in the tCRS logit for Formulation B vs A (unitless logit)") # Table S14 'Product Attributes (Formulation B vs A)' -0.0259 (SE 0.354)
    e_sbcma_logit <- 0.0983; label("Slope of the tCRS logit per unit log(soluble BCMA / 144.5 ng/mL) (unitless logit)") # Table S14 'Log(Baseline Soluble BCMA)' 0.0983 (SE 0.1409)
    e_line_4l_logit <- -0.2404; label("Shift in the tCRS logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S14 'No. of prior MM Regimens (3 vs 2)' -0.2404 (SE 0.3683)
    e_line_5l_logit <- 0.2182; label("Shift in the tCRS logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S14 'No. of prior MM Regimens (4 vs 2)' 0.2182 (SE 0.4065)
    e_ada_pos_logit <- -0.3467; label("Shift in the tCRS logit for post-infusion ADA positivity (unitless logit)") # Table S14 'Post-infusion ADA Positive (Yes vs No)' -0.3467 (SE 0.3322)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_anti_il6 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_crs_anti_il6 <- logit_ref +
      e_dose_idecel_logit * DOSE_IDECEL_MCELLS +
      e_sexm_logit * (1 - SEXF) +
      e_wt_logit * (WT - 75) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_sbcma_logit * log(SBCMA / 144.5) +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L +
      e_ada_pos_logit * ADA_POS

    prob_crs_anti_il6 <- expit(logit_crs_anti_il6)

    prob_crs_anti_il6 ~ add(addSd_prob_crs_anti_il6)
  })
}
