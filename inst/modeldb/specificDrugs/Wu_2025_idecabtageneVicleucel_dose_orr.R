Wu_2025_idecabtageneVicleucel_dose_orr <- function() {
  description <- paste0(
    "Exploratory linear logistic DOSE-response model for the probability ",
    "of overall response (partial response or better, independent review ",
    "committee) after a single infusion of the BCMA-directed CAR T-cell ",
    "therapy idecabtagene vicleucel (ide-cel) in adults with ",
    "triple-class-exposed relapsed/refractory multiple myeloma (Wu 2025; ",
    "ide-cel arm of the phase 3 KarMMa-3 trial, all 225 infused patients). ",
    "The actual infused dose enters linearly (0.0009 per 10^6 CAR+ T cells, ",
    "95% CI -0.005 to 0.007): flat over the tested range. Body weight and ",
    "extramedullary plasmacytoma lower the logit. The dose is uncentred; ",
    "body weight and log soluble BCMA are centred on 75 kg and 144.5 ng/mL ",
    "(an assumption: the paper does not print the centring of this model; ",
    "these are the reference values it states for its dose-exposure model). ",
    "Companion to the exposure-response model ",
    "Wu_2025_idecabtageneVicleucel_orr. No ODE; no random effect and no ",
    "residual error (Bernoulli likelihood)."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S12.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark dose-response regression; no time dimension)",
    dosing = "n/a (no dose events; the actual infused dose enters as the DOSE_IDECEL_MCELLS covariate column)",
    concentration = "prob_orr_central (probability of overall response by independent review, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    DOSE_IDECEL_MCELLS = list(
      description = "Actual infused dose of idecabtagene vicleucel, in millions (10^6) of CAR+ T cells.",
      units = "10^6 CAR+ T cells",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear and uncentred; the intercept-plus-dose term at the KarMMa-3 doses reproduces the observed response rate only without dose centring (vignette).",
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
    DIS_EMD = list(
      description = "1 = extramedullary plasmacytoma present at baseline, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Carried from the final ORR exposure-response model.",
      source_name = "Baseline Extramedullary Plasmacytoma (Yes vs No)"
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
    # Wu 2025 Table S12, 'Parameter Estimates of ORR Dose-Response Model'.
    logit_ref <- 1.4139; label("Overall-response logit for the reference subject at zero dose (unitless logit)") # Table S12 Intercept 1.4139 (SE 1.477)
    e_dose_idecel_logit <- 9.00e-04; label("Slope of the overall-response logit per 10^6 CAR+ T cells of actual dose (1/(10^6 cells))") # Table S12 'Actual dose (x10^6 cells)' 9.00E-04 (SE 0.003)
    e_sexm_logit <- 0.3363; label("Shift in the overall-response logit for male vs female sex (unitless logit)") # Table S12 'Sex (Male vs Female)' 0.3363 (SE 0.402)
    e_wt_logit <- -0.0244; label("Slope of the overall-response logit per kg of body weight above 75 kg (1/kg)") # Table S12 'Baseline Body Weight (kg)' -0.0244 (SE 0.011)
    e_form_idecel_b_logit <- -0.0118; label("Shift in the overall-response logit for Formulation B vs A (unitless logit)") # Table S12 'Product Attributes (Formulation B vs A)' -0.0118 (SE 0.409)
    e_sbcma_logit <- -0.2853; label("Slope of the overall-response logit per unit log(soluble BCMA / 144.5 ng/mL) (unitless logit)") # Table S12 'Log(Baseline Soluble BCMA)' -0.2853 (SE 0.176)
    e_line_4l_logit <- -0.2228; label("Shift in the overall-response logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S12 'No. of prior MM Regimens (3 vs 2)' -0.2228 (SE 0.427)
    e_line_5l_logit <- 0.1596; label("Shift in the overall-response logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S12 'No. of prior MM Regimens (4 vs 2)' 0.1596 (SE 0.466)
    e_emd_logit <- -1.1973; label("Shift in the overall-response logit for baseline extramedullary plasmacytoma (unitless logit)") # Table S12 'Baseline Extramedullary Plasmacytoma (Yes vs No)' -1.1973 (SE 0.377)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr_central <- logit_ref +
      e_dose_idecel_logit * DOSE_IDECEL_MCELLS +
      e_sexm_logit * (1 - SEXF) +
      e_wt_logit * (WT - 75) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_sbcma_logit * log(SBCMA / 144.5) +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L +
      e_emd_logit * DIS_EMD

    prob_orr_central <- expit(logit_orr_central)

    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
