Wu_2025_idecabtageneVicleucel_dose_pfs <- function() {
  description <- paste0(
    "Exploratory Cox proportional-hazards DOSE-response model for ",
    "progression-free survival (independent review committee) after a ",
    "single infusion of the BCMA-directed CAR T-cell therapy idecabtagene ",
    "vicleucel (ide-cel) in adults with triple-class-exposed ",
    "relapsed/refractory multiple myeloma (Wu 2025; ide-cel arm of the ",
    "phase 3 KarMMa-3 trial, all 225 infused patients). The actual infused ",
    "dose enters linearly (log hazard ratio -0.003 per 10^6 CAR+ T cells, ",
    "95% CI -0.007 to 0, not significant). The model returns the RELATIVE ",
    "hazard (hr, and its log lhr) against a reference subject (female, ",
    "Formulation A, two prior regimens, ADA-negative, not triple-class ",
    "refractory, no extramedullary plasmacytoma, no steroid in the last ",
    "prior regimen, dose 450 x 10^6 cells, body weight 75 kg, soluble BCMA ",
    "144.5 ng/mL, ferritin 163 ug/L, TNF-alpha 27 pg/mL); the reference ",
    "values only fix which subject has hr = 1 and do not change any hazard ",
    "ratio. Baseline TNF-alpha enters on the log2 scale (see ini()). No ",
    "baseline hazard is encoded: a Cox regression leaves h0(t) unspecified. ",
    "Companion to the exposure-response model ",
    "Wu_2025_idecabtageneVicleucel_pfs. No ODE, no random effect, no ",
    "residual error."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S13.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (semiparametric Cox relative hazard; the baseline hazard and therefore the time scale are unspecified by the fit)",
    dosing = "n/a (no dose events; the actual infused dose enters as the DOSE_IDECEL_MCELLS covariate column)",
    concentration = "hr (relative hazard of progression or death vs the reference subject, unitless; also lhr, the log relative hazard)"
  )

  covariateData <- list(
    DOSE_IDECEL_MCELLS = list(
      description = "Actual infused dose of idecabtagene vicleucel, in millions (10^6) of CAR+ T cells.",
      units = "10^6 CAR+ T cells",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear; referenced to 450 x 10^6 cells (the top of the KarMMa-3 target range) so that hr = 1 at a typical dose.",
      source_name = "Actual dose (x10^6 cells)"
    ),
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear, referenced to 75 kg.",
      source_name = "Baseline Body Weight (kg)"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen, measured on the day of the ide-cel infusion.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, referenced to the median 144.5 ng/mL.",
      source_name = "Log2(Baseline Soluble BCMA) (ng/mL)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "Male vs female with female as reference; applied to (1 - SEXF).",
      source_name = "Sex (Male vs Female)"
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
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-specified covariate.",
      source_name = "Product Attributes (Formulation B vs A)"
    ),
    ADA_POS = list(
      description = "1 = anti-CAR antidrug antibody positive at at least one post-infusion visit, 0 = negative.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "POST-BASELINE covariate; the authors caution that its protective coefficient is confounded by the longer follow-up of patients with longer PFS.",
      source_name = "Post-infusion ADA Positive (Yes vs No)"
    ),
    REFRACTORY_TRIPLE_CLASS = list(
      description = "1 = refractory to an immunomodulatory agent, a proteasome inhibitor and an anti-CD38 antibody, 0 = not triple-class refractory.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Carried from the final PFS exposure-response model.",
      source_name = "Prior Triple Refractory (Yes vs No)"
    ),
    DIS_EMD = list(
      description = "1 = extramedullary plasmacytoma present at baseline, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-selected covariate.",
      source_name = "Baseline Extramedullary Plasmacytoma (Yes vs No)"
    ),
    FERRITIN_BL = list(
      description = "Baseline serum ferritin.",
      units = "ug/L",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, referenced to the median 163 ug/L.",
      source_name = "Log2(Baseline Ferritin) (ug/L)"
    ),
    LAST_PRIOR_CORTICOSTEROID = list(
      description = "1 = the last prior multiple-myeloma regimen before ide-cel included a corticosteroid, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-selected covariate.",
      source_name = "Last Prior Med: Steroids (Yes vs No)"
    ),
    TNFA = list(
      description = "Baseline serum tumour necrosis factor alpha, measured on the day of the ide-cel infusion.",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, referenced to the median 27 pg/mL. Table S13 labels the row 'Baseline TNF-alpha (pg/mL)' without a transform; read as log2 (see ini()).",
      source_name = "Baseline TNF-alpha (pg/mL)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 225L,
    n_studies = 1L,
    n_observations = "225 right-censored progression-free-survival times (all infused patients of the ide-cel arm)",
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
    # Wu 2025 Table S13, 'Parameter Estimates of PFS Dose-Response Model'
    # (log hazard ratios).
    e_dose_idecel_haz <- -0.003; label("Log hazard ratio per 10^6 CAR+ T cells of actual dose (1/(10^6 cells))") # Table S13 'Actual dose (x10^6 cells)' -0.003 (SE 0.002)
    e_wt_haz <- 0.008; label("Log hazard ratio per kg of body weight (1/kg)") # Table S13 'Baseline Body Weight (kg)' 0.008 (SE 0.006)
    e_sbcma_haz <- 0.183; label("Log hazard ratio per doubling of baseline soluble BCMA (unitless)") # Table S13 'Log2(Baseline Soluble BCMA)' 0.183 (SE 0.072)
    e_sexm_haz <- -0.088; label("Log hazard ratio for male vs female sex (unitless)") # Table S13 'Sex (Male vs Female)' -0.088 (SE 0.225)
    e_line_4l_haz <- 0.135; label("Log hazard ratio for 3 vs 2 prior MM regimens (unitless)") # Table S13 'No. of prior MM Regimens (3 vs 2)' 0.135 (SE 0.226)
    e_line_5l_haz <- 0.161; label("Log hazard ratio for 4 vs 2 prior MM regimens (unitless)") # Table S13 'No. of prior MM Regimens (4 vs 2)' 0.161 (SE 0.242)
    e_form_idecel_b_haz <- -0.466; label("Log hazard ratio for Formulation B vs A (unitless)") # Table S13 'Product Attributes (Formulation B vs A)' -0.466 (SE 0.273)
    e_ada_pos_haz <- -0.879; label("Log hazard ratio for post-infusion ADA positivity (unitless)") # Table S13 'Post-infusion ADA Positive (Yes vs No)' -0.879 (SE 0.204)
    e_refractory_triple_class_haz <- 0.414; label("Log hazard ratio for triple-class refractory disease (unitless)") # Table S13 'Prior Triple Refractory (Yes vs No)' 0.414 (SE 0.217)
    e_emd_haz <- 0.527; label("Log hazard ratio for baseline extramedullary plasmacytoma (unitless)") # Table S13 'Baseline Extramedullary Plasmacytoma (Yes vs No)' 0.527 (SE 0.204)
    e_ferritin_bl_haz <- 0.088; label("Log hazard ratio per doubling of baseline ferritin (unitless)") # Table S13 'Log2(Baseline Ferritin)' 0.088 (SE 0.06)
    e_last_prior_corticosteroid_haz <- 0.385; label("Log hazard ratio for a corticosteroid in the last prior regimen (unitless)") # Table S13 'Last Prior Med: Steroids (Yes vs No)' 0.385 (SE 0.443)
    # Table S13 prints this row as 'Baseline TNF-alpha (pg/mL)' with no
    # transform. Per pg/mL it would give a hazard ratio of exp(0.406 * 151)
    # = 4e26 across the observed 8.8-160 pg/mL range, so it is read as
    # log2, the transform of the same covariate in Table S7 (0.33 per
    # doubling) and of every other continuous biomarker in Table S13.
    e_tnfa_haz <- 0.406; label("Log hazard ratio per doubling of baseline TNF-alpha (unitless)") # Table S13 'Baseline TNF-alpha (pg/mL)' 0.406 (SE 0.143); log2 scale assumed
  })

  model({
    lhr <- e_dose_idecel_haz * (DOSE_IDECEL_MCELLS - 450) +
      e_wt_haz * (WT - 75) +
      e_sbcma_haz * log2(SBCMA / 144.5) +
      e_sexm_haz * (1 - SEXF) +
      e_line_4l_haz * LINE_4L +
      e_line_5l_haz * LINE_5L +
      e_form_idecel_b_haz * FORM_IDECEL_B +
      e_ada_pos_haz * ADA_POS +
      e_refractory_triple_class_haz * REFRACTORY_TRIPLE_CLASS +
      e_emd_haz * DIS_EMD +
      e_ferritin_bl_haz * log2(FERRITIN_BL / 163) +
      e_last_prior_corticosteroid_haz * LAST_PRIOR_CORTICOSTEROID +
      e_tnfa_haz * log2(TNFA / 27)

    hr <- exp(lhr)
  })
}
