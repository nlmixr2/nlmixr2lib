Wu_2025_idecabtageneVicleucel_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for progression-free ",
    "survival (PFS, independent review committee, IMWG criteria) after a ",
    "single infusion of the BCMA-directed CAR T-cell therapy idecabtagene ",
    "vicleucel (ide-cel) in adults with triple-class-exposed ",
    "relapsed/refractory multiple myeloma (Wu 2025; ide-cel arm of the ",
    "phase 3 KarMMa-3 trial, N = 220). The exposure is the noncompartmental ",
    "day-0-28 transgene-level AUC (AUC_IDECEL, day*copies/ug, droplet ",
    "digital PCR) on the log2 scale: each doubling multiplies the hazard by ",
    "exp(-0.43) = 0.65. These are the Table S7 FULL-model coefficients ",
    "(13 covariate terms), the only PFS estimates the paper prints; the ",
    "Results state that backward elimination then removed last prior ",
    "anti-CD38 therapy and baseline soluble CD137, whose refitted final ",
    "coefficients are not reported. The model returns the RELATIVE hazard ",
    "(hr, and its log lhr) against the reference subject (female, ",
    "Formulation A, two prior regimens, no extramedullary plasmacytoma, ",
    "ADA-negative, not triple-class refractory, no steroid or anti-CD38 in ",
    "the last prior regimen, continuous covariates at their medians). No ",
    "baseline hazard is encoded: a Cox regression leaves h0(t) unspecified. ",
    "No ODE, no random effect, no residual error."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S7.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (semiparametric Cox relative hazard; the baseline hazard and therefore the time scale are unspecified by the fit)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IDECEL covariate column)",
    concentration = "hr (relative hazard of progression or death vs the reference subject, unitless; also lhr, the log relative hazard)"
  )

  covariateData <- list(
    AUC_IDECEL = list(
      description = "Individual ide-cel transgene-level AUC from the infusion to day 28 (AUC0-28D), by noncompartmental analysis of peripheral-blood CAR transgene copies per microgram of genomic DNA measured by droplet digital PCR.",
      units = "day*copies/ug",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as log2(AUC_IDECEL / 1.3e6). The paper centres continuous covariates on their medians without printing the AUC0-28D median; 1.3 x 10^6 is the boundary between the second and third quartiles printed in the Figure 2a legend (Q2 0.4-1.3, Q3 1.3-2.7 x 10^6). The centring only sets which subject has hr = 1; it does not change any hazard ratio. ddPCR assay (KarMMa-3); the paper did not compare these values with the KarMMa (phase 2) assay.",
      source_name = "Log2(AUC0-28D) (days*copies/ug)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "The paper codes male vs female with female as reference, so the model applies the coefficient to (1 - SEXF).",
      source_name = "Sex (Male vs Female)"
    ),
    DIS_EMD = list(
      description = "1 = extramedullary plasmacytoma present at baseline, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-selected covariate (significant in the KarMMa analysis).",
      source_name = "Baseline Extramedullary Plasmacytoma (Yes vs No)"
    ),
    FERRITIN_BL = list(
      description = "Baseline serum ferritin.",
      units = "ug/L",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, centred on the analysis-set median 163 ug/L (Table S5 note).",
      source_name = "Log2(Baseline Ferritin) (ug/L)"
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
      notes = "The Schoenfeld test (Figure S4) indicated that formulation violates the proportional-hazards assumption, and a sensitivity analysis truncated at 12 months showed no significant formulation effect (Figure S17).",
      source_name = "Product Attributes (Formulation B vs A)"
    ),
    LAST_PRIOR_CORTICOSTEROID = list(
      description = "1 = the last prior multiple-myeloma regimen before ide-cel included a corticosteroid, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "93.6% of the analysis set (Table S3); the authors caution that the effect is estimated from very unbalanced groups.",
      source_name = "Last Prior Med: Steroids (Yes vs No)"
    ),
    ADA_POS = list(
      description = "1 = anti-CAR antidrug antibody positive at at least one post-infusion visit, 0 = negative.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "POST-BASELINE covariate (61.4% positive). The authors caution that patients with longer PFS had more follow-up in which to become ADA-positive, so the protective coefficient is confounded.",
      source_name = "Post-infusion ADA Positive (Yes vs No)"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen, measured on the day of the ide-cel infusion.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, centred on the analysis-set median 144.5 ng/mL (Table S5 note).",
      source_name = "Log2(Baseline Soluble BCMA) (ng/mL)"
    ),
    REFRACTORY_TRIPLE_CLASS = list(
      description = "1 = refractory to all three classes (an immunomodulatory agent, a proteasome inhibitor and an anti-CD38 antibody), 0 = not triple-class refractory.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "61.8% of the analysis set (Table S3). Every KarMMa-3 patient was triple-class EXPOSED; this flag marks refractoriness.",
      source_name = "Prior Triple Refractory (Yes vs No)"
    ),
    TNFA = list(
      description = "Baseline serum tumour necrosis factor alpha, measured on the day of the ide-cel infusion.",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, centred on the median 27 pg/mL (Table S3).",
      source_name = "Log2(Baseline TNF-alpha) (pg/mL)"
    ),
    LAST_PRIOR_ANTICD38 = list(
      description = "1 = the last prior multiple-myeloma regimen before ide-cel included an anti-CD38 monoclonal antibody, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "72.7% of the analysis set (Table S3). In the Table S7 full model; removed from the final model by backward elimination (Results).",
      source_name = "Last Prior Med: Anti-CD38 Antibodies (Yes vs No)"
    ),
    SCD137 = list(
      description = "Baseline serum soluble CD137 (4-1BB), measured on the day of the ide-cel infusion.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "log2, centred on the median 0.105 ng/mL (Table S3). In the Table S7 full model; removed from the final model by backward elimination (Results).",
      source_name = "Log2(Baseline sCD137) (ng/mL)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_observations = "220 right-censored progression-free-survival times (one per patient)",
    age_range = "30-81 years (median 63; 43% aged 65 or older)",
    weight_range = "41.1-144 kg (median 83)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens with progression within 60 days of the last regimen; 61.8% triple-class refractory, 23.6% extramedullary plasmacytoma",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted)",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, CK data cutoff 1 March 2022; 220 subjects with evaluable exposure parameters. Kaplan-Meier PFS by AUC0-28D quartile is plotted in Figure 2a."
  )

  ini({
    # Wu 2025 Table S7, 'Parameter Estimates of PFS E-R Full Model in
    # MM-003' (log hazard ratios). Log2-transformed continuous covariates,
    # centred on the medians (Table S7 note; see covariateData for the
    # reference values). The Results state that last prior anti-CD38
    # therapy and baseline sCD137 were then removed by backward
    # elimination, but the refitted final model is not printed, so the
    # full model is carried as published.
    e_auc_idecel_haz <- -0.43; label("Log hazard ratio per doubling of AUC0-28D (unitless)") # Table S7 'Log2(AUC0-28D)' -0.43 (SE 0.056; 95% CI -0.54 to -0.32)
    e_sexm_haz <- 0.25; label("Log hazard ratio for male vs female sex (unitless)") # Table S7 'Sex (Male vs Female)' 0.25 (SE 0.203)
    e_emd_haz <- 0.17; label("Log hazard ratio for baseline extramedullary plasmacytoma (unitless)") # Table S7 'Baseline Extramedullary Plasmacytoma (Yes vs No)' 0.17 (SE 0.213)
    e_ferritin_bl_haz <- 0.1; label("Log hazard ratio per doubling of baseline ferritin (unitless)") # Table S7 'Log2(Baseline Ferritin)' 0.1 (SE 0.06)
    e_line_4l_haz <- -0.06; label("Log hazard ratio for 3 vs 2 prior MM regimens (unitless)") # Table S7 'No. of prior MM Regimens (3 vs 2)' -0.06 (SE 0.235)
    e_line_5l_haz <- -0.38; label("Log hazard ratio for 4 vs 2 prior MM regimens (unitless)") # Table S7 'No. of prior MM Regimens (4 vs 2)' -0.38 (SE 0.275)
    e_form_idecel_b_haz <- -0.67; label("Log hazard ratio for Formulation B vs A (unitless)") # Table S7 'Product Attributes (Formulation B vs A)' -0.67 (SE 0.278)
    e_last_prior_corticosteroid_haz <- 1.16; label("Log hazard ratio for a corticosteroid in the last prior regimen (unitless)") # Table S7 'Last Prior Med: Steroids (Yes vs No)' 1.16 (SE 0.474)
    e_ada_pos_haz <- -0.89; label("Log hazard ratio for post-infusion ADA positivity (unitless)") # Table S7 'Post-infusion ADA Positive (Yes vs No)' -0.89 (SE 0.22)
    e_sbcma_haz <- 0.28; label("Log hazard ratio per doubling of baseline soluble BCMA (unitless)") # Table S7 'Log2(Baseline Soluble BCMA)' 0.28 (SE 0.081)
    e_refractory_triple_class_haz <- 0.61; label("Log hazard ratio for triple-class refractory disease (unitless)") # Table S7 'Prior Triple Refractory (Yes vs No)' 0.61 (SE 0.216)
    e_tnfa_haz <- 0.33; label("Log hazard ratio per doubling of baseline TNF-alpha (unitless)") # Table S7 'Log2(Baseline TNF-alpha)' 0.33 (SE 0.197)
    e_last_prior_anticd38_haz <- -0.32; label("Log hazard ratio for an anti-CD38 antibody in the last prior regimen (unitless)") # Table S7 'Last Prior Med: Anti-CD38 Antibodies (Yes vs No)' -0.32 (SE 0.229)
    e_scd137_haz <- 0.1; label("Log hazard ratio per doubling of baseline soluble CD137 (unitless)") # Table S7 'Log2(Baseline sCD137)' 0.1 (SE 0.159)
  })

  model({
    lhr <- e_auc_idecel_haz * log2(AUC_IDECEL / 1.3e6) +
      e_sexm_haz * (1 - SEXF) +
      e_emd_haz * DIS_EMD +
      e_ferritin_bl_haz * log2(FERRITIN_BL / 163) +
      e_line_4l_haz * LINE_4L +
      e_line_5l_haz * LINE_5L +
      e_form_idecel_b_haz * FORM_IDECEL_B +
      e_last_prior_corticosteroid_haz * LAST_PRIOR_CORTICOSTEROID +
      e_ada_pos_haz * ADA_POS +
      e_sbcma_haz * log2(SBCMA / 144.5) +
      e_refractory_triple_class_haz * REFRACTORY_TRIPLE_CLASS +
      e_tnfa_haz * log2(TNFA / 27) +
      e_last_prior_anticd38_haz * LAST_PRIOR_ANTICD38 +
      e_scd137_haz * log2(SCD137 / 0.105)

    hr <- exp(lhr)
  })
}
