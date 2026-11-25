Mistry_2018_trametinib_tumor <- function() {
  description <- paste(
    "Piecewise-linear mixed-effects model of individual target-lesion longest",
    "diameter in BRAF-mutant metastatic melanoma on trametinib (METRIC phase",
    "III arm, patients with SD/PR/CR at the first on-treatment scan). Each",
    "lesion shrinks linearly at slope d until day 63 and then grows at slope g;",
    "baseline diameter is log-normal and d and g are normal, each with nested",
    "between-patient and within-patient (lesion) random effects. Drug",
    "exposure is not an input; the treatment effect is absorbed into the",
    "slopes.",
    sep = " "
  )
  reference <- paste(
    "Mistry HB, Orrell D, Eftimie R.",
    "Model based analysis of the heterogeneity in the tumour size dynamics",
    "differentiates vemurafenib, dabrafenib and trametinib in metastatic",
    "melanoma. Cancer Chemother Pharmacol. 2018;81(2):325.",
    "doi:10.1007/s00280-017-3486-3. PMID: 29222604. PMCID: PMC5778167.",
    "Parameter values from Electronic Supplementary Material 1, Table S4.",
    sep = " "
  )
  vignette <- "Mistry_2018_melanoma_lesion_dynamics"

  units <- list(
    time = "day",
    dosing = "n/a (no drug input; treatment effect is absorbed into the empirical shrinkage and regrowth slopes)",
    concentration = "mm (the observable `tumor_size` is the longest diameter of one RECIST target lesion)"
  )

  covariateData <- list(
    LESION = list(
      description = "Target-lesion identifier within a patient (1, 2, ...)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      source_name = "j (lesion subscript)",
      notes = paste(
        "Mistry 2018 indexes each RECIST target lesion j within patient i (L_ijk).",
        "Consumed by the ini() nesting syntax '~ ... | LESION', which gives every lesion of a patient its own within-patient (WP) deviation on top of the between-patient (BP) eta shared by all lesions of that patient.",
        "Lesions are measured at the same imaging visits, so records of different lesions interleave in time within one ID.",
        "Not referenced in model(); rxode2 requires an event table to contain at least two distinct LESION values (also when the random effects are zeroed).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 157L,
    n_studies = 1L,
    age_range = "not reported in this paper (see Flaherty 2012 for the METRIC demographics)",
    weight_range = "not reported in this paper",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in this paper",
    disease_state = "BRAF V600 mutation-positive metastatic melanoma (METRIC trametinib arm)",
    dose_range = "trametinib phase III regimen of METRIC (dose not stated in Mistry 2018; not a model input)",
    regions = "not reported in this paper",
    notes = paste(
      "Data obtained through clinicalstudydatarequest.com. Only patients whose first on-treatment imaging visit showed stable disease, partial or complete response were analysed (Mistry 2018 Methods, Patients).",
      "Mistry 2018 Table 1 (trametinib): N = 157; baseline SLD median 64 mm (IQR 32-106); individual longest diameter median 18 mm (IQR 12-30); week-6 ORR 46 (29%); week-6 percent change in SLD median -18% (IQR -31 to -4).",
      "Fitted with nlme in R 3.0.2; final model is Model 8 of Supplementary Table S3 (log-likelihood -5899).",
      sep = "\n"
    )
  )

  ini({
    # Supplementary Table S4, trametinib row. log(BSL) mean on the natural-log
    # scale; d and g on the linear scale (mm/day).
    lrbase <- 3.22; label("Mean log baseline lesion longest diameter log(BSL) (log mm)") # Table S4 trametinib: log(BSL) mean 3.22 (95% CI 3.16-3.29)
    slope_shrink <- -0.13; label("Mean lesion shrinkage slope d (mm/day)") # Table S4 trametinib: d mean -0.13 (95% CI -0.15 to -0.11)
    slope_regrow <- 0.05; label("Mean lesion regrowth slope g (mm/day)") # Table S4 trametinib: g mean 0.05 (95% CI 0.03-0.07)
    tregrow <- fixed(63); label("Switch time from shrinkage to regrowth sp (day)") # Table S4 SP = 63 days; chosen from the grid 63/116/179 by log-likelihood (Methods; Table S3 log-likelihoods at SP 63/116/179), not estimated

    # Between-patient (BP) random effects: one draw per patient (ID).
    # Table S4 reports standard deviations; variance = SD^2.
    etalrbase ~ 0.0676 # Table S4 trametinib log(BSL) s.d. (BP) 0.26; 0.26^2
    etaslope_shrink ~ 0.0081 # Table S4 trametinib d s.d. (BP) 0.09; 0.09^2
    etaslope_regrow ~ 0.0144 # Table S4 trametinib g s.d. (BP) 0.12; 0.12^2

    # Within-patient (WP) random effects: one draw per lesion of each patient.
    etalrbase_lesion ~ 0.1764 | LESION # Table S4 trametinib log(BSL) s.d. (WP) 0.42; 0.42^2
    etaslope_shrink_lesion ~ 0.0169 | LESION # Table S4 trametinib d s.d. (WP) 0.13; 0.13^2
    etaslope_regrow_lesion ~ 0.0009 | LESION # Table S4 trametinib g s.d. (WP) 0.03; 0.03^2

    # Residual error: Mistry 2018 states an additive residual e_ijk but does
    # not report its magnitude anywhere in the article or supplement.
    addSd <- fixed(0); label("Additive residual error SD (mm); magnitude not reported in the source")
  })

  model({
    # Individual lesion parameters (Mistry 2018 Methods: BSL log-normal; d and
    # g normal; all random effects independent).
    rbase <- exp(lrbase + etalrbase + etalrbase_lesion)
    slope_shrink_i <- slope_shrink + etaslope_shrink + etaslope_shrink_lesion
    slope_regrow_i <- slope_regrow + etaslope_regrow + etaslope_regrow_lesion

    # Mistry 2018 piecewise-linear lesion model, as printed (L_ijk equations):
    #   L = BSL + d * t                 if t <  sp
    #   L = BSL + d * sp + g * t        if t >= sp
    # The regrowth arm is g * t (not g * (t - sp)), so the printed model has a
    # step of g * sp at the switch time.
    tumor_size <- rbase + slope_shrink_i * t * (t < tregrow) +
      (slope_shrink_i * tregrow + slope_regrow_i * t) * (t >= tregrow)
    tumor_size ~ add(addSd)
  })
}
