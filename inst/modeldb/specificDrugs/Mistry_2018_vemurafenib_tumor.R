Mistry_2018_vemurafenib_tumor <- function() {
  description <- paste(
    "Piecewise-linear mixed-effects model of individual target-lesion longest",
    "diameter in BRAF-mutant metastatic melanoma on vemurafenib (BRIM-3 phase",
    "III arm, patients with SD/PR/CR at the first on-treatment scan). Each",
    "lesion shrinks linearly at slope d until day 63 and then grows at slope g;",
    "baseline diameter is log-normal and d is normal, each with nested",
    "between-patient and within-patient (lesion) random effects. No random",
    "effect on g was supported for vemurafenib. Drug exposure is not an input;",
    "the treatment effect is absorbed into the slopes.",
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
    n_subjects = 203L,
    n_studies = 1L,
    age_range = "not reported in this paper (see Chapman 2011 for the BRIM-3 demographics)",
    weight_range = "not reported in this paper",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in this paper",
    disease_state = "BRAF V600 mutation-positive metastatic melanoma (BRIM-3 vemurafenib arm)",
    dose_range = "vemurafenib phase III regimen of BRIM-3 (dose not stated in Mistry 2018; not a model input)",
    regions = "not reported in this paper",
    notes = paste(
      "Data obtained through clinicalstudydatarequest.com. Only patients whose first on-treatment imaging visit showed stable disease, partial or complete response were analysed (Mistry 2018 Methods, Patients).",
      "Mistry 2018 Table 1 (vemurafenib): N = 203; baseline SLD median 72 mm (IQR 39-122); individual longest diameter median 17 mm (IQR 10-29); week-6 ORR 121 (60%); week-6 percent change in SLD median -34% (IQR -47 to -21).",
      "Fitted with nlme in R 3.0.2; final model is Model 16 of Supplementary Table S1 (log-likelihood -9386).",
      sep = "\n"
    )
  )

  ini({
    # Supplementary Table S4, vemurafenib row. log(BSL) mean on the natural-log
    # scale; d and g on the linear scale (mm/day).
    lrbase <- 3.38; label("Mean log baseline lesion longest diameter log(BSL) (log mm)") # Table S4 vemurafenib: log(BSL) mean 3.38 (95% CI 3.3-3.46)
    slope_shrink <- -0.26; label("Mean lesion shrinkage slope d (mm/day)") # Table S4 vemurafenib: d mean -0.26 (95% CI -0.33 to -0.20)
    slope_regrow <- 0.02; label("Mean lesion regrowth slope g (mm/day)") # Table S4 vemurafenib: g mean 0.02 (95% CI 0.003-0.04)
    tregrow <- fixed(63); label("Switch time from shrinkage to regrowth sp (day)") # Table S4 SP = 63 days; chosen from the grid 63/116/179 by log-likelihood (Methods; Table S1 Models 16-18), not estimated

    # Between-patient (BP) random effects: one draw per patient (ID).
    # Table S4 reports standard deviations; variance = SD^2.
    etalrbase ~ 0.1681 # Table S4 vemurafenib log(BSL) s.d. (BP) 0.41; 0.41^2
    etaslope_shrink ~ 0.1444 # Table S4 vemurafenib d s.d. (BP) 0.38; 0.38^2

    # Within-patient (WP) random effects: one draw per lesion of each patient.
    etalrbase_lesion ~ 0.1156 | LESION # Table S4 vemurafenib log(BSL) s.d. (WP) 0.34; 0.34^2
    etaslope_shrink_lesion ~ 0.1024 | LESION # Table S4 vemurafenib d s.d. (WP) 0.32; 0.32^2
    # No random effect on g: Table S4 g s.d. (BP) and (WP) = NA (variance shrank to 0, Table S1 Model 13 vs 16).

    # Residual error: Mistry 2018 states an additive residual e_ijk but does
    # not report its magnitude anywhere in the article or supplement.
    addSd <- fixed(0); label("Additive residual error SD (mm); magnitude not reported in the source")
  })

  model({
    # Individual lesion parameters (Mistry 2018 Methods: BSL log-normal; d and
    # g normal; all random effects independent).
    rbase <- exp(lrbase + etalrbase + etalrbase_lesion)
    slope_shrink_i <- slope_shrink + etaslope_shrink + etaslope_shrink_lesion
    slope_regrow_i <- slope_regrow

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
