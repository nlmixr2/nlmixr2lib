Courlet_2023_melanoma_ici_tgi_radiomics <- function() {
  description <- paste(
    "Tumor growth inhibition (TGI) model of real-world total tumor burden",
    "(TTB) in adults with advanced cutaneous melanoma receiving first-line",
    "immune checkpoint inhibitors, with PET radiomics covariates (Courlet",
    "2023, exploratory substudy model including radiomics features, 38",
    "patients of the PET/CT cohort). Same structure as",
    "Courlet_2023_melanoma_ici_tgi: exponential growth with a fixed",
    "literature kgrowth and a log-kill ON/OFF treatment effect gated by",
    "TRT_PHASE, one kkill for the three drugs, no pharmacokinetics. TTB0",
    "depends on two baseline PET texture features (GLDM dependence entropy,",
    "GLRLM gray level non-uniformity) and kgrowth on one baseline and two",
    "delta (first follow-up minus baseline) features, each as",
    "exp(beta * feature / feature mean). Log-normal IIV on TTB0, kgrowth and",
    "kkill; additive residual error on the log scale.",
    sep = " "
  )
  reference <- paste(
    "Courlet P, Abler D, Guidi M, Girard P, Amato F, Vietti Violi N,",
    "Dietz M, Guignard N, Wicky A, Latifyan S, De Micheli R, Jreige M,",
    "Dromain C, Csajka C, Prior JO, Venkatakrishnan K, Michielin O,",
    "Cuendet MA, Terranova N. Modeling tumor size dynamics based on",
    "real-world electronic health records and image data in advanced",
    "melanoma patients receiving immunotherapy. CPT Pharmacometrics Syst",
    "Pharmacol. 2023;12(8):1170-1181. doi:10.1002/psp4.12983.",
    "Parameter values are from Supplementary Material S5 ('Substudy model",
    "including radiomics features') and its footnote equations; the radiomics",
    "covariate form is from Supplementary Material S1 (Pharmacometric",
    "covariates assessment) and the structural model from the Monolix model",
    "file of Supplementary Material S6.",
    sep = " "
  )
  vignette <- "Courlet_2023_melanoma_ici_tgi"

  units <- list(
    time = "day",
    dosing = "none (no drug input; the ON/OFF treatment effect is switched by the TRT_PHASE covariate)",
    concentration = "cm^3 (the observed variable is total tumor burden, a tumor volume, not a drug concentration)"
  )

  compartmentData <- list(
    tumor_vol = list(
      analyte = "total tumor burden (sum of the volumes of all segmented malignant lesions)",
      units = "cm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  radiomicsNote <- paste(
    "Pyradiomics feature computed on the PET standardized-uptake-value (SUV)",
    "map resampled to 1 x 1 x 1 mm3, fixed bin width 0.2; per-lesion values",
    "are aggregated over all segmented malignant lesions by a mean weighted",
    "by each lesion's share of TTB (Methods, Radiomics features as imaging",
    "covariates on a subcohort)."
  )

  covariateData <- list(
    RADIOMIC_PET_GLDM_DEPENDENCE_ENTROPY_BL = list(
      description = "Baseline PET radiomics feature: GLDM dependence entropy",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        radiomicsNote,
        "Baseline scan. Enters TTB0 as exp(3.32 * x / 6.29); 6.29 is the feature mean in the substudy (Supplementary Material S5 footnote equation)."
      ),
      source_name = "dependence entropy (GLDM, Supplementary Material S4)"
    ),
    RADIOMIC_PET_GLRLM_GRAY_LEVEL_NONUNIFORMITY_BL = list(
      description = "Baseline PET radiomics feature: GLRLM gray level non-uniformity",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        radiomicsNote,
        "Baseline scan. Enters TTB0 as exp(0.66 * x / 74.87); 74.87 is the feature mean in the substudy."
      ),
      source_name = "gray level non-uniformity (GLRLM, Supplementary Material S4)"
    ),
    RADIOMIC_PET_GLDM_SMALL_DEPENDENCE_LOW_GRAY_LEVEL_EMPHASIS_BL = list(
      description = "Baseline PET radiomics feature: GLDM small dependence low gray level emphasis",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        radiomicsNote,
        "Baseline scan. Enters kgrowth as exp(-1.04 * x / 0.013); 0.013 is the feature mean in the substudy."
      ),
      source_name = "small dependence low gray level emphasis (GLDM, Supplementary Material S4)"
    ),
    RADIOMIC_PET_GLDM_DEPENDENCE_ENTROPY_DELTA = list(
      description = "Change in PET radiomics feature GLDM dependence entropy from baseline to the first follow-up scan",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        radiomicsNote,
        "Delta feature = value at the first follow-up minus value at baseline; time-fixed per patient. Enters kgrowth as exp(0.20 * x / 0.11); 0.11 is the feature mean in the substudy."
      ),
      source_name = "delta dependence entropy (GLDM, Supplementary Material S4)"
    ),
    RADIOMIC_PET_SHAPE_ELONGATION_DELTA = list(
      description = "Change in PET radiomics shape feature elongation from baseline to the first follow-up scan",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        radiomicsNote,
        "Delta feature = value at the first follow-up minus value at baseline; time-fixed per patient. Enters kgrowth as exp(0.07 * x / -0.028); -0.028 is the (negative) feature mean in the substudy, so a more negative delta elongation gives a higher kgrowth."
      ),
      source_name = "delta shape elongation (shape, Supplementary Material S4)"
    ),
    TRT_PHASE = list(
      description = "ON/OFF immune checkpoint inhibitor treatment regressor",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-varying. 0 before the start of ICI treatment (baseline images); 1 from ICI initiation onward, INCLUDING after treatment discontinuation (Supplementary Material S6: 'trt (=0 for baseline images and =1 for on-treatment images and images collected after treatment discontinuation)'). Set TRT_PHASE = 1 for post-discontinuation records to reproduce the model. Carried forward from the preceding record, as the Monolix regressor is.",
      source_name = "trt"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38L,
    n_studies = 1L,
    disease_state = "advanced cutaneous melanoma, first-line immune checkpoint inhibitor therapy; PET/CT imaging subcohort with radiomics features at baseline and first follow-up",
    dose_range = "ipilimumab, nivolumab and/or pembrolizumab as in the full cohort; dose was not a model input (ON/OFF effect).",
    regions = "Switzerland (Lausanne University Hospital)",
    notes = "Subset of the 91-patient real-world cohort (Courlet 2023 Table 1) for whom PET radiomics features could be extracted at baseline and first follow-up (Results, Data). The paper does not tabulate demographics or feature distributions for the subcohort beyond the feature means printed in the Supplementary Material S5 equations (Supplementary Material S4 shows the feature distributions graphically). The authors present this model as an exploratory proof of concept with limited parameter precision (TTB0 RSE 117%, omega kgrowth RSE 355%)."
  )

  ini({
    # Supplementary Material S5, 'Substudy model including radiomics
    # features' column.
    lrbase <- log(0.18)
    label("TTB0 at all radiomics features equal to 0 (cm^3)") # Supplementary Table S5 TTB0 = 0.18 cm3 (RSE 117.0%)
    lkgrowth <- fixed(log(0.0081))
    label("Exponential tumor growth rate constant kgrowth at all radiomics features equal to 0 (1/day)") # Supplementary Table S5 kgrowth = 0.0081 day-1, FIXED
    lkkill <- log(0.0044)
    label("Tumor killing rate constant kkill shared by the three ICIs (1/day)") # Supplementary Table S5 kkill = 0.0044 day-1 (RSE 56.8%)

    e_gldm_de_bl_rbase <- 3.32
    label("Effect of baseline GLDM dependence entropy (scaled by its mean 6.29) on log TTB0 (unitless)") # Supplementary Table S5 dependence entropy = 3.32 (RSE 35.3%)
    e_glrlm_gln_bl_rbase <- 0.66
    label("Effect of baseline GLRLM gray level non-uniformity (scaled by its mean 74.87) on log TTB0 (unitless)") # Supplementary Table S5 gray level non-uniformity = 0.66 (RSE 12.0%)
    e_gldm_sdlgle_bl_kgrowth <- -1.04
    label("Effect of baseline GLDM small dependence low gray level emphasis (scaled by its mean 0.013) on log kgrowth (unitless)") # Supplementary Table S5 small dependence low gray level emphasis = -1.04 (RSE 42.8%)
    e_gldm_de_delta_kgrowth <- 0.20
    label("Effect of delta GLDM dependence entropy (scaled by its mean 0.11) on log kgrowth (unitless)") # Supplementary Table S5 delta dependence entropy = 0.20 (RSE 23.3%)
    e_shape_elong_delta_kgrowth <- 0.07
    label("Effect of delta shape elongation (scaled by its mean -0.028) on log kgrowth (unitless)") # Supplementary Table S5 delta shape elongation = 0.07 (RSE 41.6%)

    # IIV: omega reported as a coefficient of variation (%) of a log-normal
    # parameter, omega^2 = log(1 + CV^2), the same convention as Table 2 of
    # the main model. Under this reading the paper's statement that the
    # radiomics features explained 75% and 71% of the variability on TTB0
    # and kgrowth while the variability on kkill increased by 35% (Results,
    # Model evaluation) is reproduced as the change in omega from the
    # substudy base model. Diagonal omega matrix.
    etalrbase ~ log(1 + 0.30^2) # Supplementary Table S5 omega TTB0 = 30%
    etalkgrowth ~ log(1 + 0.11^2) # Supplementary Table S5 omega kgrowth = 11%
    etalkkill ~ log(1 + 5.05^2) # Supplementary Table S5 omega kkill = 505%

    expSd <- 0.78
    label("Additive residual error on the log scale (SD, unitless)") # Supplementary Table S5 additive residual error (log scale) = 0.78 (RSE 8.13%)
  })
  model({
    # Supplementary Material S5 footnote equations:
    # TTB0,i = 0.18 x (e^(dependence entropy,i/6.29))^3.32
    #               x (e^(gray level non-uniformity,i/74.87))^0.66
    # kgrowth,i = 0.0081 x (e^(SDLGLE,i/0.013))^-1.04
    #                    x (e^(delta dependence entropy,i/0.11))^0.20
    #                    x (e^(delta shape elongation,i/-0.028))^0.07
    # Radiomics features enter as P = Ppop x (e^(RAD/RAD_mean))^beta, without
    # log-transformation because they can be negative (Supplementary
    # Material S1, Pharmacometric covariates assessment).
    rbase <- exp(
      lrbase +
        e_gldm_de_bl_rbase * RADIOMIC_PET_GLDM_DEPENDENCE_ENTROPY_BL / 6.29 +
        e_glrlm_gln_bl_rbase * RADIOMIC_PET_GLRLM_GRAY_LEVEL_NONUNIFORMITY_BL / 74.87 +
        etalrbase
    )
    kgrowth <- exp(
      lkgrowth +
        e_gldm_sdlgle_bl_kgrowth * RADIOMIC_PET_GLDM_SMALL_DEPENDENCE_LOW_GRAY_LEVEL_EMPHASIS_BL / 0.013 +
        e_gldm_de_delta_kgrowth * RADIOMIC_PET_GLDM_DEPENDENCE_ENTROPY_DELTA / 0.11 +
        e_shape_elong_delta_kgrowth * RADIOMIC_PET_SHAPE_ELONGATION_DELTA / (-0.028) +
        etalkgrowth
    )
    kkill <- exp(lkkill + etalkkill)

    # Supplementary Material S6: exponential growth; log-kill treatment
    # effect switched on by the trt regressor; the derivative is set to 0
    # once TTB exceeds 1e12 to avoid infinite values.
    tumor_vol(0) <- rbase
    if (tumor_vol > 1e12) {
      dtumor <- 0
    } else {
      dtumor <- kgrowth * tumor_vol - TRT_PHASE * kkill * tumor_vol
    }
    d/dt(tumor_vol) <- dtumor

    tumor_vol ~ lnorm(expSd)
  })
}
