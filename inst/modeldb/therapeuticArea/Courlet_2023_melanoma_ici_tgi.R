Courlet_2023_melanoma_ici_tgi <- function() {
  description <- paste(
    "Tumor growth inhibition (TGI) model of real-world total tumor burden",
    "(TTB, summed volume of all segmented malignant lesions on PET/CT or CT)",
    "in adults with advanced cutaneous melanoma receiving first-line immune",
    "checkpoint inhibitors (ipilimumab, nivolumab, pembrolizumab, alone or in",
    "combination; Courlet 2023, final model including clinical covariates).",
    "Exponential growth with a fixed literature growth rate kgrowth and a",
    "log-kill ON/OFF treatment effect: d/dt(TTB) = kgrowth * TTB -",
    "TRT_PHASE * kkill * TTB, with a single kkill shared by the three drugs and",
    "no pharmacokinetics. Baseline TTB depends on serum albumin and the",
    "neutrophil-to-lymphocyte ratio (power functions) and on ECOG > 0",
    "(exponential); kgrowth is lower for NRAS-mutated tumors. Log-normal IIV",
    "on TTB0, kgrowth and kkill; additive residual error on the log scale.",
    sep = " "
  )
  reference <- paste(
    "Courlet P, Abler D, Guidi M, Girard P, Amato F, Vietti Violi N,",
    "Dietz M, Guignard N, Wicky A, Latifyan S, De Micheli R, Jreige M,",
    "Dromain C, Csajka C, Prior JO, Venkatakrishnan K, Michielin O,",
    "Cuendet MA, Terranova N. Modeling tumor size dynamics based on real-world electronic",
    "health records and image data in advanced melanoma patients receiving",
    "immunotherapy. CPT Pharmacometrics Syst Pharmacol. 2023;12(8):1170-1181.",
    "doi:10.1002/psp4.12983.",
    "Parameter values are from Table 2 ('Model including clinical",
    "covariates') and its footnote equations; the structural model is the",
    "Monolix model file of Supplementary Material S6 and the covariate",
    "definitions are the Monolix project of Supplementary Material S7.",
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

  covariateData <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect on TTB0 normalized to 42 g/L as printed in the Table 2 footnote equation; the Monolix project (Supplementary Material S7) normalizes to the weighted mean 42.2266 g/L. Albumin was missing for 27 of 91 patients (30%) and imputed to the median (Methods, Clinical covariates).",
      source_name = "ALB"
    ),
    NLR = list(
      description = "Baseline neutrophil-to-lymphocyte ratio",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Power effect on TTB0 normalized to 3.14 as printed in the Table 2 footnote equation (Monolix project: 3.14052). Missing values were imputed to the median.",
      source_name = "NLR"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status greater than 0",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = ECOG 1, 2 or 4 at baseline; 0 = ECOG 0. The source column ECOG was grouped as ECOG0 = {0, NA} and ECOGGT0 = {1, 2, 4} (Supplementary Material S7), so patients with missing ECOG (17 of 91) take the reference value 0. Exponential effect on TTB0.",
      source_name = "ECOG"
    ),
    TUM_NRAS_MUT = list(
      description = "Tumor NRAS mutation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = any NRAS mutation (Q61R/K/L, G12C/D, and mutations at an unknown position, which Supplementary Material S1 lumps with the codon 12/61 group); 0 = NRAS wild type. The source column NRASMUT was grouped as notmutated_NA = {0, NA} and mutated = 1 (Supplementary Material S7), so the one patient with unknown NRAS status takes the reference value 0. Exponential effect on kgrowth.",
      source_name = "NRASMUT"
    ),
    TRT_PHASE = list(
      description = "ON/OFF immune checkpoint inhibitor treatment regressor",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-varying. 0 before the start of ICI treatment (baseline images); 1 from ICI initiation onward, INCLUDING after treatment discontinuation. Supplementary Material S6: 'trt (=0 for baseline images and =1 for on-treatment images and images collected after treatment discontinuation)'. The paper found that the same parameters describe tumor dynamics during treatment and after discontinuation (Results, Base model), so killing continues after the last dose; set TRT_PHASE = 1 for post-discontinuation records to reproduce the model. Like the Monolix regressor, the value is carried forward from the preceding record.",
      source_name = "trt"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 91L,
    n_studies = 1L,
    age_range = "median 68 years (IQR 60-76)",
    sex_female_pct = 35,
    disease_state = "advanced cutaneous melanoma, first-line immune checkpoint inhibitor therapy",
    dose_range = "ipilimumab + nivolumab (n = 20), ipilimumab + nivolumab followed by nivolumab maintenance (n = 18), ipilimumab alone (n = 23), nivolumab alone (n = 3), pembrolizumab (n = 27); median doses 225 mg ipilimumab, 85 mg nivolumab, 200 mg pembrolizumab. Dose was not a model input (ON/OFF effect).",
    regions = "Switzerland (Lausanne University Hospital)",
    notes = "Retrospective real-world cohort from electronic health records and routine imaging (Courlet 2023 Table 1). 311 TTB measurements (240 PET/CT, 71 CT), median 3 per patient (range 2-10); baseline TTB missing in 27 patients. Albumin median 42.5 g/L (IQR 40.0-45.0); NLR median 3.4 (IQR 2.4-4.3); ECOG 0/1/2/4/missing = 51/18/4/1/17; NRAS mutated 34 (37%, including 3 at an unknown position), wild type 56, unknown 1. BMI median 26 kg/m^2."
  )

  ini({
    # Table 2, 'Model including clinical covariates' column.
    lrbase <- log(7.78)
    label("Typical total tumor burden at baseline TTB0 (cm^3)") # Table 2 TTB0 = 7.78 cm3 (RSE 15.6%)
    lkgrowth <- fixed(log(0.0081))
    label("Exponential tumor growth rate constant kgrowth (1/day)") # Table 2 kgrowth = 0.0081 day-1, Fixed (Chatterjee 2017 linear-SLD rate x 3, Supplementary Material S1)
    lkkill <- log(0.008)
    label("Tumor killing rate constant kkill shared by the three ICIs (1/day)") # Table 2 kkill = 0.008 day-1 (RSE 22.1%)

    e_alb_rbase <- -5.27
    label("Power exponent of albumin/42 on TTB0 (unitless)") # Table 2 beta albumin = -5.27 (RSE 27.1%)
    e_nlr_rbase <- 0.76
    label("Power exponent of NLR/3.14 on TTB0 (unitless)") # Table 2 beta neutrophil-lymphocyte ratio = 0.76 (RSE 36.5%)
    e_ecog_ge1_rbase <- 0.65
    label("Exponential effect of ECOG > 0 on TTB0 (unitless)") # Table 2 beta ECOG >0 = 0.65 (RSE 46.1%)
    e_tum_nras_mut_kgrowth <- -1.22
    label("Exponential effect of NRAS mutation on kgrowth (unitless)") # Table 2 beta NRAS mutated = -1.22 (RSE 44.9%)

    # IIV: Table 2 reports omega as a coefficient of variation (%) of a
    # log-normally distributed parameter, CV = sqrt(exp(omega^2) - 1), so
    # omega^2 = log(1 + CV^2). Under this reading the paper's statement that
    # the covariates explained 24% and 16% of the IIV on TTB0 and kgrowth
    # (Results, Model evaluation) is reproduced as the reduction in omega
    # from the base model. Diagonal omega matrix.
    etalrbase ~ log(1 + 1.23^2) # Table 2 omega TTB0 = 123%
    etalkgrowth ~ log(1 + 0.49^2) # Table 2 omega kgrowth = 49%
    etalkkill ~ log(1 + 1.85^2) # Table 2 omega kkill = 185%

    # Monolix errorModel = constant(a) on a logNormal observation
    # (Supplementary Material S7): log(TTB_obs) = log(TTB) + a * eps.
    expSd <- 1.13
    label("Additive residual error on the log scale (SD, unitless)") # Table 2 additive residual error (log scale) = 1.13 (RSE 5.41%)
  })
  model({
    # Table 2 footnote equations:
    # TTB0,i = 7.78 x (ALB,i/42)^-5.27 x (NLR,i/3.14)^0.76 x e^(0.65 x ECOGcat,i)
    # kgrowth,i = 0.0081 x e^(-1.22 x NRAScat,i)
    rbase <- exp(
      lrbase +
        e_alb_rbase * log(ALB / 42) +
        e_nlr_rbase * log(NLR / 3.14) +
        e_ecog_ge1_rbase * ECOG_GE1 +
        etalrbase
    )
    kgrowth <- exp(lkgrowth + e_tum_nras_mut_kgrowth * TUM_NRAS_MUT + etalkgrowth)
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
