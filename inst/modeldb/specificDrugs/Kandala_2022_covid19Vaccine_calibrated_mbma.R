Kandala_2022_covid19Vaccine_calibrated_mbma <- function() {
  description <- paste0(
    "MBMA (clinical efficacy, calibrated from preclinical). COVID-19 vaccine ",
    "immune correlate of protection: the calibrated rhesus-macaque (RM) MBMA, ",
    "the first of three efficacy-scale models in Kandala 2022. The RM ",
    "viral-load suppression curve (BAL/lung parameters) is translated to ",
    "clinical vaccine efficacy (VE %) by (1) a generalized-least-squares ",
    "calibration factor Fc that rescales the convalescent-normalized clinical ",
    "SN/HCS ratio onto the RM SN-titre axis, and (2) the linear VL-to-",
    "incidence-rate mapping with IRmin fixed to 0 and IRmax to 1. VE is then ",
    "1 - IR. Predicts clinical efficacy from clinical immunogenicity (SN/HCS) ",
    "using only the RM model plus the calibration factor. The 50% protective ",
    "titre is predicted at 24.3% of the geometric-mean human convalescent ",
    "titre (Figure 3). Deterministic typical-value translation (no random ",
    "effects; parameter uncertainty is propagated by sampling in the paper), ",
    "so it simulates typical vaccine efficacy and is NOT an individual-level ",
    "model. The underlying RM model is ",
    "modellib('Kandala_2022_covid19Vaccine_rm_mbma')."
  )

  reference <- paste(
    "Kandala B, Plock N, Chawla A, Largajolli A, Robey S, Watson K, et al.",
    "Accelerating model-informed decisions for COVID-19 vaccine candidates",
    "using a model-based meta-analysis approach.",
    "EBioMedicine. 2022;84:104264. doi:10.1016/j.ebiom.2022.104264.",
    sep = " "
  )

  vignette <- "Kandala_2022_covid19Vaccine"

  units <- list(
    time = paste0(
      "not applicable -- static study-level meta-analysis with no time ",
      "course. Immunogenicity is the SN titre closest to the start of the ",
      "efficacy-assessment period; VE is accrued over the trial."
    ),
    dosing = paste0(
      "not applicable -- vaccine exposure enters through the SN_HCS_RATIO ",
      "covariate, not through dosing records."
    ),
    concentration = paste0(
      "ve is the predicted clinical vaccine efficacy in percent. The output ",
      "Cc is unused."
    )
  )

  covariateData <- list(
    SN_HCS_RATIO = list(
      description = paste0(
        "Study-arm ratio of the vaccine SN titre to the same-assay ",
        "geometric-mean human convalescent serum (HCS) titre. Dimensionless ",
        "by construction; harmonizes neutralizing titres across assays."
      ),
      units = "(fold; dimensionless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Kandala 2022 Table 3 columns 'SN Titre (GM)' and 'HCS (GM)'; the ",
        "driver of the calibration (Equation 7). Scaled by the calibration ",
        "factor Fc and converted to the RM log2 SN axis before entering the ",
        "RM suppression curve: log2(SN_cal) = log2(SN/HCS) + log2(Fc), i.e. ",
        "SN_cal = SN_HCS_RATIO * Fc (Supplement Fig S6 worked example). ",
        "Observed clinical values span roughly 0.045 to 3.9."
      ),
      source_name = "SN/HCS"
    )
  )

  population <- list(
    species = "human",
    n_studies = 6L,
    age_range = paste0(
      "Adults; two candidates (Ad26.COV2.S, BNT162b2) contribute separate ",
      "younger-adult and older-adult (65-85 years) efficacy strata."
    ),
    disease_state = paste0(
      "PCR-confirmed symptomatic COVID-19 of any severity; efficacy is the ",
      "relative reduction in incidence rate versus placebo."
    ),
    dose_range = paste0(
      "Six vaccine candidates used in calibration (BNT162b2, mRNA-1273, ",
      "Ad26.COV2.S, 2019nCoV, Gam-COVID-Vac, ChAdOx1 nCoV-19); Phase 1/2 ",
      "immunogenicity matched to Phase 3 efficacy from identical regimens ",
      "(Kandala 2022 Table 3)."
    ),
    regions = "International (published Phase 1/2 and Phase 3 trials)",
    notes = paste0(
      "Model-based meta-analysis: the unit of observation is a trial-level ",
      "efficacy with its matched SN/HCS ratio. The calibration factor Fc was ",
      "estimated by generalized non-linear least squares over n = 8 data ",
      "points (six vaccines; Ad26.COV2.S and BNT162b2 each contribute two ",
      "age strata) with a root-mean-squared prediction error of 7%. ",
      "External validation (BBV152, CoronaVac) and delta-variant predictions ",
      "are in the vignette."
    )
  )

  # ==========================================================================
  # ini(): RM MBMA structural parameters (Kandala 2022 Table 2, BAL/lung set --
  # the paper uses the BAL/lung parameters for all clinical translation because
  # the clinical efficacy endpoint is symptomatic disease of any severity)
  # plus the calibration factor Fc (Figure 3 / Methods). No random effects:
  # the translation is deterministic on the maximum-likelihood estimates, and
  # the paper propagates uncertainty by sampling the variance-covariance
  # matrix rather than through model-level variability.
  # ==========================================================================
  ini({
    vlmax <- 4.86
    label("Maximum log10 post-challenge peak viral load from the RM MBMA (VLmax, log10 copies)")  # Kandala 2022 Table 2, row 'VLmax': 4.86

    hill <- 4.14
    label("Hill coefficient of the RM SN-titre viral-load suppression curve (unitless)")  # Kandala 2022 Table 2, row 'g': 4.14

    it50 <- 5.84
    label("Log2 SN titre giving 50% viral-load suppression in BAL fluid / lung tissue (IT50, log2 titre)")  # Kandala 2022 Table 2, row 'IT50': 5.84 (BAL/lung used for clinical translation)

    vlmin <- fixed(1.6)
    label("Minimum log10 post-challenge peak viral load (VLmin, log10 copies); not a directly published value (back-solved)")  # not from source as a single number -- back-solved to 1.6 from Figure 3 + Table 2 + Fc; see vignette Assumptions and deviations

    fc <- 143.9
    label("Calibration factor mapping the clinical SN/HCS ratio onto the RM SN-titre axis (unitless)")  # Kandala 2022 Figure 3 / Results: 143.9, 95% CI [90.9, 196.9]

    # Residual is a placeholder: the paper's calibration residual scale is not
    # tabulated (uncertainty is reported as prediction intervals from sampling,
    # not a residual SD). The fixed tiny additive term exists only so rxode2
    # has an error model for the typical-value VE output.
    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on VE (percent); NOT a published value")  # not from source -- calibration residual scale is unreported; see vignette Assumptions and deviations
  })

  model({
    # ======================================================================
    # Calibrate the clinical SN/HCS ratio onto the RM log2 SN-titre axis
    # (Kandala 2022 Equation 7; Supplement Fig S6):
    #   log2(SN_cal) = log2(SN/HCS) + log2(Fc) = log2(SN_HCS_RATIO * Fc)
    # ln2 is spelled out because rxode2 provides log() but not log2().
    # ======================================================================
    ln2 <- 0.693147180559945
    sn_cal_log2 <- log(SN_HCS_RATIO * fc) / ln2

    # ======================================================================
    # RM viral-load suppression curve at the calibrated titre (Kandala 2022
    # Equations 3 and 1, BAL/lung IT50).
    # ======================================================================
    suppr <- sn_cal_log2^hill / (sn_cal_log2^hill + it50^hill)
    vl <- vlmax * exp(log(vlmin / vlmax) * suppr)

    # ======================================================================
    # Linear VL-to-incidence-rate mapping with IRmin = 0, IRmax = 1 (Kandala
    # 2022 Equation 5), then VE = 100 * (1 - IR/IR_placebo). The placebo arm
    # has baseline SN near the LLOD so IR_placebo approaches IRmax = 1
    # (Equation 6).
    # ======================================================================
    ir <- (vl - vlmin) / (vlmax - vlmin)
    ve <- 100 * (1 - ir)

    ve ~ add(addSd)
  })
}
