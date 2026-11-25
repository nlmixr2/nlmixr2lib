Zhu_2019_ipatasertib_rpfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-EFFICACY model for radiographic ",
    "progression-free survival (rPFS) as a function of the ipatasertib ",
    "steady-state daily AUC based on the average ACTUAL dose up to the ",
    "event (AUCss,event; primary metric) in men with metastatic ",
    "castration-resistant prostate cancer receiving ipatasertib 200 or ",
    "400 mg once daily or placebo with abiraterone and prednisone ",
    "(A.MARTIN phase II, n = 253). The model returns the hazard ratio ",
    "versus placebo, exp(-0.0000739 * AUC_IPAT_EVT); the exposure ",
    "effect is NOT statistically significant (P = 0.232). NO BASELINE ",
    "HAZARD IS ENCODED: a Cox regression leaves h0(t) unspecified, so ",
    "only the relative hazard is available. There is no PK layer and ",
    "no ODE; exposure is supplied as a data column. No random effect ",
    "and no residual error are estimated. Companion fit on the ",
    "nominal-dose exposure metric: Zhu_2019_ipatasertib_rpfs_nominal."
  )
  reference <- paste(
    "Zhu R, Poland B, Wada R, Liu Q, Musib L, Maslyar D, Cho E, Yu W, Ma H,",
    "Jin JY, Budha N.",
    "Exposure-response-based product profile-driven clinical utility index",
    "for ipatasertib dose selection in prostate cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(4):240.",
    "doi:10.1002/psp4.12394.",
    "Coefficients are transcribed from Supporting Information Table S2",
    "(Table S1 for the Cox proportional-hazards models); exposure units are",
    "taken from the x-axis of Supporting Figure S3.",
    sep = " "
  )
  vignette <- "Zhu_2019_ipatasertib_exposure_response"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IPAT_EVT covariate column)",
    concentration = "hr (hazard ratio of radiographic progression versus placebo, unitless; also lhr, the log hazard ratio)"
  )

  covariateData <- list(
    AUC_IPAT_EVT = list(
      description = "Individual ipatasertib steady-state daily AUC based on the average ACTUAL dose up to and including the time of the event (AUCss,event): average daily dose administered up to the event (or over the whole treatment for a patient without the event) divided by the individual apparent clearance from the ipatasertib population PK model. Supplied as data: this model has no PK layer.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhu 2019 Methods, Response and exposure data. Dose interruptions and",
        "reductions are folded in, so the averaging window is endpoint- and",
        "patient-specific. Enters LINEARLY and uncentred on the linear",
        "predictor; placebo patients carry 0. Supporting Figure S3 plots",
        "the observed range out to roughly 6600 ng*h/mL. Projecting a",
        "nominal dose onto this metric needs the endpoint-specific",
        "dose-intensity model of Zhu 2019 Eqs. 1-2, reproduced in the",
        "vignette.",
        "In this model the exposure enters with slope -0.0000739 per ng*h/mL."
      ),
      source_name = "AUCss,event"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 253L,
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "Metastatic castration-resistant prostate cancer previously treated with docetaxel chemotherapy; randomised, double-blind phase II part of the phase Ib/II A.MARTIN study (NCT01485861).",
    dose_range = "Ipatasertib 200 mg or 400 mg orally once daily, or placebo, each combined with abiraterone 1000 mg and prednisone 10 mg daily.",
    arms = c(placebo = 82L, ipatasertib_200mg = 87L, ipatasertib_400mg = 84L),
    notes = paste(
      "All 253 randomised patients provided evaluable PK, rPFS and safety",
      "data (Zhu 2019 Results). Placebo patients enter the regressions with",
      "zero ipatasertib exposure. Dose modifications, mainly for adverse",
      "events, occurred in 14.6 percent (12/82) of placebo, 14.9 percent",
      "(13/87) of 200 mg and 39.3 percent (33/84) of 400 mg patients. Age,",
      "weight and race are not reported in Zhu 2019; see the primary",
      "A.MARTIN publication (de Bono 2019, Clin Cancer Res 25:928)."
    )
  )

  ini({
    # ==================================================================
    # Zhu 2019 Table S1, 'AUC ss,event as exposure metrics':
    #   lambda(t) = lambda0(t) * exp(beta * AUC)
    # (equation printed beneath Table S1). Estimate -0.0000739, SE 0.0000619,
    # 95 percent CI -0.000195 to 0.0000473, P = 0.232. Fitted in R 3.2.3.
    #
    # NO BASELINE HAZARD IS ENCODED. lambda0(t) is the nonparametric
    # baseline of a Cox regression, a quantity the fit never reports, so
    # this model returns the relative hazard only (as
    # Chan_2025_atezolizumab_pfs.R and Gupta_2022_brigatinib_pfs.R do).
    # With AUC = 0 for placebo, hr is the rPFS hazard ratio versus the
    # placebo-plus-abiraterone arm that Zhu 2019 Figure 1 projects.
    # ==================================================================
    e_auc_haz <- -0.0000739; label("Log hazard ratio of radiographic progression per 1 ng*h/mL increase in AUC_IPAT_EVT (log scale)") # Zhu 2019 Table S1, 'AUC ss,event as exposure metrics', beta = -0.0000739 (SE 0.0000619, P 0.232)
  })

  model({
    lhr <- e_auc_haz * AUC_IPAT_EVT
    hr <- exp(lhr)
  })
}
