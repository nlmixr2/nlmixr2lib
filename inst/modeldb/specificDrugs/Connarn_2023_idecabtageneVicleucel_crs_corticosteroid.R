Connarn_2023_idecabtageneVicleucel_crs_corticosteroid <- function() {
  description <- paste0(
    "Linear logistic exposure-safety model for the probability of cytokine ",
    "release syndrome (CRS) requiring corticosteroids after a single ",
    "infusion of the BCMA-directed CAR T-cell therapy idecabtagene ",
    "vicleucel (ide-cel) in adults with relapsed and refractory multiple ",
    "myeloma (Connarn 2023; phase II KarMMa, N = 125). logit(P) = -11.1 + ",
    "0.61 * log(AUC_IDECEL), with AUC_IDECEL the noncompartmental ",
    "transgene-level AUC0-28 in day*copies/ug. No covariate was retained. ",
    "There is no cellular kinetic layer and no ODE. No random effect and no ",
    "residual error are estimated (Bernoulli likelihood)."
  )
  reference <- paste(
    "Connarn JN, Witjes H, van Zutphen-van Geffen M, de Greef R,",
    "Campbell TB, Hege K, Zhou S, Lamba M.",
    "Characterizing the exposure-response relationship of idecabtagene",
    "vicleucel in patients with relapsed/refractory multiple myeloma.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1687-1697.",
    "doi:10.1002/psp4.12922.",
    sep = " "
  )
  vignette <- "Connarn_2023_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IDECEL covariate column)",
    concentration = "prob_crs_corticosteroid (probability of CRS requiring corticosteroids, 0-1; also logit_crs_corticosteroid)"
  )

  covariateData <- list(
    AUC_IDECEL = list(
      description = "Individual ide-cel area under the transgene-level curve from the infusion to day 28 (AUC0-28), by noncompartmental analysis of peripheral-blood transgene copies per microgram of genomic DNA.",
      units = "day*copies/ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters on the NATURAL LOG scale (Figure 3d slope 0.61 per unit of",
        "log(AUC0-28)); see Connarn_2023_idecabtageneVicleucel_orr for the",
        "NCA conventions and the per-dose exposure ranges."
      ),
      source_name = "AUC0-28 (log(AUC0-28) in the model)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Screened (Table 1); no covariate was retained for CRS requiring corticosteroids."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1) and not retained."
    ),
    MCPROT = list(
      description = "Baseline serum monoclonal (M) protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    ),
    SBCMA = list(
      description = "Serum soluble BCMA at baseline and before the ide-cel infusion",
      units = "(not stated)",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 125L,
    n_studies = 1L,
    n_observations = "125 binary outcomes (18 events by Figure 3d quartile counts 1/32, 3/31, 7/31, 7/31)",
    age_range = "33-78 years (median 61; all 128 treated KarMMa patients)",
    disease_state = "Relapsed and refractory multiple myeloma after at least three prior regimens (median six); 84% triple-class refractory, 26% penta-refractory, 88% received bridging therapy",
    dose_range = "Single IV infusion of 150 (n = 4), 300 (n = 70) or 450 (n = 54) x 10^6 CAR+ T cells after fludarabine 30 mg/m^2/day plus cyclophosphamide 300 mg/m^2/day lymphodepletion",
    regions = "Multinational (KarMMa, NCT03361748)",
    notes = paste(
      "Data cutoff 19 April 2019. An Emax model could not be fitted to this",
      "endpoint, so the linear model was selected although it keeps rising",
      "where the observed Q3 and Q4 rates are identical (7/31)."
    )
  )

  ini({
    # Figure 3d embedded table, 'Estimate (%RSE) on the logit scale', final
    # linear model: logit(P) = b0 + b1 * log(AUC0-28). Uncentred, so the
    # intercept is the logit at AUC0-28 = 1 day*copies/ug.
    logit_ref <- -11.1; label("Logit of the probability of CRS requiring corticosteroids at log(AUC0-28) = 0 (unitless logit)") # Figure 3d table Intercept -11.1 (37.2 %RSE)
    e_auc_idecel_logit <- 0.61; label("Slope of the logit of CRS requiring corticosteroids per unit natural log of AUC0-28 in day*copies/ug (unitless logit)") # Figure 3d table Slope 0.61 (43.5 %RSE)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_corticosteroid <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_crs_corticosteroid <- logit_ref + e_auc_idecel_logit * log(AUC_IDECEL)

    prob_crs_corticosteroid <- expit(logit_crs_corticosteroid)

    prob_crs_corticosteroid ~ add(addSd_prob_crs_corticosteroid)
  })
}
