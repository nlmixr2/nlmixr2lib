Connarn_2023_idecabtageneVicleucel_crs_tocilizumab <- function() {
  description <- paste0(
    "Sigmoid-Emax-on-the-logit exposure-safety model for the probability of ",
    "cytokine release syndrome (CRS) requiring tocilizumab after a single ",
    "infusion of the BCMA-directed CAR T-cell therapy idecabtagene ",
    "vicleucel (ide-cel) in adults with relapsed and refractory multiple ",
    "myeloma (Connarn 2023; phase II KarMMa, N = 125). The exposure is the ",
    "noncompartmental transgene-level AUC0-28 (AUC_IDECEL, day*copies/ug) on ",
    "the natural log scale; the rate plateaus at about 61% above a ",
    "log(AUC0-28) of about 14. No covariate was retained. There is no ",
    "cellular kinetic layer and no ODE. No random effect and no residual ",
    "error are estimated (Bernoulli likelihood)."
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
    concentration = "prob_crs_tocilizumab (probability of CRS requiring tocilizumab, 0-1; also logit_crs_tocilizumab)"
  )

  covariateData <- list(
    AUC_IDECEL = list(
      description = "Individual ide-cel area under the transgene-level curve from the infusion to day 28 (AUC0-28), by noncompartmental analysis of peripheral-blood transgene copies per microgram of genomic DNA.",
      units = "day*copies/ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters on the NATURAL LOG scale (Figure 3c plotted against",
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
      notes = "Screened (Table 1); no covariate was retained for CRS requiring tocilizumab."
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
    n_observations = "125 binary outcomes (65 events by Figure 3c quartile counts 8/32, 20/31, 18/31, 19/31)",
    age_range = "33-78 years (median 61; all 128 treated KarMMa patients)",
    disease_state = "Relapsed and refractory multiple myeloma after at least three prior regimens (median six); 84% triple-class refractory, 26% penta-refractory, 88% received bridging therapy",
    dose_range = "Single IV infusion of 150 (n = 4), 300 (n = 70) or 450 (n = 54) x 10^6 CAR+ T cells after fludarabine 30 mg/m^2/day plus cyclophosphamide 300 mg/m^2/day lymphodepletion",
    regions = "Multinational (KarMMa, NCT03361748)",
    notes = paste(
      "Data cutoff 19 April 2019. The sigmoid Emax model was preferred over",
      "the linear logistic model (AIC 81.0 vs 83.0) because it reproduces",
      "the observed plateau at medium and high exposure. In the phase I",
      "CRB-401 study the observed rates (11.1% at 150 and 36.8% at",
      "450 x 10^6 cells) were well below these KarMMa-based predictions",
      "(Table S1), which the authors attribute to population, dosing and",
      "evolving tocilizumab practice."
    )
  )

  ini({
    # Figure 3c embedded table, 'Estimate (%RSE) on the logit scale', final
    # model for CRS requiring tocilizumab. Same equation and the same
    # log-scale reading of 'Log(EC50)' and 'H' as
    # Connarn_2023_idecabtageneVicleucel_orr: EC50 = exp(2.61) = 13.6 in
    # log(AUC0-28) units and Hill = exp(4.67) = 106.7. The plotted Figure 3c
    # curve is flat at expit(-2.37) = 8.5% up to log(AUC) ~ 13 and at
    # expit(0.44) = 60.8% from ~14, a step that H = 4.67 on log(AUC0-28)
    # cannot produce (it would give 17% at log(AUC) = 11 and 46% at 18).
    logite0 <- -2.37; label("Logit of the probability of CRS requiring tocilizumab at zero drug effect (unitless logit)") # Figure 3c table E0 -2.37 (43.5 %RSE)
    logitemax <- 0.44; label("Logit of the probability of CRS requiring tocilizumab at maximal drug effect (unitless logit)") # Figure 3c table Emax 0.44 (47.3 %RSE)
    lec50 <- 2.61; label("Natural log of the log(AUC0-28) at half-maximal effect, with AUC0-28 in day*copies/ug (unitless)") # Figure 3c table 'Log(EC50)' 2.61 (0.71 %RSE)
    lhill <- 4.67; label("Natural log of the Hill coefficient on log(AUC0-28) (unitless)") # Figure 3c table 'H' 4.67 (27.9 %RSE); log scale, see block comment above

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_tocilizumab <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Exposure on the natural-log scale (AUC0-28 in day*copies/ug).
    lnauc <- log(AUC_IDECEL)

    # Hill fraction lnauc^hill / (ec50^hill + lnauc^hill), written as an
    # expit of hill * (log(lnauc) - lec50) so the large Hill coefficient
    # never forms lnauc^hill directly.
    fdrug <- expit(hill * (log(lnauc) - lec50))

    logit_crs_tocilizumab <- logite0 + (logitemax - logite0) * fdrug

    prob_crs_tocilizumab <- expit(logit_crs_tocilizumab)

    prob_crs_tocilizumab ~ add(addSd_prob_crs_tocilizumab)
  })
}
