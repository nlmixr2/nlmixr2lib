Connarn_2023_idecabtageneVicleucel_orr <- function() {
  description <- paste0(
    "Sigmoid-Emax-on-the-logit exposure-response model for the overall ",
    "response rate (ORR, Independent Response Committee assessment, IMWG ",
    "criteria) after a single infusion of the BCMA-directed CAR T-cell ",
    "therapy idecabtagene vicleucel (ide-cel) in adults with relapsed and ",
    "refractory multiple myeloma (Connarn 2023; phase II KarMMa, N = 125). ",
    "The exposure is the individual noncompartmental AUC of the transgene ",
    "level from day 0 to day 28 (AUC_IDECEL, day*copies/ug), entered on the ",
    "natural log scale; female sex adds to the logit. There is no cellular ",
    "kinetic layer and no ODE: the exposure is supplied as a data column. ",
    "No random effect and no residual error are estimated (Bernoulli ",
    "likelihood)."
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
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IDECEL covariate column)",
    concentration = "prob_orr_central (probability of an overall response, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    AUC_IDECEL = list(
      description = "Individual ide-cel area under the transgene-level curve from the infusion to day 28 (AUC0-28), by noncompartmental analysis of peripheral-blood transgene copies per microgram of genomic DNA.",
      units = "day*copies/ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters on the NATURAL LOG scale: the paper 'natural",
        "log-transformed' AUC0-28 for every ER evaluation and Figure 3 plots",
        "the model against log(AUC0-28) on an 11-18 axis (about 6e4 to",
        "6.6e7 day*copies/ug). Phoenix WinNonlin NCA; BLQ values set to",
        "zero; the day-28 value log-interpolated, or extrapolated with",
        "lambda-z when day 28 fell after the last numeric observation.",
        "Figure 3 per-target-dose medians (90% ranges) of log(AUC0-28):",
        "150 x 10^6 cells 14.06 (13.71-15.76, n = 4), 300 x 10^6 14.65",
        "(12.12-16.35, n = 68), 450 x 10^6 15.32 (13.47-16.47, n = 53),",
        "digitised by the maintainers."
      ),
      source_name = "AUC0-28 (log(AUC0-28) in the model)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Figure 3a table row 'Sex (female vs male)' = 1.78 on the logit scale, so female is the indicated level. Retained by stepwise covariate selection (p = 0.003).",
      source_name = "Sex (female vs male)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Screened (continuous and binned < 65 vs >= 65 years; Table 1) and not retained."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    ),
    BSA = list(
      description = "Baseline body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    ),
    SBCMA = list(
      description = "Serum soluble BCMA at baseline and before the ide-cel infusion",
      units = "(not stated)",
      type = "continuous",
      notes = "Screened (Table 1) and not retained."
    ),
    MCPROT = list(
      description = "Baseline serum monoclonal (M) protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Table 1) and not retained for ORR; it is retained, dichotomised at 10 g/L, in the complete-response model Connarn_2023_idecabtageneVicleucel_cr."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 125L,
    n_studies = 1L,
    n_observations = "125 binary overall-response outcomes (93 responders by Figure 3a quartile counts 10/32, 25/31, 29/31, 29/31)",
    age_range = "33-78 years (median 61; all 128 treated KarMMa patients)",
    disease_state = "Relapsed and refractory multiple myeloma after at least three prior regimens (median six); 84% triple-class refractory, 26% penta-refractory, 88% received bridging therapy",
    dose_range = "Single IV infusion of 150 (n = 4), 300 (n = 70) or 450 (n = 54) x 10^6 CAR+ T cells after fludarabine 30 mg/m^2/day plus cyclophosphamide 300 mg/m^2/day lymphodepletion",
    regions = "Multinational (KarMMa, NCT03361748)",
    notes = paste(
      "Data cutoff 19 April 2019 (>= 3 months PK follow-up). 127 of 128",
      "treated patients had a measurable Cmax and 125 an AUC0-28; Cmax and",
      "AUC0-28 correlated at r = 0.93, so the ER analyses used AUC0-28.",
      "The sex split is not reported in Connarn 2023. Covariates screened",
      "(Table 1): age, weight, BSA, race, sex, ethnicity, number of prior",
      "regimens, last prior therapy, bridging therapy, extramedullary",
      "disease, ECOG PS, prior HSCT, lymphodepleting regimen, refractory vs",
      "relapsed status, concomitant CRS medication, tumour burden, urine",
      "and serum M-protein, ferritin, soluble BCMA, free light chain, IL-6,",
      "IL-15, TNF-alpha, anti-drug antibody status and target dose."
    )
  )

  ini({
    # Figure 3a embedded table, 'Estimate (%RSE) on the logit scale', final
    # ORR model. Methods equation:
    #   logit(P) = E0 + (Emax - E0) * AUC^H / (EC50^H + AUC^H) + b2 * COV
    # with AUC = log(AUC0-28). Emax is the logit at full effect (not an
    # increment), so the drug term spans Emax - E0.
    #
    # The table prints 'Log(EC50)' and 'H'. Both are log-scale estimates:
    # EC50 = exp(2.64) = 14.0 in log(AUC0-28) units, at the visible midpoint
    # of the Figure 3a curve, and the Hill coefficient is exp(3.25) = 25.8.
    # Read literally as H = 3.25 on log(AUC0-28), the curve would be far too
    # shallow: about 39% ORR at log(AUC) = 11 against the plotted 14%, and
    # the CRS-tocilizumab sibling model could not reach its plotted 8.5%
    # floor anywhere on the axis. With exp(H), both models reproduce the
    # digitised Figure 3 curves to within about 1 percentage point and the
    # Table S1 dose-response predictions within about 3 points (vignette).
    logite0 <- -2.83; label("Logit of the overall-response probability at zero drug effect, male (unitless logit)") # Figure 3a table E0 -2.83 (43.0 %RSE)
    logitemax <- 2.29; label("Logit of the overall-response probability at maximal drug effect, male (unitless logit)") # Figure 3a table Emax 2.29 (37.0 %RSE)
    lec50 <- 2.64; label("Natural log of the log(AUC0-28) at half-maximal effect, with AUC0-28 in day*copies/ug (unitless)") # Figure 3a table 'Log(EC50)' 2.64 (0.96 %RSE)
    lhill <- 3.25; label("Natural log of the Hill coefficient on log(AUC0-28) (unitless)") # Figure 3a table 'H' 3.25 (19.7 %RSE); log scale, see block comment above
    e_sexf_logit <- 1.78; label("Shift in the overall-response logit for female vs male sex (unitless logit)") # Figure 3a table 'Sex (female vs male)' 1.78 (38.3 %RSE)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
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

    logit_orr_central <- logite0 + (logitemax - logite0) * fdrug +
      e_sexf_logit * SEXF

    prob_orr_central <- expit(logit_orr_central)

    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
