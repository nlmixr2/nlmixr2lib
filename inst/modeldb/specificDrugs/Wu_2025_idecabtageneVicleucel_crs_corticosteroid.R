Wu_2025_idecabtageneVicleucel_crs_corticosteroid <- function() {
  description <- paste0(
    "Linear logistic exposure-safety model for the probability of ",
    "cytokine release syndrome requiring corticosteroid treatment (sCRS) ",
    "after a single infusion of the BCMA-directed CAR T-cell therapy ",
    "idecabtagene vicleucel (ide-cel) in adults with triple-class-exposed ",
    "relapsed/refractory multiple myeloma (Wu 2025; ide-cel arm of the ",
    "phase 3 KarMMa-3 trial, N = 220). The exposure is the early expansion ",
    "rate Cmax/Tmax of the peripheral-blood transgene level ",
    "(EXPANSION_RATE_IDECEL, copies/ug/day, droplet digital PCR) on the ",
    "natural log scale, slope 0.506 per log unit. Only the pre-selected ",
    "covariates (formulation, number of prior regimens) are carried, none ",
    "significant. The covariates are not centred: the intercept is the ",
    "logit at an expansion rate of 1 copy/ug/day. No cellular kinetic layer ",
    "and no ODE; no random effect and no residual error (Bernoulli ",
    "likelihood)."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S11.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the EXPANSION_RATE_IDECEL covariate column)",
    concentration = "prob_crs_corticosteroid (probability of CRS requiring corticosteroids, 0-1; also logit_crs_corticosteroid)"
  )

  covariateData <- list(
    EXPANSION_RATE_IDECEL = list(
      description = "Individual ide-cel expansion rate, the peak peripheral-blood transgene level divided by the time to that peak (Cmax/Tmax), from noncompartmental analysis of CAR transgene copies per microgram of genomic DNA measured by droplet digital PCR.",
      units = "copies/ug/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the natural log, uncentred. KarMMa-3 quartiles (Figures S8-S9): Q1 0-0.4, Q2 0.4-1.3, Q3 1.4-3.2, Q4 3.4-11.6 x 10^4 copies/ug/day.",
      source_name = "Cmax/Tmax (Log(Cmax/Tmax) in the model)"
    ),
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table S11 labels this row 'Formulation (sLVV vs aLVV)' while every other table and the Results call it Formulation B vs A; read here as the same Formulation B indicator. Pre-selected covariate; not significant.",
      source_name = "Formulation (sLVV vs aLVV)"
    ),
    LINE_4L = list(
      description = "1 = exactly three prior multiple-myeloma regimens, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference level = 2 prior regimens. Pre-selected covariate; not significant.",
      source_name = "No. of prior MM Regimens (3 vs 2)"
    ),
    LINE_5L = list(
      description = "1 = exactly four prior multiple-myeloma regimens, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Reference level = 2 prior regimens. Pre-selected covariate; not significant.",
      source_name = "No. of prior MM Regimens (4 vs 2)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_observations = "220 binary outcomes (60 sCRS events, 27.3%; Table S9C and Figure 3b: 7, 11, 16 and 26 of 55 by expansion-rate quartile)",
    age_range = "30-81 years (median 63; 43% aged 65 or older)",
    weight_range = "41.1-144 kg (median 83)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens with progression within 60 days of the last regimen",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted)",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, CK data cutoff 1 March 2022. No covariate passed the univariable screen, so the full and final model is the pre-selected base model."
  )

  ini({
    # Wu 2025 Table S11, 'Parameter Estimates of Exposure-Response of sCRS
    # Base Model' (final = full = base, Results). Natural log of Cmax/Tmax,
    # uncentred.
    logit_ref <- -6.1186; label("sCRS logit at an expansion rate of 1 copy/ug/day, reference subject (unitless logit)") # Table S11 Intercept -6.1186 (SE 1.3524)
    e_expansion_rate_idecel_logit <- 0.506; label("Slope of the sCRS logit per unit natural log of the expansion rate in copies/ug/day (unitless logit)") # Table S11 'Log(Cmax/Tmax)' 0.506 (SE 0.1332)
    e_form_idecel_b_logit <- 0.277; label("Shift in the sCRS logit for Formulation B vs A (unitless logit)") # Table S11 'Formulation (sLVV vs aLVV)' 0.277 (SE 0.3652)
    e_line_4l_logit <- 0.4423; label("Shift in the sCRS logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S11 'No. of prior MM Regimens (3 vs 2)' 0.4423 (SE 0.3843)
    e_line_5l_logit <- 0.2028; label("Shift in the sCRS logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S11 'No. of prior MM Regimens (4 vs 2)' 0.2028 (SE 0.4109)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_corticosteroid <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_crs_corticosteroid <- logit_ref +
      e_expansion_rate_idecel_logit * log(EXPANSION_RATE_IDECEL) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L

    prob_crs_corticosteroid <- expit(logit_crs_corticosteroid)

    prob_crs_corticosteroid ~ add(addSd_prob_crs_corticosteroid)
  })
}
