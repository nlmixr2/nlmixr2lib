Wu_2025_idecabtageneVicleucel_crs_anti_il6 <- function() {
  description <- paste0(
    "Linear logistic exposure-safety model for the probability of ",
    "cytokine release syndrome requiring anti-IL-6 treatment (tocilizumab ",
    "or siltuximab; tCRS) after a single infusion of the BCMA-directed CAR ",
    "T-cell therapy idecabtagene vicleucel (ide-cel) in adults with ",
    "triple-class-exposed relapsed/refractory multiple myeloma (Wu 2025; ",
    "ide-cel arm of the phase 3 KarMMa-3 trial, N = 220). The exposure is ",
    "the early expansion rate Cmax/Tmax of the peripheral-blood transgene ",
    "level (EXPANSION_RATE_IDECEL, copies/ug/day, droplet digital PCR) on ",
    "the natural log scale, slope 0.6551 per log unit. Post-infusion ADA ",
    "positivity lowers the logit by 0.88; formulation and the number of ",
    "prior regimens were pre-selected and not significant. The covariates ",
    "are not centred: the intercept is the logit at an expansion rate of ",
    "1 copy/ug/day. No cellular kinetic layer and no ODE; no random effect ",
    "and no residual error (Bernoulli likelihood)."
  )
  reference <- paste(
    "Wu F, Zheng X, Burnett J, Masilamani M, Zhang W, Zhong X, Caia A,",
    "Cook M, Piasecki J, Kondic A, Lamba M, Zhou J.",
    "Characterizing cellular expansion of idecabtagene vicleucel and",
    "association with clinical efficacy and safety in patients with",
    "triple-class-exposed relapsed/refractory multiple myeloma.",
    "J Clin Pharmacol. 2025;65(11):1568-1576.",
    "doi:10.1002/jcph.70075.",
    "Parameter estimates in Supporting Information Table S10.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the EXPANSION_RATE_IDECEL covariate column)",
    concentration = "prob_crs_anti_il6 (probability of CRS requiring tocilizumab or siltuximab, 0-1; also logit_crs_anti_il6)"
  )

  covariateData <- list(
    EXPANSION_RATE_IDECEL = list(
      description = "Individual ide-cel expansion rate, the peak peripheral-blood transgene level divided by the time to that peak (Cmax/Tmax), from noncompartmental analysis of CAR transgene copies per microgram of genomic DNA measured by droplet digital PCR.",
      units = "copies/ug/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the natural log, uncentred. KarMMa-3 quartiles (Figures S8-S9): Q1 0-0.4, Q2 0.4-1.3, Q3 1.4-3.2, Q4 3.4-11.6 x 10^4 copies/ug/day. Chosen for the safety endpoints as an early exposure measure (median time to CRS onset 1 day).",
      source_name = "Cmax/Tmax (Log(Cmax/Tmax) in the model)"
    ),
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-selected covariate; not significant.",
      source_name = "Product Attributes (Formulation B vs A)"
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
    ),
    ADA_POS = list(
      description = "1 = anti-CAR antidrug antibody positive at at least one post-infusion visit, 0 = negative.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "POST-BASELINE covariate (61.4% positive), measured mostly from about month 3, long after CRS onset (median day 1); the authors state no clinical relevance can be inferred from the lower tCRS odds.",
      source_name = "Post-infusion ADA Positive (Yes vs No)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_observations = "220 binary outcomes (156 tCRS events, 70.9%; Table S9C and Figure 3a: 28, 37, 43 and 48 of 55 by expansion-rate quartile)",
    age_range = "30-81 years (median 63; 43% aged 65 or older)",
    weight_range = "41.1-144 kg (median 83)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens with progression within 60 days of the last regimen",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted)",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, CK data cutoff 1 March 2022. The event is defined by management (anti-IL-6 given), not by grade."
  )

  ini({
    # Wu 2025 Table S10, 'Parameter Estimates of E-R of tCRS Full Model'
    # (final = full, Results). Natural log of Cmax/Tmax, uncentred: with
    # the Figure S8 quartile mid-points this reproduces the observed
    # quartile rates (vignette).
    logit_ref <- -4.5924; label("tCRS logit at an expansion rate of 1 copy/ug/day, reference subject (unitless logit)") # Table S10 Intercept -4.5924 (SE 1.1751)
    e_expansion_rate_idecel_logit <- 0.6551; label("Slope of the tCRS logit per unit natural log of the expansion rate in copies/ug/day (unitless logit)") # Table S10 'Log(Cmax/Tmax)' 0.6551 (SE 0.13)
    e_form_idecel_b_logit <- 0.0972; label("Shift in the tCRS logit for Formulation B vs A (unitless logit)") # Table S10 'Product Attributes (Formulation B vs A)' 0.0972 (SE 0.3682)
    e_line_4l_logit <- -0.1757; label("Shift in the tCRS logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S10 'No. of prior MM Regimens (3 vs 2)' -0.1757 (SE 0.3823)
    e_line_5l_logit <- 0.3399; label("Shift in the tCRS logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S10 'No. of prior MM Regimens (4 vs 2)' 0.3399 (SE 0.4216)
    e_ada_pos_logit <- -0.8765; label("Shift in the tCRS logit for post-infusion ADA positivity (unitless logit)") # Table S10 'Post-infusion ADA Positive (Yes vs No)' -0.8765 (SE 0.3668)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_crs_anti_il6 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_crs_anti_il6 <- logit_ref +
      e_expansion_rate_idecel_logit * log(EXPANSION_RATE_IDECEL) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L +
      e_ada_pos_logit * ADA_POS

    prob_crs_anti_il6 <- expit(logit_crs_anti_il6)

    prob_crs_anti_il6 ~ add(addSd_prob_crs_anti_il6)
  })
}
