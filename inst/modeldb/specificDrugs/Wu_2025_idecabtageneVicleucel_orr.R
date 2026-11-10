Wu_2025_idecabtageneVicleucel_orr <- function() {
  description <- paste0(
    "Sigmoid-Emax-on-the-logit exposure-efficacy model for the probability ",
    "of overall response (partial response or better, independent review ",
    "committee, IMWG criteria) after a single infusion of the BCMA-directed ",
    "CAR T-cell therapy idecabtagene vicleucel (ide-cel) in adults with ",
    "triple-class-exposed relapsed/refractory multiple myeloma (Wu 2025; ",
    "ide-cel arm of the phase 3 KarMMa-3 trial, N = 220). The exposure is ",
    "the noncompartmental day-0-28 transgene-level AUC (AUC_IDECEL_DDPCR, ",
    "day*copies/ug, droplet digital PCR) on the natural log scale. The ",
    "baseline logit E0 is fixed at -5; Emax = 7.654 is the INCREMENT on the ",
    "logit, so the response probability plateaus at expit(2.654) = 93% for ",
    "the reference subject (female, Formulation A, two prior regimens, no ",
    "extramedullary plasmacytoma). The Hill coefficient of 34 on ",
    "log(AUC0-28D) makes the switch nearly a step at log(AUC0-28D) = 11.55 ",
    "(AUC0-28D about 1.0 x 10^5). Extramedullary plasmacytoma lowers the ",
    "logit by 1.34. No cellular kinetic layer and no ODE; no random effect ",
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
    "Parameter estimates in Supporting Information Table S8.",
    sep = " "
  )
  vignette <- "Wu_2025_idecabtageneVicleucel"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IDECEL_DDPCR covariate column)",
    concentration = "prob_orr_central (probability of overall response by independent review, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    AUC_IDECEL_DDPCR = list(
      description = "Individual ide-cel transgene-level AUC from the infusion to day 28 (AUC0-28D), by noncompartmental analysis of peripheral-blood CAR transgene copies per microgram of genomic DNA measured by droplet digital PCR.",
      units = "day*copies/ug",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the natural log, the argument of the sigmoid Emax. KarMMa-3 quartile boundaries (Figure 2a legend): 0.4, 1.3 and 2.7 x 10^6, maximum 17.8 x 10^6. ddPCR assay; the paper did not compare these values with the KarMMa (phase 2) assay.",
      source_name = "AUC0-28D (log-transformed in the model)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "Male vs female with female as reference; applied to (1 - SEXF). Pre-selected covariate.",
      source_name = "Sex (Male vs Female)"
    ),
    FORM_IDECEL_B = list(
      description = "1 = ide-cel manufactured as Formulation B, 0 = Formulation A (reference).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Pre-selected covariate; not significant (95% CI includes 0).",
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
    DIS_EMD = list(
      description = "1 = extramedullary plasmacytoma present at baseline, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Retained from the univariable screen; the only significant covariate in the ORR model.",
      source_name = "Baseline Extramedullary Plasmacytoma (Yes vs No)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 220L,
    n_studies = 1L,
    n_observations = "220 binary outcomes (180 responders, 81.8%; Table S6)",
    age_range = "30-81 years (median 63; 43% aged 65 or older)",
    weight_range = "41.1-144 kg (median 83)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 70.0, Black = 7.3, Asian = 2.3, Other = 0.9, Unknown = 19.5),
    disease_state = "Triple-class-exposed relapsed/refractory multiple myeloma after two to four prior regimens with progression within 60 days of the last regimen; 61.8% triple-class refractory, 23.6% extramedullary plasmacytoma",
    dose_range = "Single IV infusion; target 150-450 x 10^6 CAR+ T cells (up to 540 x 10^6 permitted)",
    regions = "International (KarMMa-3, NCT03651128)",
    notes = "Ide-cel arm of KarMMa-3, CK data cutoff 1 March 2022. Observed ORR by AUC0-28D quartile (Figure S7): about 56%, 89%, 91% and 94%."
  )

  ini({
    # Wu 2025 Table S8, 'Parameter Estimates of ORR E-R Full Model'.
    # logit(P) = E0 + Emax * X^H / (EC50^H + X^H) + covariate shifts, with
    # X = log(AUC0-28D). Emax, Hill and EC50 were estimated on the log
    # scale: the printed 'Estimate' is the back-transformed value and the
    # printed 95% CI is on the log scale (log(7.654) = 2.035 sits in
    # 1.9-2.17, log(34.006) = 3.527 in 1.86-5.2, log(11.551) = 2.447 in
    # 2.41-2.48). Emax is an INCREMENT on the logit: read as the logit at
    # full effect, the plateau would be expit(7.654) = 99.95% and could not
    # reproduce the about 91% plateau of the Figure S7 VPC (vignette).
    logite0 <- fixed(-5); label("Logit of the overall-response probability at zero drug effect, reference subject (unitless logit)") # Table S8 E0 -5 (fixed; Results: 'baseline effect (E0) fixed at -5')
    lemax <- log(7.654); label("Natural log of the maximal increase in the overall-response logit (unitless)") # Table S8 Emax 7.654 (95% CI on the log scale 1.9 to 2.17)
    lec50 <- log(11.551); label("Natural log of the log(AUC0-28D) at half-maximal effect, AUC0-28D in day*copies/ug (unitless)") # Table S8 EC50 11.551 (95% CI on the log scale 2.41 to 2.48)
    lhill <- log(34.006); label("Natural log of the Hill coefficient on log(AUC0-28D) (unitless)") # Table S8 'Hill slope' 34.006 (95% CI on the log scale 1.86 to 5.2)
    e_sexm_logit <- 0.121; label("Shift in the overall-response logit for male vs female sex (unitless logit)") # Table S8 'Sex (Male vs Female)' 0.121 (SE 0.468)
    e_form_idecel_b_logit <- -0.19; label("Shift in the overall-response logit for Formulation B vs A (unitless logit)") # Table S8 'Product Attributes (Formulation B vs A)' -0.19 (SE 0.52)
    e_line_4l_logit <- 0.03; label("Shift in the overall-response logit for 3 vs 2 prior MM regimens (unitless logit)") # Table S8 'No. of prior MM Regimens (3 vs 2)' 0.03 (SE 0.537)
    e_line_5l_logit <- -0.005; label("Shift in the overall-response logit for 4 vs 2 prior MM regimens (unitless logit)") # Table S8 'No. of prior MM Regimens (4 vs 2)' -0.005 (SE 0.554)
    e_emd_logit <- -1.343; label("Shift in the overall-response logit for baseline extramedullary plasmacytoma (unitless logit)") # Table S8 'Baseline Extramedullary Plasmacytoma (Yes vs No)' -1.343 (SE 0.483)

    # Bernoulli likelihood in the source: no between-subject variability and
    # no residual error. The tiny fixed additive residual lets rxode2 attach
    # an error model to the probability; it is not a published quantity.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    emax <- exp(lemax)
    hill <- exp(lhill)

    # Exposure on the natural-log scale (AUC0-28D in day*copies/ug).
    lnauc <- log(AUC_IDECEL_DDPCR)

    # Hill fraction lnauc^hill / (ec50^hill + lnauc^hill), written as an
    # expit of hill * (log(lnauc) - lec50) so the large Hill coefficient
    # never forms lnauc^hill directly.
    fdrug <- expit(hill * (log(lnauc) - lec50))

    logit_orr_central <- logite0 + emax * fdrug +
      e_sexm_logit * (1 - SEXF) +
      e_form_idecel_b_logit * FORM_IDECEL_B +
      e_line_4l_logit * LINE_4L +
      e_line_5l_logit * LINE_5L +
      e_emd_logit * DIS_EMD

    prob_orr_central <- expit(logit_orr_central)

    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
