Warren_2025_orismilast <- function() {
  description <- paste(
    "In vitro (human whole blood). Sigmoidal Imax model of orismilast inhibition of",
    "IL-13 production in LPS-stimulated human whole blood, driven by an externally",
    "supplied orismilast plasma concentration (CP_ORISMILAST_NGML). Warren 2025",
    "measured the orismilast and apremilast IL-13 potencies head-to-head in the same",
    "assay and compared them against predicted plasma exposure in atopic dermatitis;",
    "orismilast is about 100-fold more potent than apremilast on this endpoint.",
    "PD-only model: the orismilast population PK model referenced by Warren 2025 is a",
    "company-internal two-compartment-with-lag analysis whose parameter values are not",
    "published in the paper, its reference list, or any supplement, so no orismilast",
    "PK model is packaged and users must supply their own concentration trajectory.",
    "The companion file Warren_2025_apremilast.R carries the apremilast arm, whose PK",
    "parameters the paper does report.",
    sep = " "
  )

  reference <- paste(
    "Warren RB, Weiss A, Felding J, Sommer MOA.",
    "Population Pharmacokinetic-Pharmacodynamic (popPK/PD) Relationship of Orismilast,",
    "A Potent and Selective PDE4B/D Inhibitor, in Atopic Dermatitis.",
    "Dermatol Ther (Heidelb). 2025;15(4):831-839. doi:10.1007/s13555-025-01371-9.",
    sep = " "
  )

  vignette <- "Warren_2025_pde4_il13_atopic_dermatitis"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external orismilast plasma-concentration covariate)",
    concentration = "(observation il13Inhibition is the fraction of the maximum achievable whole-blood IL-13 inhibition, 0-1; driving covariate CP_ORISMILAST_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_ORISMILAST_NGML = list(
      description = paste(
        "Instantaneous orismilast plasma concentration at the time of each PD",
        "observation, supplied as a time-varying covariate from observed plasma samples",
        "or an upstream PK source."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the sigmoidal Imax expression",
        "il13Inhibition = imax * CP_ORISMILAST_NGML^hill / (ec50^hill + CP_ORISMILAST_NGML^hill).",
        "Warren 2025 assessed the IL-13 potency in whole blood, so the paper applies",
        "plasma concentrations to this curve directly with no protein-binding or",
        "cellular-penetration adjustment (Warren 2025 Methods, 'Whole Blood Assay').",
        "Set to 0 outside the drug-exposure window; the inhibition then collapses to 0."
      ),
      source_name = "orismilast plasma concentration (Warren 2025 Fig. 1A, 1B)"
    )
  )

  population <- list(
    species = "in vitro (human whole blood)",
    n_subjects = 8L,
    n_studies = 1,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "Healthy blood donors (in vitro); the exposures the paper compares the curve against are from atopic dermatitis patients in the ADESOS phase 2b study (NCT05469464).",
    dose_range = "Eight in vitro orismilast concentrations; the paper contextualises the curve against 20, 30, and 40 mg oral twice-daily orismilast at steady state",
    regions = NA,
    notes = paste(
      "Fresh human whole blood (n = 8 healthy donors, Tissue Solutions) was plated in",
      "duplicate and stimulated with lipopolysaccharide (1 ug/mL) for 24 h at ~37 C and",
      "5% CO2 in the presence or absence of orismilast (MW 510.29 g/mol) at eight",
      "concentrations; IL-13 in the supernatant was measured with a Procartaplex Luminex",
      "kit. Three donors had a majority of points below the limit of quantification and",
      "one donor per treatment group was a non-responder; those were excluded, so the",
      "reported curve rests on fewer than eight donors (the exact analysed n is not",
      "stated). No cytotoxicity was observed. An asymmetric five-parameter curve was fit",
      "by least-squares regression to give the relative IC50, and the IC90 was computed",
      "from the fitted hillslope coefficient (Warren 2025 Methods, 'IL-13 Concentration",
      "in Human Whole Blood'). No between-donor variability was reported, so this model",
      "is typical-value only: it carries no IIV and no residual error."
    )
  )

  ini({
    # IL-13 inhibition potency -- LPS-stimulated human whole blood, Warren 2025 Table 2.
    # Table 2 prints the orismilast IC50 as "8 nM (4 ng/ml)". The ng/mL figure is the
    # paper's rounded conversion; 8 nM x 510.29 g/mol (MW from Warren 2025 Methods)
    # = 4.0823 ng/mL, which is the exact ng/mL equivalent of the reported molar value
    # and is what this model uses. Not re-estimated here, hence fixed().
    lec50 <- fixed(log(4.0823)); label("Orismilast concentration inhibiting 50% of whole-blood IL-13 production (IC50, ng/mL)")  # Warren 2025 Table 2: IC50 = 8 nM (4 ng/mL); 8e-9 mol/L * 510.29 g/mol = 4.0823 ng/mL

    limax <- fixed(log(1)); label("Maximum fractional inhibition of IL-13 production (unitless)")  # Warren 2025 Methods: the IC50 is estimated "relative to the maximum response achieved with each drug", so the fitted asymptote is 100% of the achievable inhibition

    # hill is not printed in Warren 2025; the paper states the IC90 was calculated from
    # the fitted hillslope coefficient. Under the sigmoidal Imax form used here the
    # published IC50/IC90 pair inverts exactly to the hillslope:
    #   hill = ln(9) / ln(IC90 / IC50) = ln(9) / ln(49 nM / 8 nM) = 1.2123
    # so this value reproduces BOTH published potency anchors by construction. The
    # molar pair is used for the inversion because it is the primary reported
    # quantity; the rounded ng/mL pair (25 / 4) would give 1.1985 instead.
    lhill <- fixed(log(1.2123)); label("Hill coefficient of the whole-blood IL-13 inhibition curve (unitless); derived from the published IC50/IC90 pair, not printed in the paper")  # derived from Warren 2025 Table 2: IC50 = 8 nM, IC90 = 49 nM
  })

  model({
    ec50 <- exp(lec50)
    imax <- exp(limax)
    hill <- exp(lhill)

    # Sigmoidal Imax inhibition of LPS-stimulated whole-blood IL-13 production,
    # expressed as the fraction of the maximum achievable inhibition (0-1).
    # CP_ORISMILAST_NGML is supplied per event row (ng/mL).
    il13Inhibition <- imax * CP_ORISMILAST_NGML^hill / (ec50^hill + CP_ORISMILAST_NGML^hill)
  })
}
