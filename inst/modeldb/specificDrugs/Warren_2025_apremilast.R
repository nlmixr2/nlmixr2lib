Warren_2025_apremilast <- function() {
  description <- "One-compartment oral population PK model for apremilast with a sigmoidal Imax model of IL-13 inhibition in LPS-stimulated human whole blood, simulated in atopic dermatitis patients (Warren 2025)"
  reference <- paste(
    "Warren RB, Weiss A, Felding J, Sommer MOA.",
    "Population Pharmacokinetic-Pharmacodynamic (popPK/PD) Relationship of Orismilast,",
    "A Potent and Selective PDE4B/D Inhibitor, in Atopic Dermatitis.",
    "Dermatol Ther (Heidelb). 2025;15(4):831-839. doi:10.1007/s13555-025-01371-9.",
    "The apremilast PK parameters in Warren 2025 Table 1 are reproduced from the",
    "FDA Otezla NDA 206088Orig1s000 clinical pharmacology and biopharmaceutics review (2014).",
    sep = " "
  )
  vignette <- "Warren_2025_pde4_il13_atopic_dermatitis"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "apremilast", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apremilast", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_PSORIASIS = list(
      description = "Plaque psoriasis disease-state indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-psoriasis (atopic dermatitis, other disease, or missing disease status)",
      notes = paste(
        "Warren 2025 Table 1 reports CL/F = 9.26 L/h for the psoriasis cohort of the FDA Otezla popPK",
        "model together with a 1.09x factor 'if other disease or missing'. The atopic dermatitis cohort",
        "simulated in this paper is 'other disease', so the authors used CL/F = 9.26 * 1.09 = 10.09 L/h",
        "(Warren 2025 Table 1 footnote). Set DIS_PSORIASIS = 0 to reproduce the published atopic",
        "dermatitis simulations and DIS_PSORIASIS = 1 to recover the psoriasis reference clearance.",
        "Time-fixed per subject."
      ),
      source_name = "disease (psoriasis vs other disease or missing)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state = "Moderate-to-severe atopic dermatitis (simulated); the underlying popPK model was estimated in plaque psoriasis patients.",
    dose_range = "30 mg and 40 mg oral twice daily (11 doses to steady state), Warren 2025 Methods",
    regions = NA,
    notes = paste(
      "Warren 2025 did not fit a new apremilast PK model. Plasma concentrations in female atopic dermatitis",
      "patients were simulated with PKPDsim from the population mean parameters of the FDA-approved final",
      "apremilast (Otezla) popPK model in psoriasis patients, tabulated in Warren 2025 Table 1. The paper",
      "states explicitly that interindividual variability of those parameters has not been published, so the",
      "simulation is a typical-value (deterministic) one and this model file carries no IIV and no residual",
      "error. The IL-13 pharmacodynamic layer comes from a head-to-head LPS-stimulated human whole-blood",
      "assay in n = 8 healthy donors (3 excluded for data below the limit of quantification or non-response),",
      "reported in Warren 2025 Table 2."
    )
  )

  ini({
    # Structural PK -- Warren 2025 Table 1 (population mean parameters of the FDA
    # Otezla popPK model in psoriasis). Transferred without re-estimation, hence fixed().
    lka <- fixed(log(1.84)); label("Absorption rate constant (1/h)")                    # Warren 2025 Table 1: Ka = 1.84 1/h
    lcl <- fixed(log(9.26)); label("Apparent clearance in the psoriasis reference cohort (CL/F, L/h)")  # Warren 2025 Table 1: CL/F = 9.26 L/h
    lvc <- fixed(log(118)); label("Apparent central volume of distribution (Vc/F, L)")  # Warren 2025 Table 1: Vc/F = 118 L

    # Covariate effect -- multiplicative factor applied to CL/F when the subject is
    # NOT a psoriasis patient. Warren 2025 applied it to reach the 10.09 L/h used for
    # the atopic dermatitis simulations.
    e_dis_psoriasis_cl <- fixed(1.09); label("Multiplicative factor on CL/F for other disease or missing disease status vs the psoriasis reference (unitless)")  # Warren 2025 Table 1 row "If other disease or missing" and table footnote

    # IL-13 pharmacodynamics -- LPS-stimulated human whole blood, Warren 2025 Table 2.
    lec50 <- fixed(log(405)); label("Apremilast concentration inhibiting 50% of whole-blood IL-13 production (IC50, ng/mL)")  # Warren 2025 Table 2: 881 nM (405 ng/mL)
    limax <- fixed(log(1)); label("Maximum fractional inhibition of IL-13 production (unitless)")  # Warren 2025 Methods: IC50 estimated "relative to the maximum response achieved with each drug", i.e. the fitted asymptote is 100%
    # hill is not printed in Warren 2025; the paper states the IC90 was calculated from
    # the fitted hillslope coefficient. Under the sigmoidal Imax form used here the
    # published IC50/IC90 pair inverts exactly to the hillslope:
    #   hill = ln(9) / ln(IC90 / IC50) = ln(9) / ln(163000 nM / 881 nM) = 0.4209
    # so this value reproduces BOTH published potency anchors by construction.
    lhill <- fixed(log(0.4209)); label("Hill coefficient of the whole-blood IL-13 inhibition curve (unitless); derived from the published IC50/IC90 pair, not printed in the paper")  # derived from Warren 2025 Table 2: IC50 = 881 nM, IC90 = 163 uM
  })

  model({
    # Individual parameters. The paper's simulation is typical-value only, so no etas.
    ka <- exp(lka)
    # Warren 2025 Table 1: CL/F is 9.26 L/h in psoriasis and is multiplied by 1.09 for
    # "other disease or missing". Atopic dermatitis is "other disease", so
    # DIS_PSORIASIS = 0 gives the 10.09 L/h the authors simulated with.
    cl <- exp(lcl) * e_dis_psoriasis_cl^(1 - DIS_PSORIASIS)
    vc <- exp(lvc)

    ec50 <- exp(lec50)
    imax <- exp(limax)
    hill <- exp(lhill)

    kel <- cl / vc

    # Warren 2025 Eqs. 1-2 (one-compartment model with first-order absorption):
    #   dA[1]/dt = -Ka * A[1]
    #   dA[2]/dt =  Ka * A[1] - (CL / V) * A[2]
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 1000 converts mg/L to the ng/mL scale used for the published concentrations,
    # the IL-13 IC50, and Warren 2025 Figs. 1C-1D.
    Cc <- 1000 * central / vc

    # Sigmoidal Imax inhibition of LPS-stimulated whole-blood IL-13 production,
    # expressed as the fraction of the maximum achievable inhibition (0-1).
    il13Inhibition <- imax * Cc^hill / (ec50^hill + Cc^hill)
  })
}
