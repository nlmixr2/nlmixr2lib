Nguyen_2021_imipenem <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in 44 Vietnamese",
    "adults hospitalised for acute exacerbations of chronic obstructive",
    "pulmonary disease (Nguyen 2021). Clearance scales as a power of",
    "Cockcroft-Gault creatinine clearance, the only covariate retained;",
    "the volume of distribution carries no covariate. Inter-individual",
    "variability is exponential on clearance and volume, and residual error",
    "is proportional.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Nguyen 2021 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Nguyen TM, Ngo TH, Truong AQ, Vu DH, Le DC, Vu NB, et al.",
    "Population pharmacokinetics and dose optimization of ceftazidime and",
    "imipenem in patients with acute exacerbations of chronic obstructive",
    "pulmonary disease.",
    "Pharmaceutics. 2021;13(4):456. doi:10.3390/pharmaceutics13040456.",
    "(The paper models ceftazidime as well as imipenem; only the imipenem",
    "model is tabulated by the review and extracted here.)",
    "Parameters transcribed from Zhang P, Zhao Y, Zhu J, Yang Y, Liang G,",
    "Wang X, Yu Z. Population pharmacokinetics of imipenem in different",
    "populations for individualized dosing: a systematic review.",
    "Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055",
    "(Tables 1-3, Supplementary Table S1).",
    sep = " "
  )
  vignette <- "Zhang_2025_imipenem_model_review"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE because the primary publication is
  # not on disk.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated by the Cockcroft-Gault equation.",
        "The review does not state whether the value is BSA-normalised;",
        "Cockcroft-Gault natively returns raw mL/min and the review's",
        "abbreviation list glosses 'CLcrCG' with no normalisation",
        "mentioned, so raw mL/min is used here."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 75.54 mL/min. Enters CL as the power term",
        "(CLcr/75.54)^0.532 (Zhang 2025 Table 3). The review does not",
        "report this cohort's renal function distribution, so whether",
        "75.54 is the cohort median cannot be checked from the secondary",
        "source. Zhang 2025 Results names this study as one of five whose",
        "final model retains creatinine clearance as its SOLE covariate.",
        "Note that the cohort is small (44 patients, 84 samples) and light",
        "(median 50 kg), so the term is fitted over a narrow range of",
        "body sizes. Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts raw",
        "Cockcroft-Gault mL/min when the source does not BSA-normalise --",
        "precedent: Bai 2024 imipenem, Wang 2024 imipenem."
      ),
      source_name        = "CLcr"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    AGE  = list(description = "Age",                 units = "years",   type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort median 65 years, IQR 60-72 (Table 1)."),
    SEXF = list(description = "Female sex",          units = "(binary)", type = "binary",    notes = "Screened, not retained (Zhang 2025 Table 3). Cohort is 41 male / 3 female (Table 1), so the term is near-unidentifiable in this study."),
    WT   = list(description = "Total body weight",   units = "kg",      type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort median 50 kg, IQR 47-55 (Table 1)."),
    BMI  = list(description = "Body mass index",     units = "kg/m^2",  type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3).")
  )

  population <- list(
    species          = "human",
    n_subjects       = 44L,
    n_studies        = 1L,
    age_median       = "65 years (IQR 60-72)",
    weight_median    = "50 kg (IQR 47-55)",
    sex_female_pct   = 6.8,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Adults hospitalised for acute exacerbations of chronic obstructive",
      "pulmonary disease (AECOPD) and treated with imipenem-cilastatin.",
      "This is the only AECOPD cohort among the 18 studies in the Zhang",
      "2025 review."
    ),
    dose_range       = paste(
      "500 mg imipenem intravenously every 6, 8 or 12 h, or 1000 mg every",
      "8 or 12 h (Zhang 2025 Supplementary Table S1). The infusion",
      "duration is not reported by the review."
    ),
    regions          = "Vietnam",
    n_concentrations = 84L,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 10); 44 patients, 84",
      "samples, sex split 41 male / 3 female. Blood was sampled 0.5 h after",
      "the infusion of the third dose and 1-2 h before the fourth dose, and",
      "assayed by HPLC-UV (Zhang 2025 Supplementary Table S1). Covariates",
      "were selected by forward inclusion (dOFV > 6.635, p < 0.01) and",
      "backward elimination (dOFV > 10.828, p < 0.001). In addition to the",
      "covariates listed in covariatesDataExcluded, the review records that",
      "the Anthonisen score, MDRD-4 creatinine clearance, respiratory",
      "distress and diuretic intake were screened and not retained; the",
      "first and third of those have no canonical column in",
      "inst/references/covariate-columns.md and are noted here rather than",
      "registered, since no coefficient is reported for any of them.",
      "Fitted in Monolix and evaluated by bootstrap, VPC, goodness-of-fit",
      "plots and NPDE (Zhang 2025 Table 2).",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Nguyen 2021 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Nguyen et al. (2021) row.
    # Typical values at the reference CLcr of 75.54 mL/min. =====
    lcl <- log(7.88); label("Clearance at CRCL = 75.54 mL/min (L/h)")  # Zhang 2025 Table 2 (Nguyen 2021): CL = 7.88 L/h; leading coefficient of the Table 3 CL formula
    lvc <- log(15.1); label("Volume of distribution (L)")              # Zhang 2025 Table 2 (Nguyen 2021): V = 15.1 L

    # ===== Covariate effect -- Zhang 2025 Table 3, Nguyen et al. (2021):
    #   CL = 7.88 * (CLcr/75.54)^0.532 * e^eta_CL
    e_crcl_cl <- 0.532; label("Power exponent on (CRCL/75.54) for CL (unitless)")  # Zhang 2025 Table 3 (Nguyen 2021): (CLcr/75.54)^0.532

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 prints IIV as a bare percentage
    # per parameter without stating the convention, and the column mixes at
    # least three conventions across the review's constituent studies (the
    # audit against the four already-extracted primaries is in the
    # vignette's 'Assumptions and deviations' section). Every model
    # transcribed from this review uses the same documented reading: the
    # printed percentage is an apparent CV of a log-normal random effect,
    # so omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.294^2)  # Zhang 2025 Table 2 (Nguyen 2021): IIV CL = 29.4%, read as an apparent CV
    etalvc ~ log(1 + 0.107^2)  # Zhang 2025 Table 2 (Nguyen 2021): IIV V = 10.7%, read as an apparent CV

    # ===== Residual error =====
    propSd <- 0.233; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Nguyen 2021): Proportional = 23.3%
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (CRCL / 75.54)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central) <- -kel * central

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
