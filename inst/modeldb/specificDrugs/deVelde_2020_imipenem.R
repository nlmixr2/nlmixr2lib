deVelde_2020_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 26 Dutch",
    "critically ill adults (de Velde 2020), parameterised in",
    "micro-constants rather than clearances: an elimination rate constant,",
    "two distribution rate constants and a central volume. The elimination",
    "rate constant scales as a power of absolute (de-normalised) CKD-EPI",
    "eGFR. Inter-individual variability is exponential on the elimination",
    "rate constant only, and residual error is proportional.",
    "This file encodes the NONMEM arm of the paper. The paper also fitted",
    "the same data non-parametrically in Pmetrics, but the review reports",
    "that arm's covariate model only in terms of per-subject medians",
    "(Ke_i,med and Ke(cov)i,med), which are not population parameters and",
    "cannot be reconstructed; see the population notes.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against de Velde 2020 when the primary",
    "is obtained.",
    sep = " "
  )
  reference <- paste(
    "de Velde F, de Winter BCM, Neely MN, Yamada WM, Koch BCP,",
    "Harbarth S, et al.",
    "Population pharmacokinetics of imipenem in critically ill patients: a",
    "parametric and nonparametric model converge on CKD-EPI estimated",
    "glomerular filtration rate as an impactful covariate.",
    "Clin Pharmacokinet. 2020;59(7):885-898.",
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
    central     = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Glomerular filtration rate estimated by the CKD-EPI equation and",
        "then DE-NORMALISED to an absolute per-patient rate in mL/min by",
        "multiplying the BSA-normalised value by the patient's body surface",
        "area. Zhang 2025's abbreviation list defines the symbol used in",
        "the formula, eGFR CKD-EPI-abs, as 'absolute CKD-EPI (i.e.,",
        "CKD-EPI, multiplied by BSA)'."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNIT HAZARD. This column is an ABSOLUTE rate in mL/min, NOT the",
        "mL/min/1.73 m^2 that a CKD-EPI calculator returns by default. The",
        "reference of 119 is only interpretable on the absolute scale --",
        "it corresponds to roughly 105 mL/min/1.73 m^2 at a typical 1.96",
        "m^2 adult BSA. Supplying a BSA-normalised eGFR here understates",
        "the ratio for any patient larger than 1.73 m^2. Enters the",
        "elimination rate constant as the power term",
        "(eGFR_CKD-EPI-abs/119)^0.655 (Zhang 2025 Table 3).",
        "Zhang 2025 Table 3 records that four renal-function estimators",
        "were screened head to head for this study -- Cockcroft-Gault,",
        "MDRD-4, CKD-EPI and Jelliffe -- alongside three body-size",
        "descriptors, and that the absolute CKD-EPI form won. Stored under",
        "the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts CKD-EPI eGFR",
        "-- precedent: Bajaj 2017 nivolumab, Krens 2020 ganciclovir."
      ),
      source_name        = "eGFR CKD-EPI-abs"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    WT  = list(description = "Total body weight",  units = "kg", type = "continuous", notes = "Screened as TBW, not retained (Zhang 2025 Table 3). Cohort median 75 kg, range 66-85 (Table 1)."),
    IBW = list(description = "Ideal body weight",  units = "kg", type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    LBW = list(description = "Lean body weight",   units = "kg", type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Not a registered canonical in inst/references/covariate-columns.md; recorded here as documentation only, since it is never referenced in model().")
  )

  population <- list(
    species          = "human",
    n_subjects       = 26L,
    n_studies        = 1L,
    age_median       = "51 years (IQR 39-54)",
    weight_median    = "75 kg (IQR 66-85)",
    sex_female_pct   = 30.8,
    race_ethnicity   = NULL,
    disease_state    = "Critically ill adults receiving imipenem-cilastatin",
    dose_range       = paste(
      "500 mg imipenem intravenously every 6 h (Zhang 2025 Supplementary",
      "Table S1). The infusion duration is not reported by the review."
    ),
    regions          = "Netherlands",
    n_concentrations = 138L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 9); 26 patients, 138",
      "samples, sex split 18 male / 8 female. Blood was sampled at peak,",
      "intermediate and trough time points and assayed by HPLC-UV (Zhang",
      "2025 Supplementary Table S1). Covariates were selected by forward",
      "inclusion (dOFV > 3.84, p < 0.05) and backward elimination",
      "(dOFV > 10.83, p < 0.001). Evaluated by VPC, NPC, goodness-of-fit",
      "plots and bootstrap (Zhang 2025 Table 2).",
      "ONLY THE NONMEM ARM IS ENCODED. The study fitted the same data twice,",
      "parametrically in NONMEM and non-parametrically in Pmetrics, and",
      "Zhang 2025 Table 2 reports both. The Pmetrics arm's structural",
      "values are available (Ke 0.681 /h, Kcp 0.374 /h, Kpc 0.495 /h, Vc",
      "31.1 L; IIV Ke 34%, Kcp 81.2%, Kpc 72%, Vc 42.6%; gamma residual",
      "error 3.4), but its covariate model is printed as",
      "Ke = Ke_i,med * ((eGFR CKD-EPI-abs)/119)^Ke(cov)i,med -- an",
      "expression in the individual subject's own posterior medians, which",
      "are not population parameters and are not tabulated. A",
      "non-parametric population is in any case a discrete support-point",
      "distribution, not the log-normal random-effect structure nlmixr2's",
      "ini() expresses, so the Pmetrics arm cannot be faithfully encoded",
      "as a sibling model file even in principle. Recorded in the vignette",
      "Errata. Note also that the review's gamma of 3.4 is a Pmetrics",
      "error-polynomial multiplier, not a residual standard deviation.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from de Velde 2020 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, de Velde et al. (2020) row,
    # 'NONMEM' block. This model is parameterised in MICRO-CONSTANTS, not
    # in clearances: the review reports Ke, Kcp and Kpc directly and gives
    # only the central volume. Clearances are therefore derived quantities
    # here (cl = kel * vc etc.) rather than estimated parameters, which is
    # the reverse of every other model in this review. =====
    lkel <- log(0.637); label("Elimination rate constant at eGFR = 119 mL/min (1/h)")  # Zhang 2025 Table 2 (de Velde 2020, NONMEM): Ke = 0.637 /h; leading coefficient of the Table 3 Ke formula
    lk12 <- log(0.166); label("Central-to-peripheral rate constant Kcp (1/h)")         # Zhang 2025 Table 2 (de Velde 2020, NONMEM): Kcp = 0.166 /h
    lk21 <- log(0.195); label("Peripheral-to-central rate constant Kpc (1/h)")         # Zhang 2025 Table 2 (de Velde 2020, NONMEM): Kpc = 0.195 /h
    lvc  <- log(29.6);  label("Central volume of distribution (L)")                    # Zhang 2025 Table 2 (de Velde 2020, NONMEM): Vc = 29.6 L

    # ===== Covariate effect -- Zhang 2025 Table 3, de Velde et al. (2020)
    # row, 'NONMEN' block (the review's own typo for NONMEM):
    #   Ke = 0.637 * ((eGFR CKD-EPI-abs)/119)^0.655 * e^eta
    # Note the covariate acts on the RATE CONSTANT, not on clearance. Since
    # the central volume carries no covariate, this is observationally the
    # same as a power effect on CL, but the distinction matters if a user
    # later adds a volume covariate.
    e_crcl_kel <- 0.655; label("Power exponent on (CRCL/119) for kel (unitless)")  # Zhang 2025 Table 3 (de Velde 2020, NONMEM): ((eGFR CKD-EPI-abs)/119)^0.655

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 prints IIV as a bare percentage
    # per parameter without stating the convention, and the column mixes at
    # least three conventions across the review's constituent studies (the
    # audit against the four already-extracted primaries is in the
    # vignette's 'Assumptions and deviations' section). Every model
    # transcribed from this review uses the same documented reading: the
    # printed percentage is an apparent CV of a log-normal random effect,
    # so omega^2 = log(1 + CV^2).
    #
    # In the NONMEM arm IIV is reported on Ke only; the review gives no
    # omega for Kcp, Kpc or Vc in that arm (the 81.2% / 72% / 42.6% figures
    # in the same Table 2 cell belong to the Pmetrics arm, which is not
    # encoded here -- see population$notes).
    etalkel ~ log(1 + 0.19^2)  # Zhang 2025 Table 2 (de Velde 2020, NONMEM): IIV Ke = 19%, read as an apparent CV

    # ===== Residual error =====
    # Proportional only, from the NONMEM arm. The Pmetrics arm's 'Gamma:
    # 3.4' is a non-parametric error-polynomial multiplier and is not a
    # residual standard deviation; it is not used here.
    propSd <- 0.348; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (de Velde 2020, NONMEM): Proportional = 34.8%
  })

  model({
    # ----- Individual PK parameters -----
    # Micro-constants are the estimated quantities for this model.
    kel <- exp(lkel + etalkel) * (CRCL / 119)^e_crcl_kel
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc  <- exp(lvc)

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
