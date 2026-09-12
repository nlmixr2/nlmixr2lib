Dinh_2022_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 24 Vietnamese",
    "critically ill adults (Dinh 2022). Clearance rises exponentially with",
    "Cockcroft-Gault creatinine clearance -- a log-linear rather than a",
    "power form, and uncentred, so the tabulated clearance is the value at",
    "zero creatinine clearance. Volumes and intercompartmental clearance",
    "carry no covariates. Inter-individual variability is exponential on",
    "all four disposition parameters and residual error is combined",
    "proportional plus additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Dinh 2022 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Dinh TD, Nguyen HN, Le BH, Nguyen TT, Nguyen HL.",
    "Population-based pharmacokinetics and dose optimization of imipenem in",
    "Vietnamese critically-ill patients.",
    "Infect Drug Resist. 2022;15:4575-4583. doi:10.2147/IDR.S373348.",
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
        "Creatinine clearance calculated by the Cockcroft-Gault equation.",
        "The review's abbreviation list glosses 'CLcrCG' with no",
        "normalisation mentioned, so raw mL/min is used here."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNCENTRED EXPONENTIAL FORM. Zhang 2025 Table 3 prints",
        "CL = 4.79 * e^(0.00642 * CLcr) with no centring term, so the",
        "tabulated typical clearance of 4.79 L/h is the value at",
        "CLcr = 0 mL/min, not at any cohort-typical renal function. This",
        "matters for interpretation: at a typical adult CLcr of 100",
        "mL/min the model gives 4.79 * e^0.642 = 9.11 L/h, which is the",
        "number comparable with the other adult models in the review, and",
        "the 4.79 L/h that Zhang 2025's abstract quotes as the lowest",
        "adult clearance in the review is therefore an extrapolated",
        "intercept rather than a clearance any patient exhibits. The same",
        "log-linear uncentred idiom is registered for Peng 2025 meropenem",
        "in inst/references/covariate-columns.md; note that the",
        "coefficient carries units of per mL/min and is NOT a",
        "dimensionless exponent.",
        "The review reports no renal-function distribution for this",
        "cohort, so the fitted range of the term is unknown and",
        "extrapolation is unguarded. Stored under the canonical CRCL",
        "column."
      ),
      source_name        = "CLcr"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    AGE  = list(description = "Age",               units = "years",    type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort mean 57.5 +/- 19.9 years (Table 1)."),
    SEXF = list(description = "Female sex",        units = "(binary)", type = "binary",     notes = "Screened, not retained (Zhang 2025 Table 3). Cohort 6/24 female (Table 1)."),
    WT   = list(description = "Actual body weight", units = "kg",      type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort mean 51.3 +/- 8.6 kg (Table 1)."),
    ALB  = list(description = "Serum albumin",     units = "g/L",      type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3).")
  )

  population <- list(
    species          = "human",
    n_subjects       = 24L,
    n_studies        = 1L,
    age_mean         = "57.5 +/- 19.9 years (mean +/- SD)",
    weight_mean      = "51.3 +/- 8.6 kg (mean +/- SD)",
    sex_female_pct   = 25.0,
    race_ethnicity   = NULL,
    disease_state    = "Critically ill adults receiving imipenem-cilastatin",
    dose_range       = paste(
      "500 mg imipenem intravenously every 6 or 8 h, or 1000 mg every 8 h",
      "(Zhang 2025 Supplementary Table S1). Zhang 2025 Results notes this",
      "as one of five studies reporting a maximum daily dose of 4000 mg,",
      "above the 1000-2000 mg/day the review gives as the usual",
      "recommendation. The infusion duration is not reported."
    ),
    regions          = "Vietnam",
    n_concentrations = 139L,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 14); 24 patients, 139",
      "samples, sex split 18 male / 6 female. Blood was sampled at 0,",
      "0.25, 0.67, 1.5 and 7 h after administration and assayed by HPLC-UV",
      "(Zhang 2025 Supplementary Table S1). Covariates were selected by",
      "the COSSAC method (conditional sampling use for stepwise approach",
      "based on correlation tests) -- the only study in the review to use",
      "it. In addition to the covariates in covariatesDataExcluded, ICU",
      "versus non-ICU unit, vasopressor use and mechanical ventilation",
      "were screened and not retained. Fitted in Monolix and evaluated by",
      "goodness-of-fit plots, NPDE, VPC and bootstrap (Zhang 2025 Table",
      "2).",
      "TYPESETTING ERROR IN THE SECONDARY SOURCE. Zhang 2025 Table 2",
      "prints this study's central volume as 'V1: 11.1 h'. The unit is a",
      "typographical slip -- a volume cannot be in hours, the sibling rows",
      "V2 and Q in the same cell are given in L and L/h respectively, and",
      "litres is the only unit dimensionally consistent with a",
      "two-compartment model. Encoded as 11.1 L. Note the coincidence that",
      "Q is also 11.1 in this model, which is presumably what invited the",
      "slip; both values are transcribed as printed.",
      "This model carries an unusually complete random-effects structure",
      "for the review -- etas on all four disposition parameters, including",
      "a 110% CV on intercompartmental clearance -- which for a 24-patient",
      "cohort suggests the peripheral parameters are weakly identified.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Dinh 2022 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Dinh et al. (2022) row.
    # NOTE that the clearance intercept is the value at CLcr = 0 because
    # the covariate model is uncentred; see covariateData$CRCL. =====
    lcl <- log(4.79); label("Clearance at CRCL = 0 mL/min (L/h)")           # Zhang 2025 Table 2 (Dinh 2022): CL = 4.79 L/h; leading coefficient of the uncentred Table 3 CL formula
    lvc <- log(11.1); label("Central volume of distribution V1 (L)")        # Zhang 2025 Table 2 (Dinh 2022): V1 = 11.1, printed with the erroneous unit 'h'; see population$notes
    lvp <- log(8.82); label("Peripheral volume of distribution V2 (L)")     # Zhang 2025 Table 2 (Dinh 2022): V2 = 8.82 L
    lq  <- log(11.1); label("Intercompartmental clearance Q (L/h)")         # Zhang 2025 Table 2 (Dinh 2022): Q = 11.1 L/h

    # ===== Covariate effect -- Zhang 2025 Table 3, Dinh et al. (2022):
    #   CL = 4.79 * e^(0.00642 * CLcr)
    # A log-linear (exponential) effect with NO centring term. The
    # coefficient has units of per mL/min, not a dimensionless exponent.
    e_crcl_cl <- 0.00642; label("Exponential coefficient of CRCL on CL (per mL/min)")  # Zhang 2025 Table 3 (Dinh 2022): e^(0.00642 * CLcr)

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
    # This is the only model in the review with an eta on every disposition
    # parameter. The 110% CV on Q gives omega^2 = log(1 + 1.10^2) = 0.7949;
    # note that under the alternative reading (omega SD = 1.10, so
    # omega^2 = 1.21) this term would differ by more than 50%, so Q is the
    # parameter in this file most sensitive to the scale convention.
    etalcl ~ log(1 + 0.387^2)  # Zhang 2025 Table 2 (Dinh 2022): IIV CL = 38.7%, read as an apparent CV
    etalvc ~ log(1 + 0.563^2)  # Zhang 2025 Table 2 (Dinh 2022): IIV V1 = 56.3%, read as an apparent CV
    etalvp ~ log(1 + 0.694^2)  # Zhang 2025 Table 2 (Dinh 2022): IIV V2 = 69.4%, read as an apparent CV
    etalq  ~ log(1 + 1.10^2)   # Zhang 2025 Table 2 (Dinh 2022): IIV Q = 110%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive. The additive term is printed in
    # mg/L, i.e. on the concentration scale, so it is taken as a standard
    # deviation, which is what nlmixr2's add() expects.
    propSd <- 0.221; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Dinh 2022): Proportional = 22.1%
    addSd  <- 0.445; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Dinh 2022): Additive = 0.445 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * CRCL)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
