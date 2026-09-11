Truong_2025_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 151 Dutch",
    "critically ill and non-critically ill adults (Truong 2025). Clearance",
    "scales as a power of Cockcroft-Gault creatinine clearance, the only",
    "covariate retained; volumes and intercompartmental clearance carry",
    "none. Inter-individual variability is exponential on clearance only,",
    "and residual error is proportional.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Truong 2025 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Truong AQ, Smeets TJL, Terrier J, Li L, Dao XC, Strojil J, et al.",
    "Inadequate imipenem dosing in patients with decreased kidney function:",
    "a global clinical pharmacokinetic study.",
    "Clin Microbiol Infect. 2025;31(9):1518-1525.",
    "doi:10.1016/j.cmi.2025.05.005.",
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
        "Reference 87.6 mL/min. Enters CL as the power term",
        "(CLcr_CG/87.6)^0.462 (Zhang 2025 Table 3). The review does not",
        "report this cohort's renal function distribution, so whether",
        "87.6 is the cohort median cannot be checked from the secondary",
        "source. Zhang 2025 Table 3 records that Cockcroft-Gault was",
        "chosen over CKD-EPI and MDRD-4, which were screened head to head",
        "alongside five body-size descriptors and ICU versus non-ICU unit.",
        "Zhang 2025 Results names this study as one of five whose final",
        "model retains creatinine clearance as its SOLE covariate. Stored",
        "under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts raw",
        "Cockcroft-Gault mL/min when the source does not BSA-normalise --",
        "precedent: Bai 2024 imipenem, Wang 2024 imipenem."
      ),
      source_name        = "CLcr CG"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    AGE  = list(description = "Age",                 units = "years",    type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort median 63 years, IQR 51-72 (Table 1)."),
    SEXF = list(description = "Female sex",          units = "(binary)", type = "binary",     notes = "Screened, not retained (Zhang 2025 Table 3). Zhang 2025 Table 1 records the sex split as 151/0, i.e. an all-male modelling cohort; if that is correct the covariate is unidentifiable by construction rather than merely non-significant, and the figure should be re-checked against the primary."),
    WT   = list(description = "Total body weight",   units = "kg",       type = "continuous", notes = "Screened as TBW, not retained (Zhang 2025 Table 3). Cohort median 70 kg, IQR 61.2-82 (Table 1). Adjusted and ideal body weight were screened separately and also not retained."),
    HT   = list(description = "Height",              units = "cm",       type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    BMI  = list(description = "Body mass index",     units = "kg/m^2",   type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    BSA  = list(description = "Body surface area",   units = "m^2",      type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3).")
  )

  population <- list(
    species          = "human",
    n_subjects       = 151L,
    n_studies        = 1L,
    age_median       = "63 years (IQR 51-72)",
    weight_median    = "70 kg (IQR 61.2-82)",
    sex_female_pct   = 0.0,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Critically ill AND non-critically ill adults receiving",
      "imipenem-cilastatin. This is the only study among the 18 in the",
      "Zhang 2025 review to pool both, and the second largest cohort after",
      "Chen 2020."
    ),
    dose_range       = paste(
      "Not reported. Zhang 2025 Supplementary Table S1 records 'NR' for",
      "this study's dosage, and the review's Results text calls it out",
      "explicitly: 'the imipenem dosage was not specified in Truong et al.",
      "(2025) study.' Simulations with this model must therefore assume a",
      "regimen rather than reproduce the study's own; the Dutch sibling",
      "study de Velde 2020 used 500 mg every 6 h."
    ),
    regions          = "Netherlands",
    n_concentrations = 322L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 17); 151 patients,",
      "322 samples. Zhang 2025 Table 1 records the sex split as 151 male /",
      "0 female, which would make this an all-male cohort; the review does",
      "not comment on it and it is unusual enough to warrant re-checking",
      "against the primary. Blood was sampled at peak, intermediate and",
      "trough time points and assayed by HPLC-UV (Zhang 2025 Supplementary",
      "Table S1). Covariates were selected by forward inclusion (p < 0.01)",
      "and backward elimination (p < 0.001). Fitted in NONMEM and",
      "evaluated by bootstrap, VPC, goodness-of-fit plots and, unusually",
      "for this review, external validation; only three of the 18 studies",
      "performed external validation (Zhang 2025 Table 2, Results).",
      "Note that this model's intercompartmental clearance of 2.9 L/h is",
      "the second lowest among the review's two-compartment models, so the",
      "distribution phase is slow relative to its siblings even though the",
      "peripheral volume (21.4 L) is substantial.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Truong 2025 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Truong et al. (2025) row.
    # Typical values at the reference Cockcroft-Gault CLcr of 87.6 mL/min
    # (Zhang 2025 Table 3 formulation). =====
    lcl <- log(14.6); label("Clearance at CRCL = 87.6 mL/min (L/h)")  # Zhang 2025 Table 2 (Truong 2025): CL = 14.6 L/h; leading coefficient of the Table 3 CL formula
    lvc <- log(28.7); label("Central volume of distribution V1 (L)")  # Zhang 2025 Table 2 (Truong 2025): V1 = 28.7 L
    lvp <- log(21.4); label("Peripheral volume of distribution V2 (L)")  # Zhang 2025 Table 2 (Truong 2025): V2 = 21.4 L
    lq  <- log(2.9);  label("Intercompartmental clearance Q (L/h)")   # Zhang 2025 Table 2 (Truong 2025): Q = 2.9 L/h

    # ===== Covariate effect -- Zhang 2025 Table 3, Truong et al. (2025):
    #   CL = 14.6 * (CLcr_CG/87.6)^0.462
    # The review prints no exponential eta term in this formula, but Table
    # 2 reports an IIV on CL of 35.9%, so the random effect is applied
    # multiplicatively in the usual way.
    e_crcl_cl <- 0.462; label("Power exponent on (CRCL/87.6) for CL (unitless)")  # Zhang 2025 Table 3 (Truong 2025): (CLcr CG/87.6)^0.462

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
    # IIV is reported on CL only; the review gives no omega for V1, V2
    # or Q.
    etalcl ~ log(1 + 0.359^2)  # Zhang 2025 Table 2 (Truong 2025): IIV CL = 35.9%, read as an apparent CV

    # ===== Residual error =====
    propSd <- 0.383; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Truong 2025): Proportional = 38.30%
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (CRCL / 87.6)^e_crcl_cl
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

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
    Cc ~ prop(propSd)
  })
}
