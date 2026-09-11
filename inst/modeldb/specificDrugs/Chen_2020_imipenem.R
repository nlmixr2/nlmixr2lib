Chen_2020_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 247 Chinese",
    "critically ill adults with and without extracorporeal membrane",
    "oxygenation (Chen 2020). Clearance carries three multiplicative",
    "covariates -- Cockcroft-Gault creatinine clearance and body weight as",
    "power terms, and ECMO support as an exponential indicator that raises",
    "clearance about three-fold. Volumes and intercompartmental clearance",
    "carry no covariates. Inter-individual variability is exponential on",
    "clearance and the central volume, and residual error is combined",
    "proportional plus additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Chen 2020 when the primary is",
    "obtained -- see the population notes for two values that look like",
    "transcription defects in the secondary source.",
    sep = " "
  )
  reference <- paste(
    "Chen W, Zhang D, Lian W, Wang X, Du W, Zhang Z, et al.",
    "Imipenem population pharmacokinetics: therapeutic drug monitoring data",
    "collected in critically ill patients with or without extracorporeal",
    "membrane oxygenation.",
    "Antimicrob Agents Chemother. 2020;64(6):e00385-20.",
    "doi:10.1128/AAC.00385-20.",
    "Parameters transcribed from Zhang P, Zhao Y, Zhu J, Yang Y, Liang G,",
    "Wang X, Yu Z. Population pharmacokinetics of imipenem in different",
    "populations for individualized dosing: a systematic review.",
    "Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055",
    "(Tables 1-3, Supplementary Table S1; the review cites this study as",
    "'Chen et al. (2020b)').",
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
        "The review does not state whether the value is BSA-normalised;",
        "Cockcroft-Gault natively returns raw mL/min and the review's",
        "abbreviation list glosses 'CLcrCG' as 'clearance creatinine",
        "estimated by Cockcroft and Gault equation' with no normalisation",
        "mentioned, so raw mL/min is used here."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 59.1 mL/min. Enters CL as the power term",
        "(CLcr/59.1)^0.295 (Zhang 2025 Table 3). The review does not report",
        "this cohort's renal function distribution, so whether 59.1 is the",
        "cohort median cannot be checked from the secondary source. Stored",
        "under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts raw",
        "Cockcroft-Gault mL/min when the source does not BSA-normalise --",
        "precedent: Bai 2024 imipenem, Wang 2024 imipenem, Couffignal 2014",
        "imipenem."
      ),
      source_name        = "CLcr"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 65.0 kg, equal to the cohort median of 65 kg (range",
        "37.5-110.0; Zhang 2025 Table 1). Enters CL as the power term",
        "(BW/65.0)^0.306. Note that the exponent is 0.306, well below the",
        "theory-based allometric 0.75, and that no weight term is carried",
        "on either volume -- both are departures from the usual allometric",
        "parameterisation and are reproduced as printed."
      ),
      source_name        = "BW"
    ),
    ECMO_STATUS = list(
      description        = "Extracorporeal membrane oxygenation support during imipenem therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no ECMO support",
      notes              = paste(
        "Zhang 2025 Table 3 prints the covariate model as two branches:",
        "with ECMO, CL = 8.88 * (CLcr/59.1)^0.295 * (BW/65.0)^0.306 *",
        "e^1.16 * e^eta_CL; without ECMO, the same expression without the",
        "e^1.16 factor. That is exactly an exponential indicator",
        "exp(1.16 * ECMO_STATUS), which is how it is encoded here.",
        "MAGNITUDE WARNING: e^1.16 = 3.19, i.e. the model says ECMO support",
        "more than triples imipenem clearance. That is a very large effect",
        "for a mechanism that Zhang 2025's own Discussion attributes mainly",
        "to circuit adsorption INCREASING the apparent volume of",
        "distribution rather than clearance. The review reports no standard",
        "error or confidence interval for the coefficient, so the magnitude",
        "cannot be checked from the secondary source and should be",
        "re-verified against the primary before the model is used to",
        "support ECMO dosing. Recorded in the vignette Errata."
      ),
      source_name        = "ECMO"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    AGE  = list(description = "Age",                    units = "years",        type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort median 67 years, range 20-97 (Table 1)."),
    SEXF = list(description = "Female sex",             units = "(binary)",     type = "binary",     notes = "Screened, not retained (Zhang 2025 Table 3). Cohort 80/247 female (Table 1). Zhang 2025 Results notes that sex was screened in 10 of the review's studies and retained in none."),
    HT   = list(description = "Height",                 units = "cm",           type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    BMI  = list(description = "Body mass index",        units = "kg/m^2",       type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    CREAT = list(description = "Serum creatinine",      units = "umol/L",       type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3); the derived CLcr was retained instead. Units not stated by the review."),
    ALT  = list(description = "Alanine transaminase",   units = "U/L",          type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    AST  = list(description = "Aspartate transaminase", units = "U/L",          type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    ALB  = list(description = "Serum albumin",          units = "g/L",          type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    TBIL = list(description = "Total bilirubin",        units = "umol/L",       type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    HGB  = list(description = "Haemoglobin",            units = "g/L",          type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    PLT  = list(description = "Platelet count",         units = "10^9 cells/L", type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    RRT_CRRT_STATUS = list(description = "Continuous renal replacement therapy", units = "(binary)", type = "binary", notes = "Screened, not retained (Zhang 2025 Table 3). ECMO was retained but CRRT was not.")
  )

  population <- list(
    species          = "human",
    n_subjects       = 247L,
    n_studies        = 1L,
    age_median       = "67 years (range 20-97)",
    weight_median    = "65 kg (range 37.5-110.0)",
    sex_female_pct   = 32.4,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Critically ill adults receiving imipenem-cilastatin, with and",
      "without extracorporeal membrane oxygenation support. This is the",
      "largest cohort among the 18 studies in the Zhang 2025 review."
    ),
    dose_range       = paste(
      "250 mg every 12 h; 500 mg every 6, 8 or 12 h; 500 mg in the morning",
      "with 250 mg at night; and 1000 mg every 6, 8 or 12 h (Zhang 2025",
      "Supplementary Table S1). The infusion duration is not reported by",
      "the review."
    ),
    regions          = "China",
    n_concentrations = 580L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 8); 247 patients, 580",
      "samples, sex split 167 male / 80 female. Two sampling schemes were",
      "pooled: routine therapeutic drug monitoring at 3 h and 0.5 h before",
      "the next dose, and a richer ECMO schedule sampling before the dose",
      "and at 0.5, 1, 2, 3, 6 and 8 h after the start of infusion, taken",
      "after the fourth dose both during ECMO and after its withdrawal.",
      "Imipenem assayed by LC-MS/MS -- one of only three studies in the",
      "review not using HPLC-UV (Zhang 2025 Supplementary Table S1).",
      "Covariates were selected by forward inclusion (dOFV > 3.84,",
      "p < 0.05) and backward elimination (dOFV > 7.88, p < 0.005).",
      "Fitted in NONMEM and evaluated by bootstrap, VPC and goodness-of-fit",
      "plots (Zhang 2025 Table 2).",
      "SUSPECTED TRANSCRIPTION DEFECT IN THE SECONDARY SOURCE. The residual",
      "error reproduced by Zhang 2025 Table 2 for this study -- 6.2%",
      "proportional plus 0.003 mg/L additive -- is implausibly small for a",
      "sparse, routine-TDM dataset of 580 samples across 247 critically ill",
      "patients, and is by a wide margin the smallest residual error among",
      "the review's 18 studies (the next smallest proportional term is",
      "15.9%, and the next smallest additive term is 0.04 mg/L). Taken at",
      "face value the model predicts observed concentrations to within a",
      "few percent, which no imipenem TDM dataset supports. The values are",
      "transcribed verbatim because inventing a substitute would be worse,",
      "but a user simulating residual-error-bearing observations from this",
      "model should treat the noise level as unverified. See also the ECMO",
      "coefficient warning in covariateData$ECMO_STATUS. Both are recorded",
      "in the vignette Errata.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Chen 2020 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Chen et al. (2020b) row.
    # Typical values for the reference subject: CLcr 59.1 mL/min, 65.0 kg,
    # no ECMO (Zhang 2025 Table 3 formulations). =====
    lcl <- log(8.88); label("Clearance at the reference subject, no ECMO (L/h)")  # Zhang 2025 Table 2 (Chen 2020b): CL = 8.88 L/h; leading coefficient of the Table 3 CL formula
    lvc <- log(20.5); label("Central volume of distribution V1 (L)")              # Zhang 2025 Table 2 (Chen 2020b): V1 = 20.5 L
    lvp <- log(8.86); label("Peripheral volume of distribution V2 (L)")           # Zhang 2025 Table 2 (Chen 2020b): V2 = 8.86 L
    lq  <- log(1.74); label("Intercompartmental clearance Q (L/h)")               # Zhang 2025 Table 2 (Chen 2020b): Q = 1.74 L/h

    # ===== Covariate effects -- Zhang 2025 Table 3, Chen et al. (2020b) row
    #   With ECMO:    CL = 8.88 * (CLcr/59.1)^0.295 * (BW/65.0)^0.306 * e^1.16 * e^eta_CL
    #   Without ECMO: CL = 8.88 * (CLcr/59.1)^0.295 * (BW/65.0)^0.306 *          e^eta_CL
    # which is exp(e_ecmo_cl * ECMO_STATUS) with e_ecmo_cl = 1.16.
    e_crcl_cl <- 0.295; label("Power exponent on (CRCL/59.1) for CL (unitless)")       # Zhang 2025 Table 3 (Chen 2020b): (CLcr/59.1)^0.295
    e_wt_cl   <- 0.306; label("Power exponent on (WT/65.0) for CL (unitless)")         # Zhang 2025 Table 3 (Chen 2020b): (BW/65.0)^0.306
    e_ecmo_cl <- 1.16;  label("Exponential coefficient of ECMO support on CL (unitless)")  # Zhang 2025 Table 3 (Chen 2020b): the e^1.16 factor present only in the 'With ECMO' branch

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
    # IIV is reported on CL and V1 only; no omega is given for V2 or Q.
    etalcl ~ log(1 + 0.177^2)  # Zhang 2025 Table 2 (Chen 2020b): IIV CL = 17.7%, read as an apparent CV
    etalvc ~ log(1 + 0.148^2)  # Zhang 2025 Table 2 (Chen 2020b): IIV V1 = 14.8%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive, both transcribed verbatim from
    # the review. See the SUSPECTED TRANSCRIPTION DEFECT note in
    # population$notes -- these are implausibly small and are the least
    # trustworthy numbers in this file.
    propSd <- 0.062; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Chen 2020b): Proportional = 6.2%
    addSd  <- 0.003; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Chen 2020b): Additive = 0.003 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (CRCL / 59.1)^e_crcl_cl * (WT / 65.0)^e_wt_cl *
          exp(e_ecmo_cl * ECMO_STATUS)
    vc <- exp(lvc + etalvc)
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
    Cc ~ prop(propSd) + add(addSd)
  })
}
