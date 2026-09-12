Jaruratanasirikul_2021_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 50 critically",
    "ill Thai adults with life-threatening severe infections, treated in an",
    "intensive care unit with or without extracorporeal membrane",
    "oxygenation (Jaruratanasirikul 2021). Clearance is a linear function",
    "of CKD-EPI estimated glomerular filtration rate centred at 89",
    "mL/min/1.73 m^2, and the central volume is a linear function of",
    "adjusted body weight centred at 60 kg; ECMO support, ECMO type, ECMO",
    "flow rate and ECMO duration were all screened and none was retained.",
    "Inter-individual variability is exponential on clearance and on both",
    "volumes, and residual error is combined proportional plus additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Jaruratanasirikul 2021 when the",
    "primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Jaruratanasirikul S, Boonpeng A, Nawakitrangsan M, Samaeng M.",
    "NONMEM population pharmacokinetics and Monte Carlo dosing simulations",
    "of imipenem in critically ill patients with life-threatening severe",
    "infections during support with or without extracorporeal membrane",
    "oxygenation in an intensive care unit.",
    "Pharmacotherapy. 2021;41(7):572-597. doi:10.1002/phar.2597.",
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
    central      = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1  = list(analyte = "imipenem", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Glomerular filtration rate estimated by the CKD-EPI equation.",
        "Zhang 2025's abbreviation list distinguishes 'eGFR CKD-EPI' from",
        "'CKD-EPI-abs, absolute CKD-EPI (i.e., CKD-EPI, multiplied by",
        "BSA)', and this study's formula uses the former, so the value is",
        "the BSA-NORMALISED rate a CKD-EPI calculator returns by default,",
        "in mL/min/1.73 m^2. The reference of 89 is consistent with that",
        "reading: 89 mL/min/1.73 m^2 is the conventional lower bound of",
        "normal renal function, whereas 89 mL/min absolute would be an odd",
        "centring constant."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 89 mL/min/1.73 m^2. Enters CL as the LINEAR-DEVIATION",
        "term CL = 13.3 + 0.112 * (eGFR_CKD-EPI - 89) L/h (Zhang 2025",
        "Table 3), NOT as a power of a ratio -- the slope 0.112 therefore",
        "carries units of L/h per mL/min/1.73 m^2. The linear form is",
        "well-behaved over the whole physiological domain here because the",
        "slope is positive: CL would only reach zero at a negative eGFR",
        "(89 - 13.3/0.112 = -30 mL/min/1.73 m^2), so no upper or lower",
        "domain guard is needed on this term. Zhang 2025 Table 3 records",
        "that four renal-function estimators were screened head to head",
        "for this study -- Cockcroft-Gault, MDRD-4, MDRD-6 and CKD-EPI --",
        "and that CKD-EPI won. Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts CKD-EPI eGFR",
        "-- precedent: Bajaj 2017 nivolumab, de Velde 2020 imipenem",
        "(which uses the ABSOLUTE CKD-EPI variant instead, so the two",
        "models' CRCL columns are on different scales and their reference",
        "values are not comparable)."
      ),
      source_name        = "eGFR CKD-EPI"
    ),
    ABW = list(
      description        = paste(
        "Adjusted body weight, the obesity-dosing size descriptor",
        "conventionally computed as ABW = IBW + 0.4 * (TBW - IBW).",
        "Zhang 2025's Table 3 footnote glosses the abbreviation as",
        "'ABW, adjusted body weight' and the same table lists 'actual BW'",
        "and 'ideal BW' as separately screened covariates, so ABW here is",
        "the interpolated dosing weight and not an alias of either."
      ),
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 60 kg. Enters the central volume as the",
        "LINEAR-DEVIATION term V1 = 13.8 - 0.348 * (ABW - 60) L (Zhang",
        "2025 Table 3), so the slope -0.348 carries units of L per kg.",
        "DOMAIN CAVEAT. The slope is NEGATIVE and the form is linear, so",
        "V1 falls to zero at ABW = 60 + 13.8/0.348 = 99.7 kg and is",
        "negative above it. The model must not be simulated above that",
        "adjusted body weight; the study's own cohort is light (median",
        "total body weight 62.9 kg, IQR 52.8-70.0, Zhang 2025 Table 1),",
        "so the term was fitted far from the zero-crossing and the",
        "extrapolation is unsupported as well as unphysical. The vignette",
        "restricts its virtual cohort accordingly.",
        "FORMULA NOT REPORTED. Neither the review nor its Table 3",
        "footnote states which ideal-body-weight formula or which",
        "interpolation factor the primary used; the 0.4 factor named in",
        "the description is the usual convention, not a value read from",
        "this source. Re-verify against Jaruratanasirikul 2021 when the",
        "primary is obtained. Founding example for the ABW canonical",
        "column registered in inst/references/covariate-columns.md with",
        "this extraction."
      ),
      source_name        = "ABW"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed. Only the
  # screened covariates that have a canonical column are listed here; the
  # remainder are recorded in population$notes.
  covariatesDataExcluded <- list(
    AGE         = list(description = "Age",                              units = "years",    type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort median 56.2 years, IQR 40.95-66.6 (Table 1)."),
    SEXF        = list(description = "Female sex",                       units = "(binary)", type = "binary",     notes = "Screened as 'gender', not retained (Zhang 2025 Table 3). Cohort is 35 male / 15 female (Table 1)."),
    WT          = list(description = "Total (actual) body weight",       units = "kg",       type = "continuous", notes = "Screened as 'actual BW', not retained (Zhang 2025 Table 3); the adjusted body weight ABW won the body-size comparison and is the retained descriptor on V1. Cohort median 62.9 kg, IQR 52.8-70.0 (Table 1)."),
    IBW         = list(description = "Ideal body weight",                units = "kg",       type = "continuous", notes = "Screened as 'ideal BW', not retained (Zhang 2025 Table 3). Lost to ABW on V1."),
    BMI         = list(description = "Body mass index",                  units = "kg/m^2",   type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    ECMO_STATUS = list(description = "ECMO treatment status",            units = "(binary)", type = "binary",     notes = "Screened as 'the use of ECMO support', not retained (Zhang 2025 Table 3). This is the substantive negative result of the study: the cohort is explicitly assembled with and without ECMO, and neither ECMO status nor ECMO type, flow rate or duration reached significance on any parameter."),
    Q_ECMO      = list(description = "ECMO circuit blood flow rate",     units = "L/min",    type = "continuous", notes = "Screened as 'ECMO flow rate', not retained (Zhang 2025 Table 3). Applies only to the ECMO-supported subset."),
    T_ECMO      = list(description = "Time since ECMO cannulation",      units = "h",        type = "continuous", notes = "Screened as 'duration of ECMO', not retained (Zhang 2025 Table 3). Applies only to the ECMO-supported subset."),
    APACHE_II   = list(description = "APACHE II score at ICU admission", units = "(score)",  type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    MECH_VENT   = list(description = "Invasive mechanical ventilation",  units = "(binary)", type = "binary",     notes = "Screened as 'mechanical ventilation support', not retained (Zhang 2025 Table 3)."),
    ALB         = list(description = "Serum albumin",                    units = "g/L",      type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Contrast Por 2021, the one study in this review that does retain albumin, on both volumes.")
  )

  population <- list(
    species          = "human",
    n_subjects       = 50L,
    n_studies        = 1L,
    age_median       = "56.2 years (IQR 40.95-66.6)",
    weight_median    = "62.9 kg (IQR 52.8-70.0)",
    sex_female_pct   = 30.0,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Critically ill adults in an intensive care unit with",
      "life-threatening severe infections, supported with or without",
      "extracorporeal membrane oxygenation (ECMO). This is one of two",
      "studies in the Zhang 2025 review that enrol ECMO-supported",
      "patients; the other is Chen 2020."
    ),
    dose_range       = paste(
      "250-500 mg imipenem intravenously every 6-12 h, 500 mg every 6 h,",
      "and 1000 mg every 6 or 8 h (Zhang 2025 Supplementary Table S1).",
      "The infusion duration is not reported by the review."
    ),
    regions          = "Thailand",
    n_concentrations = 534L,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 11); 50 patients, 534",
      "samples, sex split 35 male / 15 female. Blood was sampled at 0,",
      "0-0.5, 0.5-2, 2-4 and 4-12 h after administration and assayed by",
      "HPLC-UV (Zhang 2025 Supplementary Table S1). Covariates were",
      "selected by forward inclusion (dOFV > 3.84, p < 0.05) and backward",
      "elimination (dOFV > 6.64, p < 0.01) -- the most permissive",
      "thresholds used by any study in the review, which is worth bearing",
      "in mind when reading the two retained effects. Fitted in NONMEM and",
      "evaluated by pcVPC, NPDE and bootstrap (Zhang 2025 Table 2).",
      "SCREENED COVARIATES WITHOUT A CANONICAL COLUMN. In addition to the",
      "entries in covariatesDataExcluded, Zhang 2025 Table 3 records that",
      "ECMO type, SOFA score, acute kidney injury, fluid balance, use of",
      "inotropes, septic shock and mean arterial blood pressure were",
      "screened and not retained. These have no canonical column in",
      "inst/references/covariate-columns.md and are noted here rather than",
      "registered, since no coefficient is reported for any of them.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Jaruratanasirikul 2021 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Jaruratanasirikul et al.
    # (2021) row, with the two covariate-carrying typical values taken from
    # the Table 3 formula intercepts (see the V1 note below). =====
    lcl <- log(13.3); label("Clearance at CRCL = 89 mL/min/1.73 m^2 (L/h)")  # Zhang 2025 Table 2 (Jaruratanasirikul 2021): CL = 13.3 L/h; identical to the Table 3 CL formula intercept
    # INTERNAL INCONSISTENCY IN THE SOURCE. Zhang 2025 Table 2 prints
    # V1 = 13.6 L, but the Table 3 covariate formula reads
    # V1 = 13.8 - 0.348 * (ABW - 60), i.e. an intercept of 13.8 L at the
    # reference ABW of 60 kg. The formula intercept is used here because it
    # is the quantity the covariate term is conditioned on -- pairing the
    # Table 2 value with the Table 3 slope would shift every predicted
    # volume by a constant 0.2 L and break the reference-subject identity
    # V1(ABW = 60) = intercept. The 1.5% discrepancy is recorded in the
    # vignette's Errata. Note that CL has no such conflict: Table 2 and the
    # Table 3 intercept agree exactly at 13.3 L/h.
    lvc <- log(13.8); label("Central volume at ABW = 60 kg (L)")              # Zhang 2025 Table 3 (Jaruratanasirikul 2021) formula intercept: V1 = 13.8 L at ABW = 60 kg (Table 2 prints 13.6 L; see comment)
    lvp <- log(16.9); label("Peripheral volume of distribution (L)")          # Zhang 2025 Table 2 (Jaruratanasirikul 2021): V2 = 16.9 L
    lq  <- log(24.3); label("Intercompartmental clearance (L/h)")             # Zhang 2025 Table 2 (Jaruratanasirikul 2021): Q = 24.3 L/h

    # ===== Covariate effects -- Zhang 2025 Table 3, Jaruratanasirikul et
    # al. (2021). Both are LINEAR-DEVIATION forms, not power-of-ratio
    # forms, so the coefficients carry units:
    #   CL = 13.3 + 0.112 * (eGFR_CKD-EPI - 89)   L/h
    #   V1 = 13.8 - 0.348 * (ABW - 60)            L
    # =====
    e_crcl_cl <- 0.112;  label("Linear slope of CL on (CRCL - 89) (L/h per mL/min/1.73 m^2)")  # Zhang 2025 Table 3 (Jaruratanasirikul 2021): CL = 13.3 + 0.112 x (eGFRCKD-EPI - 89)
    e_abw_vc  <- -0.348; label("Linear slope of Vc on (ABW - 60) (L per kg)")                  # Zhang 2025 Table 3 (Jaruratanasirikul 2021): V1 = 13.8 - 0.348 x (ABW - 60)

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 prints IIV as a bare percentage
    # per parameter without stating the convention, and the column mixes at
    # least three conventions across the review's constituent studies (the
    # audit against the four already-extracted primaries is in the
    # vignette's 'Assumptions and deviations' section). Every model
    # transcribed from this review uses the same documented reading: the
    # printed percentage is an apparent CV of a log-normal random effect,
    # so omega^2 = log(1 + CV^2). The three IIVs here are large (51-67%),
    # which is what makes the reading matter most for this study: at 66.9%
    # the CV and log-SD readings differ by about 12% in omega.
    etalcl ~ log(1 + 0.510^2)  # Zhang 2025 Table 2 (Jaruratanasirikul 2021): IIV CL = 51%, read as an apparent CV
    etalvc ~ log(1 + 0.669^2)  # Zhang 2025 Table 2 (Jaruratanasirikul 2021): IIV V1 = 66.9%, read as an apparent CV
    etalvp ~ log(1 + 0.560^2)  # Zhang 2025 Table 2 (Jaruratanasirikul 2021): IIV V2 = 56%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive; this is one of the few studies in
    # the review reporting both arms.
    propSd <- 0.183; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Jaruratanasirikul 2021): Proportional = 18.3%
    addSd  <- 0.216; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Jaruratanasirikul 2021): Additive = 0.216 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # Both covariate terms are linear deviations from a reference subject,
    # so the typical value is built as (intercept + slope * (cov - ref)) and
    # the exponential eta multiplies the whole covariate-adjusted typical
    # value -- the standard NONMEM idiom for an additive covariate model
    # with a log-normal ETA (same encoding as Lamoth_2009_imipenem.R).
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 89)) * exp(etalcl)
    vc <- (exp(lvc) + e_abw_vc  * (ABW  - 60)) * exp(etalvc)
    vp <- exp(lvp + etalvp)
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
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
