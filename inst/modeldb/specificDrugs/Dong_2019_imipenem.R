Dong_2019_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 56 Chinese",
    "children with haematological malignancies (Dong 2019). Clearance",
    "carries three multiplicative power covariates -- allometric body",
    "weight, age, and Schwartz-estimated creatinine clearance -- while both",
    "volumes scale linearly with body weight and intercompartmental",
    "clearance scales allometrically. Inter-individual variability is",
    "exponential on clearance and the central volume, and residual error is",
    "combined proportional plus additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Dong 2019 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Dong L, Zhai XY, Yang YL, Wang L, Zhou Y, Shi HY, et al.",
    "Population pharmacokinetics and dosing optimization of imipenem in",
    "children with hematological malignancies.",
    "Antimicrob Agents Chemother. 2019;63(6):e00006-19.",
    "doi:10.1128/AAC.00006-19.",
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
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 18 kg, close to the cohort mean of 18.65 +/- 6.90 kg",
        "(Zhang 2025 Table 1); the review does not say how the reference",
        "was chosen. Enters CL and Q with the theory-based allometric",
        "exponent 0.75 and both volumes with an exponent of 1, per the",
        "formulations printed in Zhang 2025 Table 3:",
        "CL = 8.6 * (BW/18)^0.75 * (age/4.69)^0.265 * (CLcr/214)^0.509,",
        "V1 = 7.2 * (BW/18), V2 = 6.51 * (BW/18), Q = 0.996 * (BW/18)^0.75."
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 4.69 years, close to but not equal to the cohort mean of",
        "4.86 +/- 2.33 years (Zhang 2025 Table 1); the review prints 4.69",
        "inside the Table 3 formula and does not reconcile it with Table 1,",
        "so the equation constant is used. Enters CL as the power term",
        "(age/4.69)^0.265. Zhang 2025 Results notes that this is the only",
        "one of the review's 18 models to retain age as a covariate in an",
        "adult-or-paediatric final model besides the neonatal Dao 2022, and",
        "its Discussion reads age here as a surrogate for renal and",
        "body-composition maturation over and above weight."
      ),
      source_name        = "Age"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated by the Schwartz equation, which",
        "returns a BSA-normalised value in mL/min/1.73 m^2"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Zhang 2025 Table 3 names the covariate 'CLcr Schwartz' -- the",
        "paediatric height-and-creatinine bedside equation, which is",
        "BSA-normalised. Reference 214 mL/min/1.73 m^2, a markedly",
        "supranormal value; the review does not report the cohort renal",
        "function distribution, so whether 214 is the cohort median cannot",
        "be checked from the secondary source, and the reference should be",
        "re-verified against the primary. Enters CL as the power term",
        "(CLcr/214)^0.509. Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts",
        "Schwartz-estimated paediatric CrCL -- precedent: Jung 2024",
        "vancomycin, MedellinGaribay 2015 gentamicin."
      ),
      source_name        = "CLcr"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 56L,
    n_studies        = 1L,
    age_mean         = "4.86 +/- 2.33 years (mean +/- SD)",
    weight_mean      = "18.65 +/- 6.90 kg (mean +/- SD)",
    sex_female_pct   = 46.4,
    race_ethnicity   = NULL,
    disease_state    = "Children with haematological malignancies receiving imipenem-cilastatin",
    dose_range       = paste(
      "15-25 mg/kg imipenem intravenously every 6 h (Zhang 2025",
      "Supplementary Table S1). The infusion duration is not reported by",
      "the review, but the sampling schedule implies an infusion: samples",
      "were taken 3-5 min and 0.5-1 h after the END of infusion, and",
      "2-6 h after the START of infusion."
    ),
    regions          = "China",
    n_concentrations = 136L,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 6); 56 patients, 136",
      "plasma samples, sex split 30 male / 26 female. Imipenem assayed by",
      "HPLC-UV (Zhang 2025 Supplementary Table S1). Covariates were",
      "selected by forward inclusion (dOFV > 3.84, p < 0.05) and backward",
      "elimination (dOFV > 6.635, p < 0.01); age, body weight and Schwartz",
      "creatinine clearance were the only covariates screened and all three",
      "were retained on CL, with body weight additionally retained on V1,",
      "V2 and Q (Zhang 2025 Table 3). Fitted in NONMEM and evaluated by",
      "bootstrap, goodness-of-fit plots and NPDE (Zhang 2025 Table 2).",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Dong 2019 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Dong et al. (2019) row.
    # Typical values for the reference child: 18 kg, 4.69 years, Schwartz
    # CLcr 214 mL/min/1.73 m^2 (Zhang 2025 Table 3 formulations). =====
    lcl <- log(8.6);   label("Clearance at the reference child (L/h)")            # Zhang 2025 Table 2 (Dong 2019): CL = 8.6 L/h; leading coefficient of the Table 3 CL formula
    lvc <- log(7.2);   label("Central volume of distribution at 18 kg (L)")       # Zhang 2025 Table 2 (Dong 2019): V1 = 7.2 L; leading coefficient of the Table 3 V1 formula
    lvp <- log(6.51);  label("Peripheral volume of distribution at 18 kg (L)")    # Zhang 2025 Table 2 (Dong 2019): V2 = 6.51 L; leading coefficient of the Table 3 V2 formula
    lq  <- log(0.996); label("Intercompartmental clearance at 18 kg (L/h)")       # Zhang 2025 Table 2 (Dong 2019): Q = 0.996 L/h; leading coefficient of the Table 3 Q formula

    # ===== Covariate effects -- Zhang 2025 Table 3, Dong et al. (2019) row
    #   CL = 8.6   * (BW/18)^0.75 * (age/4.69)^0.265 * (CLcr/214)^0.509
    #   V1 = 7.2   * (BW/18)
    #   V2 = 6.51  * (BW/18)
    #   Q  = 0.996 * (BW/18)^0.75
    # The review prints the two weight exponents as the theory-based
    # allometric 0.75 on the flow terms and an implicit 1 on the volumes.
    # It does not report standard errors for any of them, and does not say
    # whether 0.75 was estimated or fixed. Both weight exponents are
    # therefore left UNWRAPPED (not fixed()) rather than guessing: 0.75 and
    # 1 are the canonical allometric values and may well have been fixed,
    # but the secondary source gives no evidence either way, and encoding a
    # wrong fixed() flag would misreport provenance. This is recorded in
    # the vignette Errata.
    e_wt_cl   <- 0.75;  label("Allometric WT exponent on CL and Q (unitless)")            # Zhang 2025 Table 3 (Dong 2019): (BW/18)^0.75 in the CL and Q formulas
    e_wt_vc   <- 1;     label("Linear WT exponent on Vc and Vp (unitless)")               # Zhang 2025 Table 3 (Dong 2019): (BW/18) with no printed exponent in the V1 and V2 formulas
    e_age_cl  <- 0.265; label("Power exponent on (AGE/4.69) for CL (unitless)")           # Zhang 2025 Table 3 (Dong 2019): (age/4.69)^0.265
    e_crcl_cl <- 0.509; label("Power exponent on (CRCL/214) for CL (unitless)")           # Zhang 2025 Table 3 (Dong 2019): (CLcr/214)^0.509

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
    etalcl ~ log(1 + 0.188^2)  # Zhang 2025 Table 2 (Dong 2019): IIV CL = 18.8%, read as an apparent CV
    etalvc ~ log(1 + 0.092^2)  # Zhang 2025 Table 2 (Dong 2019): IIV V1 = 9.2%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive. The additive term is printed in
    # mg/L, i.e. on the concentration scale, so it is taken as a standard
    # deviation, which is what nlmixr2's add() expects.
    propSd <- 0.395; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Dong 2019): Proportional = 39.5%
    addSd  <- 0.205; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Dong 2019): Additive = 0.205 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # Zhang 2025 Table 3 places the exponential eta on CL and V1; Q and V2
    # carry their covariate terms but no random effect.
    cl <- exp(lcl + etalcl) * (WT / 18)^e_wt_cl * (AGE / 4.69)^e_age_cl * (CRCL / 214)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 18)^e_wt_vc
    vp <- exp(lvp)          * (WT / 18)^e_wt_vc
    q  <- exp(lq)           * (WT / 18)^e_wt_cl

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
