Por_2021_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 23 US burn",
    "patients with and without continuous venovenous haemofiltration",
    "(Por 2021). Clearance is described by two separate branches: patients",
    "not on CVVH carry power effects of Cockcroft-Gault creatinine",
    "clearance and body weight, while patients on CVVH carry an allometric",
    "body-weight effect plus an additive filter clearance. Both volumes",
    "carry a strong inverse power effect of serum albumin, and the central",
    "volume additionally scales with body weight. Inter-individual",
    "variability is exponential on clearance and the central volume, and",
    "residual error is proportional.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Por 2021 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Por ED, Akers KS, Chung KK, Livezey JR, Selig DJ.",
    "Population pharmacokinetic modeling and simulations of imipenem in",
    "burn patients with and without continuous venovenous hemofiltration",
    "in the military health system.",
    "J Clin Pharmacol. 2021;61(9):1182-1194. doi:10.1002/jcph.1865.",
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
        "Reference 99.5 kg. The cohort is heavy, as burn cohorts often are",
        "-- mean 89.6 +/- 22.38 kg with CVVH and 105.06 +/- 28.66 kg",
        "without (Zhang 2025 Table 1) -- and 99.5 sits between the two arm",
        "means. Enters the central volume with exponent 0.74, the non-CVVH",
        "clearance branch with exponent 0.33, and the CVVH clearance branch",
        "with the theory-based allometric exponent 0.75 (Zhang 2025 Table",
        "3). Note that the two clearance branches use DIFFERENT weight",
        "exponents; this is what the review prints and is reproduced",
        "verbatim."
      ),
      source_name        = "BW"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNIT CONVERSION. The source reference value is 2.7, which is",
        "serum albumin in g/dL -- the US convention, consistent with a US",
        "study, and the only reading on which 2.7 is a physiological",
        "albumin for a burn patient (2.7 g/dL = 27 g/L, markedly",
        "hypoalbuminaemic, as expected after major burns; 2.7 g/L would be",
        "incompatible with life). The canonical ALB column carries SI g/L",
        "per inst/references/covariate-columns.md, so model() converts",
        "inline with alb_gdL <- ALB * 0.1 before forming the covariate",
        "ratio, exactly as that register entry prescribes for models",
        "calibrated on g/dL.",
        "MAGNITUDE WARNING: the peripheral-volume exponent is -3.68, an",
        "extraordinarily steep dependence -- a halving of albumin would",
        "multiply the peripheral volume by 2^3.68 = 12.8. The review",
        "reports no standard error for it and the cohort is only 23",
        "patients, so the exponent is very likely poorly identified. It is",
        "transcribed verbatim but should be treated as unreliable outside",
        "the narrow albumin range actually observed; the review does not",
        "report that range. Recorded in the vignette Errata."
      ),
      source_name        = "ALB"
    ),
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
        "Reference 145.83 mL/min -- markedly supranormal, consistent with",
        "the augmented renal clearance characteristic of the hyperdynamic",
        "circulatory state after major burns, which Zhang 2025's Discussion",
        "invokes to explain why this study reports the highest clearance",
        "(15.31 L/h) of any two-compartment model in the review. Enters the",
        "NON-CVVH clearance branch only, as the power term",
        "(CLcr/145.83)^0.46; the CVVH branch carries no renal-function term",
        "at all (Zhang 2025 Table 3). Stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md."
      ),
      source_name        = "CLcr"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous venovenous haemofiltration during imipenem therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = not receiving CVVH",
      notes              = paste(
        "Zhang 2025 Table 3 prints two entirely separate clearance",
        "equations for this study rather than a single equation with a",
        "CVVH coefficient:",
        "without CVVH, CL = 15.31 * (CLcr/145.83)^0.46 * (BW/99.5)^0.33 *",
        "e^eta_CL; with CVVH, CL = 13.78 * (BW/99.5)^0.75 * e^eta_CL +",
        "1.56. The two branches differ in intercept, in weight exponent,",
        "in whether creatinine clearance enters at all, and in carrying an",
        "additive term. They are therefore encoded as a hard switch on this",
        "indicator rather than as a multiplicative covariate effect. CVVH",
        "is a continuous renal replacement modality, so the canonical",
        "RRT_CRRT_STATUS column applies per",
        "inst/references/covariate-columns.md, whose definition explicitly",
        "names 'continuous venovenous hemofiltration CVVH / CVVHF'."
      ),
      source_name        = "CVVH"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    AGE          = list(description = "Age",                        units = "years",  type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3). Cohort mean 55 +/- 19.99 years with CVVH, 51.09 +/- 19.03 without (Table 1)."),
    URINE_VOL_24H = list(description = "Urine output",              units = "mL/24h", type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3)."),
    TBSA         = list(description = "Total burned body surface area", units = "%",  type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3), along with the separately-screened total second-degree and total third-degree burn surface areas. Not a registered canonical in inst/references/covariate-columns.md; recorded here as documentation only, since it is never referenced in model() and no coefficient is reported.")
  )

  population <- list(
    species          = "human",
    n_subjects       = 23L,
    n_studies        = 1L,
    age_mean         = "55 +/- 19.99 years with CVVH; 51.09 +/- 19.03 years without CVVH (mean +/- SD)",
    weight_mean      = "89.6 +/- 22.38 kg with CVVH; 105.06 +/- 28.66 kg without CVVH (mean +/- SD)",
    sex_female_pct   = 73.9,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Adult burn patients receiving imipenem-cilastatin, some on",
      "continuous venovenous haemofiltration. This is the only burn cohort",
      "among the 18 studies in the Zhang 2025 review, and the only one in",
      "which women are the majority (17 of 23)."
    ),
    dose_range       = paste(
      "250 mg, 500 mg or 1000 mg imipenem intravenously every 6 h (Zhang",
      "2025 Supplementary Table S1). The infusion duration is not reported",
      "by the review."
    ),
    regions          = "United States of America",
    n_concentrations = 81L,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 12); 23 patients, 81",
      "samples, sex split 6 male / 17 female. Blood was sampled at trough",
      "and 0.5-8 h after administration, and assayed by HPLC-UV (Zhang",
      "2025 Supplementary Table S1). Covariates were selected by forward",
      "inclusion (dOFV > 3.84, p < 0.05) with no backward-elimination step",
      "reported. Fitted in Pumas -- the only study in the review to use it",
      "-- and evaluated by NPDE, bootstrap, VPC and, unusually for this",
      "review, external validation; only three of the 18 studies performed",
      "external validation (Zhang 2025 Table 2, Results).",
      "The review does not report how many of the 23 patients were on",
      "CVVH, only that both arms are present and that arm-specific",
      "demographics differ (Table 1). With 23 patients split across two",
      "clearance branches, every coefficient in this model rests on a very",
      "small sample; see the albumin magnitude warning in",
      "covariateData$ALB.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Por 2021 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2 and Table 3, Por et al.
    # (2021) row. Reference subject: 99.5 kg, albumin 2.7 g/dL (27 g/L),
    # Cockcroft-Gault CLcr 145.83 mL/min. =====
    #
    # Clearance has TWO branches. The names below follow the split-clearance
    # idiom already used in this package (Lamoth 2009 imipenem's
    # lcl_renal / lcl_nonren; Shekar 2014 meropenem's piecewise RRT switch):
    #   lcl      -- typical clearance of a patient NOT on CVVH
    #   lcl_cvvh -- endogenous clearance intercept of a patient ON CVVH
    #   lcl_crrt -- the additive filter clearance contributed by CVVH
    lcl      <- log(15.31); label("Clearance at the reference subject, no CVVH (L/h)")   # Zhang 2025 Table 3 (Por 2021): leading coefficient of the 'Without CVVH' CL formula; matches the Table 2 CL of 15.31 L/h
    lcl_cvvh <- log(13.78); label("Endogenous clearance intercept on CVVH at 99.5 kg (L/h)")  # Zhang 2025 Table 3 (Por 2021): leading coefficient of the 'With CVVH' CL formula
    lcl_crrt <- log(1.56);  label("Additive CVVH filter clearance (L/h)")                # Zhang 2025 Table 3 (Por 2021): the '+ 1.56' term of the 'With CVVH' CL formula
    lvc      <- log(32.67); label("Central volume of distribution at the reference subject (L)")     # Zhang 2025 Table 2 (Por 2021): V1 = 32.67 L; leading coefficient of the Table 3 Vc formula
    lvp      <- log(41.23); label("Peripheral volume of distribution at albumin 2.7 g/dL (L)")       # Zhang 2025 Table 2 (Por 2021): V2 = 41.23 L; leading coefficient of the Table 3 Vp formula
    lq       <- log(11);    label("Intercompartmental clearance Q (L/h)")                            # Zhang 2025 Table 2 (Por 2021): Q = 11 L/h

    # ===== Covariate effects -- Zhang 2025 Table 3, Por et al. (2021) =====
    #   Vc = 32.67 * (BW/99.5)^0.74 * (ALB/2.7)^(-1.17) * e^eta_Vc
    #   Vp = 41.23 * (ALB/2.7)^(-3.68)
    #   Without CVVH: CL = 15.31 * (CLcr/145.83)^0.46 * (BW/99.5)^0.33 * e^eta_CL
    #   With CVVH:    CL = 13.78 * (BW/99.5)^0.75 * e^eta_CL + 1.56
    e_crcl_cl     <- 0.46;  label("Power exponent on (CRCL/145.83) for CL, no-CVVH branch (unitless)")  # Zhang 2025 Table 3 (Por 2021): (CLcr/145.83)^0.46
    e_wt_cl       <- 0.33;  label("Power exponent on (WT/99.5) for CL, no-CVVH branch (unitless)")      # Zhang 2025 Table 3 (Por 2021): (BW/99.5)^0.33 in the 'Without CVVH' formula
    e_wt_cl_cvvh  <- 0.75;  label("Power exponent on (WT/99.5) for CL, CVVH branch (unitless)")         # Zhang 2025 Table 3 (Por 2021): (BW/99.5)^0.75 in the 'With CVVH' formula
    e_wt_vc       <- 0.74;  label("Power exponent on (WT/99.5) for Vc (unitless)")                      # Zhang 2025 Table 3 (Por 2021): (BW/99.5)^0.74
    e_alb_vc      <- -1.17; label("Power exponent on (ALB/2.7 g/dL) for Vc (unitless)")                 # Zhang 2025 Table 3 (Por 2021): (ALB/2.7)^(-1.17)
    e_alb_vp      <- -3.68; label("Power exponent on (ALB/2.7 g/dL) for Vp (unitless)")                 # Zhang 2025 Table 3 (Por 2021): (ALB/2.7)^(-3.68)

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
    etalcl ~ log(1 + 0.305^2)  # Zhang 2025 Table 2 (Por 2021): IIV CL = 30.5%, read as an apparent CV
    etalvc ~ log(1 + 0.361^2)  # Zhang 2025 Table 2 (Por 2021): IIV V1 = 36.1%, read as an apparent CV

    # ===== Residual error =====
    propSd <- 0.30; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Por 2021): Proportional = 30%
  })

  model({
    # ----- Unit conversion -----
    # The albumin covariate was calibrated on g/dL (reference 2.7); the
    # canonical ALB column carries SI g/L. See covariateData$ALB.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL

    # ----- Individual PK parameters -----
    # Clearance switches hard between the two branches printed in Zhang
    # 2025 Table 3. The exponential eta multiplies the endogenous arm of
    # whichever branch applies; in the CVVH branch the additive filter
    # clearance of 1.56 L/h is added AFTER the eta, exactly as printed.
    cl <- ((1 - RRT_CRRT_STATUS) * exp(lcl)      * (CRCL / 145.83)^e_crcl_cl * (WT / 99.5)^e_wt_cl +
                RRT_CRRT_STATUS  * exp(lcl_cvvh) *                             (WT / 99.5)^e_wt_cl_cvvh) *
          exp(etalcl) +
          RRT_CRRT_STATUS * exp(lcl_crrt)

    vc <- exp(lvc + etalvc) * (WT / 99.5)^e_wt_vc * (alb_gdL / 2.7)^e_alb_vc
    vp <- exp(lvp)                                * (alb_gdL / 2.7)^e_alb_vp
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
