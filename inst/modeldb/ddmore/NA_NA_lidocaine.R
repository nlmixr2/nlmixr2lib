# Population PK model for lidocaine and three sequential metabolites
# (monoethylglycinexylidide [MEGX], glycinexylidide [GX], and 2,6-xylidide
# [2,6-XYL]) extracted from the DDMORE Foundation Model Repository entry
# DDMODEL00000281. The bundle does not link to a journal publication, so the
# only authoritative sources for both structure and parameter values are the
# bundle's `Executable_ddmore_final_run249.ctl` (NONMEM control stream;
# structure + initial values) and `Output_real_data_original_final_run249.res`
# (NONMEM listing; final estimates after `MINIMIZATION SUCCESSFUL`).
#
# DDMORE_Model_Accomodations.txt asserts that the uploaded model and the model
# in the (un-named) reference publication do not differ. The bundle's `.ctl`
# and `.res` $PROBLEM line is `B.dat 4-cRUN249`; the License is registered to
# BAST Inc. Ltd, suggesting the run originates from a BAST-led lidocaine
# study, but no first-author / year is recoverable from the on-disk material.
# The model file therefore uses the operator-supplied placeholder filename
# `NA_NA_lidocaine.R` and reference text "DDMORE Foundation Model Repository:
# DDMODEL00000281. No linked publication identified."

NA_NA_lidocaine <- function() {
  description <- paste(
    "Population PK model for lidocaine and three sequential metabolites",
    "(monoethylglycinexylidide [MEGX], glycinexylidide [GX], and 2,6-xylidide",
    "[2,6-XYL]) using the source's NONMEM ADVAN5 / TRANS1 general-linear",
    "rate-constant parameterisation. Lidocaine in the central compartment is",
    "metabolised in parallel to MEGX (rate constant k_megx, fixed at 0.03) and",
    "to 2,6-XYL (k_xyl, fixed at 0.007); MEGX is sequentially metabolised to",
    "GX (k_gx, estimated). GX and 2,6-XYL each carry a typical-value",
    "elimination rate constant (kel_gx, kel_xyl) modulated by binary",
    "stratifications of dose-level (DLVL > 2), bilirubin (BIL > 0.53),",
    "creatinine clearance (CRCL <= 52.7), CYP1A2-modifying co-medications",
    "(S1A2 == 3), body-mass index (BMI > 27.93), serum ALT (SGPT > 11), and",
    "lactate dehydrogenase (LDH > 195). Lidocaine apparent central volume V1",
    "is also dose-level-stratified (DLVL > 2); the three metabolite",
    "compartments share a fixed apparent volume of 100. Distributed in the",
    "DDMORE Foundation Model Repository as DDMODEL00000281; no linked",
    "publication is identified in the bundle's `Model_Accommodations.txt` or",
    "in the `.ctl` / `.res` headers. Final parameter estimates come from the",
    "bundle's `Output_real_data_original_final_run249.res` listing after the",
    "`MINIMIZATION SUCCESSFUL` block (FOCE estimation, OBJ = 10219.922,",
    "covariance step succeeded). The bundle's `.ctl` does not declare time,",
    "dose, or concentration units explicitly; `units$time = 'h'`,",
    "`units$dosing = 'mg'`, and `units$concentration = 'mg/L'` are",
    "operator-default placeholders chosen so the values flow through unit-",
    "checking consistently. See the vignette Errata for the unit ambiguity",
    "and a per-time-point self-consistency check against the bundle's",
    "Simulated_Lid_B04_ddmore.csv."
  )
  reference <- "DDMORE Foundation Model Repository: DDMODEL00000281. No linked publication identified in the bundle (Model_Accommodations.txt asserts equivalence to an unspecified reference publication). Source $PROBLEM line `B.dat 4-cRUN249`; license registered to BAST Inc. Ltd; NONMEM run dated 29/11/2016."
  vignette <- "NA_NA_lidocaine"
  ddmore_id <- "DDMODEL00000281"
  replicate_of <- NULL
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    DLVL = list(
      description        = "Source-protocol integer dose-level / regimen indicator (1-4 in the bundle's simulated dataset). Used as binary `DLVL_HIGH = as.integer(DLVL > 2)` per the source `.ctl` `IF(DLVL.GT.2)P1=0` line; switches typical-value baselines for the GX rate constant k_gx_elim and the lidocaine apparent central volume vc.",
      units              = "(integer-coded categorical)",
      type               = "categorical",
      reference_category = "DLVL <= 2 (baseline regimen; `THETA(4)` for k_gx_elim base and `THETA(14)` for vc).",
      notes              = "Carried per subject (time-fixed in the bundle's simulated dataset). `DLVL > 2` selects the higher-baseline-rate / higher-volume regimen (`THETA(5)` and `THETA(15)` respectively in the source `.ctl`). The exact biological / protocol meaning of each integer level is not recoverable from the bundle; the binary threshold matches the source.",
      source_name        = "DLVL"
    ),
    TBILI = list(
      description        = "Total serum bilirubin concentration. Used as binary `BIL_HIGH = as.integer(TBILI > 0.53)` per the source `.ctl` `IF(BIL.GT.0.53)P2=0` line; adds an additive linear modifier (THETA(6) = -0.529) to the typical-value GX elimination rate constant kel_gx in the elevated-bilirubin cohort.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = "TBILI <= 0.53 (additive modifier off; the high-bilirubin effect is +0 to kel_gx).",
      notes              = "Source column name `BIL`; canonical name `TBILI` (the `inst/references/covariate-columns.md` `TBILI` entry registers `BIL` as an alias). Rename `BIL -> TBILI` before passing the dataset to `rxSolve`. Threshold 0.53 mg/dL is well within the clinical normal range (<=1.2 mg/dL) so the binarisation reflects a paper-specific cohort split, not a clinical hepatic-impairment cutoff.",
      source_name        = "BIL"
    ),
    LDH = list(
      description        = "Serum lactate dehydrogenase activity. Used as binary `LDH_HIGH = as.integer(LDH > 195)` per the source `.ctl` `IF(LDH.GT.195)P3=0` line; switches typical-value baseline of the 2,6-xylidide elimination rate constant kel_xyl between a low-LDH (THETA(11) = 0.667) and a high-LDH (THETA(12) = 0.410) regimen.",
      units              = "U/L",
      type               = "continuous",
      reference_category = "LDH <= 195 (baseline rate constant THETA(11) = 0.667).",
      notes              = "Threshold 195 U/L sits at the upper end of the clinical reference range (140-280 U/L). Binarisation rather than continuous power form, matching the source `.ctl` exactly.",
      source_name        = "LDH"
    ),
    CRCL = list(
      description        = "Creatinine-based renal function. Used as binary `CRCL_LOW = as.integer(CRCL <= 52.7)` per the source `.ctl` `IF(CRCL.LE.52.7)P4=0` line; adds an additive linear modifier (THETA(7) = -0.319) to the typical-value GX rate constant k_gx_elim in the renal-impaired cohort.",
      units              = "mL/min (BSA-normalisation method not stated in the bundle)",
      type               = "continuous",
      reference_category = "CRCL > 52.7 (additive modifier off).",
      notes              = "The source `.ctl` does not specify whether CRCL is BSA-normalised or how it was estimated (Cockcroft-Gault vs MDRD vs measured). Documented as `mL/min` (raw clinical units) here; downstream re-users should consult the linked publication if a unit-normalisation question arises.",
      source_name        = "CRCL"
    ),
    S1A2 = list(
      description        = "Source-protocol CYP1A2-modifying co-medication / phenotype categorical indicator (integer 0-3). Used as binary `S1A2_IND = as.integer(S1A2 == 3)` per the source `.ctl` `IF(S1A2.EQ.3)P5=0` line; adds an additive linear modifier (THETA(8) = +0.853) to the typical-value GX rate constant k_gx_elim in the level-3 cohort.",
      units              = "(integer-coded categorical)",
      type               = "categorical",
      reference_category = "S1A2 != 3 (values 0, 1, 2 are pooled into the reference; additive modifier off).",
      notes              = "Carried per subject (time-fixed in the bundle's simulated dataset). The natural interpretation, given the column name encodes 'CYP1A2' and the model attaches a sizeable positive K30 modifier of +0.853 to the level-3 cohort, is a CYP1A2-induction or smoking / inducer co-medication indicator. Sibling columns `D1A2` and `H1A2` are dropped in the source `.ctl`, so only the level-3 indicator is structurally identifiable. The exact biological meaning of each integer level is not fully reconstructable from the bundle.",
      source_name        = "S1A2"
    ),
    BMI = list(
      description        = "Body mass index. Used as binary `BMI_HIGH = as.integer(BMI > 27.93)` per the source `.ctl` `IF(BMI.GT.27.93)P7=0` line; adds an additive linear modifier (THETA(9) = +0.939) to the typical-value GX rate constant k_gx_elim in the high-BMI cohort.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = "BMI <= 27.93 (additive modifier off).",
      notes              = "Threshold 27.93 kg/m^2 sits between WHO overweight (>=25) and obesity (>=30) cutoffs; the value most likely reflects the source-cohort median rather than a clinical-guideline threshold. Binarisation rather than continuous power form, matching the source `.ctl` exactly.",
      source_name        = "BMI"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity. Used as binary `SGPT_HIGH = as.integer(ALT > 11)` per the source `.ctl` `IF(SGPT.GT.11)P6=0` line; adds an additive linear modifier (THETA(10) = -0.492) to the typical-value GX rate constant kel_gx AND a separate modifier (THETA(13) = +0.229) to the typical-value 2,6-xylidide rate constant kel_xyl, in the elevated-ALT cohort.",
      units              = "U/L",
      type               = "continuous",
      reference_category = "ALT <= 11 (additive modifiers off).",
      notes              = "Source column name `SGPT` (legacy serum glutamic-pyruvic transaminase label); canonical name `ALT` (the `inst/references/covariate-columns.md` `ALT` entry registers `SGPT` as an alias paralleling `SGOT` -> `AST`). Rename `SGPT -> ALT` before passing the dataset to `rxSolve`. Threshold 11 U/L is below the lower end of the clinical reference range (~7-56 U/L for adults), so the binarisation almost certainly reflects a paper-specific cohort split rather than a clinical hepatic-impairment cutoff; the linked publication is not on disk to confirm.",
      source_name        = "SGPT"
    )
  )

  population <- list(
    n_subjects     = 325L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Patient population not stated in the DDMORE bundle. The `.res` listing reports 325 subjects contributing 1989 observations; the bundle's simulated dataset (`Simulated_Lid_B04_ddmore.csv`) has subjects receiving repeated short IV infusions of lidocaine consistent with surgical / intensive-care or anti-arrhythmic dosing. The linked publication is not on disk to confirm the indication.",
    dose_range     = "Repeated IV infusions of approximately 12 time-units' duration (AMT 21600 / RATE 1800 in the bundle's simulated dataset). Mass and time units are not declared in the source `.ctl`; under the operator-chosen `units$time = 'h'` interpretation each infusion runs ~12 h.",
    regions        = NA_character_,
    notes          = "Demographics fields marked NA because the linked publication is not on disk for this extraction. n_subjects = 325 from the `.res` listing's `TOT. NO. OF INDIVIDUALS:    325` line. The DDMORE-shipped simulated dataset (`Simulated_Lid_B04_ddmore.csv`) carries 17112 records distributed over a smaller demographic-replicated cohort and is intended only as a regression-style smoke test, not a representative clinical population."
  )

  ini({
    # Lidocaine -> MEGX (k_megx) and lidocaine -> 2,6-XYL (k_xyl) rate
    # constants. Both are FIXED in the source `.ctl` `$THETA` block (lines
    # `(0.03 FIX) ;K12` and `(0.007 FIX) ;K14`) and re-emitted unchanged in
    # the `.res` `FINAL PARAMETER ESTIMATE` block as TH 1 = 3.00E-02 and
    # TH 3 = 7.00E-03. Time-unit interpretation: see vignette Errata; the
    # `units$time = "h"` choice gives lidocaine total elimination
    # k_megx + k_xyl = 0.037 1/h with apparent half-life ln(2)/0.037 ~
    # 18.7 h, slower than the textbook lidocaine IV t1/2 (~1.5-2 h), which
    # the operator notes as a deviation pending publication recovery.
    lk_megx_form <- fixed(log(0.03))
    label("Log lidocaine -> MEGX formation rate constant (1/time-unit) - FIXED")  # `.ctl` $THETA TH 1 (FIX); `.res` FINAL TH 1 = 3.00E-02
    lk_xyl_form  <- fixed(log(0.007))
    label("Log lidocaine -> 2,6-xylidide formation rate constant (1/time-unit) - FIXED")  # `.ctl` $THETA TH 3 (FIX); `.res` FINAL TH 3 = 7.00E-03

    # MEGX -> GX rate constant (k_gx_form). Estimated; `.res` FINAL TH 2 =
    # 1.93. No covariate effects in the source.
    lk_gx_form <- log(1.93)
    label("Log MEGX -> GX formation rate constant (1/time-unit)")  # `.res` FINAL TH 2 = 1.93E+00 (init `.ctl` $THETA TH 2 = (0,1))

    # GX elimination rate constant (kel_gx). Two typical-value baselines
    # depending on DLVL: THETA(4) for DLVL <= 2, THETA(5) for DLVL > 2.
    # The .ctl `T1K30 = P1*THETA(4) + (1-P1)*THETA(5)` selects the baseline
    # additively, equivalent to `kel_gx_base = exp(lkel_gx + e_dlvl_high *
    # DLVL_HIGH)` where the effect is on log scale. Final estimates from
    # the `.res` FINAL PARAMETER ESTIMATE block: TH 4 = 1.44, TH 5 = 2.07.
    lkel_gx <- log(1.44)
    label("Log typical-value GX elimination rate constant for DLVL <= 2 (1/time-unit)")  # `.res` FINAL TH 4 = 1.44E+00
    e_dlvl_high_kel_gx <- log(2.07 / 1.44)
    label("Log fold-change in typical-value GX elimination rate constant for DLVL > 2 (unitless)")  # derived from `.res` FINAL TH 5 / TH 4 = 2.07 / 1.44

    # K30 (GX elimination) additive linear modifiers. Source applies these
    # ON THE LINEAR k_gx_elim SCALE (not log scale) sequentially in the
    # `.ctl` `$PK` block:
    #   T2K30 = T1K30 + (1-P2) * THETA(6)         ; bilirubin (BIL > 0.53)
    #   T3K30 = T2K30 + (1-P4) * THETA(7)         ; CRCL <= 52.7
    #   T4K30 = T3K30 + (1-P5) * THETA(8)         ; S1A2 == 3
    #   T5K30 = T4K30 + (1-P7) * THETA(9)         ; BMI > 27.93
    #   TK30  = T5K30 + (1-P6) * THETA(10)        ; SGPT > 11
    # with `IF(...LE.0)...=0.0001` numerical floors at each intermediate
    # step. Names follow the e_<cov>_<param> naming convention; values are
    # the THETA estimates verbatim from the `.res` listing. Final values:
    # TH 6 = -0.529, TH 7 = -0.319, TH 8 = 0.853, TH 9 = 0.939, TH 10 = -0.492.
    e_bil_kel_gx  <- -0.529
    label("Additive linear modifier on typical-value GX rate constant when BIL > 0.53 (1/time-unit)")    # `.res` FINAL TH 6 = -5.29E-01
    e_crcl_kel_gx <- -0.319
    label("Additive linear modifier on typical-value GX rate constant when CRCL <= 52.7 (1/time-unit)")  # `.res` FINAL TH 7 = -3.19E-01
    e_s1a2_kel_gx <-  0.853
    label("Additive linear modifier on typical-value GX rate constant when S1A2 == 3 (1/time-unit)")     # `.res` FINAL TH 8 = +8.53E-01
    e_bmi_kel_gx  <-  0.939
    label("Additive linear modifier on typical-value GX rate constant when BMI > 27.93 (1/time-unit)")   # `.res` FINAL TH 9 = +9.39E-01
    e_sgpt_kel_gx <- -0.492
    label("Additive linear modifier on typical-value GX rate constant when SGPT > 11 (1/time-unit)")     # `.res` FINAL TH 10 = -4.92E-01

    # 2,6-XYL elimination rate constant (kel_xyl). Two typical-value
    # baselines depending on LDH: THETA(11) for LDH <= 195, THETA(12) for
    # LDH > 195. `.res` FINAL TH 11 = 0.667, TH 12 = 0.410.
    lkel_xyl <- log(0.667)
    label("Log typical-value 2,6-xylidide elimination rate constant for LDH <= 195 (1/time-unit)")  # `.res` FINAL TH 11 = +6.67E-01
    e_ldh_high_kel_xyl <- log(0.410 / 0.667)
    label("Log fold-change in typical-value 2,6-xylidide elimination rate constant for LDH > 195 (unitless)")  # derived from `.res` FINAL TH 12 / TH 11 = 0.410 / 0.667

    # K40 (2,6-XYL elimination) additive linear modifier for SGPT > 11.
    # Source applies on the linear scale: `TK40 = T1K40 + (1-P6)*THETA(13)`
    # with `IF(...LE.0)...=0.0001` floor. `.res` FINAL TH 13 = +0.229.
    e_sgpt_kel_xyl <- 0.229
    label("Additive linear modifier on typical-value 2,6-xylidide rate constant when SGPT > 11 (1/time-unit)")  # `.res` FINAL TH 13 = +2.29E-01

    # Lidocaine apparent central volume V1 (vc). Two typical-value
    # baselines depending on DLVL: THETA(14) for DLVL <= 2, THETA(15) for
    # DLVL > 2. `.res` FINAL TH 14 = 1.32E+03, TH 15 = 1.81E+03.
    lvc <- log(1320)
    label("Log lidocaine apparent central volume V1 for DLVL <= 2 (volume-unit)")  # `.res` FINAL TH 14 = +1.32E+03
    e_dlvl_high_vc <- log(1810 / 1320)
    label("Log fold-change in lidocaine apparent central volume V1 for DLVL > 2 (unitless)")  # derived from `.res` FINAL TH 15 / TH 14 = 1810 / 1320

    # Shared apparent volume of the three metabolite compartments (Vm).
    # FIXED at 100 in the source `.ctl` `$THETA` block (`(100 FIX) ;V2,3,4`)
    # and re-emitted unchanged in `.res` FINAL TH 16 = 1.00E+02. Declared
    # as three separate parameters here, one per metabolite compartment,
    # all numerically identical and all FIXED, so the per-metabolite
    # central volumes can later be relaxed if the model structure is
    # extended.
    lvc_megx <- fixed(log(100))
    label("Log MEGX apparent central volume (volume-unit) - FIXED, shared across metabolites")  # `.res` FINAL TH 16 (FIX) = 1.00E+02
    lvc_gx   <- fixed(log(100))
    label("Log GX apparent central volume (volume-unit) - FIXED, shared across metabolites")    # `.res` FINAL TH 16 (FIX) = 1.00E+02
    lvc_xyl  <- fixed(log(100))
    label("Log 2,6-xylidide apparent central volume (volume-unit) - FIXED, shared across metabolites")  # `.res` FINAL TH 16 (FIX) = 1.00E+02

    # Inter-individual variability. Source `$OMEGA` is diagonal (no BLOCK).
    # `.res` FINAL OMEGA diagonals: ETA1 = 3.91E-01 on K30 (linear-scale
    # exponential IIV: K30 = TK30 * EXP(ETA1)); ETA2 = 2.00E-01 on K40
    # (same form, K40 = TK40 * EXP(ETA2)); ETA3 = 3.11E-01 on V1
    # (V1 = TV1 * EXP(ETA3)). The eta names attach to the log-scale
    # parameters because the `* EXP(ETA)` form is identical to adding the
    # eta on the log scale.
    etalkel_gx  ~ 0.391  # `.res` FINAL OMEGA ETA1 (diagonal) = 3.91E-01
    etalkel_xyl ~ 0.200  # `.res` FINAL OMEGA ETA2 (diagonal) = 2.00E-01
    etalvc      ~ 0.311  # `.res` FINAL OMEGA ETA3 (diagonal) = 3.11E-01

    # Residual error. Source `$ERROR` block uses additive epsilon on each
    # observed concentration: Y = F + EPS(<i>) with one EPS per observation
    # type (LID = CMT 1, MEGX = CMT 2, GX = CMT 3, 2,6-XYL = CMT 4).
    # `.res` FINAL SIGMA diagonals (variances): EPS1 = 3.64E+02, EPS2 =
    # 5.33E+01, EPS3 = 4.79E+01, EPS4 = 6.39E+00. Standard deviations
    # therefore: SD_LID = sqrt(364) ~ 19.08; SD_MEGX = sqrt(53.3) ~ 7.30;
    # SD_GX = sqrt(47.9) ~ 6.92; SD_XYL = sqrt(6.39) ~ 2.53.
    addSd      <- sqrt(364)
    label("Additive residual SD on lidocaine concentration (concentration-unit)")     # `.res` FINAL SIGMA EPS1 (diagonal) = 3.64E+02 (variance) -> SD = sqrt(364)
    addSd_megx <- sqrt(53.3)
    label("Additive residual SD on MEGX concentration (concentration-unit)")          # `.res` FINAL SIGMA EPS2 (diagonal) = 5.33E+01 (variance) -> SD = sqrt(53.3)
    addSd_gx   <- sqrt(47.9)
    label("Additive residual SD on GX concentration (concentration-unit)")            # `.res` FINAL SIGMA EPS3 (diagonal) = 4.79E+01 (variance) -> SD = sqrt(47.9)
    addSd_xyl  <- sqrt(6.39)
    label("Additive residual SD on 2,6-xylidide concentration (concentration-unit)")  # `.res` FINAL SIGMA EPS4 (diagonal) = 6.39E+00 (variance) -> SD = sqrt(6.39)
  })

  model({
    # 1. Derived covariate binarisations. Each line reproduces a single
    # `IF(...).<comparator>.<value>) P<n> = 0` line from the source `.ctl`
    # `$PK` block, expressed as the `_HIGH` / `_LOW` / `_IND` indicator the
    # additive modifiers expect.
    DLVL_HIGH <- (DLVL > 2)
    BIL_HIGH  <- (TBILI > 0.53)
    LDH_HIGH  <- (LDH > 195)
    CRCL_LOW  <- (CRCL <= 52.7)
    S1A2_IND  <- (S1A2 == 3)
    BMI_HIGH  <- (BMI > 27.93)
    SGPT_HIGH <- (ALT > 11)

    # 2. Individual rate constants. Lidocaine -> MEGX and lidocaine ->
    # 2,6-XYL formation rate constants are fixed (no IIV, no covariates).
    k_megx_form <- exp(lk_megx_form)
    k_xyl_form  <- exp(lk_xyl_form)
    # MEGX -> GX rate constant (no IIV, no covariates in source).
    k_gx_form   <- exp(lk_gx_form)

    # GX elimination rate constant kel_gx. Build the typical value by
    # applying:
    #   (a) the DLVL log-scale fold-change (multiplicative on the typical),
    #   (b) the five additive linear modifiers in the source's stepwise
    #       order, with the `LE.0` floor at each intermediate step.
    # IIV is exponential around the typical value, K30 = TK30 * EXP(ETA1)
    # in the source; equivalently kel_gx = typ_kel_gx * exp(etalkel_gx).
    typ1_kel_gx <- exp(lkel_gx + e_dlvl_high_kel_gx * DLVL_HIGH)
    typ2_kel_gx <- typ1_kel_gx + e_bil_kel_gx  * BIL_HIGH
    typ2_kel_gx <- max(typ2_kel_gx, 0.0001)
    typ3_kel_gx <- typ2_kel_gx + e_crcl_kel_gx * CRCL_LOW
    typ3_kel_gx <- max(typ3_kel_gx, 0.0001)
    typ4_kel_gx <- typ3_kel_gx + e_s1a2_kel_gx * S1A2_IND
    typ4_kel_gx <- max(typ4_kel_gx, 0.0001)
    typ5_kel_gx <- typ4_kel_gx + e_bmi_kel_gx  * BMI_HIGH
    typ5_kel_gx <- max(typ5_kel_gx, 0.0001)
    typ_kel_gx  <- typ5_kel_gx + e_sgpt_kel_gx * SGPT_HIGH
    typ_kel_gx  <- max(typ_kel_gx, 0.0001)
    kel_gx      <- typ_kel_gx * exp(etalkel_gx)

    # 2,6-xylidide elimination rate constant kel_xyl. Build the typical
    # value by applying the LDH log-scale fold-change then the SGPT
    # additive modifier with the `LE.0` floor.
    typ1_kel_xyl <- exp(lkel_xyl + e_ldh_high_kel_xyl * LDH_HIGH)
    typ_kel_xyl  <- typ1_kel_xyl + e_sgpt_kel_xyl * SGPT_HIGH
    typ_kel_xyl  <- max(typ_kel_xyl, 0.0001)
    kel_xyl      <- typ_kel_xyl * exp(etalkel_xyl)

    # Lidocaine apparent central volume V1 with DLVL log-fold-change and
    # log-normal IIV.
    vc      <- exp(lvc + e_dlvl_high_vc * DLVL_HIGH + etalvc)
    vc_megx <- exp(lvc_megx)
    vc_gx   <- exp(lvc_gx)
    vc_xyl  <- exp(lvc_xyl)

    # 3. ODE system. Source `$SUBR ADVAN5 TRANS1` with `$MODEL`
    # COMP=(CENTRAL DEFDOSE DEFOBSERVATION), COMP=(MEGX), COMP=(GX),
    # COMP=(26XYL), and the implicit OUTPUT compartment 5. The rate-matrix
    # row in the `.res` listing maps:
    #   FROM 1 TO 2 = K12 (lidocaine -> MEGX)            = k_megx_form
    #   FROM 1 TO 4 = K14 (lidocaine -> 2,6-xylidide)    = k_xyl_form
    #   FROM 2 TO 3 = K23 (MEGX -> GX)                   = k_gx_form
    #   FROM 3 TO 5 = K30 (GX -> output)                 = kel_gx
    #   FROM 4 TO 5 = K40 (2,6-xylidide -> output)       = kel_xyl
    # Lidocaine has no direct elimination pathway in this model; total
    # apparent lidocaine elimination = k_megx_form + k_xyl_form.
    d/dt(central)      <- -(k_megx_form + k_xyl_form) * central
    d/dt(central_megx) <-  k_megx_form * central - k_gx_form * central_megx
    d/dt(central_gx)   <-  k_gx_form   * central_megx - kel_gx * central_gx
    d/dt(central_xyl)  <-  k_xyl_form  * central - kel_xyl * central_xyl

    # 4. Observation variables. Source `$ERROR` block uses S1=V1, S2=V2,
    # S3=V3, S4=V4 to scale each compartment's mass to a concentration
    # (compartment_mass / compartment_volume) before applying the additive
    # epsilon. `Y = F + EPS(i)` for the matching CMT switches between
    # output types via the `IF(CMT.EQ.<n>)` block; in nlmixr2 that maps to
    # one residual-error line per output, with the source columns kept on
    # the canonical `Cc` / `Cc_<metab>` names.
    Cc      <- central      / vc
    Cc_megx <- central_megx / vc_megx
    Cc_gx   <- central_gx   / vc_gx
    Cc_xyl  <- central_xyl  / vc_xyl

    Cc      ~ add(addSd)
    Cc_megx ~ add(addSd_megx)
    Cc_gx   ~ add(addSd_gx)
    Cc_xyl  ~ add(addSd_xyl)
  })
}
