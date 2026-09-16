Gao_2026_trastuzumabRezetecan <- function() {
  description <- "Sequential two-analyte population PK model for trastuzumab rezetecan (SHR-A1811, a HER2-targeting antibody-drug conjugate with a drug-to-antibody ratio of approximately 6.0; output Cc) and its released topoisomerase-I-inhibitor payload rezetecan (output Cc_rez) in adults with HER2-expressing or HER2-mutated advanced solid tumors (Gao 2026). The intact ADC is a two-compartment model with linear elimination after IV infusion. The released payload is a one-compartment model whose formation is a first-order release from the intact ADC in the central compartment and whose elimination is linear; the payload compartment does not feed back on the ADC, matching the paper's sequential two-step estimation. The release-rate constant Krel equals RAT during Cycle 1 and RAT * ALPHA (ALPHA = 0.693) from Cycle 2 onward. Covariates are body weight and baseline tumor size on ADC clearance; body weight, age and baseline albumin on ADC central volume; baseline albumin on ADC peripheral volume; body weight, baseline tumor size and cancer type on the release rate; age on payload volume; and aspartate aminotransferase on payload clearance."
  reference <- "Gao X, Zhao K, Zhao Y, Zhang Y, Zhao J, Zhao C, Djebli N. Population Pharmacokinetics of Trastuzumab Rezetecan in Patients With HER2-Expressing or Mutated Advanced Solid Tumors. CPT Pharmacometrics Syst Pharmacol. 2026. doi:10.1002/psp4.70259. PMCID PMC13274680."
  vignette <- "Gao_2026_trastuzumabRezetecan"

  # Gao 2026 reports intact-ADC concentrations in ug/mL (assay LLOQ 1.00
  # ug/mL) and released-payload concentrations in ng/mL (assay LLOQ 0.05
  # ng/mL). The ADC subsystem is encoded on the paper's own mass scale:
  # dose in mg and volumes in L give Cc directly in mg/L = ug/mL, which is
  # the scale the additive residual SD below is expressed on.
  #
  # The payload scale is NOT fully recoverable from the paper - see the
  # `mwr` note in model() and the vignette 'Assumptions and deviations'.
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "trastuzumab rezetecan (intact ADC)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "trastuzumab rezetecan (intact ADC)", units = "mg", specimen = "serum", verified = TRUE),
    central_rez = list(analyte = "released rezetecan payload", units = "mg ADC-molar-equivalent (see model() `mwr` note)", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on intact-ADC CL (exponent 0.525) and V1 (exponent 0.544), and on the payload release-rate constant RAT (exponent -0.546). Reference 60.6 kg, the value printed as the typical patient in the Gao 2026 Figure 5 and Figure 6 captions and appearing as the denominator of every body-weight term in the Section 3.2 and Section 3.3 equations. Dosing is weight-based (1.0-8.0 mg/kg), so ADC steady-state AUC scales as BW^(1 - 0.525) rather than BW^(-0.525).",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on intact-ADC V1 (exponent 0.198) and on payload V3 (exponent 0.420). Reference 56 years (Gao 2026 Section 3.2 and 3.3 equations; the Figure 5 / Figure 6 captions define the typical patient as 56 years old).",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on intact-ADC V1 (exponent -0.363) and V2 (exponent -1.44); both volumes decrease as albumin rises. Reference 42.5 g/L (Gao 2026 Section 3.2 equations and the Figure 5 caption). Gao 2026 Table 1 abbreviation list confirms the unit is g/L, so no g/dL conversion is required.",
      source_name        = "ALB"
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on released-payload clearance CL3 (exponent -0.167); payload clearance falls, and hence payload exposure rises, as AST rises. Reference 25 U/L (Gao 2026 Section 3.3 equation and the Figure 6 caption).",
      source_name        = "AST"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size: sum of diameters of target lesions at baseline (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on intact-ADC CL (exponent 0.0686) and on the payload release-rate constant RAT (exponent 0.100). Reference 52 mm (Gao 2026 Section 3.2 and 3.3 equations; the Figure 5 / Figure 6 captions give SOD 52 mm for the typical patient). Linear sum-of-diameters construct in mm, not an SPPD area. Source column SOD_B.",
      source_name        = "SOD_B"
    ),
    TUMTP_BREAST = list(
      description        = "Tumor-type indicator: 1 = breast cancer (BC), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NSCLC is the model reference when TUMTP_BREAST, TUMTP_GASTRIC and TUMTP_OTHER are all 0)",
      notes              = "Gao 2026 Table 1 note: 'CT1, CT2 and CT3 represent Breast Cancer (BC), Gastric or Gastroesophageal Junction cancer (GC/GEJ) and Other Tumor types, respectively.' The effect is ADDITIVE on the release-rate constant RAT (units 1/day), applied before the body-weight and tumor-size power terms: RAT_typical = 0.814 + 0.0562 for breast cancer. Breast cancer is the largest group in the analysis (approximately 60% of the 645 patients).",
      source_name        = "CT1"
    ),
    TUMTP_GASTRIC = list(
      description        = "Tumor-type indicator: 1 = gastric cancer or gastroesophageal junction (GEJ) adenocarcinoma, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NSCLC is the model reference when TUMTP_BREAST, TUMTP_GASTRIC and TUMTP_OTHER are all 0)",
      notes              = "Gao 2026 Table 1 CT2. ADDITIVE effect on the release-rate constant RAT (units 1/day): RAT_typical = 0.814 - 0.129 for GC/GEJ. The canonical TUMTP_GASTRIC register entry already covers gastric cancer OR adenocarcinoma of the gastroesophageal junction, which is exactly the Gao 2026 CT2 definition.",
      source_name        = "CT2"
    ),
    TUMTP_OTHER = list(
      description        = "Tumor-type indicator: 1 = other tumor types (neither NSCLC, breast, nor gastric/GEJ), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NSCLC is the model reference when TUMTP_BREAST, TUMTP_GASTRIC and TUMTP_OTHER are all 0)",
      notes              = "Gao 2026 Table 1 CT3. ADDITIVE effect on the release-rate constant RAT (units 1/day): RAT_typical = 0.814 + 0.0464 for other tumor types. NOTE the reference group here is NSCLC, NOT the 'other' pool: Gao 2026 Section 3.3 prints the NSCLC equation with the bare 0.814 and gives BC, GC/GEJ and Others their own additive shifts. This is the inverse orientation from papers that treat 'Other' as the residual reference, so the scope is paper-specific.",
      source_name        = "CT3"
    ),
    CYCLE = list(
      description        = "Treatment cycle number (1 = first 21-day cycle, 2 = second, ...; integer count, time-varying across the treatment course)",
      units              = "(count)",
      type               = "count",
      reference_category = "n/a -- used as the piecewise indicator CYCLE == 1 versus CYCLE >= 2",
      notes              = "Required for the released-payload sub-model only. Gao 2026 Section 3.3: 'Krel = RAT during Cycle 1 and Krel = RAT * ALPHA after Cycle 1' with ALPHA = 0.693, i.e. the release rate drops to 69.3% of its Cycle-1 value from Cycle 2 onward. Cycle length is 21 days (Gao 2026 Section 2.1), so CYCLE increments every 21 days on the Q3W regimen. Does not affect intact-ADC disposition or payload elimination. Gao 2026 states that a release rate varying continuously with time or cycle was evaluated and did NOT improve the fit, so the step change is the paper's selected form.",
      source_name        = "CYCLE"
    )
  )

  # Screened but not retained in the final model: Gao 2026 Section 3.4/3.5
  # and Figures S1-S7 report that race, sex, formulation, and hepatic and
  # renal function categories had no clinically relevant impact on either
  # analyte's exposure, so no point estimates are published for them.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Assessed as a covariate and compared post hoc by subgroup (Gao 2026 Figure S4); not retained in the final model and no point estimate published."
    ),
    RACE_ASIAN = list(
      description = "Race indicator (1 = Asian, 0 = otherwise)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Ethnicity subgroups compared post hoc (Gao 2026 Figure S7); 'There is no remarkable difference in predicted exposures of Trastuzumab rezetecan between ethnicity groups'. Not retained; no point estimate published."
    ),
    CRCL = list(
      description = "Creatinine clearance (renal function category driver)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Renal function categories compared post hoc (Gao 2026 Figure S1); no remarkable difference in exposure. Not retained; no point estimate published."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (NCI-ODWG)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Hepatic function categories compared post hoc (Gao 2026 Figure S2); no remarkable difference in exposure. Not retained; no point estimate published. AST was retained as a continuous covariate on payload clearance instead."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 645,
    n_studies      = 3,
    n_observations = "18,671 concentration records (9421 intact ADC, 9250 released payload) from 26,988 total records including 8317 dosing records",
    age_median     = "56 years (typical-patient value used as the AGE covariate reference)",
    weight_median  = "60.6 kg (typical-patient value used as the WT covariate reference); 5th percentile 45 kg, 95th percentile 82 kg",
    disease_state  = "HER2-expressing or HER2-mutated advanced solid tumors: breast cancer (approximately 60% of the analysis population), gastric or gastroesophageal junction adenocarcinoma, colorectal cancer, and non-small cell lung cancer with HER2 expression, amplification or mutation.",
    dose_range     = "1.0-8.0 mg/kg IV every 3 weeks (Q3W), 21-day treatment cycles.",
    studies        = "Three phase 1 studies: SHR-A1811-I-101 (HER2-expressing or mutated advanced solid tumors), SHR-A1811-I-102 (HER2-expressing advanced gastric or gastroesophageal junction adenocarcinoma and colorectal cancer), and SHR-A1811-I-103 (advanced NSCLC with HER2 expression, amplification or mutation).",
    baseline_albumin_median = "42.5 g/L (typical-patient value used as the ALB covariate reference)",
    baseline_ast_median     = "25 U/L (typical-patient value used as the AST covariate reference)",
    baseline_tumor_size_median = "52 mm sum of target-lesion diameters (typical-patient value used as the TUMSZ covariate reference); 5th percentile 15 mm, 95th percentile 149 mm",
    regions        = "China (Jiangsu Hengrui Pharmaceuticals phase 1 programme)",
    notes          = "Trastuzumab rezetecan (SHR-A1811) is a third-generation HER2-targeting ADC: anti-HER2 antibody trastuzumab, an enzyme-cleavable linker with a chiral cyclopropyl stabilising group, and the topoisomerase-I inhibitor payload rezetecan, at a drug-to-antibody ratio of approximately 6.0. Estimation used FOCE-I in NONMEM 7.5.1 with a SEQUENTIAL two-step approach: the intact-ADC fixed- and random-effect parameters were estimated first and then FIXED while the released-payload parameters were estimated. A one-step joint fit gave close estimates but ran about 6-fold longer. Gao 2026 Section 4 notes a trend toward nonlinear (TMDD-like) elimination at the 1.0 and 2.0 mg/kg dose levels (six patients each); a nonlinear component accounted for approximately 5% of total elimination and did not improve the fit, so linear clearance only was selected.",
    dosing_note    = "Dose the `central` compartment only (IV infusion; Gao 2026 Figure 2). The released-payload compartment is driven by the intact-ADC central compartment and must NOT be dosed. Supply CYCLE as a time-varying covariate column that starts at 1 and increments every 21 days."
  )

  ini({
    # ============================================================
    # Intact ADC (trastuzumab rezetecan) -- Gao 2026 Table 1 and the
    # Section 3.2 equations. Two-compartment, linear elimination,
    # IV infusion (Gao 2026 Figure 2). No absorption parameter is
    # estimated; see the vignette Errata for the paper's stray
    # 'first-order absorption' wording in Section 3.2.
    # ============================================================
    lcl <- log(0.360); label("Intact ADC clearance (L/day)")                    # Gao 2026 Table 1: theta CL = 0.360 L/day (RSE 1.20%); Section 3.2 equation CL = 0.36 * ...
    lvc <- log(2.86);  label("Intact ADC central volume V1 (L)")                # Gao 2026 Table 1: theta V1 = 2.86 L (RSE 0.829%)
    lq  <- log(0.162); label("Intact ADC intercompartmental clearance Q (L/day)") # Gao 2026 Table 1: theta Q = 0.162 L/day (RSE 4.57%)
    lvp <- log(2.88);  label("Intact ADC peripheral volume V2 (L)")             # Gao 2026 Table 1: theta V2 = 2.88 L (RSE 4.34%)

    # Covariate effects on the intact ADC -- Gao 2026 Section 3.2 equations:
    #   CL = 0.36 * (BW/60.6)^0.525 * (SOD_B/52)^0.0686 * exp(eta_CL)
    #   V1 = 2.86 * (BW/60.6)^0.544 * (AGE/56)^0.198 * (ALB/42.5)^-0.363 * exp(eta_V1)
    #   Q  = 0.162 * exp(eta_Q)
    #   V2 = 2.88 * (ALB/42.5)^-1.44 * exp(eta_V2)
    e_wt_cl    <-  0.525;  label("Power exponent of WT on intact-ADC CL (unitless)")   # Gao 2026 Table 1: theta CL_BW = 0.525 (RSE 12.3%)
    e_tumsz_cl <-  0.0686; label("Power exponent of TUMSZ on intact-ADC CL (unitless)") # Gao 2026 Table 1: theta CL_SOD_B = 0.0686 (RSE 25.1%)
    e_wt_vc    <-  0.544;  label("Power exponent of WT on intact-ADC V1 (unitless)")   # Gao 2026 Table 1: theta V1_BW = 0.544 (RSE 8.05%)
    e_age_vc   <-  0.198;  label("Power exponent of AGE on intact-ADC V1 (unitless)")  # Gao 2026 Table 1: theta V1_AGE = 0.198 (RSE 17.9%)
    e_alb_vc   <- -0.363;  label("Power exponent of ALB on intact-ADC V1 (unitless)")  # Gao 2026 Table 1: theta V1_ALB = -0.363 (RSE 23.0%)
    e_alb_vp   <- -1.44;   label("Power exponent of ALB on intact-ADC V2 (unitless)")  # Gao 2026 Table 1: theta V2_ALB = -1.44 (RSE 31.6%)

    # ============================================================
    # Released payload (rezetecan) -- Gao 2026 Table 1 and the
    # Section 3.3 equations. One-compartment, first-order release
    # from the intact ADC, linear elimination.
    # ============================================================
    lkrel    <- log(0.814); label("Cycle-1 payload release-rate constant RAT from the intact ADC (1/day)") # Gao 2026 Table 1: theta RAT = 0.814 /day (RSE 3.12%)
    lfactor1 <- log(0.693); label("Multiplicative scaling ALPHA applied to the release rate from Cycle 2 onward (unitless)") # Gao 2026 Table 1: theta ALPHA = 0.693 (RSE 1.35%); Section 3.3 'Krel = RAT ( x ALPHA, if Cycle > 1)'
    lcl_rez  <- log(392);   label("Released-payload clearance CL3 (L/day)")     # Gao 2026 Table 1: theta CL3 = 392 L/day (RSE 2.46%)
    lvc_rez  <- fixed(log(30)); label("Released-payload volume of distribution V3 (L)")  # Gao 2026 Table 1: theta V3 = 30.0 Fix. Section 3.3: fixed 'similar to previously reported value of exatecan mesylate, to avoid identifiability issues'.

    # Covariate effects on the released payload -- Gao 2026 Section 3.3:
    #   NSCLC:    RAT = 0.814            * (BW/60.6)^-0.546 * (SOD_B/52)^0.1 * exp(eta_RAT)
    #   BC:       RAT = (0.814 + 0.0562) * (BW/60.6)^-0.546 * (SOD_B/52)^0.1 * exp(eta_RAT)
    #   GC/GEJ:   RAT = (0.814 - 0.129)  * (BW/60.6)^-0.546 * (SOD_B/52)^0.1 * exp(eta_RAT)
    #   Others:   RAT = (0.814 + 0.0464) * (BW/60.6)^-0.546 * (SOD_B/52)^0.1 * exp(eta_RAT)
    #   V3  = 30  * (AGE/56)^0.420  * exp(eta_V3)
    #   CL3 = 392 * (AST/25)^-0.167 * exp(eta_CL3)
    e_wt_krel    <- -0.546; label("Power exponent of WT on the payload release-rate constant (unitless)")    # Gao 2026 Table 1: theta RAT_BW = -0.546 (RSE 14.8%)
    e_tumsz_krel <-  0.100; label("Power exponent of TUMSZ on the payload release-rate constant (unitless)") # Gao 2026 Table 1: theta RAT_SOD_B = 0.100 (RSE 23.0%)

    # Cancer-type effects are ADDITIVE on the release-rate constant
    # (units 1/day) and are applied BEFORE the WT and TUMSZ power terms.
    # NSCLC is the reference (all three indicators zero).
    e_tumtp_breast_krel  <-  0.0562; label("Additive breast-cancer shift on the payload release-rate constant (1/day; vs NSCLC reference)")        # Gao 2026 Table 1: theta RAT_CT1 = 0.0562 /day (RSE 51.6%)
    e_tumtp_gastric_krel <- -0.129;  label("Additive gastric/GEJ-cancer shift on the payload release-rate constant (1/day; vs NSCLC reference)")   # Gao 2026 Table 1: theta RAT_CT2 = -0.129 /day (RSE 29.2%)
    e_tumtp_other_krel   <-  0.0464; label("Additive other-tumor-type shift on the payload release-rate constant (1/day; vs NSCLC reference)")     # Gao 2026 Table 1: theta RAT_CT3 = 0.0464 /day (RSE 78.0%)

    e_age_vc_rez <-  0.420; label("Power exponent of AGE on released-payload V3 (unitless)")  # Gao 2026 Table 1: theta V3_AGE = 0.420 (RSE 23.5%)
    e_ast_cl_rez <- -0.167; label("Power exponent of AST on released-payload CL3 (unitless)") # Gao 2026 Table 1: theta CL3_AST = -0.167 (RSE 26.9%)

    # ============================================================
    # Inter-individual variability. Gao 2026 Table 1 reports these as
    # omega^2 rows (variances on the log scale); the IIV was
    # 'described by an exponential model' (Section 2.3) and
    # 'assumed to have a lognormal distribution' (Section 3.2).
    # ============================================================
    etalcl ~ 0.0591                # Gao 2026 Table 1 row 'omega^2 CL' = 0.0591 (RSE 8.70%, shrinkage 8.60%)
    etalvc ~ 0.0408                # Gao 2026 Table 1 row 'omega^2 V1' = 0.0408 (RSE 23.5%, shrinkage 6.20%)
    etalq ~ 0.0798                 # Gao 2026 Table 1 row 'omega^2 Q' = 0.0798 (RSE 33.6%, shrinkage 56.0%)
    etalvp ~ 0.546                 # Gao 2026 Table 1 row 'omega^2 V2' = 0.546 (RSE 8.90%, shrinkage 21.6%)
    etalkrel ~ 0.0419              # Gao 2026 Table 1 row 'omega^2 RAT' = 0.0419 (RSE 25.3%, shrinkage 44.3%)
    etalfactor1 ~ 0.0692           # Gao 2026 Table 1 row 'omega^2 ALPHA' = 0.0692 (RSE 14.6%, shrinkage 20.0%)
    etalcl_rez ~ 0.140             # Gao 2026 Table 1 row 'omega^2 CL3' = 0.140 (RSE 14.5%, shrinkage 19.7%)
    etalvc_rez ~ 0.146             # Gao 2026 Table 1 row 'omega^2 V3' = 0.146 (RSE 13.3%, shrinkage 21.6%)

    # ============================================================
    # Residual variability. Gao 2026 Section 3.2: the intact-ADC
    # residual 'was best explained by a combined proportional and
    # additive error'; Table 1 note: 'residual variability (RUV),
    # included additive and proportional error terms for
    # Trastuzumab rezetecan and a proportional error term only for
    # released toxin.'
    #
    # The Table 1 RUV rows are NONMEM $SIGMA VARIANCES, on the same
    # 'Typical value' column as the omega^2 rows, so the SDs below
    # are square roots. Two independent checks support the variance
    # reading: (1) as SDs the proportional terms would be 3.13% and
    # 7.93% CV, implausibly tight for clinical bioanalytical assays
    # and irreconcilable with the spread in Gao 2026 Figure 1;
    # (2) as a variance the additive SD is sqrt(2.14) = 1.46 ug/mL,
    # about 1.5x the stated intact-ADC assay LLOQ of 1.00 ug/mL,
    # which is the expected magnitude.
    # ============================================================
    propSd     <- 0.176918; label("Intact-ADC proportional residual SD (fraction); sqrt(0.0313)")   # Gao 2026 Table 1 row 'Trastuzumab Rezetecan Prop RUV' = 0.0313 (variance; RSE 15.2%, shrinkage 6.90%)
    addSd      <- 1.462874; label("Intact-ADC additive residual SD (ug/mL); sqrt(2.14)")            # Gao 2026 Table 1 row 'Trastuzumab Rezetecan Add RUV' = 2.14 (variance; RSE 20.4%, shrinkage 6.90%)
    propSd_rez <- 0.281603; label("Released-payload proportional residual SD (fraction); sqrt(0.0793)") # Gao 2026 Table 1 row 'Toxin Rezetecan Prop RUV' = 0.0793 (variance; RSE 3.50%, shrinkage 8.50%)
  })
  model({
    # ============================================================
    # Individual parameters -- intact ADC (Gao 2026 Section 3.2)
    # ============================================================
    cl <- exp(lcl + etalcl) * (WT / 60.6)^e_wt_cl * (TUMSZ / 52)^e_tumsz_cl
    vc <- exp(lvc + etalvc) * (WT / 60.6)^e_wt_vc * (AGE / 56)^e_age_vc *
      (ALB / 42.5)^e_alb_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (ALB / 42.5)^e_alb_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ============================================================
    # Individual parameters -- released payload (Gao 2026 Section 3.3)
    # ============================================================
    # Cancer-type shifts are additive on the release-rate constant and
    # precede the power terms, exactly as the four printed equations
    # show. NSCLC is the reference: with all three indicators zero the
    # typical release rate is the bare 0.814 /day.
    krel_typ <- exp(lkrel) +
      e_tumtp_breast_krel  * TUMTP_BREAST +
      e_tumtp_gastric_krel * TUMTP_GASTRIC +
      e_tumtp_other_krel   * TUMTP_OTHER

    # ALPHA applies from Cycle 2 onward; Cycle 1 uses the unscaled rate.
    factor1 <- exp(lfactor1 + etalfactor1)
    factor_krel <- factor1
    if (CYCLE < 2) factor_krel <- 1.0

    krel <- krel_typ * (WT / 60.6)^e_wt_krel * (TUMSZ / 52)^e_tumsz_krel *
      exp(etalkrel) * factor_krel

    cl_rez <- exp(lcl_rez + etalcl_rez) * (AST / 25)^e_ast_cl_rez
    vc_rez <- exp(lvc_rez + etalvc_rez) * (AGE / 56)^e_age_vc_rez

    kel_rez <- cl_rez / vc_rez

    # ------------------------------------------------------------
    # Molar-mass ratio for the payload formation term.
    #
    # Gao 2026 Section 3.3 states that 'the time course of intact
    # trastuzumab rezetecan concentrations, adjusted for the molar
    # mass, was the input to the released-payload model', but the
    # paper reports NO molecular weight for either the ADC or the
    # payload, and no drug-to-antibody-ratio multiplier appears in
    # the Figure 2 schematic or in any printed equation. `mwr` is
    # therefore held at 1, which makes `central_rez` an
    # ADC-molar-equivalent amount: the payload profile's SHAPE,
    # TIMING and every covariate RATIO are exact, while its absolute
    # mass concentration carries the unreported factor
    # MW_rezetecan / MW_ADC. Set `mwr` to that ratio to obtain mass
    # units. Nothing validated in the vignette depends on `mwr`,
    # because the payload residual error is purely proportional and
    # every published payload target is a ratio. This follows the
    # library precedent for an unreported payload scale constant in
    # `Lu_2022_patritumab.R` (V_DXd fixed to 1 L). See the vignette
    # 'Assumptions and deviations'.
    # ------------------------------------------------------------
    mwr <- 1.0

    # ============================================================
    # ODE system (Gao 2026 Figure 2). The ADC is dosed by IV
    # infusion into `central`. The Krel arrow into the payload
    # compartment is drawn DASHED in Figure 2 and the paper fitted
    # the two analytes sequentially with the ADC parameters fixed,
    # so payload formation does NOT deplete the ADC: the intact-ADC
    # disposition is exactly the two-compartment system above,
    # independent of Krel. This is the same forcing-function
    # structure used by Sathe_2024_sacituzumab.R.
    # ============================================================
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_rez) <-  krel * central * mwr - kel_rez * central_rez

    # ============================================================
    # Observations
    # ============================================================
    # Intact ADC: `central` in mg and `vc` in L give Cc in mg/L =
    # ug/mL, the unit Gao 2026 uses for the intact ADC throughout
    # (Figure 1A/C/E; assay LLOQ 1.00 ug/mL).
    Cc <- central / vc
    # Released payload in ADC-molar-equivalent ug/mL; multiply by
    # mwr = MW_rezetecan / MW_ADC and by 1000 to compare against the
    # ng/mL scale of Gao 2026 Figure 1B/D/F.
    Cc_rez <- central_rez / vc_rez

    Cc     ~ add(addSd) + prop(propSd)
    Cc_rez ~ prop(propSd_rez)
  })
}
