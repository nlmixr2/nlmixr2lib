Hu_2024_ipilimumab <- function() {
  description <- "Two-compartment population PK model with stationary clearance for intravenous ipilimumab (anti-CTLA-4 IgG1) in a pooled adult and pediatric (1-17 years) oncology population, with tumor-type-dependent pediatric effects on clearance beyond body size (Hu 2024)"
  reference <- "Hu Z, Liu S, Zhao Y, Du S, Hamuro L, Shen J, Roy A, Zhu L. Nivolumab and ipilimumab population pharmacokinetics in support of pediatric dose recommendations-Going beyond the body-size effect. CPT Pharmacometrics Syst Pharmacol. 2024;13(3):476-493. doi:10.1002/psp4.13098"
  vignette <- "Hu_2024_nivolumab_ipilimumab_pediatric"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "ipilimumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ipilimumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    LBM = list(
      description        = "Baseline lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 55 kg (exponent 0.789) and on VC with the same reference (exponent 0.874); Q and VP inherit the CL and VC exponents respectively (Hu 2024 Results equations for the ipilimumab full model; Figure 2 caption: 'The effect of lean body mass was also added to intercompartmental CL and peripheral volume, respectively, estimates of which were fixed to be similar to those of CL and VC'). Table 1 Note: LBM was estimated using the Boer equation (adults) and the Peter equation (children). Missing values were imputed to the reference in the source NONMEM code (Appendix File S3).",
      source_name        = "ELBMB"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only through the three age bands that define the paper's patient-population and volume covariates: young pediatric (< 12 years), adolescent (12-17 years), adult (>= 18 years). Appendix File S3 applies the splits as AGE <= 11 and 12 <= AGE <= 17. Both the numeric power effect of age on CL (THETA(9)) and the numeric power effect of age on VC (THETA(20)) were fixed to zero in the selected Full1 model, so age enters this model only categorically.",
      source_name        = "AGE"
    ),
    TUMTP_MEL = list(
      description        = "Tumor-type indicator: melanoma",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 in the reference subject. The Hu 2024 reference patient population is adult melanoma, so TUMTP_MEL = 1 together with AGE >= 18 gives the reference (no population term on CL).",
      notes              = "Part of the three-level patient-population factor (melanoma / CNS tumor / other solid tumor) that Hu 2024 crosses with an adult-versus-pediatric split. Pediatric (< 18 years) melanoma carries an exponential CL effect of exp(-0.347) = 0.707 (29% lower than adult melanoma); adult melanoma is the reference. Table 2 tumor-type counts: MEL 1284, CNST 78, other ST 65 of N = 1427, comprising adult MEL n = 1261, adult CNST n = 6, adult other ST n = 22, pediatric ST n = 43, pediatric CNST n = 72, pediatric MEL n = 23.",
      source_name        = "POPN"
    ),
    TUMTP_CNS_PRIM = list(
      description        = "Tumor-type indicator: primary central nervous system tumor",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma or other solid tumor)",
      notes              = "Exponential effects on CL: adults exp(-0.661) = 0.516 (48% lower than adult melanoma); pediatric subjects (< 18 years) exp(-0.668) = 0.513 (49% lower). The adult CNS-tumor cohort in the ipilimumab dataset is small (n = 6) and young (mean age 19.3 years), which Hu 2024 Discussion offers as the reason the adult and pediatric CNS-tumor effects coincide here but differ in the nivolumab analysis. Pediatric composition (Table 2 footnote d): diffuse intrinsic pontine glioma n = 16, medulloblastoma n = 13, high-grade glioma n = 9, ependymoma n = 8, diffuse midline glioma n = 4, and others.",
      source_name        = "POPN"
    ),
    TUMTP_OTHER = list(
      description        = "Tumor-type indicator: other solid tumors (the residual non-melanoma, non-CNS pool)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (melanoma or CNS tumor)",
      notes              = "Exponential effects on CL: adults exp(-0.698) = 0.498 (50% lower than adult melanoma); pediatric subjects (< 18 years) exp(-0.462) = 0.630 (37% lower). Unlike the nivolumab model, the pediatric solid-tumor effect is NOT split at 12 years -- Appendix File S3 applies a single POPN = 5 term for the whole < 18-year range. Per-paper composition -- adult other (n = 22, Table 2 footnote b): Ewing sarcoma n = 5, osteosarcoma n = 5, rhabdomyosarcoma n = 3, other solid tumors n = 9; pediatric other solid tumors (n = 43, Table 2 footnote c): osteosarcoma n = 8, rhabdomyosarcoma n = 6, Ewing sarcoma n = 5, renal cell carcinoma n = 2, neuroblastoma n = 1, unspecified solid tumor n = 11, others n = 10.",
      source_name        = "POPN"
    ),
    CONMED_NIVO_1Q3W = list(
      description        = "Coadministration regimen: ipilimumab + nivolumab 1 mg/kg every 3 weeks",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ipilimumab monotherapy or any other combination regimen)",
      notes              = "Exponential effect on CL (exp(0.0417) = 1.043; 95% CI -0.00318 to 0.0866). Source column COMBO level 1 (Appendix File S3). Table 2: n = 394 of 1427 received ipilimumab + nivolumab 1 mg/kg.",
      source_name        = "COMBO"
    ),
    CONMED_NIVO_3Q3W = list(
      description        = "Coadministration regimen: ipilimumab + nivolumab 3 mg/kg every 3 weeks",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ipilimumab monotherapy or any other combination regimen)",
      notes              = "Exponential effect on CL (exp(0.316) = 1.372; 95% CI -0.0390 to 0.670, RSE 57.3%). Source column COMBO level 2 (Appendix File S3). Table 2: n = 116 of 1427 received ipilimumab + nivolumab 3 mg/kg.",
      source_name        = "COMBO"
    )
  )

  # Covariate columns present in the Hu 2024 ipilimumab analysis dataset whose
  # coefficients were held at zero in the selected Full1 model and therefore
  # carry no effect here. Documented for provenance only; not referenced in
  # model().
  covariatesDataExcluded <- list(
    CONMED_BUDESONIDE = list(
      description = "Coadministration regimen: ipilimumab 10 mg/kg q3w + oral budesonide 9 mg once daily",
      units       = "(binary)",
      type        = "binary",
      notes       = "COMBO level 3 in the source NONMEM code (Appendix File S3); THETA(17) CL_chemo was declared 0 FIX and does not appear in Table 5. Table 2 records n = 58 of 1427 in this arm."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1427L,
    n_studies      = 10L,
    n_observations = 6020L,
    n_pediatric    = 138L,
    age_range      = "1-89 years (pediatric 1-17 years; adult >= 18 years)",
    age_median     = "58 years overall; pediatric solid tumors 13 y, pediatric CNS tumors 10 y, pediatric melanoma 13 y",
    weight_range   = "10.2-160 kg",
    weight_median  = "78.1 kg overall; pediatric solid tumors 45 kg, pediatric CNS tumors 37.4 kg, pediatric melanoma 56.1 kg",
    sex_female_pct = 36.9,
    race_ethnicity = c(White = 94.8, `Black/African American` = 1.1, Asian = 1.8, Other = 2.2),
    disease_state  = "Advanced melanoma (n = 1284), central nervous system tumors (n = 78), and other solid tumors (n = 65). Analysis populations: adult MEL n = 1261, adult CNST n = 6, adult other ST n = 22, pediatric ST n = 43, pediatric CNST n = 72, pediatric MEL n = 23.",
    dose_range     = "Ipilimumab 0.3-10 mg/kg intravenous infusion q3w for 4 doses (then q12w in several studies), alone or combined with nivolumab 1 or 3 mg/kg q3w, or with oral budesonide 9 mg once daily",
    regions        = "Pooled global phase I, I/II, II, and III studies (10 trials; Table S1)",
    performance_status = "ECOG PS 0 66.6%, PS 1 31.6%, PS 2 1.7%, PS 3 0.1% (screened but not retained in the final ipilimumab model)",
    renal_function = "Baseline eGFR mean 90.6 (SD 24.5) mL/min/1.73 m^2, median 90.1",
    body_composition = "Baseline lean body mass mean 55.1 (SD 12.9) kg, median 56.8 kg; estimated by the Boer equation (adults) and the Peter equation (children)",
    notes          = "Baseline demographics per Hu 2024 Table 2 (N = 1427 across 10 ipilimumab studies, of whom 138 were pediatric). Four phase I/II studies contributed the pediatric data (ADVL1412 / CA209070, CA209908, CA184070, CA184178; Table S1)."
  )

  ini({
    # Structural parameters at the reference covariates: lean body mass 55 kg,
    # adult melanoma on ipilimumab monotherapy (Hu 2024 Results, ipilimumab
    # reference-value list; Table 5 Note). CL and Q are reported in mL/h in
    # Table 5 and are converted to L/day (x 24 / 1000) because this model keeps
    # time in days.
    lcl <- log(13.5 * 24 / 1000); label("Clearance CL_REF at the reference covariates (L/day)")                   # Hu 2024 Table 5: CL_REF = 13.5 mL/h
    lvc <- log(3.90);             label("Central volume VC_REF at the reference covariates (L)")                  # Hu 2024 Table 5: VC_REF = 3.90 L
    lq  <- log(35.8 * 24 / 1000); label("Intercompartmental clearance Q_REF at the reference covariates (L/day)") # Hu 2024 Table 5: Q_REF = 35.8 mL/h
    lvp <- log(3.47);             label("Peripheral volume VP_REF at the reference covariates (L)")               # Hu 2024 Table 5: VP_REF = 3.47 L

    # Continuous covariate effects (power form (LBM / 55)^exponent).
    e_lbm_cl <- 0.789; label("Power exponent of lean body mass on CL, also applied to Q (unitless)")  # Hu 2024 Table 5: CL_LBM = 0.789
    e_lbm_vc <- 0.874; label("Power exponent of lean body mass on VC, also applied to VP (unitless)") # Hu 2024 Table 5: VC_LBM = 0.874

    # Patient-population (tumor type x adult / pediatric) effects on CL,
    # relative to the adult melanoma reference (Hu 2024 Results, ipilimumab
    # CL0TV equation; Appendix File S3 assigns each term by POPN).
    e_cnst_cl     <- -0.661; label("Exponential coefficient of adult CNS tumors on CL (unitless)")                   # Hu 2024 Table 5: CL_CNST = -0.661
    e_oth_cl      <- -0.698; label("Exponential coefficient of adult other solid tumors on CL (unitless)")           # Hu 2024 Table 5: CL_OTH = -0.698
    e_pedoth_cl   <- -0.462; label("Exponential coefficient of pediatric (< 18 y) other solid tumors on CL (unitless)") # Hu 2024 Table 5: CL_PEDOTH = -0.462
    e_pedcnst_cl  <- -0.668; label("Exponential coefficient of pediatric (< 18 y) CNS tumors on CL (unitless)")      # Hu 2024 Table 5: CL_PEDCNST = -0.668
    e_pedmel_cl   <- -0.347; label("Exponential coefficient of pediatric (< 18 y) melanoma on CL (unitless)")        # Hu 2024 Table 5: CL_PEDMEL = -0.347

    # Combination-therapy effects on CL (exponential form).
    e_nivo1q3w_cl <-  0.0417; label("Exponential coefficient of nivolumab 1 mg/kg q3w coadministration on CL (unitless)") # Hu 2024 Table 5: CL_N1Q3 = 0.0417
    e_nivo3q3w_cl <-  0.316;  label("Exponential coefficient of nivolumab 3 mg/kg q3w coadministration on CL (unitless)") # Hu 2024 Table 5: CL_N3Q3 = 0.316

    # Age-band effects on VC (exponential form). Appendix File S3 gates these
    # on AGE alone, independently of tumor type.
    e_ped_vc <- -0.296; label("Exponential coefficient of young pediatric age (< 12 y) on VC (unitless)") # Hu 2024 Table 5: VC_PED = -0.296
    e_ado_vc <- -0.217; label("Exponential coefficient of adolescent age (12-17 y) on VC (unitless)")     # Hu 2024 Table 5: VC_ADO = -0.217

    # Inter-individual variability: a log-normal 2x2 block on CL and VC only.
    # The source control stream declares a single $OMEGA BLOCK(2) and assigns
    # Q and VP without etas (Appendix File S3: "Q = TVQ", "V2 = TVV2"), which
    # is why Table 5 lists no omega^2 for Q or VP.
    etalcl + etalvc ~ c(0.147,
                        0.0258, 0.0531)  # Hu 2024 Table 5: omega^2_CL = 0.147, cov(CL,VC) = 0.0258, omega^2_VC = 0.0531

    # Residual error: combined proportional plus additive. The source control
    # stream fixes $SIGMA to 1 and forms the weight as
    # sqrt((PERR * F)^2 + AERR^2) (Appendix File S3), so both terms are
    # standard deviations on the concentration scale.
    propSd <- 0.185; label("Proportional residual error (fraction)")   # Hu 2024 Table 5: Proportional error term = 0.185
    addSd  <- 1.14;  label("Additive residual error (ug/mL)")          # Hu 2024 Table 5: Additive error term = 1.14 ug/mL
  })
  model({
    # 1. Age-band indicators. Appendix File S3 applies the pediatric splits as
    # AGE <= 11 and 12 <= AGE <= 17; adults are the complement.
    ped   <- 1 - (AGE >= 12)                  # young pediatric, < 12 years
    ado   <- (AGE >= 12) * (1 - (AGE >= 18))  # adolescent, 12-17 years
    adult <- (AGE >= 18)                      # adult, >= 18 years
    child <- 1 - adult                        # pediatric, < 18 years

    # Patient-population effect on CL, relative to adult melanoma. Unlike the
    # companion nivolumab model, every pediatric term here pools the whole
    # < 18-year range (no 12-year split on CL).
    cl_pop <-
      e_cnst_cl    * TUMTP_CNS_PRIM * adult +
      e_oth_cl     * TUMTP_OTHER    * adult +
      e_pedoth_cl  * TUMTP_OTHER    * child +
      e_pedcnst_cl * TUMTP_CNS_PRIM * child +
      e_pedmel_cl  * TUMTP_MEL      * child

    # 2. Individual parameters. Ipilimumab CL was modeled as stationary
    # (Hu 2024 Methods: over a 12-week q3w x 4 dosing period CL was not
    # expected to change materially, so the time-varying term was dropped).
    cl <- exp(lcl + etalcl) *
      (LBM / 55)^e_lbm_cl *
      exp(cl_pop) *
      exp(e_nivo1q3w_cl * CONMED_NIVO_1Q3W) *
      exp(e_nivo3q3w_cl * CONMED_NIVO_3Q3W)

    vc <- exp(lvc + etalvc) *
      (LBM / 55)^e_lbm_vc *
      exp(e_ped_vc * ped) *
      exp(e_ado_vc * ado)

    # Q inherits the lean-body-mass exponent of CL and VP that of VC
    # (Hu 2024 Figure 2 caption). Neither carries an eta in this model.
    q  <- exp(lq)  * (LBM / 55)^e_lbm_cl
    vp <- exp(lvp) * (LBM / 55)^e_lbm_vc

    # 3. Two-compartment micro-constants and ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # 4. Observation. Dose in mg and volumes in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
