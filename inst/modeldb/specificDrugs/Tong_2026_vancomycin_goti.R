Tong_2026_vancomycin_goti <- function() {
  description <- paste(
    "Modified Goti (\"Tong\") two-compartment IV population PK model for vancomycin as implemented for",
    "model-informed precision dosing (MIPD) in the InsightRX Nova clinical decision support software,",
    "and used as the PRE-intervention default model for patients with BMI < 40 kg/m2 in Tong 2026.",
    "Cockcroft-Gault creatinine clearance is computed inside the model from AGE, WT, SEXF and CREAT and",
    "capped at 150 mL/min; CL scales as (CRCL/120)^0.8, central volume as (WT/70), and intermittent",
    "hemodialysis acts as a multiplicative factor on CL (0.7) and central volume (0.5). An additive",
    "hemodialysis clearance term is carried structurally but is identically zero in the Tong 2026 cohort,",
    "which excluded dialysis patients. All population parameters are FIXED priors for MAP Bayesian",
    "estimation; none were estimated in Tong 2026.",
    sep = " "
  )
  reference <- paste(
    "Tong DMH, Brooks JT, Keizer RJ, Hughes JH. Vancomycin target attainment improved following",
    "population pharmacokinetic model switch: a large-scale quasi-experimental study of precision",
    "dosing. JAC Antimicrob Resist. 2026. doi:10.1093/jacamr/dlag016 (Supplementary data, Code section,",
    "\"Modified Goti (Tong) model\" NONMEM control stream; Table S1).",
    "Structural model and parameter estimates originate from Goti V, Chaturvedula A, Fossler MJ et al.",
    "Ther Drug Monit 2018;40:212-221. doi:10.1097/FTD.0000000000000490; the age-adjusted-creatinine",
    "modification is described in Tong DMH, Hughes JH, Keizer RJ. Ther Drug Monit 2021;43:139-140.",
    "doi:10.1097/FTD.0000000000000819. See also modellib('Goti_2018_vancomycin') for the extraction",
    "made directly from the Goti 2018 publication, which differs in three documented respects",
    "(see vignette Errata).",
    sep = " "
  )
  vignette <- "Tong_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI < 40 kg/m2 cohort: median 79.0 kg (range 20.4-173.8).",
        "Enters twice in the supplement control stream: as the numerator of the internal",
        "Cockcroft-Gault creatinine-clearance calculation (CRCLi = (140 - AGE) * WT * 0.85^SEXF /",
        "(72 * CREAT)) and as the linear size scalar on central volume (TVV = 58.4 * (WT/70)).",
        "NOTE: the peripheral volume is NOT weight-scaled in this stream (TVV2 = THETA(3) = 38.4,",
        "with no (WT/70) factor), which Table S1 confirms as 'V2 = 38.4'.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI < 40 kg/m2 cohort: median 65.8 years (range 18.0 to 90+; ages above",
        "90 are aggregated for de-identification). Used only inside the internal Cockcroft-Gault",
        "creatinine-clearance calculation.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tong 2026 Table 1, BMI < 40 kg/m2 cohort: 39961 male / 55030 female treatment courses",
        "(57.9% female). The supplement control stream uses the OPPOSITE polarity, SEX with 1 = male,",
        "written as 0.85**(1-SEX) in the Cockcroft-Gault term so that females receive the 0.85 factor.",
        "Converted to the canonical SEXF (1 = female) via SEXF = 1 - SEX, so the term becomes",
        "0.85^SEXF. Verified by the sex-specific fat-free-mass branch of the sibling Hughes stream",
        "(IF(SEX.EQ.0) selects the FEMALE Janmahasatian coefficients), which independently pins",
        "SEX = 1 as male.",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tong 2026 Table 1, BMI < 40 kg/m2 cohort: median 0.90 mg/dL (range 0.05-25.3). The",
        "supplement control stream assigns CRcalc = CR and uses it directly as the denominator of the",
        "Cockcroft-Gault equation with the 72 constant, which fixes the unit as mg/dL. NOTE: unlike",
        "the underlying Goti 2018 publication and unlike Tong 2021 (reference 19, whose age-adjusted",
        "creatinine floor DECREASED predictive performance in the elderly and was therefore NOT",
        "adopted), this stream applies NO lower truncation of serum creatinine in elderly subjects --",
        "the only truncation present is the 150 mL/min cap on the resulting creatinine clearance.",
        sep = " "
      ),
      source_name        = "CR"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = on intermittent hemodialysis, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no intermittent hemodialysis)",
      notes              = paste(
        "Enters the supplement control stream as the exponent DIAL in TVCL = ... * TH_DIAL_CL**DIAL",
        "and TVV = ... * TH_DIAL_V**DIAL, carrying the Goti 2018 dialysis factors 0.7 on CL and 0.5",
        "on Vc. IDENTICALLY ZERO throughout the Tong 2026 study population: the Methods state",
        "'Patients were excluded if they were undergoing haemodialysis at any point during treatment.'",
        "Retained here because it is part of the implemented model structure and is exercised when the",
        "model is applied outside the Tong 2026 cohort.",
        sep = " "
      ),
      source_name        = "DIAL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 87586L,
    n_studies      = 19L,
    age_range      = "18.0 to 90+ years",
    age_median     = "65.8 years",
    weight_range   = "20.4-173.8 kg",
    weight_median  = "79.0 kg",
    sex_female_pct = 57.9,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Hospitalized adults (>= 18 years) with BMI < 40 kg/m2 receiving intravenous vancomycin under",
      "routine model-informed precision dosing; at least two doses and at least one measured",
      "concentration required. Patients undergoing haemodialysis at any point during treatment were",
      "excluded, as were patients dosed with a model other than their site's default.",
      sep = " "
    ),
    dose_range     = "Intravenous vancomycin per routine clinical practice; initial doses selected a priori, subsequent doses adapted by MAP Bayesian posterior estimates",
    regions        = "United States (19 hospital systems, patients beginning treatment August 2022 to December 2024)",
    renal_function = "Serum creatinine median 0.90 mg/dL (range 0.05-25.3); haemodialysis patients excluded",
    n_concentrations = 192013L,
    notes          = paste(
      "APPLICATION population from Tong 2026 Table 1 (BMI < 40 kg/m2 cohort: 87586 patients, 94991",
      "treatment courses, 192013 samples), i.e. the cohort this model was USED to dose as the",
      "pre-intervention default -- NOT the cohort it was estimated from. The DEVELOPMENT population",
      "is Goti 2018: 1812 hospitalized adults with 2765 concentrations (Tong 2026 Table S1 rows",
      "'# patients in development population' = 1812 and '# of samples' = 2765). Every population",
      "parameter is a FIXED prior inherited from Goti 2018; Tong 2026 estimated none of them.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # ALL population parameters are FIXED priors for MAP Bayesian estimation.
    # Source: Tong 2026 Supplementary data, Code section, "Modified Goti (Tong)
    # model" control stream, which marks every $THETA and the $OMEGA BLOCK(3)
    # with the NONMEM FIX flag. Typical subject is a 70-kg non-dialysis patient
    # with Cockcroft-Gault CrCL = 120 mL/min.
    # ------------------------------------------------------------------------
    lcl <- fixed(log(4.5));  label("Clearance at CRCL = 120 mL/min, non-dialysis (L/h)")
    # $THETA(1) 4.5 FIX; Table S1 "CL = 4.5 * (CRCLi/120)^0.8"
    lvc <- fixed(log(58.4)); label("Central volume at WT = 70 kg, non-dialysis (L)")
    # $THETA(2) 58.4 FIX; Table S1 "V = 58.4 * WT/70"
    lvp <- fixed(log(38.4)); label("Peripheral volume (L; NOT weight-scaled)")
    # $THETA(3) 38.4 FIX, entered as TVV2 = THETA(3) with no (WT/70) factor; Table S1 "V2 = 38.4"
    lq  <- fixed(log(6.5));  label("Intercompartmental clearance Q (L/h)")
    # $THETA(4) 6.5 FIX; Table S1 "Q = 6.5"

    # Covariate effects, from the same control stream:
    #   TVCL = THETA(1) * (CRCLi/120)^THETA(5) * THETA(6)^DIAL
    #   TVV  = THETA(2) * (WT/70)              * THETA(7)^DIAL
    e_crcl_cl     <- fixed(0.8); label("Power exponent on (CRCL/120) for CL (unitless)")
    # $THETA(5) 0.8 FIX; Table S1 exponent in "(CRCLi/120)^0.8"
    e_hemodial_cl <- fixed(0.7); label("Multiplicative factor on CL when RRT_HEMODIAL_STATUS = 1 (unitless)")
    # $THETA(6) 0.7 FIX (TH_DIAL_CL)
    e_hemodial_vc <- fixed(0.5); label("Multiplicative factor on Vc when RRT_HEMODIAL_STATUS = 1 (unitless)")
    # $THETA(7) 0.5 FIX (TH_DIAL_V)

    # Internal Cockcroft-Gault cap, from IF(CRCLi.GT.150.0) CRCLi = 150.0 in the
    # control stream; Table S1 states it as "If CRCL > 150, then CRCLi = 150".
    crcl_cap <- fixed(150); label("Upper cap applied to the internally computed Cockcroft-Gault CrCL (mL/min)")

    # Additive hemodialysis clearance, from CLTOT = CL + CL_HEMO in the control
    # stream. Supplied per-patient by the MIPD software; identically zero in the
    # Tong 2026 cohort, which excluded dialysis patients (Methods, Patient data
    # collection). Encoded as fixed(0) so the structural term stays visible.
    cl_hemo <- fixed(0); label("Additive hemodialysis clearance (L/h; 0 in the Tong 2026 cohort)")

    # Inter-individual variability. $OMEGA BLOCK(3) FIX with all off-diagonal
    # elements 0, so the three etas are declared independently. The stream stores
    # variances on the omega^2 = CV^2 convention (NOT log(CV^2 + 1)): each
    # sqrt(variance) reproduces the Table S1 %CV row exactly.
    # NOTE: these source traces are on their own lines rather than trailing the
    # eta declarations. rxode2 rewrites a trailing comment on an `eta ~ ...` line
    # into a label() call, and a comment containing double quotes then fails to
    # parse.
    # $OMEGA(1,1) 0.1584 FIX; sqrt = 0.398 -> Table S1 "IIV on CL (%CV) 39.8"
    etalcl ~ fixed(0.1584)
    # $OMEGA(2,2) 0.6659 FIX; sqrt = 0.816 -> Table S1 "IIV on V (%CV) 81.6"
    etalvc ~ fixed(0.6659)
    # $OMEGA(3,3) 0.326 FIX; sqrt = 0.571 -> Table S1 "IIV on V2 (%CV) 57.1"
    etalvp ~ fixed(0.326)
    # No eta on Q: the stream sets Q = TVQ with no EXP(ETA), and Table S1 reports
    # "IIV on Q (%CV)" as "-" for this model.

    # Combined additive + proportional residual error, from the $ERROR block
    # W = SQRT(IPRED**2 * PROP**2 + ADD**2) with PROP and ADD hardcoded.
    addSd  <- fixed(1.42);  label("Additive residual error (mg/L)")
    # $ERROR ADD = 1.42; Table S1 "Additive error (mg/L) 1.42"
    propSd <- fixed(0.197); label("Proportional residual error (fraction)")
    # $ERROR PROP = 0.197; Table S1 "Proportional error (%) 19.7"
  })
  model({
    # 1. Derived covariate terms. Cockcroft-Gault creatinine clearance is
    #    computed INSIDE the model (the stream does not accept CRCL as a data
    #    item for this variant) and then capped:
    #      CRCLi = ((140 - AGE) * WT * 0.85^(1-SEX)) / (72 * CR)
    #      IF(CRCLi > 150) CRCLi = 150
    #    with SEXF = 1 - SEX so that 0.85^(1-SEX) becomes 0.85^SEXF.
    crcl_cg     <- ((140 - AGE) * WT * 0.85^SEXF) / (72 * CREAT)
    crcl_capped <- min(crcl_cg, crcl_cap)

    # 2. Individual PK parameters.
    cl <- exp(lcl + etalcl) * (crcl_capped / 120)^e_crcl_cl * e_hemodial_cl^RRT_HEMODIAL_STATUS
    vc <- exp(lvc + etalvc) * (WT / 70) * e_hemodial_vc^RRT_HEMODIAL_STATUS
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # Total elimination clearance, CLTOT = CL + CL_HEMO.
    cltot <- cl + cl_hemo

    # 3. Micro-constants.
    kel <- cltot / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system, transcribed from $DES:
    #      DADT(1) = -(CLTOT/V)*A(1) + (Q/V2)*A(2) - (Q/V)*A(1)
    #      DADT(2) =                 - (Q/V2)*A(2) + (Q/V)*A(1)
    #    The stream's third state, DADT(3) = A(1)/V, is a quadrature-only AUC
    #    accumulator used for MIPD reporting; it has no effect on the PK and is
    #    omitted here (AUC is obtained by NCA in the vignette instead).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation and error. S1 = V in the stream, so dose in mg over volume
    #    in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
