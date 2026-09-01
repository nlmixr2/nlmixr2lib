Tong_2026_vancomycin_thomson <- function() {
  description <- paste(
    "Capped Thomson two-compartment IV population PK model for vancomycin as implemented for",
    "model-informed precision dosing (MIPD) in the InsightRX Nova clinical decision support software,",
    "and used as the POST-intervention default model for patients with BMI < 40 kg/m2 in Tong 2026.",
    "Clearance is a linear (slope-intercept) function of creatinine clearance centred at 66 mL/min,",
    "with creatinine clearance capped at 150 mL/min; both volumes scale linearly with body weight.",
    "An additive hemodialysis clearance term is carried structurally but is identically zero in the",
    "Tong 2026 cohort, which excluded dialysis patients. All population parameters are FIXED priors",
    "for MAP Bayesian estimation; none were estimated in Tong 2026. Switching to this model from the",
    "modified Goti model improved AUC target attainment on day 1 by 4.6 percentage points.",
    sep = " "
  )
  reference <- paste(
    "Tong DMH, Brooks JT, Keizer RJ, Hughes JH. Vancomycin target attainment improved following",
    "population pharmacokinetic model switch: a large-scale quasi-experimental study of precision",
    "dosing. JAC Antimicrob Resist. 2026. doi:10.1093/jacamr/dlag016 (Supplementary data, Code",
    "section, \"Modified Thomson model\" NONMEM control stream; Table S1).",
    "Structural model and parameter estimates originate from Thomson AH, Staatz CE, Tobin CM et al.",
    "J Antimicrob Chemother 2009;63:1050-1057. doi:10.1093/jac/dkp085; the creatinine-clearance cap",
    "was evaluated in Hughes MSA, Lee T, Faldasz JD et al. Pharmacotherapy 2024;44:794-802.",
    "doi:10.1002/phar.4613.",
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
        "Tong 2026 Table 1, BMI < 40 kg/m2 cohort: median 79.0 kg (range 20.4-173.8). Scales BOTH",
        "volumes linearly and without normalisation to a reference weight: TVV = 0.675 * WT and",
        "TVV2 = 0.732 * WT, so the coefficients are volumes per kilogram (L/kg) rather than volumes",
        "at a reference weight. Table S1 confirms 'V = 0.675 * WT' and 'V2 = 0.732 * WT'.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Supplied as a data item to this variant rather than computed inside the model (unlike the",
        "sibling modified-Goti and Hughes streams, which derive it internally). The InsightRX pipeline",
        "supplies it in L/h, and the control stream converts to mL/min on its first line,",
        "EGFR = CRCL * 16.66667 (1 L/h = 16.667 mL/min); the misleading variable name EGFR denotes",
        "the unit-converted creatinine clearance, not a BSA-normalized eGFR. This model file takes",
        "CRCL directly in mL/min, which is the scale the 66 mL/min centring value and the 150 mL/min",
        "cap are expressed on, so no conversion is applied here. Stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md, which accepts raw un-normalized mL/min when",
        "the source does not BSA-normalize (same convention as Goti_2018_vancomycin.R and",
        "Delattre_2010_amikacin.R). Tong 2026 Table 1 reports serum creatinine (median 0.90 mg/dL)",
        "rather than creatinine clearance for this cohort.",
        sep = " "
      ),
      source_name        = "CRCL"
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
      "post-intervention default -- NOT the cohort it was estimated from. The DEVELOPMENT population",
      "is Thomson 2009: 398 patients with 1557 concentrations (Tong 2026 Table S1 rows '# patients in",
      "development population' = 398 and '# of samples' = 1557). Every population parameter is a",
      "FIXED prior inherited from Thomson 2009; Tong 2026 estimated none of them.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # ALL population parameters are FIXED priors for MAP Bayesian estimation.
    # Source: Tong 2026 Supplementary data, Code section, "Modified Thomson
    # model" control stream, which marks every $THETA, both $OMEGA records and
    # $SIGMA with the NONMEM FIX flag.
    # ------------------------------------------------------------------------
    lcl <- fixed(log(2.99));  label("Clearance at CRCL = 66 mL/min (L/h)")
    # $THETA(1) 2.99 FIX; Table S1 "CL = 2.99 * (1 + 0.0154 * (CRCLi-66))"
    lvc <- fixed(log(0.675)); label("Central volume per kg body weight (L/kg)")
    # $THETA(2) 0.675 FIX, entered as TVV = THETA(2) * WT; Table S1 "V = 0.675 * WT"
    lq  <- fixed(log(2.28));  label("Intercompartmental clearance Q (L/h)")
    # $THETA(4) 2.28 FIX; Table S1 "Q = 2.28"
    lvp <- fixed(log(0.732)); label("Peripheral volume per kg body weight (L/kg)")
    # $THETA(5) 0.732 FIX, entered as TVV2 = THETA(5) * WT; Table S1 "V2 = 0.732 * WT"

    # Renal-function effect on CL. Unlike the power form used by the sibling
    # modified-Goti and Hughes models, Thomson parameterises CL as a LINEAR
    # function of creatinine clearance centred at 66 mL/min:
    #   TVCL = THETA(1) * (1 + THETA(3) * (CRCLcapped - 66))
    # so the coefficient is a fractional change in CL per mL/min of CrCL.
    e_crcl_cl <- fixed(0.0154); label("Fractional change in CL per mL/min of CRCL above 66 mL/min (1/(mL/min))")
    # $THETA(3) 0.0154 FIX (TH_CRCL); Table S1 slope in "(1 + 0.0154 * (CRCLi-66))"
    crcl_ref  <- fixed(66);     label("Centring value for the linear CRCL effect on CL (mL/min)")
    # Hardcoded 66 in TVCL = THETA(1) * (1 + TH_CRCL * (CRCLcapped - 66))

    # Creatinine-clearance cap, from IF(EGFR.GT.CRCLCAP) CRCLcapped = CRCLCAP in
    # the control stream. CRCLCAP itself is not given a numeric value in the
    # supplement listing; Table S1 supplies it as "If CRCL > 150, then CRCLi = 150"
    # -- this cap is what the paper's "capped Thomson" name refers to.
    crcl_cap <- fixed(150); label("Upper cap applied to CRCL before the CL effect (mL/min)")

    # Additive hemodialysis clearance, from CLTOT = CL + CL_HEMO in the control
    # stream. Supplied per-patient by the MIPD software; identically zero in the
    # Tong 2026 cohort, which excluded dialysis patients (Methods, Patient data
    # collection). Encoded as fixed(0) so the structural term stays visible.
    cl_hemo <- fixed(0); label("Additive hemodialysis clearance (L/h; 0 in the Tong 2026 cohort)")

    # Inter-individual variability. The stream stores variances on the
    # omega^2 = CV^2 convention (NOT log(CV^2 + 1)): each sqrt(variance)
    # reproduces the Table S1 %CV row exactly.
    #   $OMEGA BLOCK(2) FIX   0.0729 / 0.01  0.0225   -> correlated CL and Vc
    #   $OMEGA 0.2401 FIX                              -> Q
    #   $OMEGA 1.69   FIX                              -> Vp
    # NOTE: these source traces are on their own lines rather than trailing the
    # eta declarations. rxode2 rewrites a trailing comment on an `eta ~ ...` line
    # into a label() call, and a comment containing double quotes then fails to
    # parse.
    # sqrt(0.0729) = 0.27 -> Table S1 "IIV on CL (%CV) 27"; sqrt(0.0225) = 0.15 ->
    # Table S1 "IIV on V (%CV) 15"; off-diagonal 0.01 gives correlation 0.247
    etalcl + etalvc ~ fixed(c(
      0.0729,
      0.01, 0.0225
    ))
    # sqrt(0.2401) = 0.49 -> Table S1 "IIV on Q (%CV) 49"
    etalq  ~ fixed(0.2401)
    # sqrt(1.69) = 1.30 -> Table S1 "IIV on V2 (%CV) 130"
    etalvp ~ fixed(1.69)

    # Combined additive + proportional residual error, from the $ERROR block
    # W = SQRT(IPRED**2 * PROP**2 + ADD**2) with PROP and ADD hardcoded.
    addSd  <- fixed(1.6);  label("Additive residual error (mg/L)")
    # $ERROR ADD = 1.6; Table S1 "Additive error (mg/L) 1.6"
    propSd <- fixed(0.15); label("Proportional residual error (fraction)")
    # $ERROR PROP = 0.15; Table S1 "Proportional error (%) 15"
  })
  model({
    # 1. Derived covariate terms. The control stream converts the incoming CRCL
    #    from L/h to mL/min (EGFR = CRCL * 16.66667) and then caps it. This file
    #    takes CRCL already in mL/min, so only the cap is applied here.
    crcl_capped <- min(CRCL, crcl_cap)

    # 2. Individual PK parameters. Note the linear (not power) CRCL effect and
    #    that both volumes are per-kg coefficients multiplied by WT directly.
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (crcl_capped - crcl_ref))
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * WT

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
