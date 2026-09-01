Tong_2026_vancomycin_carreno <- function() {
  description <- paste(
    "Carreno two-compartment IV population PK model for vancomycin in obese adults, translated from",
    "the original non-parametric fit into a parametric form for MAP Bayesian estimation and",
    "implemented for model-informed precision dosing (MIPD) in the InsightRX Nova clinical decision",
    "support software; used as the PRE-intervention default model for patients with BMI >= 40 kg/m2 in",
    "Tong 2026. Clearance is an additive slope-intercept function of creatinine clearance, splitting",
    "into a non-renal intercept and a renal arm proportional to CRCL, each carrying its own",
    "inter-individual variability. No parameter is scaled by body size, the intercompartmental rate",
    "constants k12 and k21 are primary parameters with their own variability, and the IIV distributions",
    "are wide -- Tong 2026 discusses this flexibility as a source of unstable MAP estimates. All",
    "population parameters are FIXED priors; none were estimated in Tong 2026.",
    sep = " "
  )
  reference <- paste(
    "Tong DMH, Brooks JT, Keizer RJ, Hughes JH. Vancomycin target attainment improved following",
    "population pharmacokinetic model switch: a large-scale quasi-experimental study of precision",
    "dosing. JAC Antimicrob Resist. 2026. doi:10.1093/jacamr/dlag016 (Supplementary data, Code",
    "section, \"Carreno model\" NONMEM control stream; Table S2).",
    "Structural model and parameter estimates originate from Carreno JJ, Lomaestro B, Tietjan J et al.",
    "Pilot study of a Bayesian approach to estimate vancomycin exposure in obese patients with limited",
    "pharmacokinetic sampling. Antimicrob Agents Chemother 2017;61:e02478-16. doi:10.1128/aac.02478-16.",
    "The non-parametric-to-parametric translation and the residual-error assumptions are original to",
    "the Tong 2026 implementation (Discussion: \"Translating the Carreno model into a parametric model",
    "for use with MAP Bayesian estimation also required assumptions about residual variability\").",
    sep = " "
  )
  vignette <- "Tong_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate in this model -- no parameter is scaled by body size, which Tong 2026",
        "highlights in the Discussion as unusual for an obese-patient model. Supplied as a data item;",
        "the InsightRX pipeline supplies it in L/h and the control stream converts inline within the",
        "clearance equation, CL = SCLINTER + SCLSLOPE * (CRCL*16.667), where 1 L/h = 16.667 mL/min.",
        "This model file takes CRCL directly in mL/min, which is the scale the 0.036 slope is",
        "expressed on, so no conversion is applied here; Table S2 states the equation on that scale as",
        "'CL = 0.18 + 0.036 * CRCL'. Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts raw un-normalized mL/min when the source",
        "does not BSA-normalize. Tong 2026 Table 1 reports serum creatinine (median 0.95 mg/dL, range",
        "0.10-10.9) rather than creatinine clearance for this cohort.",
        sep = " "
      ),
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2709L,
    n_studies      = 19L,
    age_range      = "18.2 to 90+ years",
    age_median     = "60.0 years",
    weight_range   = "43.8-318.0 kg",
    weight_median  = "133.0 kg",
    sex_female_pct = 57.1,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Hospitalized adults (>= 18 years) with BMI >= 40 kg/m2 receiving intravenous vancomycin under",
      "routine model-informed precision dosing; at least two doses and at least one measured",
      "concentration required. Patients undergoing haemodialysis at any point during treatment were",
      "excluded, as were patients dosed with a model other than their site's default.",
      sep = " "
    ),
    dose_range     = "Intravenous vancomycin per routine clinical practice; initial doses selected a priori, subsequent doses adapted by MAP Bayesian posterior estimates",
    regions        = "United States (19 hospital systems, patients beginning treatment August 2022 to December 2024)",
    renal_function = "Serum creatinine median 0.95 mg/dL (range 0.10-10.9); haemodialysis patients excluded",
    n_concentrations = 6572L,
    notes          = paste(
      "APPLICATION population from Tong 2026 Table 1 (BMI >= 40 kg/m2 cohort: 2709 patients, 2964",
      "treatment courses, 6572 samples), i.e. the cohort this model was USED to dose as the",
      "pre-intervention default -- NOT the cohort it was estimated from. The DEVELOPMENT population is",
      "Carreno 2017 and is very small: 12 obese patients with 71 concentrations (Tong 2026 Table S2",
      "rows '# patients in development population' = 12 and '# of samples' = 71), fit",
      "non-parametrically. Tong 2026 Discussion notes that this small, non-size-scaled, wide-IIV model",
      "gives MAP Bayesian estimation 'considerable flexibility', so 'small variations in creatinine",
      "clearance or drug concentrations can result in very different PK parameter estimates'. Every",
      "population parameter is a FIXED prior; Tong 2026 estimated none of them.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # ALL population parameters are FIXED priors for MAP Bayesian estimation.
    # Source: Tong 2026 Supplementary data, Code section, "Carreno model" control
    # stream. Every $THETA carries the NONMEM FIX flag. NOTE: the $OMEGA BLOCK(5)
    # in that listing does NOT carry an explicit FIX flag, unlike the sibling
    # modified-Goti and modified-Thomson streams; it is nevertheless a fixed
    # prior, because Tong 2026 estimated no population parameters (Methods,
    # Pharmacokinetic analysis: only individual parameters were fitted, by MAP
    # Bayesian estimation). Encoded with fixed() accordingly -- see vignette
    # Errata.
    # ------------------------------------------------------------------------

    # Clearance is additive in a non-renal intercept and a renal arm:
    #   CL = SCLINTER + SCLSLOPE * CRCL
    # mapped onto the registered multi-component clearance canonicals
    # lcl_nonren (intercept) and lcl_renal (slope), following the additive
    # renal-plus-non-renal precedent in Bulitta_2011_cefpirome.R. Each arm
    # carries its own eta, exactly as in the control stream.
    lcl_nonren <- fixed(log(0.18));  label("Non-renal clearance intercept (L/h)")
    # $THETA(5) 0.18 FIX (TVSCLINTER); Table S2 intercept in "CL = 0.18 + 0.036 * CRCL"
    lcl_renal  <- fixed(log(0.036)); label("Renal clearance slope per mL/min of CRCL (L/h per mL/min)")
    # $THETA(2) 0.036 FIX (TVSCLSLOPE); Table S2 slope in "CL = 0.18 + 0.036 * CRCL"

    lvc <- fixed(log(25.76)); label("Central volume (L; not scaled by body size)")
    # $THETA(1) 25.76 FIX (TVV); Table S2 "V = 25.76"

    # The intercompartmental exchange is parameterised directly as first-order
    # rate constants, each with its own eta:
    #   K12 = TVK12 * EXP(ETA(3));  K21 = TVK21 * EXP(ETA(4))
    # Q and Vp are derived FROM them in $PK (TVQ = K12*V1; V2 = Q/K21).
    lk12 <- fixed(log(2.29)); label("Transfer rate constant central -> peripheral1 (1/h)")
    # $THETA(3) 2.29 FIX (TVK12); Table S2 "Q = K12 * V = 2.29 * 25.76"
    lk21 <- fixed(log(1.44)); label("Transfer rate constant peripheral1 -> central (1/h)")
    # $THETA(4) 1.44 FIX (TVK21); Table S2 "V2 = Q/K21 = 59.0/1.44"

    # Inter-individual variability. $OMEGA BLOCK(5) with all off-diagonal
    # elements 0, so the five etas are declared independently. The stream stores
    # variances on the omega^2 = CV^2 convention (NOT log(CV^2 + 1)): each
    # sqrt(variance) reproduces the Table S2 %CV row exactly. Note how wide these
    # are -- 105.7% and 120.1% CV on the two rate constants.
    #
    # These source traces are on their own lines rather than trailing the eta
    # declarations: rxode2 rewrites a trailing comment on an `eta ~ ...` line into
    # a label() call, and a comment containing double quotes then fails to parse.
    # $OMEGA(1,1); sqrt = 0.4534 -> Table S2 "IIV on V (%CV) 45.3"
    etalvc        ~ fixed(0.20559)
    # $OMEGA(2,2); sqrt = 0.5556 -> Table S2 "55.6% slope"
    etalcl_renal  ~ fixed(0.30864)
    # $OMEGA(3,3); sqrt = 1.0568 -> Table S2 "K12: 105.6"
    etalk12       ~ fixed(1.11676)
    # $OMEGA(4,4); sqrt = 1.2014 -> Table S2 "K21: 120.1"
    etalk21       ~ fixed(1.4433)
    # $OMEGA(5,5); sqrt = 0.2389 -> Table S2 "23.9% intercept"
    etalcl_nonren ~ fixed(0.057068)

    # Combined additive + proportional residual error, from the $ERROR block
    # W = SQRT(IPRED**2 * PROP**2 + ADD**2) with PROP and ADD hardcoded. These
    # two values are ASSUMPTIONS ORIGINAL TO the Tong 2026 implementation, not
    # Carreno 2017 estimates: the original model was fit non-parametrically and
    # Tong 2026 Discussion states the parametric translation "required
    # assumptions about residual variability".
    addSd  <- fixed(0.1); label("Additive residual error (mg/L; assumed by the Tong 2026 implementation)")
    # $ERROR ADD = 0.1; Table S2 "Additive error (mg/L) 0.1"
    propSd <- fixed(0.1); label("Proportional residual error (fraction; assumed by the Tong 2026 implementation)")
    # $ERROR PROP = 0.1; Table S2 "Proportional error (%) 10"
  })
  model({
    # 1. Individual PK parameters. The control stream converts the incoming CRCL
    #    from L/h to mL/min inline (CRCL*16.667); this file takes CRCL already in
    #    mL/min, the scale the 0.036 slope is expressed on.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal  <- exp(lcl_renal  + etalcl_renal) * CRCL
    cl        <- cl_nonren + cl_renal

    vc <- exp(lvc + etalvc)

    # Individual transfer rate constants, exactly as $PK writes them.
    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # Intercompartmental clearance and peripheral volume, reproducing
    # TVQ = K12*V1 and V2 = Q/K21 from $PK.
    #
    # IMPORTANT: the ODEs below must be driven by q / vp, not by kcp / kpc
    # directly. rxSolve() defaults to useLinCmt = TRUE, which rewrites a
    # recognisably linear two-compartment system into a closed form; when the
    # peripheral transfer is written straight from stored micro-constants that
    # rewrite silently drops peripheral1 and solves a ONE-compartment model.
    # Measured on this model at CRCL = 90 mL/min: terminal half-life 5.22 h and
    # Cmax 52.8 mg/L (collapsed) versus the correct 13.82 h and 27.1 mg/L.
    # Routing through q and vp makes the closed-form and ODE solvers agree to
    # machine precision, which the vignette asserts.
    q  <- kcp * vc
    vp <- q / kpc

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system, transcribed from $DES:
    #      DADT(1) = -(CL/V1)*A(1) - K12*A(1) + K21*A(2)
    #      DADT(2) =                 K12*A(1) - K21*A(2)
    #    The stream's third state, DADT(3) = A(1)/V1, is a quadrature-only AUC
    #    accumulator used for MIPD reporting; it has no effect on the PK and is
    #    omitted here (AUC is obtained by NCA in the vignette instead).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observation and error. S1 = V1 in the stream, so dose in mg over volume
    #    in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
