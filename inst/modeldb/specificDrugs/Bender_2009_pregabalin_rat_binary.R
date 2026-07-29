Bender_2009_pregabalin_rat_binary <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model for",
    "pregabalin in male Sprague-Dawley rats following a 2 h",
    "intravenous infusion (4 or 10 mg/kg/h) in a chronic-",
    "constriction-injury (CCI) neuropathic-pain model, with the",
    "concomitant administration of sildenafil encoded as a",
    "BINARY presence indicator (CONMED_SILDENAFIL). Sildenafil",
    "presence reduces pregabalin clearance by a fixed fraction",
    "(theta_SLD = 0.302, i.e. 30.2% reduction) per the paper's",
    "discrete-covariate parameterisation; the alternative",
    "continuous saturable-metabolite parameterisation is",
    "encoded in the companion file",
    "Bender_2009_pregabalin_rat_smetab.R. Crossover design with",
    "two occasions per rat (Day 1 / Day 4 with a washout)",
    "carries between-occasion variability on CL and Vc",
    "multiplexed by the OCC indicator. Parameter values from",
    "Bender 2009 Table IV (Binary Sildenafil Covariate column)."
  )
  reference <- paste(
    "Bender G, Gosset J, Florian J, Tan K, Field M, Marshall S,",
    "DeJongh J, Bies R, Danhof M. (2009). Population",
    "pharmacokinetic model of the pregabalin-sildenafil",
    "interaction in rats: Application of simulation to",
    "preclinical PK-PD study design.",
    "Pharmaceutical Research 26(10):2259-2269.",
    "doi:10.1007/s11095-009-9942-y.",
    sep = " "
  )
  vignette <- "Bender_2009_pregabalin_rat"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CONMED_SILDENAFIL = list(
      description        = paste(
        "Binary per-occasion indicator for concomitant sildenafil",
        "coadministration (1 = sildenafil bolus + 6 h steady-state",
        "infusion during this experimental occasion, 0 = saline)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Per-occasion (NOT per-subject) in the Bender 2009 crossover:",
        "each rat received saline on one occasion and sildenafil on",
        "the other. Paper alias: SLDB.",
        "Sildenafil regimen on the active occasion: 2 mg/kg bolus",
        "(1 mL/kg) over 10 min followed by 12 mg/kg/h infusion for",
        "6 h, producing steady-state saturation of PDE-5 and the",
        "N-methyl sildenafil active metabolite. Values: 0 or 1."
      ),
      source_name        = "SLDB"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued occasion / period indicator for between-",
        "occasion variability multiplexing. Values 1 and 2 identify",
        "the two crossover days (Day 1 and Day 4) within each rat."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into binary indicators oc1 and",
        "oc2 that multiplex the BOV etas on log-CL and log-Vc.",
        "The two occasions are separated by a >= 3 day washout per",
        "Methods (Pharmacokinetic Study Design)."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "rat (male Sprague-Dawley)",
    n_subjects     = 28L,
    n_studies      = 1L,
    age_range      = "(not reported in the source publication; surgically prepared adult rats)",
    weight_range   = "200-250 g at receipt (Charles River Laboratories, Margate, UK)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Chronic-constriction-injury (CCI) model of neuropathic pain",
      "(Bennett & Xie). Four loose ligatures around the right",
      "sciatic nerve produce a peripheral mononeuropathy with",
      "static and dynamic allodynia resembling the human",
      "neuropathic-pain phenotype."
    ),
    dose_range     = paste(
      "Pregabalin IV constant-rate infusion through the jugular",
      "venous catheter for 2 hours at either 4 mg/kg/h or 10",
      "mg/kg/h. Each rat received both dose levels on the two",
      "experimental days. Sildenafil arm: 2 mg/kg bolus (1 mL/kg)",
      "over 10 min followed by 12 mg/kg/h infusion for 6 h.",
      "Saline arm: matching saline bolus + infusion."
    ),
    regions        = "United Kingdom (Pfizer Global Research and Development, Pfizer UK)",
    notes          = paste(
      "28 male Sprague-Dawley rats, divided evenly (n = 7) into four",
      "crossover treatment groups (1A, 1B, 2A, 2B per Table I).",
      "Jugular venous catheters (drug administration) and carotid",
      "artery catheters (blood sampling, Dilab automated sampler)",
      "were implanted. Blood samples were taken at sparse",
      "D-optimal-selected times between 0.33 and 22 h on Day 1 and",
      "Day 4 with a >= 3-day washout (Table II). All animals were",
      "housed under 12:12 h light-dark cycle with food and water ad",
      "libitum. Plasma was analysed by LC-MS/MS (LLOQ pregabalin =",
      "70 ng/mL, sildenafil and N-methyl sildenafil metabolite =",
      "1 ng/mL each; assay linear over 100-10,000 ng/mL pregabalin,",
      "3-2,000 ng/mL sildenafil, 5-2,000 ng/mL metabolite).",
      "See Bender 2009 Methods (Animals, Surgical Procedures, Drug",
      "Administration and Blood Sampling) and Tables I-II.",
      "Body weight, age, time post CCI-surgery, and time post",
      "catheterization surgery were tested as continuous PK",
      "covariates but were not retained -- the cohort was",
      "intentionally homogeneous (inbred strain, controlled",
      "conditions), so only sildenafil exposure had a significant",
      "OFV impact. Parameters are reported in absolute units (L",
      "and L/h) at the cohort body-weight scale; users simulating",
      "for a different rat weight should scale CL and Vc",
      "allometrically downstream."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Structural PK parameters -- Bender 2009 Table IV, "Binary
    # Sildenafil Covariate" column (page 2264). NONMEM ADVAN3
    # TRANS4 two-compartment with first-order elimination from
    # the central compartment. Parameter mapping:
    #   NONMEM CL -> lcl, V1 -> lvc, Q -> lq, V2 -> lvp.
    # RSE% from Table IV trail the value comments.
    # ------------------------------------------------------------
    lcl  <- log(0.052)  ; label("Clearance (L/h)")                          # Table IV: CL = 0.052 (RSE 14.59%)
    lvc  <- log(0.271)  ; label("Central volume of distribution V1 (L)")    # Table IV: V1 = 0.271 (RSE 3.38%)
    lq   <- log(0.0151) ; label("Inter-compartmental clearance Q (L/h)")    # Table IV: Q = 0.0151 (RSE 35.76%)
    lvp  <- log(2.77)   ; label("Peripheral volume of distribution V2 (L)") # Table IV: V2 = 2.77 (RSE 51.99%)

    # ------------------------------------------------------------
    # Sildenafil covariate effect on CL (binary form). Bender
    # 2009 Methods: "When SLDB was 0 (sildenafil absence = 0,
    # sildenafil presence = 1), TVCL equals theta_CL and when
    # SLDB was 1, the theta_SLD term was subtracted from the
    # population estimate of CL", i.e.
    #     CL = theta_CL * (1 - SLDB * theta_SLD)
    # Encoded as the unitless proportional reduction fraction.
    # ------------------------------------------------------------
    e_sild_cl <- 0.302 ; label("Proportional reduction in CL when sildenafil is coadministered (fraction)")  # Table IV: theta_SLD = 0.302 (RSE 23.05%)

    # ------------------------------------------------------------
    # Between-subject variability (BSV, "omega^2_1" in Bender 2009
    # Table IV). Table IV reports the NONMEM-internal variance
    # values directly; encoded as the eta variance.
    # ------------------------------------------------------------
    etalcl ~ 0.017  # Table IV: omega^2_1 CL = 0.017 (RSE 15.91%)
    etalvc ~ 0.018  # Table IV: omega^2_1 V1 = 0.018 (RSE 9.26%)
    etalq  ~ 0.439  # Table IV: omega^2_1 Q  = 0.439 (RSE 40.87%); large -- Q is poorly identified in the binary covariate model

    # ------------------------------------------------------------
    # Between-occasion variability (BOV, "omega^2_2" in Table
    # IV) on CL and Vc, multiplexed by OCC. The crossover design
    # uses 2 occasions (Day 1 / Day 4) per rat. Following the
    # NONMEM $OMEGA BLOCK(1) SAME translation pattern of
    # Aregbe_2012_alvespimycin.R: occasion-1 carries the
    # estimated variance and occasion-2 fixes the same value.
    # ------------------------------------------------------------
    etaiov_cl_1 ~ 0.061        # Table IV: omega^2_2 CL = 0.061 (RSE 19.52%)
    etaiov_cl_2 ~ fix(0.061)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)

    etaiov_vc_1 ~ 0.015        # Table IV: omega^2_2 V1 = 0.015 (RSE 8.75%)
    etaiov_vc_2 ~ fix(0.015)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)

    # ------------------------------------------------------------
    # Residual error. Bender 2009 Results: "A proportional error
    # model was selected to account for the residual variability."
    # Table IV reports the variance sigma^2 = 0.029; the
    # proportional SD on the linear concentration scale is
    # sqrt(0.029) ~ 0.1703 (~17.0% CV).
    # ------------------------------------------------------------
    propSd <- 0.1703 ; label("Proportional residual error (fraction)")  # Table IV: sigma^2 = 0.029 (RSE 18.53%); SD = sqrt(0.029)
  })

  model({
    # ------------------------------------------------------------
    # 1. Decompose OCC into per-occasion binary indicators and
    #    build the multiplexed BOV etas. Aregbe 2012 alvespimycin
    #    precedent.
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # ------------------------------------------------------------
    # 2. Individual PK parameters. Log-normal IIV on CL/Vc/Q and
    #    log-normal BOV on CL/Vc. The binary sildenafil covariate
    #    multiplicatively reduces typical CL by the fraction
    #    e_sild_cl when CONMED_SILDENAFIL = 1.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_cl) * (1 - e_sild_cl * CONMED_SILDENAFIL)
    vc <- exp(lvc + etalvc + iov_vc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # ------------------------------------------------------------
    # 3. Two-compartment IV PK with first-order elimination from
    #    the central compartment (NONMEM ADVAN3 TRANS4
    #    mass-balance ODEs). Dose lands in central; the 2 h
    #    infusion duration is encoded on the dose record.
    # ------------------------------------------------------------
    d/dt(central)     <-  q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central     -  q       / vp * peripheral1

    # ------------------------------------------------------------
    # 4. Observation: pregabalin plasma concentration. Dose in
    #    mg / Vc in L gives concentration in mg/L; multiplying
    #    by 1000 yields ng/mL to match the units declared in
    #    metadata and the LC-MS/MS assay range (pregabalin
    #    100-10,000 ng/mL). Sanity check: at 10 mg/kg/h x 0.225
    #    kg rat the predicted Cmax at 2 h is ~22 mg/L = 22,000
    #    ng/mL, matching Bender 2009 Fig. 2 (Results).
    # ------------------------------------------------------------
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
