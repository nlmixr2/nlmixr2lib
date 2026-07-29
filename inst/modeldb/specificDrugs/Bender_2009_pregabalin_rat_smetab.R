Bender_2009_pregabalin_rat_smetab <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model for",
    "pregabalin in male Sprague-Dawley rats following a 2 h",
    "intravenous infusion (4 or 10 mg/kg/h) in a chronic-",
    "constriction-injury (CCI) neuropathic-pain model, with the",
    "concomitant administration of sildenafil encoded as a",
    "CONTINUOUS saturable inhibition driven by the time-varying",
    "plasma concentration of sildenafil's active N-methyl",
    "metabolite (SLDM). Effective CL = theta_CL * (1 - SLDM /",
    "(theta_SLD + SLDM)) with theta_SLD = 1350 ng/mL acting as",
    "the IC50 of metabolite-driven inhibition. Statistically",
    "the preferred parameterisation in the paper (delta-OFV =",
    "-42.6 vs the no-covariate base; the simpler binary form is",
    "in the companion file Bender_2009_pregabalin_rat_binary.R",
    "with delta-OFV = -8.5). Crossover design with two",
    "occasions per rat (Day 1 / Day 4 with a washout) carries",
    "between-occasion variability on CL and Vc multiplexed by",
    "the OCC indicator. Parameter values from Bender 2009 Table",
    "IV (Continuous Sildenafil Metabolite Covariate column)."
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
    SLDM = list(
      description        = paste(
        "Time-varying plasma concentration of the active N-methyl",
        "metabolite of sildenafil (ng/mL), used as the saturable-",
        "inhibition driver on pregabalin CL."
      ),
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying covariate column to be supplied by the user.",
        "In the Bender 2009 study, SLDM was measured directly in",
        "plasma alongside the pregabalin concentration (LC-MS/MS;",
        "LLOQ 1 ng/mL; linear over 5-2,000 ng/mL). The N-methyl",
        "metabolite reached Cmax ~ 2,100 ng/mL at 4-7 h post bolus",
        "and was roughly steady-state during the 6 h sildenafil",
        "infusion. Set SLDM = 0 to recover the no-sildenafil",
        "typical-CL prediction. Paper alias: [SLDM] (Results,",
        "Eq. unnumbered above Discussion)."
      ),
      source_name        = "[SLDM]"
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
      "Day 4 with a >= 3-day washout (Table II). Plasma was analysed",
      "by LC-MS/MS for pregabalin, sildenafil, and the N-methyl",
      "sildenafil metabolite simultaneously; the metabolite's",
      "calibration range was 5-2,000 ng/mL with LLOQ 1 ng/mL.",
      "See Bender 2009 Methods and Tables I-II.",
      "Body weight, age, time post CCI-surgery, and time post",
      "catheterization surgery were tested as continuous PK",
      "covariates but were not retained -- the cohort was",
      "intentionally homogeneous, so only sildenafil exposure had",
      "a significant OFV impact. Parameters are reported in",
      "absolute units (L and L/h) at the cohort body-weight scale;",
      "users simulating for a different rat weight should scale",
      "CL and Vc allometrically downstream."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Structural PK parameters -- Bender 2009 Table IV,
    # "Continuous Sildenafil Metabolite Covariate" column (page
    # 2264). NONMEM ADVAN3 TRANS4 two-compartment with first-
    # order elimination from the central compartment. Parameter
    # mapping:
    #   NONMEM CL -> lcl, V1 -> lvc, Q -> lq, V2 -> lvp.
    # RSE% from Table IV trail the value comments.
    # ------------------------------------------------------------
    lcl  <- log(0.051)  ; label("Clearance (L/h)")                          # Table IV: CL = 0.051 (RSE 8.83%)
    lvc  <- log(0.272)  ; label("Central volume of distribution V1 (L)")    # Table IV: V1 = 0.272 (RSE 3.27%)
    lq   <- log(0.0229) ; label("Inter-compartmental clearance Q (L/h)")    # Table IV: Q = 0.0229 (RSE 18.08%)
    lvp  <- log(2.45)   ; label("Peripheral volume of distribution V2 (L)") # Table IV: V2 = 2.45 (RSE 38.82%)

    # ------------------------------------------------------------
    # Sildenafil-metabolite covariate effect on CL (continuous
    # saturable form). Bender 2009 Results (Eq. above
    # Discussion):
    #     CL = theta_CL * (1 - [SLDM] / (theta_SLD + [SLDM]))
    # theta_SLD is the IC50 of [SLDM] on CL (concentration at
    # which CL is reduced by 50%), in the same units as [SLDM]
    # (ng/mL).
    # ------------------------------------------------------------
    e_sldm_cl <- 1350 ; label("IC50 of [SLDM]-driven saturable inhibition on CL (ng/mL)")  # Table IV: theta_SLD = 1350 (RSE 28.37%)

    # ------------------------------------------------------------
    # Between-subject variability (BSV, "omega^2_1" in Bender
    # 2009 Table IV) on CL, Vc, and Q. Q has much smaller BSV in
    # this model (0.081) than in the binary model (0.439) since
    # the saturable-metabolite parameterisation absorbs much of
    # the apparent between-subject variability in Q. Table IV
    # reports the NONMEM-internal variance values directly;
    # encoded as the eta variance.
    # ------------------------------------------------------------
    etalcl ~ 0.013  # Table IV: omega^2_1 CL = 0.013 (RSE 14.63%)
    etalvc ~ 0.016  # Table IV: omega^2_1 V1 = 0.016 (RSE 8.72%)
    etalq  ~ 0.081  # Table IV: omega^2_1 Q  = 0.081 (RSE 15.94%)

    # ------------------------------------------------------------
    # Between-occasion variability (BOV, "omega^2_2" in Table
    # IV) on CL and Vc, multiplexed by OCC. The crossover design
    # uses 2 occasions (Day 1 / Day 4) per rat. NONMEM $OMEGA
    # BLOCK(1) SAME translation pattern (Aregbe 2012
    # alvespimycin precedent): occasion-1 carries the estimated
    # variance and occasion-2 fixes the same value.
    # ------------------------------------------------------------
    etaiov_cl_1 ~ 0.071        # Table IV: omega^2_2 CL = 0.071 (RSE 18.52%)
    etaiov_cl_2 ~ fix(0.071)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)

    etaiov_vc_1 ~ 0.016        # Table IV: omega^2_2 V1 = 0.016 (RSE 8.33%)
    etaiov_vc_2 ~ fix(0.016)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)

    # ------------------------------------------------------------
    # Residual error. Bender 2009 Results: "A proportional error
    # model was selected to account for the residual variability."
    # Table IV reports the variance sigma^2 = 0.029; the
    # proportional SD on the linear concentration scale is
    # sqrt(0.029) ~ 0.1703 (~17.0% CV).
    # ------------------------------------------------------------
    propSd <- 0.1703 ; label("Proportional residual error (fraction)")  # Table IV: sigma^2 = 0.029 (RSE 18.85%); SD = sqrt(0.029)
  })

  model({
    # ------------------------------------------------------------
    # 1. Decompose OCC into per-occasion binary indicators and
    #    build the multiplexed BOV etas.
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2

    # ------------------------------------------------------------
    # 2. Individual PK parameters. Log-normal IIV on CL/Vc/Q and
    #    log-normal BOV on CL/Vc. The continuous saturable
    #    sildenafil-metabolite covariate multiplicatively
    #    reduces typical CL by the saturating function
    #    1 - SLDM/(e_sldm_cl + SLDM); when SLDM = 0 the factor
    #    is 1 (no inhibition), when SLDM = e_sldm_cl the factor
    #    is 0.5 (50% inhibition), and when SLDM >> e_sldm_cl
    #    the factor approaches 0 (full inhibition).
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_cl) *
          (1 - SLDM / (e_sldm_cl + SLDM))
    vc <- exp(lvc + etalvc + iov_vc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # ------------------------------------------------------------
    # 3. Two-compartment IV PK with first-order elimination from
    #    the central compartment (NONMEM ADVAN3 TRANS4
    #    mass-balance ODEs).
    # ------------------------------------------------------------
    d/dt(central)     <-  q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central     -  q       / vp * peripheral1

    # ------------------------------------------------------------
    # 4. Observation: pregabalin plasma concentration. Dose in
    #    mg / Vc in L gives concentration in mg/L; multiplying
    #    by 1000 yields ng/mL to match the units declared in
    #    metadata and the LC-MS/MS assay range.
    # ------------------------------------------------------------
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
