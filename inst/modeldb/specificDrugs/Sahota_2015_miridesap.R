Sahota_2015_miridesap <- function() {
  description <- "Target-mediated drug disposition (TMDD) PK/PD model for CPHPC (miridesap, GSK2315698, Ro 63-8695) and serum amyloid P (SAP) in healthy volunteers (study CPH113776) and patients with systemic amyloidosis (study CPH114527). Two-compartment PK for CPHPC (IV plus first-order subcutaneous depot); two-compartment turnover model for SAP with first-order endogenous production and elimination; bimolecular CPHPC + free SAP -> complex binding treated as effectively irreversible (KOFF set to zero because the complex internalisation rate is much faster than the dissociation rate). Final-model covariates (Sahota 2015 Eq. 1 and Eq. 2): creatinine clearance modifies CPHPC clearance below an 80 mL/min threshold; hepatic amyloid involvement multiplies SAP intercompartmental clearance Q4; whole-body amyloid load (categorical 0-3) multiplies SAP peripheral volume V4 in two cumulative steps; biological sex multiplies baseline plasma SAP."
  reference <- "Sahota T, Berges A, Barton S, Cookson L, Zamuner S, Richards D. Target Mediated Drug Disposition Model of CPHPC in Patients With Systemic Amyloidosis. CPT Pharmacometrics Syst Pharmacol. 2015;4(2):e15. doi:10.1002/psp4.15."
  vignette <- "Sahota_2015_miridesap"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL (CPHPC plasma); ng/mL-equivalent for SAP (mg/L * 1000)"
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated by the MDRD formula at baseline (Sahota 2015 Methods, Dataset production).",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sahota 2015 Eq. 1: cl_effect = 1 + e_crcl_cl * (min(CRCL, 80) - 80); the covariate saturates at CRCL >= 80 mL/min (multiplier = 1) and linearly reduces CL with declining renal function below the 80 mL/min threshold. Reported clearance of 6.85 L/h applies at CRCL >= 80 (Sahota 2015 Table 2 footnote a). Baseline-only; not time-varying.",
      source_name        = "CRCL"
    ),
    DIS_AMYLOID_LIVER = list(
      description        = "Hepatic amyloid involvement indicator (0 = no liver amyloid; 1 = liver amyloid present).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Sahota 2015 Eq. 2: Q4_amliver = 1 + e_amliver_q4 * DIS_AMYLOID_LIVER multiplies the typical SAP intercompartmental clearance. Patients with hepatic amyloid have a roughly 5-fold higher Q4 (1 + 4.01). Baseline-only; not time-varying. Renamed from canonical AMLIVER to DIS_AMYLOID_LIVER on 2026-06-19 per the canonical-register standardization audit.",
      source_name        = "AMLIVER"
    ),
    DIS_AMYLOID_LOAD = list(
      description        = "Whole-body amyloid load ordinal score (0 = no amyloid in healthy volunteers; 1 = small; 2 = moderate; 3 = large).",
      units              = "(categorical 0-3)",
      type               = "categorical",
      reference_category = 0,
      notes              = "Sahota 2015 Methods (Dataset production) defines the score; Eq. 2 collapses categories 0 and 1 into the reference (V4 multiplier = 1), then adds e_amload2_vp_sap once at DIS_AMYLOID_LOAD>=2 and adds e_amload3_vp_sap again at DIS_AMYLOID_LOAD=3 to enforce monotonicity. Final model yields V4 multipliers of 1, 1, 7.39, and 33.78 for DIS_AMYLOID_LOAD 0, 1, 2, and 3 respectively (Sahota 2015 Results: V4 increased 7.4-fold and 33.78-fold in moderate and large amyloid loads relative to small/none). Baseline-only; not time-varying. Renamed from canonical AMLOAD to DIS_AMYLOID_LOAD on 2026-06-19 per the canonical-register standardization audit.",
      source_name        = "AMLOAD"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female; 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Sahota 2015 Eq. 2: SAP_BASE_sex = 1 + e_sexf_sap0 * SEXF, with e_sexf_sap0 = -0.30 giving ~30% lower baseline SAP in females (Sahota 2015 Results, consistent with Nelson et al. 1991 reference range of 21 mg/L in women vs 32 mg/L in men). Paper encoded SEX as 1 = male and 2 = female; convert to canonical SEXF = as.integer(SEX == 2) when ingesting paper-formatted data.",
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy volunteers (CPH113776; 21 male subjects) and patients with systemic amyloidosis (CPH114527; cohorts 1-4 spanning small-to-large whole-body amyloid load and normal-to-moderate-severe renal impairment per Sahota 2015 Table 1).",
    dose_range     = "Sahota 2015 Methods: CPH113776 used 1 h IV infusions (5-70 mg) and 24 h IV infusion regimens (induction + maintenance, total 86-960 mg). CPH114527 used 48 h IV infusions (124.8-1440 mg) followed by one or three SC doses (10-60 mg). The model also supports SC dosing via the depot compartment with KSC and F fixed at the prior point estimates.",
    regions        = NA_character_,
    notes          = "Subject count and study identifiers from Sahota 2015 Methods (Adaptive PK-PD studies). Finer demographics (precise age and weight ranges, female %) are not enumerated in the paper Table 1; Supplementary Materials S1 and S2 were not on disk for this extraction."
  )

  ini({
    # All parameter values come from Sahota 2015 Table 2 (final model parameter
    # estimates) and Eq. 1 / Eq. 2. Each parameter line carries an in-file
    # citation to the matching Table 2 row, Eq. 1 / Eq. 2 term, or paper
    # narrative. IIV variances are reported in the paper as BSV percent (CV%);
    # they are converted to log-normal omega^2 via omega^2 = log(1 + CV^2).
    # Parameters that the paper holds constant carry a 'fixed()' wrapper.

    # ---- CPHPC two-compartment PK ----------------------------------------
    lcl  <- log(6.85)  ; label("Log CPHPC clearance at CRCL >= 80 mL/min (log L/h); CL = 6.85 L/h")   # Sahota 2015 Table 2 (clearance, footnote a)
    lvc  <- log(16.15) ; label("Log CPHPC central volume of distribution (log L); V1 = 16.15 L")      # Sahota 2015 Table 2 (central volume V1)
    lq   <- log(1.72)  ; label("Log CPHPC inter-compartmental clearance (log L/h); Q2 = 1.72 L/h")    # Sahota 2015 Table 2 (intercompartmental clearance Q2)
    lvp  <- log(17.57) ; label("Log CPHPC peripheral volume (log L); V2 = 17.57 L")                   # Sahota 2015 Table 2 (peripheral volume V2)

    # ---- CRCL effect on CPHPC CL (Sahota 2015 Eq. 1) ---------------------
    e_crcl_cl <- 0.015 ; label("Linear slope of CRCL effect on CPHPC CL below the 80 mL/min threshold (per mL/min)")  # Sahota 2015 Table 2 (CRCL x clearance) and Eq. 1

    # ---- Subcutaneous depot (KSC and F fixed at prior estimates) ---------
    lka     <- fixed(log(1.5)) ; label("Log SC absorption rate constant (log 1/h); KSC = 1.5 1/h - FIXED")                  # Sahota 2015 Table 2 (KSC: 1.5 FIX); Methods
    lfdepot <- fixed(log(1))   ; label("Log SC bioavailability (log unitless); F = 1 - FIXED")                              # Sahota 2015 Methods: bioavailability assumed complete for SC doses

    # ---- SAP turnover ---------------------------------------------------
    lkout <- log(0.046) ; label("Log SAP elimination rate constant (log 1/h); KOUT = 0.046 1/h")    # Sahota 2015 Table 2 (KOUT)
    lsap0 <- log(31.10) ; label("Log baseline plasma SAP concentration in males (log mg/L); SAP_BASE = 31.10 mg/L")  # Sahota 2015 Table 2 (SAP_BASE)

    # ---- SAP distribution -----------------------------------------------
    # V3 = V1 by structural constraint (Sahota 2015 Methods: 'complex formation
    # assumed to occur in well stirred conditions in a common central
    # compartment (V1 = V3)'); not a separately estimated parameter.
    lvp_sap <- log(12.15) ; label("Log SAP peripheral volume at DIS_AMYLOID_LOAD <= 1 (log L); V4 reference = 12.15 L")  # Sahota 2015 Table 2 (peripheral volume V4)
    lq_sap  <- log(2.84)  ; label("Log SAP inter-compartmental clearance at DIS_AMYLOID_LIVER = 0 (log L/h); Q4 reference = 2.84 L/h")  # Sahota 2015 Table 2 (intercompartmental clearance Q4)

    # ---- SAP covariate effects (Sahota 2015 Eq. 2) -----------------------
    e_amliver_q4     <- 4.01  ; label("DIS_AMYLOID_LIVER multiplicative effect on SAP Q4 (Q4 = Q4_ref * (1 + e_amliver_q4 * DIS_AMYLOID_LIVER)); ~5-fold increase when liver amyloid present")  # Sahota 2015 Table 2 (Amyloid liver x Q4); Eq. 2
    e_sexf_sap0      <- -0.30 ; label("SEXF (female) multiplicative effect on SAP baseline (SAP_BASE = SAP_BASE_ref * (1 + e_sexf_sap0 * SEXF)); ~30% lower in females")    # Sahota 2015 Table 2 (Gender x SAP baseline); Eq. 2
    e_amload2_vp_sap <- 6.39  ; label("DIS_AMYLOID_LOAD step at DIS_AMYLOID_LOAD >= 2: V4 multiplier increment (V4 = V4_ref * (1 + e_amload2_vp_sap * I(DIS_AMYLOID_LOAD>=2) + e_amload3_vp_sap * I(DIS_AMYLOID_LOAD>=3)))")  # Sahota 2015 Table 2 (Amyloid load x V4, moderate); Eq. 2
    e_amload3_vp_sap <- 26.39 ; label("DIS_AMYLOID_LOAD additional step at DIS_AMYLOID_LOAD = 3: V4 multiplier increment (large amyloid load adds on top of the DIS_AMYLOID_LOAD>=2 step)")  # Sahota 2015 Table 2 (Amyloid load x V4, large); Eq. 2

    # ---- CPHPC + free SAP <-> complex binding ----------------------------
    lkon  <- log(1.94e6) ; label("Log SAP-CPHPC binding on-rate (log L/(mol h)); KON = 1.94e6 L/(mol h)")  # Sahota 2015 Table 2 (association rate KON)
    lkint <- log(5.78)   ; label("Log SAP-CPHPC complex internalisation rate (log 1/h); KINT = 5.78 1/h")  # Sahota 2015 Table 2 (complex elimination KINT)

    # ---- IIV ------------------------------------------------------------
    # Paper reports BSV percent (CV%); converted to log-normal omega^2 via
    # omega^2 = log(1 + CV^2). Parameters reported as '15% FIX' use omega^2
    # = log(1 + 0.15^2) = 0.02225.
    etalcl     ~ 0.0470          # BSV 21.93% (Sahota 2015 Table 2, clearance row)
    etalvc     ~ 0.0874          # BSV 30.23% (Sahota 2015 Table 2, central volume V1)
    etalq      ~ fixed(0.0225)   # BSV 15% FIX (Sahota 2015 Table 2, intercompartmental clearance Q2)
    etalvp     ~ fixed(0.0225)   # BSV 15% FIX (Sahota 2015 Table 2, peripheral volume V2)
    etalka     ~ fixed(0.0225)   # BSV 15% FIX (Sahota 2015 Table 2, KSC)
    etalkout   ~ 0.1560          # BSV 41.09% (Sahota 2015 Table 2, KOUT)
    etalsap0   ~ 0.0406          # BSV 20.36% (Sahota 2015 Table 2, SAP_BASE)
    etalq_sap  ~ 0.2413          # BSV 52.24% (Sahota 2015 Table 2, intercompartmental clearance Q4)
    etalvp_sap ~ 0.3117          # BSV 60.48% (Sahota 2015 Table 2, peripheral volume V4)
    etalkon    ~ fixed(0.0225)   # BSV 15% FIX (Sahota 2015 Table 2, KON)
    etalkint   ~ fixed(0.0225)   # BSV 15% FIX (Sahota 2015 Table 2, KINT)

    # ---- Residual error -------------------------------------------------
    # Both observations carry proportional residual error only; the paper
    # reports residual SD as a percent for each output stream.
    propSd     <- 0.2862 ; label("Proportional residual SD on CPHPC plasma observations (fraction)") # Sahota 2015 Table 2 (Residual error, PK: 28.62%)
    propSd_sap <- 0.2710 ; label("Proportional residual SD on SAP plasma observations (fraction)")   # Sahota 2015 Table 2 (Residual error, SAP: 27.10%)
  })

  model({
    # ------------------------------------------------------------------
    # Molecular-weight constants for the mass <-> mole conversion at the
    # binding interface. CPHPC molecular mass is 340 Da (Sahota 2015
    # Discussion); SAP pentamer molecular mass is 127,310 Da (Sahota 2015
    # Discussion).
    # ------------------------------------------------------------------
    MW_CPHPC <- 340 * 1000        # mg/mol  (= 3.40e5)
    MW_SAP   <- 127310 * 1000     # mg/mol  (= 1.2731e8)

    # ------------------------------------------------------------------
    # Covariate factors (Sahota 2015 Eq. 1 and Eq. 2)
    # ------------------------------------------------------------------
    # CRCL effect on CPHPC clearance: piecewise linear, saturating at 80
    # (Sahota 2015 Eq. 1). CRCL is truncated to the 80 mL/min threshold
    # when CRCL > 80, so the multiplier is exactly 1 in subjects with
    # normal renal function.
    crcl_capped <- (CRCL < 80) * CRCL + (CRCL >= 80) * 80
    crcl_effect <- 1 + e_crcl_cl * (crcl_capped - 80)

    # DIS_AMYLOID_LIVER effect on SAP Q4 (Eq. 2; binary indicator)
    q4_amliver_factor <- 1 + e_amliver_q4 * DIS_AMYLOID_LIVER

    # DIS_AMYLOID_LOAD effect on SAP V4 (Eq. 2; cumulative indicators enforce
    # monotonicity. Categories 0 and 1 share the reference (multiplier = 1);
    # category 2 adds e_amload2_vp_sap; category 3 adds it again plus the
    # additional step e_amload3_vp_sap.)
    amload_ge_2 <- (DIS_AMYLOID_LOAD >= 2) * 1.0
    amload_ge_3 <- (DIS_AMYLOID_LOAD >= 3) * 1.0
    v4_amload_factor <- 1 + e_amload2_vp_sap * amload_ge_2 + e_amload3_vp_sap * amload_ge_3

    # SEXF effect on SAP baseline (Eq. 2)
    sap0_sexf_factor <- 1 + e_sexf_sap0 * SEXF

    # ------------------------------------------------------------------
    # Individual parameters
    # ------------------------------------------------------------------
    cl     <- exp(lcl     + etalcl)     * crcl_effect
    vc     <- exp(lvc     + etalvc)
    q      <- exp(lq      + etalq)
    vp     <- exp(lvp     + etalvp)
    ka     <- exp(lka     + etalka)
    kout   <- exp(lkout   + etalkout)
    sap0   <- exp(lsap0   + etalsap0)   * sap0_sexf_factor   # mg/L
    kon    <- exp(lkon    + etalkon)                          # L/(mol*h)
    kint   <- exp(lkint   + etalkint)                         # 1/h
    q_sap  <- exp(lq_sap  + etalq_sap)  * q4_amliver_factor
    vp_sap <- exp(lvp_sap + etalvp_sap) * v4_amload_factor

    # Micro-constants. Note: SAP central volume V3 equals CPHPC central
    # volume V1 (Sahota 2015 Methods: complex formed in V1 = V3), so the
    # SAP central distribution uses vc rather than a separately
    # parameterised SAP central volume.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k34 <- q_sap / vc
    k43 <- q_sap / vp_sap

    # ------------------------------------------------------------------
    # Free / bound partitioning. The 'central' state holds the total
    # CPHPC mass in mg (free plus bound); the 'total_target' state holds
    # total SAP central concentration in mol/L; the 'complex' state holds
    # the bound SAP-CPHPC complex concentration in mol/L. Free quantities
    # are obtained by subtraction.
    # ------------------------------------------------------------------
    central_bound_mass <- complex * vc * MW_CPHPC      # mg
    central_free       <- central - central_bound_mass # mg
    central_free_conc  <- central_free / (vc * MW_CPHPC) # mol/L
    sap_free           <- total_target - complex       # mol/L

    # ------------------------------------------------------------------
    # Endogenous SAP production rate (mol/L/h). At baseline the central
    # SAP compartment is at steady state, so KIN per central volume
    # balances KOUT times the baseline SAP concentration. Sahota 2015
    # Results: 'At baseline, free SAP was assumed to be in equilibrium
    # between production and elimination in the central compartment.'
    # ------------------------------------------------------------------
    sap_production <- kout * sap0 / MW_SAP

    # ------------------------------------------------------------------
    # Initial conditions: total SAP starts at the per-subject baseline
    # (sap0 / MW_SAP, in mol/L); SAP peripheral starts at the same
    # baseline so the two-compartment system begins at equilibrium; the
    # complex compartment defaults to zero before dosing.
    # ------------------------------------------------------------------
    total_target(0)      <- sap0 / MW_SAP
    target_peripheral(0) <- sap0 / MW_SAP

    # ------------------------------------------------------------------
    # ODE system. Compartments and units:
    #   depot              CPHPC SC depot (mg)
    #   central            CPHPC central total mass (mg, free + bound)
    #   peripheral1        CPHPC peripheral mass (mg)
    #   total_target       SAP central total concentration (mol/L)
    #   target_peripheral  SAP peripheral free concentration (mol/L)
    #   complex            SAP-CPHPC complex concentration (mol/L)
    # Mass-balance and binding kinetics follow the TMDD schematic
    # (Sahota 2015 Figure 3) with KOFF = 0 because the complex
    # internalisation rate is much faster than dissociation.
    # ------------------------------------------------------------------
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot - kel * central_free - k12 * central_free + k21 * peripheral1 - kint * vc * MW_CPHPC * complex
    d/dt(peripheral1)       <-  k12 * central_free - k21 * peripheral1
    d/dt(total_target)      <-  sap_production - kout * sap_free - kint * complex + k43 * target_peripheral * vp_sap / vc - k34 * sap_free
    d/dt(target_peripheral) <- -k43 * target_peripheral + k34 * sap_free * vc / vp_sap
    d/dt(complex)           <-  kon * sap_free * central_free_conc - kint * complex

    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # Observations. Outputs are rescaled to the assay-reported units:
    # CPHPC plasma concentration in ng/mL (central mass / vc gives mg/L,
    # multiplied by 1000 to convert to ng/mL); SAP plasma concentration
    # in ng/mL (total_target mol/L times pentamer MW gives mg/L, times
    # 1000 to convert to ng/mL).
    # ------------------------------------------------------------------
    Cc  <- central      / vc * 1000           # ng/mL (total CPHPC plasma)
    sap <- total_target * MW_SAP * 1000       # ng/mL (total SAP plasma)

    Cc  ~ prop(propSd)
    sap ~ prop(propSd_sap)
  })
}
