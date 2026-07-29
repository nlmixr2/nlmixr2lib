NA_NA_miridesap <- function() {
  description <- "Population PK/PD model for CPHPC (miridesap, GSK2315698, Ro 63-8695) and serum amyloid P (SAP) in healthy volunteers (CPH113776) and patients with systemic amyloidosis (CPH114527). Two-compartment PK for CPHPC (with first-order subcutaneous absorption from a depot in addition to IV infusion); two-compartment turnover model for SAP with first-order endogenous production and elimination; bimolecular CPHPC + free SAP -> complex binding (treated as effectively irreversible because internalization is fast relative to dissociation) followed by complex internalization. Final-model covariates from Sahota 2015 Eq. 1 and Eq. 2: creatinine clearance (CRCL) modifies CPHPC clearance below an 80 mL/min threshold, hepatic amyloid involvement (DIS_AMYLOID_LIVER) multiplies SAP intercompartmental clearance Q4, whole-body amyloid load (DIS_AMYLOID_LOAD: 0-3) multiplies SAP peripheral volume V4 in two categorical steps, and biological sex (SEXF) multiplies baseline plasma SAP. Distributed in the DDMORE Foundation Model Repository as DDMODEL00000262."
  reference <- "Sahota T, Berges A, Barton S, Cookson L, Zamuner S, Richards D. Target Mediated Drug Disposition Model of CPHPC in Patients With Systemic Amyloidosis. CPT Pharmacometrics Syst Pharmacol. 2015;4(2):e15. doi:10.1002/psp4.15. Companion DDMORE Foundation Model Repository entry: DDMODEL00000262."
  vignette <- "NA_NA_miridesap"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL (CPHPC plasma); ng/mL-equivalent (SAP plasma; equal to mg/L * 1000)"
  )

  ddmore_id    <- "DDMODEL00000262"
  replicate_of <- NULL

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by MDRD (Sahota 2015 Methods). Time-fixed per subject.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sahota 2015 Eq. 1: cl_effect = 1 + e_crcl_cl * (min(CRCL, 80) - 80); the covariate saturates at CRCL >= 80 mL/min (multiplier = 1, no further effect) and reduces CL linearly with declining renal function below 80. Output_real_CPHPC.lst THETA(15) = 0.0153 (paper Table 2 reports 0.015, RSE 5%). Source column 'CRCL' in Simulated_CPHPC_dataset.csv.",
      source_name        = "CRCL"
    ),
    DIS_AMYLOID_LIVER = list(
      description        = "Hepatic amyloid involvement indicator (0 = no liver amyloid, 1 = liver amyloid present).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Sahota 2015 Eq. 2: Q4_amliver = 1 + e_amliver_q4 * DIS_AMYLOID_LIVER; multiplies the typical SAP intercompartmental clearance. Patients with hepatic amyloid have a ~5-fold higher Q4 (1 + 4.01). Source column 'DIS_AMYLOID_LIVER' in Simulated_CPHPC_dataset.csv. Baseline-only; not time-varying.",
      source_name        = "DIS_AMYLOID_LIVER"
    ),
    DIS_AMYLOID_LOAD = list(
      description        = "Whole-body amyloid load categorical score (0 = no amyloid in healthy volunteers; 1 = small; 2 = moderate; 3 = large).",
      units              = "(categorical 0-3)",
      type               = "categorical",
      reference_category = 0,
      notes              = "Sahota 2015 Methods: the score combines organ-by-organ amyloid presence into a single ordinal grade. Eq. 2 collapses categories 0 and 1 into the reference (V4 multiplier = 1); category 2 adds e_amload2_vp_sap; category 3 adds e_amload2_vp_sap + e_amload3_vp_sap. The cumulative additive parameterisation enforces monotonicity. Paper reports moderate = 6.39 (RSE 39%) and large = 26.39 (RSE 26%), matching Output_real_CPHPC.lst TH20 / TH21. Source column 'DIS_AMYLOID_LOAD' in Simulated_CPHPC_dataset.csv.",
      source_name        = "DIS_AMYLOID_LOAD"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male). Source data encode SEX as 1 = male, 2 = female (Sahota 2015 NONMEM convention); derive SEXF = as.integer(SEX == 2) when feeding the bundle's Simulated_CPHPC_dataset.csv into this model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Sahota 2015 Eq. 2: SAP_BASE_sex = 1 + e_sexf_sap0 * SEXF; female SAP baseline is ~30% lower than male. Paper Table 2 reports the effect at -0.30 (RSE 23%), matching Output_real_CPHPC.lst TH19. The paper's source coding is SEX = 1 (male) / 2 (female), so the bundle's source column 'SEX' must be remapped to SEXF before use; this canonical-rename mirrors the SEXF entry in inst/references/covariate-columns.md which accepts the SEX 1/2 form as a source alias.",
      source_name        = "SEX"
    )
  )

  population <- list(
    n_subjects     = 38L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Healthy volunteers (study CPH113776) and patients with systemic amyloidosis (study CPH114527; cohorts 1-4 spanning small-to-large amyloid load and normal-to-moderate-severe renal impairment, per Sahota 2015 Table 1). The two studies were analysed simultaneously in NONMEM via aggregated dataset.",
    dose_range     = "Sahota 2015 Methods: IV bolus, IV infusion (1-48 h durations), and SC dosing (60 mg q8h post-infusion). The bundle's Simulated_CPHPC_dataset.csv exercises IV regimens spanning 5 mg/h x 1 h up to 70 mg/h x 1 h followed by 10 mg/h x 23 h (overall 300 mg/d) and 40 mg/h x 1 h followed by 20 mg/h x 23 h (overall 500 mg/d). The .ctl also exposes a subcutaneous depot (CMT 5) with fixed bioavailability (FSC = 1) and fixed first-order absorption (KASC = 0.4055 on the log scale; ka = exp(0.4055) ~ 1.50 1/h), but the simulated dataset does not exercise SC dosing.",
    regions        = NA_character_,
    notes          = "Subject count (38) is from Output_real_CPHPC.lst 'TOT. NO. OF INDIVIDUALS: 38'. Study identifiers come from the RDF model-has-description and from Sahota 2015 Methods. Finer-grained demographics (age / weight ranges, female %) are not enumerated in the paper Table 1 or in the DDMORE bundle. See the vignette for the Sahota 2015 Table 2 cross-check."
  )

  ini({
    # All parameter values come from Output_real_CPHPC.lst's FINAL PARAMETER
    # ESTIMATE block (lines 384-385) which records the FINAL MODEL fit to the
    # real CPHPC + SAP data from studies CPH113776 + CPH114527. Each value is
    # cross-checked against Sahota 2015 (CPT Pharmacometrics Syst Pharmacol
    # 4:e15) Table 2; the per-parameter inline comments cite the THETA index
    # and the matching Table 2 row. The .ctl shipped in the DDMORE bundle is
    # for a simulated-dataset run with .ctl initial values offset slightly
    # from the real-data final estimates (e.g., .ctl V = 2.76 vs .lst final
    # 2.78); we follow the .lst's FINAL values which are the paper's Table 2.
    #
    # Most THETAs sit on the log scale via the mu-modelling pattern
    #   TVx = EXP(THETA(i));  MU_i = LOG(TVx);  X = EXP(MU_i + ETA(i));
    # so THETA(i) = log(<linear value>) and is assigned directly to `l<param>`.
    # The exception is KINT, where the .ctl writes `TVKINT = THETA(8)` (no EXP
    # wrapper), making THETA(8) the linear value; that one is logged here.

    # CPHPC two-compartment PK
    lcl     <- 1.92            ; label("Log CPHPC clearance at CRCL >= 80 mL/min (log L/h); CL = exp(lcl) = 6.82 L/h (Sahota 2015 Table 2: 6.85)")    # THETA(2)
    lvc     <- 2.78            ; label("Log CPHPC central volume of distribution (log L); Vc = exp(lvc) = 16.12 L (Sahota 2015 Table 2: 16.15)")      # THETA(3)
    lq      <- 0.542           ; label("Log CPHPC inter-compartmental clearance (log L/h); Q = exp(lq) = 1.72 L/h (Sahota 2015 Table 2: 1.72)")        # THETA(4)
    lvp     <- 2.87            ; label("Log CPHPC peripheral volume (log L); Vp = exp(lvp) = 17.64 L (Sahota 2015 Table 2: 17.57)")                    # THETA(5)

    # CRCL effect on CPHPC clearance -- piecewise linear, saturating at CRCL = 80
    e_crcl_cl <- 0.0153        ; label("Slope of CRCL effect on CPHPC CL for CRCL <= 80 mL/min (per mL/min); Sahota 2015 Table 2 reports 0.015")       # THETA(15)

    # Subcutaneous depot -- both fixed; bundle's simulated dataset does not
    # exercise SC dosing but the .ctl supports it via a separate compartment.
    lka     <- fixed(0.4055)   ; label("Log SC absorption rate (log 1/h); ka = exp(lka) = 1.50 1/h - FIXED (Sahota 2015 Table 2: 1.5 FIX)")            # THETA(14) FIXED
    lfdepot <- fixed(0)        ; label("Log SC bioavailability (log unitless); F = exp(lfdepot) = 1 - FIXED")                                          # THETA(13) FIXED

    # SAP turnover and distribution
    lkout    <- -3.07          ; label("Log endogenous SAP elimination rate (log 1/h); kout = exp(lkout) = 0.0464 1/h (Sahota 2015 Table 2: 0.046)")   # THETA(1)
    lsap0    <- 3.44           ; label("Log baseline plasma SAP concentration in males (log mg/L); sap0 = exp(lsap0) = 31.19 mg/L (Sahota 2015: 31.10)") # THETA(6)
    lvp_sap  <- 2.50           ; label("Log SAP peripheral volume at DIS_AMYLOID_LOAD <= 1 (log L); Vp_sap_ref = exp(lvp_sap) = 12.18 L (Sahota 2015 Table 2: 12.15)") # THETA(16)
    lq_sap   <- 1.05           ; label("Log SAP inter-compartmental clearance at DIS_AMYLOID_LIVER = 0 (log L/h); Q_sap_ref = exp(lq_sap) = 2.86 L/h (Sahota 2015: 2.84)") # THETA(17)

    # SAP covariate effects (Sahota 2015 Eq. 2)
    e_amliver_q4      <- 4.01  ; label("DIS_AMYLOID_LIVER multiplicative effect on SAP Q4 (Q4 = Q4_ref * (1 + e_amliver_q4 * DIS_AMYLOID_LIVER)); ~5x increase when liver amyloid present") # THETA(18)
    e_sexf_sap0       <- -0.30 ; label("SEXF (female) multiplicative effect on SAP baseline (SAP_BASE = SAP_BASE_ref * (1 + e_sexf_sap0 * SEXF)); ~30% lower in females") # THETA(19)
    e_amload2_vp_sap  <- 6.39  ; label("DIS_AMYLOID_LOAD step-2 multiplicative effect on SAP V4 (V4 = V4_ref * (1 + e_amload2_vp_sap)) for DIS_AMYLOID_LOAD >= 2; ~7.4x V4 in moderate amyloid load") # THETA(20)
    e_amload3_vp_sap  <- 26.4  ; label("DIS_AMYLOID_LOAD step-3 additional multiplicative effect on V4 (V4 = V4_ref * (1 + e_amload2_vp_sap + e_amload3_vp_sap)) for DIS_AMYLOID_LOAD = 3; ~33.8x V4 in large amyloid load") # THETA(21)

    # CPHPC + free SAP <-> complex binding kinetics
    lkon    <- 14.5            ; label("Log SAP-CPHPC binding on-rate (log L/(mol h)); kon = exp(lkon) = 1.98e6 L/(mol h) (Sahota 2015 Table 2: 1.94e6)") # THETA(7)
    lkint   <- log(5.78)       ; label("Log SAP-CPHPC complex internalization rate (log 1/h); kint = 5.78 1/h (Sahota 2015 Table 2: 5.78)")             # THETA(8) (linear scale in .lst)

    # IIV -- diagonal $OMEGA in the .lst (lines 394-425); no correlations.
    # Variances copied verbatim from Output_real_CPHPC.lst FINAL OMEGA block.
    etalkout    ~ 0.169          # OMEGA1 (KOUT)
    etalcl      ~ 0.0481         # OMEGA2 (CL)
    etalvc      ~ 0.0914         # OMEGA3 (V)
    etalq       ~ fixed(0.0225)  # OMEGA4 (Q) FIXED
    etalvp      ~ fixed(0.0225)  # OMEGA5 (V2) FIXED
    etalsap0    ~ 0.0415         # OMEGA6 (SAP10)
    etalkon     ~ fixed(0.0225)  # OMEGA7 (KON) FIXED
    etalkint    ~ fixed(0.0225)  # OMEGA8 (KINT) FIXED
    etalq_sap   ~ 0.273          # OMEGA9 (Q4)
    etalvp_sap  ~ 0.366          # OMEGA10 (V4)
    etalka      ~ fixed(0.0225)  # OMEGA11 (KSC) FIXED

    # Residual error - both observations are proportional only; the $ERROR
    # block multiplies IPRED by (1 + W * EPS1) for CMT.EQ.1 (CPHPC) and CMT.EQ.3
    # (SAP), with the additive components THETA(10) and THETA(12) FIXED at 0.
    propSd     <- 0.286 ; label("Proportional residual SD on CPHPC plasma observations (fraction); Sahota 2015 Table 2 reports PK = 28.62%")  # THETA(9)
    propSd_sap <- 0.271 ; label("Proportional residual SD on SAP plasma observations (fraction); Sahota 2015 Table 2 reports SAP = 27.10%")   # THETA(11)
  })

  model({
    # ------------------------------------------------------------------
    # Constants - molecular weights and binding-equilibrium constants from
    # the .ctl $PK block (lines 39-40, 101-103).
    # ------------------------------------------------------------------
    # CPHPC molecular weight: 340.37 Da; convert to mg/mol so that mass
    # (mg) divided by MW (mg/mol) gives moles.
    MOLCPH <- 340.37 * 1000          # mg/mol
    # SAP pentamer molecular weight: 5 monomers x 25 kDa each = 125 kDa.
    # The same factor of 1000 converts kDa to mg/mol.
    MOLSAP <- 5 * 25000 * 1000       # mg/mol  (= 1.25e8 mg/mol)
    # KD (= 10 nM) is FIXED in the .ctl and is used only to compute KOFF;
    # KOFF itself is computed in the .ctl but is not used in the $DES (the
    # binding is treated as effectively irreversible because KINT >> KOFF).
    # We omit it from this model for the same reason.

    # ------------------------------------------------------------------
    # Covariate factors
    # ------------------------------------------------------------------
    # CRCL effect on CPHPC clearance - piecewise linear; saturates at
    # CRCL = 80 mL/min (Sahota 2015 Eq. 1). Below 80 the multiplier reduces
    # CL; above 80 it is 1.
    crcl_capped <- (CRCL < 80) * CRCL + (CRCL >= 80) * 80
    crcl_effect <- 1 + e_crcl_cl * (crcl_capped - 80)

    # DIS_AMYLOID_LIVER effect on SAP Q4 (Sahota 2015 Eq. 2; DIS_AMYLOID_LIVER is binary 0/1)
    q4_amliver_factor <- 1 + e_amliver_q4 * DIS_AMYLOID_LIVER

    # DIS_AMYLOID_LOAD effect on SAP V4 (Sahota 2015 Eq. 2; ordinal 0-3 with
    # categories 0 and 1 sharing the reference). Encoded as cumulative
    # indicators to enforce monotonicity.
    amload_ge_2 <- (DIS_AMYLOID_LOAD >= 2) * 1.0
    amload_ge_3 <- (DIS_AMYLOID_LOAD >= 3) * 1.0
    v4_amload_factor <- 1 + e_amload2_vp_sap * amload_ge_2 + e_amload3_vp_sap * amload_ge_3

    # SEXF effect on SAP baseline (Sahota 2015 Eq. 2; SEXF is binary 0=male, 1=female)
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

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    # SAP central uses the CPHPC central volume vc (V3 = V in Sahota 2015).
    k34 <- q_sap / vc
    k43 <- q_sap / vp_sap

    # ------------------------------------------------------------------
    # Derived state variables - free / bound partitioning. The .ctl carries
    # `central` (CPHPC, mg) as TOTAL drug (free + bound) and
    # `total_target` (SAP, mol/L) as TOTAL target (free + bound). The
    # complex compartment carries the bound concentration in mol/L. Free
    # quantities are obtained by subtracting the bound mass / concentration.
    # .ctl $DES lines 152-156.
    # ------------------------------------------------------------------
    central_bound_mass <- complex * vc * MOLCPH        # mg
    central_free       <- central - central_bound_mass # mg
    central_free_conc  <- central_free / (vc * MOLCPH) # mol/L
    sap_free           <- total_target - complex       # mol/L

    # ------------------------------------------------------------------
    # Endogenous SAP production rate (mol/L/h). At steady state the
    # production / volume balances the elimination of free SAP, so
    #   KIN = KOUT * SAP10 * V_sap_central ; KIN/(V*MOLSAP) = KOUT*SAP10/MOLSAP.
    # .ctl $PK line 120.
    # ------------------------------------------------------------------
    sap_production <- kout * sap0 / MOLSAP

    # ------------------------------------------------------------------
    # Initial conditions - SAP starts at the per-subject baseline; bound
    # complex starts at zero (default). .ctl $PK lines 122-126 use
    # A_0(3) = AM3 = SAP10/MOLSAP and A_0(4) = AM4 = SAP10/MOLSAP.
    # ------------------------------------------------------------------
    total_target(0)      <- sap0 / MOLSAP   # mol/L
    target_peripheral(0) <- sap0 / MOLSAP   # mol/L

    # ------------------------------------------------------------------
    # ODE system. Compartments and units:
    #   depot              CPHPC SC depot (mg)
    #   central            CPHPC central total mass (mg, free + bound)
    #   peripheral1        CPHPC peripheral mass (mg)
    #   total_target       SAP central total concentration (mol/L)
    #   target_peripheral  SAP peripheral free concentration (mol/L)
    #   complex            SAP-CPHPC bound complex concentration (mol/L)
    # .ctl $DES lines 138-173.
    # ------------------------------------------------------------------
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot - kel * central_free - k12 * central_free + k21 * peripheral1 - kint * vc * MOLCPH * complex
    d/dt(peripheral1)       <-  k12 * central_free - k21 * peripheral1
    d/dt(total_target)      <-  sap_production - kout * sap_free - kint * complex + k43 * target_peripheral * vp_sap / vc - k34 * sap_free
    d/dt(target_peripheral) <- -k43 * target_peripheral + k34 * sap_free * vc / vp_sap
    d/dt(complex)           <-  kon * sap_free * central_free_conc - kint * complex

    # SC depot bioavailability (= 1 in the .ctl; written out for clarity)
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # Observations - rescaled per the .ctl S1 / S3 scaling factors so
    # that the predictions match the units of the bundle's DV column:
    #   S1 = V/1000  -> Cc = central / vc * 1000 (ng/mL)
    #   S3 = 1/(1000*MOLSAP) -> sap = total_target * MOLSAP * 1000 (ng/mL,
    #     equivalent to mg/L of total SAP * 1000)
    # The proportional residual error in the .ctl applies to IPRED for
    # each output (CMT.EQ.1 -> CPHPC residual; else -> SAP residual).
    # .ctl $ERROR lines 176-198, $PK lines 135-136.
    # ------------------------------------------------------------------
    Cc  <- central      / vc * 1000           # ng/mL (total CPHPC plasma)
    sap <- total_target * MOLSAP * 1000       # ng/mL (total SAP plasma; = mg/L * 1000)

    Cc  ~ prop(propSd)
    sap ~ prop(propSd_sap)
  })
}
