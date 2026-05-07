NA_NA_miridesap <- function() {
  description <- "Population PK/PD model for CPHPC (miridesap, GSK2315698, Ro 63-8695) and serum amyloid P (SAP) in healthy volunteers and patients with systemic amyloidosis. Two-compartment PK for CPHPC (with first-order subcutaneous absorption from a depot in addition to IV infusion); two-compartment turnover model for SAP with first-order endogenous production and elimination; bimolecular CPHPC + free SAP -> complex binding (treated as effectively irreversible because internalization is fast relative to dissociation) followed by complex internalization. Renal function (CRCL) modifies CPHPC clearance below a 80 mL/min threshold. Distributed in the DDMORE Foundation Model Repository as DDMODEL00000262; the bundle ships only the BASE-MODEL .ctl, so this implementation reproduces the base structural model with the CRCL covariate. The Output_real_CPHPC.lst's FINAL MODEL adds organ-specific amyloid-load covariates (AMLIVER, AMSPLEEN, AMBONE, ...) whose .mod file is NOT in the bundle and is not extracted here."
  reference <- "DDMORE Foundation Model Repository: DDMODEL00000262. No linked publication identified; the RDF metadata names two contributing GSK studies (CPH113776 in healthy volunteers and CPH114527 in amyloidosis patients) but no corresponding journal publication is on disk in this worktree and no PubMed match was returned. See vignettes/articles/NA_NA_miridesap.Rmd for the full source-trace and Errata."
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
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Used in the piecewise-linear CL effect: cl_effect = 1 + e_crcl_cl * (min(CRCL, 80) - 80). The covariate saturates at CRCL >= 80 mL/min (multiplier = 1, no further effect on CL); below 80 mL/min, CL is reduced linearly with declining renal function. .ctl $PK lines 42-49 (Executable_simulated_CPHPC_dataset.ctl).",
      source_name        = "CRCL"
    )
  )

  population <- list(
    n_subjects     = 38L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Healthy volunteers (study CPH113776) and patients with systemic amyloidosis (study CPH114527). The DDMODEL00000262.rdf model-has-description names both contributing studies. The split between studies and finer-grained demographics are not in the DDMORE bundle; no linked publication was identified.",
    dose_range     = "Bundle's Simulated_CPHPC_dataset.csv contains IV regimens spanning 5 mg/h x 1 h (single 5 mg loading) up to 70 mg/h x 1 h followed by 10 mg/h x 23 h (overall 300 mg/d) and 40 mg/h x 1 h followed by 20 mg/h x 23 h (overall 500 mg/d) given as continuous infusions to the CPHPC central compartment. The .ctl also exposes a subcutaneous depot (CMT 5) with fixed bioavailability (FSC = 1) and fixed first-order absorption (KASC = 0.4055 on the log scale; ka = exp(0.4055) ~ 1.50 1/h), but the simulated dataset does not exercise SC dosing.",
    regions        = NA_character_,
    notes          = "Subject count (38) is from Output_real_CPHPC.lst 'TOT. NO. OF INDIVIDUALS: 38' (the .lst is for the FINAL MODEL with patient organ-load covariates, not the base model extracted here, but the underlying study population is the same). Study identifiers come from the RDF model-has-description. See vignettes/articles/NA_NA_miridesap.Rmd Errata for the full set of bundle-vs-implementation caveats, including the choice to extract the base model rather than the FINAL MODEL with the patient-organ amyloid-load covariates that the bundle does not ship a re-runnable .mod for."
  )

  ini({
    # All structural-parameter values come from Executable_simulated_CPHPC_dataset.ctl
    # ($THETA block, lines 200-216). Output_simulated_CPHPC_dataset.lst confirms
    # them as final estimates: the .ctl runs $EST METHOD=IMP EONLY=1 (importance-
    # sampling evaluation only, no optimization), so the FINAL PARAMETER ESTIMATE
    # block reproduces the .ctl initial values to within rounding. The original
    # paper-derived final estimates therefore come straight from the .ctl. No
    # linked publication was identified to cross-check against.
    #
    # The .ctl is mu-referenced: most THETAs sit on the log scale via the pattern
    #   TVx = EXP(THETA(i));  MU_i = LOG(TVx);  X = EXP(MU_i + ETA(i)).
    # Under that pattern THETA(i) = log(<linear value>), so the value is assigned
    # directly to the corresponding `l<param>` in `ini()`. The exception is KINT,
    # where the .ctl writes `TVKINT = THETA(8)` (no EXP wrapper), making THETA(8)
    # the linear value; that one is logged explicitly.

    # CPHPC two-compartment PK (.ctl $PK lines 50-76)
    lcl     <- 1.92            ; label("Log CPHPC clearance at CRCL >= 80 mL/min (log L/h); CL = exp(lcl) = 6.82 L/h")     # THETA(2)
    lvc     <- 2.76            ; label("Log CPHPC central volume of distribution (log L); Vc = exp(lvc) = 15.8 L")          # THETA(3)
    lq      <- 0.595           ; label("Log CPHPC inter-compartmental clearance (log L/h); Q = exp(lq) = 1.81 L/h")          # THETA(4)
    lvp     <- 2.84            ; label("Log CPHPC peripheral volume (log L); Vp = exp(lvp) = 17.1 L")                        # THETA(5)

    # CRCL effect on CPHPC clearance -- piecewise linear, saturating at CRCL = 80
    e_crcl_cl <- 0.0152        ; label("Slope of CRCL effect on CPHPC CL for CRCL <= 80 mL/min (1 / (mL/min))")              # THETA(15)

    # Subcutaneous depot -- both fixed; the bundle's simulated dataset does not
    # exercise SC dosing but the .ctl supports it via a separate compartment.
    lka     <- fixed(0.4055)   ; label("Log SC absorption rate (log 1/h); ka = exp(lka) = 1.50 1/h - FIXED")                 # THETA(14) FIXED
    lfdepot <- fixed(0)        ; label("Log SC bioavailability (log unitless); F = exp(lfdepot) = 1 - FIXED")                # THETA(13) FIXED

    # SAP turnover and distribution (.ctl $PK lines 78-118)
    lkout   <- -3.07           ; label("Log endogenous SAP elimination rate (log 1/h); kout = exp(lkout) = 0.0464 1/h")      # THETA(1)
    lsap0   <- 3.37            ; label("Log baseline plasma SAP concentration (log mg/L); sap0 = exp(lsap0) = 29.1 mg/L")    # THETA(6)
    lvp_sap <- 3.47            ; label("Log SAP peripheral volume (log L); Vp_sap = exp(lvp_sap) = 32.1 L")                  # THETA(16)
    lq_sap  <- 1.49            ; label("Log SAP inter-compartmental clearance (log L/h); Q_sap = exp(lq_sap) = 4.44 L/h")    # THETA(17)

    # CPHPC + free SAP <-> complex binding kinetics
    lkon    <- 14.6            ; label("Log SAP-CPHPC binding on-rate (log L/(mol h)); kon = exp(lkon) = 2.19e6 L/(mol h)")  # THETA(7)
    lkint   <- log(5.71)       ; label("Log SAP-CPHPC complex internalization rate (log 1/h); kint = 5.71 1/h")              # THETA(8) (linear scale in .ctl)

    # IIV -- diagonal $OMEGA in the .ctl (lines 217-227); no correlations.
    # Variances copied verbatim from the .ctl OMEGA block, which the
    # Output_simulated lst confirms as final estimates.
    etalkout    ~ 0.125          # OMEGA1 (KOUT)
    etalcl      ~ 0.0474         # OMEGA2 (CL)
    etalvc      ~ 0.0836         # OMEGA3 (V)
    etalq       ~ fixed(0.0225)  # OMEGA4 (Q) FIXED
    etalvp      ~ fixed(0.0225)  # OMEGA5 (V2) FIXED
    etalsap0    ~ 0.0625         # OMEGA6 (SAP10)
    etalkon     ~ fixed(0.0225)  # OMEGA7 (KON) FIXED
    etalkint    ~ fixed(0.0225)  # OMEGA8 (KINT) FIXED
    etalq_sap   ~ 1.11           # OMEGA9 (Q4)
    etalvp_sap  ~ 2.05           # OMEGA10 (V4)
    etalka      ~ fixed(0.0225)  # OMEGA11 (KSC) FIXED

    # Residual error - both observations are proportional only; the $ERROR block
    # multiplies IPRED by (1 + W * EPS1) for both CMT.EQ.1 (CPHPC) and CMT.EQ.3
    # (SAP), with the additive components THETA(10) and THETA(12) FIXED at 0.
    propSd     <- 0.286 ; label("Proportional residual SD on CPHPC plasma observations (fraction)")  # THETA(9)
    propSd_sap <- 0.268 ; label("Proportional residual SD on SAP plasma observations (fraction)")    # THETA(11)
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
    # CRCL effect on CPHPC clearance - piecewise linear; saturates at
    # CRCL = 80 mL/min. Below 80 the multiplier reduces CL; above 80 it is 1.
    # .ctl $PK lines 42-49.
    # ------------------------------------------------------------------
    crcl_capped <- (CRCL < 80) * CRCL + (CRCL >= 80) * 80
    crcl_effect <- 1 + e_crcl_cl * (crcl_capped - 80)

    # ------------------------------------------------------------------
    # Individual parameters
    # ------------------------------------------------------------------
    cl     <- exp(lcl     + etalcl)     * crcl_effect
    vc     <- exp(lvc     + etalvc)
    q      <- exp(lq      + etalq)
    vp     <- exp(lvp     + etalvp)
    ka     <- exp(lka     + etalka)
    kout   <- exp(lkout   + etalkout)
    sap0   <- exp(lsap0   + etalsap0)            # mg/L
    kon    <- exp(lkon    + etalkon)             # L/(mol*h)
    kint   <- exp(lkint   + etalkint)            # 1/h
    q_sap  <- exp(lq_sap  + etalq_sap)
    vp_sap <- exp(lvp_sap + etalvp_sap)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    # SAP central uses the CPHPC central volume vc (V3 = V in the .ctl).
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
    total_target(0)     <- sap0 / MOLSAP   # mol/L
    target_peripheral(0) <- sap0 / MOLSAP  # mol/L

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
