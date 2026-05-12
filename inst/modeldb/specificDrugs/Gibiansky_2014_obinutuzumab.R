Gibiansky_2014_obinutuzumab <- function() {
  description <- "Two-compartment population PK model of obinutuzumab (GA101, glycoengineered type II anti-CD20 mAb) in adults with chronic lymphocytic leukemia (CLL) or non-Hodgkin lymphoma (NHL); clearance is the sum of a time-independent component CL_inf and a mono-exponentially decaying time-dependent component CL_T*exp(-kdes*time), with histology (CLL / BCL / DLBCL / MCL), baseline tumor size, body weight, and sex as covariates (Gibiansky 2014)."
  reference <- "Gibiansky E, Gibiansky L, Carlile DJ, Jamois C, Buchheit V, Frey N. Population Pharmacokinetics of Obinutuzumab (GA101) in Chronic Lymphocytic Leukemia (CLL) and Non-Hodgkin's Lymphoma and Exposure-Response in CLL. CPT Pharmacometrics Syst Pharmacol. 2014;3(10):e144. doi:10.1038/psp.2014.42"
  vignette <- "Gibiansky_2014_obinutuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Estimated power exponent 0.615 (shared on CL_T and CL_inf) and 0.383 (on V1) per Table 3; fixed allometric exponents 0.75 on Q and 1.0 on V2 per Methods. Reference 75 kg per NONMEM control stream (Supplementary Table S1).",
      source_name        = "BW"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source column SEX in the NONMEM control stream encodes 1 = male, 0 = female; SEXF = 1 - SEX inverts the encoding to the canonical female-indicator. Paper's typical-value parameters (CL_T, CL_inf, V1) are reported at the female reference, so effects on CL_T, CL_inf, and V1 are applied via (1 - SEXF) to preserve verbatim source values.",
      source_name        = "SEX"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size, reported as the sum of products of perpendicular diameters of target lesions (SPPD)",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Source uses BSIZ in mm^2 (SPPD convention typical for B-cell lymphoma and CLL nodal disease, where diameters are paired into products). Used in the model only via the binary stratum indicator (BSIZ <= 1750 mm^2 vs > 1750 mm^2) introduced as CODT in the NONMEM control stream; the threshold was identified by Gibiansky 2014 as providing a better fit than splitting at other values or treating BSIZ as continuous (Methods 'Covariate model development'). The canonical TUMSZ register pools SPPD with sum-of-diameters and sum-of-linear-diameters constructs; for this paper the SPPD convention is in mm^2.",
      source_name        = "BSIZ"
    ),
    TUMTP_BCL = list(
      description        = "B-cell lymphoma (BCL) histology indicator, 1 = BCL, 0 = other histology",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other histology; CLL is the typical-value reference category)",
      notes              = "Source NONMEM control stream encodes diagnosis as DIS with integer levels 1 = CLL, 2 = BCL, 3 = DLBCL, 4 = MCL (Supplementary Table S1). Decompose into TUMTP_BCL = as.integer(DIS == 2). In Gibiansky 2014 the BCL category is a pooled residual B-cell-lymphoma group that includes follicular lymphoma (all GAUDI patients were FL; Table 1). The paper applies a shared effect (theta13, ratio 0.834) to BCL and DLBCL on both CL_T and CL_inf; the BCL effect is therefore implemented as the composite (TUMTP_BCL + TUMTP_DLBCL) inside model().",
      source_name        = "DIS"
    ),
    TUMTP_DLBCL = list(
      description        = "Diffuse large B-cell lymphoma (DLBCL) histology indicator, 1 = DLBCL, 0 = other histology",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other histology; CLL is the typical-value reference category)",
      notes              = "Source NONMEM control stream encodes diagnosis as DIS with integer levels 1 = CLL, 2 = BCL, 3 = DLBCL, 4 = MCL (Supplementary Table S1). Decompose into TUMTP_DLBCL = as.integer(DIS == 3). The paper applies a shared effect (theta13, ratio 0.834) to BCL and DLBCL on both CL_T and CL_inf; see TUMTP_BCL notes.",
      source_name        = "DIS"
    ),
    TUMTP_MCL = list(
      description        = "Mantle cell lymphoma (MCL) histology indicator, 1 = MCL, 0 = other histology",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other histology; CLL is the typical-value reference category)",
      notes              = "Source NONMEM control stream encodes diagnosis as DIS with integer levels 1 = CLL, 2 = BCL, 3 = DLBCL, 4 = MCL (Supplementary Table S1). Decompose into TUMTP_MCL = as.integer(DIS == 4). The MCL effect (theta14, ratio 1.75) is applied separately to CL_T and CL_inf. Together with TUMTP_BCL and TUMTP_DLBCL, the indicators define the kdes-NHL composite ((TUMTP_BCL + TUMTP_DLBCL + TUMTP_MCL); theta12, ratio 2.08).",
      source_name        = "DIS"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 678L,
    n_studies      = 4L,
    n_observations = 12634L,
    age_range      = "22-89 years",
    age_median     = "65.7 years (mean)",
    weight_range   = "40-140 kg",
    weight_median  = "75.6 kg (mean)",
    sex_female_pct = 42.9,
    disease_state  = "CD20+ B-cell malignancies (Gibiansky 2014 Table 2): chronic lymphocytic leukemia 342 / 678 (50.4%); B-cell lymphoma (BCL; predominantly follicular lymphoma in GAUDI) 286 / 678 (42.2%); diffuse large B-cell lymphoma 30 / 678 (4.4%); mantle cell lymphoma 20 / 678 (2.9%).",
    dose_range     = "Intravenous obinutuzumab 200-2000 mg per infusion (maximum infusion rate 400 mg/h). GAUGUIN NHL: 400 or 1600/800 mg per cycle on day 1 (with day 8 of cycle 1). GAUGUIN CLL: 1000 mg cycle 1 days 1/8/15 and day 1 of subsequent cycles. GAUDI: 400 mg with FC or CHOP, or 1600/800 mg standalone; FL first-line 1000 mg q21-28d. GAUSS phase I dose escalation 200-2000 mg; phase II 1000 mg weekly induction. CLL11: 1000 mg cycle 1 days 1/8/15 and day 1 of cycles 2-6 (28-day cycles).",
    regimen        = "Intermittent IV infusion across induction and maintenance regimens; CLL11 used the 1000 mg cycle-1 loading (days 1, 8, 15) + 1000 mg q28d maintenance regimen carried forward as the obinutuzumab labelled regimen in CLL.",
    regions        = "Multi-regional phase I-III programme.",
    studies        = "GAUGUIN (BO20999, NCT00517530, phase I/II; 131 patients, 3446 samples); GAUDI (BO21000, NCT00825149, phase Ib; 134 patients, 3634 samples); GAUSS (BO21003, NCT00576758, phase I/II; 105 patients, 2327 samples); CLL11 (BO21004, NCT01010061, phase III; 308 patients, 3227 samples).",
    bsiz_summary   = "Baseline tumor size (sum of products of perpendicular diameters, SPPD) mean (SD) 5390 (19,100) mm^2 across the full cohort; per-study means range 4420-6030 mm^2 (Gibiansky 2014 Table 2). The threshold 1750 mm^2 splits the cohort into 'low' and 'high' strata for the kdes effect.",
    bcell_summary  = "Baseline B-cell count mean (SD) 37.4 (66.6) x 10^9/L across cohort; CLL11 mean 77.75 x 10^9/L vs 1.62-11.8 x 10^9/L in the NHL studies. Baseline B-cell count and diagnosis were confounded (high counts in CLL) and diagnosis was retained as the primary covariate in lieu of B-cell count.",
    excluded_obs   = "74 postdose observations (0.6%) below the assay LLOQ of 4.05 ng/mL were excluded from the analysis.",
    notes          = "Pooled phase I-III population PK dataset. The four trials enrolled previously treated CLL/NHL patients (GAUGUIN, GAUDI, GAUSS) plus previously untreated comorbid CLL patients (CLL11). Age and renal function (creatinine clearance) were tested but did not significantly affect PK parameters; race and antidrug antibody status were not retained as covariates."
  )

  ini({
    # ---- Structural PK parameters (Gibiansky 2014 Table 3, final covariate model) ----
    # Typical values are reported as exp(theta_i) for log-scale parameters; we store the
    # log of the back-transformed value. Reference subject: CLL diagnosis, female (SEX = 0),
    # body weight 75 kg, BSIZ > 1750 mm^2 (Supplementary Table S1 NONMEM control stream
    # mean-zero centring on CLL/female/75 kg/high-BSIZ).
    lkdes    <- log(0.0359);    label("Decay rate of time-dependent clearance kdes (1/day)")               # Gibiansky 2014 Table 3 exp(theta1) = 0.0359
    lcl_time <- log(0.231);     label("Initial value of time-dependent clearance CL_T (L/day)")            # Gibiansky 2014 Table 3 exp(theta2) = 0.231
    lcl_ss   <- log(0.0828);    label("Time-independent (steady-state) clearance CL_inf (L/day)")          # Gibiansky 2014 Table 3 exp(theta3) = 0.0828
    lvc      <- log(2.76);      label("Central volume of distribution V1 (L)")                              # Gibiansky 2014 Table 3 exp(theta4) = 2.76
    lvp      <- log(1.01);      label("Peripheral volume of distribution V2 (L)")                           # Gibiansky 2014 Table 3 exp(theta5) = 1.01
    lq       <- log(1.29);      label("Intercompartmental clearance Q (L/day)")                             # Gibiansky 2014 Table 3 exp(theta6) = 1.29

    # ---- Allometric / weight-scaling exponents ----
    # Estimated WT exponents on the central parameters (theta7 shared between CL_T and CL_inf,
    # theta8 on V1); fixed allometric exponents on the peripheral parameters per Methods
    # ('allometric scaling with fixed power coefficients on clearances and volumes of 0.75
    # and 1, respectively').
    e_wt_cl  <-  0.615;         label("Power exponent of (WT/75) on CL_T and CL_inf (shared; unitless)")    # Gibiansky 2014 Table 3 theta7 = 0.615
    e_wt_vc  <-  0.383;         label("Power exponent of (WT/75) on V1 (unitless)")                          # Gibiansky 2014 Table 3 theta8 = 0.383
    e_wt_q   <- fixed(0.75);    label("Power exponent of (WT/75) on Q (unitless; fixed allometric)")        # Gibiansky 2014 Methods (allometric exponent 0.75 on clearances, fixed)
    e_wt_vp  <- fixed(1.0);     label("Power exponent of (WT/75) on V2 (unitless; fixed allometric)")       # Gibiansky 2014 Methods (allometric exponent 1.0 on volumes, fixed)

    # ---- Sex effects (paper-reported values for SEX = 1 male; applied via (1 - SEXF)) ----
    e_sex_cl_time <- log(1.49); label("log(ratio) of CL_T for males vs females (applied via (1-SEXF))")     # Gibiansky 2014 Table 3 exp(theta9) = 1.49
    e_sex_cl_ss   <- log(1.22); label("log(ratio) of CL_inf for males vs females (applied via (1-SEXF))")   # Gibiansky 2014 Table 3 exp(theta10) = 1.22
    e_sex_vc      <- log(1.18); label("log(ratio) of V1 for males vs females (applied via (1-SEXF))")       # Gibiansky 2014 Table 3 exp(theta11) = 1.18

    # ---- Diagnosis effects ----
    # NHL composite (BCL + DLBCL + MCL) effect on kdes; shared BCL/DLBCL effect and a
    # separate MCL effect on CL_T and CL_inf (effects shared between the two CL components
    # per the NONMEM MU expressions; theta13 in MU_2 and MU_3 and theta14 in MU_2 and MU_3).
    e_nhl_kdes    <- log(2.08);  label("log(ratio) of kdes for NHL vs CLL (any of BCL/DLBCL/MCL)")           # Gibiansky 2014 Table 3 exp(theta12) = 2.08
    e_bcldlbcl_cl <- log(0.834); label("log(ratio) of CL_T and CL_inf for BCL or DLBCL vs CLL (shared)")     # Gibiansky 2014 Table 3 exp(theta13) = 0.834
    e_mcl_cl      <- log(1.75);  label("log(ratio) of CL_T and CL_inf for MCL vs CLL (shared)")              # Gibiansky 2014 Table 3 exp(theta14) = 1.75

    # ---- Baseline tumor size effect ----
    e_bsizlow_kdes <- log(2.65); label("log(ratio) of kdes for BSIZ <= 1750 mm^2 vs BSIZ > 1750 mm^2")       # Gibiansky 2014 Table 3 exp(theta15) = 2.65

    # ---- Inter-individual variability (omega^2 = log-scale variance; Table 3) ----
    # CV% = 100 * sqrt(exp(omega^2) - 1): kdes 201%, CL_T 122%, CL_inf 41.5%, V1 18.6%,
    # V2 65.9%, Q 120%. Correlations of random effects were tested and not retained
    # in the final model per Methods ('correlation of the random effects was not included').
    etalkdes    ~ 1.62;   # Gibiansky 2014 Table 3 Omega(1,1) (CV = 201%)
    etalcl_time ~ 0.907;  # Gibiansky 2014 Table 3 Omega(2,2) (CV = 122%)
    etalcl_ss   ~ 0.159;  # Gibiansky 2014 Table 3 Omega(3,3) (CV = 41.5%)
    etalvc      ~ 0.034;  # Gibiansky 2014 Table 3 Omega(4,4) (CV = 18.6%)
    etalvp      ~ 0.361;  # Gibiansky 2014 Table 3 Omega(5,5) (CV = 65.9%)
    etalq       ~ 0.89;   # Gibiansky 2014 Table 3 Omega(6,6) (CV = 120%)

    # ---- Residual error (Gibiansky 2014 Table 3) ----
    # Combined additive + proportional model on observed plasma concentration (ug/mL).
    # The published NONMEM error block also includes an extra IIV term ETA(7) on the
    # proportional residual SD (Omega(7,7) = 0.274, CV = 56.1%); this is a heteroscedastic
    # residual-error layer that is not directly supported by nlmixr2's standard
    # `add() + prop()` form and is therefore omitted. The model retains the typical-value
    # proportional and additive SDs from the paper. See the vignette's Assumptions and
    # deviations section for the rationale.
    propSd <- 0.1783;     label("Proportional residual error (fraction; sqrt(sigma^2_prop = 0.0318))")      # Gibiansky 2014 Table 3 Sigma(1,1) = 0.0318
    addSd  <- 0.1646;     label("Additive residual error (ug/mL; sqrt(sigma^2_add = 0.0271))")              # Gibiansky 2014 Table 3 Sigma(2,2) = 0.0271 (ug/mL)^2
  })

  model({
    # ---- Derived covariate indicators ---------------------------------------------
    # Source NONMEM control stream uses SEX = 1 for males; canonical SEXF = 1 for
    # females. The paper's typical-value parameters anchor the female reference, so
    # the male-vs-female effect is applied via (1 - SEXF).
    sex_male  <- 1 - SEXF

    # Diagnosis composite indicators per the NONMEM control stream:
    #   CODD  = 1 if DIS in {2, 3} (BCL or DLBCL) -> bcl_dlbcl
    #   CODD1 = 1 if DIS == 4 (MCL)               -> TUMTP_MCL (used directly)
    #   kdes NHL indicator = (CODD + CODD1) = any NHL histology -> nhl
    bcl_dlbcl <- TUMTP_BCL + TUMTP_DLBCL
    nhl       <- TUMTP_BCL + TUMTP_DLBCL + TUMTP_MCL

    # Baseline-tumor-size low-stratum indicator (NONMEM: CODT = 1 if BSIZ <= 1750).
    # TUMSZ holds the per-subject baseline SPPD (mm^2); the model uses only the
    # categorical threshold per Gibiansky 2014 Methods.
    bsiz_low  <- (TUMSZ <= 1750)

    # ---- Individual structural parameters ----------------------------------------
    # Allometric scaling: CL_T, CL_inf, V1 use estimated WT exponents (theta7, theta8);
    # Q and V2 retain fixed allometric coefficients (0.75 and 1.0). Sex and diagnosis
    # effects are applied as log-ratio offsets to the central-compartment parameters.
    kdes    <- exp(lkdes    + e_nhl_kdes     * nhl       + e_bsizlow_kdes * bsiz_low + etalkdes)
    cl_time <- exp(lcl_time + e_sex_cl_time  * sex_male  + e_wt_cl  * log(WT/75) +
                   e_bcldlbcl_cl * bcl_dlbcl + e_mcl_cl  * TUMTP_MCL + etalcl_time)
    cl_ss   <- exp(lcl_ss   + e_sex_cl_ss    * sex_male  + e_wt_cl  * log(WT/75) +
                   e_bcldlbcl_cl * bcl_dlbcl + e_mcl_cl  * TUMTP_MCL + etalcl_ss)
    vc      <- exp(lvc      + e_sex_vc       * sex_male  + e_wt_vc  * log(WT/75) + etalvc)
    vp      <- exp(lvp      + e_wt_vp        * log(WT/75) + etalvp)
    q       <- exp(lq       + e_wt_q         * log(WT/75) + etalq)

    # ---- Time-varying clearance --------------------------------------------------
    # Gibiansky 2014 Methods 'Base PK model development': CL(t) = CL_T * exp(-kdes * t)
    # + CL_inf. The integration time variable 'time' refers to the time since the start
    # of the simulation / first observation, matching NONMEM TIME.
    cl <- cl_time * exp(-kdes * time) + cl_ss

    # ---- Micro-constants ---------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system --------------------------------------------------------------
    # Two-compartment linear disposition with a time-varying central-compartment
    # elimination rate kel(t) = cl(t) / vc. Equivalent to Gibiansky 2014 $DES.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Observation and residual error -----------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
