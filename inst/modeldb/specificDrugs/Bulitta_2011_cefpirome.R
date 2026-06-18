Bulitta_2011_cefpirome <- function() {
  description <- paste(
    "Three-compartment population PK model for IV cefpirome with",
    "simultaneous fit of plasma concentrations and amounts excreted",
    "unchanged in urine. Built from a pooled cohort of 24 Caucasian",
    "adults: 12 cystic fibrosis (CF) patients and 12 healthy volunteers",
    "(HVs) each given a single 10-min IV infusion of 2 g cefpirome.",
    "Body size is captured by allometric scaling on lean body mass (LBM)",
    "with fixed exponents 0.75 on clearance terms and 1.0 on volumes",
    "(reference LBM = 53 kg). Total clearance is split into an estimated",
    "renal arm (CL_R, urinary recovery is tracked in the canonical urine",
    "compartment) and a non-renal arm (CL_NR). A CF / HV cohort indicator",
    "DIS_CF (1 = CF patient, 0 = HV reference) carries three disease-",
    "specific scale factors estimated by the paper: FCYF_CLR = 1.07",
    "applied to CL_R, FCYF_CLNR = 1.13 applied to CL_NR, and FCYF_VSS",
    "= 0.98 applied uniformly to V1 (central), V2 (shallow peripheral),",
    "and V3 (deep peripheral). The inter-compartmental clearances Q12",
    "(central <-> shallow) and Q23 (central <-> deep) are shared across",
    "cohorts. Typical-value clearance and volume estimates are anchored",
    "to DIS_CF = 0 (HV reference) per the DIS_CF covariate convention",
    "registered in inst/references/covariate-columns.md."
  )
  reference <- paste(
    "Bulitta JB, Kinzig M, Landersdorfer CB, Holzgrabe U, Stephan U,",
    "Sorgel F. Comparable Population Pharmacokinetics and Pharmacodynamic",
    "Breakpoints of Cefpirome in Cystic Fibrosis Patients and Healthy",
    "Volunteers. Antimicrob Agents Chemother. 2011 Jun;55(6):2927-2936.",
    "doi:10.1128/AAC.01484-10"
  )
  vignette <- "Bulitta_2011_cefpirome"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL_R, CL_NR, Q12, Q23 (fixed exponent 0.75)",
        "and V1, V2, V3 (fixed exponent 1.0) with reference LBM = 53 kg per",
        "Bulitta 2011 Population PK analysis, body size paragraph and",
        "Table 3 footnote a. LBM was computed by the Cheymol / James formula",
        "(Bulitta 2011 Table 1 footnote a citing Cheymol 1972 / James 1976);",
        "the cohort medians were 45.7 kg (CF) and 50.0 kg (HV)."
      ),
      source_name        = "LBM"
    ),
    DIS_CF = list(
      description        = paste(
        "Cystic-fibrosis cohort indicator: 1 = adult CF patient, 0 = healthy",
        "adult volunteer (reference). Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "The paper does not assign a column name to the CF/HV cohort flag;",
        "the CF-vs-HV split is encoded in the structural NONMEM/S-ADAPT",
        "model via three group-scale factors FCYF_CLR, FCYF_CLNR, and",
        "FCYF_VSS that multiply the HV typical values for CF subjects",
        "(Bulitta 2011 Methods, Population PK analysis, Between-subject-",
        "variability model paragraph). Typical lparam values in this file",
        "are anchored to DIS_CF = 0 (HV) so the FCYF factors are recovered",
        "as positive log-scale effects via exp(e_cf_<param> * DIS_CF)."
      ),
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "18-35 years (CF 18-34, median 22.5; HV 20-35, median 29)",
    weight_range   = "31.5-85.0 kg total body weight (CF 31.5-66.5, median 53.3; HV 53.0-85.0, median 63.6)",
    lbm_range      = "26.2-62.7 kg lean body mass (CF 26.2-55.9, median 45.7; HV 41.8-62.7, median 50.0)",
    height_range   = "140-183 cm (CF 140-183, median 170; HV 161-182, median 175)",
    bmi_range      = "13.2-28.4 kg/m^2 (CF 13.2-20.3, median 19.0; HV 17.7-28.4, median 20.6)",
    crcl_range     = "89.3-164 mL/min Cockcroft-Gault for a nominal 70 kg subject (CF 89.3-164, median 131; HV 92.3-150, median 116)",
    sex_female_pct = 41.7,
    race_ethnicity = "100% Caucasian (Methods, Subjects paragraph).",
    disease_state  = paste(
      "Two parallel groups of similar lean body mass: 12 adult CF patients",
      "and 12 healthy volunteers, all with normal renal function and no",
      "evidence of acute pulmonary exacerbation at the time of dosing."
    ),
    dose_range     = "Single 10-min IV infusion of 2 g (= 2000 mg) cefpirome dissolved in 20 mL water for injection.",
    regions        = "Germany (single-centre study by IBMP, Nurnberg-Heroldsberg).",
    notes          = paste(
      "Baseline demographics from Bulitta 2011 Table 1. Final estimation by",
      "NONMEM VI level 1.2 with FOCE+INTERACTION and a parallel S-ADAPT 1.57",
      "beta MC-PEM analysis; the paper used the NONMEM estimates for the",
      "Monte Carlo simulations and the values in Table 3 (NONMEM column) are",
      "the ones encoded here. The cefpirome assay (HPLC) was linear",
      "between 0.614 and 322 mg/L in plasma and between 17.1 and 8000 mg/L",
      "in urine; cefpirome plasma protein binding was 10% (paper Monte Carlo",
      "simulation paragraph, citing references 53, 63). A 10-min zero-order",
      "infusion duration TK0 was fixed structurally; simulate this model by",
      "supplying a dose with rate = amt / (10/60) (10 minutes expressed in",
      "hours) or dur = 10/60 in event records."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # STRUCTURAL CLEARANCE AND VOLUME PARAMETERS -- Bulitta 2011 Table 3
    # (NONMEM LBM-allometric column). Typical values are reported separately
    # for CF patients and healthy volunteers (HVs); the values used here are
    # the HV-column point estimates (Table 3 last column under the NONMEM
    # heading), because typicals are anchored to the DIS_CF = 0 reference per
    # the DIS_CF canonical convention registered alongside this model file.
    # CF values are recovered at DIS_CF = 1 via the e_cf_* covariate effects.
    # ------------------------------------------------------------------------
    lcl_renal  <- log(5.78);  label("Renal CL at DIS_CF = 0 (HV typical, L/h)")
    # Bulitta 2011 Table 3 NONMEM HV column: CL_R = 5.78 L/h
    lcl_nonren <- log(1.17);  label("Non-renal CL at DIS_CF = 0 (HV typical, L/h)")
    # Bulitta 2011 Table 3 NONMEM HV column: CL_NR = 1.17 L/h
    lvc        <- log(7.44);  label("Central volume V1 at DIS_CF = 0 (HV typical, L)")
    # Bulitta 2011 Table 3 NONMEM HV column: V1 = 7.44 L
    lvp        <- log(2.67);  label("Shallow peripheral volume V2 at DIS_CF = 0 (HV typical, L)")
    # Bulitta 2011 Table 3 NONMEM HV column: V2 = 2.67 L
    lvp2       <- log(4.63);  label("Deep peripheral volume V3 at DIS_CF = 0 (HV typical, L)")
    # Bulitta 2011 Table 3 NONMEM HV column: V3 = 4.63 L
    lq         <- log(25.6);  label("Inter-compartmental CL central <-> shallow peripheral, Q12 (L/h)")
    # Bulitta 2011 Table 3 NONMEM (single value, no CF/HV split): CLic_shallow = 25.6 L/h
    lq2        <- log(2.81);  label("Inter-compartmental CL central <-> deep peripheral, Q23 (L/h)")
    # Bulitta 2011 Table 3 NONMEM (single value, no CF/HV split): CLic_deep = 2.81 L/h

    # ------------------------------------------------------------------------
    # ALLOMETRIC SCALING -- Bulitta 2011 Population PK analysis, Body size
    # paragraph: "We set the allometric exponent at 0.75 for all clearance
    # terms and at 1.0 for all volumes." Reference LBM = 53 kg.
    # ------------------------------------------------------------------------
    e_lbm_cl_q  <- fixed(0.75);  label("Allometric (LBM) exponent on CL_R, CL_NR, Q12, Q23 (unitless)")
    # Bulitta 2011 Body size paragraph (allometric scaling on LBM)
    e_lbm_vc_vp <- fixed(1.00);  label("Allometric (LBM) exponent on V1, V2, V3 (unitless)")
    # Bulitta 2011 Body size paragraph (allometric scaling on LBM)

    # ------------------------------------------------------------------------
    # DISEASE-SPECIFIC SCALE FACTORS (DIS_CF effects) -- Bulitta 2011 Table 4
    # LBM-allometric row, encoded as log-scale effects on the HV typical:
    #   FCYF_CLR  = 1.07 (95% CI 0.97-1.20) -> e_cf_cl_renal  = log(1.07) =  0.06766
    #   FCYF_CLNR = 1.13 (95% CI 0.79-1.57) -> e_cf_cl_nonren = log(1.13) =  0.12222
    #   FCYF_VSS  = 0.98 (95% CI 0.90-1.06) -> e_cf_vc_vp_vp2 = log(0.98) = -0.02020
    # FCYF_VSS is applied uniformly to V1, V2, V3 (the paper estimates a
    # single FCYF for volume of distribution at steady state); the convention
    # for shared exponents across multiple volumes is e_<cov>_vc_vp_vp2 by
    # analogy with the existing two-volume e_<cov>_vc_vp form.
    # ------------------------------------------------------------------------
    e_cf_cl_renal   <-  0.06766;  label("Log-scale FCYF_CLR effect of CF status on CL_R (unitless)")
    # Bulitta 2011 Table 4 LBM-allometric row: FCYF_CLR = 1.07
    e_cf_cl_nonren  <-  0.12222;  label("Log-scale FCYF_CLNR effect of CF status on CL_NR (unitless)")
    # Bulitta 2011 Table 4 LBM-allometric row: FCYF_CLNR = 1.13
    e_cf_vc_vp_vp2  <- -0.02020;  label("Log-scale FCYF_VSS effect of CF status on V1, V2, V3 (unitless)")
    # Bulitta 2011 Table 4 LBM-allometric row: FCYF_VSS = 0.98

    # ------------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Bulitta 2011 Table 3 NONMEM column,
    # BSV section. Footnote c clarifies that the reported BSV values are
    # "apparent coefficients of variation" i.e. sqrt(omega^2) on the log
    # scale; the back-transform is omega^2 = (CV)^2 (not log(1 + CV^2)):
    #   CL_R       : 0.135 -> omega^2 = 0.135^2 = 0.018225
    #   CL_NR      : 0.322 -> omega^2 = 0.322^2 = 0.103684
    #   V1         : 0.448 -> omega^2 = 0.448^2 = 0.200704
    #   V2         : 0.735 -> omega^2 = 0.735^2 = 0.540225
    #   V3         : 0.111 -> omega^2 = 0.111^2 = 0.012321
    # The V1 / V2 random effects are correlated (Bulitta 2011 Table 3
    # footnote d: r(V1, V2) = -0.838 in NONMEM), so they enter as a block:
    #   cov(eta_V1, eta_V2) = -0.838 * 0.448 * 0.735 = -0.275893
    # Inter-compartmental CL BSVs were fixed to 0 in the NONMEM final model
    # (Bulitta 2011 Table 3 NONMEM column for BSV(CLic_shallow) and
    # BSV(CLic_deep)), so etalq and etalq2 are not declared here.
    # ------------------------------------------------------------------------
    etalcl_renal             ~ 0.018225  # Bulitta 2011 Table 3 NONMEM: BSV(CL_R) = 0.135 apparent CV
    etalcl_nonren            ~ 0.103684  # Bulitta 2011 Table 3 NONMEM: BSV(CL_NR) = 0.322 apparent CV
    etalvp2                  ~ 0.012321  # Bulitta 2011 Table 3 NONMEM: BSV(V3) = 0.111 apparent CV
    etalvc + etalvp ~ c(0.200704,
                       -0.275893, 0.540225)
    # Bulitta 2011 Table 3 NONMEM: BSV(V1) = 0.448, BSV(V2) = 0.735, r(V1,V2) = -0.838

    # ------------------------------------------------------------------------
    # RESIDUAL ERROR -- Bulitta 2011 Table 3 NONMEM column. The paper used a
    # combined additive + proportional model for plasma concentrations and a
    # combined additive + proportional model for amounts in urine.
    # ------------------------------------------------------------------------
    propSd            <- 0.0875;  label("Proportional plasma residual SD (fraction)")
    # Bulitta 2011 Table 3 NONMEM: CV_C = 8.75% for plasma concentrations
    addSd             <- 0.429;   label("Additive plasma residual SD (mg/L)")
    # Bulitta 2011 Table 3 NONMEM: SD_C = 0.429 mg/L for plasma concentrations
    propSd_urineAmt   <- 0.182;   label("Proportional urine cumulative-amount residual SD (fraction)")
    # Bulitta 2011 Table 3 NONMEM: CV_AU = 18.2% for amounts in urine
    addSd_urineAmt    <- 6.12;    label("Additive urine cumulative-amount residual SD (mg)")
    # Bulitta 2011 Table 3 NONMEM: SD_AU = 6.12 mg for amounts in urine
  })

  model({
    # ---------------------------------------------------------------------
    # Reference values.
    # ---------------------------------------------------------------------
    lbm_ref <- 53   # kg; reference LBM (Bulitta 2011 Body size paragraph)

    # ---------------------------------------------------------------------
    # Individual PK parameters with covariate effects.
    # CF effects enter on the log scale via DIS_CF (= 1 for CF, 0 for HV);
    # at DIS_CF = 0 the typical values reduce to the HV column of Table 3.
    # Allometric scaling on LBM with fixed exponents per the paper.
    # ---------------------------------------------------------------------
    cl_renal  <- exp(lcl_renal  + etalcl_renal  + e_cf_cl_renal  * DIS_CF) *
                 (LBM / lbm_ref) ^ e_lbm_cl_q
    cl_nonren <- exp(lcl_nonren + etalcl_nonren + e_cf_cl_nonren * DIS_CF) *
                 (LBM / lbm_ref) ^ e_lbm_cl_q
    vc        <- exp(lvc  + etalvc  + e_cf_vc_vp_vp2 * DIS_CF) *
                 (LBM / lbm_ref) ^ e_lbm_vc_vp
    vp        <- exp(lvp  + etalvp  + e_cf_vc_vp_vp2 * DIS_CF) *
                 (LBM / lbm_ref) ^ e_lbm_vc_vp
    vp2       <- exp(lvp2 + etalvp2 + e_cf_vc_vp_vp2 * DIS_CF) *
                 (LBM / lbm_ref) ^ e_lbm_vc_vp
    q         <- exp(lq)  * (LBM / lbm_ref) ^ e_lbm_cl_q
    q2        <- exp(lq2) * (LBM / lbm_ref) ^ e_lbm_cl_q

    # ---------------------------------------------------------------------
    # Micro-constants.
    # ---------------------------------------------------------------------
    kel_renal  <- cl_renal  / vc
    kel_nonren <- cl_nonren / vc
    k12        <- q  / vc
    k21        <- q  / vp
    k13        <- q2 / vc
    k31        <- q2 / vp2

    # ---------------------------------------------------------------------
    # 3-compartment ODE system with parallel renal + non-renal elimination
    # and cumulative urinary excretion (mass entering urine = CL_R * Cc).
    # IV doses target central directly; users supply rate (e.g. amt / (10/60))
    # or dur (10/60 h) to encode the 10-min zero-order infusion.
    # ---------------------------------------------------------------------
    d/dt(central)     <- -(kel_renal + kel_nonren) * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-   k13 * central - k31 * peripheral2
    d/dt(urine)       <-   kel_renal * central

    # ---------------------------------------------------------------------
    # Observations.
    # Cc      : plasma concentration in mg/L (dose in mg, vc in L -> mg/L).
    # urineAmt: cumulative amount excreted unchanged in urine, in mg.
    # ---------------------------------------------------------------------
    Cc       <- central / vc
    urineAmt <- urine

    Cc       ~ add(addSd)          + prop(propSd)
    urineAmt ~ add(addSd_urineAmt) + prop(propSd_urineAmt)
  })
}
