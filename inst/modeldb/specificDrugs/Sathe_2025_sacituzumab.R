Sathe_2025_sacituzumab <- function() {
  description <- "Updated coupled three-analyte population PK model for sacituzumab govitecan (SG, the ADC; output Cc), free SN-38 (released payload; output Cc_sn38), and total antibody (tAB; output Cc_tab) after pooling data from TROPiCS-02 (HR+/HER2- mBC) with IMMU-132-01 and ASCENT (mTNBC + mUC + HR+/HER2- mBC + other solid tumors) (Sathe 2025). Structure is identical to Sathe 2024 (see Sathe_2024_sacituzumab): three two-compartment models with body-weight allometric scaling; SG has IIV on CL and a baseline-albumin power covariate on CL; free SN-38 is generated from SG central by a first-order release rate KREL with apparent V1 and V2 fixed to literature values; tAB has time-dependent CL (asymptotic onset), correlated IIV on CL and V1, and covariates of baseline albumin (CL), tumor type (CL), and sex (V1). Parameter values are the pooled 3-study updated estimates from Sathe 2025 Tables 2, 3, and 4. Simulation requires dosing two compartments simultaneously (central and central_tab) for each SG infusion event."
  reference <- "Sathe AG, Jones AK, Diderichsen PM, Wang X, Chang P, Verret W, Girish S. Sacituzumab Govitecan Population Pharmacokinetics: Updated Analyses Using HR+/HER2- Metastatic Breast Cancer Data From the Phase 3 TROPiCS-02 Trial. Clin Transl Sci. 2025;18(8):e70291. doi:10.1111/cts.70291"
  vignette <- "Sathe_2025_sacituzumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight used for allometric scaling of CL/Q and V1/V2 of SG, free SN-38, and tAB; reference 70 kg (NONMEM control stream TVWT constant; pooled analysis median 69 kg).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum albumin used as a power covariate on SG CL (exponent -0.395; Sathe 2025 Table 2) and tAB CL (exponent -0.734; Sathe 2025 Table 4); reference 38 g/L (NONMEM control stream). Source paper column BALB (baseline albumin) mapped to canonical ALB.",
      source_name        = "BALB"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Sathe 2025 dataset uses SEXF directly. Effect applies on tAB V1 only: V1 is +15.3% in males (SEXF = 0) relative to females (the reference) in Sathe 2025 Table 4 (was +12.1% in Sathe 2024). Effect coded as `1 + e_sex_vc_tab * (1 - SEXF)` to match the source's male-deviation parameterization.",
      source_name        = "SEXF"
    ),
    TUMTP_OTHER = list(
      description        = "Tumor-type indicator: 1 = 'Other' epithelial cancer (small-cell and non-small-cell lung cancer, colorectal, esophageal, pancreatic ductal adenocarcinoma, etc.); 0 = mTNBC, mUC, or HR+/HER2- mBC",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (mTNBC, mUC, or HR+/HER2- mBC)",
      notes              = "Effect on tAB CL only: -11.2% CL when TUMTP_OTHER = 1 (Sathe 2025 Table 4; was -13.4% in Sathe 2024). Source column PAT2 takes integer levels (1 = mTNBC, 2 = mUC or HR+/HER2- mBC, 4 = Other) and the source NONMEM control stream collapses PAT2 = 1 and PAT2 = 2 into the reference (no effect). 'Other' composition (unchanged from Sathe 2024) = small-cell lung cancer, non-small-cell lung cancer, colorectal, esophageal, pancreatic ductal adenocarcinoma, etc.; all 184 patients from IMMU-132-01 (Sathe 2025 Table 1). Scope is paper-specific.",
      source_name        = "PAT2 (recoded: 4 -> 1, 1/2 -> 0)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age (years); tested but not retained in any of the three final models (Sathe 2025 Section 3.6, Figure 3, Figures S4-S5).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in forest-plot covariate assessment; not statistically or clinically significant."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault, mL/min); tested but not retained (mild-to-moderate renal impairment did not have a clinically relevant impact; severe RI n=4 excluded).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the covariate assessment; no dose adjustment recommended based on renal function."
    ),
    UGT1A1 = list(
      description = "UGT1A1 genotype (*1/*1, *1/*28, *28/*28, other); tested but not retained (Sathe 2025 Section 3.6).",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Model-predicted exposures for *1/*1, *1/*28, and *28/*28 genotypes were comparable."
    ),
    TROP2 = list(
      description = "Trop-2 expression (categorical staining in IMMU-132-01; H-score in ASCENT and TROPiCS-02); tested but not retained.",
      units       = "(mixed)",
      type        = "mixed",
      notes       = "Not a statistically significant covariate; evaluated as both continuous and categorical."
    )
  )

  population <- list(
    n_subjects     = 789,
    n_studies      = 3,
    age_range      = "27-88 years",
    age_median     = "58 years",
    weight_range   = "31-140 kg",
    weight_median  = "69 kg",
    sex_female_pct = 85,
    race_ethnicity = c(White = 79),
    disease_state  = "Pooled solid tumors: metastatic triple-negative breast cancer (mTNBC, n = 277), metastatic urothelial cancer (mUC, n = 36), HR+/HER2- metastatic breast cancer (n = 292), and other solid tumors (n = 184; small-cell and non-small-cell lung cancer, colorectal, esophageal, pancreatic ductal adenocarcinoma, etc.).",
    dose_range     = "IMMU-132-01: 8, 10, 12, or 18 mg/kg IV on days 1 and 8 of 21-day cycles. ASCENT and TROPiCS-02: 10 mg/kg IV on days 1 and 8 of 21-day cycles.",
    regions        = "Multinational pooled analysis (IMMU-132-01 + ASCENT + TROPiCS-02)",
    studies        = "IMMU-132-01 (NCT01631552, n = 276: mTNBC = 24, mUC = 36, HR+/HER2- mBC = 32, other = 184), ASCENT (NCT02574455, n = 253; all mTNBC), and TROPiCS-02 (NCT03901339, n = 260; all HR+/HER2- mBC).",
    baseline_albumin_median = "39 g/L (range 19-51)",
    baseline_clcr_median    = "91 mL/min (range 22-262); 51% normal (>= 90), 38% mild impairment (60 to < 90), 10% moderate impairment (30 to < 60), <1% severe impairment.",
    ecog_ps        = "0: 38%, 1: 61%, 2: 0.4% (Sathe 2025 Table 1).",
    ugt1a1_distribution = "*1/*1: 38%, *1/*28: 39%, *28/*28: 11%, Other: 1%, Missing: 10% (Sathe 2025 Table 1).",
    ada_positive_pct   = "<1% (4/789 patients from ASCENT; 0/260 from TROPiCS-02; 0/276 from IMMU-132-01).",
    notes          = "Sacituzumab govitecan (SG) is an ADC of an anti-Trop-2 humanized monoclonal antibody (hRS7) covalently linked to the topoisomerase 1 inhibitor SN-38 via a hydrolyzable linker (average drug-to-antibody ratio = 8). Structural models were not re-assessed in Sathe 2025 (previous Sathe 2024 structure carried forward; only pooled 3-study parameter re-estimation). Demographics summary: Sathe 2025 Section 3.1 and Table 1."
  )

  ini({
    # ============================================================
    # SG (sacituzumab govitecan, the ADC) - Sathe 2025 Table 2
    # (updated final PopPK model using pooled IMMU-132-01 + ASCENT
    # + TROPiCS-02 data). Values listed here are the paper-Table 2
    # point estimates; the underlying NONMEM control stream (Sathe
    # 2025 supplement Section a) uses initial values before the
    # combined-dataset re-estimation.
    # ============================================================
    lcl    <- log(0.128);   label("SG clearance (L/h)")                                  # Sathe 2025 Table 2 / Table S2 (CL_SG = 0.128)
    lvc    <- log(2.65);    label("SG central volume (L)")                                # Sathe 2025 Table 2 / Table S2 (V1_SG = 2.65)
    lq     <- log(0.00513); label("SG intercompartmental clearance (L/h)")                # Sathe 2025 Table 2 / Table S2 (Q_SG = 0.00513)
    lvp    <- log(0.929);   label("SG peripheral volume (L)")                             # Sathe 2025 Table 2 / Table S2 (V2_SG = 0.929)
    e_wt_cl_q <-  0.523; label("Body-weight allometric exponent on SG CL and Q (unitless)")   # Sathe 2025 Table 2 / Table S2 (0.523)
    e_wt_vc_vp <-  0.540; label("Body-weight allometric exponent on SG V1 and V2 (unitless)") # Sathe 2025 Table 2 / Table S2 (0.540)
    e_alb_cl <- -0.395; label("Baseline-albumin power exponent on SG CL (unitless)")          # Sathe 2025 Table 2 / Table S2 (-0.395)

    # ============================================================
    # Free SN-38 (released payload, sequential to SG) - Sathe 2025 Table 3
    # Log-scale estimates as reported; untransformed values shown
    # in labels (Table 3 "Untransformed estimate" column).
    # ============================================================
    lkrel    <- -2.37; label("First-order SG-to-free-SN-38 release rate (log 1/h); KREL = 0.0937 1/h")             # Sathe 2025 Table 3 / Table S3
    lcl_sn38  <-  5.99; label("Apparent SN-38 clearance (log L/h); CLSN38/F = 401 L/h")                             # Sathe 2025 Table 3 / Table S3
    lq_sn38   <-  5.49; label("Apparent SN-38 intercompartmental clearance (log L/h); QSN38/F = 243 L/h")           # Sathe 2025 Table 3 / Table S3
    lvc_sn38  <- fixed(log(49));   label("Apparent SN-38 central volume V1SN38/F (L); FIXED to Sathe 2024 [ref 10] literature value")     # Sathe 2025 Table 3; value carried from Sathe 2024 [ref 10] (originally Klein et al. Clin Pharmacol Ther 2002;72:638-647)
    lvp_sn38  <- fixed(log(2177)); label("Apparent SN-38 peripheral volume V2SN38/F (L); FIXED to Sathe 2024 [ref 10] literature value")  # Sathe 2025 Table 3; value carried from Sathe 2024 [ref 10] (originally Klein et al. Clin Pharmacol Ther 2002;72:638-647)
    e_wt_cl_q_sn38 <- 0.519; label("Body-weight allometric exponent on free-SN-38 CL and Q (unitless)")             # Sathe 2025 Table 3 / Table S3

    # ============================================================
    # Total antibody (tAB) - Sathe 2025 Table 4
    # NB Sathe 2025 Table S4 footnote c: CL_tAB and Q_tAB
    # parameters were reported multiplied by 100 in the estimation
    # scale (1.55 -> 0.0155 L/h; 1.05 -> 0.0105 L/h). The values
    # stored here are the actual (post-division) values.
    # ============================================================
    lcl_tab    <- log(0.0155); label("tAB clearance, baseline at t = 0 (L/h)")                                          # Sathe 2025 Table 4
    lvc_tab    <- log(2.97);   label("tAB central volume (L)")                                                          # Sathe 2025 Table 4
    lq_tab     <- log(0.0105); label("tAB intercompartmental clearance (L/h)")                                          # Sathe 2025 Table 4
    lvp_tab    <- log(1.32);   label("tAB peripheral volume (L)")                                                       # Sathe 2025 Table 4
    e_wt_cl_q_tab  <- 0.422; label("Body-weight allometric exponent on tAB CL and Q (unitless)")                        # Sathe 2025 Table 4
    e_wt_vc_vp_tab <- 0.458; label("Body-weight allometric exponent on tAB V1 and V2 (unitless)")                       # Sathe 2025 Table 4
    e_alb_cl_tab   <- -0.734; label("Baseline-albumin power exponent on tAB CL (unitless)")                             # Sathe 2025 Table 4
    e_tumor_cl_tab <- -0.112; label("Tumor-type 'Other' fractional change on tAB CL (multiplicative; unitless)")        # Sathe 2025 Table 4
    e_sex_vc_tab   <-  0.153; label("Male-sex fractional change on tAB V1 (multiplicative; applied when SEXF = 0; unitless)") # Sathe 2025 Table 4
    maxRed_tab <- 16.8;    label("tAB CL maximum relative reduction (%) at t -> infinity")                              # Sathe 2025 Table 4 (16.8%; half-time ~74 days)
    keff_tab   <- 3.91e-4; label("tAB CL time-decline rate constant (1/h); half-time ~74 days")                         # Sathe 2025 Table 4 (log rate constant = -7.85 -> exp(-7.85) = 3.91e-4; log(2)/3.91e-4 h = 74 days)

    # ============================================================
    # Inter-individual variability
    # ============================================================
    etalcl ~ 0.0136                                            # Sathe 2025 Table 2 / Table S2: IIV variance for CL_SG = 0.0136 (CV ~ 11.7%)
    etalkrel + etalcl_sn38 ~ c(0.397,
                              0.406, 0.630)                     # Sathe 2025 Table 3 / Table S3 BLOCK(2): var(KREL)=0.397, cov(KREL, CLSN38/F)=0.406, var(CLSN38/F)=0.630
    etalcl_tab + etalvc_tab ~ c(0.110,
                              0.0390, 0.0397)                   # Sathe 2025 Table 4 / Table S4 BLOCK(2): var(CL_tAB)=0.110, cov(CL_tAB, V1_tAB)=0.0390, var(V1_tAB)=0.0397

    # ============================================================
    # Residual error
    # In NONMEM: SG and free SN-38 are additive on log(DV)
    # (= proportional in linear space); tAB is linear additive +
    # proportional. Time-after-last-dose and study-indicator
    # effects on residual variance are reported but omitted from
    # the simulation model (see vignette Assumptions and deviations).
    # ============================================================
    propSd        <- 0.198; label("SG proportional residual SD on log scale (Sathe 2025 Table 2: RUV SD on log SG)")
    propSd_sn38   <- 0.344; label("Free SN-38 proportional residual SD on log scale (Sathe 2025 Table 3: exp(-1.07) = 0.344)")
    addSd_tab     <- 21.8;  label("tAB additive residual SD (ug/mL; Sathe 2025 Table 4)")
    propSd_tab    <- 0.191; label("tAB proportional residual SD as fraction (Sathe 2025 Table 4)")
  })
  model({
    # ------------------------------------------------------------
    # SG (sacituzumab govitecan, the ADC)
    # ------------------------------------------------------------
    bwt_cl_factor <- (WT / 70) ^ e_wt_cl_q
    bwt_v_factor  <- (WT / 70) ^ e_wt_vc_vp
    alb_cl_factor <- (ALB / 38) ^ e_alb_cl

    cl <- exp(lcl + etalcl) * bwt_cl_factor * alb_cl_factor
    vc <- exp(lvc)           * bwt_v_factor
    q  <- exp(lq)            * bwt_cl_factor
    vp <- exp(lvp)           * bwt_v_factor

    ke  <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------
    # Free SN-38 (sequential to SG via first-order release KREL)
    # ------------------------------------------------------------
    bwt_cl_sn38_factor <- (WT / 70) ^ e_wt_cl_q_sn38

    krel     <- exp(lkrel + etalkrel)
    cl_sn38  <- exp(lcl_sn38 + etalcl_sn38) * bwt_cl_sn38_factor
    q_sn38   <- exp(lq_sn38)                * bwt_cl_sn38_factor
    vc_sn38  <- exp(lvc_sn38)
    vp_sn38  <- exp(lvp_sn38)

    ke_sn38  <- cl_sn38 / vc_sn38
    k12_sn38 <- q_sn38  / vc_sn38
    k21_sn38 <- q_sn38  / vp_sn38

    # ------------------------------------------------------------
    # Total antibody (tAB), with time-dependent CL
    # td_factor goes from 1 at t = 0 to (1 - maxRed_tab/100) at
    # t -> infinity with rate constant keff_tab.
    # ------------------------------------------------------------
    bwt_cl_tab_factor <- (WT / 70) ^ e_wt_cl_q_tab
    bwt_v_tab_factor  <- (WT / 70) ^ e_wt_vc_vp_tab
    alb_cl_tab_factor <- (ALB / 38) ^ e_alb_cl_tab
    tumor_cl_tab_factor <- 1 + e_tumor_cl_tab * TUMTP_OTHER
    sex_vc_tab_factor   <- 1 + e_sex_vc_tab   * (1 - SEXF)
    td_cl_tab_factor    <- 1 - (maxRed_tab / 100) * (1 - exp(-keff_tab * t))

    cl_tab <- exp(lcl_tab + etalcl_tab) * bwt_cl_tab_factor * alb_cl_tab_factor *
             tumor_cl_tab_factor * td_cl_tab_factor
    vc_tab <- exp(lvc_tab + etalvc_tab) * bwt_v_tab_factor  * sex_vc_tab_factor
    q_tab  <- exp(lq_tab)               * bwt_cl_tab_factor
    vp_tab <- exp(lvp_tab)              * bwt_v_tab_factor

    ke_tab  <- cl_tab / vc_tab
    k12_tab <- q_tab  / vc_tab
    k21_tab <- q_tab  / vp_tab

    # ------------------------------------------------------------
    # ODE system (named compartments; user doses central and
    # central_tab simultaneously to simulate one SG infusion).
    # KREL drives free-SN-38 generation from SG central without a
    # back-coupling sink on SG, matching Sathe 2025 supplement
    # Section b $DES (identical to Sathe 2024).
    # ------------------------------------------------------------
    d/dt(central)       <- -ke  * central - k12  * central + k21  * peripheral1
    d/dt(peripheral1)   <-  k12 * central - k21  * peripheral1

    d/dt(central_sn38)     <-  krel     * central     - ke_sn38 * central_sn38 -
                              k12_sn38 * central_sn38 + k21_sn38 * peripheral1_sn38
    d/dt(peripheral1_sn38) <-  k12_sn38 * central_sn38 - k21_sn38 * peripheral1_sn38

    d/dt(central_tab)      <- -ke_tab   * central_tab - k12_tab  * central_tab + k21_tab * peripheral1_tab
    d/dt(peripheral1_tab)  <-  k12_tab  * central_tab - k21_tab  * peripheral1_tab

    # ------------------------------------------------------------
    # Observations
    # Concentrations expressed in ug/mL (= mg/L) given dose in mg
    # and volumes in L. Sathe 2025 reports SG and free SN-38 in
    # ng/mL and tAB in ug/mL; multiply Cc and Cc_sn38 by 1000 to
    # match the paper's ng/mL units when comparing.
    # ------------------------------------------------------------
    Cc      <- central       / vc
    Cc_sn38 <- central_sn38  / vc_sn38
    Cc_tab  <- central_tab   / vc_tab

    Cc      ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Cc_tab  ~ add(addSd_tab) + prop(propSd_tab)
  })
}
