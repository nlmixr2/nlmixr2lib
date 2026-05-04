Sathe_2024_sacituzumab <- function() {
  description <- "Coupled three-analyte population PK model for sacituzumab govitecan (SG, the ADC; output Cc), free SN-38 (released payload; output Cc_sn38), and total antibody (tAB; output Cc_tab) in adults with metastatic triple-negative breast cancer and other solid tumors (Sathe 2024). All three analytes are described by two-compartment models with body-weight allometric scaling. SG carries IIV on CL and a baseline-albumin power covariate on CL. Free SN-38 is generated from SG by a first-order release rate KREL with apparent volumes fixed to literature values (Klein 2002). tAB has time-dependent CL (asymptotic onset, max ~17% reduction at t1/2 ~48 days), correlated IIV on CL and V1, and covariates of baseline albumin (CL), tumor type (CL), and sex (V1). Simulation requires dosing two compartments simultaneously (central and central_tab) for each SG infusion event."
  reference <- "Sathe AG, Singh I, Singh P, Diderichsen PM, Wang X, Chang P, Taqui A, Phan S, Girish S, Othman AA. Population Pharmacokinetics of Sacituzumab Govitecan in Patients with Metastatic Triple-Negative Breast Cancer and Other Solid Tumors. Clin Pharmacokinet. 2024;63(5):669-681. doi:10.1007/s40262-024-01366-3"
  vignette <- "Sathe_2024_sacituzumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight used for allometric scaling of CL/Q and V1/V2 of SG, free SN-38, and tAB; reference 70 kg.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum albumin used as a power covariate on SG CL (exponent -0.355) and tAB CL (exponent -0.735); reference 38 g/L. Source paper uses 'BALB' (baseline albumin) for the column name; mapped to canonical ALB.",
      source_name        = "BALB"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Source paper analysis dataset uses SEXF directly. Effect applies on tAB V1 only: V1 is +12.1% in males (SEXF = 0) relative to females (the reference). Effect coded as `1 + e_sex_vc_tab * (1 - SEXF)` to match the source's male-deviation parameterization.",
      source_name        = "SEXF"
    ),
    TUMTP_OTH = list(
      description        = "Tumor-type indicator: 1 = 'Other' epithelial cancer (NSCLC, SCLC, colorectal, esophageal, pancreatic ductal adenocarcinoma, etc.); 0 = mTNBC, mUC, or HR+/HER2- mBC",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (mTNBC, mUC, or HR+/HER2- mBC)",
      notes              = "Effect on tAB CL only: -13.4% CL when TUMTP_OTH = 1 (Sathe 2024 Table 3). Source column PAT2 takes integer levels (1 = mTNBC, 2 = mUC or HR+/HER2- mBC, 4 = Other) and the source NONMEM control stream collapses PAT2 = 1 and PAT2 = 2 into the reference (no effect). Composition of 'Other' in this analysis = small-cell lung cancer, non-small-cell lung cancer, colorectal cancer, esophageal cancer, pancreatic ductal adenocarcinoma, etc. (Sathe 2024 Methods 2.4 and Table S2). Scope is paper-specific: a different paper's 'Other' tumor pool is not interchangeable.",
      source_name        = "PAT2 (recoded: 4 -> 1, 1/2 -> 0)"
    )
  )

  population <- list(
    n_subjects     = 529,
    n_studies      = 2,
    age_range      = "27-88 years",
    age_median     = "58 years",
    weight_range   = "37-140 kg",
    weight_median  = "70 kg",
    sex_female_pct = 78,
    race_ethnicity = c(White = 84, Black_or_African_American = 7, Asian = 2, Other = 6),
    disease_state  = "Metastatic triple-negative breast cancer (mTNBC, n = 277), metastatic urothelial cancer (mUC, n = 36), HR+/HER2- metastatic breast cancer (n = 32), and other solid tumors (n = 184; epithelial cancers including small-cell lung, non-small-cell lung, colorectal, esophageal, pancreatic ductal adenocarcinoma, etc.).",
    dose_range     = "Phase I/II (IMMU-132-01): 8, 10, 12, or 18 mg/kg IV on days 1 and 8 of 21-day cycles. Phase III (ASCENT): 10 mg/kg IV on days 1 and 8 of 21-day cycles.",
    regions        = "Multinational pooled analysis (IMMU-132-01 + ASCENT)",
    studies        = "IMMU-132-01 (NCT01631552, n = 276 of which mTNBC = 24, mUC = 36, HR+/HER2- mBC = 32, other = 184) and ASCENT (NCT02574455, n = 253; all mTNBC).",
    baseline_albumin_median = "38 g/L (range 19-50)",
    baseline_clcr_median    = "91 mL/min (range 22-262); 51% normal (>= 90), 38% mild impairment (60 to < 90), 11% moderate impairment (30 to < 60), <1% severe impairment.",
    notes          = "Sacituzumab govitecan is an ADC of an anti-Trop-2 humanized monoclonal antibody (hRS7) covalently linked to the topoisomerase 1 inhibitor SN-38 via a hydrolyzable linker (average drug-to-antibody ratio = 8). Per Sathe 2024 Methods 2.3, the tAB model was fit using the SG dose directly, then the 'true' tAB CL/V were obtained by multiplying the model estimates by 0.92 (mass ratio of naked antibody to SG). The CL/V values stored in this model are the model-fit estimates intended to be used with the SG dose, not the 0.92-scaled 'true' values. Demographics summary: Sathe 2024 Section 3.1 and Table S2.",
    dosing_note    = "Each SG infusion event generates inputs to both the SG compartment and the tAB compartment. To simulate, provide TWO dose events per infusion: one with cmt = 'central' and one with cmt = 'central_tab', sharing the same amt, rate (or dur), time, and ii/addl. The vignette demonstrates this convention."
  )

  ini({
    # ============================================================
    # SG (sacituzumab govitecan, the ADC) - Sathe 2024 Table 1
    # Final estimates listed in the supplement (Section 6, Free
    # SN-38 control stream, $THETA FIX block) carry more precision
    # than the rounded paper-table values; both agree.
    # ============================================================
    lcl    <- log(0.13333);   label("SG clearance (L/h)")                            # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(1) FIX
    lvc    <- log(2.77011);   label("SG central volume (L)")                          # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(2) FIX
    lq     <- log(0.0055141); label("SG intercompartmental clearance (L/h)")          # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(3) FIX
    lvp    <- log(0.907787);  label("SG peripheral volume (L)")                       # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(4) FIX
    allo_cl <-  0.507571; label("Body-weight allometric exponent on SG CL and Q (unitless)")  # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(6) FIX
    allo_v  <-  0.532396; label("Body-weight allometric exponent on SG V1 and V2 (unitless)") # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(7) FIX
    e_alb_cl <- -0.355253; label("Baseline-albumin power exponent on SG CL (unitless)")        # Sathe 2024 Table 1; supplement Section 6.b $THETA TH(8) FIX

    # ============================================================
    # Free SN-38 (released payload, sequential to SG) - Sathe 2024 Table 2
    # The paper reports KREL, CLSN38/F, and QSN38/F as log-scale
    # estimates; the values below are the log estimates exactly as
    # reported. Untransformed values (KREL = 0.0961 1/h,
    # CLSN38/F = 409 L/h, QSN38/F = 247 L/h) are noted in labels.
    # ============================================================
    lkrel    <- -2.34; label("First-order SG-to-free-SN-38 release rate (log 1/h); KREL = 0.0961 1/h")             # Sathe 2024 Table 2
    lcl_sn38  <-  6.02; label("Apparent SN-38 clearance (log L/h); CLSN38/F = 409 L/h")                              # Sathe 2024 Table 2
    lq_sn38   <-  5.51; label("Apparent SN-38 intercompartmental clearance (log L/h); QSN38/F = 247 L/h")            # Sathe 2024 Table 2
    lvc_sn38  <- log(49);   label("Apparent SN-38 central volume V1SN38/F (L); FIXED to Klein 2002 literature value")    # Sathe 2024 Table 2; literature ref [19] = Klein et al. Clin Pharmacol Ther 2002;72:638-647
    lvp_sn38  <- log(2177); label("Apparent SN-38 peripheral volume V2SN38/F (L); FIXED to Klein 2002 literature value") # Sathe 2024 Table 2; literature ref [19] = Klein et al. Clin Pharmacol Ther 2002;72:638-647
    allo_cl_sn38 <- 0.500; label("Body-weight allometric exponent on free-SN-38 CL and Q (unitless)") # Sathe 2024 Table 2

    # ============================================================
    # Total antibody (tAB) - Sathe 2024 Table 3
    # Values stored here are the model-fit estimates (used with
    # the SG dose, see population$notes); 'true' tAB CL/V are
    # obtained by multiplying by 0.92 per Sathe 2024 Methods 2.3.
    # ============================================================
    lcl_tab    <- log(0.016); label("tAB clearance, baseline at t = 0 (L/h)")                                          # Sathe 2024 Table 3
    lvc_tab    <- log(3.06);  label("tAB central volume (L)")                                                           # Sathe 2024 Table 3
    lq_tab     <- log(0.010); label("tAB intercompartmental clearance (L/h)")                                           # Sathe 2024 Table 3
    lvp_tab    <- log(1.20);  label("tAB peripheral volume (L)")                                                        # Sathe 2024 Table 3
    allo_cl_tab <- 0.372; label("Body-weight allometric exponent on tAB CL and Q (unitless)")                            # Sathe 2024 Table 3
    allo_v_tab  <- 0.446; label("Body-weight allometric exponent on tAB V1 and V2 (unitless)")                           # Sathe 2024 Table 3
    e_alb_cl_tab   <- -0.735; label("Baseline-albumin power exponent on tAB CL (unitless)")                              # Sathe 2024 Table 3
    e_tumor_cl_tab <- -0.134; label("Tumor-type 'Other' fractional change on tAB CL (multiplicative; unitless)")          # Sathe 2024 Table 3
    e_sex_vc_tab    <-  0.121; label("Male-sex fractional change on tAB V1 (multiplicative; applied when SEXF = 0; unitless)") # Sathe 2024 Table 3
    maxRed_tab <- 16.7;    label("tAB CL maximum relative reduction (%) at t -> infinity")                              # Sathe 2024 Table 3
    keff_tab   <- 6.08e-4; label("tAB CL time-decline rate constant (1/h); half-time ~48 days")                          # Sathe 2024 Table 3

    # ============================================================
    # Inter-individual variability
    # ============================================================
    etalcl ~ 0.011                                            # Sathe 2024 Table 1: IIV variance for CL_SG = 0.011 (CV ~ 10.5%)
    etalkrel + etalcl_sn38 ~ c(0.332,
                              0.269, 0.411)                     # Sathe 2024 Table 2 BLOCK(2): var(KREL) = 0.332, cov(KREL, CLSN38/F) = 0.269, var(CLSN38/F) = 0.411
    etalcl_tab + etalvc_tab ~ c(0.100,
                              0.045, 0.046)                     # Sathe 2024 Table 3 BLOCK(2): var(CL_tAB) = 0.100, cov(CL_tAB, V1_tAB) = 0.045, var(V1_tAB) = 0.046

    # ============================================================
    # Residual error
    # In NONMEM: SG and free SN-38 are log-additive on log(DV)
    # (= proportional in linear space); tAB is linear additive +
    # proportional. A time-after-last-dose effect on RUV variance
    # and a study indicator on RUV are reported but omitted from
    # the simulation model (see vignette Assumptions and deviations).
    # ============================================================
    propSd    <- 0.204429; label("SG proportional residual SD on log scale (Sathe 2024 Table 1; supplement TH(5))")
    propSd_sn38   <- 0.357;    label("Free SN-38 proportional residual SD on log scale (Sathe 2024 Table 2; exp(-1.03))")
    addSd_tab     <- 27.3;     label("tAB additive residual SD (ug/mL; Sathe 2024 Table 3)")
    propSd_tab    <- 0.207;    label("tAB proportional residual SD as fraction (Sathe 2024 Table 3)")
  })
  model({
    # ------------------------------------------------------------
    # SG (sacituzumab govitecan, the ADC)
    # ------------------------------------------------------------
    bwt_cl_factor <- (WT / 70) ^ allo_cl
    bwt_v_factor  <- (WT / 70) ^ allo_v
    alb_cl_factor <- (ALB / 38) ^ e_alb_cl

    cl <- exp(lcl + etalcl) * bwt_cl_factor * alb_cl_factor
    vc <- exp(lvc)            * bwt_v_factor
    q  <- exp(lq)             * bwt_cl_factor
    vp <- exp(lvp)            * bwt_v_factor

    ke  <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------
    # Free SN-38 (sequential to SG via first-order release KREL)
    # ------------------------------------------------------------
    bwt_cl_sn38_factor <- (WT / 70) ^ allo_cl_sn38

    krel    <- exp(lkrel + etalkrel)
    cl_sn38  <- exp(lcl_sn38 + etalcl_sn38) * bwt_cl_sn38_factor
    q_sn38   <- exp(lq_sn38)               * bwt_cl_sn38_factor
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
    bwt_cl_tab_factor <- (WT / 70) ^ allo_cl_tab
    bwt_v_tab_factor  <- (WT / 70) ^ allo_v_tab
    alb_cl_tab_factor <- (ALB / 38) ^ e_alb_cl_tab
    tumor_cl_tab_factor <- 1 + e_tumor_cl_tab * TUMTP_OTH
    sex_vc_tab_factor    <- 1 + e_sex_vc_tab    * (1 - SEXF)
    td_cl_tab_factor    <- 1 - (maxRed_tab / 100) * (1 - exp(-keff_tab * t))

    cl_tab <- exp(lcl_tab + etalcl_tab) * bwt_cl_tab_factor * alb_cl_tab_factor *
             tumor_cl_tab_factor * td_cl_tab_factor
    vc_tab <- exp(lvc_tab + etalvc_tab) * bwt_v_tab_factor  * sex_vc_tab_factor
    q_tab  <- exp(lq_tab)              * bwt_cl_tab_factor
    vp_tab <- exp(lvp_tab)             * bwt_v_tab_factor

    ke_tab  <- cl_tab / vc_tab
    k12_tab <- q_tab  / vc_tab
    k21_tab <- q_tab  / vp_tab

    # ------------------------------------------------------------
    # ODE system (named compartments; user doses central and
    # central_tab simultaneously to simulate one SG infusion).
    # KREL drives free-SN-38 generation from SG without a back-
    # coupling sink on SG, matching the sequential-release
    # parameterization in Sathe 2024 Methods 2.3.
    # ------------------------------------------------------------
    d/dt(central)       <- -ke  * central - k12  * central + k21  * peripheral1
    d/dt(peripheral1)   <-  k12 * central - k21  * peripheral1

    d/dt(central_sn38)     <-  krel   * central - ke_sn38 * central_sn38 -
                              k12_sn38 * central_sn38 + k21_sn38 * peripheral1_sn38
    d/dt(peripheral1_sn38) <-  k12_sn38 * central_sn38 - k21_sn38 * peripheral1_sn38

    d/dt(central_tab)      <- -ke_tab  * central_tab - k12_tab  * central_tab + k21_tab  * peripheral1_tab
    d/dt(peripheral1_tab)  <-  k12_tab * central_tab - k21_tab  * peripheral1_tab

    # ------------------------------------------------------------
    # Observations
    # Concentrations expressed in ug/mL (= mg/L) given dose in mg
    # and volumes in L. Sathe 2024 reports SG and free SN-38 in
    # ng/mL and tAB in ug/mL; multiply Cc and Cc_sn38 by 1000 to
    # match the paper's ng/mL units when comparing.
    # ------------------------------------------------------------
    Cc    <- central   / vc
    Cc_sn38 <- central_sn38 / vc_sn38
    Cc_tab  <- central_tab  / vc_tab

    Cc    ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Cc_tab  ~ add(addSd_tab) + prop(propSd_tab)
  })
}
