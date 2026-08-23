Ibrahim_2025_ibrutinib_venetoclax <- function() {
  description <- paste(
    "QSP / semi-mechanistic PK-PD. Ibrutinib + venetoclax combination model for chronic lymphocytic leukemia (CLL),",
    "the calibrated combination arm of the Ibrahim 2025 framework. Identical ODE structure and ibrutinib parameters to",
    "modellib('Ibrahim_2025_ibrutinib_cll') -- one phosphorylated-Btk turnover compartment plus four CLL cell",
    "subpopulations, with leukocyte count, lymphocyte count, SPD and spleen volume as outputs -- extended by a",
    "venetoclax (BCL-2 inhibitor) killing effect on the CLL cells. Venetoclax acts through two saturable Emax terms",
    "with DIFFERENT potencies by site, reproducing the observation that venetoclax clears circulating CLL cells far",
    "more efficiently than it clears lymph-node disease: in lymphoid tissue it adds a death rate to the three tissue",
    "subpopulations with EC50 = 2.24 ug/mL, while in peripheral blood it multiplies the natural CLL death rate with",
    "EC50 = 0.04 ug/mL and a very large Emax (3465), which is what drives the deep measurable-residual-disease (MRD)",
    "responses the paper simulates. All venetoclax terms vanish at zero venetoclax concentration, so setting",
    "CONC_VEN_MGL = 0 reduces this model exactly to the ibrutinib monotherapy model. Adds `mrd`, the peripheral-blood",
    "MRD percentage (CLL cells as a percentage of total leukocytes), which is the paper's combination-therapy",
    "endpoint. Neither drug has a PK compartment here: ibrutinib enters through the time-varying covariate AUC_IBRU",
    "and venetoclax through the time-varying covariate CONC_VEN_MGL. Sister model files from the same paper:",
    "modellib('Ibrahim_2025_ibrutinib_cll'), modellib('Ibrahim_2025_ibrutinib_bp')."
  )
  reference <- paste(
    "Ibrahim EIK, Friberg LE.",
    "Optimizing ibrutinib posology in chronic lymphocytic leukemia using a semi-mechanistic pharmacometric framework.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2186-2198.",
    "doi:10.1002/psp4.70124.",
    "Open Access under CC BY-NC.",
    "Ibrutinib structure and parameters from the authors' RxODE control stream",
    "(Data S1, PSP-2025-0220-s02.docx, `run8634_eff`) and Table 1;",
    "venetoclax effect equations and parameters from Table S2 (Data S1, PSP-2025-0220-s01.docx).",
    "Venetoclax EC50 values are literature values from Gopalakrishnan S, Wierda W, Chyla B, et al.",
    "Clin Pharmacol Ther. 2021;109(2):424-432. doi:10.1002/cpt.2005.",
    sep = " "
  )
  vignette <- "Ibrahim_2025_ibrutinib"
  paper_specific_compartments <- c(
    "pbtk", "cll_subpop1", "cll_subpop2", "cll_subpop3", "cll_bld"
  )
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; ibrutinib and venetoclax exposure enter as the time-varying covariates AUC_IBRU and CONC_VEN_MGL)",
    concentration = "leukocyte and lymphocyte counts in 10^9 cells/L; SPD in cm^2; spleen volume in cc; MRD in percent (no output is a drug concentration)"
  )

  compartmentData <- list(
    pbtk        = list(analyte = "phosphorylated Bruton tyrosine kinase (pBtk)", units = "relative quantity (1 = 100% of baseline)", specimen = "not applicable", verified = TRUE),
    cll_subpop1 = list(analyte = "CLL cells, stroma-attached proliferating colony 1 (fast detachment)", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_subpop2 = list(analyte = "CLL cells, stroma-attached proliferating colony 2 (slow detachment)", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_subpop3 = list(analyte = "CLL cells, released lymphoid-tissue pool exiting to blood", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_bld     = list(analyte = "CLL cells, resting peripheral-blood pool", units = "10^9 cells", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    AUC_IBRU = list(
      description        = "Daily 0-24 h area under the ibrutinib plasma concentration-time curve, AUC(0-24).",
      units              = "h*ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying; identical to the column used by modellib('Ibrahim_2025_ibrutinib_cll'). Tracks the patient's",
        "current daily ibrutinib dose level (420, 280 or 140 mg/day across the de-escalation schedules simulated by",
        "the paper) and drops to 0 during interruptions. Ibrahim 2025 computed it from the two-compartment ibrutinib",
        "population PK model of Marostica et al. (Cancer Chemother Pharmacol. 2015;75(1):111-121), which is NOT part",
        "of nlmixr2lib. Enters the pBtk production-inhibition Imax function as AUC_IBRU / (28.4 + AUC_IBRU)."
      ),
      source_name        = "auc"
    ),
    CONC_VEN_MGL = list(
      description        = "Venetoclax plasma concentration.",
      units              = "ug/mL (equivalently mg/L)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Venetoclax has no PK compartment in this model: exactly as ibrutinib enters through AUC_IBRU,",
        "venetoclax enters only through this concentration column. Ibrahim 2025 generated venetoclax",
        "concentration-time profiles with the two-compartment population PK model of Jones et al.",
        "(AAPS J. 2016;18(5):1192-1202), which is neither reproduced in the paper nor part of nlmixr2lib, so",
        "downstream users must supply the concentrations from that model or from observed data.",
        "Set CONC_VEN_MGL = 0 for an ibrutinib-monotherapy arm: every venetoclax term vanishes and the model reduces",
        "exactly to modellib('Ibrahim_2025_ibrutinib_cll'). The dosing schedule the paper simulated was a venetoclax",
        "ramp-up starting after ibrutinib cycle 2 or 3 at 20 mg/day with weekly escalation through 50, 100 and",
        "200 mg/day to a maintenance dose of 400 mg/day (Methods 2.4).",
        "The two EC50 values are on this scale (0.04 ug/mL in blood, 2.24 ug/mL in tissue; Table S2)."
      ),
      source_name        = "Ct,venetoclax"
    ),
    LINE_1L = list(
      description        = "First-line-therapy indicator: 1 = treatment-naive (TN) at baseline, 0 = relapsed/refractory (R/R).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (relapsed/refractory)",
      notes              = paste(
        "Time-fixed per subject; identical to the column used by modellib('Ibrahim_2025_ibrutinib_cll'), where the",
        "six places it acts on are documented in full. In this combination model it additionally scales every",
        "venetoclax effect indirectly, because all four venetoclax terms are expressed as multiples of the",
        "peripheral-blood CLL death rate kd,bld, which is itself 43% lower in R/R patients (Table 1).",
        "POLARITY WARNING: the source column `arm` is the INVERSE of this canonical; convert with LINE_1L = 1 - arm."
      ),
      source_name        = "arm"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 246L,
    n_studies      = 2L,
    age_range      = "mean 70 (SD 8.9) years",
    disease_state  = "Chronic lymphocytic leukemia (CLL); 151 (61%) treatment-naive, 95 (39%) relapsed/refractory",
    dose_range     = "ibrutinib 420 mg once daily (reduced to 280 and 140 mg/day in the de-escalation schedules); venetoclax ramped 20 -> 50 -> 100 -> 200 -> 400 mg/day weekly, started after ibrutinib cycle 2 or 3",
    regions        = "United States and international (PCYC-1102 phase Ib/II; PCYC-1115 phase III)",
    notes          = paste(
      "The ibrutinib layer was estimated on the pooled n = 246 PCYC-1102 + PCYC-1115 population (Table S1). The",
      "venetoclax layer was NOT estimated on that population: its EC50 values are literature values from",
      "Gopalakrishnan 2021 and its Emax values were fine-tuned against digitized venetoclax-monotherapy data from the",
      "same publication, then externally validated against digitized data from the CLARITY (Hillmen 2019, n = 53",
      "relapsed/refractory) and CAPTIVATE (Moreno 2023, first-line) ibrutinib + venetoclax studies (Table S2,",
      "Figure S1). Consequently the paper reports NO interindividual variability on any venetoclax parameter --",
      "Discussion: 'IIV in venetoclax response could not be considered, given that the data were derived from",
      "digitized sources.' All between-subject variability in this model therefore comes from the ibrutinib layer."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # IBRUTINIB LAYER. Identical to modellib('Ibrahim_2025_ibrutinib_cll');
    # see that file for the full per-parameter source trace, the day/month
    # unit conversion, and the Table 1 lambda_dec label transposition.
    # Values are Ibrahim 2025 Table 1 point estimates.
    # ------------------------------------------------------------------

    # --- pBtk turnover -------------------------------------------------
    pbtk_base    <- fixed(1)
    label("Baseline phosphorylated-Btk level (relative quantity; 1 = 100% active)")  # Table 1 'pBtkbaseline (%) = 100 FIX'
    lkout_pbtk   <- log(0.524)
    label("Log turnover (elimination) rate constant of pBtk in treatment-naive patients (1/day)")  # Table 1 kout,pBtk,TN = 0.524 (95% CI 0.496-0.554)
    e_rr_kout_pbtk <- log(1.76)
    label("Log fold-change in the pBtk turnover rate constant in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on kout,pBtk = 1.76' (95% CI 0.631-4.90)
    imax_pbtk    <- fixed(1)
    label("Maximum inhibition of pBtk production by ibrutinib (unitless)")  # Table 1 'IMAX = 1 FIX'
    liauc50_pbtk <- log(28.4)
    label("Log daily AUC(0-24) giving 50% of maximum pBtk-production inhibition (h*ng/mL)")  # Table 1 IAUC50,pBtk = 28.4 (95% CI 20.4-39.5)

    # --- CLL cell and biomarker baselines -------------------------------
    lrbase_cllbld_tn  <- log(208)
    label("Log baseline CLL cell number in peripheral blood, treatment-naive patients (10^9 cells)")  # Table 1 CLLbld,baseline,TN = 208 (95% CI 145-299)
    lrbase_cllbld_rr  <- log(53.2)
    label("Log baseline CLL cell number in peripheral blood, relapsed/refractory patients (10^9 cells)")  # Table 1 CLLbld,baseline,R/R = 53.2 (95% CI 32.5-87.3)
    lrbase_clltiss    <- log(3279)
    label("Log baseline total proliferating CLL cell number in lymphoid tissues (10^9 cells)")  # Table 1 CLLtiss,baseline = 3279 (95% CI 2420-4430)
    lrbase_spd        <- log(24.3)
    label("Log baseline SPD in treatment-naive patients (cm^2)")  # Table 1 SPDbaseline,TN = 24.3 (95% CI 17.6-33.6)
    e_rr_rbase_spd    <- log(1.58)
    label("Log fold-change in baseline SPD in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on SPDbaseline = 1.58' (95% CI 1.00-2.49)
    lrbase_spleen     <- log(314)
    label("Log baseline spleen volume (cc)")  # Table 1 SPLEENbaseline = 314 (95% CI 214-463)
    lleuknrm_tn       <- log(37.3)
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, treatment-naive patients (10^9 cells)")  # Table 1 LEUKnrm,TN = 37.3 (95% CI 32.6-42.6)
    lleuknrm_rr       <- log(19.9)
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, relapsed/refractory patients (10^9 cells)")  # Table 1 LEUKnrm,R/R = 19.9 (95% CI 14-28.3)
    lspdnrm           <- log(2.70)
    label("Log normal (non-leukaemic) lymph-node size (cm^2)")  # Table 1 SPDnrm = 2.70 (95% CI 1.98-3.67)
    lspleennrm        <- log(252)
    label("Log normal (non-leukaemic) spleen volume (cc)")  # Table 1 SPLEENnrm = 252 (95% CI 200-317)
    vbld              <- fixed(5)
    label("Blood volume used to convert peripheral-blood cell numbers to counts per litre (L)")  # Table 1 'Vbld (L) = 5 FIX'
    logit_flymphocyte <- log(0.310 / (1 - 0.310))
    label("Logit of the fraction of the normal lymphocyte count relative to the normal leukocyte count (unitless)")  # Table 1 f_lymphocyte = 0.310 (95% CI 0.253-0.374)

    # --- CLL cell kinetics ---------------------------------------------
    lkp        <- log(0.124 / 30)
    label("Log proliferation rate constant of CLL cells (1/day)")  # Table 1 kp = 0.124 month^-1 (95% CI 0.102-0.153); /30 days
    lkh        <- log(0.573)
    label("Log homing rate constant of CLL cells from peripheral blood back to lymphoid tissue (1/day)")  # Table 1 kh = 0.573 (95% CI 0.485-0.678)
    lkd_bld    <- log(0.0153)
    label("Log natural death rate constant of CLL cells in peripheral blood, treatment-naive patients (1/day)")  # Table 1 kd,bld,TN = 0.0153 (95% CI 0.0119-0.0196)
    e_rr_kd_bld <- log(0.571)
    label("Log fold-change in the peripheral-blood CLL death rate constant in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on kd,bld = 0.571' (95% CI 0.440-0.742)
    lkd_tiss   <- log(0.162)
    label("Log ibrutinib-induced death rate constant of CLL cells in lymphoid tissue (1/day)")  # Table 1 kd,tiss = 0.162 (95% CI 0.119-0.22)

    # --- ibrutinib effect slopes ---------------------------------------
    slp1  <- fixed(1)
    label("Slope of the ibrutinib inhibitory effect on proliferation (kp) and homing (kh) (unitless)")  # Table 1 'slp1 = 1 FIX'
    lslp2 <- log(12.0)
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop1 from stroma (unitless)")  # Table 1 slp2 = 12.0 (95% CI 9.79-14.7)
    lslp3 <- log(0.140)
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop2 from stroma (unitless)")  # Table 1 slp3 = 0.140 (95% CI 0.106-0.183)

    # --- subpopulation fractions (logit scale) --------------------------
    logit_f1 <- log(0.483 / (1 - 0.483))
    label("Logit of f1, the fraction of stroma-attached CLL cells that belong to the slow-detaching second subpopulation (unitless)")  # Table 1 f1 = 0.483 (95% CI 0.444-0.522)
    logit_f2 <- log(0.246 / (1 - 0.246))
    label("Logit of f2, the fraction of CLLtiss,baseline that belongs to the released third subpopulation (unitless)")  # Table 1 f2 = 0.246 (95% CI 0.202-0.297)

    # --- resistance -----------------------------------------------------
    # Assigned to relapsed/refractory patients; Table 1's 'TN' / 'R/R' row
    # labels for lambda_dec are transposed relative to the authors' code
    # (`resist = exp(-lmbd*iarm*t)`, iarm = 1 in R/R) and to three separate
    # prose statements. Full argument in the ini() of
    # modellib('Ibrahim_2025_ibrutinib_cll') and in the vignette Errata.
    llambda_dec <- log(0.0230 / 30)
    label("Log exponential decay constant of the ibrutinib proliferation/apoptosis effect in relapsed/refractory patients (1/day)")  # Table 1 lambda_dec = 0.0230 month^-1 (95% CI 0.0109-0.0484), labelled 'TN' but belonging to R/R; /30 days

    # --- Box-Cox shape and shared-eta scaling ---------------------------
    boxcox_cllbld <- -0.199
    label("Box-Cox shape parameter for the baseline peripheral-blood CLL count random effects (unitless)")  # Table 1 'Box-Cox for CLLbld,baseline random effects' = -0.199 (95% CI -0.329 to -0.0684)
    sc_kdb <- fixed(0.824)
    label("Scaling of the shared CLLtiss,baseline random effect onto kd,bld, reproducing their 100% correlation (unitless)")  # derived from Table 1 CV% for CLLtiss,baseline (152%) and kd,bld,TN (112%); Table 1 footnote b

    # ------------------------------------------------------------------
    # VENETOCLAX LAYER (Table S2, Data S1 PSP-2025-0220-s01.docx).
    #
    # Table S2 footnote a gives the two effect forms verbatim:
    #   CLL subpopulations in tissues:
    #     kd,bld * Emax,i * Ct,venetoclax / (EC50,tissue + Ct,venetoclax)
    #     where i denotes CLL subpopulation 1&2 or 3
    #   CLL subpopulation 4 in peripheral blood:
    #     kd,bld * (1 + Emax,4 * Ct,venetoclax / (EC50,blood + Ct,venetoclax))
    # Every venetoclax parameter is FIXED, not estimated: the two EC50s are
    # literature values from Gopalakrishnan 2021 and the three Emax values
    # were fine-tuned against digitized data from the same paper (Table S2
    # footnote b). None carries interindividual variability -- the paper's
    # Discussion states that "IIV in venetoclax response could not be
    # considered, given that the data were derived from digitized sources."
    # ------------------------------------------------------------------
    ec50_ven_bld  <- fixed(0.04)
    label("Venetoclax concentration giving 50% of the maximum killing effect on peripheral-blood CLL cells (ug/mL)")  # Table S2 'EC50, blood (ug/mL) = 0.04'; literature value from Gopalakrishnan 2021 (doi:10.1002/cpt.2005)
    ec50_ven_tiss <- fixed(2.24)
    label("Venetoclax concentration giving 50% of the maximum killing effect on lymphoid-tissue CLL cells (ug/mL)")  # Table S2 'EC50, tissues (ug/mL) = 2.24'; literature value from Gopalakrishnan 2021 (doi:10.1002/cpt.2005)
    emax_ven_tiss12 <- fixed(0.63)
    label("Maximum venetoclax killing effect on CLL subpopulations 1 and 2 in lymphoid tissue, as a multiple of kd,bld (unitless)")  # Table S2 'Emax,1&2 = 0.63'; fine-tuned on digitized Gopalakrishnan 2021 data (footnote b)
    emax_ven_tiss3  <- fixed(3.15)
    label("Maximum venetoclax killing effect on CLL subpopulation 3 in lymphoid tissue, as a multiple of kd,bld (unitless)")  # Table S2 'Emax,3 = 3.15'; fine-tuned on digitized Gopalakrishnan 2021 data (footnote b)
    emax_ven_bld    <- fixed(3465)
    label("Maximum fractional venetoclax enhancement of the peripheral-blood CLL death rate constant (unitless)")  # Table S2 'Emax,4 = 3465'; fine-tuned on digitized Gopalakrishnan 2021 data (footnote b)

    # --- IIV (ibrutinib layer only; see the sister file for the CV% arithmetic)
    etalkout_pbtk       ~ 3.6997     # Table 1 kout,pBtk CV% = 628 -> ln(1 + 6.28^2)
    etalrbase_cllbld_tn ~ 1.7689     # Table 1 CLLbld,baseline,TN CV% = 133 -> footnote a: omega^2 = 1.33^2
    etalrbase_cllbld_rr ~ 4.1616     # Table 1 CLLbld,baseline,R/R CV% = 204 -> footnote a: omega^2 = 2.04^2
    etalogit_flymphocyte ~ 0.37085   # Table 1 f_lymphocyte CV% = 67 -> ln(1 + 0.67^2)
    etalrbase_clltiss   ~ 1.19726    # Table 1 CLLtiss,baseline CV% = 152 -> ln(1 + 1.52^2); also drives kd,bld via sc_kdb
    etalrbase_spleen    ~ 0.78300    # Table 1 SPLEENbaseline CV% = 109 -> ln(1 + 1.09^2)
    etalleuknrm_tn      ~ 0.115556   # Table 1 LEUKnrm,TN CV% = 35 -> ln(1 + 0.35^2)
    etalleuknrm_rr      ~ 0.417684   # Table 1 LEUKnrm,R/R CV% = 72 -> ln(1 + 0.72^2)
    etalspleennrm       ~ 0.316354   # Table 1 SPLEENnrm CV% = 61 -> ln(1 + 0.61^2)
    etalogit_f1         ~ 2.74188    # Table 1 f1 CV% = 381 -> ln(1 + 3.81^2)
    etalogit_f2         ~ 1.09463    # Table 1 f2 CV% = 141 -> ln(1 + 1.41^2)
    etalkp              ~ 0.75308    # Table 1 kp CV% = 106 -> ln(1 + 1.06^2)
    etalkh              ~ 3.72440    # Table 1 kh CV% = 636 -> ln(1 + 6.36^2)
    etaliauc50_pbtk     ~ 2.99621    # Table 1 IAUC50,pBtk CV% = 436 -> ln(1 + 4.36^2)
    etalkd_tiss         ~ 1.64131    # Table 1 kd,tiss CV% = 204 -> ln(1 + 2.04^2)
    etallambda_dec      ~ 1.20639    # Table 1 lambda_dec CV% = 153 -> ln(1 + 1.53^2)
    etalspdnrm_rad      ~ 0.63313    # Table 1 footnote d inter-radiologist variability in SPDnrm = 94% -> ln(1 + 0.94^2)

    etalrbase_spd + etalspdnrm ~ c(0.98954,
                                   0.398635, 0.380093)  # Table 1 SPDbaseline,TN CV% = 130, SPDnrm CV% = 68, footnote b correlation 65%

    # --- residual error --------------------------------------------------
    expSd_leukocyte    <- 0.20
    label("Log-scale (log-normal) residual error SD for leukocyte count (unitless on the log scale)")  # Table 1 'RUV leukocyte count (%) = 20'
    expSd_lymphocyte   <- 0.21
    label("Log-scale (log-normal) residual error SD for lymphocyte count (unitless on the log scale)")  # Table 1 'RUV lymphocyte count (%) = 21'
    expSd_tumorSpd     <- 0.27
    label("Log-scale (log-normal) residual error SD for total SPD (unitless on the log scale)")  # Table 1 'RUV SPD (%) = 27'
    expSd_spleenVolume <- 0.13
    label("Log-scale (log-normal) residual error SD for spleen volume (unitless on the log scale)")  # Table 1 'RUV spleen volume (%) = 13'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters (ibrutinib layer).
    # ------------------------------------------------------------------
    etabox_cllbld_tn <- ((exp(etalrbase_cllbld_tn)^boxcox_cllbld) - 1) / boxcox_cllbld
    etabox_cllbld_rr <- ((exp(etalrbase_cllbld_rr)^boxcox_cllbld) - 1) / boxcox_cllbld

    cllbld0_tn <- exp(lrbase_cllbld_tn + etabox_cllbld_tn)
    cllbld0_rr <- exp(lrbase_cllbld_rr + etabox_cllbld_rr)

    clltiss0   <- exp(lrbase_clltiss + etalrbase_clltiss)
    spleenbase <- exp(lrbase_spleen + etalrbase_spleen)
    spleennrm  <- exp(lspleennrm + etalspleennrm)
    leuknrm_tn <- exp(lleuknrm_tn + etalleuknrm_tn)
    leuknrm_rr <- exp(lleuknrm_rr + etalleuknrm_rr)
    # Two random effects applied as a product of exponentials: rxode2
    # cannot mu-reference `exp(theta + eta1 + eta2)`. Algebraically identical.
    spdnrm     <- exp(lspdnrm + etalspdnrm) * exp(etalspdnrm_rad)

    kp      <- exp(lkp + etalkp)
    kh      <- exp(lkh + etalkh)
    kd_tiss <- exp(lkd_tiss + etalkd_tiss)
    slp2    <- exp(lslp2)
    slp3    <- exp(lslp3)

    # Fixed effects parenthesised together: rxode2's mu-referencing rejects
    # two bare population parameters in one expression. Arithmetic unchanged.
    kout_pbtk_tn <- exp(lkout_pbtk + etalkout_pbtk)
    kout_pbtk_rr <- exp((lkout_pbtk + e_rr_kout_pbtk) + etalkout_pbtk)
    iauc50_pbtk  <- exp(liauc50_pbtk + etaliauc50_pbtk)
    lambda_dec   <- exp(llambda_dec + etallambda_dec)

    kd_bld_tn <- exp(lkd_bld + etalrbase_clltiss * sc_kdb)
    kd_bld_rr <- exp((lkd_bld + e_rr_kd_bld) + etalrbase_clltiss * sc_kdb)

    odds_f1 <- exp(logit_f1 + etalogit_f1)
    f1      <- odds_f1 / (1 + odds_f1)
    odds_f2 <- exp(logit_f2 + etalogit_f2)
    f2      <- odds_f2 / (1 + odds_f2)

    odds_flymphocyte <- exp(logit_flymphocyte + etalogit_flymphocyte)
    f_lymphocyte     <- odds_flymphocyte / (1 + odds_flymphocyte)

    # ------------------------------------------------------------------
    # 2. Covariate model. LINE_1L = 1 - arm, so LINE_1L = 1 is treatment-naive.
    # ------------------------------------------------------------------
    cllbld0   <- cllbld0_tn * LINE_1L + cllbld0_rr * (1 - LINE_1L)
    leuknrm   <- leuknrm_tn * LINE_1L + leuknrm_rr * (1 - LINE_1L)
    kout_pbtk <- kout_pbtk_tn * LINE_1L + kout_pbtk_rr * (1 - LINE_1L)
    kd_bld    <- kd_bld_tn * LINE_1L + kd_bld_rr * (1 - LINE_1L)
    spdbase   <- exp((lrbase_spd + e_rr_rbase_spd * (1 - LINE_1L)) + etalrbase_spd)

    # ------------------------------------------------------------------
    # 3. Ibrutinib effects.
    # ------------------------------------------------------------------
    kin_pbtk <- kout_pbtk * pbtk_base
    effauc   <- imax_pbtk * AUC_IBRU / (iauc50_pbtk + AUC_IBRU)
    pbtk_eff <- pbtk_base - pbtk
    resist   <- exp(-lambda_dec * (1 - LINE_1L) * t)
    kdtch    <- kp

    eff_prol  <- slp1 * resist * pbtk_eff
    eff_death <- kd_tiss * resist * pbtk_eff
    eff_home  <- slp1 * pbtk_eff
    eff_dtch1 <- slp2 * pbtk_eff
    eff_dtch2 <- slp3 * pbtk_eff

    cll_subpop3_base <- f2 * clltiss0
    kdist            <- (kh + kd_bld) * (cllbld0 / cll_subpop3_base)

    # ------------------------------------------------------------------
    # 4. Venetoclax effects (Table S2 footnote a). Both site-specific
    #    terms are scaled by kd,bld, so they inherit its TN/RR covariate.
    #    Each vanishes at CONC_VEN_MGL = 0, and the blood term collapses to
    #    kd_bld itself, so the whole model reduces to the ibrutinib
    #    monotherapy model when venetoclax is absent.
    # ------------------------------------------------------------------
    ven_sat_tiss <- CONC_VEN_MGL / (ec50_ven_tiss + CONC_VEN_MGL)
    ven_sat_bld  <- CONC_VEN_MGL / (ec50_ven_bld + CONC_VEN_MGL)

    # additional death rate constant on the tissue subpopulations
    ven_kill_tiss12 <- kd_bld * emax_ven_tiss12 * ven_sat_tiss
    ven_kill_tiss3  <- kd_bld * emax_ven_tiss3 * ven_sat_tiss

    # multiplicative enhancement of the peripheral-blood death rate constant
    kd_bld_ven <- kd_bld * (1 + emax_ven_bld * ven_sat_bld)

    # ------------------------------------------------------------------
    # 5. ODE system: the Ibrahim 2025 monotherapy structure with the
    #    venetoclax killing terms added.
    # ------------------------------------------------------------------
    d/dt(pbtk) <- kin_pbtk * (1 - effauc) - kout_pbtk * pbtk
    pbtk(0)    <- pbtk_base

    d/dt(cll_subpop1) <- kp * (1 - eff_prol) * cll_subpop1 -
      kdtch * (1 + eff_dtch1) * cll_subpop1 -
      ven_kill_tiss12 * cll_subpop1
    cll_subpop1(0)    <- clltiss0 * (1 - f2) * (1 - f1)

    d/dt(cll_subpop2) <- kp * (1 - eff_prol) * cll_subpop2 -
      kdtch * (1 + eff_dtch2) * cll_subpop2 -
      ven_kill_tiss12 * cll_subpop2
    cll_subpop2(0)    <- clltiss0 * (1 - f2) * f1

    d/dt(cll_subpop3) <- kp * (1 - eff_prol) * cll_subpop3 +
      kdtch * (1 + eff_dtch1) * cll_subpop1 +
      kdtch * (1 + eff_dtch2) * cll_subpop2 +
      kh * (1 - eff_home) * cll_bld -
      (kdist + eff_death + ven_kill_tiss3) * cll_subpop3
    cll_subpop3(0)    <- clltiss0 * f2

    d/dt(cll_bld) <- kdist * cll_subpop3 -
      kh * (1 - eff_home) * cll_bld -
      kd_bld_ven * cll_bld
    cll_bld(0)    <- cllbld0

    # ------------------------------------------------------------------
    # 6. Observations (Ibrahim 2025 equations 2-5, plus MRD).
    # ------------------------------------------------------------------
    cll_tiss_total <- cll_subpop1 + cll_subpop2 + cll_subpop3

    tumorSpd     <- cll_tiss_total * (spdbase / clltiss0) + spdnrm
    spleenVolume <- cll_tiss_total * (spleenbase / clltiss0) + spleennrm
    leukocyte    <- (cll_bld + leuknrm) / vbld
    lymphocyte   <- (cll_bld + leuknrm * f_lymphocyte) / vbld

    totalCll     <- cll_bld + cll_tiss_total
    btkOccupancy <- (pbtk_base - pbtk) * 100

    # Peripheral-blood measurable residual disease, the paper's
    # combination-therapy endpoint: "MRD was assessed in peripheral blood as
    # the percentage of CLL cells relative to the total leukocyte count"
    # (Methods 2.4). The paper's MRD-negativity threshold is 0.001%.
    mrd <- 100 * cll_bld / (cll_bld + leuknrm)

    leukocyte    ~ lnorm(expSd_leukocyte)
    lymphocyte   ~ lnorm(expSd_lymphocyte)
    tumorSpd     ~ lnorm(expSd_tumorSpd)
    spleenVolume ~ lnorm(expSd_spleenVolume)
  })
}
