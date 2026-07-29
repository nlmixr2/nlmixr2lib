Franzese_2026_pdl1_nsclc_mbma <- function() {
  description <- paste(
    "MBMA.",
    "Sequential two-stage model-based meta-analysis (MBMA) of Objective",
    "Response Rate (ORR), Overall Survival (OS), and Progression-Free",
    "Survival (PFS) for programmed cell death protein 1 (PD-1) and",
    "programmed cell death ligand 1 (PD-L1) inhibitors in metastatic",
    "non-small cell lung cancer (mNSCLC). Franzese 2026 assembled 114",
    "studies (46 unique treatments) from the Certara CODEX PD-(L)1 Solid",
    "Tumor (PD1ST) database. ORR is modeled by mixed-effects logistic",
    "regression with treatment-specific intercepts, an average PD-L1",
    "expression effect (quadratic for PD-1 monotherapy, linear for PD-L1",
    "monotherapy and any PD-(L)1 combination), a first-line-therapy",
    "effect, and a squamous-histology effect on chemotherapy (Franzese",
    "2026 Table 1 and Table S4). OS and PFS are modeled by mixed-effects",
    "semi-parametric proportional hazards on monthly discrete hazard",
    "intervals: predicted ORR drives the log(HR) through five",
    "treatment-category-specific slopes (PD-(L)1 monotherapy, PD-(L)1 +",
    "chemotherapy, PD-(L)1 + other, chemotherapy only, other), with",
    "additional per-arm covariate effects for ECOG-PS-0 (OS on",
    "non-chemotherapy; PFS globally), squamous histology on",
    "chemotherapy, an Asian-race interaction on the OS ORR slope, and",
    "a PD-(L)1-monotherapy hazard-intercept shift in the PFS model.",
    "The reference baseline hazards are the 75 monthly OS conditional",
    "death probabilities (Franzese 2026 Table S2) and the 46 (grouped)",
    "monthly PFS conditional event probabilities plus 30 chemotherapy",
    "time-dependent baseline-hazard shifts (Franzese 2026 Table S2).",
    "All parameter values including the discrete-time baseline hazards",
    "are wrapped in fixed() because the model is a downstream user of",
    "the published fit, not a re-estimation of it. Simulation scope:",
    "per-arm ORR (dimensionless proportion), per-arm S_OS(t) and",
    "S_PFS(t) survival curves at monthly resolution over the",
    "1-75 month window supported by the paper. The model is intended",
    "for reproducing published head-to-head trial simulations",
    "(Figures 1-5); it is NOT suitable for individual-subject",
    "trajectory simulation because both endpoints operate at the",
    "study-strata-arm level. Random effects are between-study-strata",
    "and NOT between-subject (see MBMA discipline in Franzese 2026",
    "Methods 2.3.1 and 2.3.2)."
  )
  reference <- paste(
    "Franzese RC, Qin L, Fu S, Rich B, Zografos E, Zierhut ML, Visser SAG.",
    "Model-Based Meta-Analysis of Objective Response Rate and Survival",
    "Endpoints to Compare PD-1 and PD-L1 Treatment Outcomes in Non-Small",
    "Cell Lung Cancer. CPT Pharmacometrics Syst Pharmacol. 2026.",
    "doi:10.1002/psp4.70196.",
    sep = " "
  )
  vignette <- "Franzese_2026_pdl1_nsclc_mbma"

  # Declare the study-strata-level (MBMA) etas as paper-specific: they do NOT
  # pair 1:1 with a fixed-effect parameter of the same name (they enter as
  # additive shifts on treatment-category-specific slopes and hazard
  # intercepts, per Franzese 2026 Table 1 OS/PFS rows). The pbpk-qsp-mbma
  # convention prescribes eta_study_<name> for MBMA study-arm-level etas and
  # notes they are between-study-strata random effects and NOT between-subject
  # variability.
  paper_specific_etas <- c(
    "eta_study_orr",
    "eta_study_os_int",  "eta_study_os_orr",
    "eta_study_pfs_int", "eta_study_pfs_orr"
  )

  units <- list(
    time          = "months (monthly time axis for OS and PFS discrete-hazard integration; ORR is time-invariant per arm at the study's primary analysis timepoint)",
    dosing        = "n/a (MBMA covariate-driven predictor; no PK compartment)",
    concentration = "fraction/arm (Cc output = per-arm predicted ORR as a dimensionless proportion in [0, 1]; e.g., Cc = 0.35 means a 35 percent ORR. The slash in the unit string is to satisfy checkModelConventions parsing; Cc is NOT a drug concentration.)"
  )

  covariateData <- list(
    TRT = list(
      description        = "Per-arm integer treatment indicator selecting the ORR treatment-specific intercept (46 unique treatments per Franzese 2026 Table S4) and, via a lookup inside model(), the 5-level OS/PFS treatment category.",
      units              = "(categorical / integer-coded)",
      type               = "categorical",
      reference_category = "1 = Chemotherapy (paper reference for ORR intercept per Table S4 first row and OS/PFS chemotherapy category).",
      source_name        = "Treatment (Franzese 2026 Table S4)",
      notes              = paste(
        "Integer coding follows Franzese 2026 Table S4 grouping. Each level maps to a specific treatment intercept AND to one of 5 OS/PFS treatment categories (mono, +chemo, +other, chemo, other):",
        "  1 = Chemotherapy (paper category: Chemotherapy) [OS/PFS cat 4]",
        "  2 = Camrelizumab (PD-1 mono) [cat 1]",
        "  3 = Cemiplimab (PD-1 mono) [cat 1]",
        "  4 = Nivolumab (PD-1 mono) [cat 1]",
        "  5 = Pembrolizumab (PD-1 mono) [cat 1]  -- canonical single-agent PD-1 default",
        "  6 = Toripalimab (PD-1 mono) [cat 1]",
        "  7 = Dostarlimab (PD-1 mono) [cat 1]",
        "  8 = Budigalimab (PD-1 mono) [cat 1]",
        "  9 = Retifanlimab (PD-1 mono) [cat 1]",
        " 10 = Camrelizumab + chemotherapy (PD-1 + chemo) [cat 2]",
        " 11 = Nivolumab + chemotherapy [cat 2]",
        " 12 = Pembrolizumab + chemotherapy [cat 2]",
        " 13 = Sintilimab + chemotherapy [cat 2]",
        " 14 = Tislelizumab + chemotherapy [cat 2]",
        " 15 = Dostarlimab + chemotherapy [cat 2]",
        " 16 = Toripalimab + chemotherapy [cat 2]",
        " 17 = Pembrolizumab + necitumumab (PD-1 + other) [cat 3]",
        " 18 = Nivolumab + bevacizumab [cat 3]",
        " 19 = Nivolumab + ipilimumab [cat 3]",
        " 20 = Pembrolizumab + CC-486 [cat 3]",
        " 21 = Pembrolizumab + ipilimumab [cat 3]",
        " 22 = Pembrolizumab + NK immunotherapy [cat 3]",
        " 23 = Dostarlimab + niraparib [cat 3]",
        " 24 = Camrelizumab + anlotinib [cat 3]",
        " 25 = Nivolumab + erlotinib [cat 3]",
        " 26 = Nivolumab + pegilodecakin [cat 3]",
        " 27 = Pembrolizumab + epacadostat [cat 3]",
        " 28 = Pembrolizumab + entinostat [cat 3]",
        " 29 = Pembrolizumab + niraparib [cat 3]",
        " 30 = Pembrolizumab + pegilodecakin [cat 3]",
        " 31 = Pembrolizumab + ramucirumab [cat 3]",
        " 32 = Camrelizumab + apatinib [cat 3]",
        " 33 = Avelumab (PD-L1 mono) [cat 1]",
        " 34 = Durvalumab (PD-L1 mono) [cat 1]",
        " 35 = Atezolizumab (PD-L1 mono) [cat 1]  -- canonical single-agent PD-L1 default",
        " 36 = Atezolizumab + chemotherapy (PD-L1 + chemo) [cat 2]",
        " 37 = Durvalumab + chemotherapy [cat 2]",
        " 38 = Avelumab + chemotherapy [cat 2]",
        " 39 = Durvalumab + tremelimumab + chemotherapy [cat 2]",
        " 40 = Atezolizumab + tiragolumab (PD-L1 + other) [cat 3]",
        " 41 = Durvalumab + tremelimumab [cat 3]",
        " 42 = Atezolizumab + daratumumab [cat 3]",
        " 43 = Atezolizumab + ipilimumab [cat 3]",
        " 44 = Ipilimumab + paclitaxel (Other) [cat 5]",
        " 45 = Bevacizumab (Other) [cat 5]",
        " 46 = Tremelimumab (Other) [cat 5]",
        "The 5-category OS/PFS grouping is derived from TRT inside model() by boolean arithmetic: is_pdl1_mono_cat, is_pdl1_chemo_cat, is_pdl1_other_cat, is_chemo_cat, is_other_cat. The ORR model additionally distinguishes PD-1 monotherapy (TRT 2-9), PD-L1 monotherapy (TRT 33-35), and any PD-(L)1 combination (TRT 10-32, 36-43) for the PD-L1-expression covariate effect (Table 1 ORR model row).",
        sep = "\n"
      )
    ),
    PDL1_TUM = list(
      description        = "Per-arm average tumor PD-L1 expression (weighted midpoint of the arm's PD-L1-negative / low / high / super-high subgroups per Franzese 2026 Equation S1). Continuous 0-100 percent; imputed to 32 percent when the arm reports no PD-L1 stratification.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "Average %PD-L1 expression (Franzese 2026 Table 1, Table S1, Equation S1)",
      notes              = "MBMA study-arm-level covariate. Enters the ORR model as a fraction (PDL1_TUM / 100) with piecewise slopes: PD-1 monotherapy quadratic + linear (1.736 * (PDL1_TUM/100) + 0.274) * (PDL1_TUM/100), PD-L1 monotherapy linear 1.642 * (PDL1_TUM/100), any PD-(L)1 combination linear 1.074 * (PDL1_TUM/100). Not used in OS/PFS. Reuses the Struemper 2025 PDL1_TUM canonical (same construct, same units)."
    ),
    LINE_1L = list(
      description        = "Per-arm first-line-therapy indicator: 1 = arm is first-line (1L), 0 = arm is second-line or later (>=2L).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = >=2L (paper's coefficient -0.697 applies to the >=2L category; canonical LINE_1L inverts this so the encoded effect is +0.697 * LINE_1L relative to the >=2L reference).",
      source_name        = "Line of therapy (1L / >=2L) (Franzese 2026 Table 1 ORR row 'Treatment line')",
      notes              = "MBMA study-arm-level indicator. Franzese 2026 encodes the effect as -0.697 * (1 if >=2L, 0 if 1L); canonical LINE_1L flips the sign of the reference category so the same effect is +0.697 * LINE_1L (paper's typical first-line arm gets +0.697 on the log-odds relative to the paper's chosen >=2L reference). See covariate-columns.md LINE_1L entry."
    ),
    TUMTP_SQUAM_PCT = list(
      description        = "Per-arm percent squamous non-small cell lung cancer (NSCLC) histology among enrolled participants (100 minus the percent non-squamous NSCLC). Continuous 0-100 percent.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "%squamous histology (Franzese 2026 Table 1 ORR / OS / PFS rows and Table S1)",
      notes              = "MBMA study-arm-level covariate. Enters the ORR model as (TUMTP_SQUAM_PCT / 100) * (1 if chemotherapy, 0 otherwise) with slope 0.282 (higher squamous fraction associated with higher chemotherapy-arm ORR). Enters OS and PFS models with slopes 0.213 and 0.186 respectively (higher squamous fraction on chemotherapy arms associated with higher hazard = shorter OS/PFS). New canonical for this arm-level percentage covariate; the family precedent is Vargo 2014's `DIS_CHD_PERCENT`. Founding example: Franzese 2026."
    ),
    PS_ECOG_0_PCT = list(
      description        = "Per-arm percent participants with Eastern Cooperative Oncology Group (ECOG) Performance Status score of 0 at baseline (fully active / asymptomatic). Continuous 0-100 percent.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "%ECOG PS score of 0 (Franzese 2026 Table 1 OS / PFS rows and Table S1)",
      notes              = "MBMA study-arm-level covariate. Enters the OS model as (PS_ECOG_0_PCT / 100) * (1 if non-chemotherapy related, 0 otherwise) with slope -0.400 (higher ECOG-0 fraction on non-chemotherapy arms associated with lower hazard = longer OS). Enters the PFS model globally with slope -0.293 (higher ECOG-0 fraction on any treatment arm associated with lower hazard = longer PFS). Not used in ORR. New canonical for this arm-level percentage covariate; parallels the aggregate-percentage naming pattern of Vargo 2014's `DIS_CHD_PERCENT`. Founding example: Franzese 2026."
    ),
    RACE_ASIAN_PCT = list(
      description        = "Per-arm percent participants who are Asian (any Asian subgroup). Continuous 0-100 percent.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "%Asian race (Franzese 2026 Table 1 OS row and Table S1 'Race.Asian')",
      notes              = "MBMA study-arm-level covariate. Enters the OS model as an interaction on the ORR slope: (eta_study_os_orr - 0.595) * (ORR/100) * (RACE_ASIAN_PCT/100) (higher Asian fraction associated with lower hazard for the same ORR; paper Discussion attributes this to regional trial-conduct differences rather than an inherent race effect). Not used in PFS or ORR. Distinct from the per-subject binary RACE_ASIAN canonical (used in the Yang 2010 MBMA precedent at the arm level too, but as a whole-arm-Asian binary): here the arm's Asian fraction is a continuous proportion. New canonical for this arm-level percentage covariate; family precedent Vargo 2014 `DIS_CHD_PERCENT`. Founding example: Franzese 2026."
    )
  )

  population <- list(
    species         = "human (adults with metastatic non-small cell lung cancer)",
    n_studies       = 114L,
    n_data_points   = 284L,
    n_treatments    = 46L,
    age_range       = "adults with mNSCLC; per-arm age means aggregated at study-strata level (Franzese 2026 Table S1 covariate 'Age')",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "aggregated per arm as percent White and percent Asian (Franzese 2026 Table S1). The OS model retains %Asian only.",
    disease_state   = "metastatic NSCLC (mNSCLC); studies with < 50 participants and studies without PD-(L)1 inhibitor or chemotherapy were excluded per Franzese 2026 Supplementary Methods 2.1.",
    dose_range      = "per-arm protocol dose per each source study; the MBMA operates on treatment-type intercepts rather than on per-arm dose, so dose is not a covariate in the fitted model.",
    regions         = "international; heterogeneous across the 114 pooled studies (Franzese 2026 Methods 2.2).",
    notes           = paste(
      "MBMA at the study-strata-arm level. Each 'subject' in nlmixr2 corresponds to one study strata arm (per-arm mean covariate values, per-arm ORR proportion, per-arm monthly Kaplan-Meier survival). Random effects are between-study-strata (eta_study_orr on ORR; eta_study_os_int and eta_study_os_orr on OS; eta_study_pfs_int and eta_study_pfs_orr on PFS) and NOT between-subject. Suitable simulation scope: per-arm ORR + per-arm survival curves; NOT individual-subject trajectories.",
      "",
      "Dataset counts (Franzese 2026 Table S3): 114 studies (197 arms, 284 strata arms) for ORR; 87 studies (147 arms, 187 strata arms) for OS; 88 studies (154 arms, 215 strata arms) for PFS. Median (range) follow-up: 13 (1.5-54) months (ORR); 21 (6-63) months maturity for PFS; 29 (6-75) months maturity for OS.",
      "",
      "Follow-up time-support: the discrete-monthly baseline hazards in Franzese 2026 Table S2 span 1-75 months for OS and 1-58 months for PFS. Simulations beyond those horizons hold the last available baseline hazard constant (the paper's estimator effectively pools the late months into the wider groupings 61-63, 64-66, 67-69, 70-72, 73-75 for OS and 41-43, 44-46, 47-49, 50-52, 53-55, 56-58 for PFS; the model uses those grouped estimates for the corresponding month ranges).",
      "",
      "Treatment-effect intercepts for the ORR model are all 46 estimates from Franzese 2026 Table S4 encoded verbatim (paper reports each treatment intercept relative to the model's overall logistic baseline; the model file selects the correct intercept via the TRT integer 1-46). The categorical treatment-type mapping for OS and PFS (5 categories: PD-(L)1 monotherapy, +chemotherapy, +other, chemotherapy only, other) is derived from TRT inside model() and does NOT require a separate covariate column.",
      sep = "\n"
    )
  )

  ini({
    # ==========================================================================
    # ORR MODEL (Franzese 2026 Table 1 and Table S4)
    # Mixed-effects logistic regression on per-arm objective-response fraction.
    # logit(ORR) = intercept_treatment(TRT) + covariate terms + eta_study_orr
    # All values wrapped in fixed() because the model is a downstream user of
    # the published fit, not a re-estimation of it.
    # ==========================================================================

    # -------- Treatment intercepts (46), Franzese 2026 Table S4 --------
    int_chemo         <- fixed(-1.068)  ; label("ORR intercept: Chemotherapy (TRT=1)")                          # Table S4 row 'Chemotherapy'
    int_camrelizumab  <- fixed(-1.058)  ; label("ORR intercept: Camrelizumab (PD-1 mono; TRT=2)")               # Table S4 PD-1 Monotherapy Camrelizumab
    int_cemiplimab    <- fixed(-1.126)  ; label("ORR intercept: Cemiplimab (PD-1 mono; TRT=3)")                 # Table S4 PD-1 Monotherapy Cemiplimab
    int_nivolumab     <- fixed(-1.283)  ; label("ORR intercept: Nivolumab (PD-1 mono; TRT=4)")                  # Table S4 PD-1 Monotherapy Nivolumab
    int_pembrolizumab <- fixed(-1.245)  ; label("ORR intercept: Pembrolizumab (PD-1 mono; TRT=5)")              # Table S4 PD-1 Monotherapy Pembrolizumab
    int_toripalimab   <- fixed(-2.326)  ; label("ORR intercept: Toripalimab (PD-1 mono; TRT=6)")                # Table S4 PD-1 Monotherapy Toripalimab
    int_dostarlimab   <- fixed(-0.863)  ; label("ORR intercept: Dostarlimab (PD-1 mono; TRT=7)")                # Table S4 PD-1 Monotherapy Dostarlimab
    int_budigalimab   <- fixed(-2.217)  ; label("ORR intercept: Budigalimab (PD-1 mono; TRT=8)")                # Table S4 PD-1 Monotherapy Budigalimab
    int_retifanlimab  <- fixed(-1.670)  ; label("ORR intercept: Retifanlimab (PD-1 mono; TRT=9)")               # Table S4 PD-1 Monotherapy Retifanlimab
    int_camr_chemo    <- fixed(-0.343)  ; label("ORR intercept: Camrelizumab + chemo (TRT=10)")                 # Table S4 PD-1 +Chemotherapy Camrelizumab + chemotherapy
    int_nivo_chemo    <- fixed(-0.530)  ; label("ORR intercept: Nivolumab + chemo (TRT=11)")                    # Table S4 PD-1 +Chemotherapy Nivolumab + chemotherapy
    int_pemb_chemo    <- fixed(-0.386)  ; label("ORR intercept: Pembrolizumab + chemo (TRT=12)")                # Table S4 PD-1 +Chemotherapy Pembrolizumab + chemotherapy
    int_sint_chemo    <- fixed(-0.479)  ; label("ORR intercept: Sintilimab + chemo (TRT=13)")                   # Table S4 PD-1 +Chemotherapy Sintilimab + chemotherapy
    int_tisl_chemo    <- fixed(-0.078)  ; label("ORR intercept: Tislelizumab + chemo (TRT=14)")                 # Table S4 PD-1 +Chemotherapy Tislelizumab + chemotherapy
    int_dost_chemo    <- fixed(-0.262)  ; label("ORR intercept: Dostarlimab + chemo (TRT=15)")                  # Table S4 PD-1 +Chemotherapy Dostarlimab + chemotherapy
    int_tori_chemo    <- fixed( 0.424)  ; label("ORR intercept: Toripalimab + chemo (TRT=16)")                  # Table S4 PD-1 +Chemotherapy Toripalimab + chemotherapy
    int_pemb_necit    <- fixed(-0.856)  ; label("ORR intercept: Pembrolizumab + necitumumab (TRT=17)")          # Table S4 PD-1 +Other Pembrolizumab + necitumumab
    int_nivo_bev      <- fixed( 0.185)  ; label("ORR intercept: Nivolumab + bevacizumab (TRT=18)")              # Table S4 PD-1 +Other Nivolumab + bevacizumab
    int_nivo_ipi      <- fixed(-1.015)  ; label("ORR intercept: Nivolumab + ipilimumab (TRT=19)")               # Table S4 PD-1 +Other Nivolumab + ipilimumab
    int_pemb_cc486    <- fixed(-1.237)  ; label("ORR intercept: Pembrolizumab + CC-486 (TRT=20)")               # Table S4 PD-1 +Other Pembrolizumab + CC-486
    int_pemb_ipi      <- fixed(-0.879)  ; label("ORR intercept: Pembrolizumab + ipilimumab (TRT=21)")           # Table S4 PD-1 +Other Pembrolizumab + ipilimumab
    int_pemb_nk       <- fixed(-0.336)  ; label("ORR intercept: Pembrolizumab + NK immunotherapy (TRT=22)")     # Table S4 PD-1 +Other Pembrolizumab + NK immunotherapy
    int_dost_nira     <- fixed(-3.233)  ; label("ORR intercept: Dostarlimab + niraparib (TRT=23)")              # Table S4 PD-1 +Other Dostarlimab + niraparib
    int_camr_anlo     <- fixed(-2.332)  ; label("ORR intercept: Camrelizumab + anlotinib (TRT=24)")             # Table S4 PD-1 +Other Camrelizumab + anlotinib
    int_nivo_erlo     <- fixed(-1.834)  ; label("ORR intercept: Nivolumab + erlotinib (TRT=25)")                # Table S4 PD-1 +Other Nivolumab + erlotinib
    int_nivo_pegi     <- fixed(-1.648)  ; label("ORR intercept: Nivolumab + pegilodecakin (TRT=26)")            # Table S4 PD-1 +Other Nivolumab + pegilodecakin
    int_pemb_epac     <- fixed(-1.314)  ; label("ORR intercept: Pembrolizumab + epacadostat (TRT=27)")          # Table S4 PD-1 +Other Pembrolizumab + epacadostat
    int_pemb_enti     <- fixed(-1.618)  ; label("ORR intercept: Pembrolizumab + entinostat (TRT=28)")           # Table S4 PD-1 +Other Pembrolizumab + entinostat
    int_pemb_nira     <- fixed(-1.057)  ; label("ORR intercept: Pembrolizumab + niraparib (TRT=29)")            # Table S4 PD-1 +Other Pembrolizumab + niraparib
    int_pemb_pegi     <- fixed(-0.659)  ; label("ORR intercept: Pembrolizumab + pegilodecakin (TRT=30)")        # Table S4 PD-1 +Other Pembrolizumab + pegilodecakin
    int_pemb_ramu     <- fixed(-0.617)  ; label("ORR intercept: Pembrolizumab + ramucirumab (TRT=31)")          # Table S4 PD-1 +Other Pembrolizumab + ramucirumab
    int_camr_apat     <- fixed(-0.344)  ; label("ORR intercept: Camrelizumab + apatinib (TRT=32)")              # Table S4 PD-1 +Other Camrelizumab + apatinib
    int_avelumab      <- fixed(-1.418)  ; label("ORR intercept: Avelumab (PD-L1 mono; TRT=33)")                 # Table S4 PD-L1 Monotherapy Avelumab
    int_durvalumab    <- fixed(-1.441)  ; label("ORR intercept: Durvalumab (PD-L1 mono; TRT=34)")               # Table S4 PD-L1 Monotherapy Durvalumab
    int_atezolizumab  <- fixed(-1.624)  ; label("ORR intercept: Atezolizumab (PD-L1 mono; TRT=35)")             # Table S4 PD-L1 Monotherapy Atezolizumab
    int_atez_chemo    <- fixed(-0.464)  ; label("ORR intercept: Atezolizumab + chemo (TRT=36)")                 # Table S4 PD-L1 +Chemotherapy Atezolizumab + chemotherapy
    int_durv_chemo    <- fixed(-0.637)  ; label("ORR intercept: Durvalumab + chemo (TRT=37)")                   # Table S4 PD-L1 +Chemotherapy Durvalumab + chemotherapy
    int_avel_chemo    <- fixed(-0.380)  ; label("ORR intercept: Avelumab + chemo (TRT=38)")                     # Table S4 PD-L1 +Chemotherapy Avelumab + chemotherapy
    int_durv_trem_ch  <- fixed(-0.729)  ; label("ORR intercept: Durvalumab + tremelimumab + chemo (TRT=39)")    # Table S4 PD-L1 +Chemotherapy Durvalumab + tremelimumab + chemotherapy
    int_atez_tira     <- fixed(-0.706)  ; label("ORR intercept: Atezolizumab + tiragolumab (TRT=40)")           # Table S4 PD-L1 +Other Atezolizumab + tiragolumab
    int_durv_trem     <- fixed(-1.219)  ; label("ORR intercept: Durvalumab + tremelimumab (TRT=41)")            # Table S4 PD-L1 +Other Durvalumab + tremelimumab
    int_atez_dara     <- fixed(-2.770)  ; label("ORR intercept: Atezolizumab + daratumumab (TRT=42)")           # Table S4 PD-L1 +Other Atezolizumab + daratumumab
    int_atez_ipi      <- fixed(-0.610)  ; label("ORR intercept: Atezolizumab + ipilimumab (TRT=43)")            # Table S4 PD-L1 +Other Atezolizumab + ipilimumab
    int_ipi_pac       <- fixed(-0.850)  ; label("ORR intercept: Ipilimumab + paclitaxel (Other; TRT=44)")       # Table S4 Other treatments Ipilimumab + paclitaxel
    int_bevacizumab   <- fixed( 0.012)  ; label("ORR intercept: Bevacizumab (Other; TRT=45)")                   # Table S4 Other treatments Bevacizumab
    int_tremelimumab  <- fixed(-1.820)  ; label("ORR intercept: Tremelimumab (Other; TRT=46)")                  # Table S4 Other treatments Tremelimumab

    # -------- ORR covariate effects, Franzese 2026 Table 1 / Table S4 --------
    # Paper encodes -0.697 * (1 if >=2L, 0 if 1L); canonical LINE_1L flips sign of reference:
    # equivalent code below uses +e_line_1l_orr * LINE_1L so LINE_1L = 1 (first-line arm)
    # yields +0.697 relative to the paper's >=2L reference.
    e_line_1l_orr        <- fixed( 0.697)  ; label("ORR: LINE_1L effect (paper reference >=2L; sign flipped for canonical LINE_1L)")           # Table 1 ORR row: -0.697 * (1 if >=2L)
    e_pdl1_orr_pd1_quad  <- fixed( 1.736)  ; label("ORR: quadratic coefficient of (PDL1_TUM/100)^2 on PD-1 monotherapy arms (unitless)")       # Table 1 ORR row 'Average PD-L1 expression effect for PD-1 inhibitor' quadratic term
    e_pdl1_orr_pd1_lin   <- fixed( 0.274)  ; label("ORR: linear coefficient of (PDL1_TUM/100) on PD-1 monotherapy arms (unitless)")            # Table 1 ORR row 'Average PD-L1 expression effect for PD-1 inhibitor' linear term
    e_pdl1_orr_pdl1      <- fixed( 1.642)  ; label("ORR: linear coefficient of (PDL1_TUM/100) on PD-L1 monotherapy arms (unitless)")           # Table 1 ORR row 'Average PD-L1 expression effect for PD-L1 inhibitor'
    e_pdl1_orr_combo     <- fixed( 1.074)  ; label("ORR: linear coefficient of (PDL1_TUM/100) on any PD-(L)1 combination arms (unitless)")     # Table 1 ORR row 'Average PD-L1 expression effect for PD-(L)1 combination'
    e_squam_orr          <- fixed( 0.282)  ; label("ORR: linear coefficient of (TUMTP_SQUAM_PCT/100) on chemotherapy arms (unitless)")         # Table 1 ORR row '% Squamous histology'

    # -------- ORR study-strata random effect (Franzese 2026 Table 1) --------
    # eta_i ~ N(0, 0.346^2) => variance = 0.119716.
    eta_study_orr ~ 0.119716  # Table 1 ORR row: 'eta ~ (0, 0.346^2)' => 0.346^2 = 0.119716

    # ==========================================================================
    # OS MODEL (Franzese 2026 Table 1, Table 2, Table S2 monthly baseline hazards)
    # log(HR) = sum of treatment-category ORR-slope terms (with eta_study_os_orr
    #          entering only via the four PD-(L)1 or chemotherapy slopes -- the
    #          'other' slope has no random-effect term per Table 1)
    #        + squamous-on-chemotherapy shift + ECOG-0-on-non-chemotherapy shift
    #        + Asian-race-on-ORR-slope interaction
    #        + eta_study_os_int (random hazard-intercept shift per strata)
    # ==========================================================================

    # -------- OS treatment-type effects on ORR slope (Franzese 2026 Table 1, Table 2) --------
    e_orr_os_mono       <- fixed(-1.649)  ; label("OS: ORR-slope for PD-(L)1 monotherapy (Table 1, Table 2)")                          # Table 1 OS: '(eta_ORR - 1.649) * (%ORR/100) * (1 if PD-(L)1 monotherapy)'
    e_orr_os_chemo      <- fixed(-1.353)  ; label("OS: ORR-slope for PD-(L)1 + chemotherapy (Table 1, Table 2)")                       # Table 1 OS: '(eta_ORR - 1.353) * (%ORR/100) * (1 if PD-(L)1 plus chemotherapy)'
    e_orr_os_other      <- fixed(-1.507)  ; label("OS: ORR-slope for PD-(L)1 + 'other' therapy (Table 1, Table 2)")                    # Table 1 OS: '(eta_ORR - 1.507) * (%ORR/100) * (1 if PD-(L)1 plus other therapy)'
    e_orr_os_chemoonly  <- fixed(-1.651)  ; label("OS: ORR-slope for chemotherapy-only arms (Table 1, Table 2)")                       # Table 1 OS: '(eta_ORR - 1.651) * (%ORR/100) * (1 if chemotherapy)'
    e_orr_os_othertreat <- fixed(-1.285)  ; label("OS: ORR-slope for 'other' therapy arms (Table 1, Table 2)")                         # Table 1 OS: '(eta_ORR - 1.285) * (%ORR/100) * (1 if other therapy)'
    e_squam_os          <- fixed( 0.213)  ; label("OS: linear coefficient of (TUMTP_SQUAM_PCT/100) on chemotherapy arms (Table 2)")    # Table 2 OS: 'Shift in hazard for %squamous histology on chemotherapy in intercept'
    e_ecog_os           <- fixed(-0.400)  ; label("OS: linear coefficient of (PS_ECOG_0_PCT/100) on non-chemotherapy arms (Table 2)")  # Table 2 OS: 'Shift in hazard for %ECOG PS score=0 on non-chemotherapy related in intercept'
    e_asian_os          <- fixed(-0.595)  ; label("OS: additive shift to the ORR-slope for the (RACE_ASIAN_PCT/100) interaction (Table 2)")  # Table 2 OS: '%Asian race in ORR slope' (interaction: (eta_ORR - 0.595) * (%ORR/100) * (%Asian/100))

    # -------- OS study-strata random effects (Franzese 2026 Table 1) --------
    # (eta_int_os, eta_orr_os) ~ N(0, diag(0.011^2, 0.673^2)) -- diagonal (no correlation).
    eta_study_os_int + eta_study_os_orr ~ c(0.000121,
                                            0.0,      0.452929)  # Table 1 OS: 'eta_int ~ (0, 0.011^2), eta_ORR ~ (0, 0.673^2)' with 0 off-diagonal

    # -------- OS monthly baseline log(conditional death), Franzese 2026 Table S2 --------
    p_os_01 <- fixed(-2.743) ; label("OS log(P_m) month 1 (Table S2)")   # Table S2 OS month 1
    p_os_02 <- fixed(-2.502) ; label("OS log(P_m) month 2 (Table S2)")   # Table S2 OS month 2
    p_os_03 <- fixed(-2.511) ; label("OS log(P_m) month 3 (Table S2)")   # Table S2 OS month 3
    p_os_04 <- fixed(-2.449) ; label("OS log(P_m) month 4 (Table S2)")   # Table S2 OS month 4
    p_os_05 <- fixed(-2.379) ; label("OS log(P_m) month 5 (Table S2)")   # Table S2 OS month 5
    p_os_06 <- fixed(-2.402) ; label("OS log(P_m) month 6 (Table S2)")   # Table S2 OS month 6
    p_os_07 <- fixed(-2.388) ; label("OS log(P_m) month 7 (Table S2)")   # Table S2 OS month 7
    p_os_08 <- fixed(-2.410) ; label("OS log(P_m) month 8 (Table S2)")   # Table S2 OS month 8
    p_os_09 <- fixed(-2.307) ; label("OS log(P_m) month 9 (Table S2)")   # Table S2 OS month 9
    p_os_10 <- fixed(-2.408) ; label("OS log(P_m) month 10 (Table S2)")  # Table S2 OS month 10
    p_os_11 <- fixed(-2.411) ; label("OS log(P_m) month 11 (Table S2)")  # Table S2 OS month 11
    p_os_12 <- fixed(-2.426) ; label("OS log(P_m) month 12 (Table S2)")  # Table S2 OS month 12
    p_os_13 <- fixed(-2.400) ; label("OS log(P_m) month 13 (Table S2)")  # Table S2 OS month 13
    p_os_14 <- fixed(-2.409) ; label("OS log(P_m) month 14 (Table S2)")  # Table S2 OS month 14
    p_os_15 <- fixed(-2.421) ; label("OS log(P_m) month 15 (Table S2)")  # Table S2 OS month 15
    p_os_16 <- fixed(-2.482) ; label("OS log(P_m) month 16 (Table S2)")  # Table S2 OS month 16
    p_os_17 <- fixed(-2.447) ; label("OS log(P_m) month 17 (Table S2)")  # Table S2 OS month 17
    p_os_18 <- fixed(-2.554) ; label("OS log(P_m) month 18 (Table S2)")  # Table S2 OS month 18
    p_os_19 <- fixed(-2.490) ; label("OS log(P_m) month 19 (Table S2)")  # Table S2 OS month 19
    p_os_20 <- fixed(-2.578) ; label("OS log(P_m) month 20 (Table S2)")  # Table S2 OS month 20
    p_os_21 <- fixed(-2.671) ; label("OS log(P_m) month 21 (Table S2)")  # Table S2 OS month 21
    p_os_22 <- fixed(-2.677) ; label("OS log(P_m) month 22 (Table S2)")  # Table S2 OS month 22
    p_os_23 <- fixed(-2.761) ; label("OS log(P_m) month 23 (Table S2)")  # Table S2 OS month 23
    p_os_24 <- fixed(-2.814) ; label("OS log(P_m) month 24 (Table S2)")  # Table S2 OS month 24
    p_os_25 <- fixed(-2.642) ; label("OS log(P_m) month 25 (Table S2)")  # Table S2 OS month 25
    p_os_26 <- fixed(-2.743) ; label("OS log(P_m) month 26 (Table S2)")  # Table S2 OS month 26
    p_os_27 <- fixed(-2.790) ; label("OS log(P_m) month 27 (Table S2)")  # Table S2 OS month 27
    p_os_28 <- fixed(-2.912) ; label("OS log(P_m) month 28 (Table S2)")  # Table S2 OS month 28
    p_os_29 <- fixed(-2.886) ; label("OS log(P_m) month 29 (Table S2)")  # Table S2 OS month 29
    p_os_30 <- fixed(-2.884) ; label("OS log(P_m) month 30 (Table S2)")  # Table S2 OS month 30
    p_os_31 <- fixed(-2.883) ; label("OS log(P_m) month 31 (Table S2)")  # Table S2 OS month 31
    p_os_32 <- fixed(-2.392) ; label("OS log(P_m) month 32 (Table S2)")  # Table S2 OS month 32
    p_os_33 <- fixed(-3.156) ; label("OS log(P_m) month 33 (Table S2)")  # Table S2 OS month 33
    p_os_34 <- fixed(-2.923) ; label("OS log(P_m) month 34 (Table S2)")  # Table S2 OS month 34
    p_os_35 <- fixed(-3.061) ; label("OS log(P_m) month 35 (Table S2)")  # Table S2 OS month 35
    p_os_36 <- fixed(-3.444) ; label("OS log(P_m) month 36 (Table S2)")  # Table S2 OS month 36
    p_os_37 <- fixed(-3.239) ; label("OS log(P_m) month 37 (Table S2)")  # Table S2 OS month 37
    p_os_38 <- fixed(-3.266) ; label("OS log(P_m) month 38 (Table S2)")  # Table S2 OS month 38
    p_os_39 <- fixed(-3.235) ; label("OS log(P_m) month 39 (Table S2)")  # Table S2 OS month 39
    p_os_40 <- fixed(-3.471) ; label("OS log(P_m) month 40 (Table S2)")  # Table S2 OS month 40
    p_os_41 <- fixed(-3.471) ; label("OS log(P_m) month 41 (Table S2)")  # Table S2 OS month 41
    p_os_42 <- fixed(-3.471) ; label("OS log(P_m) month 42 (Table S2)")  # Table S2 OS month 42
    p_os_43 <- fixed(-2.982) ; label("OS log(P_m) month 43 (Table S2)")  # Table S2 OS month 43
    p_os_44 <- fixed(-2.982) ; label("OS log(P_m) month 44 (Table S2)")  # Table S2 OS month 44
    p_os_45 <- fixed(-2.982) ; label("OS log(P_m) month 45 (Table S2)")  # Table S2 OS month 45
    p_os_46 <- fixed(-3.497) ; label("OS log(P_m) month 46 (Table S2)")  # Table S2 OS month 46
    p_os_47 <- fixed(-3.497) ; label("OS log(P_m) month 47 (Table S2)")  # Table S2 OS month 47
    p_os_48 <- fixed(-3.497) ; label("OS log(P_m) month 48 (Table S2)")  # Table S2 OS month 48
    p_os_49 <- fixed(-3.383) ; label("OS log(P_m) month 49 (Table S2)")  # Table S2 OS month 49
    p_os_50 <- fixed(-3.383) ; label("OS log(P_m) month 50 (Table S2)")  # Table S2 OS month 50
    p_os_51 <- fixed(-3.383) ; label("OS log(P_m) month 51 (Table S2)")  # Table S2 OS month 51
    p_os_52 <- fixed(-3.469) ; label("OS log(P_m) month 52 (Table S2)")  # Table S2 OS month 52
    p_os_53 <- fixed(-3.469) ; label("OS log(P_m) month 53 (Table S2)")  # Table S2 OS month 53
    p_os_54 <- fixed(-3.469) ; label("OS log(P_m) month 54 (Table S2)")  # Table S2 OS month 54
    p_os_55 <- fixed(-3.879) ; label("OS log(P_m) month 55 (Table S2)")  # Table S2 OS month 55
    p_os_56 <- fixed(-3.879) ; label("OS log(P_m) month 56 (Table S2)")  # Table S2 OS month 56
    p_os_57 <- fixed(-3.879) ; label("OS log(P_m) month 57 (Table S2)")  # Table S2 OS month 57
    p_os_58 <- fixed(-3.451) ; label("OS log(P_m) month 58 (Table S2)")  # Table S2 OS month 58
    p_os_59 <- fixed(-3.451) ; label("OS log(P_m) month 59 (Table S2)")  # Table S2 OS month 59
    p_os_60 <- fixed(-3.451) ; label("OS log(P_m) month 60 (Table S2)")  # Table S2 OS month 60
    p_os_61_63 <- fixed(-2.957) ; label("OS log(P_m) months 61-63 grouped (Table S2)")  # Table S2 OS months 61-63
    p_os_64_66 <- fixed(-4.946) ; label("OS log(P_m) months 64-66 grouped (Table S2)")  # Table S2 OS months 64-66
    p_os_67_69 <- fixed(-2.827) ; label("OS log(P_m) months 67-69 grouped (Table S2)")  # Table S2 OS months 67-69
    p_os_70_72 <- fixed(-2.397) ; label("OS log(P_m) months 70-72 grouped (Table S2)")  # Table S2 OS months 70-72
    p_os_73_75 <- fixed(-2.252) ; label("OS log(P_m) months 73-75 grouped (Table S2)")  # Table S2 OS months 73-75

    # ==========================================================================
    # PFS MODEL (Franzese 2026 Table 1, Table 2, Table S2 monthly baseline hazards)
    # log(HR) = PD-(L)1-monotherapy shift + treatment-category ORR-slope terms
    #        + chemotherapy time-dependent shift (using Table S2 chemo column)
    #        + squamous-on-chemotherapy shift + ECOG-0 global shift
    #        + eta_study_pfs_int
    # Random effects eta_pfs_int, eta_pfs_orr are correlated (r = -0.642).
    # ==========================================================================

    # -------- PFS fixed effects (Franzese 2026 Table 1 and Table 2) --------
    e_pdl1mono_pfs      <- fixed( 0.410) ; label("PFS: PD-(L)1 monotherapy hazard-intercept shift (Table 1, Table 2)")                 # Table 1 PFS: '0.410 * (1 if PD-(L)1 monotherapy)'; Table 2 'Shift in hazard for PD-(L)1 monotherapy in intercept'
    e_orr_pfs_mono      <- fixed(-2.234) ; label("PFS: ORR-slope for PD-(L)1 monotherapy (Table 1, Table 2)")                          # Table 1 PFS: '(eta_ORR - 2.234) * (%ORR/100) * (1 if PD-(L)1 monotherapy)'
    e_orr_pfs_chemo     <- fixed(-1.278) ; label("PFS: ORR-slope for PD-(L)1 + chemotherapy (Table 1, Table 2)")                       # Table 1 PFS: '(eta_ORR - 1.278) * (%ORR/100) * (1 if PD-(L)1 plus chemotherapy)'
    e_orr_pfs_other     <- fixed(-1.192) ; label("PFS: ORR-slope for PD-(L)1 + 'other' therapy (Table 1, Table 2)")                    # Table 1 PFS: '(eta_ORR - 1.192) * (%ORR/100) * (1 if PD-(L)1 plus other therapy)'
    e_orr_pfs_chemoonly <- fixed(-2.217) ; label("PFS: ORR-slope for chemotherapy-only arms (Table 1, Table 2)")                       # Table 1 PFS: '(eta_ORR - 2.217) * (%ORR/100) * (1 if chemotherapy)'
    e_orr_pfs_othertreat<- fixed(-1.872) ; label("PFS: ORR-slope for 'other' therapy arms (Table 1)")                                  # Table 1 PFS: '(eta_ORR - 1.872) * (%ORR/100) * (1 if other therapy)'
    e_squam_pfs         <- fixed( 0.186) ; label("PFS: linear coefficient of (TUMTP_SQUAM_PCT/100) on chemotherapy arms (Table 2)")    # Table 2 PFS: 'Shift in hazard for %squamous histology on chemotherapy in intercept'
    e_ecog_pfs          <- fixed(-0.293) ; label("PFS: linear coefficient of (PS_ECOG_0_PCT/100) globally (Table 2)")                  # Table 2 PFS: 'Shift in hazard for %ECOG PS score=0 in intercept'

    # -------- PFS study-strata random effects (Franzese 2026 Table 1) --------
    # (eta_int_pfs, eta_orr_pfs) ~ N(0, [[0.199^2, r*0.199*0.251], [., 0.251^2]]) with r = -0.642.
    # Covariance = -0.642 * 0.199 * 0.251 = -0.032062 (rounded from paper's fractional r).
    eta_study_pfs_int + eta_study_pfs_orr ~ c(0.039601,
                                              -0.032063, 0.063001)  # Table 1 PFS: 'eta_int ~ (0, 0.199^2), eta_ORR ~ (0, 0.251^2), r = -0.642'; covariance = r * SD_int * SD_orr = -0.032063

    # -------- PFS monthly baseline log(conditional event), Franzese 2026 Table S2 --------
    p_pfs_01 <- fixed(-2.227) ; label("PFS log(P_m) month 1 (Table S2)")   # Table S2 PFS month 1
    p_pfs_02 <- fixed(-1.076) ; label("PFS log(P_m) month 2 (Table S2)")   # Table S2 PFS month 2
    p_pfs_03 <- fixed(-1.610) ; label("PFS log(P_m) month 3 (Table S2)")   # Table S2 PFS month 3
    p_pfs_04 <- fixed(-1.970) ; label("PFS log(P_m) month 4 (Table S2)")   # Table S2 PFS month 4
    p_pfs_05 <- fixed(-1.606) ; label("PFS log(P_m) month 5 (Table S2)")   # Table S2 PFS month 5
    p_pfs_06 <- fixed(-1.754) ; label("PFS log(P_m) month 6 (Table S2)")   # Table S2 PFS month 6
    p_pfs_07 <- fixed(-1.814) ; label("PFS log(P_m) month 7 (Table S2)")   # Table S2 PFS month 7
    p_pfs_08 <- fixed(-2.200) ; label("PFS log(P_m) month 8 (Table S2)")   # Table S2 PFS month 8
    p_pfs_09 <- fixed(-1.828) ; label("PFS log(P_m) month 9 (Table S2)")   # Table S2 PFS month 9
    p_pfs_10 <- fixed(-2.145) ; label("PFS log(P_m) month 10 (Table S2)")  # Table S2 PFS month 10
    p_pfs_11 <- fixed(-2.171) ; label("PFS log(P_m) month 11 (Table S2)")  # Table S2 PFS month 11
    p_pfs_12 <- fixed(-2.290) ; label("PFS log(P_m) month 12 (Table S2)")  # Table S2 PFS month 12
    p_pfs_13 <- fixed(-2.213) ; label("PFS log(P_m) month 13 (Table S2)")  # Table S2 PFS month 13
    p_pfs_14 <- fixed(-2.248) ; label("PFS log(P_m) month 14 (Table S2)")  # Table S2 PFS month 14
    p_pfs_15 <- fixed(-2.327) ; label("PFS log(P_m) month 15 (Table S2)")  # Table S2 PFS month 15
    p_pfs_16 <- fixed(-2.476) ; label("PFS log(P_m) month 16 (Table S2)")  # Table S2 PFS month 16
    p_pfs_17 <- fixed(-2.476) ; label("PFS log(P_m) month 17 (Table S2)")  # Table S2 PFS month 17
    p_pfs_18 <- fixed(-2.591) ; label("PFS log(P_m) month 18 (Table S2)")  # Table S2 PFS month 18
    p_pfs_19 <- fixed(-2.454) ; label("PFS log(P_m) month 19 (Table S2)")  # Table S2 PFS month 19
    p_pfs_20 <- fixed(-2.548) ; label("PFS log(P_m) month 20 (Table S2)")  # Table S2 PFS month 20
    p_pfs_21 <- fixed(-2.503) ; label("PFS log(P_m) month 21 (Table S2)")  # Table S2 PFS month 21
    p_pfs_22 <- fixed(-2.745) ; label("PFS log(P_m) month 22 (Table S2)")  # Table S2 PFS month 22
    p_pfs_23 <- fixed(-2.815) ; label("PFS log(P_m) month 23 (Table S2)")  # Table S2 PFS month 23
    p_pfs_24 <- fixed(-2.998) ; label("PFS log(P_m) month 24 (Table S2)")  # Table S2 PFS month 24
    p_pfs_25 <- fixed(-2.997) ; label("PFS log(P_m) month 25 (Table S2)")  # Table S2 PFS month 25
    p_pfs_26 <- fixed(-3.276) ; label("PFS log(P_m) month 26 (Table S2)")  # Table S2 PFS month 26
    p_pfs_27 <- fixed(-3.245) ; label("PFS log(P_m) month 27 (Table S2)")  # Table S2 PFS month 27
    p_pfs_28 <- fixed(-2.803) ; label("PFS log(P_m) month 28 (Table S2)")  # Table S2 PFS month 28
    p_pfs_29 <- fixed(-2.870) ; label("PFS log(P_m) month 29 (Table S2)")  # Table S2 PFS month 29
    p_pfs_30 <- fixed(-3.088) ; label("PFS log(P_m) month 30 (Table S2)")  # Table S2 PFS month 30
    p_pfs_31 <- fixed(-3.159) ; label("PFS log(P_m) month 31 (Table S2)")  # Table S2 PFS month 31
    p_pfs_32 <- fixed(-3.213) ; label("PFS log(P_m) month 32 (Table S2)")  # Table S2 PFS month 32
    p_pfs_33 <- fixed(-3.355) ; label("PFS log(P_m) month 33 (Table S2)")  # Table S2 PFS month 33
    p_pfs_34 <- fixed(-3.354) ; label("PFS log(P_m) month 34 (Table S2)")  # Table S2 PFS month 34
    p_pfs_35 <- fixed(-3.391) ; label("PFS log(P_m) month 35 (Table S2)")  # Table S2 PFS month 35
    p_pfs_36 <- fixed(-3.402) ; label("PFS log(P_m) month 36 (Table S2)")  # Table S2 PFS month 36
    p_pfs_37 <- fixed(-3.634) ; label("PFS log(P_m) month 37 (Table S2)")  # Table S2 PFS month 37
    p_pfs_38 <- fixed(-3.642) ; label("PFS log(P_m) month 38 (Table S2)")  # Table S2 PFS month 38
    p_pfs_39 <- fixed(-3.665) ; label("PFS log(P_m) month 39 (Table S2)")  # Table S2 PFS month 39
    p_pfs_40 <- fixed(-3.314) ; label("PFS log(P_m) month 40 (Table S2)")  # Table S2 PFS month 40
    p_pfs_41_43 <- fixed(-3.717) ; label("PFS log(P_m) months 41-43 grouped (Table S2)")  # Table S2 PFS months 41-43
    p_pfs_44_46 <- fixed(-3.245) ; label("PFS log(P_m) months 44-46 grouped (Table S2)")  # Table S2 PFS months 44-46
    p_pfs_47_49 <- fixed(-5.373) ; label("PFS log(P_m) months 47-49 grouped (Table S2)")  # Table S2 PFS months 47-49
    p_pfs_50_52 <- fixed(-2.750) ; label("PFS log(P_m) months 50-52 grouped (Table S2)")  # Table S2 PFS months 50-52
    p_pfs_53_55 <- fixed(-2.630) ; label("PFS log(P_m) months 53-55 grouped (Table S2)")  # Table S2 PFS months 53-55
    p_pfs_56_58 <- fixed(-3.057) ; label("PFS log(P_m) months 56-58 grouped (Table S2)")  # Table S2 PFS months 56-58

    # -------- PFS chemotherapy time-dependent shift log(HR), Franzese 2026 Table S2 (chemotherapy shift column) --------
    # Applied additively to log(HR) for chemotherapy arms only (per Table 1 PFS row: '+phi_k * (1 if chemotherapy)').
    # Month 1 has no shift (paper reports dash).
    ch_pfs_02 <- fixed(-0.160) ; label("PFS chemotherapy shift log(HR) month 2 (Table S2)")   # Table S2 PFS chemotherapy shift month 2
    ch_pfs_03 <- fixed( 0.279) ; label("PFS chemotherapy shift log(HR) month 3 (Table S2)")   # Table S2 PFS chemotherapy shift month 3
    ch_pfs_04 <- fixed( 0.443) ; label("PFS chemotherapy shift log(HR) month 4 (Table S2)")   # Table S2 PFS chemotherapy shift month 4
    ch_pfs_05 <- fixed( 0.774) ; label("PFS chemotherapy shift log(HR) month 5 (Table S2)")   # Table S2 PFS chemotherapy shift month 5
    ch_pfs_06 <- fixed( 0.845) ; label("PFS chemotherapy shift log(HR) month 6 (Table S2)")   # Table S2 PFS chemotherapy shift month 6
    ch_pfs_07 <- fixed( 0.923) ; label("PFS chemotherapy shift log(HR) month 7 (Table S2)")   # Table S2 PFS chemotherapy shift month 7
    ch_pfs_08 <- fixed( 1.045) ; label("PFS chemotherapy shift log(HR) month 8 (Table S2)")   # Table S2 PFS chemotherapy shift month 8
    ch_pfs_09 <- fixed( 0.929) ; label("PFS chemotherapy shift log(HR) month 9 (Table S2)")   # Table S2 PFS chemotherapy shift month 9
    ch_pfs_10 <- fixed( 0.869) ; label("PFS chemotherapy shift log(HR) month 10 (Table S2)")  # Table S2 PFS chemotherapy shift month 10
    ch_pfs_11 <- fixed( 0.820) ; label("PFS chemotherapy shift log(HR) month 11 (Table S2)")  # Table S2 PFS chemotherapy shift month 11
    ch_pfs_12 <- fixed( 0.815) ; label("PFS chemotherapy shift log(HR) month 12 (Table S2)")  # Table S2 PFS chemotherapy shift month 12
    ch_pfs_13 <- fixed( 0.732) ; label("PFS chemotherapy shift log(HR) month 13 (Table S2)")  # Table S2 PFS chemotherapy shift month 13
    ch_pfs_14 <- fixed( 0.694) ; label("PFS chemotherapy shift log(HR) month 14 (Table S2)")  # Table S2 PFS chemotherapy shift month 14
    ch_pfs_15 <- fixed( 0.763) ; label("PFS chemotherapy shift log(HR) month 15 (Table S2)")  # Table S2 PFS chemotherapy shift month 15
    ch_pfs_16 <- fixed( 0.990) ; label("PFS chemotherapy shift log(HR) month 16 (Table S2)")  # Table S2 PFS chemotherapy shift month 16
    ch_pfs_17 <- fixed( 1.021) ; label("PFS chemotherapy shift log(HR) month 17 (Table S2)")  # Table S2 PFS chemotherapy shift month 17
    ch_pfs_18 <- fixed( 1.206) ; label("PFS chemotherapy shift log(HR) month 18 (Table S2)")  # Table S2 PFS chemotherapy shift month 18
    ch_pfs_19 <- fixed( 0.765) ; label("PFS chemotherapy shift log(HR) month 19 (Table S2)")  # Table S2 PFS chemotherapy shift month 19
    ch_pfs_20 <- fixed( 0.803) ; label("PFS chemotherapy shift log(HR) month 20 (Table S2)")  # Table S2 PFS chemotherapy shift month 20
    ch_pfs_21 <- fixed( 0.920) ; label("PFS chemotherapy shift log(HR) month 21 (Table S2)")  # Table S2 PFS chemotherapy shift month 21
    ch_pfs_22 <- fixed( 0.839) ; label("PFS chemotherapy shift log(HR) month 22 (Table S2)")  # Table S2 PFS chemotherapy shift month 22
    ch_pfs_23 <- fixed( 0.963) ; label("PFS chemotherapy shift log(HR) month 23 (Table S2)")  # Table S2 PFS chemotherapy shift month 23
    ch_pfs_24 <- fixed( 1.002) ; label("PFS chemotherapy shift log(HR) month 24 (Table S2)")  # Table S2 PFS chemotherapy shift month 24
    ch_pfs_25 <- fixed( 1.354) ; label("PFS chemotherapy shift log(HR) month 25 (Table S2)")  # Table S2 PFS chemotherapy shift month 25
    ch_pfs_26 <- fixed( 1.515) ; label("PFS chemotherapy shift log(HR) month 26 (Table S2)")  # Table S2 PFS chemotherapy shift month 26
    ch_pfs_27 <- fixed( 1.390) ; label("PFS chemotherapy shift log(HR) month 27 (Table S2)")  # Table S2 PFS chemotherapy shift month 27
    ch_pfs_28 <- fixed( 0.661) ; label("PFS chemotherapy shift log(HR) month 28 (Table S2; RSE 101%)")  # Table S2 PFS chemotherapy shift month 28
    ch_pfs_29 <- fixed( 0.668) ; label("PFS chemotherapy shift log(HR) month 29 (Table S2; RSE 108%)")  # Table S2 PFS chemotherapy shift month 29
    ch_pfs_30_40 <- fixed( 0.805) ; label("PFS chemotherapy shift log(HR) months 30-40 grouped (Table S2)")  # Table S2 PFS chemotherapy shift months 30-40 (same 0.805 value across grouped months per Table S2 RSE 49%)
    ch_pfs_41_56 <- fixed( 1.234) ; label("PFS chemotherapy shift log(HR) months 41-56 grouped (Table S2)")  # Table S2 PFS chemotherapy shift months 41-56 (same 1.234 value across grouped months per Table S2 RSE 37%)

    # ==========================================================================
    # Residual error (nominal). The MBMA fits observations weighted by binomial
    # sample sizes (see Franzese 2026 Equations 3 and 5). For simulation-only
    # use we set a small nominal SD on the Cc output so nlmixr2 has an error
    # model, but no scalar residual SD is reported by the paper.
    # ==========================================================================
    addSd <- fixed(0.01) ; label("Nominal additive residual SD on Cc (dimensionless ORR proportion). Paper does not report a single scalar residual SD; observations are weighted by binomial sample sizes.")  # Franzese 2026 Methods 2.3 (residuals weighted per Equations 3 and 5)
  })

  model({
    # ============ Treatment-type indicators derived from TRT integer ============
    # ORR-model PD-(L)1-class indicators (Table 1 ORR row):
    is_pd1_mono   <- (TRT >= 2)  * (TRT <= 9)                            # TRT 2-9   are PD-1 monotherapy
    is_pdl1_mono  <- (TRT == 33) + (TRT == 34) + (TRT == 35)             # TRT 33-35 are PD-L1 monotherapy
    is_pdl1_combo <- ((TRT >= 10) * (TRT <= 16)) +                       # TRT 10-16 PD-1 + chemo
                     ((TRT >= 17) * (TRT <= 32)) +                       # TRT 17-32 PD-1 + other
                     ((TRT >= 36) * (TRT <= 39)) +                       # TRT 36-39 PD-L1 + chemo
                     ((TRT >= 40) * (TRT <= 43))                         # TRT 40-43 PD-L1 + other
    is_chemo_only <- (TRT == 1)                                          # TRT 1     Chemotherapy

    # OS/PFS 5-treatment-category indicators (Table 1 OS/PFS rows):
    is_pdl1_mono_cat  <- is_pd1_mono + is_pdl1_mono                                 # cat 1
    is_pdl1_chemo_cat <- ((TRT >= 10) * (TRT <= 16)) +                              # cat 2
                        ((TRT >= 36) * (TRT <= 39))
    is_pdl1_other_cat <- ((TRT >= 17) * (TRT <= 32)) +                              # cat 3
                        ((TRT >= 40) * (TRT <= 43))
    is_chemo_cat      <- (TRT == 1)                                                 # cat 4
    is_other_cat      <- (TRT >= 44) * (TRT <= 46)                                  # cat 5

    # Non-chemotherapy indicator (for the OS ECOG-0 effect, which applies to any non-chemotherapy arm).
    is_nonchemo <- 1 - is_chemo_cat

    # ============ ORR MODEL ============
    # Select the treatment intercept for this arm (arithmetic-indicator form; matches Struemper 2025 precedent).
    intercept_treatment <-
      (TRT ==  1) * int_chemo         + (TRT ==  2) * int_camrelizumab  +
      (TRT ==  3) * int_cemiplimab    + (TRT ==  4) * int_nivolumab     +
      (TRT ==  5) * int_pembrolizumab + (TRT ==  6) * int_toripalimab   +
      (TRT ==  7) * int_dostarlimab   + (TRT ==  8) * int_budigalimab   +
      (TRT ==  9) * int_retifanlimab  + (TRT == 10) * int_camr_chemo    +
      (TRT == 11) * int_nivo_chemo    + (TRT == 12) * int_pemb_chemo    +
      (TRT == 13) * int_sint_chemo    + (TRT == 14) * int_tisl_chemo    +
      (TRT == 15) * int_dost_chemo    + (TRT == 16) * int_tori_chemo    +
      (TRT == 17) * int_pemb_necit    + (TRT == 18) * int_nivo_bev      +
      (TRT == 19) * int_nivo_ipi      + (TRT == 20) * int_pemb_cc486    +
      (TRT == 21) * int_pemb_ipi      + (TRT == 22) * int_pemb_nk       +
      (TRT == 23) * int_dost_nira     + (TRT == 24) * int_camr_anlo     +
      (TRT == 25) * int_nivo_erlo     + (TRT == 26) * int_nivo_pegi     +
      (TRT == 27) * int_pemb_epac     + (TRT == 28) * int_pemb_enti     +
      (TRT == 29) * int_pemb_nira     + (TRT == 30) * int_pemb_pegi     +
      (TRT == 31) * int_pemb_ramu     + (TRT == 32) * int_camr_apat     +
      (TRT == 33) * int_avelumab      + (TRT == 34) * int_durvalumab    +
      (TRT == 35) * int_atezolizumab  + (TRT == 36) * int_atez_chemo    +
      (TRT == 37) * int_durv_chemo    + (TRT == 38) * int_avel_chemo    +
      (TRT == 39) * int_durv_trem_ch  + (TRT == 40) * int_atez_tira     +
      (TRT == 41) * int_durv_trem     + (TRT == 42) * int_atez_dara     +
      (TRT == 43) * int_atez_ipi      + (TRT == 44) * int_ipi_pac       +
      (TRT == 45) * int_bevacizumab   + (TRT == 46) * int_tremelimumab

    # Franzese 2026 Table 1 ORR row: full logit(ORR) equation with study-strata eta.
    # Note: paper's line-of-therapy encoding is -0.697 * (1 if >=2L), inverted to +e_line_1l_orr * LINE_1L
    # using the canonical LINE_1L (1 = first-line arm) so the same effect direction is preserved.
    pdl1_frac <- PDL1_TUM / 100
    logit_orr <- intercept_treatment +
                 e_line_1l_orr * LINE_1L +
                 (e_pdl1_orr_pd1_quad * pdl1_frac + e_pdl1_orr_pd1_lin) * pdl1_frac * is_pd1_mono +
                 e_pdl1_orr_pdl1  * pdl1_frac * is_pdl1_mono +
                 e_pdl1_orr_combo * pdl1_frac * is_pdl1_combo +
                 e_squam_orr * (TUMTP_SQUAM_PCT / 100) * is_chemo_only +
                 eta_study_orr

    # Predicted per-arm ORR (proportion) via inverse logit.
    orr_pred <- exp(logit_orr) / (1 + exp(logit_orr))

    # ============ OS MODEL log(HR) ============
    # Per Table 1 OS row: eta_orr_os enters INSIDE each of the four PD-(L)1 or chemotherapy
    # ORR-slope terms and inside the Asian interaction; the 'other therapy' slope has NO
    # eta_orr_os per Table 1. eta_int_os enters as a strata-level hazard-intercept shift.
    log_hr_os <- (eta_study_os_orr + e_orr_os_mono)       * orr_pred * is_pdl1_mono_cat +
                 (eta_study_os_orr + e_orr_os_chemo)      * orr_pred * is_pdl1_chemo_cat +
                 (eta_study_os_orr + e_orr_os_other)      * orr_pred * is_pdl1_other_cat +
                 (eta_study_os_orr + e_orr_os_chemoonly)  * orr_pred * is_chemo_cat +
                 (                    e_orr_os_othertreat) * orr_pred * is_other_cat +
                 e_squam_os * (TUMTP_SQUAM_PCT / 100) * is_chemo_cat +
                 e_ecog_os  * (PS_ECOG_0_PCT / 100)   * is_nonchemo +
                 (eta_study_os_orr + e_asian_os) * orr_pred * (RACE_ASIAN_PCT / 100) +
                 eta_study_os_int

    # ============ PFS MODEL log(HR) (excluding time-dependent chemotherapy shift) ============
    # Per Table 1 PFS row: 0.410 shift on PD-(L)1 monotherapy in intercept + treatment-category
    # ORR-slope terms (eta_orr_pfs enters inside the PD-(L)1 and chemotherapy slopes; the 'other'
    # slope has no eta_orr_pfs). eta_int_pfs enters as a strata-level hazard-intercept shift.
    # The chemotherapy time-dependent shift phi_k is added month-by-month via the lookup below.
    log_hr_pfs_static <- e_pdl1mono_pfs * is_pdl1_mono_cat +
                         (eta_study_pfs_orr + e_orr_pfs_mono)       * orr_pred * is_pdl1_mono_cat +
                         (eta_study_pfs_orr + e_orr_pfs_chemo)      * orr_pred * is_pdl1_chemo_cat +
                         (eta_study_pfs_orr + e_orr_pfs_other)      * orr_pred * is_pdl1_other_cat +
                         (eta_study_pfs_orr + e_orr_pfs_chemoonly)  * orr_pred * is_chemo_cat +
                         (                    e_orr_pfs_othertreat) * orr_pred * is_other_cat +
                         e_squam_pfs * (TUMTP_SQUAM_PCT / 100) * is_chemo_cat +
                         e_ecog_pfs  * (PS_ECOG_0_PCT / 100) +
                         eta_study_pfs_int

    # ============ Monthly baseline hazard lookup (piecewise-constant on t in months) ============
    # OS log(P_m) selector: 1 <= t <= 1 -> month 1; t <= 2 -> month 2; ...
    # Above the 75-month support, the last grouped estimate (73-75) is held.
    log_p_os_now <-
      (t <= 1) * p_os_01 +
      (t > 1)  * (t <= 2)  * p_os_02 + (t > 2)  * (t <= 3)  * p_os_03 +
      (t > 3)  * (t <= 4)  * p_os_04 + (t > 4)  * (t <= 5)  * p_os_05 +
      (t > 5)  * (t <= 6)  * p_os_06 + (t > 6)  * (t <= 7)  * p_os_07 +
      (t > 7)  * (t <= 8)  * p_os_08 + (t > 8)  * (t <= 9)  * p_os_09 +
      (t > 9)  * (t <= 10) * p_os_10 + (t > 10) * (t <= 11) * p_os_11 +
      (t > 11) * (t <= 12) * p_os_12 + (t > 12) * (t <= 13) * p_os_13 +
      (t > 13) * (t <= 14) * p_os_14 + (t > 14) * (t <= 15) * p_os_15 +
      (t > 15) * (t <= 16) * p_os_16 + (t > 16) * (t <= 17) * p_os_17 +
      (t > 17) * (t <= 18) * p_os_18 + (t > 18) * (t <= 19) * p_os_19 +
      (t > 19) * (t <= 20) * p_os_20 + (t > 20) * (t <= 21) * p_os_21 +
      (t > 21) * (t <= 22) * p_os_22 + (t > 22) * (t <= 23) * p_os_23 +
      (t > 23) * (t <= 24) * p_os_24 + (t > 24) * (t <= 25) * p_os_25 +
      (t > 25) * (t <= 26) * p_os_26 + (t > 26) * (t <= 27) * p_os_27 +
      (t > 27) * (t <= 28) * p_os_28 + (t > 28) * (t <= 29) * p_os_29 +
      (t > 29) * (t <= 30) * p_os_30 + (t > 30) * (t <= 31) * p_os_31 +
      (t > 31) * (t <= 32) * p_os_32 + (t > 32) * (t <= 33) * p_os_33 +
      (t > 33) * (t <= 34) * p_os_34 + (t > 34) * (t <= 35) * p_os_35 +
      (t > 35) * (t <= 36) * p_os_36 + (t > 36) * (t <= 37) * p_os_37 +
      (t > 37) * (t <= 38) * p_os_38 + (t > 38) * (t <= 39) * p_os_39 +
      (t > 39) * (t <= 40) * p_os_40 + (t > 40) * (t <= 41) * p_os_41 +
      (t > 41) * (t <= 42) * p_os_42 + (t > 42) * (t <= 43) * p_os_43 +
      (t > 43) * (t <= 44) * p_os_44 + (t > 44) * (t <= 45) * p_os_45 +
      (t > 45) * (t <= 46) * p_os_46 + (t > 46) * (t <= 47) * p_os_47 +
      (t > 47) * (t <= 48) * p_os_48 + (t > 48) * (t <= 49) * p_os_49 +
      (t > 49) * (t <= 50) * p_os_50 + (t > 50) * (t <= 51) * p_os_51 +
      (t > 51) * (t <= 52) * p_os_52 + (t > 52) * (t <= 53) * p_os_53 +
      (t > 53) * (t <= 54) * p_os_54 + (t > 54) * (t <= 55) * p_os_55 +
      (t > 55) * (t <= 56) * p_os_56 + (t > 56) * (t <= 57) * p_os_57 +
      (t > 57) * (t <= 58) * p_os_58 + (t > 58) * (t <= 59) * p_os_59 +
      (t > 59) * (t <= 60) * p_os_60 + (t > 60) * (t <= 63) * p_os_61_63 +
      (t > 63) * (t <= 66) * p_os_64_66 + (t > 66) * (t <= 69) * p_os_67_69 +
      (t > 69) * (t <= 72) * p_os_70_72 + (t > 72)             * p_os_73_75

    # PFS log(P_m) selector: analogous to OS but ending at month 58 with grouped tail.
    log_p_pfs_now <-
      (t <= 1) * p_pfs_01 +
      (t > 1)  * (t <= 2)  * p_pfs_02 + (t > 2)  * (t <= 3)  * p_pfs_03 +
      (t > 3)  * (t <= 4)  * p_pfs_04 + (t > 4)  * (t <= 5)  * p_pfs_05 +
      (t > 5)  * (t <= 6)  * p_pfs_06 + (t > 6)  * (t <= 7)  * p_pfs_07 +
      (t > 7)  * (t <= 8)  * p_pfs_08 + (t > 8)  * (t <= 9)  * p_pfs_09 +
      (t > 9)  * (t <= 10) * p_pfs_10 + (t > 10) * (t <= 11) * p_pfs_11 +
      (t > 11) * (t <= 12) * p_pfs_12 + (t > 12) * (t <= 13) * p_pfs_13 +
      (t > 13) * (t <= 14) * p_pfs_14 + (t > 14) * (t <= 15) * p_pfs_15 +
      (t > 15) * (t <= 16) * p_pfs_16 + (t > 16) * (t <= 17) * p_pfs_17 +
      (t > 17) * (t <= 18) * p_pfs_18 + (t > 18) * (t <= 19) * p_pfs_19 +
      (t > 19) * (t <= 20) * p_pfs_20 + (t > 20) * (t <= 21) * p_pfs_21 +
      (t > 21) * (t <= 22) * p_pfs_22 + (t > 22) * (t <= 23) * p_pfs_23 +
      (t > 23) * (t <= 24) * p_pfs_24 + (t > 24) * (t <= 25) * p_pfs_25 +
      (t > 25) * (t <= 26) * p_pfs_26 + (t > 26) * (t <= 27) * p_pfs_27 +
      (t > 27) * (t <= 28) * p_pfs_28 + (t > 28) * (t <= 29) * p_pfs_29 +
      (t > 29) * (t <= 30) * p_pfs_30 + (t > 30) * (t <= 31) * p_pfs_31 +
      (t > 31) * (t <= 32) * p_pfs_32 + (t > 32) * (t <= 33) * p_pfs_33 +
      (t > 33) * (t <= 34) * p_pfs_34 + (t > 34) * (t <= 35) * p_pfs_35 +
      (t > 35) * (t <= 36) * p_pfs_36 + (t > 36) * (t <= 37) * p_pfs_37 +
      (t > 37) * (t <= 38) * p_pfs_38 + (t > 38) * (t <= 39) * p_pfs_39 +
      (t > 39) * (t <= 40) * p_pfs_40 + (t > 40) * (t <= 43) * p_pfs_41_43 +
      (t > 43) * (t <= 46) * p_pfs_44_46 + (t > 46) * (t <= 49) * p_pfs_47_49 +
      (t > 49) * (t <= 52) * p_pfs_50_52 + (t > 52) * (t <= 55) * p_pfs_53_55 +
      (t > 55)             * p_pfs_56_58

    # PFS chemotherapy time-dependent shift phi_k (Table S2 chemo column).
    # Month 1 has no shift (Table S2 dash); grouped months 30-40 = 0.805, months 41-56 = 1.234.
    chemo_shift_now <-
      (t > 1)  * (t <= 2)  * ch_pfs_02 + (t > 2)  * (t <= 3)  * ch_pfs_03 +
      (t > 3)  * (t <= 4)  * ch_pfs_04 + (t > 4)  * (t <= 5)  * ch_pfs_05 +
      (t > 5)  * (t <= 6)  * ch_pfs_06 + (t > 6)  * (t <= 7)  * ch_pfs_07 +
      (t > 7)  * (t <= 8)  * ch_pfs_08 + (t > 8)  * (t <= 9)  * ch_pfs_09 +
      (t > 9)  * (t <= 10) * ch_pfs_10 + (t > 10) * (t <= 11) * ch_pfs_11 +
      (t > 11) * (t <= 12) * ch_pfs_12 + (t > 12) * (t <= 13) * ch_pfs_13 +
      (t > 13) * (t <= 14) * ch_pfs_14 + (t > 14) * (t <= 15) * ch_pfs_15 +
      (t > 15) * (t <= 16) * ch_pfs_16 + (t > 16) * (t <= 17) * ch_pfs_17 +
      (t > 17) * (t <= 18) * ch_pfs_18 + (t > 18) * (t <= 19) * ch_pfs_19 +
      (t > 19) * (t <= 20) * ch_pfs_20 + (t > 20) * (t <= 21) * ch_pfs_21 +
      (t > 21) * (t <= 22) * ch_pfs_22 + (t > 22) * (t <= 23) * ch_pfs_23 +
      (t > 23) * (t <= 24) * ch_pfs_24 + (t > 24) * (t <= 25) * ch_pfs_25 +
      (t > 25) * (t <= 26) * ch_pfs_26 + (t > 26) * (t <= 27) * ch_pfs_27 +
      (t > 27) * (t <= 28) * ch_pfs_28 + (t > 28) * (t <= 29) * ch_pfs_29 +
      (t > 29) * (t <= 40) * ch_pfs_30_40 + (t > 40) * (t <= 56) * ch_pfs_41_56

    # Chemotherapy shift applies only to chemotherapy arms (paper Table 1 PFS row).
    log_hr_pfs <- log_hr_pfs_static + chemo_shift_now * is_chemo_cat

    # ============ Cumulative-hazard ODEs and survival curves ============
    # Convert monthly discrete conditional-event log-probability P_m to a piecewise-constant
    # instantaneous hazard rate h_m = -log(1 - exp(P_m)) per month (paper Table 1 defines
    # log(S_0(t_k)) as the sum of log(1 - exp(P_m)) contributions, so exp(-cumhaz) reproduces
    # the paper's baseline survival). Multiply by exp(log_hr_*) for proportional-hazards HR.
    haz_os_baseline  <- -log(1 - exp(log_p_os_now)  + 1e-30)
    haz_pfs_baseline <- -log(1 - exp(log_p_pfs_now) + 1e-30)

    haz_os  <- haz_os_baseline  * exp(log_hr_os)
    haz_pfs <- haz_pfs_baseline * exp(log_hr_pfs)

    d/dt(cumhaz_os)  <- haz_os
    d/dt(cumhaz_pfs) <- haz_pfs

    S_os  <- exp(-cumhaz_os)
    S_pfs <- exp(-cumhaz_pfs)

    # ============ Observation: predicted per-arm ORR ============
    # Cc is the per-arm predicted objective response fraction (dimensionless proportion).
    # S_os and S_pfs are emitted as additional model outputs alongside Cc so the vignette
    # can plot per-arm survival curves.
    Cc <- orr_pred

    Cc ~ add(addSd)
  })
}
attr(Franzese_2026_pdl1_nsclc_mbma, "message") <-
  "MBMA of PD-(L)1 immunotherapy in metastatic non-small cell lung cancer (Franzese 2026, 114 studies, 46 unique treatments). Integrated ORR (mixed-effects logistic regression) + OS + PFS (semi-parametric proportional hazards with monthly discrete baseline hazard from Franzese 2026 Table S2). Inputs: TRT integer 1-46 (see covariateData[[TRT]]$notes for the mapping), PDL1_TUM (%), LINE_1L (binary), TUMTP_SQUAM_PCT (%), PS_ECOG_0_PCT (%), RACE_ASIAN_PCT (%). Outputs: Cc = per-arm predicted ORR proportion; S_os = per-arm S_OS(t) survival curve; S_pfs = per-arm S_PFS(t) survival curve. Simulation scope: per-arm summary; NOT suitable for individual-subject trajectories. All parameter values wrapped in fixed() because the model is a downstream user of the published fit, not a re-estimation of it."
Franzese_2026_pdl1_nsclc_mbma
