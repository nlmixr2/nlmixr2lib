Abegesah_2025_durvalumab <- function() {
  description <- "Two-compartment population PK model for durvalumab (anti-PD-L1 IgG1 kappa) with sigmoidal time-varying clearance in adults with advanced solid tumours, updating the POSEIDON model with TOPAZ-1 biliary tract cancer patients treated with durvalumab plus gemcitabine/cisplatin (Abegesah 2025)"
  reference <- paste(
    "Abegesah A, Oh D-Y, Lim K, Fan C, Chen C, Kim C, Wang J, Xynos I,",
    "Zotkiewicz M, Ren S, Phipps A, Gibbs M, Zhou D. Population",
    "pharmacokinetics and exposure-response analysis of durvalumab in",
    "combination with gemcitabine and cisplatin in patients with advanced",
    "biliary tract cancer. Cancer Chemother Pharmacol. 2025;95:23.",
    "doi:10.1007/s00280-024-04743-8",
    sep = " "
  )
  vignette <- "Abegesah_2025_durvalumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Abegesah 2025 "Bioanalytical methods" states that
  # durvalumab was measured in SERUM (ECL assay, LLOQ-ULOQ 0.05-3.2 ug/mL
  # on the assayed dilution), so the specimen is serum rather than plasma.
  compartmentData <- list(
    central     = list(analyte = "durvalumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "durvalumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on CL (exponent 0.338) and Vc (exponent 0.515), both normalized to 69.4 kg. The normalizer is printed inside the CL and Vc equations on page 3 of Abegesah 2025 and equals the 'Previous 5 studies' median of Table 2 (69.4 kg), i.e. the reference value inherited from the POSEIDON model, NOT the pooled 6-study median (69.0 kg). Do not conflate the two.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent -0.589) normalized to 39 g/L, printed in the CL_cont.cov equation on page 3. Reported and used in SI g/L (Table 2 median 39.0 g/L), so no g/dL conversion applies. The same 39 g/L reference is used by the sibling AstraZeneca model Hwang_2022_tremelimumab.R.",
      source_name        = "alb"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (raw Cockcroft-Gault style, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.136) normalized to 85.66 mL/min, printed in the CL_cont.cov equation on page 3. Table 2 reports 'Creatinine clearance (mL/min)' with no BSA normalization, so the per-model unit is mL/min rather than the register default mL/min/1.73 m^2 (same precedent as the sibling durvalumab model deVries_2025_durvalumab.R, which carries the essentially identical Baverel 2018 normalizer of 85.65 mL/min).",
      source_name        = "CrCL"
    ),
    LDH = list(
      description        = "Baseline lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.0515) normalized to 247 U/L, printed in the CL_cont.cov equation on page 3. As with body weight, 247 U/L is the 'Previous 5 studies' median of Table 2 (the inherited POSEIDON reference), not the pooled 6-study median of 242 U/L.",
      source_name        = "LDH"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects for females on CL (1 - 0.161) and Vc (1 - 0.140), printed in the CL_cat.cov and Vc equations on page 3 with male as the explicit reference (the '1_male' factor). Fig. 1 corroborates both magnitudes and directions (-16.1% on CL, -14% on Vc). Direction agrees with the sibling durvalumab models Ogasawara_2020_durvalumab.R (female CL 0.791) and deVries_2025_durvalumab.R (female CL 0.857).",
      source_name        = "sex"
    ),
    ECOG_GE1 = list(
      description        = "Baseline ECOG performance status >= 1 (1 = ECOG 1 or worse, 0 = ECOG 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes              = "Multiplicative effect on CL of (1 - 0.0501), printed in the CL_cat.cov equation on page 3 as the paper's binary 'ECOGbin' with ECOGbin = 0 the explicit reference. Table 2 records ECOG 0/1/2 as 37.91/61.83/0.10 percent, so the ordinal score collapses to this binary indicator. NOTE: the Fig. 1 tornado plot draws 'ECOG restricted activity' as +5.3% (a CL INCREASE), which contradicts the sign of both the printed equation and the Table 4 estimate (-0.0501); the equation is followed here per the standing text-vs-equation policy, and the negative direction is independently corroborated by the upstream Baverel 2018 durvalumab model as transcribed in deVries_2025_durvalumab.R (0.937^ECOG_GE1, i.e. 6.3% LOWER CL). See the vignette Errata.",
      source_name        = "ECOGbin"
    ),
    CONMED_CHEMO = list(
      description        = "Durvalumab co-administered with platinum-based chemotherapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (durvalumab monotherapy)",
      notes              = "The paper's 'comb' column is a three-level categorical (comb = 0 durvalumab monotherapy, comb = 1 durvalumab + chemotherapy, comb = 2 durvalumab + tremelimumab + chemotherapy) carrying multiplicative CL factors of 1, (1 - 0.163) and (1 - 0.0929) respectively. It is stored here as two canonical columns, CONMED_CHEMO and CONMED_TREMELIMUMAB, so that each column keeps its register meaning; the model reconstructs the comb = 1 stratum as CONMED_CHEMO * (1 - CONMED_TREMELIMUMAB). TOPAZ-1 (durvalumab + gemcitabine/cisplatin) is the comb = 1 stratum. NOTE: the Table 4 legend glosses 'Comb1' as 'durvalumab and tremelimumab without chemotherapy', which is wrong -- Fig. 1 labels the -16.3% bar 'Combo Durva + Chemo', and Table 2 assigns all 314 TOPAZ-1 durvalumab + gem/cis patients to the second combination level. Both of the paper's own data sources agree against its legend. See the vignette Errata.",
      source_name        = "comb"
    ),
    CONMED_TREMELIMUMAB = list(
      description        = "Durvalumab co-administered with tremelimumab (always on a chemotherapy backbone in this dataset)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no tremelimumab co-administration)",
      notes              = "Identifies the paper's comb = 2 stratum (durvalumab + tremelimumab + chemotherapy, n = 327, the POSEIDON triplet arm per Table 2), carrying a multiplicative CL factor of (1 - 0.0929). Corroborated by the Fig. 1 bar labelled 'Combo Durva + Treme + Chemo' at -9.3%. Because no tremelimumab-without-chemotherapy patients exist in this pooled dataset, CONMED_TREMELIMUMAB = 1 implies CONMED_CHEMO = 1 and this column alone selects the comb = 2 factor. This is the mirror image of the sibling model Hwang_2022_tremelimumab.R, where tremelimumab is the analyte and COMBO_DURVA marks durvalumab co-administration.",
      source_name        = "comb"
    ),
    TUMTP_BLADDER = list(
      description        = "Bladder (urothelial) carcinoma tumor-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types; the pooled reference stratum is non-small-cell lung cancer, the paper's tumtyp = 0)",
      notes              = "Multiplicative effect on CL of (1 + 0.0698). Table 4 names this row only 'Tumor type 2 on CL'; the tumor-type identity is recovered from the Fig. 1 tornado plot, which labels a +7% clearance bar 'TUM TYP Bladder' -- and 0.0698 is the only tumor-type coefficient that gives +7.0%. Bladder patients enter through the Study 1108 advanced-solid-tumour cohort (durvalumab's first-approved indication was urothelial carcinoma); Table 2's coarser 'Primary indication' column does not resolve them separately.",
      source_name        = "tumtyp"
    ),
    TUMTP_BTC = list(
      description        = "Biliary tract cancer tumor-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types; the pooled reference stratum is non-small-cell lung cancer, the paper's tumtyp = 0)",
      notes              = "Multiplicative effect on CL of (1 + 0.166), the largest tumor-type effect in the model. Table 4 names this row only 'Tumor type 3 on CL'; Fig. 1 labels the +16.6% bar 'TUM TYP BTC', and 0.166 gives exactly +16.6%. This is the TOPAZ-1 stratum (n = 314, 10.02% of the pooled analysis dataset) and the covariate added by this analysis -- the Results state that tumor type on clearance was 'the only new covariate that was added' to the inherited POSEIDON model.",
      source_name        = "tumtyp"
    ),
    TUMTP_OTHER = list(
      description        = "Residual 'other tumor type' indicator (the paper's unlabelled tumtyp = 1 stratum)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all other tumor types; the pooled reference stratum is non-small-cell lung cancer, the paper's tumtyp = 0)",
      notes              = "Multiplicative effect on CL of (1 - 0.0101). UNIDENTIFIED STRATUM: Abegesah 2025 never labels tumtyp = 1. Table 4 calls it only 'Tumor type 1 on CL' and the Fig. 1 tornado plot omits it entirely (only the Bladder and BTC bars are drawn), so unlike tumtyp = 2 and tumtyp = 3 there is no figure label to recover it from. Its composition is INFERRED here as the residual non-NSCLC / non-bladder / non-BTC pool -- small-cell lung cancer (n = 281 per Table 2) plus the miscellaneous advanced solid tumours of Study 1108 -- which is why the register's residual `TUMTP_OTHER` column is used rather than a named tumor-type canonical. The effect is negligible and not distinguishable from zero (Table 4 RSE 185%, bootstrap 95% CI [-0.0228; 0.0466] spanning zero), so no simulation conclusion turns on the inference. The SIGN also conflicts between the paper's own two reports of it: Table 4 prints +0.0101 while the CL_cat.cov equation on page 3 prints (1 - 0.0101). Every other categorical factor in that equation equals (1 + theta) using the Table 4 signed estimate, so tumtyp = 1 is the single term that breaks the pattern; the equation is followed here per the standing text-vs-equation policy. See the vignette Errata.",
      source_name        = "tumtyp"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the covariate analysis but not retained: 'The influence of age and race was not significant' (Results, Durvalumab PopPK). Table 2 reports a median of 63 years (range 19-96). Race was screened on the same footing and is likewise absent from the final model; no race indicator column is defined because the paper reports no coefficient for one."
    ),
    ADA_POS = list(
      description = "Treatment-emergent anti-drug antibody positive status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported and discussed but not retained as a PK covariate: 111 of 3141 patients (3.53%) were treatment-emergent ADA positive across the 6 pooled studies and 'the exposure was comparable for ADA positive and negative patients' (Abstract; Table 2). No ADA coefficient appears in Table 4 or in the CL / Vc equations, so the covariate is documented for provenance only and is deliberately absent from model()."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3141L,
    n_studies      = 6L,
    age_range      = "19-96 years",
    age_median     = "63 years",
    weight_range   = "31.0-175 kg",
    weight_median  = "69.0 kg",
    sex_female_pct = 36.61,
    race_ethnicity = c(White = 67.22, Asian = 26.16, Black = 2.48, Other = 4.19, Multiple = 0.06),
    disease_state  = "Advanced solid tumours (non-small cell lung cancer 60.22%, advanced solid tumours 21.19%, biliary tract cancer 10.02%, small cell lung cancer 8.98%)",
    dose_range     = "3-10 mg/kg every 2 weeks, 15 mg/kg every 3 weeks, 20 mg/kg every 4 weeks, and 1500 mg flat every 3 or 4 weeks, all as intravenous infusions (Table 1; doses < 3 mg/kg were excluded from the population PK analysis because they were not dose-proportional)",
    regions        = "Global",
    treatment_regimen = "Durvalumab monotherapy 60.32%, durvalumab + chemotherapy 28.15% (includes all 314 TOPAZ-1 patients), durvalumab + tremelimumab + chemotherapy 10.29%",
    ecog_distribution = "ECOG 0 37.91%, ECOG 1 61.83%, ECOG 2 0.10%, missing 0.13%",
    notes          = "Baseline demographics per Table 2 of Abegesah 2025 (pooled analysis dataset N = 3141). Studies pooled: CD-ON-MEDI4736-1108 (Study 1108, phase 1/2, advanced solid tumours, n = 1012), D4191C00003 (ATLANTIC, phase 2, NSCLC, n = 443), D4191C00001 (PACIFIC, phase 3, unresectable stage III NSCLC, n = 473), D419QC00001 (CASPIAN, phase 3, extensive-disease SCLC, n = 260), D419MC00004 (POSEIDON, phase 3, metastatic NSCLC, n = 326) and D933AC00001 (TOPAZ-1, phase 3, advanced BTC, n = 314). Other baseline medians: albumin 39.0 g/L (4.10-57.1), creatinine clearance 85.9 mL/min (25.7-363), LDH 242 U/L (18.0-15800), neutrophil-to-lymphocyte ratio 3.46 (0.006-58.8), total bilirubin 0.460 mg/dL (0.006-3.72). Treatment-emergent ADA positive 3.53%. NCI hepatic function normal 85.23% / mild 15.70% / moderate 0.89% / severe 0.03%. The covariate reference values used by the model (weight 69.4 kg, LDH 247 U/L, creatinine clearance 85.66 mL/min) are the 'Previous 5 studies' medians inherited from the POSEIDON model, not the pooled 6-study medians of this analysis."
  )

  ini({
    # ---- Structural PK ---------------------------------------------------
    # Abegesah 2025 Table 4, "Population parameter" block. CL and Q are per
    # DAY (the Table 4 Unit column reads "L/day"); the volumes are in L. The
    # typical-value CL of 0.298 L/day is the t = 0 baseline: the time-varying
    # multiplier exp(cl_tv_mult) equals exp(0) = 1 at t = 0.
    lcl <- log(0.298); label("Baseline clearance CL at the reference covariates (L/day)")   # Table 4: CL = 0.298 L/day (RSE 1.99%, bootstrap 95% CI [0.283; 0.315])
    lvc <- log(3.42); label("Central volume of distribution V1 (L)")                        # Table 4: V1 = 3.42 L (RSE 0.931%, bootstrap 95% CI [3.37; 3.47])
    lq <- log(0.452); label("Intercompartmental clearance Q (L/day)")                       # Table 4: Q = 0.452 L/day (RSE 5.76%, bootstrap 95% CI [0.374; 0.535])
    lvp <- log(1.99); label("Peripheral volume of distribution V2 (L)")                     # Table 4: V2 = 1.99 L (RSE 2.23%, bootstrap 95% CI [1.83; 2.15])

    # ---- Time-varying clearance ------------------------------------------
    # Sigmoidal (Hill) multiplier on log-CL, printed on page 3 as the
    # exp(Tmax * t / (TC50 + t)) factor of the CL_T,i equation. Written below
    # in the general Hill form used by the sibling AstraZeneca model
    # Hwang_2022_tremelimumab.R,
    #   cl_tv_mult = cl_time_max * t^lambda / (cl_t50^lambda + t^lambda),
    # which reduces to the printed equation exactly at lambda = 1. The
    # Table 4 "LAM" row reports 1.00 with no RSE, no bootstrap median and no
    # confidence interval, so it is encoded as fixed().
    #
    # Falsifier: cl_time_max = -0.498 implies an asymptotic clearance ratio of
    # exp(-0.498) = 0.608, i.e. a 39.2% decrease -- matching the Abstract's
    # "the clearance could decrease up to 39% over the time course of
    # treatment" to the printed precision.
    cl_time_max <- -0.498; label("Asymptotic log-change in CL over time (unitless)")        # Table 4: Tmax change CL = -0.498 (RSE 3.90%, bootstrap 95% CI [-0.551; -0.446])
    cl_t50 <- 61.3; label("Time at which the change in CL is 50%% of its asymptote (days)") # Table 4: TC50 change CL = 61.3 days (RSE 8.57%, bootstrap 95% CI [44.5; 84.5])
    cl_time_hill <- fixed(1.00); label("Sigmoidicity exponent of time on CL (unitless)")    # Table 4: LAM = 1.00, reported with no RSE and no confidence interval

    # ---- Covariate effects on CL -----------------------------------------
    # Continuous covariates are power terms normalized to the reference values
    # printed inside the CL_cont.cov equation on page 3 (alb / 39,
    # CrCL / 85.66, LDH / 247, WT / 69.4).
    e_alb_cl <- -0.589; label("Power exponent of baseline albumin on CL (unitless)")        # Table 4: Albumin on CL = -0.589; page 3 CL_cont.cov: (alb_i / 39)^-0.589
    e_crcl_cl <- 0.136; label("Power exponent of creatinine clearance on CL (unitless)")    # Table 4: CrCL on CL = 0.136; page 3 CL_cont.cov: (CrCL_i / 85.66)^0.136
    e_ldh_cl <- 0.0515; label("Power exponent of baseline LDH on CL (unitless)")            # Table 4: LDH on CL = 0.0515; page 3 CL_cont.cov: (LDH_i / 247)^0.0515
    e_wt_cl <- 0.338; label("Power exponent of body weight on CL (unitless)")               # Table 4: Body weight on CL = 0.338; page 3 CL_cont.cov: (WT_i / 69.4)^0.338

    # Categorical covariates enter CL_cat.cov (page 3) as multiplicative
    # factors of the form (1 + theta) using the signed Table 4 estimate, with
    # the reference level contributing a factor of exactly 1.
    e_ecog_cl <- -0.0501; label("Fractional change in CL for ECOG performance status >= 1") # Table 4: ECOG status on CL = -0.0501; page 3 CL_cat.cov: (1 - 0.0501)_{ECOGbin=1}
    e_sexf_cl <- -0.161; label("Fractional change in CL for female sex")                    # Table 4: Sex on CL = -0.161; page 3 CL_cat.cov: (1 - 0.161)_{female}
    e_chemo_cl <- -0.163; label("Fractional change in CL for durvalumab + chemotherapy")    # Table 4: COMB1 on CL = -0.163; page 3 CL_cat.cov: (1 - 0.163)_{comb=1}
    e_treme_cl <- -0.0929; label("Fractional change in CL for durvalumab + tremelimumab + chemotherapy") # Table 4: COMB2 on CL = -0.0929; page 3 CL_cat.cov: (1 - 0.0929)_{comb=2}

    # Tumor-type indicators. The tumtyp = 1 sign follows the page 3 equation,
    # which prints (1 - 0.0101) where Table 4 prints +0.0101; see the
    # covariateData[[TUMTP_OTHER]] note and the vignette Errata.
    e_tumtp_other_cl <- -0.0101; label("Fractional change in CL for the residual 'other' tumor type") # Page 3 CL_cat.cov: (1 - 0.0101)_{tumtyp=1}; Table 4 "Tumor type 1 on CL" prints +0.0101 (RSE 185%, bootstrap 95% CI [-0.0228; 0.0466])
    e_tumtp_bladder_cl <- 0.0698; label("Fractional change in CL for bladder (urothelial) carcinoma")  # Table 4: Tumor type 2 on CL = 0.0698; page 3 CL_cat.cov: (1 + 0.0698)_{tumtyp=2}; Fig. 1 "TUM TYP Bladder" +7%
    e_tumtp_btc_cl <- 0.166; label("Fractional change in CL for biliary tract cancer")                 # Table 4: Tumor type 3 on CL = 0.166; page 3 CL_cat.cov: (1 + 0.166)_{tumtyp=3}; Fig. 1 "TUM TYP BTC" +16.6%

    # ---- Covariate effects on the central volume -------------------------
    # Page 3 Vc equation: Vc,i = 3.42 * (WT_i / 69.4)^0.515 * 1_male
    #                            * (1 - 0.140)_female
    e_wt_vc <- 0.515; label("Power exponent of body weight on Vc (unitless)")               # Table 4: Body weight on V1 = 0.515; page 3 Vc equation: (WT_i / 69.4)^0.515
    e_sexf_vc <- -0.140; label("Fractional change in Vc for female sex")                    # Table 4: Sex on V1 = -0.140; page 3 Vc equation: (1 - 0.140)_{female}

    # ---- Between-subject variability -------------------------------------
    # Table 4 "Interindividual variability" block. The Methods state that
    # "Correlations between CL and V1 was estimated via omega block", and
    # Table 4 reports the off-diagonal directly as "Cov CL-V1" = 0.0390 (a
    # covariance, not a correlation). The implied correlation is
    # 0.0390 / sqrt(0.0795 * 0.0593) = 0.568.
    # Approximate log-scale CVs: sqrt(exp(0.0795) - 1) = 28.8% for CL and
    # sqrt(exp(0.0593) - 1) = 24.7% for Vc.
    etalcl + etalvc ~ c(
      0.0795,
      0.0390, 0.0593
    )                                                                                       # Table 4: ETA CL = 0.0795 (shrinkage 19.1%), Cov CL-V1 = 0.0390, ETA V1 = 0.0593 (shrinkage 25.7%)

    # ADDITIVE (not log-normal) IIV on the time-varying-CL asymptote. The
    # sibling AstraZeneca model Hwang_2022_tremelimumab.R -- same modelling
    # group, same senior author (Zhou D), and the same
    # EMPIR = Tmax * TIME^LAM / (TC50^LAM + TIME^LAM) parameterization --
    # defines Tmax_i = THETA + ETA in its published NONMEM control stream.
    # Abegesah 2025 prints only exp(eta_i) on CL itself and does not state the
    # form for Tmax, so the sibling's verified idiom is followed here; see the
    # vignette Assumptions and deviations.
    etacl_time_max ~ 0.0623                                                                 # Table 4: ETA Tmax = 0.0623 (RSE 9.01%, shrinkage 55.4%)

    # ---- Residual unexplained variability --------------------------------
    # "A combination of proportional and additive residual error model was
    # implemented" (Results). Table 4's Estimate column holds STANDARD
    # DEVIATIONS, not variances: its Bootstrap-median column reports 0.0649
    # and 23.0 for these two rows, which are 0.255^2 = 0.0650 and
    # 4.75^2 = 22.6 -- i.e. that one column slipped to the variance scale
    # while the Estimate and 95% CI columns ([0.246; 0.263] and [3.55; 6.17])
    # stayed on the SD scale. The SD reading is used.
    propSd <- 0.255; label("Proportional residual error (fraction)")                        # Table 4: Proportional component = 0.255 (RSE 0.596%, 95% CI [0.246; 0.263])
    addSd <- 4.75; label("Additive residual error (ug/mL)")                                 # Table 4: Additive component = 4.75 ug/mL (RSE 7.27%, 95% CI [3.55; 6.17])
  })

  model({
    # ---- Combination-therapy stratum indicators --------------------------
    # The paper's three-level 'comb' column is stored as two canonical
    # columns. Tremelimumab is only ever given on a chemotherapy backbone in
    # this pooled dataset, so CONMED_TREMELIMUMAB alone selects comb = 2 and
    # the chemotherapy-without-tremelimumab stratum (comb = 1, which contains
    # all of TOPAZ-1) is the product below. Both zero gives comb = 0,
    # durvalumab monotherapy, the reference level.
    comb_chemo_only <- CONMED_CHEMO * (1 - CONMED_TREMELIMUMAB)

    # ---- Categorical covariate multiplier on CL (page 3, CL_cat.cov) -----
    cl_cat_cov <- (1 + e_chemo_cl)^comb_chemo_only *
      (1 + e_treme_cl)^CONMED_TREMELIMUMAB *
      (1 + e_ecog_cl)^ECOG_GE1 *
      (1 + e_sexf_cl)^SEXF *
      (1 + e_tumtp_other_cl)^TUMTP_OTHER *
      (1 + e_tumtp_bladder_cl)^TUMTP_BLADDER *
      (1 + e_tumtp_btc_cl)^TUMTP_BTC

    # ---- Continuous covariate multiplier on CL (page 3, CL_cont.cov) -----
    cl_cont_cov <- (ALB / 39)^e_alb_cl *
      (CRCL / 85.66)^e_crcl_cl *
      (LDH / 247)^e_ldh_cl *
      (WT / 69.4)^e_wt_cl

    # ---- Time-varying clearance multiplier -------------------------------
    # cl_tv_mult(0) = 0 so cl(0) = cl_base, and
    # cl_tv_mult(t -> Inf) = cl_time_max + etacl_time_max.
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_tv_mult <- cl_time_max_i * t^cl_time_hill /
      (cl_t50^cl_time_hill + t^cl_time_hill)

    # ---- Individual structural parameters (page 3, CL_T,i and Vc,i) ------
    cl_base <- exp(lcl + etalcl) * cl_cat_cov * cl_cont_cov
    cl <- cl_base * exp(cl_tv_mult)

    vc <- exp(lvc + etalvc) * (WT / 69.4)^e_wt_vc * (1 + e_sexf_vc)^SEXF
    q <- exp(lq)
    vp <- exp(lvp)

    # ---- Micro-constants -------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system ------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Observation and combined residual error -------------------------
    # Dose in mg and volumes in L give central / vc in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
