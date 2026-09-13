Zhao_2026_durvalumab <- function() {
  description <- "Two-compartment population PK model for durvalumab (anti-PD-L1 IgG1 kappa) with sigmoidal time-varying clearance in adults with solid tumours, updating the pooled five-study model with the phase III AEGEAN cohort of resectable stage II to IIIB (N2) non-small-cell lung cancer treated with perioperative durvalumab plus platinum-based chemotherapy (Zhao 2026)"
  reference <- paste(
    "Zhao X, Ding J, Zhao J, Zhang L, Abegesah A, Zhang Y-q, O'Brien C,",
    "Doherty GJ, Chen AC, Lim K, Ren S, Ma P, Zhou D. Population",
    "pharmacokinetics and exposure-response analysis of durvalumab in",
    "patients with resectable stage II to IIIB (N2) NSCLC in the phase III",
    "AEGEAN study. Br J Clin Pharmacol. 2026;92(3):980-996.",
    "doi:10.1002/bcp.70287",
    sep = " "
  )
  vignette <- "Zhao_2026_durvalumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Zhao 2026 Methods 2.2 ("Bioanalytical methods") states
  # that "Serum samples for determination of durvalumab concentrations in the
  # AEGEAN study were analysed" by an electrochemiluminescence assay with an
  # LLOQ-ULOQ range of 0.05-3.2 ug/mL on the assayed dilution, so the specimen
  # is serum rather than plasma.
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
      notes              = "Power effects on CL (exponent 0.378) and Vc (exponent 0.503), both normalized to 69.4 kg. The normalizer 69.4 is printed inside the CL_cont.cov and Vc,i equations in Results 3.2; it is the reference weight INHERITED from the previous five-study model, not this analysis's own pooled median of 69.6 kg (Table 1, Total column). The two differ by only 0.2 kg here, but do not conflate them -- the same 69.4 kg normalizer appears in the sibling branch of this lineage, Abegesah_2025_durvalumab.R.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent -0.526) normalized to 39 g/L, printed in the CL_cont.cov equation in Results 3.2. Reported and used in SI g/L (Table 1 median 39.0 g/L, range 3.70-78.0), so no g/dL conversion applies. Albumin is the most influential covariate in the model: Results 3.2 reports a +21.3% change in steady-state CL at the 5th percentile of albumin, the largest single covariate effect in the Figure 1C tornado plot. The same 39 g/L reference is used by the sibling AstraZeneca models Abegesah_2025_durvalumab.R and Hwang_2022_tremelimumab.R.",
      source_name        = "alb"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (raw Cockcroft-Gault style, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (exponent 0.112) normalized to 85.66 mL/min, printed in the CL_cont.cov equation in Results 3.2. Table 1 reports 'Creatinine clearance (mL/min)' with no BSA normalization, so the per-model unit is mL/min rather than the register default mL/min/1.73 m^2 (same precedent as the sibling durvalumab models Abegesah_2025_durvalumab.R, which carries the identical 85.66 normalizer, and deVries_2025_durvalumab.R, whose Baverel 2018 normalizer is 85.65).",
      source_name        = "CrCL"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects for females on CL (1 - 0.166) and Vc (1 - 0.144), printed in the CL_cat.cov and Vc,i equations in Results 3.2 with male as the explicit reference (the '1_male' factor). Direction agrees with every sibling durvalumab model: Abegesah_2025_durvalumab.R (female CL 1 - 0.161, Vc 1 - 0.140), Ogasawara_2020_durvalumab.R (female CL 0.791) and deVries_2025_durvalumab.R (female CL 0.857).",
      source_name        = "sex"
    ),
    ECOG_GE1 = list(
      description        = "Baseline ECOG performance status >= 1 (1 = ECOG 1 or worse, 0 = ECOG 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status 0, 'normal activity')",
      notes              = "Multiplicative effect on CL of (1 - 0.0604), printed in the CL_cat.cov equation in Results 3.2 as the paper's binary 'ECOGbin', whose note defines ECOGbin = 0 as 'normal activity' and ECOGbin = 1 as 'restricted activity or in bed <= 50% of the time'. Table 2 records those three strata as 40.5 / 59.2 / 0.0934 percent, so the ordinal score collapses to this binary indicator. Unlike the sibling Abegesah_2025_durvalumab.R -- where the Fig. 1 tornado plot contradicts the equation's sign -- Zhao 2026 is self-consistent: Table 3 and the printed equation both give -0.0604.",
      source_name        = "ECOGbin"
    ),
    CONMED_CHEMO = list(
      description        = "Durvalumab co-administered with standard-of-care platinum-based chemotherapy (without tremelimumab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (durvalumab monotherapy)",
      notes              = "The paper's 'comb' column is a three-level categorical defined in the Results 3.2 equation note (comb = 0 durvalumab monotherapy, comb = 1 durvalumab + SOC chemotherapy, comb = 2 durvalumab + tremelimumab + SOC chemotherapy) carrying multiplicative CL factors of 1, (1 - 0.0701) and (1 - 0.0578) respectively. It is stored here as two canonical columns, CONMED_CHEMO and CONMED_TREMELIMUMAB, so that each column keeps its register meaning; the model reconstructs the comb = 1 stratum as CONMED_CHEMO * (1 - CONMED_TREMELIMUMAB). IMPORTANT -- this covariate is TIME-VARYING within AEGEAN: the neoadjuvant phase (durvalumab 1500 mg Q3W x 4 with platinum doublet) is comb = 1, whereas the post-surgical adjuvant phase (durvalumab 1500 mg Q4W x 12 as monotherapy) is comb = 0. All published AEGEAN exposure metrics (Tables S3 and S4) are derived from the neoadjuvant phase, i.e. at comb = 1.",
      source_name        = "comb"
    ),
    CONMED_TREMELIMUMAB = list(
      description        = "Durvalumab co-administered with tremelimumab (always on a chemotherapy backbone in this dataset)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no tremelimumab co-administration)",
      notes              = "Identifies the paper's comb = 2 stratum (durvalumab + tremelimumab + SOC chemotherapy, the POSEIDON triplet arm per Table S1), carrying a multiplicative CL factor of (1 - 0.0578) per Table 3 row 'COMB 2 on CL'. Because no tremelimumab-without-chemotherapy patients exist in this pooled dataset, CONMED_TREMELIMUMAB = 1 implies CONMED_CHEMO = 1 and this column alone selects the comb = 2 factor. No AEGEAN patient receives tremelimumab, so this covariate is 0 throughout the validation vignette. This is the mirror image of the sibling model Hwang_2022_tremelimumab.R, where tremelimumab is the analyte and durvalumab co-administration is the covariate.",
      source_name        = "comb"
    )
  )

  covariatesDataExcluded <- list(
    LDH = list(
      description = "Baseline lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "DELIBERATELY REMOVED BY THIS ANALYSIS. LDH on CL was a retained covariate of the previous five-study model (and is still carried by the parallel branch Abegesah_2025_durvalumab.R, exponent 0.0515 normalized to 247 U/L). Results 3.2 states that 'the significance of all covariates in the previous model except for LDH was confirmed' and that 'the final model removed the LDH covariate and retained the other covariates from the previous model'. Accordingly no LDH term appears in Table 3 or in the printed CL_cont.cov equation, and none is implemented here. Table 1 still reports the distribution (Total median 235 U/L; AEGEAN patients had notably lower LDH, median 184 U/L, than previous studies, median 247 U/L) -- which is the stated reason the covariate lost significance when AEGEAN was added."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained. Table 1 reports a median of 63.0 years (range 19.0-96.0) and Table 2 the categories < 65 / 65-75 / >= 75 at 55.4 / 34.5 / 10.0 percent. Results 3.3 states that the impact of age categories on durvalumab exposure was investigated and that 'none of the investigated covariates were considered to have a significant effect on durvalumab exposure'. No age coefficient appears in Table 3."
    ),
    ADA_POS = list(
      description = "Treatment-emergent anti-drug antibody (TEADA) positive status",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported and discussed but explicitly NOT evaluated as a PK covariate: 'The number of TEADA-positive patients was small (3.58%) in the PopPK dataset, and TEADA was not evaluated as a covariate' (Results 3.3). For the 25 TEADA-positive AEGEAN patients the empirical-Bayes exposures differed from TEADA-negative patients by < 20% (Figure S3), so durvalumab PK is unlikely to be affected by ADA status. Table 2 gives 115/3212 (3.58%) positive overall and 25/385 (6.49%) within AEGEAN."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3205L,
    n_studies      = 6L,
    n_observations = 12466L,
    age_range      = "19.0-96.0 years",
    age_median     = "63.0 years",
    weight_range   = "31.0-175 kg",
    weight_median  = "69.6 kg",
    sex_female_pct = 35.0,
    race_ethnicity = c(
      White = 67.3, Asian = 25.1, Black = 2.33, Other = 3.80,
      `American Indian/Alaskan Native` = 1.15,
      `Native Hawaiian or Other Pacific Islander` = 0.249,
      Multiple = 0.0623
    ),
    disease_state  = "Solid tumours (primary indication non-small-cell lung cancer 70.6%, advanced solid tumour 20.6%, small-cell lung cancer 8.75%; by tumour type lung 79.4%, bladder 5.95%, other 14.7%)",
    dose_range     = "Pooled across six studies (Table S1): 0.1-10 mg/kg every 2 weeks, 15 mg/kg every 3 weeks and 20 mg/kg every 4 weeks in Study 1108; 10 mg/kg every 2 weeks in ATLANTIC and PACIFIC; 1500 mg every 3 weeks in CASPIAN; 1500 mg every 3 weeks (+/- tremelimumab 75 mg) then 1500 mg every 4 weeks in POSEIDON; and in AEGEAN 1500 mg every 3 weeks for 4 neoadjuvant cycles with platinum-based chemotherapy followed by 1500 mg every 4 weeks for 12 adjuvant cycles. All intravenous. Doses < 3 mg/kg were excluded from the population PK analysis for lack of dose proportionality (Table S1 footnote a).",
    regions        = "Global (Europe 38.5%, North America 33.5%, Asia 22.6%, South America 3.61%, Africa 0.778%, other 0.996%)",
    treatment_regimen = "Durvalumab monotherapy, durvalumab + standard-of-care chemotherapy, and durvalumab + tremelimumab + standard-of-care chemotherapy (the paper's comb = 0, 1 and 2 strata). AEGEAN contributes comb = 1 during the neoadjuvant phase and comb = 0 during the adjuvant phase.",
    ecog_distribution = "Normal activity 40.5%, restricted activity 59.2%, in bed <= 50% of the time 0.0934%, missing 0.125%",
    renal_function = "By baseline creatinine clearance (Table 2 footnote b): normal (>= 90 mL/min) 43.0%, mild (60-<90) 40.6%, moderate (30-<60) 14.4%, severe (< 30) 1.96%",
    hepatic_function = "NCI-ODWG: normal 86.4%, mild 11.8%, moderate 0.436%, severe 0.0311%, missing 1.34%",
    notes          = "Baseline demographics per Tables 1 and 2 of Zhao 2026, which summarise N = 3212. The final analysis dataset is 12 466 PK samples from 3205 patients: 2827 evaluable patients from the five previous studies plus 385 from AEGEAN gives 3212, from which 7 were excluded for physiologically impossible covariate values in previous analyses, and 145 below-LLOQ samples (1.16%) were dropped (Results 3.1). Studies pooled (Table S1): CD-ON-MEDI4736-1108 (Study 1108, phase 1/2, advanced solid tumours, n = 1012), D4191C00003 (ATLANTIC, phase 2, locally advanced/metastatic NSCLC, n = 443), D4191C00001 (PACIFIC, phase 3, unresectable stage III NSCLC, n = 473), D419QC00001 (CASPIAN, phase 3, extensive-stage SCLC, n = 260), D419MC00004 (POSEIDON, phase 3, metastatic NSCLC, n = 326) and D9106C00001 (AEGEAN, phase 3, resectable stage II to IIIB (N2) NSCLC, n = 385 evaluable). Other pooled baseline medians: creatinine clearance 85.5 mL/min (25.7-718), albumin 39.0 g/L (3.70-78.0), LDH 235 U/L (4.00-15 800), AST 20.0 IU/L, ALT 18.0 IU/L, total bilirubin 0.420 mg/dL, neutrophil-to-lymphocyte ratio 3.30 (63.6% missing). Treatment-emergent ADA positive 3.58%. The covariate reference values used inside the model equations (weight 69.4 kg, creatinine clearance 85.66 mL/min, albumin 39 g/L) are inherited from the previous five-study model, not this analysis's own pooled medians (69.6 kg, 85.5 mL/min, 39.0 g/L)."
  )

  ini({
    # ---- Structural PK ---------------------------------------------------
    # Zhao 2026 Table 3, "Population parameter" block. CL and Q are per DAY
    # (the Table 3 Unit column reads "L/day"); the volumes are in L. The
    # typical-value CL of 0.285 L/day is the t = 0 baseline: the time-varying
    # multiplier exp(cl_tv_mult) equals exp(0) = 1 at t = 0. Discussion
    # corroborates: "the typical clearance and V1 were 0.285 L/day and 3.42 L".
    lcl <- log(0.285); label("Baseline clearance CL at the reference covariates (L/day)")   # Table 3: CL = 0.285 L/day (RSE 1.68%, bootstrap median 0.285, 95% CI [0.270; 0.303])
    lvc <- log(3.42); label("Central volume of distribution V1 (L)")                        # Table 3: V1 = 3.42 L (RSE 0.962%, bootstrap median 3.42, 95% CI [3.38; 3.48])
    lq <- log(0.381); label("Intercompartmental clearance Q (L/day)")                       # Table 3: Q intercompartmental = 0.381 L/day (RSE 4.90%, bootstrap median 0.381, 95% CI [0.297; 0.465])
    lvp <- log(2.30); label("Peripheral volume of distribution V2 (L)")                     # Table 3: V2 = 2.30 L (RSE 2.09%, bootstrap median 2.30, 95% CI [2.15; 2.44])

    # ---- Time-varying clearance ------------------------------------------
    # Sigmoidal (Hill) multiplier on log-CL. Results 3.2 prints the CL_T,i
    # equation as
    #   CL_T,i = 0.285 * CL_cat.cov * CL_cont.cov * exp(-0.412 * t / (48 + t))
    #            * exp(eta_i)
    # Written below in the general Hill form used by the sibling AstraZeneca
    # models Abegesah_2025_durvalumab.R and Hwang_2022_tremelimumab.R,
    #   cl_tv_mult = cl_time_max * t^lambda / (cl_t50^lambda + t^lambda),
    # which reduces to the printed equation exactly at lambda = 1. The Table 3
    # "LAM change CL" row reports 1.00 with no RSE, no bootstrap median and no
    # confidence interval, so it is encoded as fixed().
    #
    # Falsifier: cl_time_max = -0.412 implies an asymptotic clearance ratio of
    # exp(-0.412) = 0.662, i.e. a 33.8% decrease -- matching the Discussion's
    # "clearance could decrease by a maximum of 34%" to the printed precision.
    cl_time_max <- -0.412; label("Asymptotic log-change in CL over time (unitless)")        # Table 3: Tmax change CL = -0.412 (RSE 4.91%, bootstrap median -0.414, 95% CI [-0.465; -0.360]); the Table 3 Unit column prints 'L/day' for this row, which is a misprint -- it is a unitless log-scale change (see vignette Errata)
    cl_t50 <- 48.0; label("Time at which the change in CL is 50%% of its asymptote (days)") # Table 3: TC50 change CL = 48.0 days (RSE 10.9%, bootstrap median 47.4, 95% CI [32.8; 71.9])
    cl_time_hill <- fixed(1.00); label("Sigmoidicity exponent of time on CL (unitless)")    # Table 3: LAM change CL = 1.00, reported with no RSE, no bootstrap median and no confidence interval

    # ---- Covariate effects on CL -----------------------------------------
    # Continuous covariates are power terms normalized to the reference values
    # printed inside the CL_cont.cov equation in Results 3.2:
    #   CL_cont.cov = (alb_i / 39)^-0.526 * (CrCL_i / 85.66)^0.112
    #                 * (WT_i / 69.4)^0.378
    e_alb_cl <- -0.526; label("Power exponent of baseline albumin on CL (unitless)")        # Table 3: Albumin on CL = -0.526 (RSE 2.88%, bootstrap median -0.537, 95% CI [-0.776; -0.389]); Results 3.2 CL_cont.cov: (alb_i / 39)^-0.526
    e_crcl_cl <- 0.112; label("Power exponent of creatinine clearance on CL (unitless)")    # Table 3: Creatinine clearance on CL = 0.112 (RSE 19.1%, bootstrap median 0.111, 95% CI [0.0699; 0.153]); Results 3.2 CL_cont.cov: (CrCL_i / 85.66)^0.112
    e_wt_cl <- 0.378; label("Power exponent of body weight on CL (unitless)")               # Table 3: Bodyweight on CL = 0.378 (RSE 9.31%, bootstrap median 0.382, 95% CI [0.311; 0.453]); Results 3.2 CL_cont.cov: (WT_i / 69.4)^0.378

    # Categorical covariates enter CL_cat.cov (Results 3.2) as multiplicative
    # factors of the form (1 + theta) using the signed Table 3 estimate, with
    # the reference level contributing a factor of exactly 1:
    #   CL_cat.cov = 1_{comb=0} * (1 - 0.0701)_{comb=1} * (1 - 0.0578)_{comb=2}
    #                * 1_{ECOGbin=0} * (1 - 0.0604)_{ECOGbin=1}
    #                * 1_male * (1 - 0.166)_female
    e_ecog_cl <- -0.0604; label("Fractional change in CL for ECOG performance status >= 1 (fraction)")               # Table 3: ECOG status on CL = -0.0604 (RSE 19.1%, bootstrap median -0.0590, 95% CI [-0.0836; -0.0325]); Results 3.2 CL_cat.cov: (1 - 0.0604)_{ECOGbin=1}
    e_sexf_cl <- -0.166; label("Fractional change in CL for female sex (fraction)")                                  # Table 3: Sex on CL = -0.166 (RSE 7.40%, bootstrap median -0.165, 95% CI [-0.189; -0.139]); Results 3.2 CL_cat.cov: (1 - 0.166)_female
    e_chemo_cl <- -0.0701; label("Fractional change in CL for durvalumab + SOC chemotherapy (fraction)")             # Table 3: COMB1 on CL = -0.0701 (RSE 17.9%, bootstrap median -0.0694, 95% CI [-0.0984; -0.0406]); Results 3.2 CL_cat.cov: (1 - 0.0701)_{comb=1}
    e_treme_cl <- -0.0578; label("Fractional change in CL for durvalumab + tremelimumab + SOC chemotherapy (fraction)") # Table 3: COMB 2 on CL = -0.0578 (RSE 29.0%, bootstrap median -0.0561, 95% CI [-0.108; -0.00867]); Results 3.2 CL_cat.cov: (1 - 0.0578)_{comb=2}

    # ---- Covariate effects on the central volume -------------------------
    # Results 3.2 Vc equation:
    #   Vc,i = 3.42 * (WT_i / 69.4)^0.503 * 1_male * (1 - 0.144)_female
    e_wt_vc <- 0.503; label("Power exponent of body weight on Vc (unitless)")               # Table 3: Bodyweight on V1 = 0.503 (RSE 5.73%, bootstrap median 0.502, 95% CI [0.446; 0.558]); Results 3.2 Vc,i: (WT_i / 69.4)^0.503
    e_sexf_vc <- -0.144; label("Fractional change in Vc for female sex (fraction)")         # Table 3: Sex on V1 = -0.144 (RSE 8.27%, bootstrap median -0.143, 95% CI [-0.166; -0.121]); Results 3.2 Vc,i: (1 - 0.144)_female

    # ---- Between-subject variability -------------------------------------
    # Table 3 "Interindividual variability" block. The off-diagonal is reported
    # directly as "Cov CL-V1" = 0.0408 (a covariance, not a correlation); the
    # implied correlation is 0.0408 / sqrt(0.0845 * 0.0563) = 0.592.
    # Approximate log-scale CVs: sqrt(exp(0.0845) - 1) = 29.7% for CL and
    # sqrt(exp(0.0563) - 1) = 24.1% for Vc.
    etalcl + etalvc ~ c(
      0.0845,
      0.0408, 0.0563
    )                                                                                       # Table 3: ETA CL = 0.0845 (RSE 3.43%, shrinkage 19.4%), Cov CL-V1 = 0.0408 (RSE 5.77%), ETA V1 = 0.0563 (RSE 3.74%, shrinkage 28.9%)

    # ADDITIVE (not log-normal) IIV on the time-varying-CL asymptote. The
    # sibling AstraZeneca model Hwang_2022_tremelimumab.R -- same modelling
    # group, same senior author (Zhou D), and the same
    # EMPIR = Tmax * TIME^LAM / (TC50^LAM + TIME^LAM) parameterization --
    # defines Tmax_i = THETA + ETA in its published NONMEM control stream, and
    # the sibling durvalumab model Abegesah_2025_durvalumab.R follows it. Zhao
    # 2026 prints only exp(eta_i) on CL itself and does not state the form for
    # Tmax, so that verified idiom is followed here; see the vignette
    # Assumptions and deviations. The 56.5% shrinkage on this eta means the
    # published estimate is weakly informed by the data.
    etacl_time_max ~ 0.0534                                                                 # Table 3: ETA Tmax = 0.0534 (RSE 9.53%, bootstrap median 0.0537, 95% CI [0.0341; 0.0772], shrinkage 56.5%)

    # ---- Residual unexplained variability --------------------------------
    # Table 3 "Residual variability" block, reported on the STANDARD DEVIATION
    # scale: the Estimate, Bootstrap-median and 95% CI columns agree with each
    # other for both rows (0.253 / 0.253 / [0.245; 0.261] and
    # 5.38 / 5.35 / [4.16; 6.80]). Note that the sibling
    # Abegesah_2025_durvalumab.R has a bootstrap-median column that slipped to
    # the variance scale for these two rows; Zhao 2026 does NOT, so there is no
    # scale ambiguity to resolve here. Epsilon-shrinkage was 13.5%.
    propSd <- 0.253; label("Proportional residual error (fraction)")                        # Table 3: Proportional component = 0.253 (RSE 0.627%, bootstrap median 0.253, 95% CI [0.245; 0.261], shrinkage 13.5%)
    addSd <- 5.38; label("Additive residual error (ug/mL)")                                 # Table 3: Additive component = 5.38 ug/mL (RSE 6.38%, bootstrap median 5.35, 95% CI [4.16; 6.80], shrinkage 13.5%)
  })

  model({
    # ---- Combination-therapy stratum indicators --------------------------
    # The paper's three-level 'comb' column is stored as two canonical columns.
    # Tremelimumab is only ever given on a chemotherapy backbone in this pooled
    # dataset, so CONMED_TREMELIMUMAB alone selects comb = 2 and the
    # chemotherapy-without-tremelimumab stratum (comb = 1, which contains the
    # AEGEAN neoadjuvant phase) is the product below. Both zero gives comb = 0,
    # durvalumab monotherapy, the reference level.
    comb_chemo_only <- CONMED_CHEMO * (1 - CONMED_TREMELIMUMAB)

    # ---- Categorical covariate multiplier on CL (Results 3.2, CL_cat.cov) --
    cl_cat_cov <- (1 + e_chemo_cl)^comb_chemo_only *
      (1 + e_treme_cl)^CONMED_TREMELIMUMAB *
      (1 + e_ecog_cl)^ECOG_GE1 *
      (1 + e_sexf_cl)^SEXF

    # ---- Continuous covariate multiplier on CL (Results 3.2, CL_cont.cov) --
    cl_cont_cov <- (ALB / 39)^e_alb_cl *
      (CRCL / 85.66)^e_crcl_cl *
      (WT / 69.4)^e_wt_cl

    # ---- Time-varying clearance multiplier -------------------------------
    # cl_tv_mult(0) = 0 so cl(0) = cl_base, and
    # cl_tv_mult(t -> Inf) = cl_time_max + etacl_time_max.
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_tv_mult <- cl_time_max_i * t^cl_time_hill /
      (cl_t50^cl_time_hill + t^cl_time_hill)

    # ---- Individual structural parameters (Results 3.2, CL_T,i and Vc,i) ---
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
