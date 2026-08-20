Berkhout_2015_osteoporosis_placebo_qsp <- function() {
  description <- paste0(
    "QSP. Systems-pharmacology disease-progression model of the placebo ",
    "(calcium 500 mg/day) arm of postmenopausal osteoporosis. Applies the ",
    "reduced Lemaire bone-remodeling core (Schmidt 2011 / Post 2013) to ",
    "470 women from the placebo arm of the EPIC study, evaluated over 4 ",
    "years with baseline years-since-menopause (YSM) spanning 0.5-27 yr. ",
    "The mechanistic core has two dimensionless ODE states (relative ",
    "active osteoblast y, relative active osteoclast z, both = 1 at ",
    "menopause onset) driven by (i) an exponential estrogen-loss disease ",
    "progression function f(t) = exp(-k_estrogen * t) and (ii) a placebo ",
    "(calcium) inhibition factor PCa acting from t_start (subject-",
    "specific YSM at study entry) as a delayed-onset / slow-offset dip: ",
    "PCa = 1 - PLAC * (1 - exp(-k_Ca_onset * (t - t_start))) * exp(",
    "-k_Ca_offset * (t - t_start)). The core drives four biomarkers: ",
    "algebraic transducers for urine N-telopeptide (NTX = NTX0 * z^q_NTX) ",
    "and serum bone-specific alkaline phosphatase (BSAP = BSAP0 * (1 + ",
    "k_BSAP0) * y^q_BSAP; k_BSAP0 rescales the fixed tibolone-study U/L ",
    "baseline to the EPIC ng/mL assay), and two indirect-response ODEs ",
    "for lumbar-spine and total-hip BMD driven by osteoblast (production ",
    "coefficient D_AOB) and osteoclast (degradation coefficient D_AOC). ",
    "BMD baselines are BMI-corrected around the cohort median BMI of ",
    "25.34114 kg/m^2. IIV is a 2 x 2 block on (NTX0, BSAP0) plus a 2 x 2 ",
    "block on (BMD_LS_0, BMD_TH_0); residual error is log-normal on every ",
    "biomarker. The 5 system parameters (z_s, k_B, k_estrogen, D_A, b) ",
    "are fixed to the tibolone-study (Post 2013) values."
  )
  reference <- paste(
    "Berkhout J, Stone JA, Verhamme KM, Stricker BH, Sturkenboom MC,",
    "Danhof M, Post TM.",
    "Application of a systems pharmacology-based placebo population",
    "model to analyze long-term data of postmenopausal osteoporosis.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(9):516-526.",
    "doi:10.1002/psp4.12006.",
    "Structural model inherited from Post et al. J Pharmacokinet",
    "Pharmacodyn. 2013;40(2):143-156, which reduced the Lemaire et al.",
    "J Theor Biol 2004 cell-interaction model per Schmidt et al.",
    "J Pharmacokinet Pharmacodyn. 2011;38(6):873-900.",
    sep = " "
  )
  vignette <- "Berkhout_2015_osteoporosis_placebo"
  units <- list(
    time          = "day (days since onset of menopause; t = 0 at menopause onset)",
    dosing        = "n/a (placebo / disease-progression model; no drug input event)",
    concentration = paste0(
      "NTX (nmol bone collagen equivalents (bce) / mmol creatinine); ",
      "BSAP (ng/mL); BMD_LS and BMD_TH (g/cm^2); ",
      "osteoblast and osteoclast (dimensionless, relative to baseline = 1)"
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    osteoblast = list(analyte = "active osteoblast", units = NA_character_, specimen = "tissue", verified = FALSE),
    osteoclast = list(analyte = "active osteoclast", units = NA_character_, specimen = "tissue", verified = FALSE),
    BMD_LS     = list(analyte = "bone mineral density (lumbar spine)", units = NA_character_, specimen = "tissue", verified = FALSE),
    BMD_TH     = list(analyte = "bone mineral density (total hip)", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    BMI = list(
      description        = paste0(
        "Body mass index at study baseline. Enters the BMD baseline ",
        "equations as a fractional shift from the cohort median: ",
        "BMD_LS(0) = BMD_LS_0 * (1 + BMI_frac_LS * (BMI - 25.34114)) ",
        "and BMD_TH(0) = BMD_TH_0 * (1 + BMI_frac_TH * (BMI - 25.34114)). ",
        "Time-fixed per subject in the source EPIC-arm fit."
      ),
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = 25.34114,
      notes              = paste0(
        "Centred (referenced) at 25.34114 kg/m^2 in the source NONMEM code ",
        "(supplement $PK block: 'BMI-25.34114'). The paper text (Methods, ",
        "'Body composition is known to induce changes in bone morphology') ",
        "reports the median EPIC 2 cohort BMI as 25.2 kg/m^2, so the coded ",
        "25.34114 is a small refinement over the printed median; the ",
        "coded value is authoritative for the ini() coefficients ",
        "BMI_frac_LS = 0.0111 and BMI_frac_TH = 0.0154 (Table 2)."
      ),
      source_name        = "BMI"
    ),
    T_ENTRY = list(
      description        = paste0(
        "Subject-specific time (in days since onset of menopause) at ",
        "which the study (placebo / 500 mg calcium daily) treatment ",
        "starts -- i.e., the subject's YSM at baseline expressed on the ",
        "model integration axis (days). Equals YSM_baseline * 365.25 in ",
        "typical use. Enters the PCa placebo function so that PCa = 1 ",
        "for t < T_ENTRY (subject on pre-study disease trajectory) and ",
        "PCa = 1 - (1 - exp(-k_Ca_onset * (t - T_ENTRY))) * exp(",
        "-k_Ca_offset * (t - T_ENTRY)) for t >= T_ENTRY."
      ),
      units              = "day (days since menopause onset)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed per subject. In the source NONMEM data set this is ",
        "the column 'STDA' (start-day-of-treatment), which the $PK block ",
        "aliases as STAR. Subject YSM at baseline in the EPIC study ",
        "spans 0.5-27 years (Methods, Table 1) so T_ENTRY ranges roughly ",
        "180-9860 days. The complementary column 'ENDA' (end-day-of-",
        "treatment, SEND in $PK) exists in the source data but is not ",
        "used to gate PCa in the published control stream -- PCa remains ",
        "on for all t >= T_ENTRY -- so no T_END covariate is encoded; if ",
        "extending the simulation beyond the 4-year observation window, ",
        "the user should decide whether the placebo function should ",
        "continue or be switched off manually. Uses the canonical ",
        "T_ENTRY register entry (per-subject study-entry time on the ",
        "model integration axis; founding example Delor 2013)."
      ),
      source_name        = "STDA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 470L,
    n_studies       = 1L,
    age_range       = "45-59 years at baseline (inclusion criterion; Methods, Subject population)",
    age_median      = "53.3 +/- 3.7 years (EPIC 2, Table 1)",
    weight_range    = "not reported (Table 1 records BMI, not weight)",
    weight_median   = NA_character_,
    sex_female_pct  = 100,
    race_ethnicity  = paste0(
      "Not reported per subject in Table 1. Four study centers (two in ",
      "the United States and two in Europe); the source NONMEM $INPUT ",
      "declares a RACE column but the paper does not tabulate the racial ",
      "composition or fit a race effect."
    ),
    disease_state   = paste0(
      "Postmenopausal women at risk of osteoporosis (EPIC study inclusion ",
      "criteria: at least 6 months past menopause at baseline, in good ",
      "general health, no laboratory evidence of confounding systemic ",
      "disease). Only 10% of enrollees per center were allowed to have ",
      "an LS-BMD below 0.8 g/cm^2 at baseline, so the cohort is largely ",
      "non-osteoporotic at study entry."
    ),
    dose_range      = paste0(
      "Placebo tablet + at least 500 mg calcium per day (baseline dietary ",
      "calcium plus supplements as needed). All subjects were told to ",
      "achieve at least 500 mg/day; monitored by food-frequency ",
      "questionnaire at baseline and annually (Methods)."
    ),
    regions         = "United States (2 centers) and Europe (2 centers)",
    ysm_range       = "0.5-27 years since menopause at baseline (mean 5.7 +/- 5.4 yr; EPIC 2, Table 1)",
    bmi_range       = "mean 25.2 +/- 3.6 kg/m^2 (EPIC 2, Table 1)",
    ls_bmd_baseline = "0.94 +/- 0.12 g/cm^2 (EPIC 2, Table 1)",
    th_bmd_baseline = "0.85 +/- 0.12 g/cm^2 (EPIC 2, Table 1)",
    ntx_baseline    = "88.0 +/- 45.0 nmol bce/mmol cr (EPIC 2, Table 1)",
    bsap_baseline   = "11.1 +/- 4.4 ng/mL (EPIC 2, Table 1)",
    notes           = paste0(
      "EPIC (Early Postmenopausal Intervention Cohort) placebo arm. ",
      "Randomized double-blind alendronate-vs-placebo study; only the ",
      "n = 470 placebo arm is used here. Observations: LS-BMD and TH-BMD ",
      "by dual-energy X-ray absorptiometry (Hologic model 2000) twice at ",
      "baseline and annually for 4 years; NTX (urine, Osteomark) and ",
      "BSAP (serum, Ostase; n = 205 sub-sample) every 6 months. Note ",
      "that BSAP in the tibolone study (Gallagher 2001) used a different ",
      "assay reporting in U/L with baseline ~106.7 U/L; the k_BSAP0 ",
      "scaling parameter (Table 2, EPIC 2 = -0.896) converts the fixed ",
      "tibolone BSAP0 = 97.4 U/L to the EPIC ng/mL scale: BSAP0_ngml = ",
      "97.4 * (1 + (-0.896)) = 97.4 * 0.104 ~ 10.13 ng/mL, close to the ",
      "observed 11.1 ng/mL baseline. Model was fit in NONMEM 7.3.0 with ",
      "FOCE-I on log-transformed observations (LTBS) -- so the ",
      "'proportional' residual errors in Table 2 are natural-log-scale ",
      "SDs, encoded here via 'lnorm(sd)'."
    ),
    scope_note      = paste0(
      "The paper reports two model fits: EPIC 1 (subset of 222 women ",
      "with 1 <= YSM <= 5 yr, using the original tibolone zeroth-order ",
      "BMD equation) as a QUALIFICATION check, and EPIC 2 (all 470 women ",
      "with the indirect-response BMD equation) as the FINAL model. This ",
      "file encodes only the EPIC 2 final model (Table 2 column ",
      "'Value (%CV) EPIC 2 All YSM'). The paper's fully-mechanistic ",
      "predecessor is the Lemaire et al. 2004 cell-interaction model; ",
      "the reduced form used here follows Schmidt et al. 2011 and Post ",
      "et al. 2013."
    )
  )

  ini({
    # =====================================================================
    # SYSTEM-RELATED PARAMETERS (all fixed to the tibolone-study estimates;
    # Berkhout 2015 Table 2, column 'Value (%CV) tibolone study' with the
    # word 'fixed' in the EPIC columns, and NONMEM $THETA lines TH1-TH5).
    # =====================================================================

    zs        <- fixed(0.659);   label("Constant z_s in the osteoclast (z) function pi_z = z/(z+z_s) (dimensionless)")
    # Table 2, tibolone: 'fixed at 0.659'; NONMEM $THETA: '(0.659) FIX'.

    k_B       <- fixed(0.0109);  label("Osteoblast (y) apoptosis / turnover rate k_B (1/day)")
    # Table 2, tibolone: 'fixed at 0.0109'; NONMEM $THETA: '(0.0109) FIX'.

    k_estrogen <- fixed(0.00763); label("Estrogen-loss disease-progression rate k_estrogen (1/day); f(t) = exp(-k_estrogen*t)")
    # Table 2, tibolone: 'fixed at 0.00763'; NONMEM $THETA: '(0.00763) FIX'.

    D_A       <- fixed(1);       label("Osteoclast (z) apoptosis rate factor D_A (1/day)")
    # Table 2, tibolone: 'Could not be estimated, fixed at 1'; NONMEM $THETA: '(1) FIX'.

    b_baseline <- fixed(1);      label("Disease process baseline b (dimensionless; beta_0 = 1 at menopause onset)")
    # Table 2, tibolone: 'Could not be estimated, fixed at 1'; NONMEM $THETA: '(1) FIX'.

    # =====================================================================
    # PLACEBO-EFFECT (CALCIUM) PARAMETERS (Berkhout 2015 Table 2, EPIC 2).
    # =====================================================================

    k_Ca_onset  <- 0.0009;   label("Calcium (placebo) elimination-rate onset constant k_Ca_onset (1/day)")
    # Table 2, EPIC 2 k_Ca,onset = 0.0009 (%CV 12.4).

    k_Ca_offset <- 0.000226; label("Calcium (placebo) elimination-rate offset constant k_Ca_offset (1/day)")
    # Table 2, EPIC 2 K_Ca,offset = 0.000226 (%CV 21.9).

    # =====================================================================
    # BIOMARKER-TRANSDUCER PARAMETERS (Berkhout 2015 Table 2, EPIC 2).
    # =====================================================================

    NTX_0    <- 49.5;         label("Baseline urine NTX at menopause onset NTX_0 (nmol bce/mmol cr)")
    # Table 2, EPIC 2 NTX_0 = 49.5 (%CV 5.7).

    BSAP_0   <- fixed(97.4);  label("Baseline serum BSAP in tibolone-study units BSAP_0 (U/L; Post 2013 value)")
    # Table 2, all columns: 'fixed at 97.4' (footnote a: 'Fixed at the tibolone value, see main text for explanation').
    # NONMEM $THETA: '(97.4) FIX'.

    k_BSAP0  <- -0.896;       label("BSAP baseline unit-scaling parameter k_BSAP0 (dimensionless; effective BSAP baseline = BSAP_0 * (1 + k_BSAP0))")
    # Table 2, EPIC 2 k_BSAP0 = -0.896 (%CV 0.3). Rescales BSAP_0 = 97.4 U/L to the EPIC ng/mL assay: 97.4 * (1 - 0.896) ~ 10.13 ng/mL.

    q_NTX    <- 0.56;         label("NTX transducer exponent q_NTX (dimensionless; NTX = NTX_0 * z^q_NTX)")
    # Table 2, EPIC 2 q_NTX = 0.56 (%CV 15.6).

    q_BSAP   <- fixed(0.286); label("BSAP transducer exponent q_BSAP (dimensionless; BSAP = (BSAP_0*(1+k_BSAP0)) * y^q_BSAP; tibolone value)")
    # Table 2, EPIC columns: 'fixed at tibolone value' (footnote a). NONMEM $THETA: '(0, 0.286, 10) FIX'.

    # =====================================================================
    # BMD (INDIRECT-RESPONSE) PARAMETERS (Berkhout 2015 Table 2, EPIC 2).
    # k_out is derived in model() as k_in / BMD(0) per the source NONMEM
    # $PK block; only k_in is estimated / reported.
    # =====================================================================

    BMD_LS_0  <- 0.99;   label("Baseline lumbar-spine BMD at BMI = 25.34114 kg/m^2 BMD_LS_0 (g/cm^2)")
    # Table 2, EPIC 2 LS-BMD_0 = 0.99 (%CV 0.7).

    BMD_TH_0  <- 0.88;   label("Baseline total-hip BMD at BMI = 25.34114 kg/m^2 BMD_TH_0 (g/cm^2)")
    # Table 2, EPIC 2 TH-BMD_0 = 0.88 (%CV 0.7).

    k_in_LS   <- 1.13;   label("LS-BMD zero-order production rate k_in_LS (source-reported unit: mg/day; see vignette Errata for the unit discussion)")
    # Table 2, EPIC 2 k_in,ls = 1.13 (%CV 22.7).

    k_in_TH   <- 0.295;  label("TH-BMD zero-order production rate k_in_TH (source-reported unit: mg/day; see vignette Errata for the unit discussion)")
    # Table 2, EPIC 2 k_in,th = 0.295 (%CV 14.7).

    BMI_frac_LS <- 0.0111; label("BMI fractional-shift coefficient on LS-BMD baseline BMI_frac_LS (1/(kg/m^2); BMD_LS(0) = BMD_LS_0 * (1 + BMI_frac_LS * (BMI - 25.34114)))")
    # Table 2, EPIC 2 BMI-LS-BMD_0 fraction = 0.0111 (%CV 14.5).

    BMI_frac_TH <- 0.0154; label("BMI fractional-shift coefficient on TH-BMD baseline BMI_frac_TH (1/(kg/m^2); BMD_TH(0) = BMD_TH_0 * (1 + BMI_frac_TH * (BMI - 25.34114)))")
    # Table 2, EPIC 2 BMI-TH-BMD_0 fraction = 0.0154 (%CV 11.4).

    D_AOB <- -0.121;   label("Osteoblast stimulation coefficient on BMD production D_AOB (dimensionless)")
    # Table 2, EPIC 2 D_AOB = 0.121 (%CV 6.0) printed WITHOUT a sign. Encoded here as NEGATIVE
    # (-0.121) based on three independent lines of evidence: (i) the NONMEM $THETA initial value
    # in the supplementary control stream is -0.124 with bounds (-1, 1), (ii) Table 2 has a
    # demonstrated sign-drop typesetting error on the kBSAP0 row where the PDF renders "20.896"
    # for the true value -0.896 (matching NONMEM initial -0.894), and (iii) only NEGATIVE
    # values reproduce the postmenopausal LS-BMD decline shown in paper Figure 3 (VPC);
    # POSITIVE 0.121 / 0.0456 give BMD INCREASE over the 4-year study window (unphysical).
    # See vignette Errata for the full sign-inference discussion.

    D_AOC <- -0.0456;  label("Osteoclast stimulation coefficient on BMD degradation D_AOC (dimensionless)")
    # Table 2, EPIC 2 D_AOC = 0.0456 (%CV 10.6) printed WITHOUT a sign. Encoded here as NEGATIVE
    # (-0.0456) matching NONMEM $THETA initial -0.0467 and the physiological requirement above.
    # See vignette Errata for the same discussion as D_AOB.

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY (Berkhout 2015 Table 2, EPIC 2).
    # Table 2 reports IIV as %CV on the exponential model P_i = P * exp(eta).
    # Variance on the natural-log scale is Var[eta] = log(1 + CV^2).
    # NONMEM $OMEGA BLOCK(2) initial values match to 3 decimals.
    # =====================================================================

    # 2x2 block on (NTX_0, BSAP_0). Table 2 EPIC 2: IIV NTX_0 = 40%,
    # IIV BSAP_0 = 32%, correlation = 0.50. Variances:
    #   var(eta_NTX0)   = log(1 + 0.40^2) = log(1.16)   = 0.14842
    #   var(eta_BSAP0)  = log(1 + 0.32^2) = log(1.1024) = 0.09748
    #   cov             = 0.50 * sqrt(0.14842 * 0.09748) = 0.06015
    etaNTX_0 + etaBSAP_0 ~ c(
      0.14842,
      0.06015, 0.09748
    )
    # Table 2, EPIC 2: IIV NTX_0 = 40% (%CV 3.5); IIV BSAP_0 = 32% (%CV 5.6); IIV corr = 0.50 (%CV 12.1).

    # 2x2 block on (BMD_LS_0, BMD_TH_0). Table 2 EPIC 2: IIV LS-BMD_0 = 12%,
    # IIV TH-BMD_0 = 12%, correlation = 0.60. Variances:
    #   var(eta_LS)   = log(1 + 0.12^2) = log(1.0144) = 0.01430
    #   var(eta_TH)   = log(1 + 0.12^2) = 0.01430
    #   cov           = 0.60 * sqrt(0.01430 * 0.01430) = 0.008579
    etaBMD_LS_0 + etaBMD_TH_0 ~ c(
      0.01430,
      0.008579, 0.01430
    )
    # Table 2, EPIC 2: IIV LS-BMD_0 = 12% (%CV 3.4); IIV TH-BMD_0 = 12% (%CV 3.4); IIV corr = 0.60 (%CV 1.4).

    # =====================================================================
    # RESIDUAL VARIABILITY (Berkhout 2015 Table 2, EPIC 2).
    # Model was fit with LTBS (Supplementary Information Eq. 8:
    #   ln(y_obs) = ln(y_pred) + eps
    # so the values below are natural-log-scale SDs, encoded via lnorm().
    # A second residual-error term for extreme values (SW4 / SW5 in the
    # $ERROR block) that switches on when a NTX or BSAP log observation
    # falls outside the 1%-99% quantile is documented in vignette Errata
    # but not encoded here (the extreme mask depends on observation-level
    # switches SW4 / SW5 that nlmixr2's residual model does not accept).
    # =====================================================================

    expSd_NTX    <- 0.314; label("Log-normal residual SD for NTX (natural-log scale)")
    # Table 2, EPIC 2 eps_NTX (SD) = 0.314 (%CV 2.0).

    expSd_BSAP   <- 0.184; label("Log-normal residual SD for BSAP (natural-log scale)")
    # Table 2, EPIC 2 eps_BSAP (SD) = 0.184 (%CV 4.6).

    expSd_BMD_LS <- 0.022; label("Log-normal residual SD for LS-BMD (natural-log scale)")
    # Table 2, EPIC 2 eps_LS (SD) = 0.022 (%CV 2.7).

    expSd_BMD_TH <- 0.020; label("Log-normal residual SD for TH-BMD (natural-log scale)")
    # Table 2, EPIC 2 eps_TH (SD) = 0.020 (%CV 2.3).
  })

  model({
    # =====================================================================
    # 1. Placebo activity switch. PCa modulates the osteoclast-driving
    #    factor fdbf so that placebo (calcium) inhibits RANK receptor
    #    occupancy for t >= T_ENTRY. Before T_ENTRY the placebo is off
    #    (PCa = 1) and the state evolves under pure disease progression.
    #    Guard the exponentials by clamping the argument to 0 when t <
    #    T_ENTRY so PCa collapses to 1 identically -- the ifelse checks
    #    both branches to keep the compiled expression well-defined.
    # =====================================================================

    dt_placebo <- ifelse(time < T_ENTRY, 0, time - T_ENTRY)
    plac_on    <- ifelse(time < T_ENTRY, 0, 1)
    PCa <- 1 - plac_on * (1 - exp(-k_Ca_onset * dt_placebo)) * exp(-k_Ca_offset * dt_placebo)

    # =====================================================================
    # 2. Disease progression f(t) = exp(-k_estrogen * t). t = 0 at
    #    menopause onset; k_estrogen is the first-order decay rate of the
    #    estrogen-driven disease-progression signal.
    # =====================================================================

    f_dspg <- exp(-k_estrogen * time)

    # =====================================================================
    # 3. Osteoclast / osteoblast auxiliary terms (paper Eq. 1, NONMEM $DES).
    #    pi_1 = 1/(1+z_s) is the baseline value of pi_z (z = 1).
    #    piz(z) = z / (z + z_s); piz1 = piz / pi_1;
    #    pi_1z = piz1^2 (appears squared in the fdbf denominator);
    #    DFR = 1 / (1 + b) with b = 1 gives DFR = 0.5.
    # =====================================================================

    pi_1  <- 1 / (1 + zs)
    piz   <- osteoclast / (osteoclast + zs)
    piz1  <- piz / pi_1
    pi_1z <- piz1 * piz1
    DFR   <- 1 / (1 + b_baseline)

    # =====================================================================
    # 4. Osteoclast driving factor fdbf and the two mechanistic ODE
    #    right-hand sides (paper Eq. 1; NONMEM $DES lines DADT(2) / DADT(3)).
    #    Both states start at 1 (relative baseline activity) at menopause
    #    onset.
    # =====================================================================

    fdbf <- (osteoblast / (1 + f_dspg * pi_1z)) / DFR * PCa

    d/dt(osteoblast) <- k_B * (piz1 - osteoblast)
    d/dt(osteoclast) <- D_A * pi_1 * (fdbf - piz1 * osteoclast)

    osteoblast(0) <- 1
    osteoclast(0) <- 1

    # =====================================================================
    # 5. Biomarker transducers (paper Eq. 4 with the k_BSAP0 unit-scaling
    #    modification of paper Eq. 6). Both NTX_0 and BSAP_0 carry
    #    exponential IIV via the etaNTX_0 / etaBSAP_0 block.
    #      NTX  = NTX_0 * exp(eta) * z^q_NTX
    #      BSAP = (BSAP_0 * (1 + k_BSAP0)) * exp(eta) * y^q_BSAP
    # =====================================================================

    NTX_baseline_i  <- NTX_0 * exp(etaNTX_0)
    BSAP_baseline_i <- BSAP_0 * (1 + k_BSAP0) * exp(etaBSAP_0)

    NTX  <- NTX_baseline_i  * osteoclast^q_NTX
    BSAP <- BSAP_baseline_i * osteoblast^q_BSAP

    # =====================================================================
    # 6. BMD indirect-response ODEs (paper Eq. 6; NONMEM $DES DADT(6/7)).
    #    Baselines are BMI-corrected around the cohort median 25.34114
    #    kg/m^2 and carry exponential IIV. k_out is derived per NONMEM
    #    $PK block as k_in / BMD_baseline, so the steady-state at
    #    (y = z = 1) has d/dt(BMD) = k_in * (D_AOB - D_AOC), which is
    #    small but non-zero (see vignette Errata for the discussion).
    # =====================================================================

    BMD_LS_baseline_i <- BMD_LS_0 * (1 + BMI_frac_LS * (BMI - 25.34114)) * exp(etaBMD_LS_0)
    BMD_TH_baseline_i <- BMD_TH_0 * (1 + BMI_frac_TH * (BMI - 25.34114)) * exp(etaBMD_TH_0)

    k_out_LS <- k_in_LS / BMD_LS_baseline_i
    k_out_TH <- k_in_TH / BMD_TH_baseline_i

    d/dt(BMD_LS) <- k_in_LS * (1 + D_AOB * osteoblast) - k_out_LS * (1 + D_AOC * osteoclast) * BMD_LS
    d/dt(BMD_TH) <- k_in_TH * (1 + D_AOB * osteoblast) - k_out_TH * (1 + D_AOC * osteoclast) * BMD_TH

    BMD_LS(0) <- BMD_LS_baseline_i
    BMD_TH(0) <- BMD_TH_baseline_i

    # =====================================================================
    # 7. Observation model. LTBS in the source: ln(y_obs) = ln(y_pred) +
    #    eps, encoded here via '~ lnorm(sd)'. The paper's per-observation
    #    'extreme-quantile' second-tier residual variance (SW4/SW5 switch
    #    in $ERROR) is documented in vignette Errata but not encoded here.
    # =====================================================================

    NTX    ~ lnorm(expSd_NTX)
    BSAP   ~ lnorm(expSd_BSAP)
    BMD_LS ~ lnorm(expSd_BMD_LS)
    BMD_TH ~ lnorm(expSd_BMD_TH)
  })
}
