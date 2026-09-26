Choi_2016_gcc4401c <- function() {
  description <- paste(
    "Population PK/PD model for GCC-4401C, an oral direct factor Xa",
    "inhibitor, in healthy male volunteers (Choi 2016; pooled first-in-human",
    "single-ascending-dose and single-and-multiple-ascending-dose phase I",
    "studies). Plasma and urine data were fit simultaneously with a",
    "two-compartment model with sequential zero-order release into the depot",
    "(D1), an absorption lag (ALAG1) and first-order absorption; body weight",
    "enters Vc as a power function normalised to 75 kg. Elimination is split",
    "into a linear non-renal clearance and a saturable renal clearance whose",
    "magnitude is suppressed by an inhibitory Emax function of the plasma",
    "concentration. A study-specific scaling of the apparent Vc encodes the",
    "relative accuracy of the two bioanalytical assays. Eight pharmacodynamic",
    "markers are carried as direct-effect (no-delay) outputs driven by the",
    "plasma concentration: coagulation factor X activity, factor X",
    "chromogenic activity, anti-factor Xa activity, prothrombin time in INR",
    "and in seconds, activated partial thromboplastin time, antithrombin III",
    "activity and the low-molecular-weight-heparin anti-Xa assay. The PD",
    "layer was fit sequentially on individual Bayesian PK estimates. Baselines",
    "for aPTT, AT III and LMWH are not reported anywhere in the paper or its",
    "supplement, so those three outputs are the drug-induced CHANGE from",
    "baseline and are named with a _chg suffix; see the vignette Errata.",
    "Companion rivaroxaban model: modellib('Choi_2016_rivaroxaban').",
    sep = " "
  )
  reference <- paste(
    "Choi HY, Choi S, Kim YH, Lim HS.",
    "Population pharmacokinetic and pharmacodynamic modeling analysis of",
    "GCC-4401C, a novel direct factor Xa inhibitor, in healthy volunteers.",
    "CPT Pharmacometrics Syst Pharmacol. 2016 Oct;5(10):532-543.",
    "doi:10.1002/psp4.12103.",
    "Structural detail and the plasma/urine NONMEM control stream are taken",
    "from Supplementary Data file PSP4-5-532-s008 accompanying the article;",
    "the coagulation-factor-X PD control stream is PSP4-5-532-s009.",
    sep = " "
  )
  vignette <- "Choi_2016_gcc4401c_rivaroxaban"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function on Vc normalised to 75 kg (Choi 2016 Eq. 5:",
        "Typical Vc = Vc(75) * (WT/75)^H). Mean (SD) weight in the pooled",
        "analysis population was 76.6 (9.7) kg (Table 1)."
      ),
      source_name = "WT"
    ),
    STUDY_SMAD = list(
      description = paste(
        "Study indicator: 1 = the single-and-multiple-ascending-dose study",
        "(NCT01954238), 0 = the first-in-human single-ascending-dose study",
        "(NCT01651234)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single-ascending-dose study, NCT01651234)",
      notes = paste(
        "Switches three things. (i) The apparent central volume, encoding the",
        "relative accuracy of the two bioanalytical assays: the supplementary",
        "control stream sets S2 = V2/1000 and, for STUDY = 2,",
        "S2 = (V2/1000)/THETA(10) with THETA(10) = 0.82, i.e. plasma",
        "concentrations measured in the S&MAD study are about 18% lower on",
        "average than those of the SAD study (Choi 2016 Eqs. 7-8 and Table 2",
        "footnote d). (ii) The sigmoidicity of the coagulation-factor-X PD",
        "model (Table 3a). (iii) The EC50 of the two prothrombin-time PD",
        "models and the Emax of the aPTT PD model (Table 3d, 3e, 3f)."
      ),
      source_name = "STUDY"
    ),
    OCC = list(
      description = "Integer occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Choi 2016 implemented IOV on D1 and ALAG1 'in every inter-dose",
        "interval, during which plasma were drawn for PK and the period after",
        "the last dose'. The supplementary plasma/urine control stream",
        "(PSP4-5-532-s008) carries seven occasion slots (ETA(6)-ETA(12) for D1",
        "and ETA(13)-ETA(19) for ALAG1, each $OMEGA BLOCK(1) SAME), so OCC",
        "takes values 1-7 here. The PD control stream (PSP4-5-532-s009) groups",
        "the same occasion column into three IOV slots -- OCC < 6, OCC == 6 and",
        "OCC == 7 -- and that grouping is reproduced for every PD parameter",
        "the paper reports a separate IOV for."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "GCC-4401C", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "GCC-4401C", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "GCC-4401C", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "GCC-4401C", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 2L,
    n_observations = 1689L,
    age_range = "mean (SD) 30.3 (8.6) years",
    weight_range = "mean (SD) 76.6 (9.7) kg",
    height_range = "mean (SD) 176.5 (6.7) cm",
    sex_female_pct = 0,
    race_ethnicity = c(White = 43.6, Black = 42.6, Asian = 8.5, Other = 5.3),
    disease_state = "Healthy male volunteers.",
    dose_range = paste(
      "Single oral doses of 2.5, 5, 10, 20, 40 or 80 mg under overnight",
      "fasting (SAD study, 48 subjects, 6 active + 2 placebo per group);",
      "single dose on day 1 then once daily on days 3-9 at 10, 20, 40, 60 or",
      "80 mg under overnight fasting (S&MAD study, 46 subjects)."
    ),
    regions = "United States (ClinicalTrials.gov NCT01651234 and NCT01954238).",
    notes = paste(
      "Baseline demographics from Choi 2016 Table 1; the two studies did not",
      "differ materially. PK dataset: 1401 plasma GCC-4401C concentrations",
      "(576 SAD, 825 S&MAD), 120 urine concentrations and 168 rivaroxaban",
      "plasma concentrations (the last analysed in the companion model).",
      "Urine was collected over 0-4, 4-8, 8-12, 12-24 and 24-48 h in the 10,",
      "20, 40 and 80 mg SAD groups only. PD markers were measured at the PK",
      "sampling times. Estimation used NONMEM 7.2 ADVAN6 with FOCE-I.",
      "Ecarin-stimulated thrombin activity (ESTA) was measured but could not",
      "be described by any of the tested PD models -- the observed values were",
      "not monotone in dose -- so it is not part of this model."
    )
  )

  ini({
    # =====================================================================
    # PHARMACOKINETICS -- Choi 2016 Table 2(a) (plasma + urine GCC-4401C,
    # both studies). Structure verified line-by-line against the
    # supplementary NONMEM control stream PSP4-5-532-s008.
    # =====================================================================
    lka <- log(3.50)
    label("First-order absorption rate constant out of the depot, Ka (1/h)") # Table 2(a) Ka = 3.50 (RSE 19.7%; 95% CI 2.15-4.85)
    ld1 <- log(0.395)
    label("Duration of zero-order release into the depot, D1 (h)") # Table 2(a) D1 = 0.395 (RSE 10.5%; 95% CI 0.314-0.476)
    ltlag <- log(0.20)
    label("Absorption lag time, ALAG1 (h)") # Table 2(a) ALAG1 = 0.20 (RSE 11.1%; 95% CI 0.16-0.24)
    lvc <- log(55.7)
    label("Central volume of distribution in a 75 kg subject, Vc (L)") # Table 2(a) Vc = 55.7 (RSE 5.8%; 95% CI 49.4-62.0)
    e_wt_vc <- 0.67
    label("Power exponent of body weight on Vc, normalised to 75 kg (unitless)") # Table 2(a) 'H for WT and Vc' = 0.67 (RSE 50.2%); Eq. 5
    lvp <- log(27.8)
    label("Peripheral volume of distribution, Vp (L)") # Table 2(a) Vp = 27.8 (RSE 7.0%; 95% CI 24.0-31.6)
    lq <- log(2.99)
    label("Inter-compartmental clearance, Q (L/h)") # Table 2(a) Q = 2.99 (RSE 4.8%; 95% CI 2.71-3.27)
    lcl_nonren <- log(11.2)
    label("Non-renal clearance, CLNR (L/h)") # Table 2(a) CLNR = 11.2 (RSE 3.8%; 95% CI 10.4-12.0)
    lcl_renal <- log(0.78)
    label("Baseline renal clearance before saturation, BCLR (L/h)") # Table 2(a) CLR = 0.78 (RSE 9.7%; 95% CI 0.63-0.93)
    logitimax <- logit(0.88)
    label("Maximum fractional inhibition of renal clearance, IMAX (unitless)") # Table 2(a) IMAX = 0.88 (RSE 7.5%; 95% CI 0.75-1.00); NONMEM bound (0, 0.9, 1)
    lic50 <- log(166.0)
    label("Plasma concentration at half-maximal inhibition of renal clearance, IC50 (ng/mL)") # Table 2(a) IC50 = 166.0 (RSE 37.1%; 95% CI 45.3-286.7)
    e_study_smad_vc <- 0.82
    label("Ratio of measured plasma concentration in the S&MAD study to the SAD study; equivalently the apparent Vc is divided by this factor in the S&MAD study (unitless)") # Table 2(a) 'Assay' = 0.82 (RSE 4.6%; 95% CI 0.74-0.89); Eqs. 7-8

    # ---- PK inter-individual variability (Table 2(a); variances on the
    # ---- log scale, CV% in parentheses). The supplementary control stream
    # ---- fits Vc, CLNR and CLR as an $OMEGA BLOCK(3); the paper reports
    # ---- only the three diagonal variances, so the covariances are not
    # ---- reproduced here (see vignette Errata).
    etalka ~ 0.97 # Table 2(a) 'IIV Ka' = 0.97 (CV 128.0%; RSE 35.2%)
    etalvc ~ 0.09 # Table 2(a) 'IIV Vc' = 0.09 (CV 30.5%; RSE 39.9%)
    etalvp ~ 0.11 # Table 2(a) 'IIV Vp' = 0.11 (CV 33.8%; RSE 26.2%)
    etalcl_nonren ~ 0.046 # Table 2(a) 'IIV CLNR' = 0.046 (CV 21.7%; RSE 30.4%)
    etalcl_renal ~ 0.052 # Table 2(a) 'IIV CLR' = 0.052 (CV 23.2%; RSE 68.8%)

    # ---- Combined IIV + IOV on D1 and ALAG1, seven occasions. The paper
    # ---- states IIV and IOV could not be estimated separately for these
    # ---- two parameters, so a single lumped variance is carried on the
    # ---- occasion level (control stream: D1 = THETA(8)*EXP(IOV) with
    # ---- ETA(6)-ETA(12) $OMEGA BLOCK(1) SAME, and likewise ALAG1 with
    # ---- ETA(13)-ETA(19)). Occasions 2-7 are fixed to the occasion-1
    # ---- variance to encode BLOCK(1) SAME.
    etaiov_d1_1 ~ 1.49 # Table 2(a) 'IIV + IOV D1' = 1.49 (CV 185.4%; RSE 35.0%)
    etaiov_d1_2 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_d1_3 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_d1_4 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_d1_5 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_d1_6 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_d1_7 ~ fixed(1.49) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_1 ~ 0.14 # Table 2(a) 'IIV + IOV ALAG1' = 0.14 (CV 38.8%; RSE 37.2%)
    etaiov_tlag_2 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_3 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_4 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_5 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_6 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME
    etaiov_tlag_7 ~ fixed(0.14) # $OMEGA BLOCK(1) SAME

    # ---- PK residual error. The control stream writes
    # ---- W = sqrt(THETA_add^2 + THETA_prop^2 * IPRED^2) with the additive
    # ---- term fixed at 1e-6, i.e. purely proportional, and $SIGMA 1 FIX.
    propSd <- 0.20
    label("Proportional residual error on plasma concentration (fraction)") # Table 2(a) 'e (proportional), plasma' = 0.20 (SD; RSE 3.0%)
    propSd_Ae <- 0.28
    label("Proportional residual error on the amount excreted in urine (fraction)") # Table 2(a) 'e (proportional), urine' = 0.28 (SD; RSE 10.1%)

    # =====================================================================
    # PHARMACODYNAMICS -- Choi 2016 Table 3. Every marker is a direct-effect
    # model: the paper tested ligand-receptor association-dissociation
    # models and found no time delay, and the concentration-effect plots
    # showed no hysteresis. The generic inhibitory form, from the
    # supplementary CFX control stream PSP4-5-532-s009, is
    #   E = BASE - Emax * Cc^gamma / (EC50^gamma + Cc^gamma)
    # with the stimulatory markers taking the same form with a plus sign.
    #
    # BASELINES. The paper used each subject's own observed pre-dose value
    # as BASE and did not parameterise it ("Baseline PD values of each
    # endpoint were used as such in each model without parameterization").
    # Where Supplementary Figure S4 plots the simulated marker, the typical
    # baseline is read off the pre-dose median and carried here as a fixed
    # parameter; where it does not, the output is the CHANGE from baseline
    # and carries a _chg suffix. No baseline has been invented.
    # =====================================================================

    # ---- (a) Coagulation factor X activity, % (both studies; sigmoid Emax)
    lrbase_cfx <- fixed(log(100))
    label("Typical baseline coagulation factor X activity (%)") # Digitised from Supplementary Figure S4 (Coagulation Factor X panels): pre-dose median 100%; not estimated by the paper
    lemax_cfx <- log(81.2)
    label("Maximum decrease in coagulation factor X activity, Emax (%)") # Table 3(a) Emax = 81.2 (RSE 11.3%; 95% CI 63.2-99.2)
    lec50_cfx <- log(2880.0)
    label("GCC-4401C concentration at half-maximal factor X inhibition, EC50 (ng/mL)") # Table 3(a) EC50 = 2880.0 (RSE 34.1%; 95% CI 957.2-4802.8)
    lhill_cfx_sad <- log(1.07)
    label("Sigmoidicity of the factor X model in the SAD study, gamma (unitless)") # Table 3(a) gamma_SAD = 1.07 (RSE 14.3%; 95% CI 0.77-1.37)
    lhill_cfx_smad <- log(0.60)
    label("Sigmoidicity of the factor X model in the S&MAD study, gamma (unitless)") # Table 3(a) gamma_SAD/MAD = 0.60 (RSE 11.7%; 95% CI 0.46-0.74)
    etalemax_cfx ~ 0.10 # Table 3(a) 'IIV Emax' = 0.10 (CV 33.1%; RSE 48.5%)
    etalec50_cfx ~ 0.98 # Table 3(a) 'IIV EC50' = 0.98 (CV 129.4%; RSE 39.0%)
    etalhill_cfx ~ 0.24 # Table 3(a) 'IIV gamma' = 0.24 (CV 51.8%; RSE 32.7%)
    etaiov_hill_cfx_1 ~ 0.09 # Table 3(a) 'IOV gamma' = 0.09 (CV 30.1%; RSE 35.3%)
    etaiov_hill_cfx_2 ~ fixed(0.09) # $OMEGA BLOCK(1) SAME (PSP4-5-532-s009)
    etaiov_hill_cfx_3 ~ fixed(0.09) # $OMEGA BLOCK(1) SAME (PSP4-5-532-s009)
    addSd_cfx <- 3.94
    label("Additive residual error on coagulation factor X activity (%)") # Table 3(a) 'e (additive), %' = 3.94 (RSE 22.5%)
    propSd_cfx <- 0.06
    label("Proportional residual error on coagulation factor X activity (fraction)") # Table 3(a) 'e (proportional)' = 0.06 (SD; RSE 12.9%)

    # ---- (b) Factor X chromogenic activity assay, % (S&MAD study; sigmoid Emax)
    lrbase_fxcaa <- fixed(log(100))
    label("Typical baseline factor X chromogenic activity (%)") # Digitised from Supplementary Figure S4 (Factor X Chromogenic Activity Assay panels): pre-dose median 100%; not estimated by the paper
    lemax_fxcaa <- log(99.8)
    label("Maximum decrease in factor X chromogenic activity, Emax (%)") # Table 3(b) Emax = 99.8 (RSE 7.3%; 95% CI 85.5-114.1)
    lec50_fxcaa <- log(420.0)
    label("GCC-4401C concentration at half-maximal chromogenic factor X inhibition, EC50 (ng/mL)") # Table 3(b) EC50 = 420.0 (RSE 16.4%; 95% CI 284.8-555.2)
    lhill_fxcaa <- log(0.90)
    label("Sigmoidicity of the chromogenic factor X model, gamma (unitless)") # Table 3(b) gamma = 0.90 (RSE 7.5%; 95% CI 0.76-1.03)
    etalemax_fxcaa ~ 0.01 # Table 3(b) 'IIV Emax' = 0.01 (CV 11.3%; RSE 52.1%)
    etalec50_fxcaa ~ 0.08 # Table 3(b) 'IIV + IOV EC50' = 0.08 (RSE 42.1%); lumped, carried at the subject level
    etalhill_fxcaa ~ 0.06 # Table 3(b) 'IIV gamma' = 0.06 (CV 25.7%; RSE 44.9%)
    addSd_fxcaa <- 3.03
    label("Additive residual error on factor X chromogenic activity (%)") # Table 3(b) 'e (additive), %' = 3.03 (RSE 34.7%)
    propSd_fxcaa <- 0.06
    label("Proportional residual error on factor X chromogenic activity (fraction)") # Table 3(b) 'e (proportional)' = 0.06 (SD; RSE 18.3%)

    # ---- (c) Anti-factor Xa activity, IU/mL (S&MAD study; sigmoid Emax, stimulatory)
    lrbase_afx <- fixed(log(0.05))
    label("Typical baseline anti-factor Xa activity (IU/mL)") # Digitised from Supplementary Figure S4 (Anti-Factor Xa panels): pre-dose median sits on the axis, ~0.05 IU/mL; not estimated by the paper
    lemax_afx <- log(3.24)
    label("Maximum increase in anti-factor Xa activity, Emax (IU/mL)") # Table 3(c) Emax = 3.24 (RSE 21.3%; 95% CI 1.89-4.59)
    lec50_afx <- log(695.0)
    label("GCC-4401C concentration at half-maximal anti-factor Xa effect, EC50 (ng/mL)") # Table 3(c) EC50 = 695.0 (RSE 26.0%; 95% CI 340.2-1049.8)
    lhill_afx <- log(1.25)
    label("Sigmoidicity of the anti-factor Xa model, gamma (unitless)") # Table 3(c) gamma = 1.25 (RSE 6.0%; 95% CI 1.10-1.40)
    etalemax_afx ~ 0.04 # Table 3(c) 'IIV Emax' = 0.04 (CV 20.8%; RSE 62.3%)
    etalec50_afx ~ 0.005 # Table 3(c) 'IIV + IOV EC50' = 0.005 (RSE 48.8%); lumped, carried at the subject level
    etalhill_afx ~ 0.02 # Table 3(c) 'IIV gamma' = 0.02 (CV 13.4%; RSE 51.4%)
    addSd_afx <- 0.04
    label("Additive residual error on anti-factor Xa activity (IU/mL)") # Table 3(c) 'e (additive), IU/mL' = 0.04 (RSE 12.2%)
    propSd_afx <- 0.14
    label("Proportional residual error on anti-factor Xa activity (fraction)") # Table 3(c) 'e (proportional)' = 0.14 (SD; RSE 8.1%)

    # ---- (d) Prothrombin time, INR (both studies; sigmoid Emax, stimulatory)
    lrbase_ptinr <- fixed(log(1.05))
    label("Typical baseline prothrombin time (INR)") # Digitised from Supplementary Figure S4 (PT (INR) panels): pre-dose median 1.05; not estimated by the paper
    lemax_ptinr <- log(1.32)
    label("Maximum increase in prothrombin time, Emax (INR)") # Table 3(d) Emax = 1.32 (RSE 20.3%; 95% CI 0.79-1.85)
    lec50_ptinr_sad <- log(426.0)
    label("Concentration at half-maximal PT (INR) effect in the SAD study, EC50 (ng/mL)") # Table 3(d) EC50_SAD = 426.0 (RSE 24.4%; 95% CI 222.2-629.8)
    lec50_ptinr_smad <- log(1350.0)
    label("Concentration at half-maximal PT (INR) effect in the S&MAD study, EC50 (ng/mL)") # Table 3(d) EC50_SAD/MAD = 1350.0 (RSE 27.3%; 95% CI 628.7-2071.3)
    lhill_ptinr <- log(1.23)
    label("Sigmoidicity of the PT (INR) model, gamma (unitless)") # Table 3(d) gamma = 1.23 (RSE 7.5%; 95% CI 1.05-1.41)
    etalemax_ptinr ~ 0.48 # Table 3(d) 'IIV Emax' = 0.48 (CV 78.0%; RSE 43.0%)
    etalec50_ptinr ~ 0.89 # Table 3(d) 'IIV EC50' = 0.89 (CV 120.1%; RSE 52.3%)
    etalhill_ptinr ~ 0.11 # Table 3(d) 'IIV gamma' = 0.11 (CV 34.3%; RSE 74.1%)
    etaiov_hill_ptinr_1 ~ 0.04 # Table 3(d) 'IOV gamma' = 0.04 (CV 19.9%; RSE 69.6%)
    etaiov_hill_ptinr_2 ~ fixed(0.04) # $OMEGA BLOCK(1) SAME
    etaiov_hill_ptinr_3 ~ fixed(0.04) # $OMEGA BLOCK(1) SAME
    propSd_ptinr <- 0.05
    label("Proportional residual error on prothrombin time in INR (fraction)") # Table 3(d) 'e (proportional)' = 0.05 (SD; RSE 7.7%); no additive term reported
    # ---- (e) Prothrombin time, seconds (both studies; sigmoid Emax, stimulatory)
    lrbase_ptsec <- fixed(log(12.2))
    label("Typical baseline prothrombin time (s)") # Digitised from Supplementary Figure S4 (PT (sec) panels): pre-dose median 12.2 s; not estimated by the paper
    lemax_ptsec <- log(15.2)
    label("Maximum increase in prothrombin time, Emax (s)") # Table 3(e) Emax = 15.2 (RSE 21.5%; 95% CI 8.8-21.6)
    lec50_ptsec_sad <- log(563.0)
    label("Concentration at half-maximal PT (s) effect in the SAD study, EC50 (ng/mL)") # Table 3(e) EC50_SAD = 563.0 (RSE 34.3%; 95% CI 184.7-941.3)
    lec50_ptsec_smad <- log(1450.0)
    label("Concentration at half-maximal PT (s) effect in the S&MAD study, EC50 (ng/mL)") # Table 3(e) EC50_SAD/MAD = 1450.0 (RSE 30.0%; 95% CI 597.4-2302.6)
    lhill_ptsec <- log(1.16)
    label("Sigmoidicity of the PT (s) model, gamma (unitless)") # Table 3(e) gamma = 1.16 (RSE 8.0%; 95% CI 0.98-1.34)
    etalemax_ptsec ~ 0.41 # Table 3(e) 'IIV Emax' = 0.41 (CV 71.6%; RSE 64.7%)
    etalec50_ptsec ~ 0.77 # Table 3(e) 'IIV EC50' = 0.77 (CV 107.8%; RSE 64.5%)
    etalhill_ptsec ~ 0.04 # Table 3(e) 'IIV gamma' = 0.04 (CV 18.9%; RSE 65.0%)
    etaiov_hill_ptsec_1 ~ 0.03 # Table 3(e) 'IOV gamma' = 0.03 (CV 17.8%; RSE 60.1%)
    etaiov_hill_ptsec_2 ~ fixed(0.03) # $OMEGA BLOCK(1) SAME
    etaiov_hill_ptsec_3 ~ fixed(0.03) # $OMEGA BLOCK(1) SAME
    propSd_ptsec <- 0.04
    label("Proportional residual error on prothrombin time in seconds (fraction)") # Table 3(e) 'e (proportional)' = 0.04 (SD; RSE 6.9%); no additive term reported

    # ---- (f) Activated partial thromboplastin time, seconds (both studies;
    # ---- sigmoid Emax, stimulatory). No baseline aPTT is reported or
    # ---- plotted in any available source, so the output is the PROLONGATION.
    lemax_aptt_sad <- log(16.9)
    label("Maximum prolongation of aPTT in the SAD study, Emax (s)") # Table 3(f) Emax_SAD = 16.9 (RSE 13.1%; 95% CI 12.5-21.3)
    lemax_aptt_smad <- log(20.4)
    label("Maximum prolongation of aPTT in the S&MAD study, Emax (s)") # Table 3(f) Emax_SAD/MAD = 20.4 (RSE 12.5%; 95% CI 15.4-25.4)
    lec50_aptt <- log(573.0)
    label("Concentration at half-maximal aPTT prolongation, EC50 (ng/mL)") # Table 3(f) EC50 = 573.0 (RSE 22.7%; 95% CI 318.2-827.8)
    lhill_aptt <- log(1.37)
    label("Sigmoidicity of the aPTT model, gamma (unitless)") # Table 3(f) gamma = 1.37 (RSE 13.2%; 95% CI 1.02-1.72)
    etalemax_aptt ~ 0.09 # Table 3(f) 'IIV Emax' = 0.09 (CV 30.3%; RSE 86.1%)
    etalec50_aptt ~ 0.33 # Table 3(f) 'IIV EC50' = 0.33 (CV 62.4%; RSE 62.3%)
    etalhill_aptt ~ 0.19 # Table 3(f) 'IIV gamma' = 0.19 (CV 45.9%; RSE 39.2%)
    etaiov_hill_aptt_1 ~ 0.03 # Table 3(f) 'IOV gamma' = 0.03 (CV 17.8%; RSE 35.5%)
    etaiov_hill_aptt_2 ~ fixed(0.03) # $OMEGA BLOCK(1) SAME
    etaiov_hill_aptt_3 ~ fixed(0.03) # $OMEGA BLOCK(1) SAME
    propSd_aptt_chg <- 0.05
    label("Proportional residual error on the aPTT prolongation (fraction)") # Table 3(f) 'e (proportional)' = 0.05 (SD; RSE 1.0%)

    # ---- (g) Antithrombin III activity (SAD study; linear). No baseline is
    # ---- reported or plotted, so the output is the CHANGE from baseline.
    lslope_atiii <- log(0.006)
    label("Linear slope of the antithrombin III response on plasma concentration (activity units per ng/mL)") # Table 3(g) SLOPE = 0.006 (RSE 44.8%; 95% CI 0.001-0.011)
    etalslope_atiii ~ 3.73 # Table 3(g) 'IIV SLOPE' = 3.73 (CV 637.8%; RSE 38.6%)
    addSd_atiii_chg <- 9.28
    label("Additive residual error on the antithrombin III change from baseline (activity units)") # Table 3(g) 'e (additive)' = 9.28 (RSE 10.9%)

    # ---- (h) Low-molecular-weight-heparin anti-Xa assay (SAD study; simple
    # ---- Emax). No baseline is reported or plotted, so the output is the
    # ---- CHANGE from baseline.
    lemax_lmwh <- log(3.83)
    label("Maximum increase in the LMWH anti-Xa assay, Emax (assay units; unit not stated in Table 3(h))") # Table 3(h) Emax = 3.83 (RSE 11.4%; 95% CI 2.98-4.68)
    lec50_lmwh <- log(759.0)
    label("Concentration at half-maximal LMWH anti-Xa effect, EC50 (ng/mL)") # Table 3(h) EC50 = 759.0 (RSE 15.2%; 95% CI 533.6-984.4)
    etalec50_lmwh ~ 0.29 # Table 3(h) 'IIV EC50' = 0.29 (CV 57.4%; RSE 33.9%)
    addSd_lmwh_chg <- 0.08
    label("Additive residual error on the LMWH anti-Xa change from baseline (assay units)") # Table 3(h) 'e (additive)' = 0.08 (RSE 20.6%)
    propSd_lmwh_chg <- 0.17
    label("Proportional residual error on the LMWH anti-Xa change from baseline (fraction)") # Table 3(h) 'e (proportional)' = 0.17 (SD; RSE 14.3%)
  })

  model({
    # =====================================================================
    # 1. Occasion indicators. Seven PK occasions (control stream
    #    PSP4-5-532-s008); the PD control stream PSP4-5-532-s009 groups the
    #    same column into three slots.
    # =====================================================================
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    poc1 <- (OCC < 6)
    poc2 <- (OCC == 6)
    poc3 <- (OCC == 7)

    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2 + oc3 * etaiov_d1_3 +
      oc4 * etaiov_d1_4 + oc5 * etaiov_d1_5 + oc6 * etaiov_d1_6 +
      oc7 * etaiov_d1_7
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2 + oc3 * etaiov_tlag_3 +
      oc4 * etaiov_tlag_4 + oc5 * etaiov_tlag_5 + oc6 * etaiov_tlag_6 +
      oc7 * etaiov_tlag_7

    # =====================================================================
    # 2. Individual PK parameters
    # =====================================================================
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + iov_d1)
    tlag <- exp(ltlag + iov_tlag)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc # Eq. 5
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- exp(lcl_renal + etalcl_renal)
    imax <- expit(logitimax)
    ic50 <- exp(lic50)

    # Study-specific concentration scaling. The control stream sets
    # S2 = Vc/1000 and, for the S&MAD study, S2 = (Vc/1000)/0.82, i.e. the
    # predicted concentration is multiplied by 0.82 (Eqs. 7-8).
    assay_smad <- 1 - STUDY_SMAD * (1 - e_study_smad_vc)

    # =====================================================================
    # 3. Micro-constants and the saturable renal elimination. `Cc` is
    #    computed before the ODEs because it drives the inhibitory Emax
    #    function on renal clearance (Eq. 6, control-stream $DES:
    #    K24 = CLR * (1 - IMX*C2/(C2 + IC50)) / V2).
    # =====================================================================
    Cc <- 1000 * central / vc * assay_smad # ng/mL (dose in mg, vc in L)
    cl_renal_eff <- cl_renal * (1 - imax * Cc / (Cc + ic50)) # Eq. 6
    k20 <- cl_nonren / vc
    k24 <- cl_renal_eff / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 -
      (k20 + k24) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(urine) <- k24 * central

    dur(depot) <- d1 # zero-order release; dose records need rate = -2
    alag(depot) <- tlag

    Ae <- urine # cumulative amount excreted in urine (mg)

    # =====================================================================
    # 4. Direct-effect pharmacodynamics (Table 3). No effect compartment:
    #    the paper found no hysteresis and no improvement from
    #    ligand-receptor association-dissociation models.
    # =====================================================================
    # (a) Coagulation factor X activity (%)
    emax_cfx <- exp(lemax_cfx + etalemax_cfx)
    ec50_cfx <- exp(lec50_cfx + etalec50_cfx)
    hill_cfx <- exp((1 - STUDY_SMAD) * lhill_cfx_sad + STUDY_SMAD * lhill_cfx_smad +
      etalhill_cfx + poc1 * etaiov_hill_cfx_1 + poc2 * etaiov_hill_cfx_2 +
      poc3 * etaiov_hill_cfx_3)
    cfx <- exp(lrbase_cfx) -
      emax_cfx * Cc^hill_cfx / (ec50_cfx^hill_cfx + Cc^hill_cfx)

    # (b) Factor X chromogenic activity assay (%)
    emax_fxcaa <- exp(lemax_fxcaa + etalemax_fxcaa)
    ec50_fxcaa <- exp(lec50_fxcaa + etalec50_fxcaa)
    hill_fxcaa <- exp(lhill_fxcaa + etalhill_fxcaa)
    fxcaa <- exp(lrbase_fxcaa) -
      emax_fxcaa * Cc^hill_fxcaa / (ec50_fxcaa^hill_fxcaa + Cc^hill_fxcaa)

    # (c) Anti-factor Xa activity (IU/mL)
    emax_afx <- exp(lemax_afx + etalemax_afx)
    ec50_afx <- exp(lec50_afx + etalec50_afx)
    hill_afx <- exp(lhill_afx + etalhill_afx)
    afx <- exp(lrbase_afx) +
      emax_afx * Cc^hill_afx / (ec50_afx^hill_afx + Cc^hill_afx)

    # (d) Prothrombin time (INR)
    emax_ptinr <- exp(lemax_ptinr + etalemax_ptinr)
    ec50_ptinr <- exp((1 - STUDY_SMAD) * lec50_ptinr_sad +
      STUDY_SMAD * lec50_ptinr_smad + etalec50_ptinr)
    hill_ptinr <- exp(lhill_ptinr + etalhill_ptinr +
      poc1 * etaiov_hill_ptinr_1 + poc2 * etaiov_hill_ptinr_2 +
      poc3 * etaiov_hill_ptinr_3)
    ptinr <- exp(lrbase_ptinr) +
      emax_ptinr * Cc^hill_ptinr / (ec50_ptinr^hill_ptinr + Cc^hill_ptinr)

    # (e) Prothrombin time (seconds)
    emax_ptsec <- exp(lemax_ptsec + etalemax_ptsec)
    ec50_ptsec <- exp((1 - STUDY_SMAD) * lec50_ptsec_sad +
      STUDY_SMAD * lec50_ptsec_smad + etalec50_ptsec)
    hill_ptsec <- exp(lhill_ptsec + etalhill_ptsec +
      poc1 * etaiov_hill_ptsec_1 + poc2 * etaiov_hill_ptsec_2 +
      poc3 * etaiov_hill_ptsec_3)
    ptsec <- exp(lrbase_ptsec) +
      emax_ptsec * Cc^hill_ptsec / (ec50_ptsec^hill_ptsec + Cc^hill_ptsec)

    # (f) Activated partial thromboplastin time: PROLONGATION (s), because
    #     no baseline aPTT is reported in the paper or its supplement.
    emax_aptt <- exp((1 - STUDY_SMAD) * lemax_aptt_sad +
      STUDY_SMAD * lemax_aptt_smad + etalemax_aptt)
    ec50_aptt <- exp(lec50_aptt + etalec50_aptt)
    hill_aptt <- exp(lhill_aptt + etalhill_aptt +
      poc1 * etaiov_hill_aptt_1 + poc2 * etaiov_hill_aptt_2 +
      poc3 * etaiov_hill_aptt_3)
    aptt_chg <- emax_aptt * Cc^hill_aptt / (ec50_aptt^hill_aptt + Cc^hill_aptt)

    # (g) Antithrombin III activity: CHANGE from baseline, because no
    #     baseline is reported in the paper or its supplement.
    slope_atiii <- exp(lslope_atiii + etalslope_atiii)
    atiii_chg <- slope_atiii * Cc

    # (h) LMWH anti-Xa assay: CHANGE from baseline, because no baseline is
    #     reported in the paper or its supplement.
    emax_lmwh <- exp(lemax_lmwh)
    ec50_lmwh <- exp(lec50_lmwh + etalec50_lmwh)
    lmwh_chg <- emax_lmwh * Cc / (ec50_lmwh + Cc)

    # =====================================================================
    # 5. Residual error
    # =====================================================================
    Cc ~ prop(propSd)
    Ae ~ prop(propSd_Ae)
    cfx ~ add(addSd_cfx) + prop(propSd_cfx)
    fxcaa ~ add(addSd_fxcaa) + prop(propSd_fxcaa)
    afx ~ add(addSd_afx) + prop(propSd_afx)
    ptinr ~ prop(propSd_ptinr)
    ptsec ~ prop(propSd_ptsec)
    aptt_chg ~ prop(propSd_aptt_chg)
    atiii_chg ~ add(addSd_atiii_chg)
    lmwh_chg ~ add(addSd_lmwh_chg) + prop(propSd_lmwh_chg)
  })
}
