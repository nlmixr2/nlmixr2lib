Rathi_2021_belantamab <- function() {
  description <- "Sequential three-analyte population PK model for the antibody-drug conjugate belantamab mafodotin in patients with relapsed/refractory multiple myeloma (Rathi 2021; DREAMM-1 and DREAMM-2). The ADC (output Cc, ug/mL) follows a linear two-compartment model with a sigmoidal time-varying clearance and covariate effects of baseline body weight, albumin, soluble BCMA, IgG, sex, dose below 1 mg/kg and study. First-order deconjugation (kdec) converts ADC into naked antibody, which shares the ADC distribution parameters; total antibody (output Cc_tab, nM) is the sum of ADC and naked antibody times an assay correction factor. The cys-mcMMAF payload (output Cc_mcmmaf, nM) follows a linear two-compartment model fed by the ADC clearance flux times an empirical drug-to-antibody ratio that decays exponentially with time after each dose and resets at every dose, with inter-occasion variability on its clearance and central volume."
  reference <- "Rathi C, Collins J, Struemper H, Opalinska J, Jewell RC, Ferron-Brady G. Population pharmacokinetics of belantamab mafodotin, a BCMA-targeting agent in patients with relapsed/refractory multiple myeloma. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):851-863. doi:10.1002/psp4.12660"
  vignette <- "Rathi_2021_belantamab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL for the ADC (Cc); nM for total antibody (Cc_tab) and cys-mcMMAF (Cc_mcmmaf)"
  )

  # What each ODE state holds. The ADC and naked-antibody states hold mass
  # amounts in mg (the dose unit); the payload states hold molar amounts in
  # nmol. The Supplementary Methods control streams dose AMT in mg and, for the
  # total-mAb and cys-mcMMAF fits, convert it to umol with F1 = 1/152.1 (ADC
  # molecular weight 152.1 kg/mol). Plasma throughout (Methods, 'Study
  # population and PK sampling').
  compartmentData <- list(
    central = list(analyte = "belantamab mafodotin (ADC)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "belantamab mafodotin (ADC)", units = "mg", specimen = "plasma", verified = TRUE),
    central_dar0 = list(
      analyte = "naked (unconjugated) belantamab antibody",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_dar0 = list(
      analyte = "naked (unconjugated) belantamab antibody",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_mcmmaf = list(analyte = "cys-mcMMAF", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_mcmmaf = list(analyte = "cys-mcMMAF", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline, centred at 75 kg (the DREAMM-2 typical patient, Methods 'Covariate effects on Cycle 1 ADC and cys-mcMMAF exposure'). Power exponents 0.507 on ADC CL and 0.500 on ADC V1 (Table 1; control-stream THETA(17) and THETA(10)).",
      source_name = "BWT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline, centred at 40 g/L. Power exponents -0.972 on ADC CL and -0.471 on ADC V1 (Table 1).",
      source_name = "BALB"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline, centred at 100 ng/mL. Power exponents 0.126 on ADC CL (Table 1) and 0.131 on the cys-mcMMAF central volume (Table S4). Missing values were imputed at the study-specific median in the source analysis.",
      source_name = "BSBCMA"
    ),
    IGG = list(
      description = "Baseline serum immunoglobulin G",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline, centred at 10 g/L. Power exponents 0.127 on ADC CL (Table 1) and 0.154 on the cys-mcMMAF central volume (Table S4).",
      source_name = "BIGG"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The control stream carries SEX = 1 for males and applies the male factor as THETA(14)**SEX on V1, so the typical V1 of 4.47 L is the female value and males take 4.47 * 1.14. Encoded as e_sexm_vc^(1 - SEXF), i.e. SEXF = 1 - SEX.",
      source_name = "SEX"
    ),
    STUDY_DREAMM1 = list(
      description = "Study indicator: 1 = DREAMM-1 (BMA117159, phase I), 0 = DREAMM-2 (205678, phase II)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (DREAMM-2)",
      notes = "The control streams carry DREAMM = 1 or 2 and apply the study effects as THETA**(2 - DREAMM), so STUDY_DREAMM1 = 2 - DREAMM. Multiplies ADC CL by 0.608 and V1 by 0.66, total-antibody CORR by 0.848 and the cys-mcMMAF central volume by (1 + 0.507). Also switches the cys-mcMMAF observation onto the DREAMM-1 VALA assay scale (Methods equation 1) and selects the DREAMM-1 residual error.",
      source_name = "DREAMM"
    ),
    DOSE_LOW = list(
      description = "Low-dose indicator: 1 = the administered belantamab mafodotin dose is below 1 mg/kg, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dose >= 1 mg/kg)",
      notes = "The control streams set ADOSE1 = 1 unless the actual dose ADOSE (mg/kg) exceeds 1. In the studied range this is the DREAMM-1 part 1 cohorts at 0.03, 0.06, 0.12, 0.24, 0.48 and 0.96 mg/kg (Table S1); no subject received exactly 1 mg/kg. Multiplies ADC CL by 0.623 and V2 by 0.129 (Table 1), kdec by 0.0001 (held constant, Table S3) and CORR by 0.846 (Table S3).",
      source_name = "ADOSE1"
    ),
    OCC = list(
      description = "Dosing occasion (cycle) for the cys-mcMMAF inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "The cys-mcMMAF control stream builds IOV from four occasion indicators: OCC = 1, 2, 3 and OCC > 3. Rathi 2021 Results attribute the IOV to 'the difference in profiles for different cycles', so OCC is the dosing-cycle number. Not used by the ADC or total-antibody sub-models.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 291L,
    n_studies = 2L,
    age_range = "34-89 years",
    age_median = "65 years",
    weight_range = "42.4-130 kg",
    weight_median = "74.0 kg",
    sex_female_pct = 46,
    race_ethnicity = c(White = 83),
    disease_state = "Relapsed/refractory multiple myeloma, heavily pretreated (DREAMM-2: refractory to a proteasome inhibitor and an immunomodulatory agent and refractory and/or intolerant to an anti-CD38 antibody).",
    dose_range = "Belantamab mafodotin IV every 3 weeks: 0.03-4.6 mg/kg as a 1-h infusion in DREAMM-1 and 2.5 or 3.4 mg/kg as a 30-min infusion in DREAMM-2 (Table S1).",
    regions = "Multinational (DREAMM-1 and DREAMM-2).",
    biomarker_summary = "Median (range): albumin 39.0 (19.0-50.0) g/L; sBCMA 80.0 (1.95-2030) ng/mL; IgG 8.50 (0.30-99.7) g/L; eGFR 71.6 (14.0-166) mL/min/1.73 m^2 (Table S2).",
    reference_subject = "DREAMM-2 typical patient: 65-year-old man, 75 kg, sBCMA 100 ng/mL, IgG 10 g/L, albumin 40 g/L, 2.5 mg/kg (Methods).",
    notes = "DREAMM-1 (NCT02064387; n = 73) and DREAMM-2 (NCT03525678; n = 218: 95 at 2.5 mg/kg and 99 at 3.4 mg/kg frozen liquid, 24 at 3.4 mg/kg lyophilised). 2248 ADC, 2210 total-mAb, 286 DREAMM-1 (VALA assay) and 1313 DREAMM-2 (PKC assay) cys-mcMMAF concentrations. Baseline demographics in Table S2."
  )

  ini({
    # ------------------------------------------------------------
    # ADC sub-model. Rathi 2021 Table 1 prints these to two decimals; the
    # unrounded final estimates are the values held constant in the
    # Supplementary Methods total-mAb and cys-mcMMAF control streams
    # ($THETA 1-23), which round to Table 1. Reference subject: female,
    # DREAMM-2, 75 kg, albumin 40 g/L, sBCMA 100 ng/mL, IgG 10 g/L,
    # dose >= 1 mg/kg.
    # ------------------------------------------------------------
    lcl <- log(0.936); label("ADC clearance at time zero, CL0 (L/day)") # Table 1 TVCL 0.94; control stream THETA(1) 0.936; Abstract '0.936 L/day'
    lvc <- log(4.47); label("ADC central volume V1 (L)") # Table 1 TVV1 4.47; control stream THETA(2)
    lq <- log(0.781); label("ADC intercompartmental clearance Q (L/day)") # Table 1 TVQ 0.78; control stream THETA(3) 0.781
    lvp <- log(5.96); label("ADC peripheral volume V2 (L)") # Table 1 TVV2 5.96; control stream THETA(4)

    cl_time_max <- -0.327; label("Maximal log-scale change in ADC CL over time, IMAX (unitless)") # Table 1 IMAX -0.33; control stream THETA(5) -0.327
    lcl_t50 <- log(53.9); label("Time to 50 percent of the maximal ADC CL change, TI50 (day)") # Table 1 TI50 53.90; control stream THETA(6)
    lcl_time_hill <- log(4.30); label("Hill coefficient of the time-varying ADC CL, GAMMA (unitless)") # Table 1 GAMMA 4.30; control stream THETA(7)

    e_wt_cl <- 0.507; label("Power exponent of body weight on ADC CL (unitless)") # Table 1 theta_CL_BWT 0.51; control stream THETA(17) 0.507
    e_wt_vc <- 0.500; label("Power exponent of body weight on ADC V1 (unitless)") # Table 1 theta_V1_BWT 0.50; control stream THETA(10) 0.500
    e_alb_cl <- -0.972; label("Power exponent of albumin on ADC CL (unitless)") # Table 1 theta_CL_BALB -0.97; control stream THETA(12) -0.972
    e_alb_vc <- -0.471; label("Power exponent of albumin on ADC V1 (unitless)") # Table 1 theta_V1_BALB -0.47; control stream THETA(13) -0.471
    e_sexm_vc <- 1.14; label("Multiplicative factor on ADC V1 for males (unitless)") # Table 1 theta_V1_SEX 1.14; control stream THETA(14)
    e_sbcma_cl <- 0.126; label("Power exponent of sBCMA on ADC CL (unitless)") # Table 1 theta_CL_BSBCMA 0.13; control stream THETA(15) 0.126
    e_igg_cl <- 0.127; label("Power exponent of IgG on ADC CL (unitless)") # Table 1 theta_CL_BIGG 0.13; control stream THETA(16) 0.127
    e_dose_low_vp <- 0.129; label("Multiplicative factor on ADC V2 for a dose below 1 mg/kg (unitless)") # Table 1 theta_V2_ADOSE<1 0.13; control stream THETA(20) 0.129
    e_dose_low_cl <- 0.623; label("Multiplicative factor on ADC CL for a dose below 1 mg/kg (unitless)") # Table 1 theta_CL_ADOSE<1 0.62; control stream THETA(21) 0.623
    e_study_dreamm1_vc <- 0.66; label("Multiplicative factor on ADC V1 for DREAMM-1 (unitless)") # Table 1 theta_V1_DREAMM1 0.66; control stream THETA(22)
    e_study_dreamm1_cl <- 0.608; label("Multiplicative factor on ADC CL for DREAMM-1 (unitless)") # Table 1 theta_CL_DREAMM1 0.61; control stream THETA(23) 0.608

    # ------------------------------------------------------------
    # Total-antibody sub-model (Rathi 2021 Table S3). The naked antibody
    # shares V1, V2 and Q with the ADC (Results, 'Population PK model for
    # total mAb').
    # ------------------------------------------------------------
    lkdec <- log(0.0590); label("ADC deconjugation rate constant to naked antibody, Kdec (1/day)") # Table S3 TVKDEC 0.0590
    e_dose_low_kdec <- fixed(0.0001); label("Multiplicative factor on Kdec for a dose below 1 mg/kg (unitless)") # Table S3 theta_Kdec_ADOSE<1 '0.0001 FIX'
    lcal_slope_tab <- log(1.14); label("Total-antibody assay correction factor, CORR (unitless)") # Table S3 TVCORR 1.14
    e_study_dreamm1_cal_slope_tab <- 0.848; label("Multiplicative factor on CORR for DREAMM-1 (unitless)") # Table S3 theta_CORR_DREAMM1 0.848
    e_dose_low_cal_slope_tab <- 0.846; label("Multiplicative factor on CORR for a dose below 1 mg/kg (unitless)") # Table S3 theta_CORR_ADOSE<1 0.846

    # ------------------------------------------------------------
    # cys-mcMMAF sub-model (Rathi 2021 Table S4). K34 / K43 are the
    # payload central-to-peripheral and peripheral-to-central rate
    # constants.
    # ------------------------------------------------------------
    lcl_mcmmaf <- log(648); label("cys-mcMMAF clearance CLMMAF (L/day)") # Table S4 TVCLMMAF 648
    lvc_mcmmaf <- log(30.4); label("cys-mcMMAF central volume V3 (L)") # Table S4 TVV3 30.4
    lk12_mcmmaf <- log(74.9); label("cys-mcMMAF central-to-peripheral rate constant K34 (1/day)") # Table S4 K34 74.9
    lk21_mcmmaf <- log(5.62); label("cys-mcMMAF peripheral-to-central rate constant K43 (1/day)") # Table S4 K43 5.62
    lkdar <- log(0.067); label("First-order decline rate of the drug-to-antibody ratio with time after dose, RATE-DAR (1/day)") # Table S4 RATE-DAR 0.067; Results 'time to 50% reduction in DAR was 10.3 days'
    # Product DAR0 is read from a dataset column (control stream DAR0 = DAR)
    # and no value is printed anywhere in the paper or supplement. Back-solved
    # by the maintainers from the paper's typical-patient cys-mcMMAF Cavg of
    # 217 pg/mL (Results, 'Typical patient and impact of covariates on
    # exposure measures'); the model output is linear in DAR0.
    mdar <- fixed(3.93); label("Drug product drug-to-antibody ratio at the time of dosing, DAR0 (unitless)") # back-solved from Results cys-mcMMAF Cavg 217 pg/mL; not printed in the source
    e_sbcma_vc_mcmmaf <- 0.131; label("Power exponent of sBCMA on the cys-mcMMAF central volume (unitless)") # Table S4 theta_V3_BSBCMA 0.131
    e_igg_vc_mcmmaf <- 0.154; label("Power exponent of IgG on the cys-mcMMAF central volume (unitless)") # Table S4 theta_V3_BIGG 0.154
    e_study_dreamm1_vc_mcmmaf <- 0.507; label("Fractional shift in the cys-mcMMAF central volume for DREAMM-1 (unitless)") # Table S4 'Effect of DREAMM1/BMA117159 on V3' 0.507, applied as (1 + theta)
    # DREAMM-1 cys-mcMMAF was measured with the VALA assay; Methods equation 1
    # maps the PKC-scale prediction onto it, Conc_VALA = 1.486 * Conc_PKC -
    # 8.360 pg/mL. The control stream applies it in nM: SLP = 1.486 and
    # INT = (8.360/1000)/1.045 = 0.008 nM. Both were derived outside the fit.
    cal_slope_vala <- fixed(1.486); label("DREAMM-1 VALA-assay recalibration slope for cys-mcMMAF (unitless)") # Methods equation 1; control stream SLP = 1.486
    cal_int_vala <- fixed(-0.008); label("DREAMM-1 VALA-assay recalibration intercept for cys-mcMMAF (nM)") # Methods equation 1 '-8.360' pg/mL; control stream INT = 8.360/1000/1.045 nM, subtracted

    # ------------------------------------------------------------
    # IIV. Table 1 prints the ADC IIV as percent CV; the log-normal etas use
    # omega^2 = log(CV^2 + 1), the inverse of the Table S3 / S4 footnote
    # '%CV = SQRT(EXP(eta) - 1) * 100'. IMAX enters additively on the linear
    # scale (control stream IMAX = TVIMAX + ETA(5)), so its CV is SD / |theta|
    # and omega^2 = (1.23 * 0.327)^2.
    # ------------------------------------------------------------
    etalcl ~ 0.11556 # Table 1 ETACL 35.00 percent
    etalvc ~ 0.021092 # Table 1 ETAV1 14.60 percent
    etalq ~ 0.13355 # Table 1 ETAQ 37.80 percent
    etalvp ~ 0.19805 # Table 1 ETAV2 46.80 percent
    etacl_time_max ~ 0.16177 # Table 1 ETAIMAX 123.00 percent, additive eta: (1.23 * 0.327)^2
    etalcl_t50 ~ 0.086729 # Table 1 ETATI50 30.10 percent

    etalkdec ~ 0.14023 # Table S3 ETA(KDEC) 38.8 percent
    etalcal_slope_tab ~ 0.0057292 # Table S3 ETA(CORR) 7.58 percent

    # Table S4 reports the BLOCK(3) variances and the CL-V3 and V3-K43
    # covariances; the CL-K43 covariance is not reported and is set to 0.
    etalcl_mcmmaf + etalvc_mcmmaf + etalk21_mcmmaf ~ c(
      0.100,
      0.112, 0.260,
      0, 0.257, 0.956
    ) # Table S4 OM1E 0.100, COV(CLMMAF,V3) 0.112, OM2E 0.260, COV(V3,K43) 0.257, OM3E 0.956

    # IOV on cys-mcMMAF CL and V3, one eta pair per occasion (OCC = 1, 2, 3,
    # > 3) with a common variance (control stream $OMEGA BLOCK(2) SAME). The
    # within-occasion CL-V3 covariance is not reported and is set to 0.
    etaiov_cl_mcmmaf_1 ~ 0.100 # Table S4 'IOV, CLMMAF' 0.100
    etaiov_cl_mcmmaf_2 ~ fixed(0.100) # same variance as occasion 1
    etaiov_cl_mcmmaf_3 ~ fixed(0.100) # same variance as occasion 1
    etaiov_cl_mcmmaf_4 ~ fixed(0.100) # same variance as occasion 1
    etaiov_vc_mcmmaf_1 ~ 0.185 # Table S4 'IOV, V3' 0.185
    etaiov_vc_mcmmaf_2 ~ fixed(0.185) # same variance as occasion 1
    etaiov_vc_mcmmaf_3 ~ fixed(0.185) # same variance as occasion 1
    etaiov_vc_mcmmaf_4 ~ fixed(0.185) # same variance as occasion 1

    # ------------------------------------------------------------
    # Residual error. ADC and total antibody: additive on the log scale,
    # with the printed SIGMA a variance. cys-mcMMAF: combined proportional +
    # additive with separate SIGMAs per study.
    # ------------------------------------------------------------
    expSd <- 0.13784; label("ADC residual SD on the log scale (unitless)") # Table 1 'RES ERR, additive SIGMA on log scale' 0.019 (variance)
    expSd_tab <- 0.22716; label("Total-antibody residual SD on the log scale (unitless)") # Table S3 'RES ERR, additive SIGMA on log scale' 0.0516 (variance)
    propSd_mcmmaf <- 0.22583; label("cys-mcMMAF proportional residual SD, DREAMM-2 (fraction)") # Table S4 DREAMM-2 SIGMA1 0.051 (variance)
    addSd_mcmmaf <- 0.014142; label("cys-mcMMAF additive residual SD, DREAMM-2 (nM)") # Table S4 DREAMM-2 SIGMA2 0.0002 (variance)
    propSd_mcmmaf_dreamm1 <- 0.29172; label("cys-mcMMAF proportional residual SD, DREAMM-1 (fraction)") # Table S4 DREAMM-1 SIGMA3 0.0851 (variance)
    addSd_mcmmaf_dreamm1 <- 0.030; label("cys-mcMMAF additive residual SD, DREAMM-1 (nM)") # Table S4 DREAMM-1 SIGMA4 0.0009 (variance)
  })
  model({
    # Molecular weights in the control streams: F1 = 1/152.1 converts the mg
    # ADC dose to umol (152.1 kg/mol), and cys-mcMMAF ng/mL is divided by
    # 1.045 to give nM (1.045 kg/mol). 1 mg ADC = 1000/152.1 nmol.
    mw_adc <- 152.1
    mw_mcmmaf <- 1.045

    # ------------------------------------------------------------
    # ADC individual parameters (control stream $PK). Covariate factors on
    # the categorical effects are theta^indicator.
    # ------------------------------------------------------------
    cl_base <- exp(lcl + etalcl) *
      (ALB / 40)^e_alb_cl *
      (SBCMA / 100)^e_sbcma_cl *
      (IGG / 10)^e_igg_cl *
      (WT / 75)^e_wt_cl *
      e_dose_low_cl^DOSE_LOW *
      e_study_dreamm1_cl^STUDY_DREAMM1
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_t50 <- exp(lcl_t50 + etalcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    # Sigmoidal time-varying CL, CLT = exp(IMAX * TIME^GAMMA / (TI50^GAMMA +
    # TIME^GAMMA)), with TIME the time since the first dose; rxode2's t
    # matches it when the first dose is at t = 0.
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    vc <- exp(lvc + etalvc) *
      (WT / 75)^e_wt_vc *
      (ALB / 40)^e_alb_vc *
      e_sexm_vc^(1 - SEXF) *
      e_study_dreamm1_vc^STUDY_DREAMM1
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * e_dose_low_vp^DOSE_LOW

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------
    # Naked antibody and total-antibody correction (Table S3).
    # ------------------------------------------------------------
    kdec <- exp(lkdec + etalkdec) * e_dose_low_kdec^DOSE_LOW
    cal_slope_tab <- exp(lcal_slope_tab + etalcal_slope_tab) *
      e_study_dreamm1_cal_slope_tab^STUDY_DREAMM1 *
      e_dose_low_cal_slope_tab^DOSE_LOW
    # The total-mAb control stream sets K10 = CL/V1 - KDEC on the ADC and
    # K30 = K10 on the naked antibody: the ADC loses CL/V1 in total (kdec of
    # it to naked antibody), and the naked antibody is eliminated at
    # CL/V1 - KDEC.
    kel_dar0 <- kel - kdec

    # ------------------------------------------------------------
    # cys-mcMMAF individual parameters (Table S4) with IOV by occasion.
    # ------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC > 3)
    iov_cl_mcmmaf <- oc1 * etaiov_cl_mcmmaf_1 + oc2 * etaiov_cl_mcmmaf_2 + oc3 * etaiov_cl_mcmmaf_3 + oc4 * etaiov_cl_mcmmaf_4
    iov_vc_mcmmaf <- oc1 * etaiov_vc_mcmmaf_1 + oc2 * etaiov_vc_mcmmaf_2 + oc3 * etaiov_vc_mcmmaf_3 + oc4 * etaiov_vc_mcmmaf_4
    cl_mcmmaf <- exp(lcl_mcmmaf + etalcl_mcmmaf + iov_cl_mcmmaf)
    vc_mcmmaf <- exp(lvc_mcmmaf + etalvc_mcmmaf + iov_vc_mcmmaf) *
      (SBCMA / 100)^e_sbcma_vc_mcmmaf *
      (IGG / 10)^e_igg_vc_mcmmaf *
      (1 + e_study_dreamm1_vc_mcmmaf * STUDY_DREAMM1)
    k12_mcmmaf <- exp(lk12_mcmmaf)
    k21_mcmmaf <- exp(lk21_mcmmaf + etalk21_mcmmaf)
    kel_mcmmaf <- cl_mcmmaf / vc_mcmmaf

    # Empirical drug-to-antibody ratio, reset to DAR0 at every dose and
    # declining exponentially with the time after the most recent dose
    # (Results equation; control stream DARF = DAR0 * EXP(-BETA * TAD)).
    kdar <- exp(lkdar)
    tad_dar <- tad()
    dar <- mdar * exp(-kdar * tad_dar)

    # ------------------------------------------------------------
    # ODEs. Dose (mg ADC) goes into central.
    # ------------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_dar0) <- kdec * central - kel_dar0 * central_dar0 - k12 * central_dar0 + k21 * peripheral1_dar0
    d/dt(peripheral1_dar0) <- k12 * central_dar0 - k21 * peripheral1_dar0
    # cys-mcMMAF input = ADC clearance flux (umol ADC/day in the control
    # stream) times the current DAR; here in nmol/day.
    d/dt(central_mcmmaf) <- kel * central * (1000 / mw_adc) * dar - kel_mcmmaf * central_mcmmaf - k12_mcmmaf * central_mcmmaf + k21_mcmmaf * peripheral1_mcmmaf
    d/dt(peripheral1_mcmmaf) <- k12_mcmmaf * central_mcmmaf - k21_mcmmaf * peripheral1_mcmmaf

    # ------------------------------------------------------------
    # Observations.
    # ------------------------------------------------------------
    # ADC, mg/L = ug/mL.
    Cc <- central / vc
    # Total antibody in nM: (ADC + naked antibody) * CORR.
    Cc_tab <- (central + central_dar0) / vc * (1000 / mw_adc) * cal_slope_tab
    # cys-mcMMAF in nM on the PKC-assay scale, mapped onto the VALA scale for
    # DREAMM-1 records.
    Cc_mcmmaf_pkc <- central_mcmmaf / vc_mcmmaf
    Cc_mcmmaf <- (1 - STUDY_DREAMM1) * Cc_mcmmaf_pkc + STUDY_DREAMM1 * (cal_slope_vala * Cc_mcmmaf_pkc + cal_int_vala)

    propSd_mcmmaf_i <- propSd_mcmmaf * (1 - STUDY_DREAMM1) + propSd_mcmmaf_dreamm1 * STUDY_DREAMM1
    addSd_mcmmaf_i <- addSd_mcmmaf * (1 - STUDY_DREAMM1) + addSd_mcmmaf_dreamm1 * STUDY_DREAMM1

    Cc ~ lnorm(expSd)
    Cc_tab ~ lnorm(expSd_tab)
    Cc_mcmmaf ~ add(addSd_mcmmaf_i) + prop(propSd_mcmmaf_i)
  })
}
