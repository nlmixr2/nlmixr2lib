Dabour_2026_doxorubicin <- function() {
  description <- paste(
    "Integrated PK/PD model of doxorubicin in women with stage I-III breast",
    "cancer receiving adjuvant doxorubicin/cyclophosphamide (Dabour 2026):",
    "three-compartment doxorubicin PK with first-order formation of",
    "doxorubicinol (one-compartment), an NT-proBNP indirect-response model",
    "with five transit compartments driven linearly by doxorubicin (both",
    "carried unchanged from Liang 2016), and the GDF15 indirect-response model",
    "with four transit compartments and an Emax stimulation of production by",
    "doxorubicin, whose baseline depends on cumulative doxorubicin dose and BMI."
  )
  reference <- paste(
    "Dabour MS, Blaes AH, Jacobson PA, Kirstein MN, Zordoky BN, Cheng S.",
    "Evaluation of Growth Differentiation Factor 15 as an Early Biomarker of",
    "Anthracycline-Induced Cardiotoxicity Using an Integrated PK/PD Model.",
    "Clin Transl Sci. 2026;19(5):e70595.",
    "doi:10.1111/cts.70595.",
    "Doxorubicin / doxorubicinol PK and NT-proBNP PD layers from: Liang S,",
    "Brundage RC, Jacobson PA, Blaes A, Kirstein MN. Pharmacokinetic-",
    "pharmacodynamic modelling of acute N-terminal pro B-type natriuretic",
    "peptide after doxorubicin infusion in breast cancer. Br J Clin Pharmacol.",
    "2016;82(3):773-783. doi:10.1111/bcp.12989."
  )
  vignette <- "Dabour_2026_doxorubicin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Two biomarker transit chains share one model, so each chain and its
  # circulating pool carry the biomarker name.
  paper_specific_compartments <- c(
    "transit1_ntprobnp",
    "transit2_ntprobnp",
    "transit3_ntprobnp",
    "transit4_ntprobnp",
    "transit5_ntprobnp",
    "ntprobnp",
    "transit1_gdf15",
    "transit2_gdf15",
    "transit3_gdf15",
    "transit4_gdf15",
    "gdf15"
  )

  covariateData <- list(
    DOSE_DOXORUBICIN_CUM_MGM2 = list(
      description = "Cumulative doxorubicin dose received up to and including the current chemotherapy cycle",
      units = "mg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the GDF15 baseline, normalized to 60 mg/m^2:",
        "(DOSE_DOXORUBICIN_CUM_MGM2 / 60)^0.352 (Dabour 2026 Equation 4).",
        "In the source study each cycle delivered 60 mg/m^2, so the column",
        "takes the values 60, 120, 180 and 240 mg/m^2 for cycles 1-4 and must",
        "be strictly positive. It changes between the two study occasions of",
        "a subject and is constant within an occasion."
      ),
      source_name = "Cumulative DOX Dose"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on the GDF15 baseline, normalized to 25 kg/m^2: (BMI / 25)^0.956 (Dabour 2026 Equation 4).",
      source_name = "BMI"
    ),
    OCC = list(
      description = "Study occasion (chemotherapy cycle within the crossover study)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each subject contributed two consecutive chemotherapy cycles; OCC = 1",
        "for the first and OCC = 2 for the second study cycle. Selects the",
        "inter-occasion etas on doxorubicin CL, doxorubicinol CL/fm and V/fm,",
        "NT-proBNP kin0 and the GDF15 baseline. Any value other than 1 or 2",
        "zeros every inter-occasion eta."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    central = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    central_doxol = list(
      analyte = "doxorubicinol",
      units = "mg (divided by the unidentified fraction metabolized fm)",
      specimen = "plasma",
      verified = TRUE
    ),
    transit1_ntprobnp = list(
      analyte = "NT-proBNP production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_ntprobnp = list(
      analyte = "NT-proBNP production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3_ntprobnp = list(
      analyte = "NT-proBNP production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4_ntprobnp = list(
      analyte = "NT-proBNP production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5_ntprobnp = list(
      analyte = "NT-proBNP production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    ntprobnp = list(analyte = "NT-proBNP", units = "pg/mL", specimen = "plasma", verified = TRUE),
    transit1_gdf15 = list(
      analyte = "GDF15 production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2_gdf15 = list(
      analyte = "GDF15 production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3_gdf15 = list(
      analyte = "GDF15 production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4_gdf15 = list(
      analyte = "GDF15 production signal",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    gdf15 = list(analyte = "GDF15", units = "pg/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 17L,
    n_studies = 1L,
    age_range = "28-58 years",
    age_median = "45 years",
    weight_range = "54.5-140.4 kg",
    weight_median = "65.3 kg",
    bmi_range = "21.5-53.0 kg/m^2",
    bmi_median = "25.3 kg/m^2",
    bsa_median = "1.71 m^2 (range 1.60-2.50)",
    sex_female_pct = 100,
    disease_state = "Stage I-III breast cancer receiving adjuvant doxorubicin plus cyclophosphamide (dose-dense, every 14 days for four cycles)",
    dose_range = "Doxorubicin 60 mg/m^2 short IV infusion (intended 15 min) per cycle; two consecutive cycles sampled per subject, at cycles 1-4 of the regimen (cumulative 60-240 mg/m^2)",
    regions = "United States (University of Minnesota)",
    co_medication = "Cyclophosphamide 600 mg/m^2 each cycle; enalapril (5 then 10 mg/day for 7 days) before one of the two study cycles in a randomized crossover",
    n_observations = "271 GDF15 concentrations (Dabour 2026); 236 doxorubicin and 237 doxorubicinol concentrations (Liang 2016)",
    notes = paste(
      "Secondary analysis of the randomized crossover doxorubicin-enalapril",
      "interaction study (Blaes 2015). Demographics: Dabour 2026 Table 1",
      "(identical cohort to Liang 2016 Table 1). Samples by cycle: 1, 4, 16",
      "and 13 subject-cycles at cycles 1, 2, 3 and 4."
    )
  )

  ini({
    # --- Doxorubicin PK (Liang 2016 Table 2; carried unchanged by Dabour 2026) ---
    lcl <- fixed(log(54.2)); label("Doxorubicin clearance CL (L/h)") # Liang 2016 Table 2 CL_dox 54.2 L/h; Dabour 2026 Methods 'without changes in ... population parameters'
    lvc <- fixed(log(16.9)); label("Doxorubicin central volume Vc (L)") # Liang 2016 Table 2 V_C-dox 16.9 L
    lq <- fixed(log(66.1)); label("Doxorubicin intercompartmental clearance to peripheral1 (L/h); paper Q2") # Liang 2016 Table 2 Q_2-dox 66.1 L/h
    lvp <- fixed(log(1650)); label("Doxorubicin peripheral1 volume (L); paper V2") # Liang 2016 Table 2 V_2-dox 1650 L
    lq2 <- fixed(log(23.4)); label("Doxorubicin intercompartmental clearance to peripheral2 (L/h); paper Q3") # Liang 2016 Table 2 Q_3-dox 23.4 L/h
    lvp2 <- fixed(log(61.8)); label("Doxorubicin peripheral2 volume (L); paper V3") # Liang 2016 Table 2 V_3-dox 61.8 L

    # --- Doxorubicinol PK, apparent parameters scaled by 1/fm (Liang 2016 Table 2) ---
    lcle_doxol <- fixed(log(106)); label("Doxorubicinol elimination clearance divided by fraction metabolized, CL_doxol/fm (L/h)") # Liang 2016 Table 2 CL_doxo'ol/fm 106 L/h
    lvc_doxol <- fixed(log(1880)); label("Doxorubicinol volume divided by fraction metabolized, V_doxol/fm (L)") # Liang 2016 Table 2 V_C-doxo'ol/fm 1880 L

    # --- NT-proBNP indirect response, 5 transit compartments (Liang 2016 Table 3) ---
    lkin_ntprobnp <- fixed(log(110)); label("NT-proBNP baseline production rate kin0 (pg/mL/h)") # Liang 2016 Table 3 k_in_0 110; table unit 'ng/ml/h' read as pg/mL/h (Liang Figure 3C axis in pg/mL)
    lkout_ntprobnp <- fixed(log(1.66)); label("NT-proBNP first-order degradation rate constant kout (1/h)") # Liang 2016 Table 3 k_out 1.66, fixed to literature value
    lslope_ntprobnp <- fixed(log(115)); label("Linear slope of doxorubicin concentration on NT-proBNP production (L/mg)") # Liang 2016 Table 3 beta 115; applied to doxorubicin in mg/L (= ug/mL)
    lmtt_ntprobnp <- fixed(log(34.7)); label("NT-proBNP mean transit time MTT = 5 / ktr (h)") # Liang 2016 Table 3 MTT 34.7 h; Equation 2

    # --- GDF15 indirect response, 4 transit compartments (Dabour 2026 Table 2) ---
    lktr_gdf15 <- log(0.851); label("GDF15 transit rate constant Ktr (1/h)") # Dabour 2026 Table 2 K_tr 0.851 h-1 (RSE 0.792%)
    lrbase_gdf15 <- log(744.866); label("Baseline GDF15 R0 at cumulative dose 60 mg/m^2 and BMI 25 kg/m^2 (pg/mL)") # Dabour 2026 Table 2 R_0 744.866 pg/mL (RSE 15.7%)
    lkout_gdf15 <- log(1.012); label("GDF15 first-order degradation rate constant Kout (1/h)") # Dabour 2026 Table 2 K_out 1.012 h-1 (RSE 24.6%)
    lemax_gdf15 <- log(16.837); label("Maximal fractional stimulation of GDF15 production by doxorubicin Emax (unitless)") # Dabour 2026 Table 2 E_max 16.837 (RSE 11.5%)
    lec50_gdf15 <- log(0.0385); label("Doxorubicin concentration giving half-maximal GDF15 stimulation EC50 (mg/L = ug/mL)") # Dabour 2026 Table 2 EC50 0.0385 ug/mL (RSE 21.3%)
    e_dose_doxorubicin_cum_mgm2_rbase_gdf15 <- 0.352; label("Power exponent of cumulative doxorubicin dose (/60 mg/m^2) on GDF15 baseline (unitless)") # Dabour 2026 Table 2 beta1 0.352 (RSE 32.8%); Equation 4
    e_bmi_rbase_gdf15 <- 0.956; label("Power exponent of BMI (/25 kg/m^2) on GDF15 baseline (unitless)") # Dabour 2026 Table 2 beta2 0.956 (RSE 38.2%); Equation 4

    # --- Between-subject variability; omega^2 = log(CV^2 + 1) ---
    etalcl ~ fixed(0.006915) # Liang 2016 Table 2 BSV CL_dox 8.33% CV
    etalvc ~ fixed(0.01096) # Liang 2016 Table 2 BSV V_C-dox 10.5% CV
    etalq ~ fixed(0.02344) # Liang 2016 Table 2 BSV Q_2-dox 15.4% CV
    etalvp ~ fixed(0.06886) # Liang 2016 Table 2 BSV V_2-dox 26.7% CV
    etalcle_doxol ~ fixed(0.03365) # Liang 2016 Table 2 BSV CL_doxo'ol/fm 18.5% CV
    etalvc_doxol ~ fixed(0.1303) # Liang 2016 Table 2 BSV V_doxo'ol/fm 37.3% CV
    etalkin_ntprobnp ~ fixed(0.2113) # Liang 2016 Table 3 BSV k_in_0 48.5% CV
    etalslope_ntprobnp ~ fixed(0.4357) # Liang 2016 Table 3 BSV beta 73.9% CV
    etalmtt_ntprobnp ~ fixed(0.1402) # Liang 2016 Table 3 BSV MTT 38.8% CV
    etalrbase_gdf15 ~ 0.09587 # Dabour 2026 Table 2 IIV R_0 31.72% CV
    etalkout_gdf15 + etalec50_gdf15 ~ c(0.2141, -0.2100, 0.3056) # Dabour 2026 Table 2 IIV K_out 48.86% CV, EC50 59.79% CV, corr -0.821

    # --- Inter-occasion variability (two study cycles); shared variance per parameter ---
    etaiov_cl_1 ~ fixed(0.001489) # Liang 2016 Table 2 IOV CL_dox 3.86% CV
    etaiov_cl_2 ~ fixed(0.001489) # Liang 2016 Table 2 IOV CL_dox 3.86% CV; same variance
    etaiov_cle_doxol_1 ~ fixed(0.01406) # Liang 2016 Table 2 IOV CL_doxo'ol/fm 11.9% CV
    etaiov_cle_doxol_2 ~ fixed(0.01406) # Liang 2016 Table 2 IOV CL_doxo'ol/fm 11.9% CV; same variance
    etaiov_vc_doxol_1 ~ fixed(0.008779) # Liang 2016 Table 2 IOV V_doxo'ol/fm 9.39% CV
    etaiov_vc_doxol_2 ~ fixed(0.008779) # Liang 2016 Table 2 IOV V_doxo'ol/fm 9.39% CV; same variance
    etaiov_kin_ntprobnp_1 ~ fixed(0.07444) # Liang 2016 Table 3 IOV k_in_0 27.8% CV
    etaiov_kin_ntprobnp_2 ~ fixed(0.07444) # Liang 2016 Table 3 IOV k_in_0 27.8% CV; same variance
    etaiov_rbase_gdf15_1 ~ 0.01211 # Dabour 2026 Table 2 IOV R_0 11.04% CV
    etaiov_rbase_gdf15_2 ~ fixed(0.01211) # Dabour 2026 Table 2 IOV R_0 11.04% CV; same variance

    # --- Residual error ---
    propSd <- fixed(0.0873); label("Doxorubicin proportional residual SD (fraction)") # Liang 2016 Table 2 RUV-proportional-dox 8.73% CV
    addSd <- fixed(0.60); label("Doxorubicin additive residual SD (ng/mL)") # Liang 2016 Table 2 RUV-additive-dox 0.60 SD
    propSd_doxol <- fixed(0.149); label("Doxorubicinol proportional residual SD (fraction)") # Liang 2016 Table 2 RUV-proportional-doxo'ol 14.9% CV
    propSd_ntprobnp <- fixed(0.238); label("NT-proBNP proportional residual SD (fraction)") # Liang 2016 Table 3 RUV-BNP 23.8% CV
    propSd_gdf15 <- 0.156; label("GDF15 proportional residual SD (fraction)") # Dabour 2026 Table 2 RUV proportional 0.156
  })

  model({
    # Occasion indicators for the inter-occasion etas
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2
    iov_cle_doxol <- oc1 * etaiov_cle_doxol_1 + oc2 * etaiov_cle_doxol_2
    iov_vc_doxol <- oc1 * etaiov_vc_doxol_1 + oc2 * etaiov_vc_doxol_2
    iov_kin_ntprobnp <- oc1 * etaiov_kin_ntprobnp_1 + oc2 * etaiov_kin_ntprobnp_2
    iov_rbase_gdf15 <- oc1 * etaiov_rbase_gdf15_1 + oc2 * etaiov_rbase_gdf15_2

    # Doxorubicin and doxorubicinol PK
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    cle_doxol <- exp(lcle_doxol + etalcle_doxol + iov_cle_doxol)
    vc_doxol <- exp(lvc_doxol + etalvc_doxol + iov_vc_doxol)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_doxol <- cle_doxol / vc_doxol

    # Doxorubicin concentration driving both biomarkers (mg/L = ug/mL)
    cdox <- central / vc

    # NT-proBNP: kin = kin0 * (1 + slope * Cdox), five transit compartments
    kin_ntprobnp <- exp(lkin_ntprobnp + etalkin_ntprobnp + iov_kin_ntprobnp)
    kout_ntprobnp <- exp(lkout_ntprobnp)
    slope_ntprobnp <- exp(lslope_ntprobnp + etalslope_ntprobnp)
    mtt_ntprobnp <- exp(lmtt_ntprobnp + etalmtt_ntprobnp)
    ktr_ntprobnp <- 5 / mtt_ntprobnp

    # GDF15: kin = R0 * kout, Emax stimulation, four transit compartments
    rbase_gdf15 <- exp(lrbase_gdf15 + etalrbase_gdf15 + iov_rbase_gdf15) *
      (DOSE_DOXORUBICIN_CUM_MGM2 / 60)^e_dose_doxorubicin_cum_mgm2_rbase_gdf15 *
      (BMI / 25)^e_bmi_rbase_gdf15
    ktr_gdf15 <- exp(lktr_gdf15)
    kout_gdf15 <- exp(lkout_gdf15 + etalkout_gdf15)
    emax_gdf15 <- exp(lemax_gdf15)
    ec50_gdf15 <- exp(lec50_gdf15 + etalec50_gdf15)
    kin_gdf15 <- rbase_gdf15 * kout_gdf15

    # Steady-state initial conditions (Liang 2016; Dabour 2026 Equation 2)
    transit1_ntprobnp(0) <- kin_ntprobnp / ktr_ntprobnp
    transit2_ntprobnp(0) <- kin_ntprobnp / ktr_ntprobnp
    transit3_ntprobnp(0) <- kin_ntprobnp / ktr_ntprobnp
    transit4_ntprobnp(0) <- kin_ntprobnp / ktr_ntprobnp
    transit5_ntprobnp(0) <- kin_ntprobnp / ktr_ntprobnp
    ntprobnp(0) <- kin_ntprobnp / kout_ntprobnp
    transit1_gdf15(0) <- kin_gdf15 / ktr_gdf15
    transit2_gdf15(0) <- kin_gdf15 / ktr_gdf15
    transit3_gdf15(0) <- kin_gdf15 / ktr_gdf15
    transit4_gdf15(0) <- kin_gdf15 / ktr_gdf15
    gdf15(0) <- rbase_gdf15

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    # All doxorubicin clearance enters the fm-scaled doxorubicinol amount:
    # fm * CL * Cdox formation divided by fm, eliminated at (CL/fm) / (V/fm).
    d/dt(central_doxol) <- kel * central - kel_doxol * central_doxol

    d/dt(transit1_ntprobnp) <- kin_ntprobnp * (1 + slope_ntprobnp * cdox) -
      ktr_ntprobnp * transit1_ntprobnp
    d/dt(transit2_ntprobnp) <- ktr_ntprobnp * (transit1_ntprobnp - transit2_ntprobnp)
    d/dt(transit3_ntprobnp) <- ktr_ntprobnp * (transit2_ntprobnp - transit3_ntprobnp)
    d/dt(transit4_ntprobnp) <- ktr_ntprobnp * (transit3_ntprobnp - transit4_ntprobnp)
    d/dt(transit5_ntprobnp) <- ktr_ntprobnp * (transit4_ntprobnp - transit5_ntprobnp)
    d/dt(ntprobnp) <- ktr_ntprobnp * transit5_ntprobnp - kout_ntprobnp * ntprobnp

    d/dt(transit1_gdf15) <- kin_gdf15 * (1 + emax_gdf15 * cdox / (ec50_gdf15 + cdox)) -
      ktr_gdf15 * transit1_gdf15
    d/dt(transit2_gdf15) <- ktr_gdf15 * (transit1_gdf15 - transit2_gdf15)
    d/dt(transit3_gdf15) <- ktr_gdf15 * (transit2_gdf15 - transit3_gdf15)
    d/dt(transit4_gdf15) <- ktr_gdf15 * (transit3_gdf15 - transit4_gdf15)
    d/dt(gdf15) <- ktr_gdf15 * transit4_gdf15 - kout_gdf15 * gdf15

    # Observations: doxorubicin and doxorubicinol in ng/mL; biomarkers in pg/mL
    Cc <- cdox * 1000
    Cc_doxol <- central_doxol / vc_doxol * 1000

    Cc ~ add(addSd) + prop(propSd)
    Cc_doxol ~ prop(propSd_doxol)
    ntprobnp ~ prop(propSd_ntprobnp)
    gdf15 ~ prop(propSd_gdf15)
  })
}
