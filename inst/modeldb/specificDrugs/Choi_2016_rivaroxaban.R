Choi_2016_rivaroxaban <- function() {
  description <- paste(
    "Population PK/PD model for rivaroxaban 20 mg once daily taken with food",
    "in healthy male volunteers, fit as the active comparator arm of the",
    "GCC-4401C single-and-multiple-ascending-dose phase I study (Choi 2016).",
    "Disposition is two-compartment with linear clearance; absorption is a",
    "Weibull-type process (shape gamma) combined with zero-order release into",
    "the depot over a duration D1, with no absorption lag and no retained",
    "covariates. Six pharmacodynamic markers are carried as direct-effect",
    "(no-delay) outputs driven by the plasma concentration: coagulation",
    "factor X activity, factor X chromogenic activity, anti-factor Xa",
    "activity, prothrombin time in INR and in seconds, and activated partial",
    "thromboplastin time. The PD layer was fit sequentially on individual",
    "Bayesian PK estimates. The baseline aPTT is not reported anywhere in the",
    "paper or its supplement, so that output is the drug-induced CHANGE from",
    "baseline and is named with a _chg suffix; see the vignette Errata.",
    "Companion GCC-4401C model: modellib('Choi_2016_gcc4401c').",
    sep = " "
  )
  reference <- paste(
    "Choi HY, Choi S, Kim YH, Lim HS.",
    "Population pharmacokinetic and pharmacodynamic modeling analysis of",
    "GCC-4401C, a novel direct factor Xa inhibitor, in healthy volunteers.",
    "CPT Pharmacometrics Syst Pharmacol. 2016 Oct;5(10):532-543.",
    "doi:10.1002/psp4.12103.",
    "No NONMEM control stream for the rivaroxaban arm is deposited; the",
    "absorption structure here is reconstructed from the GCC-4401C stream",
    "(Supplementary Data PSP4-5-532-s008) and adjudicated against Figure 2c",
    "-- see the vignette Errata.",
    sep = " "
  )
  vignette <- "Choi_2016_gcc4401c_rivaroxaban"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    OCC = list(
      description = "Integer occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "The rivaroxaban arm received 20 mg on day 1 and once daily on days",
        "3-9, so the same seven occasions as the GCC-4401C arm apply (Choi",
        "2016 Methods; the GCC-4401C control stream PSP4-5-532-s008 carries",
        "seven $OMEGA BLOCK(1) SAME occasion slots). OCC takes values 1-7.",
        "For the two PD parameters with a separately-reported IOV",
        "(the coagulation-factor-X and chromogenic-factor-X EC50), the PD",
        "control stream PSP4-5-532-s009 groups the same column into three",
        "slots -- OCC < 6, OCC == 6 and OCC == 7 -- and that grouping is",
        "reproduced here."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight (kg)",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened and statistically significant on Vc, but NOT retained:",
        "'the final model was chosen as a base model without any covariate",
        "because the PK model for rivaroxaban containing WT as a covariate",
        "was not successful, with frequent rounding errors, which is likely",
        "because of less data for rivaroxaban in the current study' (Choi",
        "2016 Discussion). Retained on Vc in the companion GCC-4401C model."
      )
    ),
    AGE = list(
      description = "Age (years)",
      units = "years",
      type = "continuous",
      notes = "Screened as a candidate covariate (Choi 2016 Methods) and not retained."
    ),
    HT = list(
      description = "Height (cm)",
      units = "cm",
      type = "continuous",
      notes = "Screened as a candidate covariate (Choi 2016 Methods) and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    n_observations = 168L,
    age_range = "mean (SD) 30.5 (6.5) years for the whole S&MAD cohort",
    weight_range = "mean (SD) 77.6 (9.6) kg for the whole S&MAD cohort",
    height_range = "mean (SD) 178.7 (6.8) cm for the whole S&MAD cohort",
    sex_female_pct = 0,
    race_ethnicity = "Not reported separately for the rivaroxaban sub-group; the whole S&MAD cohort (n = 46) was 16 White, 21 African American, 5 Asian, 4 Other (Choi 2016 Table 1).",
    disease_state = "Healthy male volunteers.",
    dose_range = "Rivaroxaban 20 mg orally, 30 minutes after a standard breakfast, on day 1 and once daily on days 3-9.",
    regions = "United States (ClinicalTrials.gov NCT01954238).",
    notes = paste(
      "Six additional subjects in the 20 mg group of the GCC-4401C S&MAD",
      "study received rivaroxaban as an active comparator, contributing 168",
      "plasma rivaroxaban concentrations. Rivaroxaban was given fed, as",
      "recommended in its label, whereas GCC-4401C was given fasted -- the",
      "food difference is a deliberate part of the comparison and must be",
      "kept in mind when comparing the two models. Estimation used NONMEM",
      "7.2 ADVAN6 with FOCE-I. Baseline demographics from Choi 2016 Table 1."
    )
  )

  ini({
    # =====================================================================
    # PHARMACOKINETICS -- Choi 2016 Table 2(b) (plasma rivaroxaban, S&MAD
    # study). "The absorption process was described by a Weibull model with
    # mixed first and zero order absorption" (Results). The depot is filled
    # by a zero-order input of duration D1 and emptied by the Weibull
    # hazard h(t) = ka * gamma * (ka * t)^(gamma - 1), which is the
    # gamma = 1 generalisation of the first-order depot in the GCC-4401C
    # control stream (PSP4-5-532-s008) and matches the form used elsewhere
    # in the library (see Lindauer_2017_lacosamide_seizure.R).
    #
    # UNITS OF Ka. Table 2(b) prints "Ka, 1/h  2.24", but 2.24 behaves as
    # the Weibull SCALE in hours, not as a rate: the scale reading
    # reproduces Figure 2c (typical Tmax 3.6 h, Cmax 192 ng/mL at 20 mg
    # versus roughly 3 h and 185 ng/mL read off the figure), whereas the
    # rate reading gives Tmax 1.5 h and Cmax 280 ng/mL, which is
    # incompatible with it. The equivalent Weibull rate constant carried
    # here is therefore 1/2.24 = 0.446 /h. See the vignette Errata.
    # =====================================================================
    lka <- log(1 / 2.24)
    label("Weibull absorption rate constant, the reciprocal of the printed scale (1/h)") # Table 2(b) Ka = 2.24 printed as 1/h but behaving as the Weibull scale in h (RSE 79.9%; 95% CI -1.27 to 5.75); adjudicated against Figure 2c
    lgamma_abs <- log(1.34)
    label("Shape parameter of the Weibull absorption, gamma (unitless)") # Table 2(b) gamma = 1.34 (RSE 25.2%; 95% CI 0.68-2.00)
    ld1 <- log(1.07)
    label("Duration of zero-order release into the depot, D1 (h)") # Table 2(b) D1 = 1.07 (RSE 40.9%; 95% CI 0.21-1.93)
    lvc <- log(58.4)
    label("Central volume of distribution, Vc (L)") # Table 2(b) Vc = 58.4 (RSE 9.1%; 95% CI 48.0-68.8)
    lvp <- log(24.4)
    label("Peripheral volume of distribution, Vp (L)") # Table 2(b) Vp = 24.4 (RSE 14.2%; 95% CI 17.6-31.2)
    lq <- log(3.49)
    label("Inter-compartmental clearance, Q (L/h)") # Table 2(b) Q = 3.49 (RSE 27.8%; 95% CI 1.59-5.39)
    lcl <- log(10.1)
    label("Apparent total clearance, CL (L/h)") # Table 2(b) CL = 10.1 (RSE 5.5%; 95% CI 9.0-11.2)

    # ---- PK inter-individual variability (Table 2(b)). The eta on the
    # ---- Weibull parameter is carried with a plus sign on the reciprocal
    # ---- scale; because it is an independent zero-mean normal, the induced
    # ---- distribution of the absorption parameter is unchanged by the sign.
    etalka ~ 3.3 # Table 2(b) 'IIV Ka' = 3.3 (CV 503.1%; RSE 45.3%)
    etalvc ~ 0.01 # Table 2(b) 'IIV Vc' = 0.01 (CV 10.0%; RSE 104.1%)
    etalvp ~ 0.06 # Table 2(b) 'IIV Vp' = 0.06 (CV 25.7%; RSE 37.1%)
    etalcl ~ 0.015 # Table 2(b) 'IIV CL' = 0.015 (CV 12.3%; RSE 41.7%)

    # ---- Combined IIV + IOV on D1 over seven occasions. "IOV and IIV could
    # ---- not be separately estimated, and the lumped variances of IOV and
    # ---- IIV for D1 were estimated" (Results). Occasions 2-7 are fixed to
    # ---- the occasion-1 variance to encode $OMEGA BLOCK(1) SAME.
    etaiov_d1_1 ~ 0.13 # Table 2(b) 'IIV + IOV D1' = 0.13 (CV 144.8%; RSE 48.6%)
    etaiov_d1_2 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME
    etaiov_d1_3 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME
    etaiov_d1_4 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME
    etaiov_d1_5 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME
    etaiov_d1_6 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME
    etaiov_d1_7 ~ fixed(0.13) # $OMEGA BLOCK(1) SAME

    propSd <- 0.28
    label("Proportional residual error on plasma concentration (fraction)") # Table 2(b) 'e (proportional)' = 0.28 (SD; RSE 5.4%)

    # =====================================================================
    # PHARMACODYNAMICS -- Choi 2016 Table 3, rivaroxaban rows. Direct-effect
    # models throughout (no hysteresis, no delay). Baselines are read off
    # Supplementary Figure S4, which plots GCC-4401C and rivaroxaban on the
    # same axes; where the figure does not plot the marker the output is the
    # CHANGE from baseline and carries a _chg suffix.
    # =====================================================================

    # ---- (a) Coagulation factor X activity, % (sigmoid Emax)
    lrbase_cfx <- fixed(log(100))
    label("Typical baseline coagulation factor X activity (%)") # Digitised from Supplementary Figure S4 (Coagulation Factor X panels): pre-dose median 100%; not estimated by the paper
    lemax_cfx <- log(42.9)
    label("Maximum decrease in coagulation factor X activity, Emax (%)") # Table 3(a) rivaroxaban Emax = 42.9 (RSE 19.0%; 95% CI 26.9-58.9)
    lec50_cfx <- log(194.0)
    label("Rivaroxaban concentration at half-maximal factor X inhibition, EC50 (ng/mL)") # Table 3(a) rivaroxaban EC50 = 194.0 (RSE 27.4%; 95% CI 99.9-288.1)
    lhill_cfx <- log(0.90)
    label("Sigmoidicity of the factor X model, gamma (unitless)") # Table 3(a) rivaroxaban gamma = 0.90 (RSE 21.3%; 95% CI 0.53-1.28)
    etalec50_cfx ~ 0.33 # Table 3(a) rivaroxaban 'IIV EC50' = 0.33 (CV 62.1%; RSE 117.8%)
    etaiov_ec50_cfx_1 ~ 0.25 # Table 3(a) rivaroxaban 'IOV EC50' = 0.25 (CV 52.9%; RSE 155.5%)
    etaiov_ec50_cfx_2 ~ fixed(0.25) # $OMEGA BLOCK(1) SAME
    etaiov_ec50_cfx_3 ~ fixed(0.25) # $OMEGA BLOCK(1) SAME
    propSd_cfx <- 0.06
    label("Proportional residual error on coagulation factor X activity (fraction)") # Table 3(a) rivaroxaban 'e (proportional)' = 0.06 (SD; RSE 6.9%); no additive term reported

    # ---- (b) Factor X chromogenic activity assay, % (simple Emax)
    lrbase_fxcaa <- fixed(log(100))
    label("Typical baseline factor X chromogenic activity (%)") # Digitised from Supplementary Figure S4 (Factor X Chromogenic Activity Assay panels): pre-dose median 100%; not estimated by the paper
    lemax_fxcaa <- log(112.0)
    label("Maximum decrease in factor X chromogenic activity, Emax (%)") # Table 3(b) rivaroxaban Emax = 112.0 (RSE 3.2%; 95% CI 104.9-119.1)
    lec50_fxcaa <- log(126.0)
    label("Rivaroxaban concentration at half-maximal chromogenic factor X inhibition, EC50 (ng/mL)") # Table 3(b) rivaroxaban EC50 = 126.0 (RSE 13.1%; 95% CI 93.7-158.3)
    etalec50_fxcaa ~ 0.07 # Table 3(b) rivaroxaban 'IIV EC50' = 0.07 (CV 26.9%; RSE 62.9%)
    etaiov_ec50_fxcaa_1 ~ 0.01 # Table 3(b) rivaroxaban 'IOV EC50' = 0.01 (CV 8.2%; RSE 90.1%)
    etaiov_ec50_fxcaa_2 ~ fixed(0.01) # $OMEGA BLOCK(1) SAME
    etaiov_ec50_fxcaa_3 ~ fixed(0.01) # $OMEGA BLOCK(1) SAME
    addSd_fxcaa <- 3.89
    label("Additive residual error on factor X chromogenic activity (%)") # Table 3(b) rivaroxaban 'e (additive), %' = 3.89 (RSE 15.0%)
    propSd_fxcaa <- 0.05
    label("Proportional residual error on factor X chromogenic activity (fraction)") # Table 3(b) rivaroxaban 'e (proportional)' = 0.05 (SD; RSE 33.2%)

    # ---- (c) Anti-factor Xa activity, IU/mL (linear)
    lrbase_afx <- fixed(log(0.05))
    label("Typical baseline anti-factor Xa activity (IU/mL)") # Digitised from Supplementary Figure S4 (Anti-Factor Xa panels): pre-dose median sits on the axis, ~0.05 IU/mL; not estimated by the paper
    lslope_afx <- log(0.005)
    label("Linear slope of anti-factor Xa activity on plasma concentration (IU/mL per ng/mL)") # Table 3(c) rivaroxaban SLOPE = 0.005 (RSE 6.4%; 95% CI 0.004-0.006)
    etalslope_afx ~ 0.01 # Table 3(c) rivaroxaban 'IIV SLOPE' = 0.01 (CV 10.8%; RSE 43.7%)
    addSd_afx <- 0.04
    label("Additive residual error on anti-factor Xa activity (IU/mL)") # Table 3(c) rivaroxaban 'e (additive), IU/mL' = 0.04 (RSE 11.7%)
    propSd_afx <- 0.29
    label("Proportional residual error on anti-factor Xa activity (fraction)") # Table 3(c) rivaroxaban 'e (proportional)' = 0.29 (SD; RSE 13.0%)

    # ---- (d) Prothrombin time, INR (simple Emax)
    lrbase_ptinr <- fixed(log(1.05))
    label("Typical baseline prothrombin time (INR)") # Digitised from Supplementary Figure S4 (PT (INR) panels): pre-dose median 1.05; not estimated by the paper
    lemax_ptinr <- log(0.71)
    label("Maximum increase in prothrombin time, Emax (INR)") # Table 3(d) rivaroxaban Emax = 0.71 (RSE 26.5%); the printed 95% CI '104.9-119.1' is a duplicate of the Table 3(b) row and cannot belong to this parameter
    lec50_ptinr <- log(434.0)
    label("Rivaroxaban concentration at half-maximal PT (INR) effect, EC50 (ng/mL)") # Table 3(d) rivaroxaban EC50 = 434.0 (RSE 37.3%; 95% CI 116.5-751.5)
    etalec50_ptinr ~ 0.06 # Table 3(d) rivaroxaban 'IIV + IOV EC50' = 0.06 (CV 25.8%; RSE 53.9%); lumped, carried at the subject level
    addSd_ptinr <- 0.04
    label("Additive residual error on prothrombin time in INR (INR)") # Table 3(d) rivaroxaban 'e (additive), INR' = 0.04 (RSE 6.2%)

    # ---- (e) Prothrombin time, seconds (simple Emax)
    lrbase_ptsec <- fixed(log(12.2))
    label("Typical baseline prothrombin time (s)") # Digitised from Supplementary Figure S4 (PT (sec) panels): pre-dose median 12.2 s; not estimated by the paper
    lemax_ptsec <- log(6.86)
    label("Maximum increase in prothrombin time, Emax (s)") # Table 3(e) rivaroxaban Emax = 6.86 (RSE 39.1%; 95% CI 1.61-12.11)
    lec50_ptsec <- log(418.0)
    label("Rivaroxaban concentration at half-maximal PT (s) effect, EC50 (ng/mL)") # Table 3(e) rivaroxaban EC50 = 418.0 (RSE 60.3%; 95% CI -75.9 to 911.9)
    etalec50_ptsec ~ 0.07 # Table 3(e) rivaroxaban 'IIV + IOV EC50' = 0.07 (CV 26.6%; RSE 71.5%); lumped, carried at the subject level
    addSd_ptsec <- 0.04
    label("Additive residual error on prothrombin time in seconds (s)") # Table 3(e) rivaroxaban 'e (additive), sec' = 0.04 (RSE 9.3%)

    # ---- (f) Activated partial thromboplastin time (sigmoid Emax). No
    # ---- baseline aPTT is reported or plotted, so the output is the
    # ---- PROLONGATION.
    lemax_aptt <- log(12.4)
    label("Maximum prolongation of aPTT, Emax (s)") # Table 3(f) rivaroxaban Emax = 12.4 (RSE 10.3%; 95% CI 9.9-14.9)
    lec50_aptt <- log(135.0)
    label("Rivaroxaban concentration at half-maximal aPTT prolongation, EC50 (ng/mL)") # Table 3(f) rivaroxaban EC50 = 135.0 (RSE 14.3%; 95% CI 97.2-172.8)
    lhill_aptt <- log(0.94)
    label("Sigmoidicity of the aPTT model, gamma (unitless)") # Table 3(f) rivaroxaban gamma = 0.94 (RSE 18.1%; 95% CI 0.61-1.28)
    etalhill_aptt ~ 0.27 # Table 3(f) rivaroxaban 'IIV + IOV gamma' = 0.27 (CV 55.7%; RSE 38.5%); lumped, carried at the subject level
    addSd_aptt_chg <- 0.04
    label("Additive residual error on the aPTT prolongation (s)") # Table 3(f) rivaroxaban 'e (additive), sec' = 0.04 (RSE 9.3%)
  })

  model({
    # =====================================================================
    # 1. Occasion indicators (seven PK occasions; three PD slots)
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

    # =====================================================================
    # 2. Individual PK parameters
    # =====================================================================
    ka <- exp(lka + etalka)
    gamma_abs <- exp(lgamma_abs)
    d1 <- exp(ld1 + iov_d1)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    cl <- exp(lcl + etalcl)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # =====================================================================
    # 3. Weibull absorption. The hazard is evaluated at time after the most
    #    recent dose, so it restarts with each dose; the epsilon keeps the
    #    power finite if the model is re-fit with gamma < 1.
    # =====================================================================
    tad_depot <- tad(depot)
    haz_abs <- ka * gamma_abs * (ka * tad_depot + 1e-10)^(gamma_abs - 1)

    d/dt(depot) <- -haz_abs * depot
    d/dt(central) <- haz_abs * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    dur(depot) <- d1 # zero-order release; dose records need rate = -2

    Cc <- 1000 * central / vc # ng/mL (dose in mg, vc in L)

    # =====================================================================
    # 4. Direct-effect pharmacodynamics (Table 3, rivaroxaban rows)
    # =====================================================================
    # (a) Coagulation factor X activity (%)
    emax_cfx <- exp(lemax_cfx)
    ec50_cfx <- exp(lec50_cfx + etalec50_cfx + poc1 * etaiov_ec50_cfx_1 +
      poc2 * etaiov_ec50_cfx_2 + poc3 * etaiov_ec50_cfx_3)
    hill_cfx <- exp(lhill_cfx)
    cfx <- exp(lrbase_cfx) -
      emax_cfx * Cc^hill_cfx / (ec50_cfx^hill_cfx + Cc^hill_cfx)

    # (b) Factor X chromogenic activity assay (%). Emax (112%) exceeds the
    #     100% baseline, so the fitted curve can fall below zero at
    #     concentrations far above those studied -- see the vignette Errata.
    emax_fxcaa <- exp(lemax_fxcaa)
    ec50_fxcaa <- exp(lec50_fxcaa + etalec50_fxcaa + poc1 * etaiov_ec50_fxcaa_1 +
      poc2 * etaiov_ec50_fxcaa_2 + poc3 * etaiov_ec50_fxcaa_3)
    fxcaa <- exp(lrbase_fxcaa) - emax_fxcaa * Cc / (ec50_fxcaa + Cc)

    # (c) Anti-factor Xa activity (IU/mL), linear
    slope_afx <- exp(lslope_afx + etalslope_afx)
    afx <- exp(lrbase_afx) + slope_afx * Cc

    # (d) Prothrombin time (INR)
    emax_ptinr <- exp(lemax_ptinr)
    ec50_ptinr <- exp(lec50_ptinr + etalec50_ptinr)
    ptinr <- exp(lrbase_ptinr) + emax_ptinr * Cc / (ec50_ptinr + Cc)

    # (e) Prothrombin time (seconds)
    emax_ptsec <- exp(lemax_ptsec)
    ec50_ptsec <- exp(lec50_ptsec + etalec50_ptsec)
    ptsec <- exp(lrbase_ptsec) + emax_ptsec * Cc / (ec50_ptsec + Cc)

    # (f) Activated partial thromboplastin time: PROLONGATION (s), because
    #     no baseline aPTT is reported in the paper or its supplement.
    emax_aptt <- exp(lemax_aptt)
    ec50_aptt <- exp(lec50_aptt)
    hill_aptt <- exp(lhill_aptt + etalhill_aptt)
    aptt_chg <- emax_aptt * Cc^hill_aptt / (ec50_aptt^hill_aptt + Cc^hill_aptt)

    # =====================================================================
    # 5. Residual error
    # =====================================================================
    Cc ~ prop(propSd)
    cfx ~ prop(propSd_cfx)
    fxcaa ~ add(addSd_fxcaa) + prop(propSd_fxcaa)
    afx ~ add(addSd_afx) + prop(propSd_afx)
    ptinr ~ add(addSd_ptinr)
    ptsec ~ add(addSd_ptsec)
    aptt_chg ~ add(addSd_aptt_chg)
  })
}
