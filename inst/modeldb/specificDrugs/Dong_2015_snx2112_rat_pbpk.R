Dong_2015_snx2112_rat_pbpk <- function() {
  description <- paste0(
    "Preclinical (rat). PBPK (whole-body, perfusion-limited, 8 compartments, ",
    "coded in MATLAB). Nonparticulate model of SNX-2112 (Hsp90 inhibitor) ",
    "disposition after an intravenous bolus of a 45% propylene-glycol ",
    "cosolvent solution (Dong et al. 2015, Int J Nanomedicine, Eqs 1-8). ",
    "Plasma, heart, spleen, intestine, kidney, liver, lung and a lumped ",
    "remainder ('other'); the spleen and intestine drain into the liver and ",
    "the liver, heart, kidney and remainder drain into the lung, which returns ",
    "to plasma at the plasma cardiac output. Elimination is hepatic only, as ",
    "CLH times the liver tissue concentration. Deterministic: no ",
    "inter-individual or residual variance is reported. WARNING: the ",
    "published parameter set does not reproduce the source's own tissue data ",
    "-- the lung:plasma ratio is predicted about 200-fold below the measured ",
    "one because the fitted Kp,lu (0.0213) pins it, and systemic clearance ",
    "is 1.5- to 1.7-fold below the observed value; see the vignette before ",
    "using lung or absolute plasma predictions."
  )
  reference <- paste0(
    "Dong D, Wang X, Wang H, Zhang X, Wang Y, Wu B. Elucidating the in vivo ",
    "fate of nanocrystals using a physiologically based pharmacokinetic ",
    "model: a case study with the anticancer agent SNX-2112. Int J ",
    "Nanomedicine. 2015;10:2521-2535. doi:10.2147/IJN.S79734"
  )
  vignette <- "Dong_2015_snx2112_nanocrystals"
  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L"
  )

  # Every ODE state is an AMOUNT of SNX-2112 in umol; model() divides by the
  # Table 1 volume to obtain the concentration in uM.
  compartmentData <- list(
    plasma = list(analyte = "SNX-2112", units = "umol", specimen = "plasma", verified = TRUE),
    heart = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    intestine = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    other = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = NULL,
    weight_range = "193.0 +/- 8.9 g (mean +/- SD); the model physiology is for a 250 g rat (Table 1 caption)",
    sex_female_pct = 0,
    disease_state = "Healthy rats.",
    dose_range = "2 mg/kg SNX-2112 intravenous bolus (tail vein) in 45% propylene-glycol cosolvent.",
    regions = "Jinan University, Guangzhou, China",
    notes = paste0(
      "Destructive sampling: n = 5 rats at each of 10 times (10, 20 and 30 min; ",
      "1, 2, 4, 6, 8, 12 and 18 h), with plasma, heart, liver, spleen, lung, ",
      "kidney and intestine collected from each rat (Methods, 'Pharmacokinetic ",
      "study of SNX-2112'). The model was fitted to the arm means with 1/y^2 ",
      "weighting in MATLAB 8.0. The nanocrystal arm is a separate model, ",
      "Dong_2015_snx2112_nanocrystal_rat_pbpk."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Plasma flows (L/h) -- Table 1, 'Plasma flow (L/h)' column, for a 250 g
    # rat. Table 1 cites reference 42, a bioanalytical paper; the values are
    # those of reference 43 (Kagan 2011 AAPS J, Table I, rat), which gives
    # them as percent of the plasma cardiac output: liver 18.3, kidney 14.1,
    # spleen 1.0, intestine 14.3, heart 4.9. Fixed during model fitting
    # (Methods, 'Nonparticulate PBPK model').
    # ----------------------------------------------------------------
    qc <- fixed(2.69)
    label("Plasma cardiac output Qco (L/h)") # Methods 'Nonparticulate PBPK model': 'Plasma cardiac output (Qco) for rats was set at 2.69 L/h'; Table 1 row 'Lung' = 2.6900
    q_liver <- fixed(0.492)
    label("Liver plasma flow (L/h); used for both Qha and Qli of Eqs 7-8")
    # Table 1 row 'Liver' = 0.4920. Eqs 7 and 8 use two symbols, Qha (hepatic artery, Figure 1 legend) and Qli (liver outflow), but Table 1 prints one liver flow; the single printed value is used for both. See the vignette 'Hepatic flow' section.
    q_kidney <- fixed(0.379)
    label("Kidney plasma flow (L/h)") # Table 1 row 'Kidney' = 0.3790
    q_spleen <- fixed(0.0269)
    label("Spleen plasma flow (L/h)") # Table 1 row 'Spleen' = 0.0269
    q_intestine <- fixed(0.385)
    label("Intestine plasma flow (L/h)") # Table 1 row 'Intestine' = 0.3850
    q_heart <- fixed(0.132)
    label("Heart plasma flow (L/h)") # Table 1 row 'Heart' = 0.1320
    q_other <- fixed(1.28)
    label("Remaining-tissue plasma flow (L/h)") # Table 1 row 'Remaining' = 1.2800

    # ----------------------------------------------------------------
    # Volumes (L) -- Table 1, 'Organ volume (L)' column, 250 g rat
    # (Kagan 2011 Table I percent of body weight x 0.250 kg, density 1).
    # ----------------------------------------------------------------
    v_plasma <- fixed(0.008695)
    label("Plasma volume Vpl (L)")
    # NOT printed in Dong 2015 or in Kagan 2011. Recovered from Table 1, whose 'Remaining' row is the residual of a 250 g (0.250 L) rat: 0.250 - 0.221000 (remaining) - 0.020305 (liver + kidney + spleen + intestine + lung + heart) = 0.008695 L (34.8 mL/kg). Figure 9A's nanoparticulate plasma curve starts at Dose/Vpl = 1.0765 umol / 0.008695 L = 124 uM, consistent with this value.
    v_liver <- fixed(0.00915)
    label("Liver volume (L)") # Table 1 row 'Liver' = 0.009150
    v_kidney <- fixed(0.00183)
    label("Kidney volume (L)") # Table 1 row 'Kidney' = 0.001830
    v_spleen <- fixed(0.0005)
    label("Spleen volume (L)") # Table 1 row 'Spleen' = 0.000500
    v_intestine <- fixed(0.00675)
    label("Intestine volume (L)") # Table 1 row 'Intestine' = 0.006750
    v_lung <- fixed(0.00125)
    label("Lung volume (L)") # Table 1 row 'Lung' = 0.001250
    v_heart <- fixed(0.000825)
    label("Heart volume (L)") # Table 1 row 'Heart' = 0.000825
    v_other <- fixed(0.221)
    label("Remaining-tissue volume (L)") # Table 1 row 'Remaining' = 0.221000

    # ----------------------------------------------------------------
    # Tissue-to-plasma partition coefficients -- Table 4. The kidney,
    # spleen, intestine and heart values were fixed to the observed
    # tissue:plasma AUC(0-18 h) ratios (Table 4 footnote b, Figure 6); the
    # liver, lung and remainder values were estimated (footnote a).
    # ----------------------------------------------------------------
    lkp_kidney <- fixed(log(5.92))
    label("Log kidney-to-plasma partition coefficient Kp,kd (unitless)") # Table 4 row 'Kp,kd' = 5.9200 (fixed, footnote b)
    lkp_spleen <- fixed(log(6.78))
    label("Log spleen-to-plasma partition coefficient Kp,sp (unitless)") # Table 4 row 'Kp,sp' = 6.7800 (fixed, footnote b)
    lkp_intestine <- fixed(log(12.2))
    label("Log intestine-to-plasma partition coefficient Kp,in (unitless)") # Table 4 row 'Kp,in' = 12.2000 (fixed, footnote b)
    lkp_heart <- fixed(log(2.24))
    label("Log heart-to-plasma partition coefficient Kp,ht (unitless)") # Table 4 row 'Kp,ht' = 2.2400 (fixed, footnote b)
    lkp_liver <- log(4.03)
    label("Log liver-to-plasma partition coefficient Kp,li (unitless)") # Table 4 row 'Kp,li' = 4.0300 (%CV 8.1, footnote a)
    lkp_lung <- log(0.0213)
    label("Log lung-to-plasma partition coefficient Kp,lu (unitless)") # Table 4 row 'Kp,lu' = 0.0213 (%CV 9.5, footnote a)
    lkp_other <- log(5.15)
    label("Log remaining-tissue-to-plasma partition coefficient Kp,rm (unitless)") # Table 4 row 'Kp,rm' = 5.1500 (%CV 4.3, footnote a)

    # ----------------------------------------------------------------
    # Hepatic clearance -- Table 4, estimated (footnote a). Eq 7 applies it
    # to the liver TISSUE concentration Cli, not to the liver outflow
    # concentration Cli / Kp,li.
    # ----------------------------------------------------------------
    lclh <- log(0.0246)
    label("Log apparent hepatic clearance CLH acting on the liver tissue concentration (L/h)") # Table 4 row 'CLH' = 0.0246 L/h (%CV 7.6, footnote a)

    # ----------------------------------------------------------------
    # Residual error. The model was fitted by weighted least squares with
    # 1/y^2 weights (Methods, 'Data analysis'), i.e. a proportional error,
    # but its magnitude is not reported. Encoded as zero rather than
    # invented; the model is for typical-value simulation.
    # ----------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - magnitude not reported in source)")
  })

  model({
    kp_kidney <- exp(lkp_kidney)
    kp_spleen <- exp(lkp_spleen)
    kp_intestine <- exp(lkp_intestine)
    kp_heart <- exp(lkp_heart)
    kp_liver <- exp(lkp_liver)
    kp_lung <- exp(lkp_lung)
    kp_other <- exp(lkp_other)
    clh <- exp(lclh)

    # The single Table 1 liver flow enters Eq 7 as the hepatic-artery inflow
    # Qha and Eqs 7-8 as the liver outflow Qli.
    q_hepatic_artery <- q_liver

    # Concentrations (uM) = amount (umol) / volume (L).
    c_plasma <- plasma / v_plasma
    c_heart <- heart / v_heart
    c_spleen <- spleen / v_spleen
    c_intestine <- intestine / v_intestine
    c_other <- other / v_other
    c_kidney <- kidney / v_kidney
    c_liver <- liver / v_liver
    c_lung <- lung / v_lung

    # Eq 1: Vpl dCpl/dt = Qco (Clu / Kp,lu - Cpl). The intravenous bolus is
    # dosed into this compartment.
    d/dt(plasma) <- qc * (c_lung / kp_lung - c_plasma)
    # Eq 2
    d/dt(heart) <- q_heart * (c_plasma - c_heart / kp_heart)
    # Eq 3
    d/dt(spleen) <- q_spleen * (c_plasma - c_spleen / kp_spleen)
    # Eq 4
    d/dt(intestine) <- q_intestine * (c_plasma - c_intestine / kp_intestine)
    # Eq 5
    d/dt(other) <- q_other * (c_plasma - c_other / kp_other)
    # Eq 6
    d/dt(kidney) <- q_kidney * (c_plasma - c_kidney / kp_kidney)
    # Eq 7: Vli dCli/dt = Qha Cpl + Qin Cin / Kp,in + Qsp Csp / Kp,sp
    #   - Qli Cli / Kp,li - CLH Cli
    d/dt(liver) <- q_hepatic_artery * c_plasma +
      q_intestine * c_intestine / kp_intestine +
      q_spleen * c_spleen / kp_spleen -
      q_liver * c_liver / kp_liver -
      clh * c_liver
    # Eq 8: Vlu dClu/dt = Qli Cli / Kp,li + Qht Cht / Kp,ht + Qrm Crm / Kp,rm
    #   + Qkd Ckd / Kp,kd - Qco Clu / Kp,lu
    d/dt(lung) <- q_liver * c_liver / kp_liver +
      q_heart * c_heart / kp_heart +
      q_other * c_other / kp_other +
      q_kidney * c_kidney / kp_kidney -
      qc * c_lung / kp_lung

    # Observation: plasma SNX-2112 (uM), the UPLC-QTOF/MS plasma assay.
    Cc <- c_plasma
    Cc ~ prop(propSd)
  })
}
