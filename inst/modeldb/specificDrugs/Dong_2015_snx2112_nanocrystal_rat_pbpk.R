Dong_2015_snx2112_nanocrystal_rat_pbpk <- function() {
  description <- paste0(
    "Preclinical (rat). PBPK (whole-body, perfusion-limited, coded in ",
    "MATLAB). Nanoparticulate model of SNX-2112 (Hsp90 inhibitor) after an ",
    "intravenous bolus of 203 nm poloxamer-188-stabilised nanocrystals ",
    "(Dong et al. 2015, Int J Nanomedicine, Eqs 9-19). Extends the ",
    "nonparticulate model (Dong_2015_snx2112_rat_pbpk, whose parameters are ",
    "fixed here) with a nanoparticulate species in plasma, liver and spleen: ",
    "nanocrystals are taken up from plasma into the liver and spleen by ",
    "apparent uptake clearances and release dissolved drug by a first-order ",
    "release constant Krel in each of the three compartments. The observed ",
    "plasma concentration is total (dissolved + nanoparticulate) drug. ",
    "Deterministic: no inter-individual or residual variance is reported. ",
    "WARNING: inherits the nonparticulate model's failure to reproduce the ",
    "source's lung and clearance data, and the printed release term ",
    "(Krel x V x C, Krel in 1/h) does not reproduce the source's own ",
    "Figure 9 simulation, in which particulate drug leaves the liver and ",
    "spleen within about 2 h; as printed, the particles stay in the spleen ",
    "for days. See the vignette before using tissue predictions."
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

  # Every ODE state is an AMOUNT of SNX-2112 in umol. A bare name holds the
  # dissolved (nonparticulate) drug and the _np suffix holds drug still in
  # nanocrystal form.
  compartmentData <- list(
    plasma = list(analyte = "SNX-2112", units = "umol", specimen = "plasma", verified = TRUE),
    plasma_np = list(analyte = "SNX-2112 nanocrystal", units = "umol", specimen = "plasma", verified = TRUE),
    heart = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    spleen_np = list(analyte = "SNX-2112 nanocrystal", units = "umol", specimen = "tissue", verified = TRUE),
    intestine = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    other = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "SNX-2112", units = "umol", specimen = "tissue", verified = TRUE),
    liver_np = list(analyte = "SNX-2112 nanocrystal", units = "umol", specimen = "tissue", verified = TRUE),
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
    dose_range = "2 mg/kg SNX-2112 intravenous bolus (tail vein) as poloxamer-188-stabilised nanocrystals (203 nm, zeta potential -11.6 mV).",
    regions = "Jinan University, Guangzhou, China",
    notes = paste0(
      "Destructive sampling: n = 5 rats at each of 10 times (10, 20 and 30 min; ",
      "1, 2, 4, 6, 8, 12 and 18 h), plasma and six tissues per rat. Only the ",
      "three nanoparticle-specific parameters (Upli, Upsp, Krel) were estimated ",
      "from this arm; every nonparticulate parameter was fixed to its ",
      "cosolvent-arm estimate (Results, 'Nanoparticulate PBPK modeling')."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Physiology -- Table 1 (250 g rat) and Methods. Identical to
    # Dong_2015_snx2112_rat_pbpk; see that file for the derivations.
    # ----------------------------------------------------------------
    qc <- fixed(2.69)
    label("Plasma cardiac output Qco (L/h)") # Methods 'Nonparticulate PBPK model': Qco = 2.69 L/h; Table 1 row 'Lung' = 2.6900
    q_liver <- fixed(0.492)
    label("Liver plasma flow (L/h); used for both Qha and Qli of Eqs 17 and 19") # Table 1 row 'Liver' = 0.4920; single printed value used for both symbols (see vignette)
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
    v_plasma <- fixed(0.008695)
    label("Plasma volume Vpl (L)") # NOT printed; Table 1 residual 0.250 - 0.221000 - 0.020305 = 0.008695 L (see Dong_2015_snx2112_rat_pbpk)
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
    # Nonparticulate drug parameters -- Table 4, all FIXED in this arm
    # (Results, 'Nanoparticulate PBPK modeling': 'Parameters associated the
    # nonparticulate form were taken as the corresponding estimates from
    # nonparticulate PBPK modeling and were fixed during model fitting').
    # ----------------------------------------------------------------
    lkp_kidney <- fixed(log(5.92))
    label("Log kidney-to-plasma partition coefficient Kp,kd (unitless)") # Table 4 row 'Kp,kd' = 5.9200
    lkp_spleen <- fixed(log(6.78))
    label("Log spleen-to-plasma partition coefficient Kp,sp (unitless)") # Table 4 row 'Kp,sp' = 6.7800
    lkp_intestine <- fixed(log(12.2))
    label("Log intestine-to-plasma partition coefficient Kp,in (unitless)") # Table 4 row 'Kp,in' = 12.2000
    lkp_heart <- fixed(log(2.24))
    label("Log heart-to-plasma partition coefficient Kp,ht (unitless)") # Table 4 row 'Kp,ht' = 2.2400
    lkp_liver <- fixed(log(4.03))
    label("Log liver-to-plasma partition coefficient Kp,li (unitless)") # Table 4 row 'Kp,li' = 4.0300 (cosolvent-arm estimate)
    lkp_lung <- fixed(log(0.0213))
    label("Log lung-to-plasma partition coefficient Kp,lu (unitless)") # Table 4 row 'Kp,lu' = 0.0213 (cosolvent-arm estimate)
    lkp_other <- fixed(log(5.15))
    label("Log remaining-tissue-to-plasma partition coefficient Kp,rm (unitless)") # Table 4 row 'Kp,rm' = 5.1500 (cosolvent-arm estimate)
    lclh <- fixed(log(0.0246))
    label("Log apparent hepatic clearance CLH acting on the liver tissue concentration (L/h)") # Table 4 row 'CLH' = 0.0246 L/h (cosolvent-arm estimate)

    # ----------------------------------------------------------------
    # Nanoparticle-specific parameters -- Table 4, estimated from the
    # nanocrystal arm (no footnote).
    # ----------------------------------------------------------------
    lclup_liver <- log(0.0906)
    label("Log apparent nanoparticle uptake clearance from plasma into liver Upli (L/h)") # Table 4 row 'Upli' = 0.0906 L/h (%CV 36.0)
    lclup_spleen <- log(0.182)
    label("Log apparent nanoparticle uptake clearance from plasma into spleen Upsp (L/h)") # Table 4 row 'Upsp' = 0.1820 L/h (%CV 26.0)
    lkrel <- log(0.0257)
    label("Log first-order in vivo release rate constant of dissolved drug from nanocrystals Krel (1/h)") # Table 4 row 'Krel' = 0.0257 1/h (%CV 24.0)

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
    clup_liver <- exp(lclup_liver)
    clup_spleen <- exp(lclup_spleen)
    krel <- exp(lkrel)

    # The single Table 1 liver flow enters Eq 17 as Qha and Eqs 17 and 19
    # as Qli.
    q_hepatic_artery <- q_liver

    c_plasma <- plasma / v_plasma
    c_plasma_np <- plasma_np / v_plasma
    c_heart <- heart / v_heart
    c_spleen <- spleen / v_spleen
    c_spleen_np <- spleen_np / v_spleen
    c_intestine <- intestine / v_intestine
    c_other <- other / v_other
    c_kidney <- kidney / v_kidney
    c_liver <- liver / v_liver
    c_liver_np <- liver_np / v_liver
    c_lung <- lung / v_lung

    # Release fluxes (umol/h) exactly as printed: Krel Vti Cti^Nano, i.e.
    # Krel times the nanoparticulate AMOUNT.
    rel_plasma <- krel * v_plasma * c_plasma_np
    rel_spleen <- krel * v_spleen * c_spleen_np
    rel_liver <- krel * v_liver * c_liver_np

    # Eq 9
    d/dt(plasma) <- qc * (c_lung / kp_lung - c_plasma) + rel_plasma
    # Eq 10: the nanocrystal bolus is dosed into this compartment.
    d/dt(plasma_np) <- -rel_plasma - clup_liver * c_plasma_np -
      clup_spleen * c_plasma_np
    # Eq 11
    d/dt(heart) <- q_heart * (c_plasma - c_heart / kp_heart)
    # Eq 12
    d/dt(spleen) <- q_spleen * (c_plasma - c_spleen / kp_spleen) + rel_spleen
    # Eq 13
    d/dt(spleen_np) <- clup_spleen * c_plasma_np - rel_spleen
    # Eq 14
    d/dt(intestine) <- q_intestine * (c_plasma - c_intestine / kp_intestine)
    # Eq 15
    d/dt(other) <- q_other * (c_plasma - c_other / kp_other)
    # Eq 16
    d/dt(kidney) <- q_kidney * (c_plasma - c_kidney / kp_kidney)
    # Eq 17
    d/dt(liver) <- q_hepatic_artery * c_plasma +
      q_intestine * c_intestine / kp_intestine +
      q_spleen * c_spleen / kp_spleen -
      q_liver * c_liver / kp_liver -
      clh * c_liver +
      rel_liver
    # Eq 18
    d/dt(liver_np) <- clup_liver * c_plasma_np - rel_liver
    # Eq 19
    d/dt(lung) <- q_liver * c_liver / kp_liver +
      q_heart * c_heart / kp_heart +
      q_other * c_other / kp_other +
      q_kidney * c_kidney / kp_kidney -
      qc * c_lung / kp_lung

    # Total-drug concentrations, the quantities the assay measures (it does
    # not distinguish dissolved from nanocrystalline drug; Figure 9 'Total
    # drug = nonparticulate + nanoparticulate drug').
    c_liver_total <- c_liver + c_liver_np
    c_spleen_total <- c_spleen + c_spleen_np

    # Observation: total plasma SNX-2112 (uM).
    Cc <- c_plasma + c_plasma_np
    Cc ~ prop(propSd)
  })
}
