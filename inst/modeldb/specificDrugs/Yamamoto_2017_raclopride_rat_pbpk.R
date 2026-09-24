Yamamoto_2017_raclopride_rat_pbpk <- function() {
  description <- paste(
    "PBPK (LeiCNS-PK1.0 comprehensive CNS physiologically-based model).",
    "Preclinical (rat, male Wistar). Raclopride disposition in plasma and",
    "nine CNS compartments: brain microvasculature, brain extracellular",
    "fluid, brain intracellular fluid, lysosomes, and the four",
    "cerebrospinal-fluid spaces (lateral ventricle, third + fourth ventricle,",
    "cisterna magna, subarachnoid space) draining in series back to the brain",
    "microvasculature. Transport across the blood-brain barrier and the",
    "blood-CSF barrier is the sum of a paracellular clearance (aqueous",
    "diffusivity over barrier width, on the 0.006-0.016 percent of surface",
    "area left open by tight junctions) and a transcellular clearance",
    "(transmembrane permeability on the remaining 99.8 percent), with net",
    "active transport carried by asymmetry factors and pH partitioning by",
    "pH-dependent factors. The plasma side is the empirical three-compartment",
    "model of Table 2, fitted by the authors in NONMEM 7.3; every CNS",
    "parameter is fixed to rat physiology (Table 3) or to the compound's",
    "physicochemical properties (Tables 4 and 5) and none was fitted to the",
    "CNS data, so the CNS profiles are genuine predictions."
  )
  reference <- paste(
    "Yamamoto Y, Valitalo PA, Huntjens DR, Proost JH, Vermeulen A,",
    "Krauwinkel W, Beukers MW, van den Berg DJ, Hartman R, Wong YC,",
    "Danhof M, van Hasselt JGC, de Lange ECM. Predicting Drug",
    "Concentration-Time Profiles in Multiple CNS Compartments Using a",
    "Comprehensive Physiologically-Based Pharmacokinetic Model. CPT",
    "Pharmacometrics Syst Pharmacol. 2017;6(11):765-777.",
    "doi:10.1002/psp4.12250.",
    "Model structure and all parameter values are from the main text",
    "(Tables 1-5) and the supplementary NONMEM control streams",
    "(Supplementary Material S3-S13); the asymmetry-factor and",
    "binding-factor closed forms are Supplementary Material S1."
  )
  vignette <- "Yamamoto_2017_cns_pbpk_rat"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "raclopride", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "raclopride", units = "ng", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "raclopride", units = "ng", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "raclopride", units = "ng", specimen = "plasma", verified = TRUE),
    brain_ecf = list(analyte = "raclopride", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_icf = list(analyte = "raclopride", units = "ng", specimen = "tissue", verified = TRUE),
    brain_lysosome = list(analyte = "raclopride", units = "ng", specimen = "tissue", verified = TRUE),
    brain_csf_lv = list(analyte = "raclopride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "raclopride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "raclopride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "raclopride", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 19,
    weight_range = "225-275 g (Supplementary Material S2)",
    disease_state = "healthy",
    dose_range = "0.56 mg/kg (10 min infusion), intravenous",
    notes = paste(
      "Study design from Table 1: 19 animals, 0.56 mg/kg (10 min infusion),",
      "data from Table 1 (no external reference given). Plasma and brain",
      "microdialysis (brain extracellular fluid and CSF) samples were collected",
      "after femoral-vein administration; see Supplementary Material S2 for",
      "surgery and bioanalysis. Concentrations in the CNS compartments are",
      "unbound concentrations."
    )
  )

  ini({
    # ---- Empirical plasma pharmacokinetics (Table 2) ----
    lcl <- log(46.4); label("Plasma clearance (mL/min)")  # Table 2 CL_PL
    lvc <- log(48.9); label("Central volume of distribution (mL)")  # Table 2 V_PL
    lq <- log(13.4); label("Intercompartmental clearance to peripheral 1 (mL/min)")  # Table 2 Q_PL_PER1
    lvp <- log(684); label("Peripheral 1 volume of distribution (mL)")  # Table 2 V_PER1
    lq2 <- log(69.2); label("Intercompartmental clearance to peripheral 2 (mL/min)")  # Table 2 Q_PL_PER2
    lvp2 <- log(493); label("Peripheral 2 volume of distribution (mL)")  # Table 2 V_PER2

    # ---- Interindividual variability (Table 2) ----
    # Table 2 reports these as coefficients of variation; the control
    # streams carry the corresponding variances, which are the squares of
    # the tabulated values (for example 0.239^2 = 0.0571).
    etalcl ~ 0.0208  # Table 2 CV 14.4 percent

    # ---- Residual unexplained variability (Table 2) ----
    # Only plasma carries residual error: the control streams fix EPS(2)-EPS(4),
    # the brain-ECF and CSF residuals, to zero because the CNS model predicts
    # rather than fits those data.
    propSd <- 0.141; label("Proportional residual error (fraction)")  # Table 2 sigma_plasma proportional

    # ---- Drug-specific CNS parameters (Table 4) ----
    # Transmembrane permeability and aqueous diffusivity are the in silico
    # predictions of Eqs. 2 and 1 from logP = 1.3 and molecular weight = 347;
    # the tabulated values are used directly rather than re-deriving them.
    P0 <- fixed(6.6e-4); label("Transmembrane permeability (cm/min)")  # Table 4 transmembrane permeability
    Daq <- fixed(3.1e-4); label("Aqueous diffusivity coefficient (cm^2/min)")  # Table 4 aqueous diffusivity coefficient

    # Asymmetry factors carry the net effect of active transport. They were
    # back-solved from the measured Kp,uu values with the closed forms of
    # Supplementary Material S1: influx factors are estimated and efflux
    # factors held at 1 when Kp,uu > 1, and the other way round when Kp,uu < 1.
    AFin1 <- fixed(1); label("Asymmetry factor into brain extracellular fluid (unitless)")  # Table 4 AFin1
    AFin2 <- fixed(1); label("Asymmetry factor into lateral-ventricle CSF (unitless)")  # Table 4 AFin2
    AFin3 <- fixed(1); label("Asymmetry factor into third + fourth ventricle CSF (unitless)")  # Table 4 AFin3
    AFout1 <- fixed(1.4); label("Asymmetry factor out of brain extracellular fluid (unitless)")  # Table 4 AFout1
    AFout2 <- fixed(1.1); label("Asymmetry factor out of lateral-ventricle CSF (unitless)")  # Table 4 AFout2
    AFout3 <- fixed(1.9); label("Asymmetry factor out of third + fourth ventricle CSF (unitless)")  # Table 4 AFout3

    # ---- pH-dependent factors (Table 5) ----
    # Each factor is the ratio of the uncharged fraction of the compound in
    # the donor compartment to the uncharged fraction in plasma (pH 7.4),
    # from the Henderson-Hasselbalch equations 10-17. Because brain ECF and
    # all CSF sites share pH 7.3, PHF1 to PHF4 are equal; brain ICF is pH 7.0
    # so PHF5 and PHF6 are equal; lysosomes are pH 5.0.
    PHF1 <- fixed(0.80); label("pH-dependent factor, brain extracellular fluid to microvasculature (unitless)")  # Table 5 PHF1
    PHF2 <- fixed(0.80); label("pH-dependent factor, lateral-ventricle CSF to microvasculature (unitless)")  # Table 5 PHF2
    PHF3 <- fixed(0.80); label("pH-dependent factor, third + fourth ventricle CSF to microvasculature (unitless)")  # Table 5 PHF3
    PHF4 <- fixed(0.80); label("pH-dependent factor, brain extracellular to intracellular fluid (unitless)")  # Table 5 PHF4
    PHF5 <- fixed(0.40); label("pH-dependent factor, brain intracellular to extracellular fluid (unitless)")  # Table 5 PHF5
    PHF6 <- fixed(0.40); label("pH-dependent factor, brain intracellular fluid to lysosome (unitless)")  # Table 5 PHF6
    PHF7 <- fixed(0.0041); label("pH-dependent factor, lysosome to brain intracellular fluid (unitless)")  # Table 5 PHF7

    # ---- Brain-tissue binding factor (Table 5) ----
    # Ratio of tissue-bound to unbound drug in the brain extracellular
    # space, back-solved from the total brain-to-plasma ratio Kp
    # (Supplementary Material S1). It enters only the total-brain output:
    # in the control streams it multiplies and then divides out of every
    # differential equation, so it does not affect the ODE solution.
    BF <- fixed(8.5); label("Brain-tissue binding factor (unitless)")  # Table 5 BF

    # ---- Rat CNS physiology (Table 3); identical for all ten drugs ----
    # Table 3 reports volumes in uL; converted to mL to match the mL/min
    # flows and the ng / (ng/mL) dosing and concentration units.
    V_TOT <- fixed(1880 / 1000); label("Total brain volume (mL)")  # Table 3: 1880 uL
    V_ECF <- fixed(290 / 1000);  label("Brain extracellular fluid volume (mL)")  # Table 3: 290 uL
    V_ICF <- fixed(1440 / 1000); label("Brain intracellular fluid volume (mL)")  # Table 3: 1440 uL
    V_LYS <- fixed(18 / 1000);   label("Total lysosomal volume (mL)")  # Table 3: 18 uL (1:80 of brain ICF)
    V_LV  <- fixed(50 / 1000);   label("Lateral-ventricle CSF volume (mL)")  # Table 3: 50 uL
    V_TFV <- fixed(50 / 1000);   label("Third + fourth ventricle CSF volume (mL)")  # Table 3: 50 uL
    V_CM  <- fixed(17 / 1000);   label("Cisterna-magna CSF volume (mL)")  # Table 3: 17 uL
    V_SAS <- fixed(180 / 1000);  label("Subarachnoid-space CSF volume (mL)")  # Table 3: 180 uL
    V_MV  <- fixed(60 / 1000);   label("Brain microvascular volume (mL)")  # Table 3: 60 uL

    Q_CBF <- fixed(1.2);    label("Cerebral blood flow (mL/min)")  # Table 3
    Q_ECF <- fixed(0.0002); label("Brain extracellular-fluid bulk flow (mL/min)")  # Table 3
    Q_CSF <- fixed(0.0022); label("Cerebrospinal-fluid flow (mL/min)")  # Table 3

    SA_BBB    <- fixed(263);  label("Blood-brain barrier surface area (cm^2)")  # Table 3
    SA_BCSFB1 <- fixed(12.5); label("BCSFB surface area around the lateral ventricle (cm^2)")  # Table 3 footnote d
    SA_BCSFB2 <- fixed(12.5); label("BCSFB surface area around the third + fourth ventricle (cm^2)")  # Table 3 footnote d
    SA_BCM    <- fixed(3000); label("Total brain cell-membrane surface area (cm^2)")  # Table 3
    SA_LYSO   <- fixed(1440); label("Total lysosomal membrane surface area (cm^2)")  # Table 3

    # Effective surface-area fractions. Table 3 footnotes b and c give these as
    # PERCENTAGES: 99.8 percent of each barrier is available for transcellular
    # diffusion, 0.006 percent of the BBB and 0.016 percent of the BCSFB for
    # paracellular diffusion. The control streams (S3-S13) hard-code the same
    # numbers as the bare fractions 0.998, 0.00006 and 0.00016.
    f_trans      <- fixed(0.998);   label("Transcellular effective surface-area fraction (unitless)")  # Table 3 footnotes b, c
    f_para_BBB   <- fixed(0.00006); label("Paracellular effective surface-area fraction, BBB (unitless)")  # Table 3 footnote b
    f_para_BCSFB <- fixed(0.00016); label("Paracellular effective surface-area fraction, BCSFB (unitless)")  # Table 3 footnote c

    # Table 3 gives the BBB width as 0.3-0.5 um and states that 0.5 um was used.
    # The same width is applied at the BCSFB: every control stream uses one
    # THETA for the ratio Daq / width across both barriers.
    w_BBB <- fixed(0.5 * 1e-4); label("Barrier width (cm)")  # Table 3: 0.5 um used in the model
  })

  model({
    # ---- 1. Empirical plasma pharmacokinetics ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    # ---- 2. Passive diffusion clearances (Eqs. 4, 5, 8, 9) ----
    # Paracellular: Qp = (Daq / width) * paracellular surface area.
    # Transcellular: Qt = 1/2 * P0 * transcellular surface area, where the
    # factor 1/2 corrects for passage over two membranes rather than one.
    Qp_BBB    <- (Daq / w_BBB) * SA_BBB * f_para_BBB
    Qt_BBB    <- 0.5 * P0 * SA_BBB * f_trans
    Qp_BCSFB1 <- (Daq / w_BBB) * SA_BCSFB1 * f_para_BCSFB
    Qt_BCSFB1 <- 0.5 * P0 * SA_BCSFB1 * f_trans
    Qp_BCSFB2 <- (Daq / w_BBB) * SA_BCSFB2 * f_para_BCSFB
    Qt_BCSFB2 <- 0.5 * P0 * SA_BCSFB2 * f_trans
    Q_BCM     <- P0 * SA_BCM
    Q_LYSO    <- P0 * SA_LYSO

    # ---- 3. Directional clearances (Eqs. 6, 7 and 18-24; Table 5 footnote) ----
    # Active transport enters as an asymmetry factor on the transcellular
    # term only; the pH-dependent factor scales the whole efflux clearance.
    Q_BBB_in     <- Qp_BBB + Qt_BBB * AFin1
    Q_BBB_out    <- (Qp_BBB + Qt_BBB * AFout1) * PHF1
    Q_BCSFB1_in  <- Qp_BCSFB1 + Qt_BCSFB1 * AFin2
    Q_BCSFB1_out <- (Qp_BCSFB1 + Qt_BCSFB1 * AFout2) * PHF2
    Q_BCSFB2_in  <- Qp_BCSFB2 + Qt_BCSFB2 * AFin3
    Q_BCSFB2_out <- (Qp_BCSFB2 + Qt_BCSFB2 * AFout3) * PHF3

    # ---- 4. Compartment concentrations ----
    C_PL  <- central / vc
    C_MV  <- brain_vascular / V_MV
    C_ECF <- brain_ecf / V_ECF
    C_ICF <- brain_icf / V_ICF
    C_LYS <- brain_lysosome / V_LYS
    C_LV  <- brain_csf_lv / V_LV
    C_TFV <- brain_csf_tfv / V_TFV
    C_CM  <- brain_csf_cm / V_CM
    C_SAS <- brain_csf_sas / V_SAS

    # ---- 5. Plasma disposition (control stream $DES, DADT(1)-DADT(3)) ----
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      (q2 / vc) * central + (q2 / vp2) * peripheral2 -
      Q_CBF * C_PL + Q_CBF * C_MV
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(peripheral2) <- (q2 / vc) * central - (q2 / vp2) * peripheral2

    # ---- 6. CNS PBPK system (control stream $DES, DADT(4)-DADT(12)) ----
    # Brain microvasculature is the donor compartment for BOTH barriers and
    # receives the subarachnoid CSF back-flow, closing the CSF loop.
    d/dt(brain_vascular) <- Q_CBF * C_PL - Q_CBF * C_MV -
      Q_BBB_in * C_MV + Q_BBB_out * C_ECF -
      Q_BCSFB1_in * C_MV + Q_BCSFB1_out * C_LV -
      Q_BCSFB2_in * C_MV + Q_BCSFB2_out * C_TFV +
      Q_CSF * C_SAS
    d/dt(brain_ecf) <- Q_BBB_in * C_MV - Q_BBB_out * C_ECF -
      Q_BCM * PHF4 * C_ECF + Q_BCM * PHF5 * C_ICF -
      Q_ECF * C_ECF
    d/dt(brain_icf) <- Q_BCM * PHF4 * C_ECF - Q_BCM * PHF5 * C_ICF -
      Q_LYSO * PHF6 * C_ICF + Q_LYSO * PHF7 * C_LYS
    d/dt(brain_lysosome) <- Q_LYSO * PHF6 * C_ICF - Q_LYSO * PHF7 * C_LYS
    d/dt(brain_csf_lv) <- Q_BCSFB1_in * C_MV - Q_BCSFB1_out * C_LV +
      Q_ECF * C_ECF - Q_CSF * C_LV
    d/dt(brain_csf_tfv) <- Q_BCSFB2_in * C_MV - Q_BCSFB2_out * C_TFV +
      Q_CSF * C_LV - Q_CSF * C_TFV
    d/dt(brain_csf_cm) <- Q_CSF * C_TFV - Q_CSF * C_CM
    d/dt(brain_csf_sas) <- Q_CSF * C_CM - Q_CSF * C_SAS

    # ---- 7. Observations ----
    # Cc is total plasma concentration, the only endpoint that was fitted.
    # The CNS outputs are unbound concentrations and are predictions.
    Cc <- C_PL
    Cbrain_ecf <- C_ECF
    Cbrain_icf <- C_ICF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_cm <- C_CM
    # Total brain tissue: bound drug in the extracellular space plus the
    # free extracellular, intracellular and lysosomal amounts, over total
    # brain volume (control stream $ERROR, DBR).
    Cbrain_total <- (brain_ecf * BF + brain_ecf + brain_icf + brain_lysosome) / V_TOT
    Cc ~ prop(propSd)
  })
}
