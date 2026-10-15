Hirasawa_2022_inulin_rat_pbpk <- function() {
  description <- paste(
    "PBPK (LeiCNS-PK3.1 CNS physiologically-based model). Preclinical (rat).",
    "Inulin (Unlabelled or [14c]-Labelled) disposition after intra-CSF administration in healthy rats.",
    "LeiCNS-PK3.1 replaces the unidirectional constant CSF flow of",
    "LeiCNS-PK3.0 by bidirectional, site-dependent CSF movement (downward",
    "and upward rates between lateral ventricles, third + fourth ventricles",
    "and cisterna magna, and between cisterna magna and subarachnoid space),",
    "a drug-specific transfer clearance from the subarachnoid space to",
    "plasma, a choroid-plexus microvessel compartment that is the donor for",
    "both blood-CSF barriers, and a correction factor on paracellular",
    "permeability. The brain parenchyma (microvessels, ECF, cell membrane,",
    "ICF, lysosomes) is LeiCNS-PK3.0 with rat physiology. Asymmetry factors",
    "for active transport are back-solved from the reported Kp,uu values.",
    "The three drug-independent CSF movement parameters are fixed to the",
    "sucrose estimates; sasQCSF and CFPPA are the drug-specific estimates.",
    "Plasma side: the empirical 1-compartment unbound plasma PK model of",
    "Supplementary Table S1 (from ref [19], Reed and Woodbury 1963), mass-balanced with the CNS.",
    "No active transport: all asymmetry factors are 1."
  )
  reference <- paste(
    "Hirasawa M, de Lange ECM. Revisiting Cerebrospinal Fluid Flow",
    "Direction and Rate in Physiologically Based Pharmacokinetic Model.",
    "Pharmaceutics. 2022;14(9):1764. doi:10.3390/pharmaceutics14091764.",
    "The rat LeiCNS-PK3.0 physiology (volumes, flows, surface areas, widths,",
    "pH) and the Daq / P0 / Qp / Qt equations are from Saleh MAA, Bloemberg",
    "JM, Elassaiss-Schaap J, de Lange ECM. Drug Distribution in Brain and",
    "Cerebrospinal Fluids in Relation to IC50 Values in Aging and",
    "Alzheimer's Disease, Using the Physiologically Based LeiCNS-PK3.0",
    "Model. J Pharmacokinet Pharmacodyn. 2021;48(5):725-741.",
    "doi:10.1007/s10928-021-09768-7 (Supplementary table 2 and",
    "Supplementary equations; reference [6] of Hirasawa 2022)."
  )
  vignette <- "Hirasawa_2022_leicns_pk31_csf_flow"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  # Choroid-plexus microvessels, separated from the whole-brain microvessels
  # in LeiCNS-PK3.1 as the donor for both blood-CSF barriers.
  paper_specific_compartments <- c("brain_choroid_plexus_vascular")

  compartmentData <- list(
    central = list(analyte = "inulin", units = "ng", specimen = "plasma", verified = TRUE),
    brain_vascular = list(analyte = "inulin", units = "ng", specimen = "plasma", verified = TRUE),
    brain_choroid_plexus_vascular = list(analyte = "inulin", units = "ng", specimen = "plasma", verified = TRUE),
    brain_ecf = list(analyte = "inulin", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_cell_membrane = list(analyte = "inulin", units = "ng", specimen = "tissue", verified = TRUE),
    brain_icf = list(analyte = "inulin", units = "ng", specimen = "tissue", verified = TRUE),
    brain_lysosome = list(analyte = "inulin", units = "ng", specimen = "tissue", verified = TRUE),
    brain_csf_lv = list(analyte = "inulin", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "inulin", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "inulin", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "inulin", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat",
    disease_state = "healthy",
    weight = "250 g assumed (Section 2.5)",
    dose_range = "2.5 mg/kg inulin ICV (12.5 uL, right lateral ventricle) and 2 uCi [14C]inulin IC (100 uL)",
    notes = paste(
      "Data: inulin ICV (Noguchi 2017, ref [20]; 2.5 mg/kg in 12.5 uL) and [14C]inulin IC (Reed and Woodbury 1963, ref [19]; 2 uCi in 100 uL), fitted simultaneously.",
      "Published group-mean profiles digitised by the authors (WebPlotDigitizer);",
      "subject counts are those of the cited source studies and are not",
      "restated. Rats are assumed 250 g (Section 2.5) when converting per-kg",
      "doses. Intra-CSF doses are split across the CSF compartments by the",
      "dosing-volume assumptions of Supplementary Table S3."
    )
  )

  ini({
    # ---- Empirical unbound plasma PK model: Supplementary Table S1, ref [19], Reed and Woodbury 1963 ----
    lcl <- log(0.00463); label("Central clearance, CLcen (mL/min)")  # Table S1 CLcen
    lvc <- log(9.6); label("Central volume, Vcen (mL)")  # Table S1 Vcen

    # ---- Drug physicochemical properties: Supplementary Table S2 ----
    mwt <- fixed(6179.4); label("Molecular weight (g/mol)")  # Table S2
    logp <- fixed(-62); label("Octanol/water partition coefficient, logP (unitless)")  # Table S2
    pka <- fixed(11.27); label("Acidic ionization constant (unitless)")  # Table S2
    pkb <- fixed(-4); label("Basic ionization constant (unitless)")  # Table S2

    # ---- Drug-independent CSF movement: Table 2, [3H]/[14C]sucrose fit ----
    # Reported in uL/min; divided by 1000 to give mL/min.
    # Table 2 footnote a: fixed to the values estimated from the sucrose data.
    Q_CSF_ven_D <- fixed(0.251 / 1000); label("Downward ventricular CSF movement rate, venQCSF,D (mL/min)")  # Table 2: 0.251 uL/min
    Q_CSF_cis_D <- fixed(1.07 / 1000); label("Downward cisternal CSF movement rate CM to SAS, cisQCSF,D (mL/min)")  # Table 2: 1.07 uL/min
    ud_ratio <- fixed(1.59); label("Ratio of upward to downward CSF movement rates (unitless)")  # Table 2: 1.59

    # ---- Drug-specific CNS parameters ----
    Q_SAS <- 1.22 / 1000; label("Transfer clearance from SAS to plasma, sasQCSF (mL/min)")  # Table 2: 1.22 uL/min
    cf_ppa <- 0.262; label("Correction factor on paracellular permeability, CFPPA (unitless)")  # Table 2: 0.262

    # ---- Choroid-plexus microvessels: Section 2.7 ----
    # CP weight 4 mg per 250 g rat; CP blood volume 196 uL/g; CP blood flow
    # 4 mL/min/g tissue.
    V_CPMV <- fixed(0.004 * 0.196); label("Choroid-plexus microvessel volume (mL)")  # Section 2.7: 4 mg x 196 uL/g
    Q_CPBF <- fixed(0.004 * 4); label("Choroid-plexus blood flow (mL/min)")  # Section 2.7: 4 mg x 4 mL/min/g

    # ---- Rat CNS physiology: Saleh 2021 Supplementary table 2 ----
    V_MV <- fixed(0.054); label("Brain microvessel volume (mL)")  # Saleh 2021 Suppl. table 2
    V_ECF <- fixed(0.36); label("Brain ECF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_ICF <- fixed(1.44); label("Brain ICF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_LYS <- fixed(0.018); label("Brain lysosome volume (mL)")  # Saleh 2021 Suppl. table 2
    # Brain phospholipid volume = 5% of total brain volume (1.8 mL).
    V_BCM <- fixed(0.05 * 1.8); label("Brain cell-membrane (phospholipid) volume (mL)")  # Saleh 2021 main text and Suppl. table 2
    # Both lateral ventricles (3.75 uL each, Hirasawa 2022 Section 2.5).
    V_LV <- fixed(0.0075); label("Lateral-ventricle CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_TFV <- fixed(0.0075); label("Third + fourth ventricle CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_CM <- fixed(0.017); label("Cisterna-magna CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_SAS <- fixed(0.135); label("Subarachnoid-space CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    Q_CBF <- fixed(2.87); label("Cerebral blood flow (mL/min)")  # Saleh 2021 Suppl. table 2
    Q_ECF <- fixed(0.0002); label("Brain ECF bulk flow (mL/min)")  # Saleh 2021 Suppl. table 2
    SA_BBB <- fixed(155); label("BBB surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    # Total BCSFB area; half each at the lateral and third + fourth ventricles.
    SA_BCSFB <- fixed(25); label("BCSFB total surface area (cm^2)")  # Saleh 2021 Suppl. table 2 (footnote 13)
    SA_BCM <- fixed(4250); label("Brain cell-membrane surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    SA_LYSO <- fixed(2700); label("Lysosomal membrane surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    # Effective surface areas are percentages in the source table.
    f_trans <- fixed(0.998); label("Transcellular effective surface-area fraction (unitless)")  # Saleh 2021 Suppl. table 2: 99.8%
    f_para_BBB <- fixed(0.006 / 100); label("Paracellular effective surface-area fraction, BBB (unitless)")  # Saleh 2021 Suppl. table 2: 0.006%
    f_para_BCSFB <- fixed(0.05 / 100); label("Paracellular effective surface-area fraction, BCSFB (unitless)")  # Saleh 2021 Suppl. table 2: 0.05%
    # Widths in um, converted to cm; the BCSFB row shares the BBB value
    # (merged table cell).
    w_BBB <- fixed(0.5 * 1e-4); label("BBB width (cm)")  # Saleh 2021 Suppl. table 2: 0.5 um
    w_BCSFB <- fixed(0.5 * 1e-4); label("BCSFB width (cm)")  # Saleh 2021 Suppl. table 2: 0.5 um
    pH_PL <- fixed(7.4); label("Plasma and microvessel pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_ECF <- fixed(7.3); label("Brain ECF pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_CSF <- fixed(7.3); label("CSF pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_ICF <- fixed(7.0); label("Brain ICF pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_LYS <- fixed(5.0); label("Lysosomal pH (unitless)")  # Saleh 2021 Suppl. table 2

    # ---- Residual unexplained variability on plasma: Table S1 ----
    propSd <- 0.0579; label("Proportional residual error (fraction)")  # Table S1: 5.79%
  })

  model({
    # ---- 1. Empirical unbound plasma PK (Table S1) ----
    cl <- exp(lcl)
    vc <- exp(lvc)

    # ---- 2. pH-dependent neutral fractions (Saleh 2021 equations) ----
    PHF_MV <- (1 / (1 + 10^(pH_PL - pka))) * (1 / (1 + 10^(pkb - pH_PL)))
    PHF_ECF <- (1 / (1 + 10^(pH_ECF - pka))) * (1 / (1 + 10^(pkb - pH_ECF)))
    PHF_CSF <- (1 / (1 + 10^(pH_CSF - pka))) * (1 / (1 + 10^(pkb - pH_CSF)))
    PHF_ICF <- (1 / (1 + 10^(pH_ICF - pka))) * (1 / (1 + 10^(pkb - pH_ICF)))
    PHF_LYS <- (1 / (1 + 10^(pH_LYS - pka))) * (1 / (1 + 10^(pkb - pH_LYS)))

    # ---- 3. Passive diffusion (Saleh 2021 Supplementary equations) ----
    # Both regressions are per-second; x 60 gives per-minute.
    Daq <- 10^(-4.113 - 0.4609 * log10(mwt)) * 60
    P0 <- 10^(0.939 * logp - 6.21) * 60
    # CFPPA scales the paracellular clearance at the BBB and both BCSFBs.
    Qp_BBB <- cf_ppa * (Daq / w_BBB) * SA_BBB * f_para_BBB
    Qt_BBB <- 0.5 * P0 * SA_BBB * f_trans
    Qp_BCSFB <- cf_ppa * (Daq / w_BCSFB) * (SA_BCSFB / 2) * f_para_BCSFB
    Qt_BCSFB <- 0.5 * P0 * (SA_BCSFB / 2) * f_trans
    CL_wo <- P0 * SA_BCM
    CL_ow <- CL_wo / 10^logp
    Q_LYSO <- P0 * SA_LYSO

    # ---- 4. Bidirectional CSF movement (Figure 1B) ----
    Q_CSF_ven_U <- ud_ratio * Q_CSF_ven_D
    Q_CSF_cis_U <- ud_ratio * Q_CSF_cis_D

    # ---- 5. Asymmetry factors ----
    # No active transport or CNS metabolism (Section 2.6): all factors are 1.
    AF_BBB_in <- 1
    AF_BBB_ef <- 1
    AF_LV_in <- 1
    AF_LV_ef <- 1
    AF_TFV_in <- 1
    AF_TFV_ef <- 1

    # ---- 6. Concentrations (all unbound) ----
    C_PL <- central / vc
    C_MV <- brain_vascular / V_MV
    C_CP <- brain_choroid_plexus_vascular / V_CPMV
    C_ECF <- brain_ecf / V_ECF
    C_BCM <- brain_cell_membrane / V_BCM
    C_ICF <- brain_icf / V_ICF
    C_LYS <- brain_lysosome / V_LYS
    C_LV <- brain_csf_lv / V_LV
    C_TFV <- brain_csf_tfv / V_TFV
    C_CM <- brain_csf_cm / V_CM
    C_SAS <- brain_csf_sas / V_SAS

    # ---- 7. Barrier fluxes ----
    # Paracellular transport carries all unbound drug; transcellular only the
    # neutral fraction, scaled by the asymmetry factor.
    inBBB <- (Qp_BBB + AF_BBB_in * Qt_BBB * PHF_MV) * C_MV
    efBBB <- (Qp_BBB + AF_BBB_ef * Qt_BBB * PHF_ECF) * C_ECF
    inLV <- (Qp_BCSFB + AF_LV_in * Qt_BCSFB * PHF_MV) * C_CP
    efLV <- (Qp_BCSFB + AF_LV_ef * Qt_BCSFB * PHF_CSF) * C_LV
    inTFV <- (Qp_BCSFB + AF_TFV_in * Qt_BCSFB * PHF_MV) * C_CP
    efTFV <- (Qp_BCSFB + AF_TFV_ef * Qt_BCSFB * PHF_CSF) * C_TFV

    # ---- 8. Plasma, mass-balanced with the CNS (Figure 1B) ----
    d/dt(central) <- -cl * C_PL -
      Q_CBF * (C_PL - C_MV) - Q_CPBF * (C_PL - C_CP) +
      Q_SAS * C_SAS

    # ---- 9. CNS (Figure 1B) ----
    d/dt(brain_vascular) <- Q_CBF * (C_PL - C_MV) - inBBB + efBBB
    d/dt(brain_choroid_plexus_vascular) <- Q_CPBF * (C_PL - C_CP) - inLV + efLV - inTFV + efTFV
    d/dt(brain_ecf) <- inBBB - efBBB - Q_ECF * C_ECF -
      CL_wo * PHF_ECF * C_ECF + CL_ow * C_BCM
    d/dt(brain_cell_membrane) <- CL_wo * PHF_ECF * C_ECF + CL_wo * PHF_ICF * C_ICF -
      2 * CL_ow * C_BCM
    d/dt(brain_icf) <- CL_ow * C_BCM - CL_wo * PHF_ICF * C_ICF -
      Q_LYSO * PHF_ICF * C_ICF + Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_lysosome) <- Q_LYSO * PHF_ICF * C_ICF - Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_csf_lv) <- inLV - efLV + Q_ECF * C_ECF -
      Q_CSF_ven_D * C_LV + Q_CSF_ven_U * C_TFV
    d/dt(brain_csf_tfv) <- inTFV - efTFV + Q_CSF_ven_D * C_LV -
      (Q_CSF_ven_U + Q_CSF_ven_D) * C_TFV + Q_CSF_ven_U * C_CM
    d/dt(brain_csf_cm) <- Q_CSF_ven_D * C_TFV - (Q_CSF_ven_U + Q_CSF_cis_D) * C_CM +
      Q_CSF_cis_U * C_SAS
    d/dt(brain_csf_sas) <- Q_CSF_cis_D * C_CM - (Q_CSF_cis_U + Q_SAS) * C_SAS

    # ---- 10. Observations ----
    # Cc is the unbound plasma concentration (Table S1 models are fitted to
    # unbound plasma data); the CSF outputs are unbound by assumption
    # (Section 2.1: unbound fraction in CSF = 1).
    Cc <- C_PL
    Cbrain_ecf <- C_ECF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_tfv <- C_TFV
    Cbrain_csf_cm <- C_CM
    Cbrain_csf_sas <- C_SAS
    Cc ~ prop(propSd)
  })
}
