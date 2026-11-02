Gulave_2025_morphine_pbpk <- function() {
  description <- paste(
    "PBPK (LeiCNS-PK3.0 CNS physiologically-based model, human version).",
    "Unbound morphine distribution into human brain extracellular fluid",
    "(brain ECF) after intravenous dosing, with nonlinear blood-brain",
    "barrier (BBB) transport and an optional P-glycoprotein (P-gp)",
    "drug-drug interaction at the BBB. The plasma input is the",
    "one-compartment population PK model of Oosten 2017 (70 kg adult)",
    "acting as a forcing function. The CNS structure is LeiCNS-PK3.0:",
    "brain microvasculature, brain ECF, a phospholipid brain-cell-membrane",
    "compartment, brain intracellular fluid, lysosomes and four CSF",
    "compartments. The extent of BBB transport Kp,uu,BBB is a decreasing",
    "function of the instantaneous unbound plasma concentration (above 1",
    "below about 109 ng/mL, below 1 above it); the BBB asymmetry factor",
    "back-calculated from it is scaled by the human-to-rat P-gp expression",
    "ratio fAF_BBB = 0.218. A P-gp inhibitor enters through the",
    "Rowland-Matin term 1 + fu,p*[I]/Ki on the transcellular efflux",
    "clearance, supplied as the covariate INH_PGP_RATIO (0 = no inhibitor)."
  )
  reference <- paste(
    "Gulave B, Lesmana A, de Lange ECM, van Hasselt JGC. Do",
    "P-glycoprotein-mediated drug-drug interactions at the blood-brain",
    "barrier impact morphine brain distribution? J Pharmacokinet",
    "Pharmacodyn. 2025;52(1):10. doi:10.1007/s10928-024-09957-0.",
    "Human LeiCNS-PK3.0 physiology, the Daq / P0 / Qp / Qt regressions",
    "and the asymmetry-factor closed forms: Saleh MAA, Loo CF,",
    "Elassaiss-Schaap J, de Lange ECM. J Pharmacokinet Pharmacodyn.",
    "2021;48(5):725-741. doi:10.1007/s10928-021-09768-7 (reference [19]).",
    "Morphine LeiCNS-PK3.0 extension and the fAF_BBB translation:",
    "Gulave B, Budda D, Saleh MAA, van Hasselt JGC, de Lange ECM. Eur J",
    "Pharm Sci. 2023;187:106482. doi:10.1016/j.ejps.2023.106482",
    "(reference [18]). Plasma PK: Oosten AW et al. Clin Pharmacokinet.",
    "2017;56(7):733-746. doi:10.1007/s40262-016-0471-7 (reference [20])."
  )
  vignette <- "Gulave_2025_morphine_pbpk"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  covariateData <- list(
    INH_PGP_RATIO = list(
      description = paste(
        "Unbound P-gp inhibitor concentration divided by its inhibitory",
        "constant, fu,p * [I] / Ki (or IC50). Divides the BBB transcellular",
        "efflux clearance by (1 + INH_PGP_RATIO) (Rowland-Matin static",
        "model, Gulave 2025 Eqs. 2-4). The inhibitor is assumed at steady",
        "state, so the value is constant within a subject."
      ),
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = "0 (no inhibitor; the uninhibited model)",
      notes = paste(
        "Supplementary Table 1 tabulates the term per inhibitor in its 'DDI'",
        "column. Those printed values are 1000-fold larger than fu,p * Css / Ki",
        "computed in consistent units (Css in ng/mL converted to uM with the",
        "inhibitor molecular weight, Ki in uM); the paper's Figures 2-3 are",
        "reproduced only with the printed values, so the published",
        "simulations used them. See the vignette Assumptions and deviations."
      ),
      source_name = "fu,p,inhibitor * C_inhibitor / Ki (Eq. 2); 'DDI' column of Supplementary Table 1"
    )
  )

  compartmentData <- list(
    central = list(analyte = "morphine", units = "ng", specimen = "plasma", verified = TRUE),
    brain_vascular = list(analyte = "morphine", units = "ng", specimen = "plasma", verified = TRUE),
    brain_ecf = list(analyte = "morphine", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_cell_membrane = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_icf = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_lysosome = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_csf_lv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2000,
    n_studies = 0,
    weight_range = "70 kg (single typical adult body weight; no covariates)",
    disease_state = paste(
      "Virtual adult population for a simulation study; the plasma PK model",
      "was developed in adult cancer patients (Oosten 2017)"
    ),
    dose_range = paste(
      "Morphine 10 mg IV every 4 h for 7 days (clinical regimen); 50 and",
      "100 mg every 4 h in the sensitivity analyses"
    ),
    co_medication = "34 P-gp inhibitors at their clinical maintenance doses (Supplementary Table 1)",
    notes = paste(
      "No clinical data were fitted. 2000 virtual subjects with interindividual",
      "variability in plasma PK only (Methods 'Simulation scenarios'). All CNS",
      "parameters are fixed to human LeiCNS-PK3.0 physiology (Saleh 2021",
      "Supplementary Table 2) or derived from morphine physicochemical",
      "properties (Table 1)."
    )
  )

  ini({
    # ---- Plasma PK: Oosten 2017 one-compartment model (Table 1) ----
    lcl <- log(91.9 * 1000 / 60); label("Morphine clearance (mL/min)")  # Table 1: 91.9 L/hr (= 1531.67 mL/min, Gulave 2023 Table 1)
    lvc <- log(278 * 1000); label("Morphine central volume (mL)")  # Table 1: 278 L

    # ---- Morphine physicochemical properties: Table 1 (DrugBank) ----
    mwt <- fixed(285.34); label("Molecular weight (g/mol)")  # Table 1
    logp <- fixed(0.99); label("Lipophilicity, logP (unitless)")  # Table 1
    pka <- fixed(10.26); label("Acidic ionization constant (unitless)")  # Table 1
    pkb <- fixed(9.12); label("Basic ionization constant (unitless)")  # Table 1
    fu_plasma <- fixed(0.65); label("Unbound fraction in plasma (unitless)")  # Table 1

    # ---- Nonlinear Kp,uu,BBB: Table 1 ----
    # Kp,uu,BBB = (14.543 * Cu^-0.895) / (5.902^-0.895 + Cu^-0.895)
    #           = kpuu_max / (1 + (Cu / c50_kpuu)^hill_kpuu), Cu in ng/mL.
    kpuu_max <- fixed(14.543); label("Kp,uu,BBB at vanishing unbound plasma concentration (unitless)")  # Table 1 Kp,uu,BBB function
    c50_kpuu <- fixed(5.902); label("Unbound plasma concentration halving Kp,uu,BBB (ng/mL)")  # Table 1 Kp,uu,BBB function
    hill_kpuu <- fixed(0.895); label("Hill exponent of the Kp,uu,BBB function (unitless)")  # Table 1 Kp,uu,BBB function

    # ---- Rat-to-human translation of the BBB asymmetry factor: Table 1 ----
    faf_bbb <- fixed(0.21815); label("Human-to-rat BBB P-gp expression ratio scaling the asymmetry factor (unitless)")  # Table 1 fAFBBB (4.21 / 19.28 fmol/ug)

    # ---- Human CNS physiology: Saleh 2021 Supplementary Table 2 ----
    V_MV <- fixed(45); label("Brain microvasculature volume (mL)")  # Saleh 2021 Suppl. table 2
    V_ECF <- fixed(253); label("Brain ECF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_ICF <- fixed(1000); label("Brain ICF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_LYS <- fixed(12.5); label("Total lysosome volume (mL)")  # Saleh 2021 Suppl. table 2
    V_LV <- fixed(20); label("Lateral ventricle CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_TFV <- fixed(3); label("Third + fourth ventricle CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_CM <- fixed(1); label("Cisterna magna CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    V_SAS <- fixed(116); label("Subarachnoid space CSF volume (mL)")  # Saleh 2021 Suppl. table 2
    # Brain cell membrane (phospholipid) volume is 5% of total brain volume
    # (Saleh 2021 main text); total brain volume 1250 mL.
    V_BCM <- fixed(0.05 * 1250); label("Brain phospholipid volume (mL)")  # Saleh 2021: 5% of V_tot = 1250 mL

    Q_CBF <- fixed(689); label("Cerebral blood flow (mL/min)")  # Saleh 2021 Suppl. table 2
    Q_ECF <- fixed(0.2); label("Brain ECF bulk flow (mL/min)")  # Saleh 2021 Suppl. table 2
    Q_CSF <- fixed(0.42); label("CSF flow (mL/min)")  # Saleh 2021 Suppl. table 2

    SA_BBB <- fixed(150000); label("BBB surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    SA_BCSFB <- fixed(15000); label("BCSFB total surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    SA_BCM <- fixed(2666520); label("Brain cell membrane surface area (cm^2)")  # Saleh 2021 Suppl. table 2
    SA_LYSO <- fixed(1980260); label("Lysosomal membrane surface area (cm^2)")  # Saleh 2021 Suppl. table 2

    # Effective surface areas are printed as percentages.
    f_trans <- fixed(99.8 / 100); label("Transcellular effective surface-area fraction (unitless)")  # Saleh 2021 Suppl. table 2: 99.8%
    f_para_BBB <- fixed(0.004 / 100); label("Paracellular effective surface-area fraction, BBB (unitless)")  # Saleh 2021 Suppl. table 2: 0.004%
    f_para_BCSFB <- fixed(0.016 / 100); label("Paracellular effective surface-area fraction, BCSFB (unitless)")  # Saleh 2021 Suppl. table 2: 0.016%

    # Barrier widths in um converted to cm. The BCSFB width cell is merged
    # with the BBB row in the source table (0.5 um).
    w_BBB <- fixed(0.5 * 1e-4); label("BBB width (cm)")  # Saleh 2021 Suppl. table 2: 0.5 um
    w_BCSFB <- fixed(0.5 * 1e-4); label("BCSFB width (cm)")  # Saleh 2021 Suppl. table 2: 0.5 um (merged cell)

    pH_PL <- fixed(7.4); label("Plasma and brain microvasculature pH (unitless)")  # Saleh 2021 Suppl. table 2 (merged cell)
    pH_ECF <- fixed(7.3); label("Brain ECF pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_ICF <- fixed(7); label("Brain ICF pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_LYS <- fixed(5); label("Lysosomal pH (unitless)")  # Saleh 2021 Suppl. table 2
    pH_CSF <- fixed(7.3); label("CSF pH (unitless)")  # Saleh 2021 Suppl. table 2

    # ---- Interindividual variability: Table 1 (reported as variances) ----
    etalcl ~ 0.222  # Table 1: IIV central clearance (as variance)
    etalvc ~ 0.747  # Table 1: IIV central compartment volume (as variance)

    # ---- Residual error: Table 1 (reported as a variance) ----
    propSd <- sqrt(0.286); label("Proportional residual error on plasma morphine (fraction)")  # Table 1: 0.286 as variance
  })

  model({
    # ---- 1. Plasma PK (Oosten 2017; 70 kg) ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # ---- 2. pH-dependent neutral (uncharged) fractions ----
    # Acidic term 1/(1 + 10^(pH - pKa)); basic term 1/(1 + 10^(pKb - pH)).
    PHF_MV <- (1 / (1 + 10^(pH_PL - pka))) * (1 / (1 + 10^(pkb - pH_PL)))
    PHF_ECF <- (1 / (1 + 10^(pH_ECF - pka))) * (1 / (1 + 10^(pkb - pH_ECF)))
    PHF_ICF <- (1 / (1 + 10^(pH_ICF - pka))) * (1 / (1 + 10^(pkb - pH_ICF)))
    PHF_LYS <- (1 / (1 + 10^(pH_LYS - pka))) * (1 / (1 + 10^(pkb - pH_LYS)))
    PHF_CSF <- (1 / (1 + 10^(pH_CSF - pka))) * (1 / (1 + 10^(pkb - pH_CSF)))

    # ---- 3. Passive diffusivity and permeability ----
    # Both regressions return per-second units; x 60 converts to per-minute.
    # Saleh 2021 supplementary equations; Gulave 2025 Supplementary Eqs. 1-2.
    Daq <- 10^(-4.113 - 0.4609 * log10(mwt)) * 60
    P0 <- 10^(0.939 * logp - 6.21) * 60

    # ---- 4. Barrier and intracellular clearances ----
    # The total BCSFB surface area is split 50/50 between the lateral-
    # ventricle and third + fourth-ventricle barriers.
    Qp_BBB <- (Daq / w_BBB) * (SA_BBB * f_para_BBB)
    Qt_BBB <- 0.5 * P0 * (SA_BBB * f_trans)
    Qp_BCSFB <- (Daq / w_BCSFB) * ((SA_BCSFB / 2) * f_para_BCSFB)
    Qt_BCSFB <- 0.5 * P0 * ((SA_BCSFB / 2) * f_trans)
    CL_wo <- P0 * SA_BCM
    CL_ow <- CL_wo / (10^logp)
    Q_LYSO <- P0 * SA_LYSO

    # ---- 5. Concentration-dependent Kp,uu (Table 1) ----
    # Evaluated at the instantaneous unbound plasma concentration (ng/mL).
    # Kp,uu at the BCSFB is not reported for morphine and is set equal to
    # Kp,uu,BBB; the CSF chain is downstream of brain ECF.
    Cu_plasma <- fu_plasma * (central / vc)
    kpuu_ecf <- kpuu_max / (1 + (Cu_plasma / c50_kpuu)^hill_kpuu)
    kpuu_lv <- kpuu_ecf
    kpuu_cm <- kpuu_ecf

    # ---- 6. Asymmetry factors (Saleh 2021 closed forms) ----
    # Kp,uu >= 1: efflux AF = 1, influx AF back-solved; Kp,uu < 1: influx
    # AF = 1, efflux AF back-solved (Gulave 2025 Supplementary Eqs. 3-4).
    # Only the back-solved BBB factor is scaled by fAF_BBB (Gulave 2023
    # Methods 2.2: rat-derived Kp,uu translated by the P-gp expression ratio).
    if (kpuu_ecf >= 1) {
      AF_BBB_ef <- 1
      AF_BBB_in <- faf_bbb * (-(kpuu_ecf * Q_CBF * Q_ECF + kpuu_cm * Qp_BBB * Q_CSF +
        Q_CBF * (kpuu_ecf * (PHF_ECF * Qt_BBB + Qp_BBB) - Qp_BBB)) /
        (PHF_MV * Qt_BBB * (kpuu_cm * Q_CSF - Q_CBF)))
      AF_LV_ef <- 1
      AF_TFV_ef <- 1
      AF_LV_in <- (kpuu_ecf * Q_CBF * Q_ECF - Q_CSF * (kpuu_lv * Q_CBF + kpuu_cm * Qp_BCSFB) -
        Q_CBF * (kpuu_lv * (PHF_CSF * Qt_BCSFB + Qp_BCSFB) - Qp_BCSFB)) /
        (PHF_MV * Qt_BCSFB * (kpuu_cm * Q_CSF - Q_CBF))
      AF_TFV_in <- (Q_CSF * (Q_CBF * (kpuu_lv - kpuu_cm) - kpuu_cm * Qp_BCSFB) +
        Q_CBF * (-kpuu_cm * (PHF_CSF * Qt_BCSFB + Qp_BCSFB) + Qp_BCSFB)) /
        (PHF_MV * Qt_BCSFB * (kpuu_cm * Q_CSF - Q_CBF))
    } else {
      AF_BBB_in <- 1
      AF_BBB_ef <- faf_bbb * (-(kpuu_ecf * Q_CBF * Q_ECF +
        kpuu_cm * Q_CSF * (PHF_MV * Qt_BBB + Qp_BBB) +
        Q_CBF * (kpuu_ecf * Qp_BBB - PHF_MV * Qt_BBB - Qp_BBB)) /
        (kpuu_ecf * PHF_ECF * Qt_BBB * Q_CBF))
      AF_LV_in <- 1
      AF_TFV_in <- 1
      AF_LV_ef <- (kpuu_ecf * Q_CBF * Q_ECF -
        Q_CSF * (kpuu_cm * (PHF_MV * Qt_BCSFB + Qp_BCSFB) + kpuu_lv * Q_CBF) +
        Q_CBF * (-kpuu_lv * Qp_BCSFB + PHF_MV * Qt_BCSFB + Qp_BCSFB)) /
        (kpuu_lv * PHF_CSF * Qt_BCSFB * Q_CBF)
      AF_TFV_ef <- (Q_CSF * (Q_CBF * (kpuu_lv - kpuu_cm) -
        kpuu_cm * (PHF_MV * Qt_BCSFB + Qp_BCSFB)) +
        Q_CBF * (-kpuu_cm * Qp_BCSFB + PHF_MV * Qt_BCSFB + Qp_BCSFB)) /
        (kpuu_cm * PHF_CSF * Qt_BCSFB * Q_CBF)
    }

    # ---- 7. P-gp DDI at the BBB (Eqs. 2-4) ----
    # QtBBB_DDI = QtBBB * AFeBBB / (1 + fu,p * [I] / Ki), applied to the
    # transcellular efflux term only.
    ddi_effect <- 1 + INH_PGP_RATIO

    # ---- 8. CNS compartment concentrations ----
    C_MVu <- fu_plasma * (brain_vascular / V_MV)
    C_ECF <- brain_ecf / V_ECF
    C_BCM <- brain_cell_membrane / V_BCM
    C_ICF <- brain_icf / V_ICF
    C_LYS <- brain_lysosome / V_LYS
    C_LV <- brain_csf_lv / V_LV
    C_TFV <- brain_csf_tfv / V_TFV
    C_CM <- brain_csf_cm / V_CM
    C_SAS <- brain_csf_sas / V_SAS

    # ---- 9. Barrier fluxes ----
    # Neutral drug crosses by both routes; charged drug only paracellularly.
    inBBB <- Qp_BBB * C_MVu + Qt_BBB * AF_BBB_in * PHF_MV * C_MVu
    efBBB <- Qp_BBB * C_ECF + Qt_BBB * (AF_BBB_ef / ddi_effect) * PHF_ECF * C_ECF
    inLV <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_LV_in * PHF_MV * C_MVu
    efLV <- Qp_BCSFB * C_LV + Qt_BCSFB * AF_LV_ef * PHF_CSF * C_LV
    inTFV <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_TFV_in * PHF_MV * C_MVu
    efTFV <- Qp_BCSFB * C_TFV + Qt_BCSFB * AF_TFV_ef * PHF_CSF * C_TFV

    # ---- 10. Plasma model (forcing function) ----
    # CNS uptake does not deplete plasma (LeiCNS-PK3.0 plasma input).
    d/dt(central) <- -(cl / vc) * central

    # ---- 11. CNS PBPK system ----
    d/dt(brain_vascular) <- Q_CBF * (central / vc - brain_vascular / V_MV) -
      inBBB + efBBB - inLV + efLV - inTFV + efTFV
    d/dt(brain_ecf) <- inBBB - efBBB - Q_ECF * C_ECF -
      CL_wo * PHF_ECF * C_ECF + CL_ow * C_BCM
    d/dt(brain_cell_membrane) <- CL_wo * PHF_ECF * C_ECF + CL_wo * PHF_ICF * C_ICF -
      2 * CL_ow * C_BCM
    d/dt(brain_icf) <- CL_ow * C_BCM - CL_wo * PHF_ICF * C_ICF -
      Q_LYSO * PHF_ICF * C_ICF + Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_lysosome) <- Q_LYSO * PHF_ICF * C_ICF - Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_csf_lv) <- inLV - efLV + Q_ECF * C_ECF - Q_CSF * C_LV
    d/dt(brain_csf_tfv) <- inTFV - efTFV + Q_CSF * C_LV - Q_CSF * C_TFV
    d/dt(brain_csf_cm) <- Q_CSF * C_TFV - Q_CSF * C_CM
    d/dt(brain_csf_sas) <- Q_CSF * C_CM - Q_CSF * C_SAS

    # ---- 12. Observations ----
    # Cc: total plasma morphine; Ccu: unbound plasma; Cbrain_ecf: unbound
    # brain ECF morphine, the paper's endpoint.
    Cc <- central / vc
    Ccu <- Cu_plasma
    Cbrain_ecf <- C_ECF
    Cbrain_icf <- C_ICF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_cm <- C_CM
    Cc ~ prop(propSd)
  })
}
