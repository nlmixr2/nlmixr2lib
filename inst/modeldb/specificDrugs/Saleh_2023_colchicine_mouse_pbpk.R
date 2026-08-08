Saleh_2023_colchicine_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (LeiCNS-PK3.0 CNS physiologically-based model, mouse version).",
    "Preclinical (mouse, NMRI).",
    "Nine-compartment CNS PBPK model predicting unbound colchicine",
    "concentrations in brain extracellular fluid (brain ECF) after a",
    "single 1.5 mg/kg intravenous dose. The CNS structure is the mouse",
    "re-parameterisation of LeiCNS-PK3.0: brain microvasculature",
    "(brain_vascular), brain ECF, a phospholipid brain-cell-membrane",
    "binding compartment (brain_cell_membrane), brain intracellular",
    "fluid (brain_icf), lysosomes (brain_lysosome), and the four CSF",
    "compartments (lateral ventricles, third + fourth ventricles,",
    "cisterna magna, subarachnoid space) draining in series back to",
    "plasma. Transport across the BBB and BCSFB is the sum of a",
    "paracellular clearance (Qp, charged + neutral drug) and a",
    "transcellular clearance (Qt, neutral drug only, scaled by",
    "pH-dependent neutral fractions PHF and by asymmetry factors AF",
    "that encode active transport). Every CNS parameter is fixed to",
    "mouse physiology (Table III) or derived from the drug's",
    "physicochemical properties (Table I); none was fitted here. The",
    "BBB asymmetry factors are back-calculated from Kp,uu,BBB =",
    "0.14 (taken from ref [20]; Table V). The plasma PK model is the",
    "empirical two-compartment model of Table IV (in-house NONMEM fit) and acts",
    "purely as a forcing function: brain uptake does not deplete",
    "plasma, exactly as published."
  )
  reference <- paste(
    "Saleh MAA, Gulave B, Campagne O, Stewart CF, Elassaiss-Schaap J,",
    "de Lange ECM. Using the LeiCNS-PK3.0 Physiologically-Based",
    "Pharmacokinetic Model to Predict Brain Extracellular Fluid",
    "Pharmacokinetics in Mice. Pharm Res. 2023;40(11):2555-2566.",
    "doi:10.1007/s11095-023-03554-5.",
    "CNS model structure, the supplementary equations for Daq / P0 /",
    "Qp / Qt and the asymmetry-factor closed forms are defined in",
    "Saleh MAA et al. J Pharmacokinet Pharmacodyn. 2021;48(5):725-741.",
    "doi:10.1007/s10928-021-09768-7 (reference [13] of the 2023 paper)."
  )
  vignette <- "Saleh_2023_leicns_pk30_mouse_brain"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  compartmentData <- list(
    central              = list(analyte = "colchicine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "colchicine", units = "ng", specimen = "plasma", verified = TRUE),
    brain_vascular       = list(analyte = "colchicine", units = "ng", specimen = "plasma", verified = TRUE),
    brain_ecf            = list(analyte = "colchicine", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_cell_membrane  = list(analyte = "colchicine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_icf            = list(analyte = "colchicine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_lysosome       = list(analyte = "colchicine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_csf_lv         = list(analyte = "colchicine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv        = list(analyte = "colchicine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm         = list(analyte = "colchicine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas        = list(analyte = "colchicine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species      = "mouse (NMRI)",
    disease_state = "healthy / tumour-bearing laboratory mice; CNS physiology parameterised for the healthy mouse",
    dose_range   = "1.5 mg/kg IV (single dose; Table IV)",
    notes        = paste(
      "Unbound plasma and microdialysis brain-ECF concentration-time",
      "data for colchicine came from ref [20]",
      "(Table II). Subject counts per drug are not reported. Doses in",
      "Table IV are mg/kg while all volumes in Table III / IV are",
      "absolute (mL); converting a mg/kg dose to the ng amount this",
      "model expects therefore needs a body weight, which the paper",
      "never states. The validation vignette uses 0.025 kg (standard",
      "adult laboratory mouse)."
    )
  )

  ini({
    # ---- Empirical plasma PK model (forcing function): Table IV ----
    lcl <- log(0.487); label("Central clearance (mL/min)")  # Table IV
    lvc <- log(70.1); label("Central compartment volume (mL)")  # Table IV
    lq  <- log(2.11); label("Intercompartmental clearance (mL/min)")  # Table IV
    lvp <- log(705); label("Peripheral compartment volume (mL)")  # Table IV

    # ---- Drug physicochemical properties: Table I ----
    mwt  <- fixed(399); label("Molecular weight (g/mol)")  # Table I
    logp <- fixed(1.59); label("Lipophilicity, logP (unitless)")  # Table I (ALOGPS)
    pka  <- fixed(15.06); label("Acidic ionization constant (unitless)")  # Table I (Chemaxon)
    pkb  <- fixed(-1.2); label("Basic ionization constant (unitless)")  # Table I (Chemaxon)

    # ---- Drug biological properties: Tables II and V ----
    fu_plasma <- fixed(0.61); label("Unbound fraction in plasma (unitless)")  # Table II
    kpuu_ecf  <- fixed(0.14); label("Kp,uu,BBB: unbound brain-ECF-to-plasma ratio (unitless)")  # Table V (ref [20])
    # Kp,uu at the two BCSFB barriers is not reported for mouse. Set equal to
    # Kp,uu,BBB following the authors' own LeiCNS-Source drug_pk.csv default.
    # The CSF chain is strictly downstream of brain ECF, so these do not
    # affect the brain-ECF prediction (see vignette Assumptions).
    kpuu_lv <- fixed(0.14); label("Kp,uu at the BCSFB, lateral ventricles (unitless)")  # assumed = Kp,uu,BBB
    kpuu_cm <- fixed(0.14); label("Kp,uu at the BCSFB, cisterna magna (unitless)")  # assumed = Kp,uu,BBB
    # All 10 drugs have Kp,uu,BBB < 1, so the influx asymmetry factor is 1 and
    # the efflux asymmetry factors are back-solved in model() (Saleh 2021).
    AF_in <- fixed(1); label("Influx asymmetry factor (unitless)")  # Saleh 2021 AF selection rule

    # ---- Mouse CNS physiology: Table III ----
    # Volumes are reported in uL; converted to mL for consistency with the
    # mL/min flows and the ng / (ng/mL) dosing and concentration units.
    V_MV  <- fixed(5 / 1000);      label("Brain microvasculature volume (mL)")  # Table III: 5 uL
    V_ECF <- fixed(67 / 1000);     label("Brain ECF volume (mL)")  # Table III: 67 uL
    V_ICF <- fixed(288 / 1000);    label("Brain ICF volume (mL)")  # Table III: 288 uL (80% of total brain)
    V_LYS <- fixed(3.6 / 1000);    label("Total lysosome volume (mL)")  # Table III: 3.6 uL (1.25% of brain ICF)
    V_LV  <- fixed(1.0275 / 1000); label("Lateral ventricle CSF volume (mL)")  # Table III: 1.0275 uL
    V_TFV <- fixed(2.5 / 1000);    label("Third + fourth ventricle CSF volume (mL)")  # Table III: 2.5 uL
    V_CM  <- fixed(2.13 / 1000);   label("Cisterna magna CSF volume (mL)")  # Table III: 2.13 uL
    V_SAS <- fixed(16.88 / 1000);  label("Subarachnoid space CSF volume (mL)")  # Table III: 16.88 uL
    # Brain phospholipid (cell membrane) volume = 0.05 x total brain volume
    V_BCM <- fixed(0.05 * 360 / 1000); label("Brain phospholipid volume (mL)")  # Table III: 0.05 volume fraction x 360 uL total brain

    Q_CBF <- fixed(0.46134);   label("Cerebral blood flow (mL/min)")  # Table III
    Q_ECF <- fixed(0.0003744); label("Brain ECF bulk flow (mL/min)")  # Table III
    Q_CSF <- fixed(0.000343);  label("CSF flow (mL/min)")  # Table III

    SA_BBB   <- fixed(19.76);  label("BBB surface area (cm^2)")  # Table III
    SA_BCSFB <- fixed(9.88);   label("BCSFB total surface area (cm^2)")  # Table III (50% of BBB)
    SA_BCM   <- fixed(1006.5); label("Brain cell membrane surface area (cm^2)")  # Table III
    SA_LYSO  <- fixed(540);    label("Lysosomal surface area (cm^2)")  # Table III

    # Table III lists the effective surface areas under a '(unitless)' header,
    # but they are PERCENTAGES: Yamamoto 2017 (LeiCNS-PK1.0, reference [43] of
    # this table) states verbatim that '99.8% of total SA_BBB ... is used for
    # transcellular diffusion, whereas 0.006% of total SA_BBB ... [is] used for
    # paracellular diffusion'. Reading 0.006 as a fraction would inflate
    # paracellular clearance 100-fold.
    f_trans       <- fixed(0.998 / 100); label("Effective transcellular surface-area fraction (unitless)")  # Table III: 0.998%
    f_para_BBB    <- fixed(0.006 / 100); label("Effective paracellular surface-area fraction, BBB (unitless)")  # Table III: 0.006%
    f_para_BCSFB  <- fixed(0.05 / 100);  label("Effective paracellular surface-area fraction, BCSFB (unitless)")  # Table III: 0.05%

    # Barrier widths reported in um; converted to cm.
    w_BBB   <- fixed(0.7 * 1e-4); label("BBB width (cm)")  # Table III: 0.7 um
    w_BCSFB <- fixed(1.7 * 1e-4); label("BCSFB width (cm)")  # Table III: 1.7 um

    pH_PL  <- fixed(7.4); label("Plasma pH (unitless)")  # Table III
    pH_ECF <- fixed(7.4); label("Brain ECF pH (unitless)")  # Table III
    pH_ICF <- fixed(7.2); label("Brain ICF pH (unitless)")  # Table III
    pH_LYS <- fixed(5.5); label("Lysosomal pH (unitless)")  # Table III
    pH_CSF <- fixed(7.2); label("CSF pH (unitless)")  # Table III

    # ---- Between-subject variability: Table IV (reported as VARIANCES) ----
    # Table IV reports all interindividual variances as 0 for colchicine;
    # the plasma model is therefore deterministic (no etas), as published.

    # ---- Residual unexplained variability: Table IV (reported as VARIANCES) ----
    propSd <- sqrt(0.0224); label("Proportional residual error (fraction)")  # Table IV: 0.0224 as variance
  })

  model({
    # ---- 1. Empirical plasma PK parameters (Table IV) ----
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ---- 2. pH-dependent neutral (uncharged) fractions ----
    # PHF = fraction of unbound drug that is uncharged at the local pH.
    # Acidic term 1/(1 + 10^(pH - pKa)); basic term 1/(1 + 10^(pKb - pH)).
    # Saleh 2021 supplementary equations.
    PHF_MV  <- (1 / (1 + 10^(pH_PL - pka))) * (1 / (1 + 10^(pkb - pH_PL)))
    PHF_ECF <- (1 / (1 + 10^(pH_ECF - pka))) * (1 / (1 + 10^(pkb - pH_ECF)))
    PHF_ICF <- (1 / (1 + 10^(pH_ICF - pka))) * (1 / (1 + 10^(pkb - pH_ICF)))
    PHF_LYS <- (1 / (1 + 10^(pH_LYS - pka))) * (1 / (1 + 10^(pkb - pH_LYS)))
    PHF_CSF <- (1 / (1 + 10^(pH_CSF - pka))) * (1 / (1 + 10^(pkb - pH_CSF)))

    # ---- 3. Passive diffusivity and permeability ----
    # Both regressions return per-second units; x 60 converts to per-minute,
    # pinned numerically against the published rat values in Yamamoto 2017
    # Tables 4-5. Saleh 2021 supplementary equations.
    Daq <- 10^(-4.113 - 0.4609 * log10(mwt)) * 60
    P0  <- 10^(0.939 * logp - 6.21) * 60

    # ---- 4. Barrier and intracellular clearances ----
    # The total BCSFB surface area is split 50/50 between the lateral-
    # ventricle and third+fourth-ventricle barriers.
    Qp_BBB   <- (Daq / w_BBB) * (SA_BBB * f_para_BBB)
    Qt_BBB   <- 0.5 * P0 * (SA_BBB * f_trans)
    Qp_BCSFB <- (Daq / w_BCSFB) * ((SA_BCSFB / 2) * f_para_BCSFB)
    Qt_BCSFB <- 0.5 * P0 * ((SA_BCSFB / 2) * f_trans)
    # Non-specific phospholipid binding: water-to-oil and oil-to-water
    # clearances share the octanol/water partition ratio 10^logP.
    CL_wo  <- P0 * SA_BCM
    CL_ow  <- CL_wo / (10^logp)
    Q_LYSO <- P0 * SA_LYSO

    # ---- 5. Asymmetry factors (back-solved from Kp,uu) ----
    # Closed forms from the Saleh 2021 supplementary equations, solved for the
    # efflux factor with the influx factor fixed at 1 (valid when Kp,uu < 1).
    AF_BBB_ef <- -(kpuu_ecf * Q_CBF * Q_ECF +
                   kpuu_cm * Q_CSF * (AF_in * PHF_MV * Qt_BBB + Qp_BBB) +
                   Q_CBF * (kpuu_ecf * Qp_BBB - AF_in * PHF_MV * Qt_BBB - Qp_BBB)) /
                 (kpuu_ecf * PHF_ECF * Qt_BBB * Q_CBF)
    AF_LV_ef <- (kpuu_ecf * Q_CBF * Q_ECF -
                 Q_CSF * (kpuu_cm * (AF_in * PHF_MV * Qt_BCSFB + Qp_BCSFB) + kpuu_lv * Q_CBF) +
                 Q_CBF * (-kpuu_lv * Qp_BCSFB + AF_in * PHF_MV * Qt_BCSFB + Qp_BCSFB)) /
               (kpuu_lv * PHF_CSF * Qt_BCSFB * Q_CBF)
    AF_TFV_ef <- (Q_CSF * (Q_CBF * (kpuu_lv - kpuu_cm) -
                           kpuu_cm * (AF_in * PHF_MV * Qt_BCSFB + Qp_BCSFB)) +
                  Q_CBF * (-kpuu_cm * Qp_BCSFB + AF_in * PHF_MV * Qt_BCSFB + Qp_BCSFB)) /
                (kpuu_cm * PHF_CSF * Qt_BCSFB * Q_CBF)

    # ---- 6. CNS compartment concentrations ----
    # Every CNS state holds an amount; each compartment carries the TOTAL
    # unbound concentration, of which the fraction PHF is uncharged.
    C_MVu <- fu_plasma * (brain_vascular / V_MV)
    C_ECF <- brain_ecf / V_ECF
    C_BCM <- brain_cell_membrane / V_BCM
    C_ICF <- brain_icf / V_ICF
    C_LYS <- brain_lysosome / V_LYS
    C_LV  <- brain_csf_lv / V_LV
    C_TFV <- brain_csf_tfv / V_TFV
    C_CM  <- brain_csf_cm / V_CM
    C_SAS <- brain_csf_sas / V_SAS

    # ---- 7. Barrier fluxes ----
    # Neutral drug crosses by both paracellular and transcellular routes;
    # charged drug is restricted to the paracellular route.
    inBBB <- Qp_BBB * C_MVu + Qt_BBB * AF_in * PHF_MV * C_MVu
    efBBB <- Qp_BBB * C_ECF + Qt_BBB * AF_BBB_ef * PHF_ECF * C_ECF
    inLV  <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_in * PHF_MV * C_MVu
    efLV  <- Qp_BCSFB * C_LV + Qt_BCSFB * AF_LV_ef * PHF_CSF * C_LV
    inTFV <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_in * PHF_MV * C_MVu
    efTFV <- Qp_BCSFB * C_TFV + Qt_BCSFB * AF_TFV_ef * PHF_CSF * C_TFV

    # ---- 8. Empirical plasma model (forcing function) ----
    # Published as a forcing function: CNS uptake does NOT deplete plasma,
    # so the plasma compartment is deliberately not mass-balanced with brain.
    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # ---- 9. CNS PBPK system ----
    d/dt(brain_vascular)      <- Q_CBF * (central / vc - brain_vascular / V_MV) -
                                 inBBB + efBBB - inLV + efLV - inTFV + efTFV
    d/dt(brain_ecf)           <- inBBB - efBBB - Q_ECF * C_ECF -
                                 CL_wo * PHF_ECF * C_ECF + CL_ow * C_BCM
    d/dt(brain_cell_membrane) <- CL_wo * PHF_ECF * C_ECF + CL_wo * PHF_ICF * C_ICF -
                                 2 * CL_ow * C_BCM
    d/dt(brain_icf)           <- CL_ow * C_BCM - CL_wo * PHF_ICF * C_ICF -
                                 Q_LYSO * PHF_ICF * C_ICF + Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_lysosome)      <- Q_LYSO * PHF_ICF * C_ICF - Q_LYSO * PHF_LYS * C_LYS
    d/dt(brain_csf_lv)        <- inLV - efLV + Q_ECF * C_ECF - Q_CSF * C_LV
    d/dt(brain_csf_tfv)       <- inTFV - efTFV + Q_CSF * C_LV - Q_CSF * C_TFV
    d/dt(brain_csf_cm)        <- Q_CSF * C_TFV - Q_CSF * C_CM
    d/dt(brain_csf_sas)       <- Q_CSF * C_CM - Q_CSF * C_SAS

    # ---- 10. Observations ----
    # Cc is the TOTAL plasma concentration described by the empirical model;
    # Ccu is the unbound plasma concentration that drives CNS uptake.
    # Cbrain_ecf is the unbound brain-ECF concentration, the paper's endpoint.
    Cc          <- central / vc
    Ccu         <- fu_plasma * (central / vc)
    Cbrain_ecf  <- C_ECF
    Cbrain_icf  <- C_ICF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_cm <- C_CM
    Cc ~ prop(propSd)
  })
}

