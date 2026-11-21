Sun_2025_acyclovir_pbpk <- function() {
  description <- paste(
    "PBPK (LeiCNS3.0 CNS physiologically-based model, human version).",
    "Acyclovir distribution into brain extracellular fluid (brain ECF) and",
    "the four CSF compartments of adults, fitted by full Bayesian analysis",
    "to aggregate plasma and lumbar CSF (subarachnoid space, SAS) data",
    "after oral valacyclovir. Plasma is a one-compartment model with",
    "first-order absorption (oral valacyclovir) that acts as a forcing",
    "function for the LeiCNS3.0 CNS system: brain microvasculature, brain",
    "ECF, a phospholipid brain-cell-membrane compartment, brain",
    "intracellular fluid, lysosomes and the lateral-ventricle, third +",
    "fourth-ventricle, cisterna-magna and SAS CSF compartments. Active",
    "transport enters as efflux asymmetry factors back-calculated from",
    "Kp,uu,ECF = 0.32 and Kp,uu,CSF = 0.23. Parameter values are posterior",
    "means; there is no interindividual variability. The published brain",
    "ECF and plasma simulations are not reproduced by the printed equations",
    "(see the vignette)."
  )
  reference <- paste(
    "Sun M, Manson ML, Martson A-G, Bodilsen J, de Lange ECM, Guo T.",
    "Revisiting Acyclovir Dosing for Adult Viral Encephalitis Using a Full",
    "Bayesian LeiCNS PBPK Modeling Approach. Clin Pharmacokinet.",
    "2025;64(9):1413. doi:10.1007/s40262-025-01545-w.",
    "LeiCNS-PK3.0 structure, the Daq / P0 / Qp / Qt regressions, the",
    "asymmetry-factor closed forms and the brain phospholipid volume:",
    "Saleh MAA, Loo CF, Elassaiss-Schaap J, de Lange ECM. J Pharmacokinet",
    "Pharmacodyn. 2021;48(5):725-741. doi:10.1007/s10928-021-09768-7",
    "(reference [19])."
  )
  vignette <- "Sun_2025_acyclovir_pbpk"

  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "valacyclovir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "acyclovir", units = "mg", specimen = "plasma", verified = TRUE),
    brain_vascular = list(analyte = "acyclovir", units = "mg", specimen = "plasma", verified = TRUE),
    brain_ecf = list(analyte = "acyclovir", units = "mg", specimen = "brain ISF", verified = TRUE),
    brain_cell_membrane = list(analyte = "acyclovir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_icf = list(analyte = "acyclovir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_lysosome = list(analyte = "acyclovir", units = "mg", specimen = "tissue", verified = TRUE),
    brain_csf_lv = list(analyte = "acyclovir", units = "mg", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "acyclovir", units = "mg", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "acyclovir", units = "mg", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "acyclovir", units = "mg", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_studies = 2,
    disease_state = paste(
      "Fitted data: adults with normal renal function (Smith 2010; the",
      "chronic-kidney-disease arm was excluded) and adults with multiple",
      "sclerosis (Lycke 2003). Simulation target: adults with HSV or VZV",
      "viral encephalitis, represented through inflamed-CNS priors on SAS",
      "volume, BBB pore size and CSF flow (Table S2)."
    ),
    dose_range = paste(
      "Fitted: oral valacyclovir 1000 mg TID (Lycke 2003) and 2000 mg QID",
      "(Smith 2010). Simulated: IV acyclovir 8-20 mg/kg every 8 h or every",
      "6 h as 1-, 1.5- or 2-h infusions."
    ),
    notes = paste(
      "Only aggregate (mean) concentrations digitised from the two",
      "publications were available, so the number of subjects is not",
      "reported and no interindividual variability was estimated",
      "(Discussion). Body weight for the mg/kg simulations is not reported."
    )
  )

  ini({
    # ---- Plasma PK: Table S3 posterior means ----
    lcl <- log(567); label("Acyclovir clearance CLe (mL/min)")  # Table S3: CLe mean 567 mL/min
    lvc <- log(62481); label("Acyclovir central volume VCEN (mL)")  # Table S3: VCEN mean 62481 mL
    lka <- log(0.0074); label("Valacyclovir absorption rate constant Ka (1/min)")  # Table S3: Ka mean 0.0074 1/min
    lfdepot <- fixed(log(0.545)); label("Bioavailability of acyclovir from oral valacyclovir (fraction)")  # Table S1: F 54.5% (fixed; not used in the IV simulations)

    # ---- Acyclovir physicochemical properties: Table S1 (DrugBank) ----
    mwt <- fixed(225.2); label("Molecular weight (g/mol)")  # Table S1
    logp <- fixed(-0.95); label("Lipophilicity, logP (unitless)")  # Table S1
    pka <- fixed(11.98); label("Acidic ionization constant (unitless)")  # Table S1
    pkb <- fixed(3.02); label("Basic ionization constant (unitless)")  # Table S1
    # Plasma protein binding is not listed in Table S1 and the CNS
    # simulations follow the total plasma concentration: SAS/plasma AUC in
    # Figure 3 (5.9 / 26.1 = 0.23) equals Kp,uu,CSF only when fu = 1.
    fu_plasma <- fixed(1); label("Unbound fraction in plasma (unitless)")  # Not reported; set to 1 (see comment)

    # ---- Extent of CNS transport: Table S1 ----
    kpuu_ecf <- fixed(0.32); label("Kp,uu,ECF, brain ECF-to-plasma unbound ratio (unitless)")  # Table S1: Kpuu,ECF 0.32
    kpuu_csf <- fixed(0.23); label("Kp,uu,CSF, CSF-to-plasma unbound ratio (unitless)")  # Table S1: Kpuu,CSF 0.23

    # ---- CNS parameters estimated in the Bayesian analysis: Table S3 ----
    V_SAS <- 278; label("Subarachnoid space CSF volume (mL)")  # Table S3: VSAS mean 278 mL
    PZ_BBB <- 0.00145; label("BBB pore size (um)")  # Table S3: PZ_BBB mean 0.00145 um
    SA_BBB <- 152997; label("BBB surface area (cm^2)")  # Table S3: SABBB mean 152997 cm^2
    w_BBB <- 0.0000515; label("BBB width (cm)")  # Table S3: Width_BBB mean 0.0000515 cm
    Q_CSF <- 1.096; label("CSF flow (mL/min)")  # Table S3: QCSF mean 1.096 mL/min

    # Reference BBB pore size against which PZ_BBB scales the paracellular
    # BBB clearance (healthy prior mean).
    PZ_BBB_ref <- fixed(0.0007); label("Reference (healthy) BBB pore size (um)")  # Table S2: PZ_BBB healthy mean 0.0007 um

    # ---- Fixed human CNS physiology: Table S1 (from Saleh 2021) ----
    V_MV <- fixed(47); label("Brain microvasculature volume (mL)")  # Table S1: VMV 47 mL
    V_ECF <- fixed(253); label("Brain ECF volume (mL)")  # Table S1: VECF 253 mL
    V_ICF <- fixed(1000); label("Brain ICF volume (mL)")  # Table S1: VICF 1000 mL
    V_LYS <- fixed(12.5); label("Total lysosome volume (mL)")  # Table S1: VLYS 12.5 mL
    V_LV <- fixed(20); label("Lateral ventricle CSF volume (mL)")  # Table S1: VLV 20 mL
    V_TFV <- fixed(3); label("Third + fourth ventricle CSF volume (mL)")  # Table S1: VTFV 3 mL
    V_CM <- fixed(1); label("Cisterna magna CSF volume (mL)")  # Table S1: VCM 1 mL
    # Brain phospholipid (cell-membrane) volume: 5% of total brain volume
    # (Saleh 2021); total brain volume 1250 mL from Table S1.
    V_BCM <- fixed(0.05 * 1250); label("Brain phospholipid volume (mL)")  # Saleh 2021: 5% of Vtot; Table S1: Vtot 1250 mL

    Q_CBF <- fixed(689); label("Cerebral blood flow (mL/min)")  # Table S1: QCBF 689 mL/min
    Q_ECF <- fixed(0.2); label("Brain ECF bulk flow (mL/min)")  # Table S1: QECF 0.2 mL/min

    SA_BCSFB <- fixed(15000); label("BCSFB total surface area (cm^2)")  # Table S1: SABCSFB 15000 cm^2
    SA_BCM <- fixed(2666520); label("Brain cell membrane surface area (cm^2)")  # Table S1: SABCM 2666520 cm^2
    SA_LYSO <- fixed(1980260); label("Lysosomal membrane surface area (cm^2)")  # Table S1: SALYS 1980260 cm^2

    # Effective surface areas are printed as percentages.
    f_trans <- fixed(99.8 / 100); label("Transcellular effective surface-area fraction (unitless)")  # Table S1: SABBB,T and SABCSFB,T 99.8%
    f_para_BBB <- fixed(0.004 / 100); label("Paracellular effective surface-area fraction, BBB (unitless)")  # Table S1: SABBB,P 0.004%
    f_para_BCSFB <- fixed(0.016 / 100); label("Paracellular effective surface-area fraction, BCSFB (unitless)")  # Table S1: SABCSFB,P 0.016%

    w_BCSFB <- fixed(0.5 * 1e-4); label("BCSFB width (cm)")  # Table S1: WBCSFB 0.5 um

    pH_PL <- fixed(7.4); label("Plasma and brain microvasculature pH (unitless)")  # Table S1: pHPL and pHMV 7.4
    pH_ECF <- fixed(7.3); label("Brain ECF pH (unitless)")  # Table S1: pHECF 7.3
    pH_ICF <- fixed(7); label("Brain ICF pH (unitless)")  # Table S1: pHICF 7
    pH_LYS <- fixed(5); label("Lysosomal pH (unitless)")  # Table S1: pHLYS 5
    pH_CSF <- fixed(7.3); label("CSF pH (unitless)")  # Table S1: pHCSF 7.3

    # ---- Residual error: Table S3 (log-normal likelihood) ----
    expSd <- 0.301; label("Log-normal residual SD, plasma acyclovir (log scale)")  # Table S3: Sigma_Plasma mean 0.301
    expSd_Cbrain_csf_sas <- 0.137; label("Log-normal residual SD, lumbar CSF (SAS) acyclovir (log scale)")  # Table S3: Sigma_SAS mean 0.137
  })

  model({
    # ---- 1. Plasma PK ----
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- cl / vc

    # Both Kp,uu,CSF sites (lateral ventricle and cisterna magna) take the
    # single reported Kp,uu,CSF.
    kpuu_lv <- kpuu_csf
    kpuu_cm <- kpuu_csf

    # ---- 2. pH-dependent neutral (uncharged) fractions ----
    # Acidic term 1/(1 + 10^(pH - pKa)); basic term 1/(1 + 10^(pKb - pH)).
    PHF_MV <- (1 / (1 + 10^(pH_PL - pka))) * (1 / (1 + 10^(pkb - pH_PL)))
    PHF_ECF <- (1 / (1 + 10^(pH_ECF - pka))) * (1 / (1 + 10^(pkb - pH_ECF)))
    PHF_ICF <- (1 / (1 + 10^(pH_ICF - pka))) * (1 / (1 + 10^(pkb - pH_ICF)))
    PHF_LYS <- (1 / (1 + 10^(pH_LYS - pka))) * (1 / (1 + 10^(pkb - pH_LYS)))
    PHF_CSF <- (1 / (1 + 10^(pH_CSF - pka))) * (1 / (1 + 10^(pkb - pH_CSF)))

    # ---- 3. Passive diffusivity and permeability (Saleh 2021) ----
    # Both regressions return per-second units; x 60 converts to per-minute.
    Daq <- 10^(-4.113 - 0.4609 * log10(mwt)) * 60
    P0 <- 10^(0.939 * logp - 6.21) * 60

    # ---- 4. Barrier and intracellular clearances ----
    # The paracellular BBB clearance is scaled linearly by the estimated pore
    # size relative to its healthy value; the printed LeiCNS-PK3.0 equation
    # has no pore-size term. The BCSFB pore size is fixed at its Table S1
    # value (0.027 um), so its scaling factor is 1. The total BCSFB surface
    # area is split 50/50 between the lateral-ventricle and third +
    # fourth-ventricle barriers.
    Qp_BBB <- (Daq / w_BBB) * (SA_BBB * f_para_BBB) * (PZ_BBB / PZ_BBB_ref)
    Qt_BBB <- 0.5 * P0 * (SA_BBB * f_trans)
    Qp_BCSFB <- (Daq / w_BCSFB) * ((SA_BCSFB / 2) * f_para_BCSFB)
    Qt_BCSFB <- 0.5 * P0 * ((SA_BCSFB / 2) * f_trans)
    CL_wo <- P0 * SA_BCM
    CL_ow <- CL_wo / (10^logp)
    Q_LYSO <- P0 * SA_LYSO

    # ---- 5. Asymmetry factors (Saleh 2021 closed forms) ----
    # All Kp,uu values are below 1, so the influx factors are 1 and the
    # efflux factors are back-calculated from Kp,uu.
    AF_BBB_in <- 1
    AF_BBB_ef <- -(kpuu_ecf * Q_CBF * Q_ECF +
      kpuu_cm * Q_CSF * (PHF_MV * Qt_BBB + Qp_BBB) +
      Q_CBF * (kpuu_ecf * Qp_BBB - PHF_MV * Qt_BBB - Qp_BBB)) /
      (kpuu_ecf * PHF_ECF * Qt_BBB * Q_CBF)
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

    # ---- 6. CNS compartment concentrations (mg/mL) ----
    C_MVu <- fu_plasma * (brain_vascular / V_MV)
    C_ECF <- brain_ecf / V_ECF
    C_BCM <- brain_cell_membrane / V_BCM
    C_ICF <- brain_icf / V_ICF
    C_LYS <- brain_lysosome / V_LYS
    C_LV <- brain_csf_lv / V_LV
    C_TFV <- brain_csf_tfv / V_TFV
    C_CM <- brain_csf_cm / V_CM
    C_SAS <- brain_csf_sas / V_SAS

    # ---- 7. Barrier fluxes ----
    # Neutral drug crosses by both routes; charged drug only paracellularly.
    inBBB <- Qp_BBB * C_MVu + Qt_BBB * AF_BBB_in * PHF_MV * C_MVu
    efBBB <- Qp_BBB * C_ECF + Qt_BBB * AF_BBB_ef * PHF_ECF * C_ECF
    inLV <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_LV_in * PHF_MV * C_MVu
    efLV <- Qp_BCSFB * C_LV + Qt_BCSFB * AF_LV_ef * PHF_CSF * C_LV
    inTFV <- Qp_BCSFB * C_MVu + Qt_BCSFB * AF_TFV_in * PHF_MV * C_MVu
    efTFV <- Qp_BCSFB * C_TFV + Qt_BCSFB * AF_TFV_ef * PHF_CSF * C_TFV

    # ---- 8. Plasma model (forcing function) ----
    # CNS uptake does not deplete plasma (LeiCNS-PK3.0 plasma input).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    # ---- 9. CNS PBPK system ----
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

    # ---- 10. Observations (mg/L) ----
    # Cc: total plasma acyclovir; Cbrain_csf_sas: lumbar CSF (SAS), the two
    # fitted endpoints. Cbrain_ecf: brain ECF, the paper's efficacy site.
    Cc <- 1000 * central / vc
    Cbrain_ecf <- 1000 * C_ECF
    Cbrain_csf_sas <- 1000 * C_SAS
    Cc ~ lnorm(expSd)
    Cbrain_csf_sas ~ lnorm(expSd_Cbrain_csf_sas)
  })
}
