Fornari_2019_carboplatin_rat <- function() {
  description <- paste(
    "Preclinical (rat). QSP. Haematopoiesis systems pharmacology model of",
    "carboplatin-induced bone marrow toxicity in male RccHan:WIST rats",
    "(Fornari 2019). A two-compartment free-platinum PK model drives a DNA",
    "damage state (adduct formation proportional to free platinum, first-order",
    "repair). Damage kills the proliferating fraction of the bone marrow",
    "multipotent progenitors (MPP), common myeloid progenitors (CMP) and",
    "megakaryocyte-erythrocyte progenitors (MEP), and of the first two",
    "proliferative transit compartments of the reticulocyte and platelet",
    "lineages. MPP feed CMP and MEP; CMP feed three-compartment non-proliferative",
    "neutrophil and monocyte maturation chains; MEP feed the reticulocyte and",
    "platelet chains (two proliferative plus one non-proliferative transit",
    "compartment); circulating neutrophils, monocytes, reticulocytes, platelets",
    "and red blood cells close the system. Negative feedback from CMP/MEP on",
    "the stem-cell inflow, from circulating cells on MPP differentiation into",
    "CMP and MEP, and from circulating reticulocytes, red cells and platelets",
    "on transit-cell proliferation regulate recovery. All differentiation and",
    "proliferation rates are derived from the untreated steady state, so every",
    "state starts and stays at its baseline without drug. Naive-pooled fit:",
    "no IIV; log-additive residual error per endpoint.",
    sep = " "
  )
  reference <- paste(
    "Fornari C, Oplustil O'Connor L, Pin C, Smith A, Yates JWT, Cheung SYA,",
    "Jodrell DI, Mettetal JT, Collins TA. Quantifying Drug-Induced Bone",
    "Marrow Toxicity Using a Novel Haematopoiesis Systems Pharmacology Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(11):858-868.",
    "doi:10.1002/psp4.12459. PMID: 31508894. PMCID: PMC6875710.",
    "Equations 1-13 and Tables 1-2 of the main text; Supplementary Text S5",
    "(rat PK model), S6-S7 (steady-state derivation and reparameterisation);",
    "Tables S3 (rat PK) and S4 (residual error); Supplementary Model Code",
    "(MATLAB parameter block).",
    sep = " "
  )
  vignette <- "Fornari_2019_carboplatin_myelosuppression"
  units <- list(
    time = "h",
    dosing = "mg/kg (carboplatin IV bolus; the PK volumes and clearances act per kg of body weight, see vignette)",
    concentration = "uM (free platinum, Cc); 10^9 cells/L (all cell-count outputs)"
  )

  dosing <- c("central")

  # Bone marrow progenitors, lineage-specific transit chains and circulating
  # pools of the Fornari 2019 haematopoiesis tree (Figure 1) plus the DNA
  # damage state of Equation 3. None is a reusable library compartment.
  paper_specific_compartments <- c(
    "damage",
    "mpp",
    "cmp",
    "mep",
    "precursor1_mono",
    "precursor2_mono",
    "precursor3_mono",
    "precursor1_ret",
    "precursor2_ret",
    "precursor3_ret",
    "circ_mono",
    "circ_ret",
    "circ_rbc"
  )

  compartmentData <- list(
    central = list(analyte = "carboplatin (free platinum)", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "carboplatin (free platinum)", units = "mg/kg", specimen = "plasma", verified = TRUE),
    damage = list(
      analyte = "carboplatin-induced DNA damage (platinum-DNA adducts)",
      units = "uM-equivalent",
      specimen = "not applicable",
      verified = TRUE
    ),
    mpp = list(
      analyte = "multipotent progenitors (CD90+/Lin- cells)",
      units = "10^9 cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    cmp = list(
      analyte = "common myeloid progenitors (CD45+ cells)",
      units = "10^9 cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    mep = list(
      analyte = "megakaryocyte-erythrocyte progenitors (CD71+ cells)",
      units = "10^9 cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor1_anc = list(
      analyte = "neutrophil transit cells (T1Neut)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2_anc = list(
      analyte = "neutrophil transit cells (T2Neut)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor3_anc = list(
      analyte = "neutrophil transit cells (T3Neut)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor1_mono = list(
      analyte = "monocyte transit cells (T1Mono)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2_mono = list(
      analyte = "monocyte transit cells (T2Mono)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor3_mono = list(
      analyte = "monocyte transit cells (T3Mono)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor1_ret = list(
      analyte = "proliferative reticulocyte transit cells (T1Ret)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2_ret = list(
      analyte = "proliferative reticulocyte transit cells (T2Ret)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor3_ret = list(
      analyte = "reticulocyte transit cells (T3Ret)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor1_plt = list(
      analyte = "proliferative platelet transit cells (T1Plt)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor2_plt = list(
      analyte = "proliferative platelet transit cells (T2Plt)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    precursor3_plt = list(
      analyte = "platelet transit cells (T3Plt)",
      units = "10^9 cells/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ_anc = list(
      analyte = "circulating neutrophils",
      units = "10^9 cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    circ_mono = list(
      analyte = "circulating monocytes",
      units = "10^9 cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    circ_ret = list(
      analyte = "circulating reticulocytes",
      units = "10^9 cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    circ_plt = list(
      analyte = "circulating platelets",
      units = "10^9 cells/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    circ_rbc = list(
      analyte = "circulating red blood cells",
      units = "10^9 cells/L",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "rat (RccHan:WIST, male)",
    n_subjects = 178L,
    n_studies = 6L,
    age_range = "approximately 12 weeks at the start of dosing",
    weight_range = "minimum 320 g at the start of dosing",
    sex_female_pct = 0,
    disease_state = "Healthy rats (haematotoxicity study).",
    dose_range = "Carboplatin 3, 10, 30 or 40 mg/kg IV (1-min tail-vein injection), single dose or every 14 days for up to 7 cycles (Table S1 studies A-F); vehicle controls.",
    regions = "AstraZeneca, United Kingdom",
    notes = "Six studies (A-F, Table S1) pooled: haematology on the Siemens Advia 2120i; bone marrow MPP (CD90+/Lin-), CMP (CD45+) and MEP (CD71+) by flow cytometry of femur marrow; free platinum by UPLC-MS/MS of plasma ultrafiltrate (30 mg/kg in-house study) plus literature data (20 mg/kg, Siddik 1987). PK and PD were fitted sequentially by naive pooling (Supplementary Text S8)."
  )

  ini({
    # Free-platinum two-compartment PK (Supplementary Text S5, Equation S1).
    # Table S3 prints CL and Q in '1/h' and volumes in 'l'; the
    # Supplementary Model Code labels CL and Q 'liter/hour'. With doses in
    # mg/kg the volumes act per kg of body weight (see vignette).
    lcl <- log(0.42); label("Free-platinum clearance (L/h, per kg)") # Table S3 'CL' = 0.42 (CV 6%); model code Cl_carbo = 0.42
    lvc <- log(0.26); label("Free-platinum central volume (L, per kg)") # Table S3 'V Cen' = 0.26 (CV 14%)
    lvp <- log(0.30); label("Free-platinum peripheral volume (L, per kg)") # Table S3 'V Per' = 0.30 (CV 9%)
    lq <- log(0.06); label("Free-platinum intercompartmental clearance (L/h, per kg)") # Table S3 'Q C' = 0.06 (CV 9%)

    # DNA damage (Equation 3)
    lkdam <- fixed(log(0.017)); label("DNA-damage formation rate constant kdam (1/h)") # Table 1 'k dam' = 0.017, Fixed
    lkrep <- fixed(log(0.017)); label("DNA-damage repair rate constant krep (1/h)") # Table 1 'k rep' = 0.017, Fixed

    # Circulating (elimination) rates of mature blood cells, literature values
    lkcirc_anc <- fixed(log(0.17)); label("Neutrophil circulating rate kcircNeut (1/h)") # Table 1 'kcircNeut' = 0.17, literature
    lkcirc_mono <- fixed(log(0.06)); label("Monocyte circulating rate kcircMono (1/h)") # Table 1 'kcircMono' = 0.06, literature
    lkcirc_plt <- fixed(log(0.01)); label("Platelet circulating rate kcircPlt (1/h)") # Table 1 'kcircPlt' = 0.01, literature
    lkcirc_rbc <- fixed(log(0.0007)); label("Red blood cell circulating rate kcircRBC (1/h)") # Table 1 'kcircRBC' = 0.0007, literature

    # Baseline (homeostatic) cell counts
    lbl_anc <- log(1.32); label("Neutrophil baseline Neut0 (10^9 cells/L)") # Table 1 'Neut0' = 1.32 (CV 2.6%)
    lbl_mono <- log(0.16); label("Monocyte baseline Mono0 (10^9 cells/L)") # Table 1 'Mono0' = 0.16 (CV 2.6%)
    lbl_plt <- log(967); label("Platelet baseline Plt0 (10^9 cells/L)") # Table 1 'Plt0' = 967 (CV 2.5%)
    lbl_ret <- log(235); label("Reticulocyte baseline Ret0 (10^9 cells/L)") # Table 1 'Ret0' = 235 (CV 2.7%)
    lbl_rbc <- log(8161); label("Red blood cell baseline RBC0 (10^9 cells/L)") # Table 1 'RBC0' = 8,161 (CV 2.2%)
    lbl_mpp <- log(69); label("Multipotent progenitor baseline MPP0 (10^9 cells/L)") # Table 1 'MPP0' = 69 (CV 4.6%)
    lbl_cmp <- log(214); label("Common myeloid progenitor baseline CMP0 (10^9 cells/L)") # Table 1 'CMP0' = 214 (CV 4.1%)
    lbl_mep <- log(156); label("Megakaryocyte-erythrocyte progenitor baseline MEP0 (10^9 cells/L)") # Table 1 'MEP0' = 156 (CV 5.2%)

    # Mean transit (maturation) times of the three-compartment chains
    lmtt_anc <- log(61.5); label("Neutrophil mean transit time MTTNeut (h)") # Table 1 'MTTNeut' = 61.5 (CV 56%)
    lmtt_mono <- log(62.2); label("Monocyte mean transit time MTTMono (h)") # Table 1 'MTTMono' = 62.2 (CV 55%)
    lmtt_plt <- log(74); label("Platelet mean transit time MTTPlt (h)") # Table 1 'MTTPlt' = 74 (CV 5.6%)
    lmtt_ret <- log(75); label("Reticulocyte mean transit time MTTRet (h)") # Table 1 'MTTRet' = 75 (CV 4.8%)

    # Steady-state reparameterisation ratios (Supplementary Text S7, Equation S7)
    lambda_t2t1_ret <- fixed(2); label("Ratio T2Ret0/T1Ret0 of the reticulocyte transit compartments, lambda1 (dimensionless)") # Table 1 'lambda1' = 2, Fixed
    lambda_t2t1_plt <- fixed(2); label("Ratio T2Plt0/T1Plt0 of the platelet transit compartments, lambda2 (dimensionless)") # Table 1 'lambda2' = 2, Fixed
    lambda_cmp <- fixed(1.8); label("Ratio of CMP differentiation to proliferation rate, lambda3 (dimensionless)") # Table 1 'lambda3' = 1.8, Fixed
    lambda_mep <- fixed(1.8); label("Ratio of MEP differentiation to proliferation rate, lambda4 (dimensionless)") # Table 1 'lambda4' = 1.8, Fixed
    lambda_mpp <- fixed(1.8); label("Ratio of MPP differentiation to proliferation rate, lambda5 (dimensionless)") # Table 1 'lambda5' = 1.8, Fixed

    # Feedback powers (Equations 11-13)
    gamma_stem <- fixed(0.07); label("Feedback power on the stem-cell inflow into MPP (dimensionless)") # Table 1 'gamma stem' = 0.07, Fixed
    gamma_mat_cmp <- fixed(0.6); label("Feedback power on MPP differentiation into CMP (dimensionless)") # Table 1 'gamma MatCMP' = 0.6, Fixed
    gamma_mat_mep <- fixed(0.3); label("Feedback power on MPP differentiation into MEP (dimensionless)") # Table 1 'gamma MatMEP' = 0.3, Fixed
    gamma_prol_trans <- fixed(0.7); label("Feedback power on reticulocyte and platelet transit-cell proliferation (dimensionless)") # Table 1 'gamma prolTrans' = 0.7, Fixed

    # Carboplatin cell-kill parameters (Equation 2, Table 2)
    slope_mpp <- 2.05; label("Carboplatin damage-driven kill slope on MPP (1/uM)") # Table 2 'SlopeMPP' = 2.05 (CV 11.8%)
    slope_cmp <- 1.47; label("Carboplatin damage-driven kill slope on CMP (1/uM)") # Table 2 'SlopeCMP' = 1.47 (CV 8.6%)
    slope_mep <- 2.19; label("Carboplatin damage-driven kill slope on MEP (1/uM)") # Table 2 'SlopeMEP' = 2.19 (CV 10.3%)
    delta_ret <- 2.8; label("Multiplier of SlopeMEP for the proliferative reticulocyte transit cells (dimensionless)") # Table 2 'delta Ret' = 2.8 (CV 21%)
    delta_plt <- 0.54; label("Multiplier of SlopeMEP for the proliferative platelet transit cells (dimensionless)") # Table 2 'delta Plt' = 0.54 (CV 30.4%)

    # Log-additive residual error (Table S4)
    expSd_MPP <- 0.33; label("Log-additive residual SD, MPP") # Table S4 'sigma MPP' = 0.33
    expSd_CMP <- 0.19; label("Log-additive residual SD, CMP") # Table S4 'sigma CMP' = 0.19
    expSd_MEP <- 0.33; label("Log-additive residual SD, MEP") # Table S4 'sigma MEP' = 0.33
    expSd_ANC <- 0.17; label("Log-additive residual SD, neutrophils") # Table S4 'sigma Neut' = 0.17
    expSd_MONO <- 0.17; label("Log-additive residual SD, monocytes") # Table S4 'sigma Mono' = 0.17
    expSd_PLT <- 0.17; label("Log-additive residual SD, platelets") # Table S4 'sigma Ptl' = 0.17
    expSd_RET <- 0.36; label("Log-additive residual SD, reticulocytes") # Table S4 'sigma Ret' = 0.36
    expSd_RBC <- 0.06; label("Log-additive residual SD, red blood cells") # Table S4 'sigma RBC' = 0.06
  })

  model({
    # Individual parameters (naive-pooled fit, no IIV)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    kdam <- exp(lkdam)
    krep <- exp(lkrep)
    kcirc_anc <- exp(lkcirc_anc)
    kcirc_mono <- exp(lkcirc_mono)
    kcirc_plt <- exp(lkcirc_plt)
    kcirc_rbc <- exp(lkcirc_rbc)
    bl_anc <- exp(lbl_anc)
    bl_mono <- exp(lbl_mono)
    bl_plt <- exp(lbl_plt)
    bl_ret <- exp(lbl_ret)
    bl_rbc <- exp(lbl_rbc)
    bl_mpp <- exp(lbl_mpp)
    bl_cmp <- exp(lbl_cmp)
    bl_mep <- exp(lbl_mep)

    # Maturation rates within each transit chain: a = 3 / MTT (Equations 6, 8)
    a_anc <- 3 / exp(lmtt_anc)
    a_mono <- 3 / exp(lmtt_mono)
    a_ret <- 3 / exp(lmtt_ret)
    a_plt <- 3 / exp(lmtt_plt)

    # Steady-state relations of the untreated system (Supplementary Text S6,
    # Equations S3-S4, with the lambda reparameterisation of Equation S7).
    # Reticulocytes leave circulation at the rate that keeps the red-cell
    # pool at RBC0.
    kcirc_ret <- kcirc_rbc * bl_rbc / bl_ret
    t_anc0 <- kcirc_anc * bl_anc / a_anc
    t_mono0 <- kcirc_mono * bl_mono / a_mono
    t_ret0 <- kcirc_ret * bl_ret / a_ret
    t1_ret0 <- t_ret0 / lambda_t2t1_ret
    t_plt0 <- kcirc_plt * bl_plt / a_plt
    t1_plt0 <- t_plt0 / lambda_t2t1_plt
    ktr_anc <- kcirc_anc * bl_anc / bl_cmp
    ktr_mono <- kcirc_mono * bl_mono / bl_cmp
    kprol_cmp <- (ktr_anc + ktr_mono) / lambda_cmp
    ktr_ret <- a_ret * t1_ret0^2 / (bl_mep * t_ret0)
    ktr_plt <- a_plt * t1_plt0^2 / (bl_mep * t_plt0)
    kprol_mep <- (ktr_ret + ktr_plt) / lambda_mep
    ktr_cmp <- (ktr_anc + ktr_mono - kprol_cmp) * bl_cmp / bl_mpp
    ktr_mep <- (ktr_ret + ktr_plt - kprol_mep) * bl_mep / bl_mpp
    kprol_mpp <- (ktr_cmp + ktr_mep) / lambda_mpp
    kprol_ret <- a_ret * (1 - t1_ret0 / t_ret0)
    kprol_plt <- a_plt * (1 - t1_plt0 / t_plt0)
    kstem <- (ktr_cmp + ktr_mep - kprol_mpp) * bl_mpp

    # Free platinum in uM: dose in mg/kg over volume in L/kg gives mg/L,
    # converted with the carboplatin molar mass (371.25 g/mol, one platinum
    # atom per molecule).
    Cc <- central / vc * 1000 / 371.25

    # Feedback functions (Equations 11-13). The 1/2 weights are 'alphaWBC =
    # 0.5' in the Supplementary Model Code.
    fdbk_stem <- (0.5 * bl_cmp / cmp + 0.5 * bl_mep / mep)^gamma_stem
    fdbk_mat_cmp <- (0.5 * bl_anc / circ_anc + 0.5 * bl_mono / circ_mono)^gamma_mat_cmp
    fdbk_mat_mep <- (bl_rbc / circ_rbc)^gamma_mat_mep
    fdbk_prol_ret <- (0.5 * bl_ret / circ_ret + 0.5 * bl_rbc / circ_rbc)^gamma_prol_trans
    fdbk_prol_plt <- (bl_plt / circ_plt)^gamma_prol_trans

    # Free-platinum PK (Equation S1)
    d/dt(central) <- -cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # DNA damage (Equation 3)
    d/dt(damage) <- kdam * Cc - krep * damage

    # Progenitors (Equations 1 and 4)
    d/dt(mpp) <- kstem * fdbk_stem + kprol_mpp * (1 - slope_mpp * damage) * mpp -
      ktr_cmp * fdbk_mat_cmp * mpp - ktr_mep * fdbk_mat_mep * mpp
    d/dt(cmp) <- kprol_cmp * (1 - slope_cmp * damage) * cmp + ktr_cmp * fdbk_mat_cmp * mpp -
      ktr_anc * cmp - ktr_mono * cmp
    d/dt(mep) <- kprol_mep * (1 - slope_mep * damage) * mep + ktr_mep * fdbk_mat_mep * mpp -
      ktr_ret * mep - ktr_plt * mep

    # Non-proliferative neutrophil and monocyte transit chains (Equation 5)
    d/dt(precursor1_anc) <- ktr_anc * cmp - a_anc * precursor1_anc
    d/dt(precursor2_anc) <- a_anc * precursor1_anc - a_anc * precursor2_anc
    d/dt(precursor3_anc) <- a_anc * precursor2_anc - a_anc * precursor3_anc
    d/dt(precursor1_mono) <- ktr_mono * cmp - a_mono * precursor1_mono
    d/dt(precursor2_mono) <- a_mono * precursor1_mono - a_mono * precursor2_mono
    d/dt(precursor3_mono) <- a_mono * precursor2_mono - a_mono * precursor3_mono

    # Reticulocyte and platelet chains: two proliferative, carboplatin-
    # sensitive transit compartments then one non-proliferative (Equation 7).
    # The T1Plt outflow uses aPlt (Equation 7 prints aRet; see vignette).
    d/dt(precursor1_ret) <- kprol_ret * fdbk_prol_ret * precursor1_ret -
      delta_ret * slope_mep * kprol_ret * damage * precursor1_ret +
      ktr_ret * mep - a_ret * precursor1_ret
    d/dt(precursor2_ret) <- kprol_ret * fdbk_prol_ret * precursor2_ret -
      delta_ret * slope_mep * kprol_ret * damage * precursor2_ret +
      a_ret * precursor1_ret - a_ret * precursor2_ret
    d/dt(precursor3_ret) <- a_ret * precursor2_ret - a_ret * precursor3_ret
    d/dt(precursor1_plt) <- kprol_plt * fdbk_prol_plt * precursor1_plt -
      delta_plt * slope_mep * kprol_plt * damage * precursor1_plt +
      ktr_plt * mep - a_plt * precursor1_plt
    d/dt(precursor2_plt) <- kprol_plt * fdbk_prol_plt * precursor2_plt -
      delta_plt * slope_mep * kprol_plt * damage * precursor2_plt +
      a_plt * precursor1_plt - a_plt * precursor2_plt
    d/dt(precursor3_plt) <- a_plt * precursor2_plt - a_plt * precursor3_plt

    # Circulating cells (Equation 9). Equation 9 omits the red-cell balance;
    # reticulocytes mature into red cells at kcircRet, the rate defined by
    # the RBC steady state in Equation S4.
    d/dt(circ_anc) <- a_anc * precursor3_anc - kcirc_anc * circ_anc
    d/dt(circ_mono) <- a_mono * precursor3_mono - kcirc_mono * circ_mono
    d/dt(circ_ret) <- a_ret * precursor3_ret - kcirc_ret * circ_ret
    d/dt(circ_plt) <- a_plt * precursor3_plt - kcirc_plt * circ_plt
    d/dt(circ_rbc) <- kcirc_ret * circ_ret - kcirc_rbc * circ_rbc

    # Homeostatic initial conditions (Equations S2-S3)
    mpp(0) <- bl_mpp
    cmp(0) <- bl_cmp
    mep(0) <- bl_mep
    precursor1_anc(0) <- t_anc0
    precursor2_anc(0) <- t_anc0
    precursor3_anc(0) <- t_anc0
    precursor1_mono(0) <- t_mono0
    precursor2_mono(0) <- t_mono0
    precursor3_mono(0) <- t_mono0
    precursor1_ret(0) <- t1_ret0
    precursor2_ret(0) <- t_ret0
    precursor3_ret(0) <- t_ret0
    precursor1_plt(0) <- t1_plt0
    precursor2_plt(0) <- t_plt0
    precursor3_plt(0) <- t_plt0
    circ_anc(0) <- bl_anc
    circ_mono(0) <- bl_mono
    circ_ret(0) <- bl_ret
    circ_plt(0) <- bl_plt
    circ_rbc(0) <- bl_rbc

    # Observations (Table S4 log-additive residual error)
    MPP <- mpp
    CMP <- cmp
    MEP <- mep
    ANC <- circ_anc
    MONO <- circ_mono
    PLT <- circ_plt
    RET <- circ_ret
    RBC <- circ_rbc
    MPP ~ lnorm(expSd_MPP)
    CMP ~ lnorm(expSd_CMP)
    MEP ~ lnorm(expSd_MEP)
    ANC ~ lnorm(expSd_ANC)
    MONO ~ lnorm(expSd_MONO)
    PLT ~ lnorm(expSd_PLT)
    RET ~ lnorm(expSd_RET)
    RBC ~ lnorm(expSd_RBC)
  })
}
