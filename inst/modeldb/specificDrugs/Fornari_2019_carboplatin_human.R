Fornari_2019_carboplatin_human <- function() {
  description <- paste(
    "QSP. Human translation of the Fornari 2019 haematopoiesis systems",
    "pharmacology model of carboplatin-induced bone marrow toxicity",
    "(prediction, not a fit to human data). Same structure as",
    "Fornari_2019_carboplatin_rat: a two-compartment free-carboplatin PK model",
    "drives a DNA damage state that kills the proliferating fraction of bone",
    "marrow multipotent (MPP), common myeloid (CMP) and megakaryocyte-",
    "erythrocyte (MEP) progenitors and of the proliferative reticulocyte and",
    "platelet transit cells, with neutrophil, monocyte, reticulocyte,",
    "platelet and red-cell lineages regulated by cross-lineage negative",
    "feedback. Physiology-dependent parameters (circulating rates, baseline",
    "counts, maturation times, allometrically scaled progenitor pool sizes)",
    "are the human values of Table 1; drug-sensitivity slopes are the rat",
    "estimates scaled by in vitro colony-forming-unit IC50 ratios (Table 2);",
    "feedback powers, steady-state ratios and DNA damage rates are as in the",
    "rat. Human free-carboplatin PK is the Zandvliet 2008 typical-value model",
    "at a creatinine clearance of 125 mL/min. No IIV.",
    sep = " "
  )
  reference <- paste(
    "Fornari C, Oplustil O'Connor L, Pin C, Smith A, Yates JWT, Cheung SYA,",
    "Jodrell DI, Mettetal JT, Collins TA. Quantifying Drug-Induced Bone",
    "Marrow Toxicity Using a Novel Haematopoiesis Systems Pharmacology Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(11):858-868.",
    "doi:10.1002/psp4.12459. PMID: 31508894. PMCID: PMC6875710.",
    "Human columns of Tables 1-2; Equation 10; Supplementary Text S9-S11 and",
    "Table S2. Human carboplatin PK from Zandvliet AS et al. Br J Clin",
    "Pharmacol. 2008;66(4):485-497. doi:10.1111/j.1365-2125.2008.03230.x",
    "(Table 2; packaged as Zandvliet_2008_carboplatin).",
    sep = " "
  )
  vignette <- "Fornari_2019_carboplatin_myelosuppression"
  units <- list(
    time = "h",
    dosing = "mg (carboplatin IV)",
    concentration = "uM (free carboplatin, Cc); 10^9 cells/L (all cell-count outputs)"
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
    central = list(analyte = "carboplatin (free, ultrafilterable)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(
      analyte = "carboplatin (free, ultrafilterable)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
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
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "Adult cancer patients receiving carboplatin (simulated translation; no human data were fitted).",
    dose_range = "Carboplatin dosed by the Calvert formula to a target AUC(0-24) of 5 mg*min/mL with GFR 125 mL/min (750 mg), two 21-day cycles (Figure 4, Supplementary Text S11).",
    renal_function = "GFR / creatinine clearance 125 mL/min (normal renal function) assumed for dosing and PK.",
    notes = "Projection from the rat model (Fornari_2019_carboplatin_rat) by replacing physiology-dependent parameters with human literature or allometrically scaled values and scaling drug-effect slopes by species-specific in vitro IC50 ratios. Predictions were compared in Figure 4 with simulations of the Schmitt 2010 carboplatin myelosuppression model (Supplementary Text S12)."
  )

  ini({
    # Free-carboplatin two-compartment PK: Zandvliet 2008 Table 2 typical
    # values (Supplementary Text S11 cites this model). Zandvliet CL =
    # 0.76 * CLcr + 1.5 L/h; with the CLcr of 125 mL/min (7.5 L/h) that
    # Supplementary Text S11 assumes for the Calvert dose, CL = 7.2 L/h.
    lcl <- fixed(log(7.2)); label("Free-carboplatin clearance (L/h)") # Zandvliet 2008 Table 2 theta1 = 0.76 x CLcr 7.5 L/h (Fornari S11 GFR 125 mL/min) + theta2 1.5 L/h
    lvc <- fixed(log(15.5)); label("Free-carboplatin central volume (L)") # Zandvliet 2008 Table 2 V1 = 15.5 L
    lvp <- fixed(log(9.86)); label("Free-carboplatin peripheral volume (L)") # Zandvliet 2008 Table 2 V2 = 9.86 L
    lq <- fixed(log(3.46)); label("Free-carboplatin intercompartmental clearance (L/h)") # Zandvliet 2008 Table 2 Q = 3.46 L/h

    # DNA damage (Equation 3), as in the rat
    lkdam <- fixed(log(0.017)); label("DNA-damage formation rate constant kdam (1/h)") # Table 1 human 'k dam' = 0.017, as in the rat
    lkrep <- fixed(log(0.017)); label("DNA-damage repair rate constant krep (1/h)") # Table 1 human 'k rep' = 0.017, as in the rat

    # Circulating (elimination) rates of mature blood cells, literature values
    lkcirc_anc <- fixed(log(0.1)); label("Neutrophil circulating rate kcircNeut (1/h)") # Table 1 human 'kcircNeut' = 0.1
    lkcirc_mono <- fixed(log(0.04)); label("Monocyte circulating rate kcircMono (1/h)") # Table 1 human 'kcircMono' = 0.04
    lkcirc_plt <- fixed(log(0.0052)); label("Platelet circulating rate kcircPlt (1/h)") # Table 1 human 'kcircPlt' = 0.0052
    lkcirc_rbc <- fixed(log(0.00037)); label("Red blood cell circulating rate kcircRBC (1/h)") # Table 1 human 'kcircRBC' = 0.00037

    # Baseline (homeostatic) cell counts. Table 1 prints a normal range for
    # four of them; the midpoint is used (see vignette).
    lbl_anc <- fixed(log(4.5)); label("Neutrophil baseline Neut0 (10^9 cells/L)") # Table 1 human 'Neut0' = 2-7; midpoint 4.5 (Figure 4 median starts at about 4.5)
    lbl_mono <- fixed(log(0.6)); label("Monocyte baseline Mono0 (10^9 cells/L)") # Table 1 human 'Mono0' printed 0.2-10, read as 0.2-1.0; midpoint 0.6
    lbl_plt <- fixed(log(345)); label("Platelet baseline Plt0 (10^9 cells/L)") # Table 1 human 'Plt0' = 345
    lbl_ret <- fixed(log(77.5)); label("Reticulocyte baseline Ret0 (10^9 cells/L)") # Table 1 human 'Ret0' = 40-115; midpoint 77.5
    lbl_rbc <- fixed(log(5000)); label("Red blood cell baseline RBC0 (10^9 cells/L)") # Table 1 human 'RBC0' = 4,100-5,900; midpoint 5000
    lbl_mpp <- fixed(log(1.3)); label("Multipotent progenitor baseline MPP0 (10^9 cells/L)") # Table 1 human 'MPP0' = 1.3, allometric scaling (Equation S8)
    lbl_cmp <- fixed(log(20.9)); label("Common myeloid progenitor baseline CMP0 (10^9 cells/L)") # Table 1 human 'CMP0' = 20.9, proportion to MPP0 (Equation S9)
    lbl_mep <- fixed(log(15)); label("Megakaryocyte-erythrocyte progenitor baseline MEP0 (10^9 cells/L)") # Table 1 human 'MEP0' = 15, proportion to MPP0 (Equation S9)

    # Mean transit (maturation) times of the three-compartment chains
    lmtt_anc <- fixed(log(210)); label("Neutrophil mean transit time MTTNeut (h)") # Table 1 human 'MTTNeut' = 210
    lmtt_mono <- fixed(log(121.5)); label("Monocyte mean transit time MTTMono (h)") # Table 1 human 'MTTMono' = 121.5
    lmtt_plt <- fixed(log(168)); label("Platelet mean transit time MTTPlt (h)") # Table 1 human 'MTTPlt' = 168
    lmtt_ret <- fixed(log(66)); label("Reticulocyte mean transit time MTTRet (h)") # Table 1 human 'MTTRet' = 66

    # Steady-state reparameterisation ratios (Equation S7), as in the rat
    lambda_t2t1_ret <- fixed(2); label("Ratio T2Ret0/T1Ret0 of the reticulocyte transit compartments, lambda1 (dimensionless)") # Table 1 human 'lambda1' = 2, as in the rat
    lambda_t2t1_plt <- fixed(2); label("Ratio T2Plt0/T1Plt0 of the platelet transit compartments, lambda2 (dimensionless)") # Table 1 human 'lambda2' = 2, as in the rat
    lambda_cmp <- fixed(1.8); label("Ratio of CMP differentiation to proliferation rate, lambda3 (dimensionless)") # Table 1 human 'lambda3' = 1.8, as in the rat
    lambda_mep <- fixed(1.8); label("Ratio of MEP differentiation to proliferation rate, lambda4 (dimensionless)") # Table 1 human 'lambda4' = 1.8, as in the rat
    lambda_mpp <- fixed(1.8); label("Ratio of MPP differentiation to proliferation rate, lambda5 (dimensionless)") # Table 1 human 'lambda5' = 1.8, as in the rat

    # Feedback powers (Equations 11-13), as in the rat
    gamma_stem <- fixed(0.07); label("Feedback power on the stem-cell inflow into MPP (dimensionless)") # Table 1 human 'gamma stem' = 0.07, as in the rat
    gamma_mat_cmp <- fixed(0.6); label("Feedback power on MPP differentiation into CMP (dimensionless)") # Table 1 human 'gamma MatCMP' = 0.6, as in the rat
    gamma_mat_mep <- fixed(0.3); label("Feedback power on MPP differentiation into MEP (dimensionless)") # Table 1 human 'gamma MatMEP' = 0.3, as in the rat
    gamma_prol_trans <- fixed(0.7); label("Feedback power on reticulocyte and platelet transit-cell proliferation (dimensionless)") # Table 1 human 'gamma prolTrans' = 0.7, as in the rat

    # Carboplatin cell-kill parameters: rat slopes scaled by IC50 ratios
    # (Equation 10, Table 2 'IC50 scaled value')
    slope_mpp <- fixed(0.79); label("Carboplatin damage-driven kill slope on MPP (1/uM)") # Table 2 human 'SlopeMPP' IC50-scaled = 0.79
    slope_cmp <- fixed(0.57); label("Carboplatin damage-driven kill slope on CMP (1/uM)") # Table 2 human 'SlopeCMP' IC50-scaled = 0.57
    slope_mep <- fixed(0.66); label("Carboplatin damage-driven kill slope on MEP (1/uM)") # Table 2 human 'SlopeMEP' IC50-scaled = 0.66
    delta_ret <- fixed(2.8); label("Multiplier of SlopeMEP for the proliferative reticulocyte transit cells (dimensionless)") # Table 2 human 'delta Ret' = 2.8, as in rat
    delta_plt <- fixed(0.54); label("Multiplier of SlopeMEP for the proliferative platelet transit cells (dimensionless)") # Table 2 human 'delta Plt' = 0.54, as in rat

    # Log-additive residual error: no human values are reported; the rat
    # Table S4 values are carried over (see vignette).
    expSd_MPP <- fixed(0.33); label("Log-additive residual SD, MPP") # Table S4 rat 'sigma MPP' = 0.33, carried to human
    expSd_CMP <- fixed(0.19); label("Log-additive residual SD, CMP") # Table S4 rat 'sigma CMP' = 0.19, carried to human
    expSd_MEP <- fixed(0.33); label("Log-additive residual SD, MEP") # Table S4 rat 'sigma MEP' = 0.33, carried to human
    expSd_ANC <- fixed(0.17); label("Log-additive residual SD, neutrophils") # Table S4 rat 'sigma Neut' = 0.17, carried to human
    expSd_MONO <- fixed(0.17); label("Log-additive residual SD, monocytes") # Table S4 rat 'sigma Mono' = 0.17, carried to human
    expSd_PLT <- fixed(0.17); label("Log-additive residual SD, platelets") # Table S4 rat 'sigma Ptl' = 0.17, carried to human
    expSd_RET <- fixed(0.36); label("Log-additive residual SD, reticulocytes") # Table S4 rat 'sigma Ret' = 0.36, carried to human
    expSd_RBC <- fixed(0.06); label("Log-additive residual SD, red blood cells") # Table S4 rat 'sigma RBC' = 0.06, carried to human
  })

  model({
    # Individual parameters (typical-value projection, no IIV)
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

    # Free carboplatin in uM: dose in mg over volume in L gives mg/L,
    # converted with the carboplatin molar mass (371.25 g/mol), so that the
    # damage drive is on the same uM scale as the rat Slope estimates.
    Cc <- central / vc * 1000 / 371.25

    # Feedback functions (Equations 11-13). The 1/2 weights are 'alphaWBC =
    # 0.5' in the Supplementary Model Code.
    fdbk_stem <- (0.5 * bl_cmp / cmp + 0.5 * bl_mep / mep)^gamma_stem
    fdbk_mat_cmp <- (0.5 * bl_anc / circ_anc + 0.5 * bl_mono / circ_mono)^gamma_mat_cmp
    fdbk_mat_mep <- (bl_rbc / circ_rbc)^gamma_mat_mep
    fdbk_prol_ret <- (0.5 * bl_ret / circ_ret + 0.5 * bl_rbc / circ_rbc)^gamma_prol_trans
    fdbk_prol_plt <- (bl_plt / circ_plt)^gamma_prol_trans

    # Free-carboplatin PK (Zandvliet 2008, Supplementary Text S11)
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
