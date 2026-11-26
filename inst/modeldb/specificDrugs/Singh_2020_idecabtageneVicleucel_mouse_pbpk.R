# Whole-body PBPK-PD model of anti-BCMA CAR T cells (bb2121, idecabtagene
# vicleucel) in RPMI-8226 xenograft-bearing NSG mice, published by Singh et
# al. 2020 (MAbs 12:1688616), Figure 1C / Figure 5A.
#
# SOURCES. Structure: the authors' Mlxtran model code in the Supplementary
# Material ("A MLXTRAN model code file on the developed PBPK-PD model, as
# described in Fig 1C"), which is the as-run anti-BCMA model and outranks the
# printed main-text Eqs. 13-33 where they differ (see the vignette Errata).
# Physiology: Supplementary Table 1 (25 g mouse). Transmigration rates and the
# liver elimination rate: main-text Table 2 (PBPK block) and the code.
# Anti-BCMA PD parameters: Table 2, PBPK-PD block, BCMA column (identical to
# the code's <PARAMETER> block). Binding constants and densities: the code
# (kon, koff, CAR density) and Table 1 (BCMA density 12,590/cell).
#
# STATE UNITS follow the code: organ states are CAR T-cell CONCENTRATIONS
# (cells/mL); is_tumor is the AMOUNT of CAR T cells in the tumour
# extravascular space (cells); complex is CAR-BCMA complexes per mL of tumour
# extravascular space; tumor is the total tumour volume in mL; effect1-4 carry
# the delayed killing rate (1/h). An intravenous dose (cells) enters blood
# through f(blood) = 1 / v_blood, which is exactly the code's dose of
# Dose / V_Blood cells/mL. An intratumoural dose (cells) goes to is_tumor.

Singh_2020_idecabtageneVicleucel_mouse_pbpk <- function() {
  description <- "Preclinical (mouse, RPMI-8226 xenograft). PBPK-PD (whole-body cell PBPK with CAR-target engagement, CAR T-cell expansion and tumour growth inhibition). Anti-BCMA CAR T cells (bb2121, idecabtagene vicleucel) in female NSG mice bearing subcutaneous BCMA-expressing RPMI-8226 multiple-myeloma tumours. Blood, lung, spleen, kidney, brain, GI tract, liver, 'other' and tumour are each split into vascular and extravascular spaces joined by first-order transmigration and returned through organ lymph flows to a lymph-node pool; CAR T cells are eliminated first-order from the liver extravascular space. In the tumour extravascular space CARs bind BCMA (second-order kon, first-order koff); the number of CAR-target complexes per tumour cell drives Emax expansion of the total intratumoural CAR T-cell pool and, through a four-compartment signal transit chain, Emax killing of an exponentially growing tumour. Deterministic typical-value model; outputs are blood CAR T cells (cells/uL), tumour volume (mm3) and CAR-target complexes per tumour cell."
  reference <- paste(
    "Singh AP, Zheng X, Lin-Schmidt X, Chen W, Carpenter TJ, Zong A, Wang W,",
    "Heald DL (2020). Development of a quantitative relationship between",
    "CAR-affinity, antigen abundance, tumor cell depletion and CAR-T cell",
    "expansion using a multiscale systems PK-PD model. MAbs 12(1):1688616.",
    "doi:10.1080/19420862.2019.1688616.",
    "Model code and physiology in the Supplementary Material (PMC6927769).",
    "Anti-BCMA xenograft data from Friedman et al. 2018 (Hum Gene Ther",
    "29:585-601).",
    sep = " "
  )
  vignette <- "Singh_2020_cart_multiscale"

  # The GI-tract vascular space has no registered vp_ organ token (the
  # register resolves the small and large intestine separately, which this
  # model lumps); its extravascular space uses the registered is_gut form.
  # Same choice as Tsai_2022_cart_liver_tumor_mouse.R.
  paper_specific_compartments <- c("vp_gut")

  units <- list(
    time = "h",
    dosing = "cells (CAR+ T cells; intravenous dose on blood, intratumoural dose on is_tumor)",
    concentration = "cells/uL (Cblood); mm3 (tumor_vol); complexes per tumour cell (cplx_per_tumor)"
  )

  compartmentData <- list(
    blood = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "whole blood", verified = TRUE),
    vp_lung = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_spleen = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_kidney = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_kidney = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_brain = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_brain = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_other = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_other = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_gut = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_gut = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tissue", verified = TRUE),
    lnode = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "lymph", verified = TRUE),
    vp_tumor = list(analyte = "anti-BCMA CAR T cells", units = "cells/mL", specimen = "tumor", verified = TRUE),
    is_tumor = list(
      analyte = "anti-BCMA CAR T cells (total, unbound + tumour-bound)",
      units = "cells",
      specimen = "tumor",
      verified = TRUE
    ),
    complex = list(analyte = "CAR-BCMA complexes", units = "complexes/mL", specimen = "tumor", verified = TRUE),
    tumor = list(analyte = "total tumour volume", units = "mL", specimen = "not applicable", verified = TRUE),
    effect1 = list(
      analyte = "delayed tumour-cell killing rate",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect2 = list(
      analyte = "delayed tumour-cell killing rate",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect3 = list(
      analyte = "delayed tumour-cell killing rate",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect4 = list(
      analyte = "delayed tumour-cell killing rate",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female NSG, subcutaneous RPMI-8226 xenograft; human CAR T cells)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "25 g reference mouse (organ volumes and flows, Supplementary Table 1)",
    disease_state = "BCMA-expressing RPMI-8226 multiple-myeloma xenograft (12,590 BCMA per cell)",
    dose_range = "Single intravenous dose of 5 x 10^6 anti-BCMA (bb2121) CAR+ T cells per mouse (fitted); 0.1-50 x 10^6 simulated",
    notes = paste(
      "Tumour-volume and blood CAR T-cell profiles from Friedman et al. 2018",
      "(Table 1 row 'Anti-BCMA CAR-T'). The PBPK",
      "transmigration and liver-elimination rates were estimated beforehand by",
      "simultaneous fitting of untransduced T-cell, anti-EGFR and anti-CD19 CAR",
      "T-cell biodistribution data (Table 1, In Vivo Biodistribution Studies)",
      "and were fixed in the PBPK-PD fit. Fitting used Monolix 8 (SAEM); no",
      "between-animal variability was retained (all omegas 0 in the deposited",
      "project) and residual-error estimates were not reported."
    )
  )

  ini({
    # ---- Physiology, 25 g mouse (Supplementary Table 1; code 'Physiological parameters') ----
    q_lung <- fixed(373)
    label("Lung blood flow (cardiac output) (Q_Lung, mL/h)") # Supp Table 1 Lung Q = 373 mL/h
    q_spleen <- fixed(8.18)
    label("Spleen blood flow (Q_Spleen, mL/h)") # Supp Table 1 Spleen Q = 8.18 mL/h
    q_gi <- fixed(75.4)
    label("GI-tract blood flow (Q_GI, mL/h)") # Supp Table 1 GI Q = 75.4 mL/h
    q_liver <- fixed(10.3)
    label("Hepatic-artery blood flow (Q_Liver, mL/h)") # Supp Table 1 Liver Q = 10.3 mL/h
    q_kidney <- fixed(68.5)
    label("Kidney blood flow (Q_Kidney, mL/h)") # Supp Table 1 Kidney Q = 68.5 mL/h
    q_brain <- fixed(11.8)
    label("Brain blood flow (Q_Brain, mL/h)") # Supp Table 1 Brain Q = 11.8 mL/h
    q_tumor <- fixed(6.0)
    label("Tumour blood flow (Q_Tumor, mL/h)") # Supp Table 1 Tumor Q = 6.0 mL/h
    q_other <- fixed(192.82)
    label("Blood flow to other tissues (Q_Other, mL/h)") # code Q_Other = 192.82 (373 minus all other organs); Supp Table 1 prints 192.8
    l_lung <- fixed(0.0746)
    label("Lung lymph flow (L_Lung, mL/h)") # Supp Table 1 Lung L = 0.0746 mL/h
    l_spleen <- fixed(0.00164)
    label("Spleen lymph flow (L_Spleen, mL/h)") # Supp Table 1 Spleen L = 0.00164 mL/h
    l_gi <- fixed(0.0151)
    label("GI-tract lymph flow (L_GI, mL/h)") # Supp Table 1 GI L = 0.0151 mL/h
    l_liver <- fixed(0.00206)
    label("Liver lymph flow (L_Liver, mL/h)") # Supp Table 1 Liver L = 0.00206 mL/h
    l_kidney <- fixed(0.0137)
    label("Kidney lymph flow (L_Kidney, mL/h)") # Supp Table 1 Kidney L = 0.0137 mL/h
    l_brain <- fixed(0.00236)
    label("Brain lymph flow (L_Brain, mL/h)") # Supp Table 1 Brain L = 0.00236 mL/h
    l_tumor <- fixed(0.0012)
    label("Tumour lymph flow (L_Tumor, mL/h)") # Supp Table 1 Tumor L = 0.0012 mL/h
    l_other <- fixed(0.0386)
    label("Lymph flow from other tissues (L_Other, mL/h)") # Supp Table 1 Other L = 0.0386 mL/h
    l_lymph <- fixed(0.14926)
    label("Total lymph flow from the lymph node to blood (L_LN, mL/h)") # code L_LN = 0.14926 (sum of organ lymph flows); Supp Table 1 LN L = 0.149
    v_blood <- fixed(0.944)
    label("Whole-blood volume (V_Blood, mL)") # Supp Table 1 Whole Blood = 0.944 mL
    vv_lung <- fixed(0.0295)
    label("Lung vascular volume (mL)") # Supp Table 1 Lung V^V = 0.0295 mL
    vi_lung <- fixed(0.0384)
    label("Lung extravascular volume (mL)") # Supp Table 1 Lung V^EV = 0.0384 mL
    vv_spleen <- fixed(0.0154)
    label("Spleen vascular volume (mL)") # Supp Table 1 Spleen V^V = 0.0154 mL
    vi_spleen <- fixed(0.0254)
    label("Spleen extravascular volume (mL)") # Supp Table 1 Spleen V^EV = 0.0254 mL
    vv_gi <- fixed(0.0166)
    label("GI-tract vascular volume (mL)") # Supp Table 1 GI V^V = 0.0166 mL
    vi_gi <- fixed(0.1815)
    label("GI-tract extravascular volume (mL)") # Supp Table 1 GI V^EV = 0.1815 mL
    vv_liver <- fixed(0.164)
    label("Liver vascular volume (mL)") # Supp Table 1 Liver V^V = 0.164 mL
    vi_liver <- fixed(0.385)
    label("Liver extravascular volume (mL)") # Supp Table 1 Liver V^EV = 0.385 mL
    vv_kidney <- fixed(0.0289)
    label("Kidney vascular volume (mL)") # Supp Table 1 Kidney V^V = 0.0289 mL
    vi_kidney <- fixed(0.0788)
    label("Kidney extravascular volume (mL)") # Supp Table 1 Kidney V^EV = 0.0788 mL
    vv_brain <- fixed(0.0107)
    label("Brain vascular volume (mL)") # Supp Table 1 Brain V^V = 0.0107 mL
    vi_brain <- fixed(0.0873)
    label("Brain extravascular volume (mL)") # Supp Table 1 Brain V^EV = 0.0873 mL
    vv_other <- fixed(0.55209)
    label("Vascular volume of other tissues (mL)") # code V_Other_v = 0.55209; Supp Table 1 prints 0.552
    vi_other <- fixed(4.11183)
    label("Extravascular volume of other tissues (mL)") # code V_Other_e = 4.11183; Supp Table 1 prints 4.11
    lv_lnode <- fixed(log(0.113))
    label("Lymph-node volume (V_LN, mL)") # Supp Table 1 LN total volume = 0.113 mL; code V_LN = 0.113
    f_vasc_tumor <- fixed(0.113)
    label("Vascular fraction of total tumour volume (unitless)") # Methods: VTumorV = 11.3% of total tumour volume (Eq. 33 block); code V_Tumor_v = 0.113 * V_Tumor

    # ---- Transmigration and elimination (Table 2, PBPK block; fixed in the PBPK-PD fit) ----
    j_lung <- fixed(4705.13)
    label("Vascular-to-extravascular transmigration rate, lung (J_Lung, 1/h)") # Table 2 JCAR-T Lung = 4705.13 (24%), estimated in the PBPK step
    j_spleen <- fixed(1799.76)
    label("Vascular-to-extravascular transmigration rate, spleen (J_Spleen, 1/h)") # Table 2 JCAR-T Spleen = 1799.76 (47.1%)
    j_liver <- fixed(716.49)
    label("Vascular-to-extravascular transmigration rate, liver (J_Liver, 1/h)") # Table 2 JCAR-T Liver = 716.49 (46.1%)
    j_kidney <- fixed(296.75)
    label("Vascular-to-extravascular transmigration rate, kidney (J_Kidney, 1/h)") # Table 2 JCAR-T Kidney = 296.75 (31.1%)
    j_tumor <- fixed(10.35)
    label("Vascular-to-extravascular transmigration rate, tumour (J_Tumor, 1/h)") # Table 2 JCAR-T Tumor = 10.35 (32.2%)
    j_other <- fixed(86.8)
    label("Vascular-to-extravascular transmigration rate, other tissues (J_Others, 1/h)") # Table 2 JCAR-T Others = 86.8 (Fixed), from Khot 2019
    j_brain <- fixed(1.4)
    label("Vascular-to-extravascular transmigration rate, brain (J_Brain, 1/h)") # Table 2 JCAR-T Brain = 1.4 (Fixed), from Khot 2019
    j_gi <- fixed(9.05)
    label("Vascular-to-extravascular transmigration rate, GI tract (J_GI, 1/h)") # Table 2 JCAR-T GI = 9.05 (Fixed), from Khot 2019
    lkel <- fixed(log(0.0288))
    label("First-order CAR T-cell elimination rate from the liver extravascular space (Kel_Liver, 1/h)") # code Kdep_Liver = 0.0288; Table 2 KelLiver = 0.029 (15.1%) is the rounded value

    # ---- CAR-BCMA binding and system constants (code; Table 1; Methods) ----
    kon <- fixed(1.08e-12)
    label("CAR-BCMA second-order association rate constant (1/((cells/mL) h))") # code adm2 amount 1.08e-12 into kon, '1/(#/ml)/hr'; with koff gives Kd = 10 nM (Table 1 anti-BCMA Kd = 10 nM)
    koff <- fixed(6.48)
    label("CAR-BCMA first-order dissociation rate constant (1/h)") # code adm3 amount 6.48 into koff, '1/hr'
    ag_tumor <- fixed(12590)
    label("BCMA density on RPMI-8226 tumour cells (receptors/cell)") # Table 1 'RPMI-8226 MM cells (12,590/cell)'; code adm4 amount 12590 into Density_TAA
    ag_car <- fixed(5000)
    label("CAR density on CAR T cells (receptors/cell)") # code Density_CAR = 5000
    tumor_cell_density <- fixed(1e8)
    label("Tumour-cell packing density in tumour tissue (cells/mL)") # Methods after Eq. 20: CTumorCell set to 10^8 tumour cells/mL; code TargetCell = 10^8
    tumor0 <- fixed(0.15)
    label("Initial total tumour volume at CAR T-cell dosing (mL)") # Methods (model-based simulations): 150 mm3 initial tumour burden; code adm5 doses 0.05-0.5 mL into V_Tumor. Scenario input, override per study

    # ---- Anti-BCMA PBPK-PD parameters (Table 2, PBPK-PD block, BCMA column; code <PARAMETER>) ----
    lkg <- log(0.00305)
    label("First-order exponential tumour growth rate (KgExp, 1/h)") # Table 2 KgExpTumor BCMA = 0.00305 (12.7%); code Kg_Exp_Tumor_pop = 0.00305
    lkkill_max <- log(0.0486)
    label("Maximum first-order tumour-cell killing rate (Kkillmax, 1/h)") # Table 2 KkillmaxCAR-T BCMA = 0.0486 (5.63%); code Kkill_max_pop = 0.0486
    lkc50_kill <- log(9.94)
    label("CAR-target complexes per tumour cell giving half-maximal killing (KC50, complexes/cell)") # Table 2 KC50CAR-T BCMA = 9.94 (5.93%); code IC50_pop = 9.94
    gam_kill <- fixed(1)
    label("Hill coefficient of the killing function (unitless)") # Table 2 gamma TGI CAR-T BCMA = 1 (fixed); code gamma_kill_pop = 1
    ltau_kill <- log(21.6)
    label("Transit time of each killing-signal transit compartment (tau, h)") # Table 2 tau BCMA = 21.6 (5.24%); code tau_kill_pop = 21.6
    lkexp_max <- log(0.0906)
    label("Maximum first-order in vivo CAR T-cell expansion rate (TMaxAct, 1/h)") # Table 2 TMaxAct BCMA = 0.0906 (2.98%); code Kexp_CART_max_pop = 0.0906
    lec50_exp <- log(5.18)
    label("CAR-target complexes per tumour cell giving half-maximal expansion (EC50Act, complexes/cell)") # Table 2 EC50Act BCMA = 5.18 (6.01%); code EC50_pop = 5.18
    gam_exp <- fixed(1)
    label("Hill coefficient of the expansion function (unitless)") # code gamma_act_pop = 1 (omega 0); Eq. 26 prints the hyperbolic form
  })

  model({
    kel <- exp(lkel)
    v_lnode <- exp(lv_lnode)
    kg <- exp(lkg)
    kkill_max <- exp(lkkill_max)
    kc50_kill <- exp(lkc50_kill)
    tau_kill <- exp(ltau_kill)
    kexp_max <- exp(lkexp_max)
    ec50_exp <- exp(lec50_exp)

    # Tumour volume (mL) is split 11.3 / 88.7 percent vascular / extravascular.
    vv_tumor <- f_vasc_tumor * tumor
    vi_tumor <- (1 - f_vasc_tumor) * tumor

    # Complexes per tumour cell: complex is per mL of tumour extravascular
    # space, which holds tumor_cell_density tumour cells per mL.
    cplx_per_tumor <- complex / tumor_cell_density
    kexp <- kexp_max * cplx_per_tumor^gam_exp / (cplx_per_tumor^gam_exp + ec50_exp^gam_exp)
    kkill <- kkill_max * cplx_per_tumor^gam_kill / (cplx_per_tumor^gam_kill + kc50_kill^gam_kill)

    # Unbound CAR T-cell concentration in the tumour extravascular space; only
    # unbound cells leave the tumour by lymph.
    c_tumor_free <- is_tumor / vi_tumor - complex / ag_car

    d/dt(blood) <- (-(q_lung + l_lung) * blood +
      (q_gi - l_gi + q_spleen - l_spleen + q_liver - l_liver) * vp_liver +
      (q_kidney - l_kidney) * vp_kidney + (q_brain - l_brain) * vp_brain +
      l_lymph * lnode + (q_other - l_other) * vp_other +
      (q_tumor - l_tumor) * vp_tumor) / v_blood
    d/dt(vp_lung) <- ((q_lung + l_lung) * blood - j_lung * vp_lung * vv_lung - q_lung * vp_lung) / vv_lung
    d/dt(is_lung) <- (j_lung * vp_lung * vv_lung - l_lung * is_lung) / vi_lung
    d/dt(vp_spleen) <- (q_spleen * vp_lung - j_spleen * vp_spleen * vv_spleen - (q_spleen - l_spleen) * vp_spleen) / vv_spleen
    d/dt(is_spleen) <- (j_spleen * vp_spleen * vv_spleen - l_spleen * is_spleen) / vi_spleen
    d/dt(vp_kidney) <- (q_kidney * vp_lung - j_kidney * vp_kidney * vv_kidney - (q_kidney - l_kidney) * vp_kidney) / vv_kidney
    d/dt(is_kidney) <- (j_kidney * vp_kidney * vv_kidney - l_kidney * is_kidney) / vi_kidney
    d/dt(vp_brain) <- (q_brain * vp_lung - j_brain * vp_brain * vv_brain - (q_brain - l_brain) * vp_brain) / vv_brain
    d/dt(is_brain) <- (j_brain * vp_brain * vv_brain - l_brain * is_brain) / vi_brain
    d/dt(vp_other) <- (q_other * vp_lung - j_other * vp_other * vv_other - (q_other - l_other) * vp_other) / vv_other
    d/dt(is_other) <- (j_other * vp_other * vv_other - l_other * is_other) / vi_other
    d/dt(vp_gut) <- (q_gi * vp_lung - j_gi * vp_gut * vv_gi - (q_gi - l_gi) * vp_gut) / vv_gi
    d/dt(is_gut) <- (j_gi * vp_gut * vv_gi - l_gi * is_gut) / vi_gi
    d/dt(vp_liver) <- (q_liver * vp_lung - j_liver * vp_liver * vv_liver -
      (q_gi - l_gi + q_spleen - l_spleen + q_liver - l_liver) * vp_liver +
      (q_gi - l_gi) * vp_gut + (q_spleen - l_spleen) * vp_spleen) / vv_liver
    d/dt(is_liver) <- (j_liver * vp_liver * vv_liver - l_liver * is_liver) / vi_liver - kel * is_liver
    d/dt(lnode) <- (l_liver * is_liver + l_gi * is_gut + l_tumor * c_tumor_free +
      l_other * is_other + l_brain * is_brain + l_kidney * is_kidney +
      l_spleen * is_spleen + l_lung * is_lung - l_lymph * lnode) / v_lnode
    d/dt(vp_tumor) <- (q_tumor * vp_lung - j_tumor * vp_tumor * vv_tumor - (q_tumor - l_tumor) * vp_tumor) / vv_tumor
    d/dt(is_tumor) <- j_tumor * vp_tumor * vv_tumor - l_tumor * c_tumor_free + kexp * is_tumor
    d/dt(complex) <- kon * ((is_tumor / vi_tumor) * ag_car - complex) *
      (tumor_cell_density * ag_tumor - complex) - koff * complex
    d/dt(tumor) <- kg * tumor - effect4 * tumor
    d/dt(effect1) <- (kkill - effect1) / tau_kill
    d/dt(effect2) <- (effect1 - effect2) / tau_kill
    d/dt(effect3) <- (effect2 - effect3) / tau_kill
    d/dt(effect4) <- (effect3 - effect4) / tau_kill

    tumor(0) <- tumor0
    # An intravenous dose in cells becomes a blood concentration (cells/mL).
    f(blood) <- 1 / v_blood

    # Outputs: blood CAR T cells per uL and tumour volume in mm3.
    Cblood <- blood / 1000
    tumor_vol <- max(tumor, 0) * 1000
  })
}
