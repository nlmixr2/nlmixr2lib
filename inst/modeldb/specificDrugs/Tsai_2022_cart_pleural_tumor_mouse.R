Tsai_2022_cart_pleural_tumor_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Minimal PBPK-PD (mPBPK-PD) model of anti-mesothelin",
    "CAR T cell cellular kinetics and pleural tumor killing after intravenous or",
    "intrapleural delivery (Tsai 2022). Ten ODE states in cell numbers: arterial",
    "and venous blood, lung vasculature and interstitium, a lumped 'other",
    "tissues' vasculature and interstitium, lymph node, and a pleural space",
    "inside the lungs that holds the tumor. CAR T cells cross from vasculature",
    "to interstitium by first-order transmigration, return to blood through",
    "lymph (with a retention factor in the pleural space) and are eliminated in",
    "the lung interstitium. In the pleural space free CAR T cells bind tumor",
    "cells through a cell-level macroscopic association / dissociation pair;",
    "each CAR T - tumor cell complex kills its tumor cell (releasing the CAR T",
    "cell) and drives CAR T proliferation inhibited by total tumor burden. The",
    "tumor grows exponentially and its volume adds to the pleural distribution",
    "volume. Deterministic typical-value model (naive-pooled fit, no IIV).",
    sep = " "
  )

  reference <- paste(
    "Tsai CH, Singh AP, Xia CQ, Wang H. Development of minimal",
    "physiologically-based pharmacokinetic-pharmacodynamic models for",
    "characterizing cellular kinetics of CAR T cells following local deliveries",
    "in mice. J Pharmacokinet Pharmacodyn. 2022;49(5):525-538.",
    "doi:10.1007/s10928-022-09818-8",
    sep = " "
  )

  vignette <- "Tsai_2022_cart_local_delivery"

  units <- list(
    time = "h",
    dosing = "cells",
    concentration = "cells/mL"
  )

  # The pleural space is a paper-specific anatomical space (Tsai 2022 Fig. 2a
  # and Supplementary 'Model Equations', Compartment 6) with no registered
  # canonical; every other state uses a registered PBPK name.
  paper_specific_compartments <- c("pleural_space")

  compartmentData <- list(
    a_arterial = list(analyte = "CAR T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    a_venous = list(analyte = "CAR T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    vp_lung = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_other = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_other = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    lnode = list(analyte = "CAR T cells", units = "cells", specimen = "lymph", verified = TRUE),
    pleural_space = list(analyte = "free CAR T cells", units = "cells", specimen = "tumor", verified = TRUE),
    tumor = list(
      analyte = "mesothelin-positive tumor cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex = list(
      analyte = "CAR T cell - tumor cell complexes",
      units = "complexes",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (NSG, tumor-bearing; human T cells)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    weight_range = "28 g reference mouse (organ volumes and flows from Shah and Betts 2012)",
    disease_state = paste(
      "Orthotopic pleural mesothelioma (mesothelin-positive) xenograft;",
      "the PK layer was checked against radiolabeled exogenous human T cells",
      "in non-tumor-bearing mice."
    ),
    dose_range = paste(
      "0.1, 0.3, 1 and 3 x 10^6 anti-mesothelin CAR T cells as a single",
      "intravenous or intrapleural injection (fitted data); 0.03-3 x 10^6",
      "cells simulated."
    ),
    notes = paste(
      "Naive-pooled fit in Phoenix WinNonlin (proportional error) to",
      "bioluminescence data digitized from Adusumilli 2014 (Sci Transl Med",
      "6:261ra151; CAR T and tumor BLI after i.v. and intrapleural dosing).",
      "Lung, blood, other-tissue and lymph physiology from Shah and Betts 2012",
      "(J Pharmacokinet Pharmacodyn 39:67-86) and T cell transmigration rates",
      "from Khot 2019 (J Pharmacol Exp Ther 368:503), as tabulated in Tsai 2022",
      "Table S2. The number of digitized mice is not reported."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Physiology (Tsai 2022 Table S2; source Shah and Betts 2012 unless noted)
    # ---------------------------------------------------------------------
    q_blood <- fixed(678)
    label("Cardiac blood flow through the lungs (Q_Blood, mL/h)") # Table S2 Q_Blood = 678 mL/hr
    v_arterial <- fixed(0.8585)
    label("Arterial blood volume (V_ab, mL)") # Table S2 V_ab = 0.8585 mL (50 percent of total blood volume)
    v_venous <- fixed(0.8585)
    label("Venous blood volume (V_vb, mL)") # Table S2 V_vb = 0.8585 mL (50 percent of total blood volume)
    l_lung <- fixed(0.746)
    label("Lung lymph flow (L_Lungs, mL/h)") # Table S2 L_Lungs = 0.746 mL/hr
    vv_lung <- fixed(0.0536)
    label("Lung vascular volume (Vv_Lungs, mL)") # Table S2 Vv_Lungs = 0.0536 mL
    vi_lung <- fixed(0.0384)
    label("Lung interstitial volume (Vi_Lungs, mL)") # Table S2 Vi_Lungs = 0.0384 mL
    j_lung <- fixed(1843)
    label("CAR T transmigration rate constant into lung interstitium (J_Lungs, 1/h)") # Table S2 J_Lungs = 1843 1/hr (Khot 2019)
    lkel <- fixed(log(0.84))
    label("T cell elimination rate constant in the lung interstitium (k_el, 1/h)") # Table S2 k_eli = 0.84 1/hr (Zhu 1996)
    q_other <- fixed(678)
    label("Blood flow into other tissues (Q_ot, mL/h)") # Table S2 Q_ot = 678 mL/hr
    l_other <- fixed(0.904)
    label("Lymph flow of other tissues (L_ot, mL/h)") # Table S2 L_ot = 0.904 mL/hr
    vv_other <- fixed(1.432)
    label("Vascular volume of other tissues (Vv_ot, mL)") # Table S2 Vv_ot = 1.432 mL
    vi_other <- fixed(4.870)
    label("Interstitial volume of other tissues (Vi_ot, mL)") # Table S2 Vi_ot = 4.870 mL
    j_other <- fixed(41.4)
    label("CAR T transmigration rate constant into other tissues (J_ot, 1/h)") # Table S2 J_ot = 41.4 1/hr, derived by Eq. 1 from Khot 2019
    l_lymph <- fixed(1.65)
    label("Total lymph flow through the lymph nodes (L_total_Lymph, mL/h)") # Table S2 L_total_Lymph = 1.65 mL/hr
    lv_lnode <- fixed(log(0.113))
    label("Total lymph node volume (Vtotal_LN, mL)") # Table S2 Vtotal_LN = 0.113 mL
    l_pleural <- fixed(0.033)
    label("Pleural space lymph flow (L_PS, mL/h)") # Table S2 L_PS = 0.033 mL/hr (2 percent of lung lymph, as in dogs; Miserocchi 1997)
    vi_pleural <- fixed(0.00728)
    label("Pleural cavity volume (Vi_PS, mL)") # Table S2 Vi_PS = 0.00728 mL (human 0.26 mL/kg; Noppen 2000)

    # ---------------------------------------------------------------------
    # Estimated parameters (Tsai 2022 Table 1 = Table S2 'Estimated' rows)
    # ---------------------------------------------------------------------
    lj_pleural <- log(0.119)
    label("CAR T transmigration rate constant from lung interstitium into the pleural space (J_PS, 1/h)") # Table 1 J_PS = 0.119 1/h (RSE 25.6)
    lr_pleural <- log(2.88)
    label("Retention factor of CAR T cells in the pleural space (R_PS, unitless)") # Table 1 R_PS = 2.88 (RSE 2.60)
    lkg <- log(0.00385)
    label("Exponential tumor growth rate constant (k_g, 1/h)") # Table 1 k_g = 0.00385 1/h (RSE 8.19); doubling time 7.5 days
    lkpro <- log(0.115)
    label("Maximum CAR T proliferation rate constant (k_pro, 1/h)") # Table 1 k_pro = 0.115 1/h (RSE 14.9)
    lki <- log(3.84e7)
    label("Tumor cell number giving 50 percent inhibition of k_pro (KI, cells)") # Table 1 KI = 3.84 x 10^7 cells (RSE 15.4)
    lkkill <- log(0.0733)
    label("Tumor cell killing rate constant per CAR T - tumor cell complex (k_kill, 1/h)") # Table 1 k_kill = 0.0733 1/h (RSE 3.60)
    lkoff <- log(6.85e-8)
    label("Cell-level macroscopic dissociation rate constant (k_off,mac, 1/h)") # Table 1 k_off,mac = 6.85 x 10^-8 1/h (RSE 3.74)

    # ---------------------------------------------------------------------
    # Binding and tumor constants held at assumed / literature values
    # ---------------------------------------------------------------------
    kon <- fixed(1e6)
    label("Cell-level macroscopic association rate constant (k_on,mac, 1/M/s)") # Table S2 k_on,mac = 10^6 M^-1 s^-1 (Faro 2017); Methods: held for identifiability
    r_taa <- fixed(1e4)
    label("Tumor-associated antigen copies per tumor cell (R_TAA, copies/cell)") # Table S2 R_TAA = 10^4 copy/cell, assumed
    v_cell_tumor <- fixed(518.3e-12)
    label("Volume of one tumor cell (mL/cell)") # Supplementary Model Equations: 518.3 fL per cell (Phillips 2012)
    f_isf_tumor <- fixed(0.2)
    label("Interstitial fraction of total tumor volume (unitless)") # Supplementary Model Equations: V_tumor = 518.3e-12 * TB / (1 - 0.2)
    tb0 <- fixed(1e8)
    label("Initial number of tumor cells in the pleural space (TB0, cells)") # Fig. 4, Table 2 and Fig. 6 simulations; per-arm fitted TB0 are in Table S4

    # ---------------------------------------------------------------------
    # Bioluminescence readout (Tsai 2022 Eqs. 3-4 and Table S4)
    # ---------------------------------------------------------------------
    s_bli_cart <- fixed(0.311398)
    label("Bioluminescence per CAR T cell (S1, photons/s/cell)") # Table S4 S1 = 311398 photons/s at the 1 x 10^6-cell first time point, i.e. per 10^6 cells; see vignette
    s_bli_tumor <- fixed(1)
    label("Bioluminescence per tumor cell (S2, photons/s/cell)") # Table S4 S2 = 1 photons/s-cell
    b_bli_tumor <- fixed(754800)
    label("Baseline tumor bioluminescence noise (B2, photons/s)") # Table S4 B2 = 754800 photons/s
  })

  model({
    # k_on,mac (1/M/s) to the per-copy rate used with copies/mL:
    # 3600 s/h and 1000 mL/L over Avogadro's number gives mL/(copy * h).
    avogadro <- 6.02214076e23
    kel <- exp(lkel)
    v_lnode <- exp(lv_lnode)
    kon_cell <- kon * 3600 * 1000 / avogadro

    j_pleural <- exp(lj_pleural)
    r_pleural <- exp(lr_pleural)
    kg <- exp(lkg)
    kpro <- exp(lkpro)
    ki <- exp(lki)
    kkill <- exp(lkkill)
    koff <- exp(lkoff)

    tumor(0) <- tb0

    # Tumor volume adds to the pleural distribution volume (Supplementary
    # Model Equations, Compartment 6).
    v_tumor <- v_cell_tumor * tumor / (1 - f_isf_tumor)
    v_pleural <- vi_pleural + v_tumor

    c_arterial <- a_arterial / v_arterial
    c_venous <- a_venous / v_venous
    c_lung_vas <- vp_lung / vv_lung
    c_lung_int <- is_lung / vi_lung
    c_other_vas <- vp_other / vv_other
    c_other_int <- is_other / vi_other
    c_lnode <- lnode / v_lnode
    c_pleural <- pleural_space / v_pleural

    # CAR T - tumor cell binding on free antigen copies in the pleural space
    binding <- kon_cell * pleural_space * (tumor - complex) * r_taa / v_pleural
    proliferation <- kpro / (1 + tumor / ki) * complex

    d/dt(a_arterial) <- c_lung_vas * q_blood - c_arterial * q_other
    d/dt(a_venous) <- c_other_vas * (q_other - l_other) + c_lnode * l_lymph -
      c_venous * (q_blood + l_lung)
    d/dt(vp_lung) <- c_venous * (q_blood + l_lung) - j_lung * vp_lung - c_lung_vas * q_blood
    d/dt(is_lung) <- j_lung * vp_lung - c_lung_int * (l_lung - l_pleural) - kel * is_lung -
      j_pleural * is_lung
    d/dt(vp_other) <- c_arterial * q_other - j_other * vp_other - c_other_vas * (q_other - l_other)
    d/dt(is_other) <- j_other * vp_other - c_other_int * l_other
    d/dt(lnode) <- c_lung_int * (l_lung - l_pleural) + l_pleural / r_pleural * c_pleural +
      c_other_int * l_other - c_lnode * l_lymph
    d/dt(pleural_space) <- j_pleural * is_lung - l_pleural / r_pleural * c_pleural - binding +
      koff * complex + kkill * complex + proliferation
    d/dt(tumor) <- kg * (tumor - complex) - kkill * complex
    d/dt(complex) <- binding - koff * complex - kkill * complex

    # Outputs. Concentrations are cells/mL (equal to cells/g at unit density).
    # Blood is the venous pool, which reproduces the Table 2 blood exposure.
    Cblood <- c_venous
    Clung <- (vp_lung + is_lung) / (vv_lung + vi_lung)
    Cpleural <- (pleural_space + complex) / v_pleural
    cart_total <- a_arterial + a_venous + vp_lung + is_lung + vp_other + is_other +
      lnode + pleural_space + complex
    bli_cart <- s_bli_cart * cart_total
    bli_tumor <- s_bli_tumor * tumor + b_bli_tumor
  })
}
