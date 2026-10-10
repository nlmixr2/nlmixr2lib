Tsai_2022_cart_liver_tumor_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Minimal PBPK-PD (mPBPK-PD) model of CAR T cell",
    "cellular kinetics and killing of a liver tumor after intravenous, portal",
    "vein or local hepatic artery delivery (Tsai 2022, theoretical simulation",
    "model). Seventeen ODE states in cell numbers: arterial and venous blood,",
    "lymph node, and vasculature plus interstitium of the lungs, liver,",
    "gastrointestinal tract, spleen, lumped 'other tissues' and a liver tumor.",
    "The tumor shares hepatic-artery and portal-vein (GI and splenic) inflow",
    "with the liver, drains to venous blood, and has vascular and interstitial",
    "volumes proportional to its cell number. CAR T cells transmigrate from",
    "vasculature to interstitium, return through lymph (with retention factors",
    "in liver, tumor and spleen), and are eliminated in the lung interstitium.",
    "In the tumor interstitium CAR T cells bind tumor cells, kill them, and",
    "proliferate with tumor-burden inhibition, using the pleural-tumor PD",
    "parameters with a twofold higher killing rate. Deterministic",
    "typical-value model.",
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

  # The gastrointestinal-tract vascular space has no registered vp_ organ
  # entry (vp_small_intestine / vp_large_intestine resolve the two segments,
  # which Tsai 2022 lumps); its interstitium uses the registered is_gut form.
  paper_specific_compartments <- c("vp_gut")

  compartmentData <- list(
    a_arterial = list(analyte = "CAR T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    a_venous = list(analyte = "CAR T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    vp_lung = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_gut = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_gut = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_spleen = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_other = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_other = list(analyte = "CAR T cells", units = "cells", specimen = "tissue", verified = TRUE),
    lnode = list(analyte = "CAR T cells", units = "cells", specimen = "lymph", verified = TRUE),
    vp_tumor = list(analyte = "CAR T cells", units = "cells", specimen = "tumor", verified = TRUE),
    is_tumor = list(analyte = "free CAR T cells", units = "cells", specimen = "tumor", verified = TRUE),
    tumor = list(analyte = "tumor cells", units = "cells", specimen = "not applicable", verified = TRUE),
    complex = list(
      analyte = "CAR T cell - tumor cell complexes",
      units = "complexes",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (tumor-bearing; human T cells)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "28 g reference mouse (organ volumes and flows from Shah and Betts 2012)",
    disease_state = paste(
      "Hypothetical liver tumor (2.5 x 10^8 cells, about 0.2 mL, about 10",
      "percent of liver volume) receiving about 10 percent of total hepatic",
      "blood flow; theoretical simulation, not fitted to tumor data."
    ),
    dose_range = paste(
      "0.3 and 3 x 10^6 CAR T cells as a single intravenous, portal vein or",
      "local hepatic artery dose (simulated)."
    ),
    notes = paste(
      "The non-tumor mPBPK layer (blood, lungs, liver, GI tract, spleen) was",
      "compared with and fitted to radiolabeled exogenous human T cell",
      "biodistribution digitized from Khot 2019 (J Pharmacol Exp Ther",
      "368:503). Physiology from Shah and Betts 2012 (J Pharmacokinet",
      "Pharmacodyn 39:67-86) and transmigration / retention constants from",
      "Khot 2019, as tabulated in Tsai 2022 Table S3. PD parameters were taken",
      "from the anti-mesothelin pleural tumor fit (Tsai_2022_cart_pleural_tumor_mouse)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Blood, lungs, lymph node (Tsai 2022 Table S3; Shah and Betts 2012)
    # ---------------------------------------------------------------------
    q_blood <- fixed(678)
    label("Cardiac blood flow through the lungs (Q_Blood, mL/h)") # Table S3 Q_Blood = 678 mL/hr
    v_arterial <- fixed(0.8585)
    label("Arterial blood volume (V_ab, mL)") # Table S3 V_ab = 0.8585 mL
    v_venous <- fixed(0.8585)
    label("Venous blood volume (V_vb, mL)") # Table S3 V_vb = 0.8585 mL
    l_lung <- fixed(0.746)
    label("Lung lymph flow (L_Lungs, mL/h)") # Table S3 L_Lungs = 0.746 mL/hr
    vv_lung <- fixed(0.0536)
    label("Lung vascular volume (Vv_Lungs, mL)") # Table S3 Vv_Lungs = 0.0536 mL
    vi_lung <- fixed(0.0384)
    label("Lung interstitial volume (Vi_Lungs, mL)") # Table S3 Vi_Lungs = 0.0384 mL
    j_lung <- fixed(1843)
    label("CAR T transmigration rate constant into lung interstitium (J_Lungs, 1/h)") # Table S3 J_Lungs = 1843 1/hr (Khot 2019)
    lkel <- fixed(log(0.84))
    label("T cell elimination rate constant in the lung interstitium (k_el, 1/h)") # Table S3 k_eli = 0.84 1/hr (Zhu 1996)
    l_lymph <- fixed(1.65)
    label("Total lymph flow through the lymph nodes (L_total_Lymph, mL/h)") # Table S3 L_total_Lymph = 1.65 mL/hr
    lv_lnode <- fixed(log(0.113))
    label("Total lymph node volume (Vtotal_LN, mL)") # Table S3 Vtotal_LN = 0.113 mL

    # ---------------------------------------------------------------------
    # Liver, GI tract, spleen (Tsai 2022 Table S3)
    # ---------------------------------------------------------------------
    q_liver_ha <- fixed(18.7)
    label("Hepatic artery blood flow into the liver (Q_Liver,HA, mL/h)") # Table S3 Q_Liver,HA = 18.7 mL/hr
    l_liver <- fixed(0.1874)
    label("Total liver lymph flow (L_total,Liver, mL/h)") # Table S3 L_total,Liver = 0.1874 mL/hr
    vv_liver <- fixed(0.298)
    label("Liver vascular volume (Vv_Liver, mL)") # Table S3 Vv_Liver = 0.298 mL
    vi_liver <- fixed(0.385)
    label("Liver interstitial volume (Vi_Liver, mL)") # Table S3 Vi_Liver = 0.385 mL
    j_liver <- fixed(126.9)
    label("CAR T transmigration rate constant into liver interstitium (J_Liver, 1/h)") # Table S3 J_Liver = 126.9 1/hr (Khot 2019)
    r_liver <- fixed(2.5)
    label("Retention factor of CAR T cells in liver and liver tumor (R_Liver, unitless)") # Table S3 R_Liver = 2.5 (Khot 2019)
    q_gi <- fixed(137)
    label("Blood flow into the GI tract (Q_GI, mL/h)") # Table S3 Q_GI = 137 mL/hr
    l_gi <- fixed(0.1508)
    label("GI tract lymph flow (L_GI, mL/h)") # Table S3 L_GI = 0.1508 mL/hr
    vv_gi <- fixed(0.03019)
    label("GI tract vascular volume (Vv_GI, mL)") # Table S3 Vv_GI = 0.03019 mL
    vi_gi <- fixed(0.1815)
    label("GI tract interstitial volume (Vi_GI, mL)") # Table S3 Vi_GI = 0.1815 mL
    j_gi <- fixed(18.1)
    label("CAR T transmigration rate constant into GI interstitium (J_GI, 1/h)") # Table S3 J_GI = 18.1 1/hr (Khot 2019)
    q_spleen <- fixed(14.88)
    label("Blood flow into the spleen (Q_Spleen, mL/h)") # Table S3 Q_Spleen = 14.88 mL/hr
    l_spleen <- fixed(0.01636)
    label("Spleen lymph flow (L_Spleen, mL/h)") # Table S3 L_Spleen = 0.01636 mL/hr
    vv_spleen <- fixed(0.028)
    label("Spleen vascular volume (Vv_Spleen, mL)") # Table S3 Vv_Spleen = 0.028 mL
    vi_spleen <- fixed(0.0254)
    label("Spleen interstitial volume (Vi_Spleen, mL)") # Table S3 Vi_Spleen = 0.0254 mL
    j_spleen <- fixed(114)
    label("CAR T transmigration rate constant into spleen interstitium (J_Spleen, 1/h)") # Table S3 J_Spleen = 114 1/hr (Khot 2019)
    r_spleen <- fixed(9.8)
    label("Retention factor of CAR T cells in the spleen (R_Spleen, unitless)") # Table S3 R_Spleen = 9.8 (Khot 2019)

    # ---------------------------------------------------------------------
    # Other tissues (Tsai 2022 Table S3)
    # ---------------------------------------------------------------------
    q_other <- fixed(507.42)
    label("Blood flow into other tissues (Q_ot, mL/h)") # Table S3 Q_ot = 507.42 mL/hr
    l_other <- fixed(0.549)
    label("Lymph flow of other tissues (L_ot, mL/h)") # Table S3 L_ot = 0.549 mL/hr
    vv_other <- fixed(1.076)
    label("Vascular volume of other tissues (Vv_ot, mL)") # Table S3 Vv_ot = 1.076 mL
    vi_other <- fixed(4.278)
    label("Interstitial volume of other tissues (Vi_ot, mL)") # Table S3 Vi_ot = 4.278 mL
    j_other <- fixed(16.7)
    label("CAR T transmigration rate constant into other tissues (J_ot, 1/h)") # Table S3 J_ot = 16.7 1/hr derived by Eq. 1 (used for Figs. 7-8); the fitted alternative is 94.3 1/hr (RSE 20.0), used for Figs. S5, S6, S8

    # ---------------------------------------------------------------------
    # Liver tumor (Tsai 2022 Table S3 'Tumor' rows, all assumed)
    # ---------------------------------------------------------------------
    q_tumor_ha <- fixed(2)
    label("Hepatic artery blood flow into the tumor (Q_Tumor,HA, mL/h)") # Table S3 Q_Tumor,HA = 2 mL/hr, assumed
    q_tumor_pv_gi <- fixed(14.5)
    label("Portal vein blood flow of GI origin into the tumor (Q_Tumor,PV,GI, mL/h)") # Table S3 Q_Tumor,PV,GI = 14.5 mL/hr, assumed
    q_tumor_pv_spleen <- fixed(1.5)
    label("Portal vein blood flow of splenic origin into the tumor (Q_Tumor,PV,Spleen, mL/h)") # Table S3 Q_Tumor,PV,Spleen = 1.5 mL/hr, assumed
    l_tumor_ha <- fixed(0.004)
    label("Tumor lymph flow attributed to the hepatic artery supply (L_Tumor,HA, mL/h)") # Table S3 L_Tumor,HA = 0.004 mL/hr, assumed
    l_tumor_pv <- fixed(0.032)
    label("Tumor lymph flow attributed to the portal vein supply (L_Tumor,PV, mL/h)") # Table S3 L_Tumor,PV = 0.032 mL/hr, assumed
    j_tumor <- fixed(126.9)
    label("CAR T transmigration rate constant into tumor interstitium (J_Tumor, 1/h)") # Table S3 J_Tumor = 126.9 1/hr, assumed equal to liver
    f_vasc_tumor <- fixed(0.15)
    label("Vascular fraction of total tumor volume (unitless)") # Table S3 Vv_Tumor = 15 percent of V_total Tumor, liver proportion
    f_isf_tumor <- fixed(0.2)
    label("Interstitial fraction of total tumor volume (unitless)") # Table S3 Vi_Tumor = 20 percent of V_total Tumor, liver proportion
    v_cell_tumor <- fixed(518.3e-12)
    label("Volume of one tumor cell (mL/cell)") # Supplementary Model Equations ii: 518.3 fL per cell (Phillips 2012)
    tb0 <- fixed(2.5e8)
    label("Initial number of liver tumor cells (TB0, cells)") # Table S3 TB0 = 2.5 x 10^8 cells, assumed

    # ---------------------------------------------------------------------
    # PD parameters transferred from the pleural tumor fit (Table S3)
    # ---------------------------------------------------------------------
    kg <- fixed(0.00385)
    label("Exponential tumor growth rate constant (k_g, 1/h)") # Table S3 k_g = 0.00385 1/hr, pleural tumor parameter
    kpro <- fixed(0.115)
    label("Maximum CAR T proliferation rate constant (k_pro, 1/h)") # Table S3 k_pro = 0.115 1/hr, pleural tumor parameter
    ki <- fixed(3.84e7)
    label("Tumor cell number giving 50 percent inhibition of k_pro (KI, cells)") # Table S3 KI = 3.84 x 10^7 cells, pleural tumor parameter
    kkill <- fixed(0.1466)
    label("Tumor cell killing rate constant per CAR T - tumor cell complex (k_kill, 1/h)") # Table S3 'assumed twice as pleural tumor' and Methods (twofold higher): 2 x 0.0733; Table S3 prints the pleural 0.0733
    r_taa <- fixed(1e4)
    label("Tumor-associated antigen copies per tumor cell (R_TAA, copies/cell)") # Table S3 R_TAA = 10^4 copy/cell, pleural tumor parameter
    kon <- fixed(1e6)
    label("Cell-level macroscopic association rate constant (k_on,mac, 1/M/s)") # Table S3 k_on,mac = 10^6 M^-1 s^-1, pleural tumor parameter
    koff <- fixed(6.85e-8)
    label("Cell-level macroscopic dissociation rate constant (k_off,mac, 1/h)") # Table S3 k_off,mac = 6.85 x 10^-8 1/hr, pleural tumor parameter
  })

  model({
    # k_on,mac (1/M/s) to the per-copy rate used with copies/mL:
    # 3600 s/h and 1000 mL/L over Avogadro's number gives mL/(copy * h).
    avogadro <- 6.02214076e23
    kel <- exp(lkel)
    v_lnode <- exp(lv_lnode)
    kon_cell <- kon * 3600 * 1000 / avogadro

    tumor(0) <- tb0

    # Tumor volumes scale with the tumor cell number (Supplementary Model
    # Equations ii: tumor cells occupy 1 - 0.2 - 0.15 of the tumor volume).
    v_tumor <- v_cell_tumor * tumor / (1 - f_isf_tumor - f_vasc_tumor)
    vv_tumor <- f_vasc_tumor * v_tumor
    vi_tumor <- f_isf_tumor * v_tumor

    # Flows routed through the tumor are subtracted from the normal liver.
    q_tumor <- q_tumor_ha + q_tumor_pv_gi + q_tumor_pv_spleen
    l_tumor <- l_tumor_ha + l_tumor_pv
    l_liver_normal <- l_liver - l_tumor_ha - l_tumor_pv
    q_liver_out <- (q_liver_ha - q_tumor_ha) + (q_gi - l_gi - q_tumor_pv_gi) +
      (q_spleen - l_spleen - q_tumor_pv_spleen) - l_liver_normal

    c_arterial <- a_arterial / v_arterial
    c_venous <- a_venous / v_venous
    c_lung_vas <- vp_lung / vv_lung
    c_lung_int <- is_lung / vi_lung
    c_liver_vas <- vp_liver / vv_liver
    c_liver_int <- is_liver / vi_liver
    c_gi_vas <- vp_gut / vv_gi
    c_gi_int <- is_gut / vi_gi
    c_spleen_vas <- vp_spleen / vv_spleen
    c_spleen_int <- is_spleen / vi_spleen
    c_other_vas <- vp_other / vv_other
    c_other_int <- is_other / vi_other
    c_lnode <- lnode / v_lnode
    c_tumor_vas <- vp_tumor / vv_tumor
    c_tumor_int <- is_tumor / vi_tumor

    # CAR T - tumor cell binding on free antigen copies in the tumor interstitium
    binding <- kon_cell * is_tumor * (tumor - complex) * r_taa / vi_tumor
    proliferation <- kpro / (1 + tumor / ki) * complex

    d/dt(a_arterial) <- c_lung_vas * q_blood - c_arterial * (q_liver_ha + q_gi + q_spleen + q_other)
    d/dt(a_venous) <- c_other_vas * (q_other - l_other) + c_lnode * l_lymph +
      c_liver_vas * q_liver_out + c_tumor_vas * (q_tumor - l_tumor) -
      c_venous * (q_blood + l_lung)
    d/dt(vp_lung) <- c_venous * (q_blood + l_lung) - j_lung * vp_lung - c_lung_vas * q_blood
    d/dt(is_lung) <- j_lung * vp_lung - c_lung_int * l_lung - kel * is_lung
    d/dt(vp_liver) <- c_arterial * (q_liver_ha - q_tumor_ha) +
      c_gi_vas * (q_gi - l_gi - q_tumor_pv_gi) +
      c_spleen_vas * (q_spleen - l_spleen - q_tumor_pv_spleen) -
      j_liver * vp_liver - c_liver_vas * q_liver_out
    d/dt(is_liver) <- j_liver * vp_liver - c_liver_int * l_liver_normal / r_liver
    d/dt(vp_gut) <- c_arterial * q_gi - j_gi * vp_gut - c_gi_vas * (q_gi - l_gi)
    d/dt(is_gut) <- j_gi * vp_gut - c_gi_int * l_gi
    d/dt(vp_spleen) <- c_arterial * q_spleen - j_spleen * vp_spleen - c_spleen_vas * (q_spleen - l_spleen)
    d/dt(is_spleen) <- j_spleen * vp_spleen - c_spleen_int * l_spleen / r_spleen
    d/dt(vp_other) <- c_arterial * q_other - j_other * vp_other - c_other_vas * (q_other - l_other)
    d/dt(is_other) <- j_other * vp_other - c_other_int * l_other
    d/dt(lnode) <- c_lung_int * l_lung + c_gi_int * l_gi + c_spleen_int * l_spleen / r_spleen +
      c_liver_int * l_liver_normal / r_liver + l_tumor * c_tumor_int / r_liver +
      c_other_int * l_other - c_lnode * l_lymph
    d/dt(vp_tumor) <- c_arterial * q_tumor_ha + c_gi_vas * q_tumor_pv_gi +
      c_spleen_vas * q_tumor_pv_spleen - j_tumor * vp_tumor -
      c_tumor_vas * (q_tumor - l_tumor)
    d/dt(is_tumor) <- j_tumor * vp_tumor - l_tumor * c_tumor_int / r_liver - binding +
      koff * complex + kkill * complex + proliferation
    d/dt(tumor) <- kg * (tumor - complex) - kkill * complex
    d/dt(complex) <- binding - koff * complex - kkill * complex

    # Fraction of a portal vein dose entering the tumor vasculature (the rest
    # enters the normal liver vasculature), by portal blood flow.
    f_pv_tumor <- (q_tumor_pv_gi + q_tumor_pv_spleen) / (q_gi - l_gi + q_spleen - l_spleen)

    # Outputs. Concentrations are cells/mL (equal to cells/g at unit density).
    Cblood <- a_venous / v_venous
    Ctumor <- (vp_tumor + is_tumor + complex) / v_tumor
    Vtumor <- v_tumor
  })
}
