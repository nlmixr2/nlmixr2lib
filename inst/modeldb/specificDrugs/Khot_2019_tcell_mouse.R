Khot_2019_tcell_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, ADAPT 5). Disposition of",
    "intravenously administered, ex vivo activated and chromium-51-labelled",
    "T cells in B16-BL6 melanoma-bearing C57BL/6 mice (Khot 2019). Whole",
    "blood, 13 organs and an other-tissues compartment, each split into a",
    "vascular and an extravascular subcompartment, plus a subcutaneous tumor,",
    "a whole-body lymph node and two inguinal lymph nodes (tumor-draining and",
    "non-draining). T cells leave each vascular space by irreversible",
    "first-order transmigration, return from the extravascular space with",
    "the lymph flow (slowed by retention factors in spleen, liver and",
    "kidney), and are eliminated only from the lung vascular space. Bone and",
    "muscle lymph is split equally between the two inguinal nodes and tumor",
    "lymph drains to the tumor-draining node. Physiology for a 28 g mouse",
    "from Shah and Betts 2012. Deterministic typical-value model;",
    "concentrations in cells/mL, equal to percent of injected dose per gram",
    "of tissue when the dose is 100.",
    sep = " "
  )

  reference <- paste(
    "Khot A, Matsueda S, Thomas VA, Koya RC, Shah DK. Measurement and",
    "quantitative characterization of whole-body pharmacokinetics of",
    "exogenously administered T cells in mice. J Pharmacol Exp Ther.",
    "2019;368(3):503-513. doi:10.1124/jpet.118.252858. Physiology from",
    "Shah DK, Betts AM. Towards a platform PBPK model to characterize the",
    "plasma and tissue disposition of monoclonal antibodies in preclinical",
    "species and human. J Pharmacokinet Pharmacodyn. 2012;39(1):67-86.",
    "doi:10.1007/s10928-011-9232-2 (Electronic Supplementary Material 1,",
    "the ADAPT model code).",
    sep = " "
  )

  vignette <- "Khot_2019_tcell_mouse"

  units <- list(
    time = "h",
    dosing = "cells",
    concentration = "cells/mL"
  )

  # The two inguinal lymph nodes are transit nodes specific to this model
  # (left inguinal = non-draining, right inguinal = tumor-draining); the
  # whole-body lymph node keeps the registered `lnode` name.
  paper_specific_compartments <- c("lnode_inguinal", "lnode_tumor_draining")

  compartmentData <- list(
    blood = list(analyte = "T cells", units = "cells", specimen = "whole blood", verified = TRUE),
    vp_lung = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_lung = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_heart = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_heart = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_kidney = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_kidney = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_brain = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_brain = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_muscle = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_muscle = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_bone = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_bone = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_skin = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_skin = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_adipose = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_adipose = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_small_intestine = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_small_intestine = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_large_intestine = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_large_intestine = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_spleen = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_spleen = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_pancreas = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_pancreas = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_liver = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_other = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    is_other = list(analyte = "T cells", units = "cells", specimen = "tissue", verified = TRUE),
    vp_tumor = list(analyte = "T cells", units = "cells", specimen = "tumor", verified = TRUE),
    is_tumor = list(analyte = "T cells", units = "cells", specimen = "tumor", verified = TRUE),
    lnode = list(analyte = "T cells", units = "cells", specimen = "lymph", verified = TRUE),
    lnode_inguinal = list(analyte = "T cells", units = "cells", specimen = "lymph", verified = TRUE),
    lnode_tumor_draining = list(analyte = "T cells", units = "cells", specimen = "lymph", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6, male, B16-BL6 melanoma-bearing)",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "8 weeks at purchase; tumors grown to about 250 mm^3 before dosing",
    weight_range = "28 g reference mouse (physiology from Shah and Betts 2012); body weights not reported",
    sex_female_pct = 0,
    disease_state = paste(
      "Subcutaneous B16-BL6 melanoma in the right flank (0.5 million cells",
      "inoculated; about 250 mm^3 at dosing; model tumor volume 0.45 mL)."
    ),
    dose_range = paste(
      "Single intravenous (penile vein) dose of 10 uCi of Cr-51-labelled",
      "activated splenic T cells, about 2 x 10^7 cells per mouse."
    ),
    notes = paste(
      "Three mice sacrificed at each of 1 h and 1, 3, 7 and 14 days; 19",
      "tissues counted and expressed as percent of injected dose per gram.",
      "T cells were isolated from the spleens of tumor-bearing mice and",
      "activated with anti-CD3, anti-CD28 and IL-2 (more than 90 percent",
      "CD3+, of which more than 85 percent CD8+). The model was fit by",
      "maximum likelihood (naive pooled) in ADAPT V."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Transmigration rate constants (Khot 2019 Table 2), first order on the
    # vascular amount of each tissue
    # ---------------------------------------------------------------------
    j_lung <- 1843
    label("Lung transmigration rate constant (1/h)") # Table 2 J_lung = 1843.0 (CV 15.4)
    j_heart <- 34.9
    label("Heart transmigration rate constant (1/h)") # Table 2 J_heart = 34.9 (CV 17.2)
    j_kidney <- 87.3
    label("Kidney transmigration rate constant (1/h)") # Table 2 J_kidney = 87.3 (CV 37.9)
    j_brain <- 1.4
    label("Brain transmigration rate constant (1/h)") # Table 2 J_brain = 1.4 (CV 21.2)
    j_muscle <- 0.5
    label("Muscle transmigration rate constant (1/h)") # Table 2 J_muscle = 0.5 (CV 19.7)
    j_bone <- 82.6
    label("Bone transmigration rate constant (1/h)") # Table 2 J_bone = 82.6 (CV 16.7)
    j_tumor <- 0.6
    label("Tumor transmigration rate constant (1/h)") # Table 2 J_tumor = 0.6 (CV 18.6)
    j_skin <- 0.6
    label("Skin transmigration rate constant (1/h)") # Table 2 J_skin = 0.6 (CV 18.7)
    j_adipose <- 1.6
    label("Adipose transmigration rate constant (1/h)") # Table 2 J_fat = 1.6 (CV 18.5)
    j_small_intestine <- 12.9
    label("Small-intestine transmigration rate constant (1/h)") # Table 2 J_SI = 12.9 (CV 17.1)
    j_large_intestine <- 5.2
    label("Large-intestine transmigration rate constant (1/h)") # Table 2 J_LI = 5.2 (CV 18.1)
    j_spleen <- 114
    label("Spleen transmigration rate constant (1/h)") # Table 2 J_spleen = 114.0 (CV 33.8)
    j_liver <- 126.9
    label("Liver transmigration rate constant (1/h)") # Table 2 J_liver = 126.9 (CV 18.9)
    j_pancreas <- 10
    label("Pancreas transmigration rate constant (1/h)") # Table 2 J_pancreas = 10.0 (CV 20.0)
    j_other <- 86.8
    label("Other-tissues transmigration rate constant (1/h)") # Table 2 J_other = 86.8 (CV 17.7)

    # ---------------------------------------------------------------------
    # Retention factors (Table 2): divide the lymphatic return of T cells
    # from the extravascular space of spleen, liver and kidney
    # ---------------------------------------------------------------------
    r_kidney <- 3.9
    label("Kidney retention factor (unitless)") # Table 2 R_kidney = 3.9 (CV 37.9)
    r_spleen <- 9.8
    label("Spleen retention factor (unitless)") # Table 2 R_spleen = 9.8 (CV 34.8)
    r_liver <- 2.5
    label("Liver retention factor (unitless)") # Table 2 R_liver = 2.5 (CV 18.4)

    # ---------------------------------------------------------------------
    # Values held constant by the authors (Table 2)
    # ---------------------------------------------------------------------
    lkel <- fixed(log(0.84))
    label("Elimination rate constant from the lung vascular space (1/h)") # Table 2 E_lung = 0.84, taken from Zhu et al. 1996
    f_recirc <- fixed(1)
    label("Fraction of T cells leaving each extravascular space that recirculates to a lymph node (fraction)") # Table 2 F_tissue = 1, an assumption of the authors
    q_tumor <- fixed(6)
    label("Tumor blood flow (mL/h)") # Table 2 Q_tumor = 6, taken from Zhu et al. 1996
    v_tumor <- fixed(0.45)
    label("Tumor volume (mL)") # Table 2 V_tumor = 0.45 mL or g (radius about 3.85 mm)

    # ---------------------------------------------------------------------
    # Tumor sub-volumes and lymph flow: not reported; derived by the
    # maintainers within the Shah and Betts 2012 framework (vignette
    # Assumptions)
    # ---------------------------------------------------------------------
    f_vasc_tumor <- fixed(0.068)
    label("Vascular fraction of the tumor volume (fraction)") # median (VV + BV) / total over the 13 Shah 2012 organs used here; not reported by Khot 2019
    f_isf_tumor <- fixed(0.174)
    label("Extravascular (interstitial) fraction of the tumor volume (fraction)") # median IV / total over the 13 Shah 2012 organs used here; not reported by Khot 2019

    # ---------------------------------------------------------------------
    # Inguinal lymph-node volumes: not reported anywhere; assumed by the
    # maintainers (vignette Assumptions). The node concentrations are
    # insensitive to them because each node turns over in well under an hour.
    # ---------------------------------------------------------------------
    v_lnode_inguinal <- fixed(0.005)
    label("Non-draining (left) inguinal lymph-node volume (mL)") # assumed 5 mg node; not reported by Khot 2019
    v_lnode_tumor_draining <- fixed(0.005)
    label("Tumor-draining (right) inguinal lymph-node volume (mL)") # assumed 5 mg node; not reported by Khot 2019

    # ---------------------------------------------------------------------
    # Physiology of a 28 g mouse (Shah and Betts 2012 ESM1, ADAPT code;
    # values in L and L/h there, mL and mL/h here). Blood flow = plasma +
    # blood-cell flow; vascular volume = vascular plasma + vascular
    # blood-cell volume; extravascular volume = endosomal + interstitial +
    # cellular volume (the whole tissue outside the vasculature).
    # ---------------------------------------------------------------------
    v_blood <- fixed(1.717)
    label("Whole-blood volume (mL)") # Shah 2012 ESM1: PL_V_M 0.00094435 + BC_V_M 0.00077265 L
    lv_lnode <- fixed(log(0.113))
    label("Whole-body lymph-node volume (mL)") # Shah 2012 ESM1: LN_V_M 0.000113 L
    f_lymph <- fixed(0.002)
    label("Lymph flow as a fraction of tissue plasma flow (fraction)") # Shah 2012 ESM1: <tissue>_LF_M = <tissue>_PLQ_M*0.002
    f_plasma <- fixed(0.55)
    label("Plasma fraction of tissue blood flow, one minus the hematocrit (fraction)") # Shah 2012 ESM1: PLQ / (PLQ + BCQ) = 0.55 for every tissue
    q_lung <- fixed(678.42)
    label("Lung blood flow (mL/h)") # Shah 2012 ESM1: LU_PLQ_M 0.373131 + LU_BCQ_M 0.305289 L/h (plasma + blood-cell flow)
    vv_lung <- fixed(0.0535547)
    label("Lung vascular volume (mL)") # Shah 2012 ESM1: LU_VV_M 0.00002945509193 + LU_BV_M 0.00002409962067 L
    vi_lung <- fixed(0.150853)
    label("Lung extravascular volume (mL)") # Shah 2012 ESM1: LU_EV_M 0.0000010220365 + LU_IV_M 0.0000384285724 + LU_CV_M 0.000111402 L (endosomal + interstitial + cellular)
    vt_lung <- fixed(0.204407)
    label("Lung total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of LU_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_heart <- fixed(66.36)
    label("Heart blood flow (mL/h)") # Shah 2012 ESM1: HE_PLQ_M 0.036498 + HE_BCQ_M 0.029862 L/h (plasma + blood-cell flow)
    vv_heart <- fixed(0.01064)
    label("Heart vascular volume (mL)") # Shah 2012 ESM1: HE_VV_M 0.000005852 + HE_BV_M 0.000004788 L
    vi_heart <- fixed(0.14136)
    label("Heart extravascular volume (mL)") # Shah 2012 ESM1: HE_EV_M 0.00000076 + HE_IV_M 0.000021736 + HE_CV_M 0.000118864 L (endosomal + interstitial + cellular)
    vt_heart <- fixed(0.152)
    label("Heart total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of HE_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_kidney <- fixed(124.56)
    label("Kidney blood flow (mL/h)") # Shah 2012 ESM1: KI_PLQ_M 0.068508 + KI_BCQ_M 0.056052 L/h (plasma + blood-cell flow)
    vv_kidney <- fixed(0.0525)
    label("Kidney vascular volume (mL)") # Shah 2012 ESM1: KI_VV_M 0.000028875 + KI_BV_M 0.000023625 L
    vi_kidney <- fixed(0.4725)
    label("Kidney extravascular volume (mL)") # Shah 2012 ESM1: KI_EV_M 0.000002625 + KI_IV_M 0.00007875 + KI_CV_M 0.000391125 L (endosomal + interstitial + cellular)
    vt_kidney <- fixed(0.525)
    label("Kidney total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of KI_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_brain <- fixed(21.42)
    label("Brain blood flow (mL/h)") # Shah 2012 ESM1: BR_PLQ_M 0.011781 + BR_BCQ_M 0.009639 L/h (plasma + blood-cell flow)
    vv_brain <- fixed(0.0194)
    label("Brain vascular volume (mL)") # Shah 2012 ESM1: BR_VV_M 0.00001067 + BR_BV_M 0.00000873 L
    vi_brain <- fixed(0.4656)
    label("Brain extravascular volume (mL)") # Shah 2012 ESM1: BR_EV_M 0.000002425 + BR_IV_M 0.0000873 + BR_CV_M 0.000375875 L (endosomal + interstitial + cellular)
    vt_brain <- fixed(0.485)
    label("Brain total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of BR_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_muscle <- fixed(156.6)
    label("Muscle blood flow (mL/h)") # Shah 2012 ESM1: MU_PLQ_M 0.08613 + MU_BCQ_M 0.07047 L/h (plasma + blood-cell flow)
    vv_muscle <- fixed(0.45276)
    label("Muscle vascular volume (mL)") # Shah 2012 ESM1: MU_VV_M 0.000249018 + MU_BV_M 0.000203742 L
    vi_muscle <- fixed(10.8662)
    label("Muscle extravascular volume (mL)") # Shah 2012 ESM1: MU_EV_M 0.000056595 + MU_IV_M 0.00147147 + MU_CV_M 0.009338175 L (endosomal + interstitial + cellular)
    vt_muscle <- fixed(11.319)
    label("Muscle total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of MU_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_bone <- fixed(27.6)
    label("Bone blood flow (mL/h)") # Shah 2012 ESM1: BO_PLQ_M 0.01518 + BO_BCQ_M 0.01242 L/h (plasma + blood-cell flow)
    vv_bone <- fixed(0.11296)
    label("Bone vascular volume (mL)") # Shah 2012 ESM1: BO_VV_M 0.000062128 + BO_BV_M 0.000050832 L
    vi_bone <- fixed(2.71104)
    label("Bone extravascular volume (mL)") # Shah 2012 ESM1: BO_EV_M 0.00001412 + BO_IV_M 0.000525264 + BO_CV_M 0.002171656 L (endosomal + interstitial + cellular)
    vt_bone <- fixed(2.824)
    label("Bone total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of BO_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_skin <- fixed(50.58)
    label("Skin blood flow (mL/h)") # Shah 2012 ESM1: SK_PLQ_M 0.027819 + SK_BCQ_M 0.022761 L/h (plasma + blood-cell flow)
    vv_skin <- fixed(0.341292)
    label("Skin vascular volume (mL)") # Shah 2012 ESM1: SK_VV_M 0.0001877106 + SK_BV_M 0.0001535814 L
    vi_skin <- fixed(4.67771)
    label("Skin extravascular volume (mL)") # Shah 2012 ESM1: SK_EV_M 0.000025095 + SK_IV_M 0.00165627 + SK_CV_M 0.002996343 L (endosomal + interstitial + cellular)
    vt_skin <- fixed(5.019)
    label("Skin total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of SK_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_adipose <- fixed(24.42)
    label("Adipose blood flow (mL/h)") # Shah 2012 ESM1: AD_PLQ_M 0.013431 + AD_BCQ_M 0.010989 L/h (plasma + blood-cell flow)
    vv_adipose <- fixed(0.03964)
    label("Adipose vascular volume (mL)") # Shah 2012 ESM1: AD_VV_M 0.000021802 + AD_BV_M 0.000017838 L
    vi_adipose <- fixed(1.94236)
    label("Adipose extravascular volume (mL)") # Shah 2012 ESM1: AD_EV_M 0.00000991 + AD_IV_M 0.00033694 + AD_CV_M 0.00159551 L (endosomal + interstitial + cellular)
    vt_adipose <- fixed(1.982)
    label("Adipose total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of AD_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_small_intestine <- fixed(105.6)
    label("Small-intestine blood flow (mL/h)") # Shah 2012 ESM1: SI_PLQ_M 0.05808 + SI_BCQ_M 0.04752 L/h (plasma + blood-cell flow)
    vv_small_intestine <- fixed(0.021112)
    label("Small-intestine vascular volume (mL)") # Shah 2012 ESM1: SI_VV_M 0.0000116116 + SI_BV_M 0.0000095004 L
    vi_small_intestine <- fixed(0.706888)
    label("Small-intestine extravascular volume (mL)") # Shah 2012 ESM1: SI_EV_M 0.00000364 + SI_IV_M 0.000126672 + SI_CV_M 0.000576576 L (endosomal + interstitial + cellular)
    vt_small_intestine <- fixed(0.728)
    label("Small-intestine total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of SI_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_large_intestine <- fixed(31.44)
    label("Large-intestine blood flow (mL/h)") # Shah 2012 ESM1: LI_PLQ_M 0.017292 + LI_BCQ_M 0.014148 L/h (plasma + blood-cell flow)
    vv_large_intestine <- fixed(0.0090915)
    label("Large-intestine vascular volume (mL)") # Shah 2012 ESM1: LI_VV_M 0.000005000325 + LI_BV_M 0.000004091175 L
    vi_large_intestine <- fixed(0.304408)
    label("Large-intestine extravascular volume (mL)") # Shah 2012 ESM1: LI_EV_M 0.0000015675 + LI_IV_M 0.000054549 + LI_CV_M 0.000248292 L (endosomal + interstitial + cellular)
    vt_large_intestine <- fixed(0.3135)
    label("Large-intestine total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of LI_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_spleen <- fixed(14.88)
    label("Spleen blood flow (mL/h)") # Shah 2012 ESM1: SP_PLQ_M 0.008184 + SP_BCQ_M 0.006696 L/h (plasma + blood-cell flow)
    vv_spleen <- fixed(0.02794)
    label("Spleen vascular volume (mL)") # Shah 2012 ESM1: SP_VV_M 0.000015367 + SP_BV_M 0.000012573 L
    vi_spleen <- fixed(0.09906)
    label("Spleen extravascular volume (mL)") # Shah 2012 ESM1: SP_EV_M 0.000000635 + SP_IV_M 0.0000254 + SP_CV_M 0.000073025 L (endosomal + interstitial + cellular)
    vt_spleen <- fixed(0.127)
    label("Spleen total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of SP_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_pancreas <- fixed(11.34)
    label("Pancreas blood flow (mL/h)") # Shah 2012 ESM1: PA_PLQ_M 0.006237 + PA_BCQ_M 0.005103 L/h (plasma + blood-cell flow)
    vv_pancreas <- fixed(0.0097)
    label("Pancreas vascular volume (mL)") # Shah 2012 ESM1: PA_VV_M 0.000005335 + PA_BV_M 0.000004365 L
    vi_pancreas <- fixed(0.0873)
    label("Pancreas extravascular volume (mL)") # Shah 2012 ESM1: PA_EV_M 0.000000485 + PA_IV_M 0.000016878 + PA_CV_M 0.000069937 L (endosomal + interstitial + cellular)
    vt_pancreas <- fixed(0.097)
    label("Pancreas total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of PA_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_liver <- fixed(18.66)
    label("Hepatic-artery blood flow to the liver (mL/h)") # Shah 2012 ESM1: LIV_PLQ_M 0.010263 + LIV_BCQ_M 0.008397 L/h (plasma + blood-cell flow)
    vv_liver <- fixed(0.298375)
    label("Liver vascular volume (mL)") # Shah 2012 ESM1: LIV_VV_M 0.00016410625 + LIV_BV_M 0.00013426875 L
    vi_liver <- fixed(1.62663)
    label("Liver extravascular volume (mL)") # Shah 2012 ESM1: LIV_EV_M 0.000009625 + LIV_IV_M 0.000384595163 + LIV_CV_M 0.001232405 L (endosomal + interstitial + cellular)
    vt_liver <- fixed(1.925)
    label("Liver total tissue volume, the denominator of the reported tissue concentration (mL)") # Shah 2012 ESM1: sum of LIV_VV_M, _BV_M, _EV_M, _IV_M, _CV_M
    q_other <- fixed(24.96)
    label("Blood flow to the other-tissues compartment before the tumor flow is drawn from it (mL/h)") # Shah 2012 ESM1: OT_PLQ_M 0.01254 + OT_BCQ_M 0.01026 + TH_PLQ_M 0.001188 + TH_BCQ_M 0.000972 L/h (thymus folded in)
    vv_other <- fixed(0.0364193)
    label("Other-tissues vascular volume (mL)") # Shah 2012 ESM1: OT_VV_M 0.000019535609097 + OT_BV_M 0.0000159837 + TH_VV_M 0.000000495 + TH_BV_M 0.000000405 L
    vi_other <- fixed(0.437899)
    label("Other-tissues extravascular volume (mL)") # Shah 2012 ESM1: OT_EV_M + OT_IV_M + OT_CV_M + TH_EV_M + TH_IV_M + TH_CV_M = 0.000002326591 + 0.0000797245183 + 0.000347748 + 0.000000045 + 0.00000153 + 0.000006525 L
  })

  model({
    kel <- exp(lkel)
    v_lnode <- exp(lv_lnode)

    # Lymph flows: f_lymph of each tissue's plasma flow (Shah 2012 ESM1). The
    # liver lymph also carries the portal plasma arriving from the gut organs
    # and spleen, exactly as LIV_LF_M is defined there.
    l_lung <- f_lymph * f_plasma * q_lung
    l_heart <- f_lymph * f_plasma * q_heart
    l_kidney <- f_lymph * f_plasma * q_kidney
    l_brain <- f_lymph * f_plasma * q_brain
    l_muscle <- f_lymph * f_plasma * q_muscle
    l_bone <- f_lymph * f_plasma * q_bone
    l_skin <- f_lymph * f_plasma * q_skin
    l_adipose <- f_lymph * f_plasma * q_adipose
    l_small_intestine <- f_lymph * f_plasma * q_small_intestine
    l_large_intestine <- f_lymph * f_plasma * q_large_intestine
    l_spleen <- f_lymph * f_plasma * q_spleen
    l_pancreas <- f_lymph * f_plasma * q_pancreas
    q_portal <- q_small_intestine + q_large_intestine + q_spleen + q_pancreas
    l_liver <- f_lymph * f_plasma * (q_liver + (1 - f_lymph) * q_portal)

    # Tumor: blood flow drawn from the other-tissues compartment so that the
    # arterial outflow of the lung equals the summed tissue inflows.
    q_other_net <- q_other - q_tumor
    l_other <- f_lymph * f_plasma * q_other_net
    l_tumor <- f_lymph * f_plasma * q_tumor
    vv_tumor <- f_vasc_tumor * v_tumor
    vi_tumor <- f_isf_tumor * v_tumor

    # Inguinal nodes pass on the lymph they receive (Eqs 11-12); the
    # whole-body node returns all lymph to blood (Eqs 1 and 10).
    l_lnode_inguinal <- (l_bone + l_muscle) / 2
    l_lnode_tumor_draining <- (l_bone + l_muscle) / 2 + l_tumor
    l_lnode <- l_lung + l_heart + l_kidney + l_brain + l_skin + l_adipose +
      l_small_intestine + l_large_intestine + l_spleen + l_pancreas + l_liver +
      l_other + l_lnode_inguinal + l_lnode_tumor_draining

    # Liver vascular outflow to blood: portal inflow net of the gut lymph plus
    # the hepatic-artery inflow net of the liver lymph (Eqs 1 and 8).
    q_liver_out <- (q_small_intestine - l_small_intestine) +
      (q_large_intestine - l_large_intestine) + (q_spleen - l_spleen) +
      (q_pancreas - l_pancreas) + (q_liver - l_liver)

    c_blood <- blood / v_blood
    c_vp_lung <- vp_lung / vv_lung
    c_is_lung <- is_lung / vi_lung
    c_vp_heart <- vp_heart / vv_heart
    c_is_heart <- is_heart / vi_heart
    c_vp_kidney <- vp_kidney / vv_kidney
    c_is_kidney <- is_kidney / vi_kidney
    c_vp_brain <- vp_brain / vv_brain
    c_is_brain <- is_brain / vi_brain
    c_vp_muscle <- vp_muscle / vv_muscle
    c_is_muscle <- is_muscle / vi_muscle
    c_vp_bone <- vp_bone / vv_bone
    c_is_bone <- is_bone / vi_bone
    c_vp_skin <- vp_skin / vv_skin
    c_is_skin <- is_skin / vi_skin
    c_vp_adipose <- vp_adipose / vv_adipose
    c_is_adipose <- is_adipose / vi_adipose
    c_vp_small_intestine <- vp_small_intestine / vv_small_intestine
    c_is_small_intestine <- is_small_intestine / vi_small_intestine
    c_vp_large_intestine <- vp_large_intestine / vv_large_intestine
    c_is_large_intestine <- is_large_intestine / vi_large_intestine
    c_vp_spleen <- vp_spleen / vv_spleen
    c_is_spleen <- is_spleen / vi_spleen
    c_vp_pancreas <- vp_pancreas / vv_pancreas
    c_is_pancreas <- is_pancreas / vi_pancreas
    c_vp_liver <- vp_liver / vv_liver
    c_is_liver <- is_liver / vi_liver
    c_vp_other <- vp_other / vv_other
    c_is_other <- is_other / vi_other
    c_vp_tumor <- vp_tumor / vv_tumor
    c_is_tumor <- is_tumor / vi_tumor
    c_lnode <- lnode / v_lnode
    c_lnode_inguinal <- lnode_inguinal / v_lnode_inguinal
    c_lnode_tumor_draining <- lnode_tumor_draining / v_lnode_tumor_draining

    # Eq 1, whole blood. The tumor venous return is not printed in Eq 1 but
    # is required for mass balance (vignette Assumptions).
    d/dt(blood) <- -(q_lung + l_lung) * c_blood +
      (q_heart - l_heart) * c_vp_heart + (q_kidney - l_kidney) * c_vp_kidney +
      (q_brain - l_brain) * c_vp_brain + (q_muscle - l_muscle) * c_vp_muscle +
      (q_bone - l_bone) * c_vp_bone + (q_skin - l_skin) * c_vp_skin +
      (q_adipose - l_adipose) * c_vp_adipose + q_liver_out * c_vp_liver +
      (q_other_net - l_other) * c_vp_other + (q_tumor - l_tumor) * c_vp_tumor +
      l_lnode * c_lnode

    # Eqs 2-3, lungs: elimination from the vascular space only
    d/dt(vp_lung) <- (q_lung + l_lung) * c_blood - q_lung * c_vp_lung -
      j_lung * vp_lung - kel * vp_lung
    d/dt(is_lung) <- j_lung * vp_lung - l_lung * f_recirc * c_is_lung

    # Eqs 4-5, typical tissues (arterial supply at the lung vascular
    # concentration)
    d/dt(vp_heart) <- q_heart * c_vp_lung - (q_heart - l_heart) * c_vp_heart - j_heart * vp_heart
    d/dt(is_heart) <- j_heart * vp_heart - l_heart * f_recirc * c_is_heart
    d/dt(vp_brain) <- q_brain * c_vp_lung - (q_brain - l_brain) * c_vp_brain - j_brain * vp_brain
    d/dt(is_brain) <- j_brain * vp_brain - l_brain * f_recirc * c_is_brain
    d/dt(vp_muscle) <- q_muscle * c_vp_lung - (q_muscle - l_muscle) * c_vp_muscle - j_muscle * vp_muscle
    d/dt(is_muscle) <- j_muscle * vp_muscle - l_muscle * f_recirc * c_is_muscle
    d/dt(vp_bone) <- q_bone * c_vp_lung - (q_bone - l_bone) * c_vp_bone - j_bone * vp_bone
    d/dt(is_bone) <- j_bone * vp_bone - l_bone * f_recirc * c_is_bone
    d/dt(vp_skin) <- q_skin * c_vp_lung - (q_skin - l_skin) * c_vp_skin - j_skin * vp_skin
    d/dt(is_skin) <- j_skin * vp_skin - l_skin * f_recirc * c_is_skin
    d/dt(vp_adipose) <- q_adipose * c_vp_lung - (q_adipose - l_adipose) * c_vp_adipose - j_adipose * vp_adipose
    d/dt(is_adipose) <- j_adipose * vp_adipose - l_adipose * f_recirc * c_is_adipose
    d/dt(vp_small_intestine) <- q_small_intestine * c_vp_lung -
      (q_small_intestine - l_small_intestine) * c_vp_small_intestine -
      j_small_intestine * vp_small_intestine
    d/dt(is_small_intestine) <- j_small_intestine * vp_small_intestine -
      l_small_intestine * f_recirc * c_is_small_intestine
    d/dt(vp_large_intestine) <- q_large_intestine * c_vp_lung -
      (q_large_intestine - l_large_intestine) * c_vp_large_intestine -
      j_large_intestine * vp_large_intestine
    d/dt(is_large_intestine) <- j_large_intestine * vp_large_intestine -
      l_large_intestine * f_recirc * c_is_large_intestine
    d/dt(vp_pancreas) <- q_pancreas * c_vp_lung - (q_pancreas - l_pancreas) * c_vp_pancreas - j_pancreas * vp_pancreas
    d/dt(is_pancreas) <- j_pancreas * vp_pancreas - l_pancreas * f_recirc * c_is_pancreas
    d/dt(vp_other) <- q_other_net * c_vp_lung - (q_other_net - l_other) * c_vp_other - j_other * vp_other
    d/dt(is_other) <- j_other * vp_other - l_other * f_recirc * c_is_other
    d/dt(vp_tumor) <- q_tumor * c_vp_lung - (q_tumor - l_tumor) * c_vp_tumor - j_tumor * vp_tumor
    d/dt(is_tumor) <- j_tumor * vp_tumor - l_tumor * f_recirc * c_is_tumor

    # Eqs 6-7, tissues with a retention factor
    d/dt(vp_kidney) <- q_kidney * c_vp_lung - (q_kidney - l_kidney) * c_vp_kidney - j_kidney * vp_kidney
    d/dt(is_kidney) <- j_kidney * vp_kidney - l_kidney * c_is_kidney / r_kidney * f_recirc
    d/dt(vp_spleen) <- q_spleen * c_vp_lung - (q_spleen - l_spleen) * c_vp_spleen - j_spleen * vp_spleen
    d/dt(is_spleen) <- j_spleen * vp_spleen - l_spleen * c_is_spleen / r_spleen * f_recirc

    # Eqs 8-9, liver: hepatic artery plus portal inflow
    d/dt(vp_liver) <- q_liver * c_vp_lung +
      (q_small_intestine - l_small_intestine) * c_vp_small_intestine +
      (q_large_intestine - l_large_intestine) * c_vp_large_intestine +
      (q_spleen - l_spleen) * c_vp_spleen + (q_pancreas - l_pancreas) * c_vp_pancreas -
      q_liver_out * c_vp_liver - j_liver * vp_liver
    d/dt(is_liver) <- j_liver * vp_liver - l_liver * c_is_liver / r_liver * f_recirc

    # Eq 10, whole-body lymph node
    d/dt(lnode) <- f_recirc * (l_lung * c_is_lung + l_heart * c_is_heart +
      l_kidney * c_is_kidney / r_kidney + l_brain * c_is_brain +
      l_skin * c_is_skin + l_adipose * c_is_adipose +
      l_small_intestine * c_is_small_intestine +
      l_large_intestine * c_is_large_intestine + l_spleen * c_is_spleen / r_spleen +
      l_pancreas * c_is_pancreas + l_liver * c_is_liver / r_liver +
      l_other * c_is_other) +
      l_lnode_inguinal * c_lnode_inguinal +
      l_lnode_tumor_draining * c_lnode_tumor_draining - l_lnode * c_lnode

    # Eqs 11-12, inguinal lymph nodes: half of the bone and muscle lymph each;
    # the tumor lymph drains to the right (tumor-draining) node only
    d/dt(lnode_inguinal) <- f_recirc * (l_bone / 2 * c_is_bone + l_muscle / 2 * c_is_muscle) -
      l_lnode_inguinal * c_lnode_inguinal
    d/dt(lnode_tumor_draining) <- f_recirc * (l_bone / 2 * c_is_bone + l_muscle / 2 * c_is_muscle +
      l_tumor * c_is_tumor) - l_lnode_tumor_draining * c_lnode_tumor_draining

    # Outputs: T cells per mL of whole tissue (cells/g at unit density). Each
    # organ's subcompartment amounts are divided by the whole organ volume,
    # matching the per-gram tissue counting of the study.
    Cblood <- c_blood
    Clung <- (vp_lung + is_lung) / vt_lung
    Cheart <- (vp_heart + is_heart) / vt_heart
    Ckidney <- (vp_kidney + is_kidney) / vt_kidney
    Cbrain <- (vp_brain + is_brain) / vt_brain
    Cmuscle <- (vp_muscle + is_muscle) / vt_muscle
    Cbone <- (vp_bone + is_bone) / vt_bone
    Cskin <- (vp_skin + is_skin) / vt_skin
    Cadipose <- (vp_adipose + is_adipose) / vt_adipose
    Csmall_intestine <- (vp_small_intestine + is_small_intestine) / vt_small_intestine
    Clarge_intestine <- (vp_large_intestine + is_large_intestine) / vt_large_intestine
    Cspleen <- (vp_spleen + is_spleen) / vt_spleen
    Cpancreas <- (vp_pancreas + is_pancreas) / vt_pancreas
    Cliver <- (vp_liver + is_liver) / vt_liver
    Ctumor <- (vp_tumor + is_tumor) / v_tumor
    Clnode_inguinal <- c_lnode_inguinal
    Clnode_tumor_draining <- c_lnode_tumor_draining
  })
}
