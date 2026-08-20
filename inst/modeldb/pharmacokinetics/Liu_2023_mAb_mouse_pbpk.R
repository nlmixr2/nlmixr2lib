Liu_2023_mAb_mouse_pbpk <- function() {
  description <- "PBPK (whole-body, 15-tissue Shah & Betts 2012 platform). Preclinical (mouse, female BALB/c, 28 g). Inter-antibody variability in monoclonal antibody plasma PK, with the antibody-specific pinocytosis-uptake coefficient F1 driven by a sigmoidal function of the antibody's heparin-chromatography retention time (HEPARIN_RT) and the convective-transport coefficient F2 fixed at the panel median. Fitted simultaneously to IV plasma profiles of 53 aglycosylated human IgG1z SEFL2.2 antibodies dosed at 2 mg/kg."
  reference <- "Liu S, Humphreys SC, Cook KD, Conner KP, Correia AR, Jacobitz AW, Yang M, Primack R, Soto M, Padaki R, Lubomirski M, Smith R, Mock M, Thomas VA. Utility of physiologically based pharmacokinetic modeling to predict inter-antibody variability in monoclonal antibody pharmacokinetics in mice. MAbs. 2023;15(1):2263926. doi:10.1080/19420862.2023.2263926"
  vignette <- "Liu_2023_mAb_mouse_pbpk"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  compartmentData <- list(
    # Tissue vascular (plasma) space: mAb amount
    vp_heart           = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_lung            = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_muscle          = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_skin            = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_adipose         = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_bone            = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_brain           = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_kidney          = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_liver           = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_thymus          = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_spleen          = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    vp_other           = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    # Tissue interstitial space: mAb amount
    is_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_thymus          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    is_other           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = TRUE),
    # Endosomal space: FcRn-unbound mAb amount
    eu_heart           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_lung            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_muscle          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_skin            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_adipose         = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_bone            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_brain           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_kidney          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_liver           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_thymus          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_spleen          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eu_other           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    # Endosomal space: FcRn-bound mAb amount
    eb_heart           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_lung            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_muscle          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_skin            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_adipose         = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_bone            = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_brain           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_kidney          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_liver           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_thymus          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_spleen          = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    eb_other           = list(analyte = "mAb", units = "nmol", specimen = "endosome", verified = TRUE),
    # Endosomal space: free (unoccupied) FcRn amount
    fr_heart           = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_lung            = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_muscle          = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_skin            = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_adipose         = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_bone            = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_brain           = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_kidney          = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_liver           = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_small_intestine = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_large_intestine = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_pancreas        = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_thymus          = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_spleen          = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    fr_other           = list(analyte = "FcRn", units = "nmol", specimen = "endosome", verified = TRUE),
    # Central plasma and lymph node
    plasma             = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode              = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = TRUE)
  )

  covariateData <- list(
    HEPARIN_RT = list(
      description        = "Heparin-chromatography retention time of the administered monoclonal antibody",
      units              = "min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Molecule-level in-vitro physicochemical attribute, not a subject-level covariate:",
        "every record for a given animal carries the retention time of the antibody that animal received.",
        "Measured on a HiTrap Heparin High Performance 1 mL column (50 mM Tris pH 7.6, 5 mM NaCl loading",
        "buffer; linear 5-400 mM NaCl gradient over 20 column volumes), recorded at the centre of the",
        "elution peak (Liu 2023 Supplementary Methods 1, Heparin Chromatography Assay).",
        "Training-set range 2.92-31.5 min, validation-set range 15.1-32.6 min",
        "(Liu 2023 Supplementary Table S1). The covariate enters the model only through the sigmoidal",
        "F1 relationship; the paper's threshold for flagging a fast-clearing antibody is 16.5 min."
      ),
      source_name        = "Heparin_RT"
    )
  )

  population <- list(
    species        = "mouse (female BALB/c, 28 g)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "6-8 weeks at dosing",
    weight_range   = "28 g (Shah & Betts 2012 mouse physiological parameter set)",
    sex_female_pct = 100,
    disease_state  = "Healthy wild-type mice; no mouse-cross-reactive target, so no target-mediated drug disposition.",
    dose_range     = "2 mg/kg IV bolus per antibody (lateral tail vein), administered as cassettes of up to 5 antibodies per animal with a combined antibody load of 10 mg/kg.",
    notes          = paste(
      "Panel of 83 aglycosylated human IgG1z stable-effector-functionless (SEFL) 2.2 antibodies",
      "(N297G plus R292C/V302C engineered disulfide). 19 antibodies with mouse cross-reactivity and",
      "8 with observed immunogenicity were excluded, leaving 56 for individual F1/F2 estimation.",
      "Ab76 (atypical profile), Ab37 and Ab30 (nonspecific-binding assay false positives) were further",
      "removed, leaving 53 antibodies in the simultaneous covariate fit reported in Table 2.",
      "An independent 14-antibody validation panel (Ab85-Ab98) was predicted a priori from HEPARIN_RT.",
      "Serum concentrations were quantified by multiplex electrochemiluminescent immunoassay",
      "(n = 3 animals per time point, mean profiles used for fitting)."
    ),
    scope_note     = paste(
      "Mechanistic platform model: no between-subject random effects. Inter-antibody variability is",
      "carried entirely by the HEPARIN_RT covariate on F1. Intended for typical-value simulation of a",
      "named antibody's mouse plasma profile given its heparin retention time."
    )
  )

  ini({
    # --- Inter-antibody variability coefficients (Liu 2023) -------------------
    # Named for the quantity each one scales rather than for the paper's symbol.
    # `clup_scale` is the paper's F1: it multiplies pinocytotic uptake (CLup).
    # `lymph_scale` is the paper's F2: it multiplies organ lymphatic flow (L).
    # clup_scale = clup_scale_min +
    #   (clup_scale_max - clup_scale_min) /
    #   (1 + (rt50_clup_scale / HEPARIN_RT) ^ hill_clup_scale)
    # (Liu 2023 Eq 1 / Eq 2; coefficients a, b, c, d of Table 2).
    lclup_scale_min  <- log(0.67);  label("Pinocytotic-uptake scaling, lower plateau (low-heparin-retention antibodies; paper's F1 coefficient a, unitless)")   # Liu 2023 Table 2, coefficient a = 0.67 (CV% 1.28)
    lclup_scale_max  <- log(3.65);  label("Pinocytotic-uptake scaling, upper plateau (high-heparin-retention antibodies; paper's F1 coefficient b, unitless)")  # Liu 2023 Table 2, coefficient b = 3.65 (CV% 3.37)
    lrt50_clup_scale <- log(17.7);  label("HEPARIN_RT giving half-maximal pinocytotic-uptake scaling (paper's F1 coefficient c, min)")                          # Liu 2023 Table 2, coefficient c = 17.7 (CV% 0.543)
    lhill_clup_scale <- log(22.7);  label("Sigmoidicity of the HEPARIN_RT to pinocytotic-uptake-scaling relationship (paper's F1 coefficient d, unitless)")     # Liu 2023 Table 2, coefficient d = 22.7 (CV% 7.89)
    llymph_scale     <- fixed(log(1.45)); label("Organ-lymphatic-flow scaling coefficient (paper's F2, unitless)")                                              # Liu 2023 Results: "we fixed F2 to 1.45, the median F2 value in this panel"

    # --- System parameters -----------------------------------------------------
    # CLup is the one platform parameter Liu 2023 changed relative to Shah & Betts
    # 2012 ("The parameter glossary and values are the same as in Shah and Betts
    # except that the CLup value is set to 0.24 L/h/L").
    lclup  <- fixed(log(0.24));    label("Pinocytosis rate per unit endosomal volume (L/h/L)")   # Liu 2023 Materials and methods, base PBPK model

    # Remaining system parameters inherited from the on-disk upstream platform paper
    # (Shah DK, Betts AM. J Pharmacokinet Pharmacodyn. 2012;39(1):67-86, Table 6 and text p.73),
    # which Liu 2023 cites as reference 36 and adopts unchanged.
    lfcrn  <- fixed(log(4.98e-5)); label("Endosomal FcRn concentration (mol/L)")                 # inherited from Shah & Betts 2012 Table 6 (CV% 11.1)
    lkdeg  <- fixed(log(42.9));    label("Endosomal degradation rate of FcRn-unbound mAb (1/h)") # inherited from Shah & Betts 2012 Table 6 (CV% 15.7)
    lclnlf <- fixed(log(9.1));     label("Lymph-node-to-plasma transfer scaling C_LNLF (unitless)")         # inherited from Shah & Betts 2012 Table 6 (CV% > 50)
    lkon   <- fixed(log(8.06e7));  label("FcRn-IgG association rate constant (1/M/h, mouse)")    # inherited from Shah & Betts 2012 text p.73, mouse value
    lkoff  <- fixed(log(6.55));    label("FcRn-IgG dissociation rate constant (1/h, mouse)")     # inherited from Shah & Betts 2012 text p.73, mouse value

    # Liu 2023 states only that a proportional error variance model was assumed
    # ("Estimation and simulation"); no residual-error magnitude is reported
    # anywhere in the paper or its supplement, so it is carried as fixed(0).
    propSd <- fixed(0); label("Proportional residual error (fraction; magnitude not reported)")
  })

  model({
    clup_scale_min  <- exp(lclup_scale_min)
    clup_scale_max  <- exp(lclup_scale_max)
    rt50_clup_scale <- exp(lrt50_clup_scale)
    hill_clup_scale <- exp(lhill_clup_scale)
    lymph_scale     <- exp(llymph_scale)

    clup_pv <- exp(lclup)
    fcrn_M  <- exp(lfcrn)
    kdeg    <- exp(lkdeg)
    c_lnlf  <- exp(lclnlf)
    kon_M   <- exp(lkon)
    koff    <- exp(lkoff)

    # Molar-to-nanomolar conversions: mAb amounts are nmol and concentrations nM,
    # so FcRn is carried in nM and Kon in 1/(nM h).
    fcrn <- fcrn_M * 1e9   # 4.98e-5 M = 49800 nM
    kon  <- kon_M * 1e-9   # 8.06e7 1/M/h = 0.0806 1/(nM h)

    # Antibody-specific pinocytosis coefficient (Liu 2023 Eq 1 / Eq 2). Sigmoidal
    # in the heparin-chromatography retention time of the administered antibody.
    clup_scale <- clup_scale_min + (clup_scale_max - clup_scale_min) / (1.0 + (rt50_clup_scale / HEPARIN_RT) ^ hill_clup_scale)

    # Fixed system constants (Shah & Betts 2012 text p.73)
    fr    <- 0.715   # fraction of FcRn-bound mAb recycled to the vascular space
    sigis <- 0.2     # lymphatic reflection coefficient, all tissues

    # === Mouse (28 g male) physiology - Shah & Betts 2012 Table 1 ==============
    # Table 1 reports volumes in mL and plasma flows in mL/h; both are divided by
    # 1000 here so the model runs in L, L/h, nmol and nM. Blood-cell
    # sub-compartments are omitted: mAb never enters them, so they are inert
    # (they do not appear in the Liu 2023 Supplementary Methods 2 equations).
    # Lymph flow is 0.2% of the tissue's plasma flow (Supplementary Methods 2,
    # "Lymph Flows": L_ORG = 0.002 * PLQ_ORG).

    v_plasma <- 0.944 / 1000        # Table 1 "Plasma" total volume 0.944 mL
    v_lnode  <- 0.113 / 1000        # Table 1 "Ly. Node" total volume 0.113 mL

    # Heart
    v_ht_vp <- 0.00585 / 1000; v_ht_e <- 0.000760 / 1000; v_ht_is <- 0.0217 / 1000
    q_ht    <- 36.5 / 1000;    sv_ht  <- 0.95
    l_ht    <- 0.002 * q_ht;   clu_ht <- clup_pv * v_ht_e

    # Lung. Volumes and reflection coefficient from Table 1 / text p.73. The lung
    # plasma flow is derived below as the sum of the arterial supplies to the
    # other 14 tissues (see sum_q_X) rather than taken as Table 1's 373 mL/h:
    # 373 mL/h includes the 1.65 mL/h perfusing the lymph node, which is not a
    # mAb-transport path in this model, and using it directly opens a 0.4%
    # arterial leak that is an order of magnitude larger than the antibodies'
    # true clearance (2.6-29 uL/h, Supplementary Table S1).
    v_lu_vp <- 0.0295 / 1000; v_lu_e <- 0.00102 / 1000; v_lu_is <- 0.0384 / 1000
    sv_lu   <- 0.95

    # Muscle
    v_mu_vp <- 0.249 / 1000;  v_mu_e <- 0.0566 / 1000;  v_mu_is <- 1.47 / 1000
    q_mu    <- 86.1 / 1000;   sv_mu  <- 0.95
    l_mu    <- 0.002 * q_mu;  clu_mu <- clup_pv * v_mu_e

    # Skin
    v_sk_vp <- 0.188 / 1000;  v_sk_e <- 0.0251 / 1000;  v_sk_is <- 1.66 / 1000
    q_sk    <- 27.8 / 1000;   sv_sk  <- 0.95
    l_sk    <- 0.002 * q_sk;  clu_sk <- clup_pv * v_sk_e

    # Adipose
    v_ad_vp <- 0.0218 / 1000; v_ad_e <- 0.00991 / 1000; v_ad_is <- 0.337 / 1000
    q_ad    <- 13.4 / 1000;   sv_ad  <- 0.95
    l_ad    <- 0.002 * q_ad;  clu_ad <- clup_pv * v_ad_e

    # Bone
    v_bo_vp <- 0.0621 / 1000; v_bo_e <- 0.0141 / 1000;  v_bo_is <- 0.525 / 1000
    q_bo    <- 15.2 / 1000;   sv_bo  <- 0.85
    l_bo    <- 0.002 * q_bo;  clu_bo <- clup_pv * v_bo_e

    # Brain
    v_br_vp <- 0.0107 / 1000; v_br_e <- 0.00243 / 1000; v_br_is <- 0.0873 / 1000
    q_br    <- 11.8 / 1000;   sv_br  <- 0.99
    l_br    <- 0.002 * q_br;  clu_br <- clup_pv * v_br_e

    # Kidney
    v_ki_vp <- 0.0289 / 1000; v_ki_e <- 0.00263 / 1000; v_ki_is <- 0.0788 / 1000
    q_ki    <- 68.5 / 1000;   sv_ki  <- 0.90
    l_ki    <- 0.002 * q_ki;  clu_ki <- clup_pv * v_ki_e

    # Liver (hepatic-artery supply plus portal return from spleen, pancreas,
    # small intestine and large intestine)
    v_li_vp <- 0.164 / 1000;  v_li_e <- 0.00963 / 1000; v_li_is <- 0.385 / 1000
    q_li    <- 10.3 / 1000;   sv_li  <- 0.85
    clu_li  <- clup_pv * v_li_e

    # Small intestine (drains to liver)
    v_si_vp <- 0.0116 / 1000; v_si_e <- 0.00364 / 1000; v_si_is <- 0.127 / 1000
    q_si    <- 58.1 / 1000;   sv_si  <- 0.90
    l_si    <- 0.002 * q_si;  clu_si <- clup_pv * v_si_e

    # Large intestine (drains to liver)
    v_lr_vp <- 0.0050 / 1000; v_lr_e <- 0.00157 / 1000; v_lr_is <- 0.0545 / 1000
    q_lr    <- 17.3 / 1000;   sv_lr  <- 0.95
    l_lr    <- 0.002 * q_lr;  clu_lr <- clup_pv * v_lr_e

    # Pancreas (drains to liver)
    v_pa_vp <- 0.00534 / 1000; v_pa_e <- 0.000485 / 1000; v_pa_is <- 0.0169 / 1000
    q_pa    <- 6.24 / 1000;    sv_pa  <- 0.90
    l_pa    <- 0.002 * q_pa;   clu_pa <- clup_pv * v_pa_e

    # Thymus
    v_th_vp <- 0.0005 / 1000; v_th_e <- 0.00005 / 1000; v_th_is <- 0.00153 / 1000
    q_th    <- 1.19 / 1000;   sv_th  <- 0.90
    l_th    <- 0.002 * q_th;  clu_th <- clup_pv * v_th_e

    # Spleen (drains to liver)
    v_sp_vp <- 0.0154 / 1000; v_sp_e <- 0.000635 / 1000; v_sp_is <- 0.0254 / 1000
    q_sp    <- 8.18 / 1000;   sv_sp  <- 0.85
    l_sp    <- 0.002 * q_sp;  clu_sp <- clup_pv * v_sp_e

    # Other (composite carcass)
    v_ot_vp <- 0.0195 / 1000; v_ot_e <- 0.00233 / 1000; v_ot_is <- 0.0797 / 1000
    q_ot    <- 10.9 / 1000;   sv_ot  <- 0.95
    l_ot    <- 0.002 * q_ot;  clu_ot <- clup_pv * v_ot_e

    # Liver lymph flow is 0.2% of its TOTAL (arterial + portal) plasma inflow, so
    # that the liver's vascular outflow QVOUT closes at 0.998 of that inflow like
    # every other organ (Supplementary Methods 2, "Lymph Flows").
    q_li_upstream <- (q_si - l_si) + (q_pa - l_pa) + (q_lr - l_lr) + (q_sp - l_sp)
    l_li          <- 0.002 * (q_li + q_li_upstream)

    # Derived lung plasma flow and lymph flow (arterial mass-balance closure).
    sum_q_X <- q_ht + q_mu + q_sk + q_ad + q_bo + q_br + q_ki + q_li +
               q_si + q_lr + q_pa + q_th + q_sp + q_ot
    q_lu    <- sum_q_X
    l_lu    <- 0.002 * q_lu
    clu_lu  <- clup_pv * v_lu_e

    # Lymph-node-to-plasma flow (Supplementary Methods 2: L_LYMPH = PLQ_LUNG * C_LNLF)
    l_lnode <- c_lnlf * q_lu

    # === Concentration aliases =================================================
    cv_p  <- plasma / v_plasma
    cv_ln <- lnode  / v_lnode

    cv_ht <- vp_heart / v_ht_vp; ceu_ht <- eu_heart / v_ht_e
    ceb_ht <- eb_heart / v_ht_e; cfr_ht <- fr_heart / v_ht_e; cis_ht <- is_heart / v_ht_is

    cv_lu <- vp_lung / v_lu_vp; ceu_lu <- eu_lung / v_lu_e
    ceb_lu <- eb_lung / v_lu_e; cfr_lu <- fr_lung / v_lu_e; cis_lu <- is_lung / v_lu_is

    cv_mu <- vp_muscle / v_mu_vp; ceu_mu <- eu_muscle / v_mu_e
    ceb_mu <- eb_muscle / v_mu_e; cfr_mu <- fr_muscle / v_mu_e; cis_mu <- is_muscle / v_mu_is

    cv_sk <- vp_skin / v_sk_vp; ceu_sk <- eu_skin / v_sk_e
    ceb_sk <- eb_skin / v_sk_e; cfr_sk <- fr_skin / v_sk_e; cis_sk <- is_skin / v_sk_is

    cv_ad <- vp_adipose / v_ad_vp; ceu_ad <- eu_adipose / v_ad_e
    ceb_ad <- eb_adipose / v_ad_e; cfr_ad <- fr_adipose / v_ad_e; cis_ad <- is_adipose / v_ad_is

    cv_bo <- vp_bone / v_bo_vp; ceu_bo <- eu_bone / v_bo_e
    ceb_bo <- eb_bone / v_bo_e; cfr_bo <- fr_bone / v_bo_e; cis_bo <- is_bone / v_bo_is

    cv_br <- vp_brain / v_br_vp; ceu_br <- eu_brain / v_br_e
    ceb_br <- eb_brain / v_br_e; cfr_br <- fr_brain / v_br_e; cis_br <- is_brain / v_br_is

    cv_ki <- vp_kidney / v_ki_vp; ceu_ki <- eu_kidney / v_ki_e
    ceb_ki <- eb_kidney / v_ki_e; cfr_ki <- fr_kidney / v_ki_e; cis_ki <- is_kidney / v_ki_is

    cv_li <- vp_liver / v_li_vp; ceu_li <- eu_liver / v_li_e
    ceb_li <- eb_liver / v_li_e; cfr_li <- fr_liver / v_li_e; cis_li <- is_liver / v_li_is

    cv_si <- vp_small_intestine / v_si_vp; ceu_si <- eu_small_intestine / v_si_e
    ceb_si <- eb_small_intestine / v_si_e; cfr_si <- fr_small_intestine / v_si_e
    cis_si <- is_small_intestine / v_si_is

    cv_lr <- vp_large_intestine / v_lr_vp; ceu_lr <- eu_large_intestine / v_lr_e
    ceb_lr <- eb_large_intestine / v_lr_e; cfr_lr <- fr_large_intestine / v_lr_e
    cis_lr <- is_large_intestine / v_lr_is

    cv_pa <- vp_pancreas / v_pa_vp; ceu_pa <- eu_pancreas / v_pa_e
    ceb_pa <- eb_pancreas / v_pa_e; cfr_pa <- fr_pancreas / v_pa_e; cis_pa <- is_pancreas / v_pa_is

    cv_th <- vp_thymus / v_th_vp; ceu_th <- eu_thymus / v_th_e
    ceb_th <- eb_thymus / v_th_e; cfr_th <- fr_thymus / v_th_e; cis_th <- is_thymus / v_th_is

    cv_sp <- vp_spleen / v_sp_vp; ceu_sp <- eu_spleen / v_sp_e
    ceb_sp <- eb_spleen / v_sp_e; cfr_sp <- fr_spleen / v_sp_e; cis_sp <- is_spleen / v_sp_is

    cv_ot <- vp_other / v_ot_vp; ceu_ot <- eu_other / v_ot_e
    ceb_ot <- eb_other / v_ot_e; cfr_ot <- fr_other / v_ot_e; cis_ot <- is_other / v_ot_is

    # === Organ vascular spaces ================================================
    # dA_V/dt = M_VI - QVOUT * C_V - (1 - sigV) * F2 * L * C_V
    #           - F1 * CLup_ORG * C_V + FR * CLup_ORG * C_E,B
    # (Liu 2023 Supplementary Methods 2, "Organ Vascular space"). The arterial
    # supply of every organ except the lung and the liver is the lung vascular
    # concentration; the lung is supplied from central plasma.

    d/dt(vp_heart) <- q_ht * cv_lu - (q_ht - l_ht) * cv_ht - (1 - sv_ht) * lymph_scale * l_ht * cv_ht -
      clup_scale * clu_ht * cv_ht + fr * clu_ht * ceb_ht
    d/dt(vp_lung) <- (q_lu + l_lu) * cv_p - q_lu * cv_lu - (1 - sv_lu) * lymph_scale * l_lu * cv_lu -
      clup_scale * clu_lu * cv_lu + fr * clu_lu * ceb_lu
    d/dt(vp_muscle) <- q_mu * cv_lu - (q_mu - l_mu) * cv_mu - (1 - sv_mu) * lymph_scale * l_mu * cv_mu -
      clup_scale * clu_mu * cv_mu + fr * clu_mu * ceb_mu
    d/dt(vp_skin) <- q_sk * cv_lu - (q_sk - l_sk) * cv_sk - (1 - sv_sk) * lymph_scale * l_sk * cv_sk -
      clup_scale * clu_sk * cv_sk + fr * clu_sk * ceb_sk
    d/dt(vp_adipose) <- q_ad * cv_lu - (q_ad - l_ad) * cv_ad - (1 - sv_ad) * lymph_scale * l_ad * cv_ad -
      clup_scale * clu_ad * cv_ad + fr * clu_ad * ceb_ad
    d/dt(vp_bone) <- q_bo * cv_lu - (q_bo - l_bo) * cv_bo - (1 - sv_bo) * lymph_scale * l_bo * cv_bo -
      clup_scale * clu_bo * cv_bo + fr * clu_bo * ceb_bo
    d/dt(vp_brain) <- q_br * cv_lu - (q_br - l_br) * cv_br - (1 - sv_br) * lymph_scale * l_br * cv_br -
      clup_scale * clu_br * cv_br + fr * clu_br * ceb_br
    d/dt(vp_kidney) <- q_ki * cv_lu - (q_ki - l_ki) * cv_ki - (1 - sv_ki) * lymph_scale * l_ki * cv_ki -
      clup_scale * clu_ki * cv_ki + fr * clu_ki * ceb_ki
    d/dt(vp_small_intestine) <- q_si * cv_lu - (q_si - l_si) * cv_si -
      (1 - sv_si) * lymph_scale * l_si * cv_si - clup_scale * clu_si * cv_si + fr * clu_si * ceb_si
    d/dt(vp_large_intestine) <- q_lr * cv_lu - (q_lr - l_lr) * cv_lr -
      (1 - sv_lr) * lymph_scale * l_lr * cv_lr - clup_scale * clu_lr * cv_lr + fr * clu_lr * ceb_lr
    d/dt(vp_pancreas) <- q_pa * cv_lu - (q_pa - l_pa) * cv_pa - (1 - sv_pa) * lymph_scale * l_pa * cv_pa -
      clup_scale * clu_pa * cv_pa + fr * clu_pa * ceb_pa
    d/dt(vp_thymus) <- q_th * cv_lu - (q_th - l_th) * cv_th - (1 - sv_th) * lymph_scale * l_th * cv_th -
      clup_scale * clu_th * cv_th + fr * clu_th * ceb_th
    d/dt(vp_spleen) <- q_sp * cv_lu - (q_sp - l_sp) * cv_sp - (1 - sv_sp) * lymph_scale * l_sp * cv_sp -
      clup_scale * clu_sp * cv_sp + fr * clu_sp * ceb_sp
    d/dt(vp_other) <- q_ot * cv_lu - (q_ot - l_ot) * cv_ot - (1 - sv_ot) * lymph_scale * l_ot * cv_ot -
      clup_scale * clu_ot * cv_ot + fr * clu_ot * ceb_ot

    # Liver: hepatic-artery inflow plus the four portal tributaries; the total
    # outflow QVOUT_LIVER = q_li - l_li + q_li_upstream returns to central plasma.
    qvout_li <- q_li - l_li + q_li_upstream
    d/dt(vp_liver) <- q_li * cv_lu +
      (q_si - l_si) * cv_si + (q_pa - l_pa) * cv_pa +
      (q_lr - l_lr) * cv_lr + (q_sp - l_sp) * cv_sp -
      qvout_li * cv_li - (1 - sv_li) * lymph_scale * l_li * cv_li -
      clup_scale * clu_li * cv_li + fr * clu_li * ceb_li

    # === Endosomal spaces ======================================================
    # Uptake into the endosome is scaled by F1; FcRn-complex efflux (recycling to
    # vascular and interstitial space) is not (Liu 2023 Supplementary Methods 2,
    # "Endosomal Space"). Free-plus-bound FcRn is conserved within each organ.

    bind_ht <- kon * ceu_ht * cfr_ht * v_ht_e; ub_ht <- koff * eb_heart
    d/dt(eu_heart) <- clup_scale * clu_ht * (cv_ht + cis_ht) - bind_ht + ub_ht - kdeg * eu_heart
    d/dt(eb_heart) <- bind_ht - ub_ht - clup_pv * eb_heart
    d/dt(fr_heart) <- ub_ht - bind_ht + clup_pv * eb_heart

    bind_lu <- kon * ceu_lu * cfr_lu * v_lu_e; ub_lu <- koff * eb_lung
    d/dt(eu_lung) <- clup_scale * clu_lu * (cv_lu + cis_lu) - bind_lu + ub_lu - kdeg * eu_lung
    d/dt(eb_lung) <- bind_lu - ub_lu - clup_pv * eb_lung
    d/dt(fr_lung) <- ub_lu - bind_lu + clup_pv * eb_lung

    bind_mu <- kon * ceu_mu * cfr_mu * v_mu_e; ub_mu <- koff * eb_muscle
    d/dt(eu_muscle) <- clup_scale * clu_mu * (cv_mu + cis_mu) - bind_mu + ub_mu - kdeg * eu_muscle
    d/dt(eb_muscle) <- bind_mu - ub_mu - clup_pv * eb_muscle
    d/dt(fr_muscle) <- ub_mu - bind_mu + clup_pv * eb_muscle

    bind_sk <- kon * ceu_sk * cfr_sk * v_sk_e; ub_sk <- koff * eb_skin
    d/dt(eu_skin) <- clup_scale * clu_sk * (cv_sk + cis_sk) - bind_sk + ub_sk - kdeg * eu_skin
    d/dt(eb_skin) <- bind_sk - ub_sk - clup_pv * eb_skin
    d/dt(fr_skin) <- ub_sk - bind_sk + clup_pv * eb_skin

    bind_ad <- kon * ceu_ad * cfr_ad * v_ad_e; ub_ad <- koff * eb_adipose
    d/dt(eu_adipose) <- clup_scale * clu_ad * (cv_ad + cis_ad) - bind_ad + ub_ad - kdeg * eu_adipose
    d/dt(eb_adipose) <- bind_ad - ub_ad - clup_pv * eb_adipose
    d/dt(fr_adipose) <- ub_ad - bind_ad + clup_pv * eb_adipose

    bind_bo <- kon * ceu_bo * cfr_bo * v_bo_e; ub_bo <- koff * eb_bone
    d/dt(eu_bone) <- clup_scale * clu_bo * (cv_bo + cis_bo) - bind_bo + ub_bo - kdeg * eu_bone
    d/dt(eb_bone) <- bind_bo - ub_bo - clup_pv * eb_bone
    d/dt(fr_bone) <- ub_bo - bind_bo + clup_pv * eb_bone

    bind_br <- kon * ceu_br * cfr_br * v_br_e; ub_br <- koff * eb_brain
    d/dt(eu_brain) <- clup_scale * clu_br * (cv_br + cis_br) - bind_br + ub_br - kdeg * eu_brain
    d/dt(eb_brain) <- bind_br - ub_br - clup_pv * eb_brain
    d/dt(fr_brain) <- ub_br - bind_br + clup_pv * eb_brain

    bind_ki <- kon * ceu_ki * cfr_ki * v_ki_e; ub_ki <- koff * eb_kidney
    d/dt(eu_kidney) <- clup_scale * clu_ki * (cv_ki + cis_ki) - bind_ki + ub_ki - kdeg * eu_kidney
    d/dt(eb_kidney) <- bind_ki - ub_ki - clup_pv * eb_kidney
    d/dt(fr_kidney) <- ub_ki - bind_ki + clup_pv * eb_kidney

    bind_li <- kon * ceu_li * cfr_li * v_li_e; ub_li <- koff * eb_liver
    d/dt(eu_liver) <- clup_scale * clu_li * (cv_li + cis_li) - bind_li + ub_li - kdeg * eu_liver
    d/dt(eb_liver) <- bind_li - ub_li - clup_pv * eb_liver
    d/dt(fr_liver) <- ub_li - bind_li + clup_pv * eb_liver

    bind_si <- kon * ceu_si * cfr_si * v_si_e; ub_si <- koff * eb_small_intestine
    d/dt(eu_small_intestine) <- clup_scale * clu_si * (cv_si + cis_si) - bind_si + ub_si -
      kdeg * eu_small_intestine
    d/dt(eb_small_intestine) <- bind_si - ub_si - clup_pv * eb_small_intestine
    d/dt(fr_small_intestine) <- ub_si - bind_si + clup_pv * eb_small_intestine

    bind_lr <- kon * ceu_lr * cfr_lr * v_lr_e; ub_lr <- koff * eb_large_intestine
    d/dt(eu_large_intestine) <- clup_scale * clu_lr * (cv_lr + cis_lr) - bind_lr + ub_lr -
      kdeg * eu_large_intestine
    d/dt(eb_large_intestine) <- bind_lr - ub_lr - clup_pv * eb_large_intestine
    d/dt(fr_large_intestine) <- ub_lr - bind_lr + clup_pv * eb_large_intestine

    bind_pa <- kon * ceu_pa * cfr_pa * v_pa_e; ub_pa <- koff * eb_pancreas
    d/dt(eu_pancreas) <- clup_scale * clu_pa * (cv_pa + cis_pa) - bind_pa + ub_pa - kdeg * eu_pancreas
    d/dt(eb_pancreas) <- bind_pa - ub_pa - clup_pv * eb_pancreas
    d/dt(fr_pancreas) <- ub_pa - bind_pa + clup_pv * eb_pancreas

    bind_th <- kon * ceu_th * cfr_th * v_th_e; ub_th <- koff * eb_thymus
    d/dt(eu_thymus) <- clup_scale * clu_th * (cv_th + cis_th) - bind_th + ub_th - kdeg * eu_thymus
    d/dt(eb_thymus) <- bind_th - ub_th - clup_pv * eb_thymus
    d/dt(fr_thymus) <- ub_th - bind_th + clup_pv * eb_thymus

    bind_sp <- kon * ceu_sp * cfr_sp * v_sp_e; ub_sp <- koff * eb_spleen
    d/dt(eu_spleen) <- clup_scale * clu_sp * (cv_sp + cis_sp) - bind_sp + ub_sp - kdeg * eu_spleen
    d/dt(eb_spleen) <- bind_sp - ub_sp - clup_pv * eb_spleen
    d/dt(fr_spleen) <- ub_sp - bind_sp + clup_pv * eb_spleen

    bind_ot <- kon * ceu_ot * cfr_ot * v_ot_e; ub_ot <- koff * eb_other
    d/dt(eu_other) <- clup_scale * clu_ot * (cv_ot + cis_ot) - bind_ot + ub_ot - kdeg * eu_other
    d/dt(eb_other) <- bind_ot - ub_ot - clup_pv * eb_other
    d/dt(fr_other) <- ub_ot - bind_ot + clup_pv * eb_other

    # === Interstitial spaces ===================================================
    # Convective entry and lymphatic exit both scale with F2; endosomal uptake
    # from the interstitium scales with F1, and the (1 - FR) fraction of recycled
    # FcRn-complex is delivered here (Supplementary Methods 2, "Interstitial Space").

    d/dt(is_heart) <- (1 - sv_ht) * lymph_scale * l_ht * cv_ht - (1 - sigis) * lymph_scale * l_ht * cis_ht -
      clup_scale * clu_ht * cis_ht + clu_ht * (1 - fr) * ceb_ht
    d/dt(is_lung) <- (1 - sv_lu) * lymph_scale * l_lu * cv_lu - (1 - sigis) * lymph_scale * l_lu * cis_lu -
      clup_scale * clu_lu * cis_lu + clu_lu * (1 - fr) * ceb_lu
    d/dt(is_muscle) <- (1 - sv_mu) * lymph_scale * l_mu * cv_mu - (1 - sigis) * lymph_scale * l_mu * cis_mu -
      clup_scale * clu_mu * cis_mu + clu_mu * (1 - fr) * ceb_mu
    d/dt(is_skin) <- (1 - sv_sk) * lymph_scale * l_sk * cv_sk - (1 - sigis) * lymph_scale * l_sk * cis_sk -
      clup_scale * clu_sk * cis_sk + clu_sk * (1 - fr) * ceb_sk
    d/dt(is_adipose) <- (1 - sv_ad) * lymph_scale * l_ad * cv_ad - (1 - sigis) * lymph_scale * l_ad * cis_ad -
      clup_scale * clu_ad * cis_ad + clu_ad * (1 - fr) * ceb_ad
    d/dt(is_bone) <- (1 - sv_bo) * lymph_scale * l_bo * cv_bo - (1 - sigis) * lymph_scale * l_bo * cis_bo -
      clup_scale * clu_bo * cis_bo + clu_bo * (1 - fr) * ceb_bo
    d/dt(is_brain) <- (1 - sv_br) * lymph_scale * l_br * cv_br - (1 - sigis) * lymph_scale * l_br * cis_br -
      clup_scale * clu_br * cis_br + clu_br * (1 - fr) * ceb_br
    d/dt(is_kidney) <- (1 - sv_ki) * lymph_scale * l_ki * cv_ki - (1 - sigis) * lymph_scale * l_ki * cis_ki -
      clup_scale * clu_ki * cis_ki + clu_ki * (1 - fr) * ceb_ki
    d/dt(is_liver) <- (1 - sv_li) * lymph_scale * l_li * cv_li - (1 - sigis) * lymph_scale * l_li * cis_li -
      clup_scale * clu_li * cis_li + clu_li * (1 - fr) * ceb_li
    d/dt(is_small_intestine) <- (1 - sv_si) * lymph_scale * l_si * cv_si -
      (1 - sigis) * lymph_scale * l_si * cis_si - clup_scale * clu_si * cis_si + clu_si * (1 - fr) * ceb_si
    d/dt(is_large_intestine) <- (1 - sv_lr) * lymph_scale * l_lr * cv_lr -
      (1 - sigis) * lymph_scale * l_lr * cis_lr - clup_scale * clu_lr * cis_lr + clu_lr * (1 - fr) * ceb_lr
    d/dt(is_pancreas) <- (1 - sv_pa) * lymph_scale * l_pa * cv_pa - (1 - sigis) * lymph_scale * l_pa * cis_pa -
      clup_scale * clu_pa * cis_pa + clu_pa * (1 - fr) * ceb_pa
    d/dt(is_thymus) <- (1 - sv_th) * lymph_scale * l_th * cv_th - (1 - sigis) * lymph_scale * l_th * cis_th -
      clup_scale * clu_th * cis_th + clu_th * (1 - fr) * ceb_th
    d/dt(is_spleen) <- (1 - sv_sp) * lymph_scale * l_sp * cv_sp - (1 - sigis) * lymph_scale * l_sp * cis_sp -
      clup_scale * clu_sp * cis_sp + clu_sp * (1 - fr) * ceb_sp
    d/dt(is_other) <- (1 - sv_ot) * lymph_scale * l_ot * cv_ot - (1 - sigis) * lymph_scale * l_ot * cis_ot -
      clup_scale * clu_ot * cis_ot + clu_ot * (1 - fr) * ceb_ot

    # === Central plasma (Supplementary Methods 2, "Plasma") ====================
    d/dt(plasma) <- (q_ht - l_ht) * cv_ht +
      (q_mu - l_mu) * cv_mu +
      (q_sk - l_sk) * cv_sk +
      (q_ad - l_ad) * cv_ad +
      (q_bo - l_bo) * cv_bo +
      (q_br - l_br) * cv_br +
      (q_ki - l_ki) * cv_ki +
      qvout_li * cv_li +
      (q_th - l_th) * cv_th +
      (q_ot - l_ot) * cv_ot +
      l_lnode * cv_ln -
      (q_lu + l_lu) * cv_p

    # === Lymph node (Supplementary Methods 2, "Lymph Nodes") ==================
    d/dt(lnode) <- (1 - sigis) * lymph_scale * (
      l_ht * cis_ht + l_lu * cis_lu + l_mu * cis_mu + l_sk * cis_sk + l_ad * cis_ad +
      l_bo * cis_bo + l_br * cis_br + l_ki * cis_ki + l_li * cis_li + l_si * cis_si +
      l_lr * cis_lr + l_pa * cis_pa + l_th * cis_th + l_sp * cis_sp + l_ot * cis_ot
    ) - l_lnode * cv_ln

    # === Initial conditions ====================================================
    # All mAb states start empty (the IV bolus enters `plasma` from the event
    # table). Every organ endosome starts with its full FcRn complement free.
    fr_heart(0) <- fcrn * v_ht_e
    fr_lung(0) <- fcrn * v_lu_e
    fr_muscle(0) <- fcrn * v_mu_e
    fr_skin(0) <- fcrn * v_sk_e
    fr_adipose(0) <- fcrn * v_ad_e
    fr_bone(0) <- fcrn * v_bo_e
    fr_brain(0) <- fcrn * v_br_e
    fr_kidney(0) <- fcrn * v_ki_e
    fr_liver(0) <- fcrn * v_li_e
    fr_small_intestine(0) <- fcrn * v_si_e
    fr_large_intestine(0) <- fcrn * v_lr_e
    fr_pancreas(0) <- fcrn * v_pa_e
    fr_thymus(0) <- fcrn * v_th_e
    fr_spleen(0) <- fcrn * v_sp_e
    fr_other(0) <- fcrn * v_ot_e

    # === Observation ===========================================================
    # Cc is the central plasma concentration in nM. Users dosing in mass units
    # convert with dose_nmol = dose_mg / MW * 1e9 (MW = antibody molar mass).
    Cc <- cv_p
    Cc ~ prop(propSd)
  })
}
