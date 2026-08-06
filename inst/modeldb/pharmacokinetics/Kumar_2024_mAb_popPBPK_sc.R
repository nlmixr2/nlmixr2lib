Kumar_2024_mAb_popPBPK_sc <- function() {
  description <- paste(
    "popPBPK (whole-body). Kumar 2024 population physiologically based",
    "pharmacokinetic model characterizing INTER-ANTIBODY variability in the",
    "clinical PK of monoclonal antibodies with linear disposition after",
    "SUBCUTANEOUS administration. This is the intravenous model",
    "(Kumar_2024_mAb_popPBPK_iv.R) with the skin bifurcated into a small",
    "injection-site 'subcutaneous' compartment (2.25 mL interstitial space,",
    "sized to a typical ~2 mL SC injection) and a 'rest of skin' compartment,",
    "giving 16 tissues x 6 sub-compartments plus central plasma, central",
    "blood cells and a lymph node - 99 ODE states. The SC dose is",
    "administered into the interstitial space of the subcutaneous",
    "compartment (cmt = 'is_subcutaneous'), NOT into an empirical first-order",
    "depot. Two absorption parameters are added there: a first-order local",
    "degradation rate kSC representing pre-systemic catabolism at the",
    "injection site, and a scaling factor S_LU multiplying the convective",
    "lymphatic flow that carries mAb from the SC interstitium to the lymph",
    "node. Both were estimated in Monolix 2021R1 by fitting to the SC",
    "profiles of the 16 mAbs that also had IV data, with the antibody-specific",
    "CLup and kdeg values fixed as regressors from the IV fit. All four",
    "random effects (CLup, kdeg, kSC, S_LU) are INTER-ANTIBODY variability,",
    "not between-subject variability: sampling them draws a new hypothetical",
    "antibody, not a new patient. Physiology is the Shah and Betts 2012 human",
    "(71 kg male) parameter set, with the skin and subcutaneous rows taken",
    "from Kumar 2024 Table 1."
  )
  reference <- paste(
    "Kumar M, Lanke S, Yadav A, Ette M, Mager DE, Shah DK (2024).",
    "Inter-Antibody Variability in the Clinical Pharmacokinetics of",
    "Monoclonal Antibodies Characterized Using Population Physiologically",
    "Based Pharmacokinetic Modeling. Antibodies 13(3):54.",
    "doi:10.3390/antib13030054.",
    "Structure, equations and human physiology inherited from",
    "Shah DK, Betts AM (2012). Towards a platform PBPK model to characterize",
    "the plasma and tissue disposition of monoclonal antibodies in",
    "preclinical species and human. J Pharmacokinet Pharmacodyn 39(1):67-86.",
    "doi:10.1007/s10928-011-9232-2.",
    sep = " "
  )
  vignette <- "Kumar_2024_mAb_popPBPK"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    vp_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "administration site", verified = FALSE),
    vp_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    vp_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    vp_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_thymus          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_other           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    vp_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "administration site", verified = FALSE),
    bc_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    bc_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    bc_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_thymus          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_other           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    bc_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_thymus          = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_thymus          = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_thymus          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eb_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    fr_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    eu_other           = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE),
    eb_other           = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    fr_other           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_heart           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_lung            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_muscle          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_skin            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_subcutaneous    = list(analyte = "mAb", units = "nmol", specimen = "administration site", verified = FALSE),
    is_adipose         = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_bone            = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_brain           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_kidney          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_liver           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_small_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    is_large_intestine = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    is_pancreas        = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_thymus          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_spleen          = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    is_other           = list(analyte = "mAb", units = "nmol", specimen = "tissue", verified = FALSE),
    plasma             = list(analyte = "mAb", units = "nmol", specimen = "plasma", verified = FALSE),
    bcc                = list(analyte = "mAb", units = "nmol", specimen = "blood cell", verified = FALSE),
    lnode              = list(analyte = "mAb", units = "nmol", specimen = "lymph", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    FDA_APPROVED = list(
      description = paste(
        "FDA approval status of the antibody (1 = FDA approved,",
        "0 = clinically tested but not approved)."
      ),
      units = "binary",
      type = "categorical",
      notes = paste(
        "Kumar 2024 Figure 6 compares the distributions of the individual",
        "(empirical Bayes) CLup and kdeg estimates between FDA-approved and",
        "clinically tested mAbs and reports medians of 0.39 vs 0.28 L/h/L",
        "(CLup) and 25.5 vs 29.3 1/h (kdeg). The paper concludes that 'no",
        "clear distinction was inferred in the disposition of FDA-approved",
        "and clinically tested mAbs', and approval status was NEVER entered",
        "into the model as a covariate - no coefficient is estimated or",
        "reported. Documented here for provenance only; it is deliberately",
        "not referenced in model()."
      )
    )
  )

  population <- list(
    n_subjects     = 16,
    n_studies      = 16,
    age_range      = "adults (predominantly healthy-volunteer phase 1 cohorts)",
    weight_range   = "71 kg reference adult male (Shah and Betts 2012 Table 4 physiology)",
    sex_female_pct = NA_real_,
    species        = "human",
    disease_state  = paste(
      "Mixed. Predominantly phase 1 single-dose studies in healthy",
      "volunteers, with some patient cohorts (asthma, psoriasis, rheumatoid",
      "arthritis and others). See Kumar 2024 Supplementary Table S3 for the",
      "per-antibody study list."
    ),
    dose_range     = "multiple SC dose levels per antibody (digitized phase 1 dose-escalation profiles)",
    regions        = "multi-national (USA, Canada, Europe, Japan, Korea, China and others; Supplementary Table S3)",
    notes          = paste(
      "The 'individuals' of this population model are ANTIBODIES, not human",
      "subjects. The SC model was built from the 16 mAbs that had BOTH IV and",
      "SC data (Table 4: adalimumab, canakinumab, guselkumab, risankizumab,",
      "secukinumab, tezepelumab, tildrakizumab, mepolizumab, daclizumab,",
      "benralizumab, belimumab, enokizumab, gevokizumab, fulranumab,",
      "olokizumab, ralpancizumab; all median %PE below 30%). A further nine",
      "mAbs with SC-only data (emicizumab, etrolizumab, fremanezumab,",
      "galcanezumab, ixekizumab, lanadelumab, omalizumab, quilizumab,",
      "tralokinumab) were held out and used only for the a priori",
      "prediction-window validation in Figure 8."
    ),
    scope_note     = paste(
      "Simulating with the etas active generates a range of plausible",
      "ANTIBODY-level typical profiles (the paper's Monte Carlo prediction",
      "window, Figure 8), NOT between-patient variability within one",
      "antibody. There is no within-antibody residual error model in this",
      "file; Kumar 2024 used a combined error model during fitting but does",
      "not report its coefficients (see vignette Errata)."
    )
  )

  ini({
    # --- Re-estimated by Kumar 2024 (Supplementary Table S2; RSE% in comments) ---
    lclup <- log(0.32)
    label("Pinocytosis rate per unit endosomal volume CLup (L/h/L)") # Suppl Table S2, 0.32 (RSE 5.6%)
    lkdeg <- log(26.1)
    label("Endosomal degradation rate of FcRn-unbound mAb kdeg (1/h)") # Suppl Table S2, 26.1 (RSE 11%)
    lksc <- log(0.0015)
    label("First-order local degradation rate at the SC injection site kSC (1/h)") # Suppl Table S2, 0.0015 (RSE 66%)
    lslu <- log(0.54)
    label("Scaling factor on lymphatic uptake from SC interstitium S_LU (unitless)") # Suppl Table S2, 0.54 (RSE 14%)

    # --- Inherited and FIXED from Shah and Betts 2012 ---
    # Kumar 2024 Section 2.5: "All model parameters apart from pinocytosis
    # rate (CLup) and non-specific degradation rate of unbound antibody
    # (kdeg) are the same as those used for humans by Shah and Betts".
    lfcrn <- fixed(log(4.98e-5))
    label("Endosomal free-FcRn concentration (mol/L)") # Shah and Betts 2012 Table 6
    lkon <- fixed(log(5.59e8))
    label("FcRn-IgG association rate constant (1/M/h, human)") # Shah and Betts 2012 text p.73
    lkoff <- fixed(log(23.9))
    label("FcRn-IgG dissociation rate constant (1/h, human)") # Shah and Betts 2012 text p.73
    llymphnode <- fixed(log(3.670))
    label("Lymph-node-to-plasma flow (L/h, human)") # Shah and Betts 2012 Table 4 row Ly. Node

    # --- Inter-antibody variability (Supplementary Table S2) ---
    # omega values are the standard deviations of the log-normal random
    # effects; nlmixr2 takes the VARIANCE, so each is squared here.
    etalclup ~ 0.5329 # omega_CLup = 73% (RSE 3.2%); 0.73^2
    etalkdeg ~ 0.2116 # omega_kdeg = 46% (RSE 17.6%); 0.46^2
    etalksc ~ 3.7249 # omega_kSC = 193% (RSE 24.7%); 1.93^2
    etalslu ~ 0.2401 # omega_S_LU = 49% (RSE 21.1%); 0.49^2
  })

  model({
    clup_pv <- exp(lclup + etalclup)
    kdeg <- exp(lkdeg + etalkdeg)
    ksc <- exp(lksc + etalksc)
    s_lu <- exp(lslu + etalslu)
    fcrn_M <- exp(lfcrn)
    kon_M <- exp(lkon)
    koff <- exp(lkoff)
    l_lnode <- exp(llymphnode)

    # Convert FcRn molar concentration to nM and kon from 1/M/h to 1/(nM h)
    # so that all states can be carried as nmol.
    fcrn <- fcrn_M * 1e9
    kon <- kon_M * 1e-9

    # Fixed system constants (Shah and Betts 2012 text p.73)
    fr <- 0.715 # Fraction of FcRn-bound mAb recycled to the vascular space
    sigis <- 0.2 # Lymphatic (interstitial) reflection coefficient, all tissues

    # === Human (71 kg male) physiology - Shah and Betts 2012 Table 4 ===
    # Volumes mL -> L; flows mL/h -> L/h. Lymph flow is plasma flow / 500
    # (Shah and Betts 2012 "Model parameters"). Vascular reflection
    # coefficients sv_* are the a priori pore-size values from the same
    # section (0.95 / 0.90 / 0.85, brain 0.99).

    # Central pools. Kumar 2024 Section 2.5: central plasma and blood cell
    # volumes are total blood volume MINUS the vascular volume held in the
    # tissues, giving 1412 and 1155 mL. This is the single physiological
    # difference from Shah_2012_mAb_PBPK.R, which used the Table 4 whole-body
    # plasma (3126 mL) and blood-cell (2558 mL) volumes directly.
    v_plasma <- 1.412 # Kumar 2024 Section 2.5
    v_bcc <- 1.155 # Kumar 2024 Section 2.5

    # Lymph node - Table 4 total volume 274 mL, single mixing compartment.
    v_lnode <- 0.274

    # Heart
    v_ht_vp <- 0.0131; v_ht_bc <- 0.0108; v_ht_e <- 0.00171; v_ht_is <- 0.0488
    q_ht <- 7.752; bcq_ht <- 6.342; sv_ht <- 0.95
    l_ht <- q_ht / 500
    clu_ht <- clup_pv * v_ht_e

    # Lung. Plasma and blood-cell flows are derived below from the sum of the
    # non-lung tissue inflows so the lung vascular mass balance closes
    # exactly; the Table 4 values (181913 and 148838 mL/h) are ~1.8% higher
    # than that sum and would leak mass at the lung-arterial junction.
    v_lu_vp <- 0.0550; v_lu_bc <- 0.0450; v_lu_e <- 0.00500; v_lu_is <- 0.300
    sv_lu <- 0.95

    # Muscle
    v_mu_vp <- 0.662; v_mu_bc <- 0.541; v_mu_e <- 0.150; v_mu_is <- 3.910
    q_mu <- 33.469; bcq_mu <- 27.383; sv_mu <- 0.95
    l_mu <- q_mu / 500
    clu_mu <- clup_pv * v_mu_e

    # Skin ("rest of skin" after the SC injection site is carved out).
    # Kumar 2024 Table 1 column "Skin (SC)". These are the Shah and Betts
    # Table 4 skin values minus the subcutaneous compartment below, e.g.
    # 3408 - 6.82 = 3401 mL total and 1125 - 2.25 = 1123 mL interstitial.
    v_sk_vp <- 0.1272; v_sk_bc <- 0.1041; v_sk_e <- 0.01701; v_sk_is <- 1.123
    q_sk <- 11.600; bcq_sk <- 9.493; sv_sk <- 0.95
    l_sk <- q_sk / 500
    clu_sk <- clup_pv * v_sk_e

    # Subcutaneous injection site. Kumar 2024 Table 1 column "SC (SC)":
    # the interstitial sub-compartment was set to 2.25 mL (a typical ~2 mL
    # SC injection) and the other sub-volumes and the flows were scaled from
    # the skin compartment keeping the skin sub-compartment ratios.
    # The vascular reflection coefficient is the skin value (Section 2.3:
    # the SC tissue "has a similar structure to other tissues").
    v_sc_vp <- 0.00025; v_sc_bc <- 0.00021; v_sc_e <- 0.00003; v_sc_is <- 0.00225
    q_sc <- 0.02325; bcq_sc <- 0.01902; sv_sc <- 0.95
    l_sc <- q_sc / 500
    clu_sc <- clup_pv * v_sc_e

    # Adipose
    v_ad_vp <- 0.148; v_ad_bc <- 0.121; v_ad_e <- 0.0673; v_ad_is <- 2.289
    q_ad <- 11.233; bcq_ad <- 9.191; sv_ad <- 0.95
    l_ad <- q_ad / 500
    clu_ad <- clup_pv * v_ad_e

    # Bone
    v_bo_vp <- 0.224; v_bo_bc <- 0.183; v_bo_e <- 0.0508; v_bo_is <- 1.891
    q_bo <- 2.591; bcq_bo <- 2.120; sv_bo <- 0.85
    l_bo <- q_bo / 500
    clu_bo <- clup_pv * v_bo_e

    # Brain. Present in Kumar 2024 Figure 1 (and required by the central
    # plasma volume arithmetic in Section 2.5) even though the Section 2.2
    # prose list of tissues omits it - see vignette Errata.
    v_br_vp <- 0.0319; v_br_bc <- 0.0261; v_br_e <- 0.00725; v_br_is <- 0.261
    q_br <- 21.453; bcq_br <- 17.553; sv_br <- 0.99
    l_br <- q_br / 500
    clu_br <- clup_pv * v_br_e

    # Kidney
    v_ki_vp <- 0.0182; v_ki_bc <- 0.0149; v_ki_e <- 0.00166; v_ki_is <- 0.0498
    q_ki <- 36.402; bcq_ki <- 29.784; sv_ki <- 0.90
    l_ki <- q_ki / 500
    clu_ki <- clup_pv * v_ki_e

    # Liver (receives splanchnic portal return - Shah and Betts 2012 Eq 11)
    v_li_vp <- 0.183; v_li_bc <- 0.149; v_li_e <- 0.0107; v_li_is <- 0.429
    q_li <- 13.210; bcq_li <- 10.808; sv_li <- 0.85
    l_li <- q_li / 500
    clu_li <- clup_pv * v_li_e

    # Small intestine (drains to liver)
    v_si_vp <- 0.00615; v_si_bc <- 0.00503; v_si_e <- 0.00193; v_si_is <- 0.0671
    q_si <- 12.368; bcq_si <- 10.120; sv_si <- 0.90
    l_si <- q_si / 500
    clu_si <- clup_pv * v_si_e

    # Large intestine (drains to liver)
    v_lr_vp <- 0.00874; v_lr_bc <- 0.00715; v_lr_e <- 0.00274; v_lr_is <- 0.0953
    q_lr <- 12.867; bcq_lr <- 10.527; sv_lr <- 0.95
    l_lr <- q_lr / 500
    clu_lr <- clup_pv * v_lr_e

    # Pancreas (drains to liver)
    v_pa_vp <- 0.00570; v_pa_bc <- 0.00466; v_pa_e <- 0.000518; v_pa_is <- 0.0180
    q_pa <- 3.056; bcq_pa <- 2.500; sv_pa <- 0.90
    l_pa <- q_pa / 500
    clu_pa <- clup_pv * v_pa_e

    # Thymus
    v_th_vp <- 0.000353; v_th_bc <- 0.000288; v_th_e <- 0.0000321; v_th_is <- 0.00109
    q_th <- 0.353; bcq_th <- 0.289; sv_th <- 0.90
    l_th <- q_th / 500
    clu_th <- clup_pv * v_th_e

    # Spleen (drains to liver)
    v_sp_vp <- 0.0268; v_sp_bc <- 0.0219; v_sp_e <- 0.00111; v_sp_is <- 0.0443
    q_sp <- 6.343; bcq_sp <- 5.189; sv_sp <- 0.85
    l_sp <- q_sp / 500
    clu_sp <- clup_pv * v_sp_e

    # Other (carcass: stomach, bladder, gallbladder, thyroid, etc.)
    v_ot_vp <- 0.204; v_ot_bc <- 0.167; v_ot_e <- 0.0243; v_ot_is <- 0.831
    q_ot <- 5.521; bcq_ot <- 4.517; sv_ot <- 0.95
    l_ot <- q_ot / 500
    clu_ot <- clup_pv * v_ot_e

    # Lung flows derived for mass-balance closure: q_lu * (1 - 1/500) = sum q_X
    sum_q_X <- q_ht + q_mu + q_sk + q_sc + q_ad + q_bo + q_br + q_ki + q_li +
      q_si + q_lr + q_pa + q_th + q_sp + q_ot
    sum_bcq_X <- bcq_ht + bcq_mu + bcq_sk + bcq_sc + bcq_ad + bcq_bo + bcq_br +
      bcq_ki + bcq_li + bcq_si + bcq_lr + bcq_pa + bcq_th + bcq_sp + bcq_ot
    q_lu <- sum_q_X / (1 - 1 / 500)
    l_lu <- q_lu / 500
    bcq_lu <- sum_bcq_X
    clu_lu <- clup_pv * v_lu_e

    # === Concentration aliases (states are amounts in nmol) ===
    cv_p <- plasma / v_plasma
    cv_bcc <- bcc / v_bcc
    cv_ln <- lnode / v_lnode

    cv_ht <- vp_heart / v_ht_vp; cb_ht <- bc_heart / v_ht_bc; ceu_ht <- eu_heart / v_ht_e
    ceb_ht <- eb_heart / v_ht_e; cfr_ht <- fr_heart / v_ht_e; cis_ht <- is_heart / v_ht_is

    cv_lu <- vp_lung / v_lu_vp; cb_lu <- bc_lung / v_lu_bc; ceu_lu <- eu_lung / v_lu_e
    ceb_lu <- eb_lung / v_lu_e; cfr_lu <- fr_lung / v_lu_e; cis_lu <- is_lung / v_lu_is

    cv_mu <- vp_muscle / v_mu_vp; cb_mu <- bc_muscle / v_mu_bc; ceu_mu <- eu_muscle / v_mu_e
    ceb_mu <- eb_muscle / v_mu_e; cfr_mu <- fr_muscle / v_mu_e; cis_mu <- is_muscle / v_mu_is

    cv_sk <- vp_skin / v_sk_vp; cb_sk <- bc_skin / v_sk_bc; ceu_sk <- eu_skin / v_sk_e
    ceb_sk <- eb_skin / v_sk_e; cfr_sk <- fr_skin / v_sk_e; cis_sk <- is_skin / v_sk_is

    cv_sc <- vp_subcutaneous / v_sc_vp; cb_sc <- bc_subcutaneous / v_sc_bc
    ceu_sc <- eu_subcutaneous / v_sc_e; ceb_sc <- eb_subcutaneous / v_sc_e
    cfr_sc <- fr_subcutaneous / v_sc_e; cis_sc <- is_subcutaneous / v_sc_is

    cv_ad <- vp_adipose / v_ad_vp; cb_ad <- bc_adipose / v_ad_bc; ceu_ad <- eu_adipose / v_ad_e
    ceb_ad <- eb_adipose / v_ad_e; cfr_ad <- fr_adipose / v_ad_e; cis_ad <- is_adipose / v_ad_is

    cv_bo <- vp_bone / v_bo_vp; cb_bo <- bc_bone / v_bo_bc; ceu_bo <- eu_bone / v_bo_e
    ceb_bo <- eb_bone / v_bo_e; cfr_bo <- fr_bone / v_bo_e; cis_bo <- is_bone / v_bo_is

    cv_br <- vp_brain / v_br_vp; cb_br <- bc_brain / v_br_bc; ceu_br <- eu_brain / v_br_e
    ceb_br <- eb_brain / v_br_e; cfr_br <- fr_brain / v_br_e; cis_br <- is_brain / v_br_is

    cv_ki <- vp_kidney / v_ki_vp; cb_ki <- bc_kidney / v_ki_bc; ceu_ki <- eu_kidney / v_ki_e
    ceb_ki <- eb_kidney / v_ki_e; cfr_ki <- fr_kidney / v_ki_e; cis_ki <- is_kidney / v_ki_is

    cv_li <- vp_liver / v_li_vp; cb_li <- bc_liver / v_li_bc; ceu_li <- eu_liver / v_li_e
    ceb_li <- eb_liver / v_li_e; cfr_li <- fr_liver / v_li_e; cis_li <- is_liver / v_li_is

    cv_si <- vp_small_intestine / v_si_vp; cb_si <- bc_small_intestine / v_si_bc
    ceu_si <- eu_small_intestine / v_si_e; ceb_si <- eb_small_intestine / v_si_e
    cfr_si <- fr_small_intestine / v_si_e; cis_si <- is_small_intestine / v_si_is

    cv_lr <- vp_large_intestine / v_lr_vp; cb_lr <- bc_large_intestine / v_lr_bc
    ceu_lr <- eu_large_intestine / v_lr_e; ceb_lr <- eb_large_intestine / v_lr_e
    cfr_lr <- fr_large_intestine / v_lr_e; cis_lr <- is_large_intestine / v_lr_is

    cv_pa <- vp_pancreas / v_pa_vp; cb_pa <- bc_pancreas / v_pa_bc; ceu_pa <- eu_pancreas / v_pa_e
    ceb_pa <- eb_pancreas / v_pa_e; cfr_pa <- fr_pancreas / v_pa_e; cis_pa <- is_pancreas / v_pa_is

    cv_th <- vp_thymus / v_th_vp; cb_th <- bc_thymus / v_th_bc; ceu_th <- eu_thymus / v_th_e
    ceb_th <- eb_thymus / v_th_e; cfr_th <- fr_thymus / v_th_e; cis_th <- is_thymus / v_th_is

    cv_sp <- vp_spleen / v_sp_vp; cb_sp <- bc_spleen / v_sp_bc; ceu_sp <- eu_spleen / v_sp_e
    ceb_sp <- eb_spleen / v_sp_e; cfr_sp <- fr_spleen / v_sp_e; cis_sp <- is_spleen / v_sp_is

    cv_ot <- vp_other / v_ot_vp; cb_ot <- bc_other / v_ot_bc; ceu_ot <- eu_other / v_ot_e
    ceb_ot <- eb_other / v_ot_e; cfr_ot <- fr_other / v_ot_e; cis_ot <- is_other / v_ot_is

    # === Vascular plasma (Shah and Betts 2012 Eq 4; liver Eq 11) ===
    d/dt(vp_heart) <- q_ht * cv_lu - (q_ht - l_ht) * cv_ht - (1 - sv_ht) * l_ht * cv_ht - clu_ht * cv_ht + clu_ht * fr * ceb_ht
    d/dt(vp_lung) <- q_lu * cv_p - (q_lu - l_lu) * cv_lu - (1 - sv_lu) * l_lu * cv_lu - clu_lu * cv_lu + clu_lu * fr * ceb_lu
    d/dt(vp_muscle) <- q_mu * cv_lu - (q_mu - l_mu) * cv_mu - (1 - sv_mu) * l_mu * cv_mu - clu_mu * cv_mu + clu_mu * fr * ceb_mu
    d/dt(vp_skin) <- q_sk * cv_lu - (q_sk - l_sk) * cv_sk - (1 - sv_sk) * l_sk * cv_sk - clu_sk * cv_sk + clu_sk * fr * ceb_sk
    d/dt(vp_subcutaneous) <- q_sc * cv_lu - (q_sc - l_sc) * cv_sc - (1 - sv_sc) * l_sc * cv_sc - clu_sc * cv_sc + clu_sc * fr * ceb_sc
    d/dt(vp_adipose) <- q_ad * cv_lu - (q_ad - l_ad) * cv_ad - (1 - sv_ad) * l_ad * cv_ad - clu_ad * cv_ad + clu_ad * fr * ceb_ad
    d/dt(vp_bone) <- q_bo * cv_lu - (q_bo - l_bo) * cv_bo - (1 - sv_bo) * l_bo * cv_bo - clu_bo * cv_bo + clu_bo * fr * ceb_bo
    d/dt(vp_brain) <- q_br * cv_lu - (q_br - l_br) * cv_br - (1 - sv_br) * l_br * cv_br - clu_br * cv_br + clu_br * fr * ceb_br
    d/dt(vp_kidney) <- q_ki * cv_lu - (q_ki - l_ki) * cv_ki - (1 - sv_ki) * l_ki * cv_ki - clu_ki * cv_ki + clu_ki * fr * ceb_ki
    d/dt(vp_small_intestine) <- q_si * cv_lu - (q_si - l_si) * cv_si - (1 - sv_si) * l_si * cv_si - clu_si * cv_si + clu_si * fr * ceb_si
    d/dt(vp_large_intestine) <- q_lr * cv_lu - (q_lr - l_lr) * cv_lr - (1 - sv_lr) * l_lr * cv_lr - clu_lr * cv_lr + clu_lr * fr * ceb_lr
    d/dt(vp_pancreas) <- q_pa * cv_lu - (q_pa - l_pa) * cv_pa - (1 - sv_pa) * l_pa * cv_pa - clu_pa * cv_pa + clu_pa * fr * ceb_pa
    d/dt(vp_thymus) <- q_th * cv_lu - (q_th - l_th) * cv_th - (1 - sv_th) * l_th * cv_th - clu_th * cv_th + clu_th * fr * ceb_th
    d/dt(vp_spleen) <- q_sp * cv_lu - (q_sp - l_sp) * cv_sp - (1 - sv_sp) * l_sp * cv_sp - clu_sp * cv_sp + clu_sp * fr * ceb_sp
    d/dt(vp_other) <- q_ot * cv_lu - (q_ot - l_ot) * cv_ot - (1 - sv_ot) * l_ot * cv_ot - clu_ot * cv_ot + clu_ot * fr * ceb_ot

    qsum_li <- (q_li - l_li) + (q_sp - l_sp) + (q_pa - l_pa) + (q_si - l_si) + (q_lr - l_lr)
    d/dt(vp_liver) <- q_li * cv_lu +
      (q_sp - l_sp) * cv_sp +
      (q_pa - l_pa) * cv_pa +
      (q_si - l_si) * cv_si +
      (q_lr - l_lr) * cv_lr -
      qsum_li * cv_li -
      (1 - sv_li) * l_li * cv_li -
      clu_li * cv_li +
      clu_li * fr * ceb_li

    # === Vascular blood cells (Eq 5; liver Eq 12) ===
    d/dt(bc_heart) <- bcq_ht * (cb_lu - cb_ht)
    d/dt(bc_lung) <- bcq_lu * (cv_bcc - cb_lu)
    d/dt(bc_muscle) <- bcq_mu * (cb_lu - cb_mu)
    d/dt(bc_skin) <- bcq_sk * (cb_lu - cb_sk)
    d/dt(bc_subcutaneous) <- bcq_sc * (cb_lu - cb_sc)
    d/dt(bc_adipose) <- bcq_ad * (cb_lu - cb_ad)
    d/dt(bc_bone) <- bcq_bo * (cb_lu - cb_bo)
    d/dt(bc_brain) <- bcq_br * (cb_lu - cb_br)
    d/dt(bc_kidney) <- bcq_ki * (cb_lu - cb_ki)
    d/dt(bc_small_intestine) <- bcq_si * (cb_lu - cb_si)
    d/dt(bc_large_intestine) <- bcq_lr * (cb_lu - cb_lr)
    d/dt(bc_pancreas) <- bcq_pa * (cb_lu - cb_pa)
    d/dt(bc_thymus) <- bcq_th * (cb_lu - cb_th)
    d/dt(bc_spleen) <- bcq_sp * (cb_lu - cb_sp)
    d/dt(bc_other) <- bcq_ot * (cb_lu - cb_ot)

    bcqsum_li <- bcq_li + bcq_sp + bcq_pa + bcq_si + bcq_lr
    d/dt(bc_liver) <- bcq_li * cb_lu +
      bcq_sp * cb_sp +
      bcq_pa * cb_pa +
      bcq_si * cb_si +
      bcq_lr * cb_lr -
      bcqsum_li * cb_li

    # === Endosomal space (Eq 6, Eq 7, plus the free-FcRn balance) ===
    bind_ht <- kon * ceu_ht * cfr_ht * v_ht_e; ub_ht <- koff * eb_heart
    d/dt(eu_heart) <- clu_ht * (cv_ht + cis_ht) - bind_ht + ub_ht - kdeg * eu_heart
    d/dt(eb_heart) <- bind_ht - ub_ht - clu_ht * ceb_ht
    d/dt(fr_heart) <- ub_ht - bind_ht + clu_ht * ceb_ht

    bind_lu <- kon * ceu_lu * cfr_lu * v_lu_e; ub_lu <- koff * eb_lung
    d/dt(eu_lung) <- clu_lu * (cv_lu + cis_lu) - bind_lu + ub_lu - kdeg * eu_lung
    d/dt(eb_lung) <- bind_lu - ub_lu - clu_lu * ceb_lu
    d/dt(fr_lung) <- ub_lu - bind_lu + clu_lu * ceb_lu

    bind_mu <- kon * ceu_mu * cfr_mu * v_mu_e; ub_mu <- koff * eb_muscle
    d/dt(eu_muscle) <- clu_mu * (cv_mu + cis_mu) - bind_mu + ub_mu - kdeg * eu_muscle
    d/dt(eb_muscle) <- bind_mu - ub_mu - clu_mu * ceb_mu
    d/dt(fr_muscle) <- ub_mu - bind_mu + clu_mu * ceb_mu

    bind_sk <- kon * ceu_sk * cfr_sk * v_sk_e; ub_sk <- koff * eb_skin
    d/dt(eu_skin) <- clu_sk * (cv_sk + cis_sk) - bind_sk + ub_sk - kdeg * eu_skin
    d/dt(eb_skin) <- bind_sk - ub_sk - clu_sk * ceb_sk
    d/dt(fr_skin) <- ub_sk - bind_sk + clu_sk * ceb_sk

    bind_sc <- kon * ceu_sc * cfr_sc * v_sc_e; ub_sc <- koff * eb_subcutaneous
    d/dt(eu_subcutaneous) <- clu_sc * (cv_sc + cis_sc) - bind_sc + ub_sc - kdeg * eu_subcutaneous
    d/dt(eb_subcutaneous) <- bind_sc - ub_sc - clu_sc * ceb_sc
    d/dt(fr_subcutaneous) <- ub_sc - bind_sc + clu_sc * ceb_sc

    bind_ad <- kon * ceu_ad * cfr_ad * v_ad_e; ub_ad <- koff * eb_adipose
    d/dt(eu_adipose) <- clu_ad * (cv_ad + cis_ad) - bind_ad + ub_ad - kdeg * eu_adipose
    d/dt(eb_adipose) <- bind_ad - ub_ad - clu_ad * ceb_ad
    d/dt(fr_adipose) <- ub_ad - bind_ad + clu_ad * ceb_ad

    bind_bo <- kon * ceu_bo * cfr_bo * v_bo_e; ub_bo <- koff * eb_bone
    d/dt(eu_bone) <- clu_bo * (cv_bo + cis_bo) - bind_bo + ub_bo - kdeg * eu_bone
    d/dt(eb_bone) <- bind_bo - ub_bo - clu_bo * ceb_bo
    d/dt(fr_bone) <- ub_bo - bind_bo + clu_bo * ceb_bo

    bind_br <- kon * ceu_br * cfr_br * v_br_e; ub_br <- koff * eb_brain
    d/dt(eu_brain) <- clu_br * (cv_br + cis_br) - bind_br + ub_br - kdeg * eu_brain
    d/dt(eb_brain) <- bind_br - ub_br - clu_br * ceb_br
    d/dt(fr_brain) <- ub_br - bind_br + clu_br * ceb_br

    bind_ki <- kon * ceu_ki * cfr_ki * v_ki_e; ub_ki <- koff * eb_kidney
    d/dt(eu_kidney) <- clu_ki * (cv_ki + cis_ki) - bind_ki + ub_ki - kdeg * eu_kidney
    d/dt(eb_kidney) <- bind_ki - ub_ki - clu_ki * ceb_ki
    d/dt(fr_kidney) <- ub_ki - bind_ki + clu_ki * ceb_ki

    bind_li <- kon * ceu_li * cfr_li * v_li_e; ub_li <- koff * eb_liver
    d/dt(eu_liver) <- clu_li * (cv_li + cis_li) - bind_li + ub_li - kdeg * eu_liver
    d/dt(eb_liver) <- bind_li - ub_li - clu_li * ceb_li
    d/dt(fr_liver) <- ub_li - bind_li + clu_li * ceb_li

    bind_si <- kon * ceu_si * cfr_si * v_si_e; ub_si <- koff * eb_small_intestine
    d/dt(eu_small_intestine) <- clu_si * (cv_si + cis_si) - bind_si + ub_si - kdeg * eu_small_intestine
    d/dt(eb_small_intestine) <- bind_si - ub_si - clu_si * ceb_si
    d/dt(fr_small_intestine) <- ub_si - bind_si + clu_si * ceb_si

    bind_lr <- kon * ceu_lr * cfr_lr * v_lr_e; ub_lr <- koff * eb_large_intestine
    d/dt(eu_large_intestine) <- clu_lr * (cv_lr + cis_lr) - bind_lr + ub_lr - kdeg * eu_large_intestine
    d/dt(eb_large_intestine) <- bind_lr - ub_lr - clu_lr * ceb_lr
    d/dt(fr_large_intestine) <- ub_lr - bind_lr + clu_lr * ceb_lr

    bind_pa <- kon * ceu_pa * cfr_pa * v_pa_e; ub_pa <- koff * eb_pancreas
    d/dt(eu_pancreas) <- clu_pa * (cv_pa + cis_pa) - bind_pa + ub_pa - kdeg * eu_pancreas
    d/dt(eb_pancreas) <- bind_pa - ub_pa - clu_pa * ceb_pa
    d/dt(fr_pancreas) <- ub_pa - bind_pa + clu_pa * ceb_pa

    bind_th <- kon * ceu_th * cfr_th * v_th_e; ub_th <- koff * eb_thymus
    d/dt(eu_thymus) <- clu_th * (cv_th + cis_th) - bind_th + ub_th - kdeg * eu_thymus
    d/dt(eb_thymus) <- bind_th - ub_th - clu_th * ceb_th
    d/dt(fr_thymus) <- ub_th - bind_th + clu_th * ceb_th

    bind_sp <- kon * ceu_sp * cfr_sp * v_sp_e; ub_sp <- koff * eb_spleen
    d/dt(eu_spleen) <- clu_sp * (cv_sp + cis_sp) - bind_sp + ub_sp - kdeg * eu_spleen
    d/dt(eb_spleen) <- bind_sp - ub_sp - clu_sp * ceb_sp
    d/dt(fr_spleen) <- ub_sp - bind_sp + clu_sp * ceb_sp

    bind_ot <- kon * ceu_ot * cfr_ot * v_ot_e; ub_ot <- koff * eb_other
    d/dt(eu_other) <- clu_ot * (cv_ot + cis_ot) - bind_ot + ub_ot - kdeg * eu_other
    d/dt(eb_other) <- bind_ot - ub_ot - clu_ot * ceb_ot
    d/dt(fr_other) <- ub_ot - bind_ot + clu_ot * ceb_ot

    # === Interstitial space (Eq 9 with antigen-binding terms zero) ===
    d/dt(is_heart) <- (1 - sv_ht) * l_ht * cv_ht - (1 - sigis) * l_ht * cis_ht - clu_ht * cis_ht + clu_ht * (1 - fr) * ceb_ht
    d/dt(is_lung) <- (1 - sv_lu) * l_lu * cv_lu - (1 - sigis) * l_lu * cis_lu - clu_lu * cis_lu + clu_lu * (1 - fr) * ceb_lu
    d/dt(is_muscle) <- (1 - sv_mu) * l_mu * cv_mu - (1 - sigis) * l_mu * cis_mu - clu_mu * cis_mu + clu_mu * (1 - fr) * ceb_mu
    d/dt(is_skin) <- (1 - sv_sk) * l_sk * cv_sk - (1 - sigis) * l_sk * cis_sk - clu_sk * cis_sk + clu_sk * (1 - fr) * ceb_sk

    # Subcutaneous interstitium - Kumar 2024 Section 2.4, the one equation the
    # paper modifies from the IV model. Relative to a generic tissue it adds
    # (a) the S_LU scaling on the convective lymphatic efflux and (b) the
    # first-order local degradation kSC. The SC dose lands in this state.
    sc_to_lymph <- s_lu * (1 - sigis) * l_sc * cis_sc
    d/dt(is_subcutaneous) <- (1 - sv_sc) * l_sc * cv_sc -
      sc_to_lymph -
      clu_sc * cis_sc +
      clu_sc * (1 - fr) * ceb_sc -
      ksc * is_subcutaneous
    d/dt(is_adipose) <- (1 - sv_ad) * l_ad * cv_ad - (1 - sigis) * l_ad * cis_ad - clu_ad * cis_ad + clu_ad * (1 - fr) * ceb_ad
    d/dt(is_bone) <- (1 - sv_bo) * l_bo * cv_bo - (1 - sigis) * l_bo * cis_bo - clu_bo * cis_bo + clu_bo * (1 - fr) * ceb_bo
    d/dt(is_brain) <- (1 - sv_br) * l_br * cv_br - (1 - sigis) * l_br * cis_br - clu_br * cis_br + clu_br * (1 - fr) * ceb_br
    d/dt(is_kidney) <- (1 - sv_ki) * l_ki * cv_ki - (1 - sigis) * l_ki * cis_ki - clu_ki * cis_ki + clu_ki * (1 - fr) * ceb_ki
    d/dt(is_liver) <- (1 - sv_li) * l_li * cv_li - (1 - sigis) * l_li * cis_li - clu_li * cis_li + clu_li * (1 - fr) * ceb_li
    d/dt(is_small_intestine) <- (1 - sv_si) * l_si * cv_si - (1 - sigis) * l_si * cis_si - clu_si * cis_si + clu_si * (1 - fr) * ceb_si
    d/dt(is_large_intestine) <- (1 - sv_lr) * l_lr * cv_lr - (1 - sigis) * l_lr * cis_lr - clu_lr * cis_lr + clu_lr * (1 - fr) * ceb_lr
    d/dt(is_pancreas) <- (1 - sv_pa) * l_pa * cv_pa - (1 - sigis) * l_pa * cis_pa - clu_pa * cis_pa + clu_pa * (1 - fr) * ceb_pa
    d/dt(is_thymus) <- (1 - sv_th) * l_th * cv_th - (1 - sigis) * l_th * cis_th - clu_th * cis_th + clu_th * (1 - fr) * ceb_th
    d/dt(is_spleen) <- (1 - sv_sp) * l_sp * cv_sp - (1 - sigis) * l_sp * cis_sp - clu_sp * cis_sp + clu_sp * (1 - fr) * ceb_sp
    d/dt(is_other) <- (1 - sv_ot) * l_ot * cv_ot - (1 - sigis) * l_ot * cis_ot - clu_ot * cis_ot + clu_ot * (1 - fr) * ceb_ot

    # === Central plasma (Eq 1) ===
    d/dt(plasma) <- (q_ht - l_ht) * cv_ht +
      (q_ki - l_ki) * cv_ki +
      (q_mu - l_mu) * cv_mu +
      (q_sk - l_sk) * cv_sk +
      (q_sc - l_sc) * cv_sc +
      (q_br - l_br) * cv_br +
      (q_ad - l_ad) * cv_ad +
      (q_th - l_th) * cv_th +
      qsum_li * cv_li +
      (q_bo - l_bo) * cv_bo +
      (q_ot - l_ot) * cv_ot +
      l_lnode * cv_ln -
      q_lu * cv_p

    # === Central blood cells (Eq 2) ===
    d/dt(bcc) <- bcq_ht * cb_ht +
      bcq_ki * cb_ki +
      bcq_mu * cb_mu +
      bcq_sk * cb_sk +
      bcq_sc * cb_sc +
      bcq_br * cb_br +
      bcq_ad * cb_ad +
      bcq_th * cb_th +
      bcqsum_li * cb_li +
      bcq_bo * cb_bo +
      bcq_ot * cb_ot -
      bcq_lu * cv_bcc

    # === Lymph node (Eq 3) ===
    d/dt(lnode) <- (1 - sigis) * l_ht * cis_ht +
      (1 - sigis) * l_lu * cis_lu +
      (1 - sigis) * l_mu * cis_mu +
      (1 - sigis) * l_sk * cis_sk +
      sc_to_lymph +
      (1 - sigis) * l_ad * cis_ad +
      (1 - sigis) * l_bo * cis_bo +
      (1 - sigis) * l_br * cis_br +
      (1 - sigis) * l_ki * cis_ki +
      (1 - sigis) * l_li * cis_li +
      (1 - sigis) * l_si * cis_si +
      (1 - sigis) * l_lr * cis_lr +
      (1 - sigis) * l_pa * cis_pa +
      (1 - sigis) * l_th * cis_th +
      (1 - sigis) * l_sp * cis_sp +
      (1 - sigis) * l_ot * cis_ot -
      l_lnode * cv_ln

    # === Initial conditions ===
    # All mAb states start at zero. The SC dose enters `is_subcutaneous` from
    # the event table (cmt = "is_subcutaneous"); an IV dose may still be given
    # into `plasma`. Free FcRn starts at 4.98e-5 M in every tissue endosome.
    fr_heart(0) <- fcrn * v_ht_e
    fr_lung(0) <- fcrn * v_lu_e
    fr_muscle(0) <- fcrn * v_mu_e
    fr_skin(0) <- fcrn * v_sk_e
    fr_subcutaneous(0) <- fcrn * v_sc_e
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

    # === Observation ===
    # Cc is the central plasma concentration in nM. Doses are in nmol;
    # convert from mass via dose_nmol = dose_mg / MW * 1e9 with MW the
    # antibody molar mass (~1.45e5 g/mol for a typical IgG1).
    Cc <- cv_p
  })
}
