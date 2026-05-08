Shah_2012_mAb_PBPK <- function() {
  description <- "Shah and Betts 2012 platform PBPK model for monoclonal antibody plasma and tissue disposition - human (71 kg) parameter set, non-target-binding (nonspecific) mAb. 15 tissues x 6 sub-compartments (vascular plasma, vascular blood cells, endosomal unbound mAb, endosomal FcRn-bound mAb, endosomal free FcRn, interstitial) plus central plasma, central blood cells, and lymph node, totalling 93 ODE states. FcRn-mediated recycling is implemented in every tissue endosomal space."
  reference <- "Shah DK, Betts AM. Towards a platform PBPK model to characterize the plasma and tissue disposition of monoclonal antibodies in preclinical species and human. J Pharmacokinet Pharmacodyn. 2012;39(1):67-86. doi:10.1007/s10928-011-9232-2"
  vignette <- "Shah_2012_mAb_PBPK"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 6,
    n_studies      = 1,
    age_range      = "adult",
    weight_range   = "70 kg reference adult",
    sex_female_pct = NA_real_,
    species        = "Human (71 kg male reference); the paper additionally reports parameter sets for mouse (28 g), rat (280 g), and rhesus monkey (6.2 kg).",
    disease_state  = "Adults with rheumatoid arthritis (validation cohort, Weisman 2003).",
    dose_range     = "Validation: 5 mg/kg IV single dose adalimumab (Weisman 2003 cohort 4).",
    regions        = "Multi-source (digitized mean profiles from Weisman 2003 for human; Bazin-Redureau 1997 for rat; Hinton 2004 for rhesus monkey; Garg/Balthasar 2007, Baxter 1994, Ferl 2005, control IgG/IVIG mouse studies for mouse).",
    notes          = "Platform PBPK model fit simultaneously to 52 plasma and tissue concentration profiles spanning four species (mouse, rat, monkey, human) and several mAbs (control IgG, 7E3, MOPC21, cT84.66, OST577, adalimumab). Four system parameters (FcRn, CLup, Kdeg, C_LNLF) were estimated from pooled data; physiology and FcRn association/dissociation rates were fixed from literature. This model file packages the human (71 kg) parameterization without antigen binding (non-target-mediated disposition). For tumor-targeting or antigen-specific applications, additional terms on cell-membrane binding (Kon_Ag, Koff_Ag, Kint) would be required - see paper Eqs 9-10.",
    scope_note     = "Mechanistic platform model: no IIV, no residual error - intended for typical-value simulation. See vignette for steady-state, mass-balance, and adalimumab plasma replication checks."
  )

  ini({
    # Four system parameters estimated by simultaneous fit to pooled multi-species
    # data (Shah & Betts 2012 Table 6).
    lfcrn   <- log(4.98e-5);  label("Endosomal FcRn concentration (mol/L)")        # Table 6 (CV% 11.1)
    lclup   <- log(0.0366);   label("Pinocytosis rate per unit endosomal volume (L/h/L)") # Table 6 (CV% 3.48)
    lkdeg   <- log(42.9);     label("Endosomal degradation rate of FcRn-unbound mAb (1/h)") # Table 6 (CV% 15.7)
    lclnlf  <- log(9.1);      label("Lymph-node-to-plasma transfer scaling C_LNLF (unitless ratio)") # Table 6 (CV% > 50)

    # Species-specific FcRn-IgG binding (text p.73, human values).
    lkon    <- log(5.59e8);   label("FcRn-IgG association rate constant (1/M/h, human)") # text p.73
    lkoff   <- log(23.9);     label("FcRn-IgG dissociation rate constant (1/h, human)")  # text p.73

    # Lymph-node return flow L_LymphNode (L/h, human). Reported in Table 4 row
    # "Ly. Node" plasma flow column = 3670 mL/h. Treated as a fixed, species-
    # specific physiological flow in this model file; the C_LNLF system
    # parameter above is retained for traceability against Table 6 but does
    # not directly enter the ODEs.
    llymphnode <- log(3.670); label("Lymph-node-to-plasma flow L_LymphNode (L/h, human)") # Table 4 row Ly.Node
  })

  model({
    fcrn_M  <- exp(lfcrn)
    clup_pv <- exp(lclup)
    kdeg    <- exp(lkdeg)
    c_lnlf  <- exp(lclnlf)
    kon_M   <- exp(lkon)
    koff    <- exp(lkoff)
    l_lnode <- exp(llymphnode)

    # Convert FcRn molar concentration to nM and Kon from 1/M/h to 1/(nM h)
    fcrn   <- fcrn_M * 1e9     # nM (49800 nM at 4.98e-5 M)
    kon    <- kon_M * 1e-9     # 1/(nM h)

    # Fixed system constants (Shah & Betts 2012 text p.73, p.77)
    fr     <- 0.715            # Fraction of FcRn-bound mAb recycled to vascular space (Garg & Balthasar 2007)
    sigis  <- 0.2              # Lymphatic reflection coefficient, all tissues (text p.73)

    # === Human (71 kg male) physiological parameters - Shah & Betts 2012 Table 4 ===
    # Volumes are converted from mL to L (mL / 1000); flows from mL/h to L/h (mL/h / 1000).

    # Central pools
    v_plasma <- 3.126
    v_bcc    <- 2.558

    # Lymph node - Total V = Cell V = 274 mL; treated as a single mixing
    # compartment. l_lnode (above) sets the outflow rate.
    v_lnode  <- 0.274

    # Heart
    v_ht_vp <- 0.0131; v_ht_bc <- 0.0108; v_ht_e <- 0.00171; v_ht_is <- 0.0488
    q_ht    <- 7.752;  bcq_ht  <- 6.342;  sv_ht  <- 0.95
    l_ht    <- q_ht / 500
    clu_ht  <- clup_pv * v_ht_e

    # Lung. Volumes from Table 4. Plasma and blood-cell flows are derived from
    # the sum of non-lung tissue arterial inflows so that the lung's vascular
    # mass balance closes exactly. The Shah & Betts 2012 Table 4 reference
    # values (q_lu = 181.913 L/h, bcq_lu = 148.838 L/h) come from a different
    # physiological reference than the tissue rows and do not perfectly equal
    # sum_X q_X + l_lu (they are ~1.8% high), which would cause a slow but
    # non-trivial leak at the lung-arterial junction. Using the sum-derived
    # values keeps total mAb mass conserved (apart from the explicit Kdeg
    # endosomal degradation).
    v_lu_vp <- 0.0550; v_lu_bc <- 0.0450; v_lu_e <- 0.00500; v_lu_is <- 0.300
    sv_lu   <- 0.95

    # Muscle
    v_mu_vp <- 0.662;  v_mu_bc <- 0.541;  v_mu_e <- 0.150;   v_mu_is <- 3.910
    q_mu    <- 33.469; bcq_mu  <- 27.383; sv_mu  <- 0.95
    l_mu    <- q_mu / 500
    clu_mu  <- clup_pv * v_mu_e

    # Skin
    v_sk_vp <- 0.127;  v_sk_bc <- 0.104;  v_sk_e <- 0.0170;  v_sk_is <- 1.125
    q_sk    <- 11.626; bcq_sk  <- 9.512;  sv_sk  <- 0.95
    l_sk    <- q_sk / 500
    clu_sk  <- clup_pv * v_sk_e

    # Adipose
    v_ad_vp <- 0.148;  v_ad_bc <- 0.121;  v_ad_e <- 0.0673;  v_ad_is <- 2.289
    q_ad    <- 11.233; bcq_ad  <- 9.191;  sv_ad  <- 0.95
    l_ad    <- q_ad / 500
    clu_ad  <- clup_pv * v_ad_e

    # Bone
    v_bo_vp <- 0.224;  v_bo_bc <- 0.183;  v_bo_e <- 0.0508;  v_bo_is <- 1.891
    q_bo    <- 2.591;  bcq_bo  <- 2.120;  sv_bo  <- 0.85
    l_bo    <- q_bo / 500
    clu_bo  <- clup_pv * v_bo_e

    # Brain
    v_br_vp <- 0.0319; v_br_bc <- 0.0261; v_br_e <- 0.00725; v_br_is <- 0.261
    q_br    <- 21.453; bcq_br  <- 17.553; sv_br  <- 0.99
    l_br    <- q_br / 500
    clu_br  <- clup_pv * v_br_e

    # Kidney
    v_ki_vp <- 0.0182; v_ki_bc <- 0.0149; v_ki_e <- 0.00166; v_ki_is <- 0.0498
    q_ki    <- 36.402; bcq_ki  <- 29.784; sv_ki  <- 0.90
    l_ki    <- q_ki / 500
    clu_ki  <- clup_pv * v_ki_e

    # Liver (special - receives portal flow from spleen, pancreas, S.Int, L.Int)
    v_li_vp <- 0.183;  v_li_bc <- 0.149;  v_li_e <- 0.0107;  v_li_is <- 0.429
    q_li    <- 13.210; bcq_li  <- 10.808; sv_li  <- 0.85
    l_li    <- q_li / 500
    clu_li  <- clup_pv * v_li_e

    # Small intestine (feeds liver portal)
    v_si_vp <- 0.00615; v_si_bc <- 0.00503; v_si_e <- 0.00193; v_si_is <- 0.0671
    q_si    <- 12.368; bcq_si  <- 10.120; sv_si  <- 0.90
    l_si    <- q_si / 500
    clu_si  <- clup_pv * v_si_e

    # Large intestine (feeds liver portal)
    v_lr_vp <- 0.00874; v_lr_bc <- 0.00715; v_lr_e <- 0.00274; v_lr_is <- 0.0953
    q_lr    <- 12.867; bcq_lr  <- 10.527; sv_lr  <- 0.95
    l_lr    <- q_lr / 500
    clu_lr  <- clup_pv * v_lr_e

    # Pancreas (feeds liver portal)
    v_pa_vp <- 0.00570; v_pa_bc <- 0.00466; v_pa_e <- 0.000518; v_pa_is <- 0.0180
    q_pa    <- 3.056;  bcq_pa  <- 2.500;  sv_pa  <- 0.90
    l_pa    <- q_pa / 500
    clu_pa  <- clup_pv * v_pa_e

    # Thymus
    v_th_vp <- 0.000353; v_th_bc <- 0.000288; v_th_e <- 0.0000321; v_th_is <- 0.00109
    q_th    <- 0.353;  bcq_th  <- 0.289;  sv_th  <- 0.90
    l_th    <- q_th / 500
    clu_th  <- clup_pv * v_th_e

    # Spleen (feeds liver portal)
    v_sp_vp <- 0.0268; v_sp_bc <- 0.0219; v_sp_e <- 0.00111; v_sp_is <- 0.0443
    q_sp    <- 6.343;  bcq_sp  <- 5.189;  sv_sp  <- 0.85
    l_sp    <- q_sp / 500
    clu_sp  <- clup_pv * v_sp_e

    # Other (composite carcass: stomach, bladder, gallbladder, thyroid, etc.)
    v_ot_vp <- 0.204;  v_ot_bc <- 0.167;  v_ot_e <- 0.0243;  v_ot_is <- 0.831
    q_ot    <- 5.521;  bcq_ot  <- 4.517;  sv_ot  <- 0.95
    l_ot    <- q_ot / 500
    clu_ot  <- clup_pv * v_ot_e

    # Derived lung plasma and blood-cell flows (mass-balance closure).
    # q_lu * (1 - 1/500) = sum_X q_X (where X iterates non-lung tissues).
    sum_q_X   <- q_ht + q_mu + q_sk + q_ad + q_bo + q_br + q_ki + q_li +
                 q_si + q_lr + q_pa + q_th + q_sp + q_ot
    sum_bcq_X <- bcq_ht + bcq_mu + bcq_sk + bcq_ad + bcq_bo + bcq_br + bcq_ki +
                 bcq_li + bcq_si + bcq_lr + bcq_pa + bcq_th + bcq_sp + bcq_ot
    q_lu      <- sum_q_X / (1 - 1/500)
    l_lu      <- q_lu / 500
    bcq_lu    <- sum_bcq_X
    clu_lu    <- clup_pv * v_lu_e

    # === Concentration aliases (used in the ODE right-hand sides) ===

    cv_p   <- plasma  / v_plasma
    cv_bcc <- bcc     / v_bcc
    cv_ln  <- lnode   / v_lnode

    cv_ht <- vp_ht / v_ht_vp; cb_ht <- bc_ht / v_ht_bc; ceu_ht <- eu_ht / v_ht_e
    ceb_ht <- eb_ht / v_ht_e; cfr_ht <- fr_ht / v_ht_e; cis_ht <- is_ht / v_ht_is

    cv_lu <- vp_lu / v_lu_vp; cb_lu <- bc_lu / v_lu_bc; ceu_lu <- eu_lu / v_lu_e
    ceb_lu <- eb_lu / v_lu_e; cfr_lu <- fr_lu / v_lu_e; cis_lu <- is_lu / v_lu_is

    cv_mu <- vp_mu / v_mu_vp; cb_mu <- bc_mu / v_mu_bc; ceu_mu <- eu_mu / v_mu_e
    ceb_mu <- eb_mu / v_mu_e; cfr_mu <- fr_mu / v_mu_e; cis_mu <- is_mu / v_mu_is

    cv_sk <- vp_sk / v_sk_vp; cb_sk <- bc_sk / v_sk_bc; ceu_sk <- eu_sk / v_sk_e
    ceb_sk <- eb_sk / v_sk_e; cfr_sk <- fr_sk / v_sk_e; cis_sk <- is_sk / v_sk_is

    cv_ad <- vp_ad / v_ad_vp; cb_ad <- bc_ad / v_ad_bc; ceu_ad <- eu_ad / v_ad_e
    ceb_ad <- eb_ad / v_ad_e; cfr_ad <- fr_ad / v_ad_e; cis_ad <- is_ad / v_ad_is

    cv_bo <- vp_bo / v_bo_vp; cb_bo <- bc_bo / v_bo_bc; ceu_bo <- eu_bo / v_bo_e
    ceb_bo <- eb_bo / v_bo_e; cfr_bo <- fr_bo / v_bo_e; cis_bo <- is_bo / v_bo_is

    cv_br <- vp_br / v_br_vp; cb_br <- bc_br / v_br_bc; ceu_br <- eu_br / v_br_e
    ceb_br <- eb_br / v_br_e; cfr_br <- fr_br / v_br_e; cis_br <- is_br / v_br_is

    cv_ki <- vp_ki / v_ki_vp; cb_ki <- bc_ki / v_ki_bc; ceu_ki <- eu_ki / v_ki_e
    ceb_ki <- eb_ki / v_ki_e; cfr_ki <- fr_ki / v_ki_e; cis_ki <- is_ki / v_ki_is

    cv_li <- vp_li / v_li_vp; cb_li <- bc_li / v_li_bc; ceu_li <- eu_li / v_li_e
    ceb_li <- eb_li / v_li_e; cfr_li <- fr_li / v_li_e; cis_li <- is_li / v_li_is

    cv_si <- vp_si / v_si_vp; cb_si <- bc_si / v_si_bc; ceu_si <- eu_si / v_si_e
    ceb_si <- eb_si / v_si_e; cfr_si <- fr_si / v_si_e; cis_si <- is_si / v_si_is

    cv_lr <- vp_lr / v_lr_vp; cb_lr <- bc_lr / v_lr_bc; ceu_lr <- eu_lr / v_lr_e
    ceb_lr <- eb_lr / v_lr_e; cfr_lr <- fr_lr / v_lr_e; cis_lr <- is_lr / v_lr_is

    cv_pa <- vp_pa / v_pa_vp; cb_pa <- bc_pa / v_pa_bc; ceu_pa <- eu_pa / v_pa_e
    ceb_pa <- eb_pa / v_pa_e; cfr_pa <- fr_pa / v_pa_e; cis_pa <- is_pa / v_pa_is

    cv_th <- vp_th / v_th_vp; cb_th <- bc_th / v_th_bc; ceu_th <- eu_th / v_th_e
    ceb_th <- eb_th / v_th_e; cfr_th <- fr_th / v_th_e; cis_th <- is_th / v_th_is

    cv_sp <- vp_sp / v_sp_vp; cb_sp <- bc_sp / v_sp_bc; ceu_sp <- eu_sp / v_sp_e
    ceb_sp <- eb_sp / v_sp_e; cfr_sp <- fr_sp / v_sp_e; cis_sp <- is_sp / v_sp_is

    cv_ot <- vp_ot / v_ot_vp; cb_ot <- bc_ot / v_ot_bc; ceu_ot <- eu_ot / v_ot_e
    ceb_ot <- eb_ot / v_ot_e; cfr_ot <- fr_ot / v_ot_e; cis_ot <- is_ot / v_ot_is

    # === Vascular plasma ODEs (Eq 4 / Eq 11) ===
    # Arterial input: lung vascular plasma feeds every tissue except itself; the lung
    # is fed from central plasma. Liver also receives portal venous outflow from
    # spleen, pancreas, small intestine, and large intestine (Eq 11).

    d/dt(vp_ht) <- q_ht * cv_lu - (q_ht - l_ht) * cv_ht - (1 - sv_ht) * l_ht * cv_ht - clu_ht * cv_ht + clu_ht * fr * ceb_ht
    d/dt(vp_lu) <- q_lu * cv_p  - (q_lu - l_lu) * cv_lu - (1 - sv_lu) * l_lu * cv_lu - clu_lu * cv_lu + clu_lu * fr * ceb_lu
    d/dt(vp_mu) <- q_mu * cv_lu - (q_mu - l_mu) * cv_mu - (1 - sv_mu) * l_mu * cv_mu - clu_mu * cv_mu + clu_mu * fr * ceb_mu
    d/dt(vp_sk) <- q_sk * cv_lu - (q_sk - l_sk) * cv_sk - (1 - sv_sk) * l_sk * cv_sk - clu_sk * cv_sk + clu_sk * fr * ceb_sk
    d/dt(vp_ad) <- q_ad * cv_lu - (q_ad - l_ad) * cv_ad - (1 - sv_ad) * l_ad * cv_ad - clu_ad * cv_ad + clu_ad * fr * ceb_ad
    d/dt(vp_bo) <- q_bo * cv_lu - (q_bo - l_bo) * cv_bo - (1 - sv_bo) * l_bo * cv_bo - clu_bo * cv_bo + clu_bo * fr * ceb_bo
    d/dt(vp_br) <- q_br * cv_lu - (q_br - l_br) * cv_br - (1 - sv_br) * l_br * cv_br - clu_br * cv_br + clu_br * fr * ceb_br
    d/dt(vp_ki) <- q_ki * cv_lu - (q_ki - l_ki) * cv_ki - (1 - sv_ki) * l_ki * cv_ki - clu_ki * cv_ki + clu_ki * fr * ceb_ki
    d/dt(vp_si) <- q_si * cv_lu - (q_si - l_si) * cv_si - (1 - sv_si) * l_si * cv_si - clu_si * cv_si + clu_si * fr * ceb_si
    d/dt(vp_lr) <- q_lr * cv_lu - (q_lr - l_lr) * cv_lr - (1 - sv_lr) * l_lr * cv_lr - clu_lr * cv_lr + clu_lr * fr * ceb_lr
    d/dt(vp_pa) <- q_pa * cv_lu - (q_pa - l_pa) * cv_pa - (1 - sv_pa) * l_pa * cv_pa - clu_pa * cv_pa + clu_pa * fr * ceb_pa
    d/dt(vp_th) <- q_th * cv_lu - (q_th - l_th) * cv_th - (1 - sv_th) * l_th * cv_th - clu_th * cv_th + clu_th * fr * ceb_th
    d/dt(vp_sp) <- q_sp * cv_lu - (q_sp - l_sp) * cv_sp - (1 - sv_sp) * l_sp * cv_sp - clu_sp * cv_sp + clu_sp * fr * ceb_sp
    d/dt(vp_ot) <- q_ot * cv_lu - (q_ot - l_ot) * cv_ot - (1 - sv_ot) * l_ot * cv_ot - clu_ot * cv_ot + clu_ot * fr * ceb_ot

    # Liver - Eq 11 (special): hepatic-artery inflow + splanchnic portal returns,
    # outflow proportional to combined plasma flow leaving liver to central plasma.
    qsum_li <- (q_li - l_li) + (q_sp - l_sp) + (q_pa - l_pa) + (q_si - l_si) + (q_lr - l_lr)
    d/dt(vp_li) <- q_li * cv_lu +
                   (q_sp - l_sp) * cv_sp +
                   (q_pa - l_pa) * cv_pa +
                   (q_si - l_si) * cv_si +
                   (q_lr - l_lr) * cv_lr -
                   qsum_li * cv_li -
                   (1 - sv_li) * l_li * cv_li -
                   clu_li * cv_li +
                   clu_li * fr * ceb_li

    # === Vascular blood-cell ODEs (Eq 5 / Eq 12) ===
    d/dt(bc_ht) <- bcq_ht * (cb_lu - cb_ht)
    d/dt(bc_lu) <- bcq_lu * (cv_bcc - cb_lu)
    d/dt(bc_mu) <- bcq_mu * (cb_lu - cb_mu)
    d/dt(bc_sk) <- bcq_sk * (cb_lu - cb_sk)
    d/dt(bc_ad) <- bcq_ad * (cb_lu - cb_ad)
    d/dt(bc_bo) <- bcq_bo * (cb_lu - cb_bo)
    d/dt(bc_br) <- bcq_br * (cb_lu - cb_br)
    d/dt(bc_ki) <- bcq_ki * (cb_lu - cb_ki)
    d/dt(bc_si) <- bcq_si * (cb_lu - cb_si)
    d/dt(bc_lr) <- bcq_lr * (cb_lu - cb_lr)
    d/dt(bc_pa) <- bcq_pa * (cb_lu - cb_pa)
    d/dt(bc_th) <- bcq_th * (cb_lu - cb_th)
    d/dt(bc_sp) <- bcq_sp * (cb_lu - cb_sp)
    d/dt(bc_ot) <- bcq_ot * (cb_lu - cb_ot)

    # Liver blood cells - Eq 12 (portal returns from spleen, pancreas, S.Int, L.Int)
    bcqsum_li <- bcq_li + bcq_sp + bcq_pa + bcq_si + bcq_lr
    d/dt(bc_li) <- bcq_li * cb_lu +
                   bcq_sp * cb_sp +
                   bcq_pa * cb_pa +
                   bcq_si * cb_si +
                   bcq_lr * cb_lr -
                   bcqsum_li * cb_li

    # === Endosomal compartment ODEs (Eq 6 / Eq 7 / Eq 8) ===
    # bind_X = Kon * C_eu * C_fr * V_e gives the rate (nmol/h) at which mAb binds
    # FcRn within tissue X's endosomal space.

    bind_ht <- kon * ceu_ht * cfr_ht * v_ht_e; ub_ht <- koff * eb_ht
    d/dt(eu_ht) <- clu_ht * (cv_ht + cis_ht) - bind_ht + ub_ht - kdeg * eu_ht
    d/dt(eb_ht) <- bind_ht - ub_ht - clu_ht * ceb_ht
    d/dt(fr_ht) <- ub_ht - bind_ht + clu_ht * ceb_ht

    bind_lu <- kon * ceu_lu * cfr_lu * v_lu_e; ub_lu <- koff * eb_lu
    d/dt(eu_lu) <- clu_lu * (cv_lu + cis_lu) - bind_lu + ub_lu - kdeg * eu_lu
    d/dt(eb_lu) <- bind_lu - ub_lu - clu_lu * ceb_lu
    d/dt(fr_lu) <- ub_lu - bind_lu + clu_lu * ceb_lu

    bind_mu <- kon * ceu_mu * cfr_mu * v_mu_e; ub_mu <- koff * eb_mu
    d/dt(eu_mu) <- clu_mu * (cv_mu + cis_mu) - bind_mu + ub_mu - kdeg * eu_mu
    d/dt(eb_mu) <- bind_mu - ub_mu - clu_mu * ceb_mu
    d/dt(fr_mu) <- ub_mu - bind_mu + clu_mu * ceb_mu

    bind_sk <- kon * ceu_sk * cfr_sk * v_sk_e; ub_sk <- koff * eb_sk
    d/dt(eu_sk) <- clu_sk * (cv_sk + cis_sk) - bind_sk + ub_sk - kdeg * eu_sk
    d/dt(eb_sk) <- bind_sk - ub_sk - clu_sk * ceb_sk
    d/dt(fr_sk) <- ub_sk - bind_sk + clu_sk * ceb_sk

    bind_ad <- kon * ceu_ad * cfr_ad * v_ad_e; ub_ad <- koff * eb_ad
    d/dt(eu_ad) <- clu_ad * (cv_ad + cis_ad) - bind_ad + ub_ad - kdeg * eu_ad
    d/dt(eb_ad) <- bind_ad - ub_ad - clu_ad * ceb_ad
    d/dt(fr_ad) <- ub_ad - bind_ad + clu_ad * ceb_ad

    bind_bo <- kon * ceu_bo * cfr_bo * v_bo_e; ub_bo <- koff * eb_bo
    d/dt(eu_bo) <- clu_bo * (cv_bo + cis_bo) - bind_bo + ub_bo - kdeg * eu_bo
    d/dt(eb_bo) <- bind_bo - ub_bo - clu_bo * ceb_bo
    d/dt(fr_bo) <- ub_bo - bind_bo + clu_bo * ceb_bo

    bind_br <- kon * ceu_br * cfr_br * v_br_e; ub_br <- koff * eb_br
    d/dt(eu_br) <- clu_br * (cv_br + cis_br) - bind_br + ub_br - kdeg * eu_br
    d/dt(eb_br) <- bind_br - ub_br - clu_br * ceb_br
    d/dt(fr_br) <- ub_br - bind_br + clu_br * ceb_br

    bind_ki <- kon * ceu_ki * cfr_ki * v_ki_e; ub_ki <- koff * eb_ki
    d/dt(eu_ki) <- clu_ki * (cv_ki + cis_ki) - bind_ki + ub_ki - kdeg * eu_ki
    d/dt(eb_ki) <- bind_ki - ub_ki - clu_ki * ceb_ki
    d/dt(fr_ki) <- ub_ki - bind_ki + clu_ki * ceb_ki

    bind_li <- kon * ceu_li * cfr_li * v_li_e; ub_li <- koff * eb_li
    d/dt(eu_li) <- clu_li * (cv_li + cis_li) - bind_li + ub_li - kdeg * eu_li
    d/dt(eb_li) <- bind_li - ub_li - clu_li * ceb_li
    d/dt(fr_li) <- ub_li - bind_li + clu_li * ceb_li

    bind_si <- kon * ceu_si * cfr_si * v_si_e; ub_si <- koff * eb_si
    d/dt(eu_si) <- clu_si * (cv_si + cis_si) - bind_si + ub_si - kdeg * eu_si
    d/dt(eb_si) <- bind_si - ub_si - clu_si * ceb_si
    d/dt(fr_si) <- ub_si - bind_si + clu_si * ceb_si

    bind_lr <- kon * ceu_lr * cfr_lr * v_lr_e; ub_lr <- koff * eb_lr
    d/dt(eu_lr) <- clu_lr * (cv_lr + cis_lr) - bind_lr + ub_lr - kdeg * eu_lr
    d/dt(eb_lr) <- bind_lr - ub_lr - clu_lr * ceb_lr
    d/dt(fr_lr) <- ub_lr - bind_lr + clu_lr * ceb_lr

    bind_pa <- kon * ceu_pa * cfr_pa * v_pa_e; ub_pa <- koff * eb_pa
    d/dt(eu_pa) <- clu_pa * (cv_pa + cis_pa) - bind_pa + ub_pa - kdeg * eu_pa
    d/dt(eb_pa) <- bind_pa - ub_pa - clu_pa * ceb_pa
    d/dt(fr_pa) <- ub_pa - bind_pa + clu_pa * ceb_pa

    bind_th <- kon * ceu_th * cfr_th * v_th_e; ub_th <- koff * eb_th
    d/dt(eu_th) <- clu_th * (cv_th + cis_th) - bind_th + ub_th - kdeg * eu_th
    d/dt(eb_th) <- bind_th - ub_th - clu_th * ceb_th
    d/dt(fr_th) <- ub_th - bind_th + clu_th * ceb_th

    bind_sp <- kon * ceu_sp * cfr_sp * v_sp_e; ub_sp <- koff * eb_sp
    d/dt(eu_sp) <- clu_sp * (cv_sp + cis_sp) - bind_sp + ub_sp - kdeg * eu_sp
    d/dt(eb_sp) <- bind_sp - ub_sp - clu_sp * ceb_sp
    d/dt(fr_sp) <- ub_sp - bind_sp + clu_sp * ceb_sp

    bind_ot <- kon * ceu_ot * cfr_ot * v_ot_e; ub_ot <- koff * eb_ot
    d/dt(eu_ot) <- clu_ot * (cv_ot + cis_ot) - bind_ot + ub_ot - kdeg * eu_ot
    d/dt(eb_ot) <- bind_ot - ub_ot - clu_ot * ceb_ot
    d/dt(fr_ot) <- ub_ot - bind_ot + clu_ot * ceb_ot

    # === Interstitial ODEs (Eq 9, antigen-binding terms zero for nonspecific mAb) ===
    d/dt(is_ht) <- (1 - sv_ht) * l_ht * cv_ht - (1 - sigis) * l_ht * cis_ht - clu_ht * cis_ht + clu_ht * (1 - fr) * ceb_ht
    d/dt(is_lu) <- (1 - sv_lu) * l_lu * cv_lu - (1 - sigis) * l_lu * cis_lu - clu_lu * cis_lu + clu_lu * (1 - fr) * ceb_lu
    d/dt(is_mu) <- (1 - sv_mu) * l_mu * cv_mu - (1 - sigis) * l_mu * cis_mu - clu_mu * cis_mu + clu_mu * (1 - fr) * ceb_mu
    d/dt(is_sk) <- (1 - sv_sk) * l_sk * cv_sk - (1 - sigis) * l_sk * cis_sk - clu_sk * cis_sk + clu_sk * (1 - fr) * ceb_sk
    d/dt(is_ad) <- (1 - sv_ad) * l_ad * cv_ad - (1 - sigis) * l_ad * cis_ad - clu_ad * cis_ad + clu_ad * (1 - fr) * ceb_ad
    d/dt(is_bo) <- (1 - sv_bo) * l_bo * cv_bo - (1 - sigis) * l_bo * cis_bo - clu_bo * cis_bo + clu_bo * (1 - fr) * ceb_bo
    d/dt(is_br) <- (1 - sv_br) * l_br * cv_br - (1 - sigis) * l_br * cis_br - clu_br * cis_br + clu_br * (1 - fr) * ceb_br
    d/dt(is_ki) <- (1 - sv_ki) * l_ki * cv_ki - (1 - sigis) * l_ki * cis_ki - clu_ki * cis_ki + clu_ki * (1 - fr) * ceb_ki
    d/dt(is_li) <- (1 - sv_li) * l_li * cv_li - (1 - sigis) * l_li * cis_li - clu_li * cis_li + clu_li * (1 - fr) * ceb_li
    d/dt(is_si) <- (1 - sv_si) * l_si * cv_si - (1 - sigis) * l_si * cis_si - clu_si * cis_si + clu_si * (1 - fr) * ceb_si
    d/dt(is_lr) <- (1 - sv_lr) * l_lr * cv_lr - (1 - sigis) * l_lr * cis_lr - clu_lr * cis_lr + clu_lr * (1 - fr) * ceb_lr
    d/dt(is_pa) <- (1 - sv_pa) * l_pa * cv_pa - (1 - sigis) * l_pa * cis_pa - clu_pa * cis_pa + clu_pa * (1 - fr) * ceb_pa
    d/dt(is_th) <- (1 - sv_th) * l_th * cv_th - (1 - sigis) * l_th * cis_th - clu_th * cis_th + clu_th * (1 - fr) * ceb_th
    d/dt(is_sp) <- (1 - sv_sp) * l_sp * cv_sp - (1 - sigis) * l_sp * cis_sp - clu_sp * cis_sp + clu_sp * (1 - fr) * ceb_sp
    d/dt(is_ot) <- (1 - sv_ot) * l_ot * cv_ot - (1 - sigis) * l_ot * cis_ot - clu_ot * cis_ot + clu_ot * (1 - fr) * ceb_ot

    # === Central plasma (Eq 1) ===
    # Sums venous returns from tissues whose blood drains directly to plasma
    # (heart, kidney, muscle, skin, brain, adipose, thymus, bone, other) plus
    # the liver venous outflow that bundles its own plus splanchnic portal
    # tributaries, plus return from the lymph node, minus arterial outflow to
    # the lung.
    d/dt(plasma) <- (q_ht - l_ht) * cv_ht +
                    (q_ki - l_ki) * cv_ki +
                    (q_mu - l_mu) * cv_mu +
                    (q_sk - l_sk) * cv_sk +
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
    # All mAb states start at zero (IV bolus dose enters plasma via the data
    # event table). Free FcRn at fcrn = 4.98e-5 M = fcrn nM in every tissue
    # endosomal space; bound FcRn-mAb starts at zero.
    fr_ht(0) <- fcrn * v_ht_e
    fr_lu(0) <- fcrn * v_lu_e
    fr_mu(0) <- fcrn * v_mu_e
    fr_sk(0) <- fcrn * v_sk_e
    fr_ad(0) <- fcrn * v_ad_e
    fr_bo(0) <- fcrn * v_bo_e
    fr_br(0) <- fcrn * v_br_e
    fr_ki(0) <- fcrn * v_ki_e
    fr_li(0) <- fcrn * v_li_e
    fr_si(0) <- fcrn * v_si_e
    fr_lr(0) <- fcrn * v_lr_e
    fr_pa(0) <- fcrn * v_pa_e
    fr_th(0) <- fcrn * v_th_e
    fr_sp(0) <- fcrn * v_sp_e
    fr_ot(0) <- fcrn * v_ot_e

    # === Observation ===
    # Cc reports the central plasma concentration (nM). The state amounts are in
    # nmol; users dosing in mg should convert via dose_nmol = dose_mg / MW * 1e9
    # where MW is the antibody molar mass (e.g., 144,190 g/mol for adalimumab).
    Cc <- cv_p
  })
}
