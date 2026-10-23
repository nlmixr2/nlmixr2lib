Zaid_2021_pb212dotamtate_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, SAAM II / SimBiology). Whole-body",
    "physiologically based pharmacokinetic model of the in vivo alpha-particle",
    "generator [212Pb]Pb-DOTAMTATE, a somatostatin analogue targeting sst2, and",
    "of its radioactive decay products in athymic nude mice bearing AR42J",
    "neuroendocrine tumour xenografts. 308 ODE states, amounts in nmol and time",
    "in min. The labelled system carries 212Pb-peptide, 212Bi-peptide (from",
    "beta decay of 212Pb to the 212Bi ground state), free 212Bi (released when",
    "212Pb decays to an excited 212Bi state, dissociation fraction 0.16),",
    "free 212Bi bound to high-molecular-weight plasma proteins and to red",
    "blood cells, and free 208Tl (from the alpha branch of 212Bi); a second,",
    "unlabelled-peptide system competes for the same sst2 sites and is fed by",
    "the decays that strip the radionuclide from the peptide. Seven",
    "sst2-expressing tissues (tumour, liver, spleen, pancreas, GI tract,",
    "adrenal, lung) carry vascular, interstitial, receptor-bound, early-",
    "endosomal and sorting compartments; the kidneys bind the peptide from",
    "their vascular space (glomerular podocytes) and add a proximal tubule,",
    "a non-specific tubular uptake compartment and urinary excretion; six",
    "sst2-negative tissues (muscle, skin, heart, bone, fat, remainder) carry",
    "vascular, interstitial and a free-212Bi cellular compartment; the brain",
    "is vascular-only. Spleen, pancreas and GI tract drain into the liver.",
    "Six sst2 densities, tumour perfusion, the sorting rate and the pancreas",
    "and kidney release rates were fitted to %ID/g biodistribution data; no",
    "interindividual or residual variability was estimated. The blood pools",
    "and the unlabelled system are not written out in the source and are",
    "reconstructed by flow balance and by the cited two-system framework.",
    sep = " "
  )
  reference <- paste(
    "Zaid NRR, Kletting P, Winter G, Prasad V, Beer AJ, Glatting G.",
    "A Physiologically Based Pharmacokinetic Model for In Vivo Alpha Particle",
    "Generators Targeting Neuroendocrine Tumors in Mice. Pharmaceutics.",
    "2021;13(12):2132. doi:10.3390/pharmaceutics13122132. Model equations",
    "(Eqs 1-64) and parameters (Table S1) are in Supplement A. Remainder-of-body",
    "flow and red-cell volume definitions follow the same authors' 212Bi rat",
    "PBPK (Zaid NRR et al. Cancer Biother Radiopharm 2021,",
    "doi:10.1089/cbr.2021.0028, Table A1); the two-system labelled/unlabelled",
    "structure follows Kletting P et al.",
    "J Nucl Med 2016;57(4):503-8, doi:10.2967/jnumed.115.164699.",
    sep = " "
  )
  vignette <- "Zaid_2021_pb212dotamtate_mouse_pbpk"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/g")
  # The intravenous dose enters the venous serum pool, as in the source
  # (Methods 2.3). Give the radiolabelled [212Pb]Pb-DOTAMTATE amount to
  # pb_ven and the co-injected unlabelled peptide to cold_ven; the
  # non-generator comparator [212Bi]Bi-DOTAMTATE is dosed to bis_ven.
  dosing <- c("pb_ven", "cold_ven", "bis_ven")

  paper_specific_compartment_pattern <- "^(pb|bis|bi|bihwpp|birbc|tl|cold)_"

  population <- list(
    species = "mouse (athymic nude, AR42J rat pancreatic tumour xenograft)",
    n_subjects = NA,
    n_studies = 1,
    disease_state = "AR42J (sst2-expressing) subcutaneous xenograft, 300 uL tumour volume assumed in the model.",
    dose_range = "0.0013 nmol [212Pb]Pb-DOTAMTATE (185 kBq, 5 uCi at 2.4 uCi/ng) intravenous bolus for the fitted study; 0.017e-3 to 0.079 nmol simulated for dosimetry.",
    regions = "Fitted to the published biodistribution of Stallons et al. (Mol Cancer Ther 2019), United States.",
    notes = paste(
      "A standard 30 g mouse; every physiological volume and flow is a",
      "literature value scaled to body weight (Table S1). The fitted data are",
      "the non-decay-corrected %ID/g in kidneys, liver, spleen, pancreas,",
      "tumour and lung at 1, 4 and 24 h (Table S2), one mean per time point.",
      "The number of animals per time point is not reported in the source.",
      sep = " "
    )
  )

  # Every state is an amount in nmol. The specimen follows the
  # sub-compartment: vascular and blood-pool states are serum (or the
  # red-cell fraction for red-cell-bound 212Bi), early-endosomal and sorting
  # states are endosomes, urinary excretion is urine, and everything else is
  # tissue (tumour for the tumour).
  compartmentData <- local({
    vsp <- c("pb", "bis", "bi", "bihwpp", "birbc", "tl", "cold")
    isp <- c("pb", "bis", "bi", "tl", "cold")
    bsp <- c("pb", "bis", "cold")
    sstr2 <- c("tumor", "liver", "spleen", "pancreas", "gi", "adrenal", "lung")
    nonsstr2 <- c("muscle", "skin", "heart", "bone", "fat", "rb")
    nm <- c(
      paste0(vsp, "_art"),
      paste0(vsp, "_ven"),
      as.vector(outer(vsp, c(sstr2, "kidney", nonsstr2), paste, sep = "_v_")),
      as.vector(outer(isp, c(sstr2, nonsstr2), paste, sep = "_int_")),
      as.vector(outer(bsp, c(sstr2, "kidney"), paste, sep = "_bnd_")),
      as.vector(outer(isp, c(sstr2, "kidney"), paste, sep = "_endo_")),
      as.vector(outer(isp, c(sstr2, "kidney"), paste, sep = "_sort_")),
      paste0(isp, "_tub_kidney"),
      paste0(isp, "_nsp_kidney"),
      paste0(isp, "_exc_kidney"),
      as.vector(outer(c("bi", "tl"), nonsstr2, paste, sep = "_cell_"))
    )
    analyteMap <- c(
      pb = "[212Pb]Pb-DOTAMTATE",
      bis = "[212Bi]Bi-DOTAMTATE",
      bi = "free 212Bi",
      bihwpp = "212Bi bound to high-molecular-weight plasma proteins",
      birbc = "212Bi bound to red blood cells",
      tl = "free 208Tl",
      cold = "unlabelled DOTAMTATE"
    )
    specimen <- ifelse(
      grepl("^birbc_", nm),
      "blood cell",
      ifelse(
        grepl("(_v_|_art$|_ven$)", nm),
        "serum",
        ifelse(
          grepl("_(endo|sort)_", nm),
          "endosome",
          ifelse(
            grepl("_exc_", nm),
            "urine",
            ifelse(grepl("_tumor$", nm), "tumor", "tissue")
          )
        )
      )
    )
    stats::setNames(
      Map(
        list,
        analyte = unname(analyteMap[sub("_.*$", "", nm)]),
        units = "nmol",
        specimen = specimen,
        verified = TRUE
      ),
      nm
    )
  })

  ini({
    # =====================================================================
    # Fitted parameters, Table 1. SAAM II and SimBiology agree to the
    # printed precision except the liver sst2 density (0.17 SAAM II, 0.16
    # SimBiology); the SimBiology value is used because the dosimetry of
    # Table 2 and Figures 5-6 was computed in SimBiology.
    # =====================================================================
    lrdens_kidney   <- log(3.52)     ; label("Log sst2 density, kidneys (nmol/L)")             # Table 1, 3.52 +/- 0.07
    lrdens_liver    <- log(0.16)     ; label("Log sst2 density, liver (nmol/L)")               # Table 1, SimBiology 0.16 (SAAM II 0.17 +/- 0.04)
    lrdens_spleen   <- log(0.73)     ; label("Log sst2 density, spleen (nmol/L)")              # Table 1, 0.73 +/- 0.07
    lrdens_lung     <- log(1.77)     ; label("Log sst2 density, lung (nmol/L)")                # Table 1, 1.77 +/- 0.14
    lrdens_pancreas <- log(6.18)     ; label("Log sst2 density, pancreas (nmol/L)")            # Table 1, 6.18 +/- 0.94
    lrdens_tumor    <- log(11.73)    ; label("Log sst2 density, tumour (nmol/L)")              # Table 1, 11.73 +/- 2.97
    lfperf_tumor    <- log(0.09)     ; label("Log tumour serum perfusion (mL/min/g)")          # Table 1, 0.09 +/- 0.03
    lksort          <- log(0.0076)   ; label("Log sorting rate, sst2-expressing tissues (1/min)") # Table 1, 0.0076 +/- 0.0003
    lkrel_pancreas  <- log(0.0011)   ; label("Log release rate, pancreas (1/min)")             # Table 1, 0.0011 +/- 0.0002
    lkrel_kidney    <- log(0.000380) ; label("Log release rate, kidneys (1/min)")              # Table 1, 0.000380 +/- 0.000001

    # =====================================================================
    # Peptide- and radionuclide-specific parameters, Table S1.
    # =====================================================================
    koff     <- fixed(0.04)    ; label("sst2 dissociation rate (1/min)")                          # Table S1, k_off = 0.04
    kd       <- fixed(1)       ; label("sst2 dissociation constant (nmol/L)")                     # Table S1, KD = 1; k_on = k_off / KD
    lam_pb   <- fixed(0.00109) ; label("Physical decay rate of 212Pb (1/min)")                    # Table S1, lambda_phy(212Pb) = 0.00109
    f_diss   <- fixed(0.16)    ; label("Fraction of 212Pb decays to excited 212Bi that free the 212Bi from the chelator") # Table S1, f = 0.16
    lam_bib  <- fixed(0.00732) ; label("Beta decay rate of 212Bi (1/min)")                        # Table S1, lambda_phy,beta(212Bi) = 0.00732
    lam_bia  <- fixed(0.00412) ; label("Alpha decay rate of 212Bi (1/min)")                       # Table S1, lambda_phy,alpha(212Bi) = 0.00412
    lam_tl   <- fixed(0.22360) ; label("Beta decay rate of 208Tl (1/min)")                        # Table S1, lambda_phy,beta(208Tl) = 0.22360
    kint     <- fixed(0.17)    ; label("Internalisation rate of the bound peptide (1/min)")       # Table S1, lambda_intern = 0.17; Methods 2.1
    krec     <- fixed(0.05)    ; label("Recycling rate from early endosomes (1/min)")             # Table S1, lambda_rec = 0.05; Methods 2.1
    krel_nt  <- fixed(0)       ; label("Release rate, normal tissues other than pancreas and kidneys (1/min)") # Table S1, lambda_rel,NT = 0.0 (preliminary fit)
    krel_tum <- fixed(0)       ; label("Release rate, tumour (1/min)")                            # Table S1, lambda_rel,T = 0.0 (preliminary fit)
    rratio_adrenal <- fixed(0.4)  ; label("Adrenal sst2 density relative to pancreas")            # Table S1, [R_Ad] = 0.4 * [R_pan]
    rratio_gi      <- fixed(0.75) ; label("GI sst2 density relative to pancreas")                 # Table S1, [R_GI] = 0.75 * [R_pan]

    # ---- Free 212Bi binding in blood (Table S1, from the rat 212Bi PBPK) --
    kass_rbc   <- fixed(0.01)  ; label("Association rate of free 212Bi with red blood cells (1/min)") # Table S1, Bi-RBC-ass = 0.01
    kdiss_rbc  <- fixed(0.017) ; label("Dissociation rate of red-cell-bound 212Bi (1/min)")           # Table S1, Bi-RBC-diss = 0.017
    kass_hwpp  <- fixed(0.20)  ; label("Association rate of free 212Bi with HWPP (1/min)")            # Table S1, Bi-HWPP-ass = 0.20
    kdiss_hwpp <- fixed(0.20)  ; label("Dissociation rate of HWPP-bound 212Bi (1/min)")               # Table S1, Bi-HWPP-diss = 0.20 (cited rat source prints 0.11; see vignette)

    # ---- Free 212Bi cellular uptake (Table S1) -----------------------------
    kup_liver  <- fixed(0.86)  ; label("Free 212Bi uptake rate, liver cells (1/min)")                 # Table S1, lambda_up,liv = 0.86
    kup_gi     <- fixed(0.27)  ; label("Free 212Bi uptake rate, GI cells; spleen, pancreas and lung share it (1/min)") # Table S1, lambda_up,GI = 0.27
    kup_bone   <- fixed(3.85)  ; label("Free 212Bi uptake rate, bone cells (1/min)")                  # Table S1, lambda_up,Bone = 3.85
    kup_skin   <- fixed(0.04)  ; label("Free 212Bi uptake rate, skin cells (1/min)")                  # Table S1, lambda_up,Skin = 0.04
    kup_muscle <- fixed(0.007) ; label("Free 212Bi uptake rate, muscle cells; fat, heart and remainder share it (1/min)") # Table S1, lambda_up,Mus = 0.007
    kup_kidney <- fixed(1)     ; label("Free 212Bi uptake rate, proximal tubule epithelium (1/min)")  # Table S1, lambda_up,K = 1
    kup_tumor  <- fixed(0)     ; label("Free 212Bi uptake rate, tumour cells (1/min)")                # not tabulated in Table S1; set to 0 by the maintainers (see vignette)
    kup_adrenal <- fixed(0)    ; label("Free 212Bi uptake rate, adrenal cells (1/min)")               # not tabulated in Table S1; set to 0 by the maintainers (see vignette)

    # ---- Permeability-surface densities (Table S1) -------------------------
    kps_muscle_p <- fixed(0.02) ; label("Peptide permeability-surface density, muscle (mL/min/g)")          # Table S1, k_Mus,p = 0.02
    kps_muscle_f <- fixed(0.22) ; label("Free-nuclide permeability-surface density, muscle (mL/min/g)")     # Table S1, k_Mus,f = 0.22
    kps_tumor_p  <- fixed(0.31) ; label("Peptide permeability-surface density, tumour (mL/min/g)")          # Table S1, k_T = 0.31

    # ---- Kidney (Table S1) -------------------------------------------------
    gfr          <- fixed(0.0001)    ; label("Glomerular filtration rate (L/min)")                  # Table S1, GFR = 0.0001
    theta_p      <- fixed(0.97)      ; label("Sieving coefficient ratio, peptide")                  # Table S1, Theta_p = 0.97
    fexc_p       <- fixed(0.98)      ; label("Excretion fraction, peptide")                         # Table S1, f_exc,p = 0.98
    theta_f      <- fixed(1)         ; label("Sieving coefficient ratio, free radionuclides")       # Table S1, Theta_f = 1
    fexc_f       <- fixed(0)         ; label("Excretion fraction, free radionuclides")              # Table S1, f_exc,f = 0.0
    v_tub        <- fixed(0.0001449) ; label("Volume of the kidney proximal tubules (L)")           # Table S1, V_K_Prox_Tub = 0.0001449
    kendo_kidney <- fixed(0.03)      ; label("Renal non-specific tubular uptake rate (1/min)")      # Table S1, lambda_Endo,k = 0.03; Methods 2.1

    # ---- Mouse physiology (Table S1) ---------------------------------------
    bw              <- fixed(0.03)  ; label("Body weight of the standard mouse (kg)")         # Table S1, BW = 0.03
    hct             <- fixed(0.45)  ; label("Hematocrit (fraction)")                          # Table S1, H = 0.45
    fsc             <- fixed(0.86)  ; label("Red-cell volume scaling factor")                 # Table S1, f_sc = 0.86
    vtot_tumor      <- fixed(300e-6); label("Total tumour volume (L)")                        # Table S1, V_T,total = 300e-6
    vfrac_tumor_int <- fixed(0.4)   ; label("Tumour interstitial volume fraction")            # Table S1, v_t,int = 0.4
    vfrac_tumor_v   <- fixed(0.1)   ; label("Tumour vascular volume fraction")                # Table S1, v_t,v = 0.1
  })

  model({
    # =====================================================================
    # 1. Fitted and derived kinetic parameters
    # =====================================================================
    rdens_kidney   <- exp(lrdens_kidney)
    rdens_liver    <- exp(lrdens_liver)
    rdens_spleen   <- exp(lrdens_spleen)
    rdens_lung     <- exp(lrdens_lung)
    rdens_pancreas <- exp(lrdens_pancreas)
    rdens_tumor    <- exp(lrdens_tumor)
    rdens_adrenal  <- rratio_adrenal * rdens_pancreas   # Table S1
    rdens_gi       <- rratio_gi * rdens_pancreas        # Table S1
    fperf_tumor    <- exp(lfperf_tumor)
    ksort          <- exp(lksort)
    krel_pancreas  <- exp(lkrel_pancreas)
    krel_kidney    <- exp(lkrel_kidney)
    krel_tumor     <- krel_tum
    krel_liver     <- krel_nt
    krel_spleen    <- krel_nt
    krel_gi        <- krel_nt
    krel_adrenal   <- krel_nt
    krel_lung      <- krel_nt

    kon    <- koff / kd                 # Table S1, k_on = k_off / KD (L/nmol/min)
    lam_ex <- f_diss * lam_pb           # Table S1, decay of 212Pb to excited 212Bi states (Bi leaves the chelator)
    lam_gr <- (1 - f_diss) * lam_pb     # Table S1, decay of 212Pb to the 212Bi ground state (Bi stays chelated)
    lam_bi <- lam_bib + lam_bia         # total 212Bi decay rate

    kup_spleen   <- kup_gi              # Table S1, lambda_up,Sp = lambda_up,GI
    kup_pancreas <- kup_gi              # Table S1, lambda_up,Panc = lambda_up,GI
    kup_lung     <- kup_gi              # Table S1, lambda_up,Lung = lambda_up,GI
    kup_fat      <- kup_muscle          # Table S1, lambda_up,Fat = lambda_up,Mus
    kup_heart    <- kup_muscle          # Table S1, lambda_up,Hrt = lambda_up,Mus
    kup_rb       <- kup_muscle          # Table S1, lambda_up,RB = lambda_up,Mus

    # =====================================================================
    # 2. Physiology, Table S1. Volumes in L, flows in L/min. Table S1 gives
    #    perfusion and permeability densities in mL/min/g with 1 mL = 1 g,
    #    so a density times a volume in L is already a flow in L/min.
    # =====================================================================
    vp     <- 0.059 * bw                # total body serum
    co     <- 0.275 * bw^0.75           # cardiac output (blood)
    qp_tot <- (1 - hct) * co            # total serum flow
    qr_tot <- hct * co                  # total red-cell flow
    vrbc   <- 0.0276 * bw               # total red-cell volume

    # ---- Total tissue volumes -------------------------------------------
    vtot_kidney   <- 0.019 * bw
    vtot_liver    <- 0.069 * bw
    vtot_spleen   <- 0.005 * bw
    vtot_pancreas <- 0.003 * bw
    vtot_gi       <- 0.037 * bw
    vtot_adrenal  <- 0.0005 * bw
    vtot_lung     <- 0.007 * bw
    vtot_muscle   <- 0.403 * bw
    vtot_skin     <- 0.188 * bw
    vtot_heart    <- 0.005 * bw
    vtot_brain    <- 0.017 * bw
    vtot_bone     <- (0.054 / 1.92) * bw    # marrow-free bone density 1.92 g/cm^3
    vtot_fat      <- (0.071 / 0.92) * bw    # fat density 0.92 g/cm^3
    # Remainder: body volume (1 mL = 1 g) minus every tissue except the tumour.
    vtot_rb <- bw - vtot_kidney - vtot_liver - vtot_spleen - vtot_pancreas -
      vtot_gi - vtot_adrenal - vtot_lung - vtot_muscle - vtot_skin -
      vtot_heart - vtot_brain - vtot_bone - vtot_fat

    # ---- Vascular serum volumes -------------------------------------------
    vv_tumor    <- vfrac_tumor_v * vtot_tumor * (1 - hct)
    vv_kidney   <- 0.055 * vtot_kidney
    vv_liver    <- 0.085 * vtot_liver
    vv_spleen   <- 0.121 * vtot_spleen
    vv_pancreas <- 0.055 * vtot_pancreas
    vv_gi       <- 0.032 * vtot_gi
    vv_adrenal  <- 0.03 * (1 - hct) * vtot_adrenal
    vv_lung     <- 0.145 * vtot_lung
    vv_muscle   <- 0.022 * vtot_muscle
    vv_skin     <- 0.037 * vtot_skin
    vv_heart    <- 0.038 * vtot_heart
    vv_brain    <- 0.022 * vtot_brain
    vv_bone     <- 0.11 * (1 - hct) * vtot_bone
    vv_fat      <- 0.011 * vtot_fat
    vv_art      <- 0.06 * vp + 0.045 * vp   # arterial serum plus half the cardiac serum
    vv_ven      <- 0.18 * vp + 0.045 * vp   # venous serum plus half the cardiac serum
    # Remainder serum: total serum less the blood pools and every tissue
    # except the tumour (Kletting 2016 framework; see vignette).
    vv_rb <- vp - vv_art - vv_ven - vv_kidney - vv_liver - vv_spleen -
      vv_pancreas - vv_gi - vv_adrenal - vv_lung - vv_muscle - vv_skin -
      vv_heart - vv_brain - vv_bone - vv_fat

    # ---- Vascular red-cell volumes --------------------------------------
    vr_tumor    <- vfrac_tumor_v * vtot_tumor * (hct / (1 - hct))   # as printed in Table S1
    vr_kidney   <- fsc * 0.045 * vtot_kidney
    vr_liver    <- fsc * 0.069 * vtot_liver
    vr_spleen   <- fsc * 0.099 * vtot_spleen
    vr_pancreas <- fsc * 0.045 * vtot_pancreas
    vr_gi       <- fsc * 0.01 * vtot_gi
    vr_adrenal  <- fsc * 0.03 * hct * vtot_adrenal
    vr_lung     <- fsc * 0.118 * vtot_lung
    vr_muscle   <- fsc * 0.018 * vtot_muscle
    vr_skin     <- fsc * 0.031 * vtot_skin
    vr_heart    <- fsc * 0.032 * vtot_heart
    vr_brain    <- fsc * 0.018 * vtot_brain
    vr_bone     <- fsc * 0.01 * hct * vtot_bone
    vr_fat      <- fsc * 0.008 * vtot_fat
    vr_art      <- fsc * (0.06 * vp + 0.045 * vp)
    vr_ven      <- fsc * (0.18 * vp + 0.045 * vp)
    # Remainder red cells: total red-cell volume less every tissue except
    # the tumour (rat Table A1 definition; the blood pools are not
    # subtracted -- see vignette).
    vr_rb <- vrbc - vr_kidney - vr_liver - vr_spleen - vr_pancreas - vr_gi -
      vr_adrenal - vr_lung - vr_muscle - vr_skin - vr_heart - vr_brain -
      vr_bone - vr_fat

    # ---- Interstitial volumes ---------------------------------------------
    vi_tumor    <- vfrac_tumor_int * vtot_tumor
    vi_liver    <- 0.2 * vtot_liver
    vi_spleen   <- 0.2 * vtot_spleen
    vi_pancreas <- 0.17 * vtot_pancreas
    vi_gi       <- 0.35 * vtot_gi
    vi_adrenal  <- 0.2 * vtot_adrenal
    vi_lung     <- 0.19 * vtot_lung
    vi_muscle   <- 0.13 * vtot_muscle
    vi_skin     <- 0.33 * vtot_skin
    vi_heart    <- 0.14 * vtot_heart
    vi_bone     <- 0.19 * vtot_bone
    vi_fat      <- 0.17 * vtot_fat
    vi_rb       <- 0.17 * vtot_rb

    # ---- Serum and red-cell flows (F_i,p = (1-H) CO f_i; F_i,RBC = H CO f_i)
    qp_kidney   <- 0.184 * qp_tot
    qp_liver    <- 0.028 * qp_tot    # hepatic artery
    qp_spleen   <- 0.022 * qp_tot
    qp_pancreas <- 0.017 * qp_tot
    qp_gi       <- 0.135 * qp_tot
    qp_adrenal  <- 0.0027 * qp_tot
    qp_muscle   <- 0.231 * qp_tot
    qp_skin     <- 0.075 * qp_tot
    qp_heart    <- 0.098 * qp_tot
    qp_bone     <- 0.041 * qp_tot
    qp_fat      <- 0.036 * qp_tot
    qp_lung     <- qp_tot             # f_Lun = 1
    qp_tumor    <- fperf_tumor * vtot_tumor   # F_T = f_T * V_T,total (serum)
    qp_rb <- qp_tot - qp_kidney - qp_liver - qp_spleen - qp_pancreas - qp_gi -
      qp_adrenal - qp_muscle - qp_skin - qp_heart - qp_bone - qp_fat - qp_tumor
    qp_livout <- qp_liver + qp_spleen + qp_pancreas + qp_gi   # portal drainage

    qr_kidney   <- 0.184 * qr_tot
    qr_liver    <- 0.028 * qr_tot
    qr_spleen   <- 0.022 * qr_tot
    qr_pancreas <- 0.017 * qr_tot
    qr_gi       <- 0.135 * qr_tot
    qr_adrenal  <- 0.0027 * qr_tot
    qr_muscle   <- 0.231 * qr_tot
    qr_skin     <- 0.075 * qr_tot
    qr_heart    <- 0.098 * qr_tot
    qr_bone     <- 0.041 * qr_tot
    qr_fat      <- 0.036 * qr_tot
    qr_lung     <- qr_tot
    qr_tumor    <- qp_tumor * hct / (1 - hct)
    qr_rb <- qr_tot - qr_kidney - qr_liver - qr_spleen - qr_pancreas - qr_gi -
      qr_adrenal - qr_muscle - qr_skin - qr_heart - qr_bone - qr_fat - qr_tumor
    qr_livout <- qr_liver + qr_spleen + qr_pancreas + qr_gi

    # ---- Permeability-surface products, PS_i = k_i * V_i,total (L/min) ----
    psp_tumor    <- kps_tumor_p * vtot_tumor
    psp_liver    <- 500 * kps_muscle_p * vtot_liver
    psp_spleen   <- 500 * kps_muscle_p * vtot_spleen
    psp_pancreas <- 100 * kps_muscle_p * vtot_pancreas
    psp_gi       <- kps_muscle_p * vtot_gi
    psp_adrenal  <- 100 * kps_muscle_p * vtot_adrenal
    psp_lung     <- 100 * kps_muscle_p * vtot_lung
    psp_muscle   <- kps_muscle_p * vtot_muscle
    psp_skin     <- kps_muscle_p * vtot_skin
    psp_heart    <- kps_muscle_p * vtot_heart
    psp_bone     <- kps_muscle_p * vtot_bone
    psp_fat      <- kps_muscle_p * vtot_fat
    psp_rb       <- kps_muscle_p * vtot_rb
    # Table S1 gives no free-nuclide value for the tumour; every tabulated
    # tissue has k_i,f / k_i,p = k_Mus,f / k_Mus,p, which is applied to k_T.
    psf_tumor    <- kps_tumor_p * (kps_muscle_f / kps_muscle_p) * vtot_tumor
    psf_liver    <- 500 * kps_muscle_f * vtot_liver
    psf_spleen   <- 500 * kps_muscle_f * vtot_spleen
    psf_pancreas <- 100 * kps_muscle_f * vtot_pancreas
    psf_gi       <- kps_muscle_f * vtot_gi
    psf_adrenal  <- 100 * kps_muscle_f * vtot_adrenal
    psf_lung     <- 100 * kps_muscle_f * vtot_lung
    psf_muscle   <- kps_muscle_f * vtot_muscle
    psf_skin     <- kps_muscle_f * vtot_skin
    psf_heart    <- kps_muscle_f * vtot_heart
    psf_bone     <- kps_muscle_f * vtot_bone
    psf_fat      <- kps_muscle_f * vtot_fat
    psf_rb       <- kps_muscle_f * vtot_rb

    # ---- Free sst2 sites, shared by the labelled and unlabelled systems:
    #      total minus the surface-bound and the early-endosomal peptide.
    #      Methods 2.1 has the receptor travel with the internalised complex
    #      and recycle to the surface after sorting, so early-endosomal
    #      complexes still occupy a site (Kletting 2016 Eq 1 subtracts only
    #      the surface-bound peptide; see vignette). ------------------------
    rfree_tumor    <- rdens_tumor * vtot_tumor - pb_bnd_tumor - bis_bnd_tumor - cold_bnd_tumor -
      pb_endo_tumor - bis_endo_tumor - cold_endo_tumor
    rfree_liver    <- rdens_liver * vtot_liver - pb_bnd_liver - bis_bnd_liver - cold_bnd_liver -
      pb_endo_liver - bis_endo_liver - cold_endo_liver
    rfree_spleen   <- rdens_spleen * vtot_spleen - pb_bnd_spleen - bis_bnd_spleen - cold_bnd_spleen -
      pb_endo_spleen - bis_endo_spleen - cold_endo_spleen
    rfree_pancreas <- rdens_pancreas * vtot_pancreas - pb_bnd_pancreas - bis_bnd_pancreas - cold_bnd_pancreas -
      pb_endo_pancreas - bis_endo_pancreas - cold_endo_pancreas
    rfree_gi       <- rdens_gi * vtot_gi - pb_bnd_gi - bis_bnd_gi - cold_bnd_gi -
      pb_endo_gi - bis_endo_gi - cold_endo_gi
    rfree_adrenal  <- rdens_adrenal * vtot_adrenal - pb_bnd_adrenal - bis_bnd_adrenal - cold_bnd_adrenal -
      pb_endo_adrenal - bis_endo_adrenal - cold_endo_adrenal
    rfree_lung     <- rdens_lung * vtot_lung - pb_bnd_lung - bis_bnd_lung - cold_bnd_lung -
      pb_endo_lung - bis_endo_lung - cold_endo_lung
    rfree_kidney   <- rdens_kidney * vtot_kidney - pb_bnd_kidney - bis_bnd_kidney - cold_bnd_kidney -
      pb_endo_kidney - bis_endo_kidney - cold_endo_kidney

    # =====================================================================
    # 3. ODEs (Supplement A, Eqs 1-64). State names are
    #    <species>_<sub-compartment>_<tissue>: species pb = [212Pb]Pb-peptide,
    #    bis = [212Bi]Bi-peptide, bi = free 212Bi, bihwpp / birbc = free 212Bi
    #    bound to plasma proteins / red cells, tl = free 208Tl, cold =
    #    unlabelled peptide; sub-compartments v (vascular), int
    #    (interstitial), bnd (sst2-bound), endo (early endosome), sort
    #    (sorting), cell (free-212Bi cellular), tub / nsp / exc (kidney
    #    tubule, non-specific uptake, urine). Stable 208Pb is not carried: it
    #    feeds back into nothing. Each tissue repeats the same equations;
    #    corrections to the printed equations are listed in the vignette.
    # =====================================================================
    # ---- Arterial pool: fed by the lung, drains to every perfused tissue --
    d/dt(pb_art) <- qp_lung * pb_v_lung / vv_lung - qp_tot * pb_art / vv_art -
      lam_pb * pb_art
    d/dt(bis_art) <- qp_lung * bis_v_lung / vv_lung - qp_tot * bis_art / vv_art +
      lam_gr * pb_art - lam_bi * bis_art
    d/dt(bi_art) <- qp_lung * bi_v_lung / vv_lung - qp_tot * bi_art / vv_art +
      kdiss_hwpp * bihwpp_art - kass_hwpp * bi_art +
      kdiss_rbc * birbc_art - kass_rbc * bi_art -
      lam_bi * bi_art + lam_ex * pb_art
    d/dt(bihwpp_art) <- qp_lung * bihwpp_v_lung / vv_lung - qp_tot * bihwpp_art / vv_art +
      kass_hwpp * bi_art - kdiss_hwpp * bihwpp_art - lam_bi * bihwpp_art
    d/dt(birbc_art) <- qr_lung * birbc_v_lung / vr_lung - qr_tot * birbc_art / vr_art +
      kass_rbc * bi_art - kdiss_rbc * birbc_art - lam_bi * birbc_art
    d/dt(tl_art) <- qp_lung * tl_v_lung / vv_lung - qp_tot * tl_art / vv_art +
      lam_bia * (bis_art + bi_art + bihwpp_art + birbc_art) - lam_tl * tl_art
    d/dt(cold_art) <- qp_lung * cold_v_lung / vv_lung - qp_tot * cold_art / vv_art +
      lam_ex * pb_art + lam_bi * bis_art
    # ---- Venous pool: collects every tissue outflow (portal organs via the
    #      liver), drains to the lung; the intravenous dose enters here -----
    d/dt(pb_ven) <- qp_tumor * pb_v_tumor / vv_tumor +
      qp_kidney * pb_v_kidney / vv_kidney +
      qp_livout * pb_v_liver / vv_liver +
      qp_adrenal * pb_v_adrenal / vv_adrenal +
      qp_muscle * pb_v_muscle / vv_muscle +
      qp_skin * pb_v_skin / vv_skin +
      qp_heart * pb_v_heart / vv_heart +
      qp_bone * pb_v_bone / vv_bone +
      qp_fat * pb_v_fat / vv_fat +
      qp_rb * pb_v_rb / vv_rb -
      qp_tot * pb_ven / vv_ven -
      lam_pb * pb_ven
    d/dt(bis_ven) <- qp_tumor * bis_v_tumor / vv_tumor +
      qp_kidney * bis_v_kidney / vv_kidney +
      qp_livout * bis_v_liver / vv_liver +
      qp_adrenal * bis_v_adrenal / vv_adrenal +
      qp_muscle * bis_v_muscle / vv_muscle +
      qp_skin * bis_v_skin / vv_skin +
      qp_heart * bis_v_heart / vv_heart +
      qp_bone * bis_v_bone / vv_bone +
      qp_fat * bis_v_fat / vv_fat +
      qp_rb * bis_v_rb / vv_rb -
      qp_tot * bis_ven / vv_ven +
      lam_gr * pb_ven - lam_bi * bis_ven
    d/dt(bi_ven) <- qp_tumor * bi_v_tumor / vv_tumor +
      qp_kidney * bi_v_kidney / vv_kidney +
      qp_livout * bi_v_liver / vv_liver +
      qp_adrenal * bi_v_adrenal / vv_adrenal +
      qp_muscle * bi_v_muscle / vv_muscle +
      qp_skin * bi_v_skin / vv_skin +
      qp_heart * bi_v_heart / vv_heart +
      qp_bone * bi_v_bone / vv_bone +
      qp_fat * bi_v_fat / vv_fat +
      qp_rb * bi_v_rb / vv_rb -
      qp_tot * bi_ven / vv_ven +
      kdiss_hwpp * bihwpp_ven - kass_hwpp * bi_ven +
      kdiss_rbc * birbc_ven - kass_rbc * bi_ven -
      lam_bi * bi_ven + lam_ex * pb_ven
    d/dt(bihwpp_ven) <- qp_tumor * bihwpp_v_tumor / vv_tumor +
      qp_kidney * bihwpp_v_kidney / vv_kidney +
      qp_livout * bihwpp_v_liver / vv_liver +
      qp_adrenal * bihwpp_v_adrenal / vv_adrenal +
      qp_muscle * bihwpp_v_muscle / vv_muscle +
      qp_skin * bihwpp_v_skin / vv_skin +
      qp_heart * bihwpp_v_heart / vv_heart +
      qp_bone * bihwpp_v_bone / vv_bone +
      qp_fat * bihwpp_v_fat / vv_fat +
      qp_rb * bihwpp_v_rb / vv_rb -
      qp_tot * bihwpp_ven / vv_ven +
      kass_hwpp * bi_ven - kdiss_hwpp * bihwpp_ven - lam_bi * bihwpp_ven
    d/dt(birbc_ven) <- qr_tumor * birbc_v_tumor / vr_tumor +
      qr_kidney * birbc_v_kidney / vr_kidney +
      qr_livout * birbc_v_liver / vr_liver +
      qr_adrenal * birbc_v_adrenal / vr_adrenal +
      qr_muscle * birbc_v_muscle / vr_muscle +
      qr_skin * birbc_v_skin / vr_skin +
      qr_heart * birbc_v_heart / vr_heart +
      qr_bone * birbc_v_bone / vr_bone +
      qr_fat * birbc_v_fat / vr_fat +
      qr_rb * birbc_v_rb / vr_rb -
      qr_tot * birbc_ven / vr_ven +
      kass_rbc * bi_ven - kdiss_rbc * birbc_ven - lam_bi * birbc_ven
    d/dt(tl_ven) <- qp_tumor * tl_v_tumor / vv_tumor +
      qp_kidney * tl_v_kidney / vv_kidney +
      qp_livout * tl_v_liver / vv_liver +
      qp_adrenal * tl_v_adrenal / vv_adrenal +
      qp_muscle * tl_v_muscle / vv_muscle +
      qp_skin * tl_v_skin / vv_skin +
      qp_heart * tl_v_heart / vv_heart +
      qp_bone * tl_v_bone / vv_bone +
      qp_fat * tl_v_fat / vv_fat +
      qp_rb * tl_v_rb / vv_rb -
      qp_tot * tl_ven / vv_ven +
      lam_bia * (bis_ven + bi_ven + bihwpp_ven + birbc_ven) - lam_tl * tl_ven
    d/dt(cold_ven) <- qp_tumor * cold_v_tumor / vv_tumor +
      qp_kidney * cold_v_kidney / vv_kidney +
      qp_livout * cold_v_liver / vv_liver +
      qp_adrenal * cold_v_adrenal / vv_adrenal +
      qp_muscle * cold_v_muscle / vv_muscle +
      qp_skin * cold_v_skin / vv_skin +
      qp_heart * cold_v_heart / vv_heart +
      qp_bone * cold_v_bone / vv_bone +
      qp_fat * cold_v_fat / vv_fat +
      qp_rb * cold_v_rb / vv_rb -
      qp_tot * cold_ven / vv_ven +
      lam_ex * pb_ven + lam_bi * bis_ven
    # ---- tumor (sst2-expressing) ----
    d/dt(pb_v_tumor) <- qp_tumor * pb_art / vv_art -
      qp_tumor * pb_v_tumor / vv_tumor +
      psp_tumor * (pb_int_tumor / vi_tumor - pb_v_tumor / vv_tumor) -
      lam_pb * pb_v_tumor
    d/dt(bis_v_tumor) <- qp_tumor * bis_art / vv_art -
      qp_tumor * bis_v_tumor / vv_tumor +
      psp_tumor * (bis_int_tumor / vi_tumor - bis_v_tumor / vv_tumor) +
      lam_gr * pb_v_tumor - lam_bi * bis_v_tumor
    d/dt(bi_v_tumor) <- qp_tumor * bi_art / vv_art -
      qp_tumor * bi_v_tumor / vv_tumor +
      psf_tumor * (bi_int_tumor / vi_tumor - bi_v_tumor / vv_tumor) +
      kdiss_hwpp * bihwpp_v_tumor - kass_hwpp * bi_v_tumor +
      kdiss_rbc * birbc_v_tumor - kass_rbc * bi_v_tumor -
      lam_bi * bi_v_tumor + lam_ex * pb_v_tumor
    d/dt(bihwpp_v_tumor) <- qp_tumor * bihwpp_art / vv_art -
      qp_tumor * bihwpp_v_tumor / vv_tumor +
      kass_hwpp * bi_v_tumor - kdiss_hwpp * bihwpp_v_tumor - lam_bi * bihwpp_v_tumor
    d/dt(birbc_v_tumor) <- qr_tumor * birbc_art / vr_art -
      qr_tumor * birbc_v_tumor / vr_tumor +
      kass_rbc * bi_v_tumor - kdiss_rbc * birbc_v_tumor - lam_bi * birbc_v_tumor
    d/dt(tl_v_tumor) <- qp_tumor * tl_art / vv_art -
      qp_tumor * tl_v_tumor / vv_tumor +
      psf_tumor * (tl_int_tumor / vi_tumor - tl_v_tumor / vv_tumor) +
      lam_bia * (bis_v_tumor + bi_v_tumor + bihwpp_v_tumor + birbc_v_tumor) - lam_tl * tl_v_tumor
    d/dt(cold_v_tumor) <- qp_tumor * cold_art / vv_art -
      qp_tumor * cold_v_tumor / vv_tumor +
      psp_tumor * (cold_int_tumor / vi_tumor - cold_v_tumor / vv_tumor) +
      lam_ex * pb_v_tumor + lam_bi * bis_v_tumor
    d/dt(pb_int_tumor) <- psp_tumor * (pb_v_tumor / vv_tumor - pb_int_tumor / vi_tumor) +
      krec * pb_endo_tumor + koff * pb_bnd_tumor -
      kon / vi_tumor * rfree_tumor * pb_int_tumor -
      lam_pb * pb_int_tumor
    d/dt(bis_int_tumor) <- psp_tumor * (bis_v_tumor / vv_tumor - bis_int_tumor / vi_tumor) +
      krec * bis_endo_tumor + koff * bis_bnd_tumor -
      kon / vi_tumor * rfree_tumor * bis_int_tumor +
      lam_gr * pb_int_tumor - lam_bi * bis_int_tumor
    d/dt(bi_int_tumor) <- psf_tumor * (bi_v_tumor / vv_tumor - bi_int_tumor / vi_tumor) -
      kup_tumor * bi_int_tumor - lam_bi * bi_int_tumor + lam_ex * (pb_int_tumor + pb_bnd_tumor)
    d/dt(tl_int_tumor) <- psf_tumor * (tl_v_tumor / vv_tumor - tl_int_tumor / vi_tumor) +
      lam_bia * (bis_int_tumor + bi_int_tumor + bis_bnd_tumor) - lam_tl * tl_int_tumor
    d/dt(cold_int_tumor) <- psp_tumor * (cold_v_tumor / vv_tumor - cold_int_tumor / vi_tumor) +
      krec * cold_endo_tumor + koff * cold_bnd_tumor -
      kon / vi_tumor * rfree_tumor * cold_int_tumor +
      lam_ex * pb_int_tumor + lam_bi * bis_int_tumor
    d/dt(pb_bnd_tumor) <- -koff * pb_bnd_tumor + kon / vi_tumor * rfree_tumor * pb_int_tumor -
      kint * pb_bnd_tumor - lam_pb * pb_bnd_tumor
    d/dt(bis_bnd_tumor) <- -koff * bis_bnd_tumor + kon / vi_tumor * rfree_tumor * bis_int_tumor -
      kint * bis_bnd_tumor +
      lam_gr * pb_bnd_tumor - lam_bi * bis_bnd_tumor
    d/dt(cold_bnd_tumor) <- -koff * cold_bnd_tumor + kon / vi_tumor * rfree_tumor * cold_int_tumor -
      kint * cold_bnd_tumor +
      lam_ex * pb_bnd_tumor + lam_bi * bis_bnd_tumor
    d/dt(pb_endo_tumor) <- kint * pb_bnd_tumor - krec * pb_endo_tumor - ksort * pb_endo_tumor -
      lam_pb * pb_endo_tumor
    d/dt(bis_endo_tumor) <- kint * bis_bnd_tumor - krec * bis_endo_tumor - ksort * bis_endo_tumor +
      lam_gr * pb_endo_tumor - lam_bi * bis_endo_tumor
    d/dt(bi_endo_tumor) <- kup_tumor * bi_int_tumor + lam_ex * pb_endo_tumor - lam_bi * bi_endo_tumor
    d/dt(tl_endo_tumor) <- lam_bia * (bis_endo_tumor + bi_endo_tumor) - lam_tl * tl_endo_tumor
    d/dt(cold_endo_tumor) <- kint * cold_bnd_tumor - krec * cold_endo_tumor - ksort * cold_endo_tumor +
      lam_ex * pb_endo_tumor + lam_bi * bis_endo_tumor
    d/dt(pb_sort_tumor) <- ksort * pb_endo_tumor - krel_tumor * pb_sort_tumor - lam_pb * pb_sort_tumor
    d/dt(bis_sort_tumor) <- ksort * bis_endo_tumor - krel_tumor * bis_sort_tumor +
      lam_gr * pb_sort_tumor - lam_bi * bis_sort_tumor
    d/dt(bi_sort_tumor) <- lam_ex * pb_sort_tumor - lam_bi * bi_sort_tumor
    d/dt(tl_sort_tumor) <- lam_bia * (bis_sort_tumor + bi_sort_tumor) - lam_tl * tl_sort_tumor
    d/dt(cold_sort_tumor) <- ksort * cold_endo_tumor - krel_tumor * cold_sort_tumor +
      lam_ex * pb_sort_tumor + lam_bi * bis_sort_tumor
    # ---- liver (sst2-expressing) ----
    d/dt(pb_v_liver) <- qp_liver * pb_art / vv_art +
      qp_spleen * pb_v_spleen / vv_spleen +
      qp_pancreas * pb_v_pancreas / vv_pancreas +
      qp_gi * pb_v_gi / vv_gi -
      qp_livout * pb_v_liver / vv_liver +
      psp_liver * (pb_int_liver / vi_liver - pb_v_liver / vv_liver) -
      lam_pb * pb_v_liver
    d/dt(bis_v_liver) <- qp_liver * bis_art / vv_art +
      qp_spleen * bis_v_spleen / vv_spleen +
      qp_pancreas * bis_v_pancreas / vv_pancreas +
      qp_gi * bis_v_gi / vv_gi -
      qp_livout * bis_v_liver / vv_liver +
      psp_liver * (bis_int_liver / vi_liver - bis_v_liver / vv_liver) +
      lam_gr * pb_v_liver - lam_bi * bis_v_liver
    d/dt(bi_v_liver) <- qp_liver * bi_art / vv_art +
      qp_spleen * bi_v_spleen / vv_spleen +
      qp_pancreas * bi_v_pancreas / vv_pancreas +
      qp_gi * bi_v_gi / vv_gi -
      qp_livout * bi_v_liver / vv_liver +
      psf_liver * (bi_int_liver / vi_liver - bi_v_liver / vv_liver) +
      kdiss_hwpp * bihwpp_v_liver - kass_hwpp * bi_v_liver +
      kdiss_rbc * birbc_v_liver - kass_rbc * bi_v_liver -
      lam_bi * bi_v_liver + lam_ex * pb_v_liver
    d/dt(bihwpp_v_liver) <- qp_liver * bihwpp_art / vv_art +
      qp_spleen * bihwpp_v_spleen / vv_spleen +
      qp_pancreas * bihwpp_v_pancreas / vv_pancreas +
      qp_gi * bihwpp_v_gi / vv_gi -
      qp_livout * bihwpp_v_liver / vv_liver +
      kass_hwpp * bi_v_liver - kdiss_hwpp * bihwpp_v_liver - lam_bi * bihwpp_v_liver
    d/dt(birbc_v_liver) <- qr_liver * birbc_art / vr_art +
      qr_spleen * birbc_v_spleen / vr_spleen +
      qr_pancreas * birbc_v_pancreas / vr_pancreas +
      qr_gi * birbc_v_gi / vr_gi -
      qr_livout * birbc_v_liver / vr_liver +
      kass_rbc * bi_v_liver - kdiss_rbc * birbc_v_liver - lam_bi * birbc_v_liver
    d/dt(tl_v_liver) <- qp_liver * tl_art / vv_art +
      qp_spleen * tl_v_spleen / vv_spleen +
      qp_pancreas * tl_v_pancreas / vv_pancreas +
      qp_gi * tl_v_gi / vv_gi -
      qp_livout * tl_v_liver / vv_liver +
      psf_liver * (tl_int_liver / vi_liver - tl_v_liver / vv_liver) +
      lam_bia * (bis_v_liver + bi_v_liver + bihwpp_v_liver + birbc_v_liver) - lam_tl * tl_v_liver
    d/dt(cold_v_liver) <- qp_liver * cold_art / vv_art +
      qp_spleen * cold_v_spleen / vv_spleen +
      qp_pancreas * cold_v_pancreas / vv_pancreas +
      qp_gi * cold_v_gi / vv_gi -
      qp_livout * cold_v_liver / vv_liver +
      psp_liver * (cold_int_liver / vi_liver - cold_v_liver / vv_liver) +
      lam_ex * pb_v_liver + lam_bi * bis_v_liver
    d/dt(pb_int_liver) <- psp_liver * (pb_v_liver / vv_liver - pb_int_liver / vi_liver) +
      krec * pb_endo_liver + koff * pb_bnd_liver -
      kon / vi_liver * rfree_liver * pb_int_liver -
      lam_pb * pb_int_liver
    d/dt(bis_int_liver) <- psp_liver * (bis_v_liver / vv_liver - bis_int_liver / vi_liver) +
      krec * bis_endo_liver + koff * bis_bnd_liver -
      kon / vi_liver * rfree_liver * bis_int_liver +
      lam_gr * pb_int_liver - lam_bi * bis_int_liver
    d/dt(bi_int_liver) <- psf_liver * (bi_v_liver / vv_liver - bi_int_liver / vi_liver) -
      kup_liver * bi_int_liver - lam_bi * bi_int_liver + lam_ex * (pb_int_liver + pb_bnd_liver)
    d/dt(tl_int_liver) <- psf_liver * (tl_v_liver / vv_liver - tl_int_liver / vi_liver) +
      lam_bia * (bis_int_liver + bi_int_liver + bis_bnd_liver) - lam_tl * tl_int_liver
    d/dt(cold_int_liver) <- psp_liver * (cold_v_liver / vv_liver - cold_int_liver / vi_liver) +
      krec * cold_endo_liver + koff * cold_bnd_liver -
      kon / vi_liver * rfree_liver * cold_int_liver +
      lam_ex * pb_int_liver + lam_bi * bis_int_liver
    d/dt(pb_bnd_liver) <- -koff * pb_bnd_liver + kon / vi_liver * rfree_liver * pb_int_liver -
      kint * pb_bnd_liver - lam_pb * pb_bnd_liver
    d/dt(bis_bnd_liver) <- -koff * bis_bnd_liver + kon / vi_liver * rfree_liver * bis_int_liver -
      kint * bis_bnd_liver +
      lam_gr * pb_bnd_liver - lam_bi * bis_bnd_liver
    d/dt(cold_bnd_liver) <- -koff * cold_bnd_liver + kon / vi_liver * rfree_liver * cold_int_liver -
      kint * cold_bnd_liver +
      lam_ex * pb_bnd_liver + lam_bi * bis_bnd_liver
    d/dt(pb_endo_liver) <- kint * pb_bnd_liver - krec * pb_endo_liver - ksort * pb_endo_liver -
      lam_pb * pb_endo_liver
    d/dt(bis_endo_liver) <- kint * bis_bnd_liver - krec * bis_endo_liver - ksort * bis_endo_liver +
      lam_gr * pb_endo_liver - lam_bi * bis_endo_liver
    d/dt(bi_endo_liver) <- kup_liver * bi_int_liver + lam_ex * pb_endo_liver - lam_bi * bi_endo_liver
    d/dt(tl_endo_liver) <- lam_bia * (bis_endo_liver + bi_endo_liver) - lam_tl * tl_endo_liver
    d/dt(cold_endo_liver) <- kint * cold_bnd_liver - krec * cold_endo_liver - ksort * cold_endo_liver +
      lam_ex * pb_endo_liver + lam_bi * bis_endo_liver
    d/dt(pb_sort_liver) <- ksort * pb_endo_liver - krel_liver * pb_sort_liver - lam_pb * pb_sort_liver
    d/dt(bis_sort_liver) <- ksort * bis_endo_liver - krel_liver * bis_sort_liver +
      lam_gr * pb_sort_liver - lam_bi * bis_sort_liver
    d/dt(bi_sort_liver) <- lam_ex * pb_sort_liver - lam_bi * bi_sort_liver
    d/dt(tl_sort_liver) <- lam_bia * (bis_sort_liver + bi_sort_liver) - lam_tl * tl_sort_liver
    d/dt(cold_sort_liver) <- ksort * cold_endo_liver - krel_liver * cold_sort_liver +
      lam_ex * pb_sort_liver + lam_bi * bis_sort_liver
    # ---- spleen (sst2-expressing) ----
    d/dt(pb_v_spleen) <- qp_spleen * pb_art / vv_art -
      qp_spleen * pb_v_spleen / vv_spleen +
      psp_spleen * (pb_int_spleen / vi_spleen - pb_v_spleen / vv_spleen) -
      lam_pb * pb_v_spleen
    d/dt(bis_v_spleen) <- qp_spleen * bis_art / vv_art -
      qp_spleen * bis_v_spleen / vv_spleen +
      psp_spleen * (bis_int_spleen / vi_spleen - bis_v_spleen / vv_spleen) +
      lam_gr * pb_v_spleen - lam_bi * bis_v_spleen
    d/dt(bi_v_spleen) <- qp_spleen * bi_art / vv_art -
      qp_spleen * bi_v_spleen / vv_spleen +
      psf_spleen * (bi_int_spleen / vi_spleen - bi_v_spleen / vv_spleen) +
      kdiss_hwpp * bihwpp_v_spleen - kass_hwpp * bi_v_spleen +
      kdiss_rbc * birbc_v_spleen - kass_rbc * bi_v_spleen -
      lam_bi * bi_v_spleen + lam_ex * pb_v_spleen
    d/dt(bihwpp_v_spleen) <- qp_spleen * bihwpp_art / vv_art -
      qp_spleen * bihwpp_v_spleen / vv_spleen +
      kass_hwpp * bi_v_spleen - kdiss_hwpp * bihwpp_v_spleen - lam_bi * bihwpp_v_spleen
    d/dt(birbc_v_spleen) <- qr_spleen * birbc_art / vr_art -
      qr_spleen * birbc_v_spleen / vr_spleen +
      kass_rbc * bi_v_spleen - kdiss_rbc * birbc_v_spleen - lam_bi * birbc_v_spleen
    d/dt(tl_v_spleen) <- qp_spleen * tl_art / vv_art -
      qp_spleen * tl_v_spleen / vv_spleen +
      psf_spleen * (tl_int_spleen / vi_spleen - tl_v_spleen / vv_spleen) +
      lam_bia * (bis_v_spleen + bi_v_spleen + bihwpp_v_spleen + birbc_v_spleen) - lam_tl * tl_v_spleen
    d/dt(cold_v_spleen) <- qp_spleen * cold_art / vv_art -
      qp_spleen * cold_v_spleen / vv_spleen +
      psp_spleen * (cold_int_spleen / vi_spleen - cold_v_spleen / vv_spleen) +
      lam_ex * pb_v_spleen + lam_bi * bis_v_spleen
    d/dt(pb_int_spleen) <- psp_spleen * (pb_v_spleen / vv_spleen - pb_int_spleen / vi_spleen) +
      krec * pb_endo_spleen + koff * pb_bnd_spleen -
      kon / vi_spleen * rfree_spleen * pb_int_spleen -
      lam_pb * pb_int_spleen
    d/dt(bis_int_spleen) <- psp_spleen * (bis_v_spleen / vv_spleen - bis_int_spleen / vi_spleen) +
      krec * bis_endo_spleen + koff * bis_bnd_spleen -
      kon / vi_spleen * rfree_spleen * bis_int_spleen +
      lam_gr * pb_int_spleen - lam_bi * bis_int_spleen
    d/dt(bi_int_spleen) <- psf_spleen * (bi_v_spleen / vv_spleen - bi_int_spleen / vi_spleen) -
      kup_spleen * bi_int_spleen - lam_bi * bi_int_spleen + lam_ex * (pb_int_spleen + pb_bnd_spleen)
    d/dt(tl_int_spleen) <- psf_spleen * (tl_v_spleen / vv_spleen - tl_int_spleen / vi_spleen) +
      lam_bia * (bis_int_spleen + bi_int_spleen + bis_bnd_spleen) - lam_tl * tl_int_spleen
    d/dt(cold_int_spleen) <- psp_spleen * (cold_v_spleen / vv_spleen - cold_int_spleen / vi_spleen) +
      krec * cold_endo_spleen + koff * cold_bnd_spleen -
      kon / vi_spleen * rfree_spleen * cold_int_spleen +
      lam_ex * pb_int_spleen + lam_bi * bis_int_spleen
    d/dt(pb_bnd_spleen) <- -koff * pb_bnd_spleen + kon / vi_spleen * rfree_spleen * pb_int_spleen -
      kint * pb_bnd_spleen - lam_pb * pb_bnd_spleen
    d/dt(bis_bnd_spleen) <- -koff * bis_bnd_spleen + kon / vi_spleen * rfree_spleen * bis_int_spleen -
      kint * bis_bnd_spleen +
      lam_gr * pb_bnd_spleen - lam_bi * bis_bnd_spleen
    d/dt(cold_bnd_spleen) <- -koff * cold_bnd_spleen + kon / vi_spleen * rfree_spleen * cold_int_spleen -
      kint * cold_bnd_spleen +
      lam_ex * pb_bnd_spleen + lam_bi * bis_bnd_spleen
    d/dt(pb_endo_spleen) <- kint * pb_bnd_spleen - krec * pb_endo_spleen - ksort * pb_endo_spleen -
      lam_pb * pb_endo_spleen
    d/dt(bis_endo_spleen) <- kint * bis_bnd_spleen - krec * bis_endo_spleen - ksort * bis_endo_spleen +
      lam_gr * pb_endo_spleen - lam_bi * bis_endo_spleen
    d/dt(bi_endo_spleen) <- kup_spleen * bi_int_spleen + lam_ex * pb_endo_spleen - lam_bi * bi_endo_spleen
    d/dt(tl_endo_spleen) <- lam_bia * (bis_endo_spleen + bi_endo_spleen) - lam_tl * tl_endo_spleen
    d/dt(cold_endo_spleen) <- kint * cold_bnd_spleen - krec * cold_endo_spleen - ksort * cold_endo_spleen +
      lam_ex * pb_endo_spleen + lam_bi * bis_endo_spleen
    d/dt(pb_sort_spleen) <- ksort * pb_endo_spleen - krel_spleen * pb_sort_spleen - lam_pb * pb_sort_spleen
    d/dt(bis_sort_spleen) <- ksort * bis_endo_spleen - krel_spleen * bis_sort_spleen +
      lam_gr * pb_sort_spleen - lam_bi * bis_sort_spleen
    d/dt(bi_sort_spleen) <- lam_ex * pb_sort_spleen - lam_bi * bi_sort_spleen
    d/dt(tl_sort_spleen) <- lam_bia * (bis_sort_spleen + bi_sort_spleen) - lam_tl * tl_sort_spleen
    d/dt(cold_sort_spleen) <- ksort * cold_endo_spleen - krel_spleen * cold_sort_spleen +
      lam_ex * pb_sort_spleen + lam_bi * bis_sort_spleen
    # ---- pancreas (sst2-expressing) ----
    d/dt(pb_v_pancreas) <- qp_pancreas * pb_art / vv_art -
      qp_pancreas * pb_v_pancreas / vv_pancreas +
      psp_pancreas * (pb_int_pancreas / vi_pancreas - pb_v_pancreas / vv_pancreas) -
      lam_pb * pb_v_pancreas
    d/dt(bis_v_pancreas) <- qp_pancreas * bis_art / vv_art -
      qp_pancreas * bis_v_pancreas / vv_pancreas +
      psp_pancreas * (bis_int_pancreas / vi_pancreas - bis_v_pancreas / vv_pancreas) +
      lam_gr * pb_v_pancreas - lam_bi * bis_v_pancreas
    d/dt(bi_v_pancreas) <- qp_pancreas * bi_art / vv_art -
      qp_pancreas * bi_v_pancreas / vv_pancreas +
      psf_pancreas * (bi_int_pancreas / vi_pancreas - bi_v_pancreas / vv_pancreas) +
      kdiss_hwpp * bihwpp_v_pancreas - kass_hwpp * bi_v_pancreas +
      kdiss_rbc * birbc_v_pancreas - kass_rbc * bi_v_pancreas -
      lam_bi * bi_v_pancreas + lam_ex * pb_v_pancreas
    d/dt(bihwpp_v_pancreas) <- qp_pancreas * bihwpp_art / vv_art -
      qp_pancreas * bihwpp_v_pancreas / vv_pancreas +
      kass_hwpp * bi_v_pancreas - kdiss_hwpp * bihwpp_v_pancreas - lam_bi * bihwpp_v_pancreas
    d/dt(birbc_v_pancreas) <- qr_pancreas * birbc_art / vr_art -
      qr_pancreas * birbc_v_pancreas / vr_pancreas +
      kass_rbc * bi_v_pancreas - kdiss_rbc * birbc_v_pancreas - lam_bi * birbc_v_pancreas
    d/dt(tl_v_pancreas) <- qp_pancreas * tl_art / vv_art -
      qp_pancreas * tl_v_pancreas / vv_pancreas +
      psf_pancreas * (tl_int_pancreas / vi_pancreas - tl_v_pancreas / vv_pancreas) +
      lam_bia * (bis_v_pancreas + bi_v_pancreas + bihwpp_v_pancreas + birbc_v_pancreas) - lam_tl * tl_v_pancreas
    d/dt(cold_v_pancreas) <- qp_pancreas * cold_art / vv_art -
      qp_pancreas * cold_v_pancreas / vv_pancreas +
      psp_pancreas * (cold_int_pancreas / vi_pancreas - cold_v_pancreas / vv_pancreas) +
      lam_ex * pb_v_pancreas + lam_bi * bis_v_pancreas
    d/dt(pb_int_pancreas) <- psp_pancreas * (pb_v_pancreas / vv_pancreas - pb_int_pancreas / vi_pancreas) +
      krec * pb_endo_pancreas + koff * pb_bnd_pancreas -
      kon / vi_pancreas * rfree_pancreas * pb_int_pancreas -
      lam_pb * pb_int_pancreas
    d/dt(bis_int_pancreas) <- psp_pancreas * (bis_v_pancreas / vv_pancreas - bis_int_pancreas / vi_pancreas) +
      krec * bis_endo_pancreas + koff * bis_bnd_pancreas -
      kon / vi_pancreas * rfree_pancreas * bis_int_pancreas +
      lam_gr * pb_int_pancreas - lam_bi * bis_int_pancreas
    d/dt(bi_int_pancreas) <- psf_pancreas * (bi_v_pancreas / vv_pancreas - bi_int_pancreas / vi_pancreas) -
      kup_pancreas * bi_int_pancreas - lam_bi * bi_int_pancreas + lam_ex * (pb_int_pancreas + pb_bnd_pancreas)
    d/dt(tl_int_pancreas) <- psf_pancreas * (tl_v_pancreas / vv_pancreas - tl_int_pancreas / vi_pancreas) +
      lam_bia * (bis_int_pancreas + bi_int_pancreas + bis_bnd_pancreas) - lam_tl * tl_int_pancreas
    d/dt(cold_int_pancreas) <- psp_pancreas * (cold_v_pancreas / vv_pancreas - cold_int_pancreas / vi_pancreas) +
      krec * cold_endo_pancreas + koff * cold_bnd_pancreas -
      kon / vi_pancreas * rfree_pancreas * cold_int_pancreas +
      lam_ex * pb_int_pancreas + lam_bi * bis_int_pancreas
    d/dt(pb_bnd_pancreas) <- -koff * pb_bnd_pancreas + kon / vi_pancreas * rfree_pancreas * pb_int_pancreas -
      kint * pb_bnd_pancreas - lam_pb * pb_bnd_pancreas
    d/dt(bis_bnd_pancreas) <- -koff * bis_bnd_pancreas + kon / vi_pancreas * rfree_pancreas * bis_int_pancreas -
      kint * bis_bnd_pancreas +
      lam_gr * pb_bnd_pancreas - lam_bi * bis_bnd_pancreas
    d/dt(cold_bnd_pancreas) <- -koff * cold_bnd_pancreas + kon / vi_pancreas * rfree_pancreas * cold_int_pancreas -
      kint * cold_bnd_pancreas +
      lam_ex * pb_bnd_pancreas + lam_bi * bis_bnd_pancreas
    d/dt(pb_endo_pancreas) <- kint * pb_bnd_pancreas - krec * pb_endo_pancreas - ksort * pb_endo_pancreas -
      lam_pb * pb_endo_pancreas
    d/dt(bis_endo_pancreas) <- kint * bis_bnd_pancreas - krec * bis_endo_pancreas - ksort * bis_endo_pancreas +
      lam_gr * pb_endo_pancreas - lam_bi * bis_endo_pancreas
    d/dt(bi_endo_pancreas) <- kup_pancreas * bi_int_pancreas + lam_ex * pb_endo_pancreas - lam_bi * bi_endo_pancreas
    d/dt(tl_endo_pancreas) <- lam_bia * (bis_endo_pancreas + bi_endo_pancreas) - lam_tl * tl_endo_pancreas
    d/dt(cold_endo_pancreas) <- kint * cold_bnd_pancreas - krec * cold_endo_pancreas - ksort * cold_endo_pancreas +
      lam_ex * pb_endo_pancreas + lam_bi * bis_endo_pancreas
    d/dt(pb_sort_pancreas) <- ksort * pb_endo_pancreas - krel_pancreas * pb_sort_pancreas - lam_pb * pb_sort_pancreas
    d/dt(bis_sort_pancreas) <- ksort * bis_endo_pancreas - krel_pancreas * bis_sort_pancreas +
      lam_gr * pb_sort_pancreas - lam_bi * bis_sort_pancreas
    d/dt(bi_sort_pancreas) <- lam_ex * pb_sort_pancreas - lam_bi * bi_sort_pancreas
    d/dt(tl_sort_pancreas) <- lam_bia * (bis_sort_pancreas + bi_sort_pancreas) - lam_tl * tl_sort_pancreas
    d/dt(cold_sort_pancreas) <- ksort * cold_endo_pancreas - krel_pancreas * cold_sort_pancreas +
      lam_ex * pb_sort_pancreas + lam_bi * bis_sort_pancreas
    # ---- gi (sst2-expressing) ----
    d/dt(pb_v_gi) <- qp_gi * pb_art / vv_art -
      qp_gi * pb_v_gi / vv_gi +
      psp_gi * (pb_int_gi / vi_gi - pb_v_gi / vv_gi) -
      lam_pb * pb_v_gi
    d/dt(bis_v_gi) <- qp_gi * bis_art / vv_art -
      qp_gi * bis_v_gi / vv_gi +
      psp_gi * (bis_int_gi / vi_gi - bis_v_gi / vv_gi) +
      lam_gr * pb_v_gi - lam_bi * bis_v_gi
    d/dt(bi_v_gi) <- qp_gi * bi_art / vv_art -
      qp_gi * bi_v_gi / vv_gi +
      psf_gi * (bi_int_gi / vi_gi - bi_v_gi / vv_gi) +
      kdiss_hwpp * bihwpp_v_gi - kass_hwpp * bi_v_gi +
      kdiss_rbc * birbc_v_gi - kass_rbc * bi_v_gi -
      lam_bi * bi_v_gi + lam_ex * pb_v_gi
    d/dt(bihwpp_v_gi) <- qp_gi * bihwpp_art / vv_art -
      qp_gi * bihwpp_v_gi / vv_gi +
      kass_hwpp * bi_v_gi - kdiss_hwpp * bihwpp_v_gi - lam_bi * bihwpp_v_gi
    d/dt(birbc_v_gi) <- qr_gi * birbc_art / vr_art -
      qr_gi * birbc_v_gi / vr_gi +
      kass_rbc * bi_v_gi - kdiss_rbc * birbc_v_gi - lam_bi * birbc_v_gi
    d/dt(tl_v_gi) <- qp_gi * tl_art / vv_art -
      qp_gi * tl_v_gi / vv_gi +
      psf_gi * (tl_int_gi / vi_gi - tl_v_gi / vv_gi) +
      lam_bia * (bis_v_gi + bi_v_gi + bihwpp_v_gi + birbc_v_gi) - lam_tl * tl_v_gi
    d/dt(cold_v_gi) <- qp_gi * cold_art / vv_art -
      qp_gi * cold_v_gi / vv_gi +
      psp_gi * (cold_int_gi / vi_gi - cold_v_gi / vv_gi) +
      lam_ex * pb_v_gi + lam_bi * bis_v_gi
    d/dt(pb_int_gi) <- psp_gi * (pb_v_gi / vv_gi - pb_int_gi / vi_gi) +
      krec * pb_endo_gi + koff * pb_bnd_gi -
      kon / vi_gi * rfree_gi * pb_int_gi -
      lam_pb * pb_int_gi
    d/dt(bis_int_gi) <- psp_gi * (bis_v_gi / vv_gi - bis_int_gi / vi_gi) +
      krec * bis_endo_gi + koff * bis_bnd_gi -
      kon / vi_gi * rfree_gi * bis_int_gi +
      lam_gr * pb_int_gi - lam_bi * bis_int_gi
    d/dt(bi_int_gi) <- psf_gi * (bi_v_gi / vv_gi - bi_int_gi / vi_gi) -
      kup_gi * bi_int_gi - lam_bi * bi_int_gi + lam_ex * (pb_int_gi + pb_bnd_gi)
    d/dt(tl_int_gi) <- psf_gi * (tl_v_gi / vv_gi - tl_int_gi / vi_gi) +
      lam_bia * (bis_int_gi + bi_int_gi + bis_bnd_gi) - lam_tl * tl_int_gi
    d/dt(cold_int_gi) <- psp_gi * (cold_v_gi / vv_gi - cold_int_gi / vi_gi) +
      krec * cold_endo_gi + koff * cold_bnd_gi -
      kon / vi_gi * rfree_gi * cold_int_gi +
      lam_ex * pb_int_gi + lam_bi * bis_int_gi
    d/dt(pb_bnd_gi) <- -koff * pb_bnd_gi + kon / vi_gi * rfree_gi * pb_int_gi -
      kint * pb_bnd_gi - lam_pb * pb_bnd_gi
    d/dt(bis_bnd_gi) <- -koff * bis_bnd_gi + kon / vi_gi * rfree_gi * bis_int_gi -
      kint * bis_bnd_gi +
      lam_gr * pb_bnd_gi - lam_bi * bis_bnd_gi
    d/dt(cold_bnd_gi) <- -koff * cold_bnd_gi + kon / vi_gi * rfree_gi * cold_int_gi -
      kint * cold_bnd_gi +
      lam_ex * pb_bnd_gi + lam_bi * bis_bnd_gi
    d/dt(pb_endo_gi) <- kint * pb_bnd_gi - krec * pb_endo_gi - ksort * pb_endo_gi -
      lam_pb * pb_endo_gi
    d/dt(bis_endo_gi) <- kint * bis_bnd_gi - krec * bis_endo_gi - ksort * bis_endo_gi +
      lam_gr * pb_endo_gi - lam_bi * bis_endo_gi
    d/dt(bi_endo_gi) <- kup_gi * bi_int_gi + lam_ex * pb_endo_gi - lam_bi * bi_endo_gi
    d/dt(tl_endo_gi) <- lam_bia * (bis_endo_gi + bi_endo_gi) - lam_tl * tl_endo_gi
    d/dt(cold_endo_gi) <- kint * cold_bnd_gi - krec * cold_endo_gi - ksort * cold_endo_gi +
      lam_ex * pb_endo_gi + lam_bi * bis_endo_gi
    d/dt(pb_sort_gi) <- ksort * pb_endo_gi - krel_gi * pb_sort_gi - lam_pb * pb_sort_gi
    d/dt(bis_sort_gi) <- ksort * bis_endo_gi - krel_gi * bis_sort_gi +
      lam_gr * pb_sort_gi - lam_bi * bis_sort_gi
    d/dt(bi_sort_gi) <- lam_ex * pb_sort_gi - lam_bi * bi_sort_gi
    d/dt(tl_sort_gi) <- lam_bia * (bis_sort_gi + bi_sort_gi) - lam_tl * tl_sort_gi
    d/dt(cold_sort_gi) <- ksort * cold_endo_gi - krel_gi * cold_sort_gi +
      lam_ex * pb_sort_gi + lam_bi * bis_sort_gi
    # ---- adrenal (sst2-expressing) ----
    d/dt(pb_v_adrenal) <- qp_adrenal * pb_art / vv_art -
      qp_adrenal * pb_v_adrenal / vv_adrenal +
      psp_adrenal * (pb_int_adrenal / vi_adrenal - pb_v_adrenal / vv_adrenal) -
      lam_pb * pb_v_adrenal
    d/dt(bis_v_adrenal) <- qp_adrenal * bis_art / vv_art -
      qp_adrenal * bis_v_adrenal / vv_adrenal +
      psp_adrenal * (bis_int_adrenal / vi_adrenal - bis_v_adrenal / vv_adrenal) +
      lam_gr * pb_v_adrenal - lam_bi * bis_v_adrenal
    d/dt(bi_v_adrenal) <- qp_adrenal * bi_art / vv_art -
      qp_adrenal * bi_v_adrenal / vv_adrenal +
      psf_adrenal * (bi_int_adrenal / vi_adrenal - bi_v_adrenal / vv_adrenal) +
      kdiss_hwpp * bihwpp_v_adrenal - kass_hwpp * bi_v_adrenal +
      kdiss_rbc * birbc_v_adrenal - kass_rbc * bi_v_adrenal -
      lam_bi * bi_v_adrenal + lam_ex * pb_v_adrenal
    d/dt(bihwpp_v_adrenal) <- qp_adrenal * bihwpp_art / vv_art -
      qp_adrenal * bihwpp_v_adrenal / vv_adrenal +
      kass_hwpp * bi_v_adrenal - kdiss_hwpp * bihwpp_v_adrenal - lam_bi * bihwpp_v_adrenal
    d/dt(birbc_v_adrenal) <- qr_adrenal * birbc_art / vr_art -
      qr_adrenal * birbc_v_adrenal / vr_adrenal +
      kass_rbc * bi_v_adrenal - kdiss_rbc * birbc_v_adrenal - lam_bi * birbc_v_adrenal
    d/dt(tl_v_adrenal) <- qp_adrenal * tl_art / vv_art -
      qp_adrenal * tl_v_adrenal / vv_adrenal +
      psf_adrenal * (tl_int_adrenal / vi_adrenal - tl_v_adrenal / vv_adrenal) +
      lam_bia * (bis_v_adrenal + bi_v_adrenal + bihwpp_v_adrenal + birbc_v_adrenal) - lam_tl * tl_v_adrenal
    d/dt(cold_v_adrenal) <- qp_adrenal * cold_art / vv_art -
      qp_adrenal * cold_v_adrenal / vv_adrenal +
      psp_adrenal * (cold_int_adrenal / vi_adrenal - cold_v_adrenal / vv_adrenal) +
      lam_ex * pb_v_adrenal + lam_bi * bis_v_adrenal
    d/dt(pb_int_adrenal) <- psp_adrenal * (pb_v_adrenal / vv_adrenal - pb_int_adrenal / vi_adrenal) +
      krec * pb_endo_adrenal + koff * pb_bnd_adrenal -
      kon / vi_adrenal * rfree_adrenal * pb_int_adrenal -
      lam_pb * pb_int_adrenal
    d/dt(bis_int_adrenal) <- psp_adrenal * (bis_v_adrenal / vv_adrenal - bis_int_adrenal / vi_adrenal) +
      krec * bis_endo_adrenal + koff * bis_bnd_adrenal -
      kon / vi_adrenal * rfree_adrenal * bis_int_adrenal +
      lam_gr * pb_int_adrenal - lam_bi * bis_int_adrenal
    d/dt(bi_int_adrenal) <- psf_adrenal * (bi_v_adrenal / vv_adrenal - bi_int_adrenal / vi_adrenal) -
      kup_adrenal * bi_int_adrenal - lam_bi * bi_int_adrenal + lam_ex * (pb_int_adrenal + pb_bnd_adrenal)
    d/dt(tl_int_adrenal) <- psf_adrenal * (tl_v_adrenal / vv_adrenal - tl_int_adrenal / vi_adrenal) +
      lam_bia * (bis_int_adrenal + bi_int_adrenal + bis_bnd_adrenal) - lam_tl * tl_int_adrenal
    d/dt(cold_int_adrenal) <- psp_adrenal * (cold_v_adrenal / vv_adrenal - cold_int_adrenal / vi_adrenal) +
      krec * cold_endo_adrenal + koff * cold_bnd_adrenal -
      kon / vi_adrenal * rfree_adrenal * cold_int_adrenal +
      lam_ex * pb_int_adrenal + lam_bi * bis_int_adrenal
    d/dt(pb_bnd_adrenal) <- -koff * pb_bnd_adrenal + kon / vi_adrenal * rfree_adrenal * pb_int_adrenal -
      kint * pb_bnd_adrenal - lam_pb * pb_bnd_adrenal
    d/dt(bis_bnd_adrenal) <- -koff * bis_bnd_adrenal + kon / vi_adrenal * rfree_adrenal * bis_int_adrenal -
      kint * bis_bnd_adrenal +
      lam_gr * pb_bnd_adrenal - lam_bi * bis_bnd_adrenal
    d/dt(cold_bnd_adrenal) <- -koff * cold_bnd_adrenal + kon / vi_adrenal * rfree_adrenal * cold_int_adrenal -
      kint * cold_bnd_adrenal +
      lam_ex * pb_bnd_adrenal + lam_bi * bis_bnd_adrenal
    d/dt(pb_endo_adrenal) <- kint * pb_bnd_adrenal - krec * pb_endo_adrenal - ksort * pb_endo_adrenal -
      lam_pb * pb_endo_adrenal
    d/dt(bis_endo_adrenal) <- kint * bis_bnd_adrenal - krec * bis_endo_adrenal - ksort * bis_endo_adrenal +
      lam_gr * pb_endo_adrenal - lam_bi * bis_endo_adrenal
    d/dt(bi_endo_adrenal) <- kup_adrenal * bi_int_adrenal + lam_ex * pb_endo_adrenal - lam_bi * bi_endo_adrenal
    d/dt(tl_endo_adrenal) <- lam_bia * (bis_endo_adrenal + bi_endo_adrenal) - lam_tl * tl_endo_adrenal
    d/dt(cold_endo_adrenal) <- kint * cold_bnd_adrenal - krec * cold_endo_adrenal - ksort * cold_endo_adrenal +
      lam_ex * pb_endo_adrenal + lam_bi * bis_endo_adrenal
    d/dt(pb_sort_adrenal) <- ksort * pb_endo_adrenal - krel_adrenal * pb_sort_adrenal - lam_pb * pb_sort_adrenal
    d/dt(bis_sort_adrenal) <- ksort * bis_endo_adrenal - krel_adrenal * bis_sort_adrenal +
      lam_gr * pb_sort_adrenal - lam_bi * bis_sort_adrenal
    d/dt(bi_sort_adrenal) <- lam_ex * pb_sort_adrenal - lam_bi * bi_sort_adrenal
    d/dt(tl_sort_adrenal) <- lam_bia * (bis_sort_adrenal + bi_sort_adrenal) - lam_tl * tl_sort_adrenal
    d/dt(cold_sort_adrenal) <- ksort * cold_endo_adrenal - krel_adrenal * cold_sort_adrenal +
      lam_ex * pb_sort_adrenal + lam_bi * bis_sort_adrenal
    # ---- lung (sst2-expressing) ----
    d/dt(pb_v_lung) <- qp_lung * pb_ven / vv_ven -
      qp_lung * pb_v_lung / vv_lung +
      psp_lung * (pb_int_lung / vi_lung - pb_v_lung / vv_lung) -
      lam_pb * pb_v_lung
    d/dt(bis_v_lung) <- qp_lung * bis_ven / vv_ven -
      qp_lung * bis_v_lung / vv_lung +
      psp_lung * (bis_int_lung / vi_lung - bis_v_lung / vv_lung) +
      lam_gr * pb_v_lung - lam_bi * bis_v_lung
    d/dt(bi_v_lung) <- qp_lung * bi_ven / vv_ven -
      qp_lung * bi_v_lung / vv_lung +
      psf_lung * (bi_int_lung / vi_lung - bi_v_lung / vv_lung) +
      kdiss_hwpp * bihwpp_v_lung - kass_hwpp * bi_v_lung +
      kdiss_rbc * birbc_v_lung - kass_rbc * bi_v_lung -
      lam_bi * bi_v_lung + lam_ex * pb_v_lung
    d/dt(bihwpp_v_lung) <- qp_lung * bihwpp_ven / vv_ven -
      qp_lung * bihwpp_v_lung / vv_lung +
      kass_hwpp * bi_v_lung - kdiss_hwpp * bihwpp_v_lung - lam_bi * bihwpp_v_lung
    d/dt(birbc_v_lung) <- qr_lung * birbc_ven / vr_ven -
      qr_lung * birbc_v_lung / vr_lung +
      kass_rbc * bi_v_lung - kdiss_rbc * birbc_v_lung - lam_bi * birbc_v_lung
    d/dt(tl_v_lung) <- qp_lung * tl_ven / vv_ven -
      qp_lung * tl_v_lung / vv_lung +
      psf_lung * (tl_int_lung / vi_lung - tl_v_lung / vv_lung) +
      lam_bia * (bis_v_lung + bi_v_lung + bihwpp_v_lung + birbc_v_lung) - lam_tl * tl_v_lung
    d/dt(cold_v_lung) <- qp_lung * cold_ven / vv_ven -
      qp_lung * cold_v_lung / vv_lung +
      psp_lung * (cold_int_lung / vi_lung - cold_v_lung / vv_lung) +
      lam_ex * pb_v_lung + lam_bi * bis_v_lung
    d/dt(pb_int_lung) <- psp_lung * (pb_v_lung / vv_lung - pb_int_lung / vi_lung) +
      krec * pb_endo_lung + koff * pb_bnd_lung -
      kon / vi_lung * rfree_lung * pb_int_lung -
      lam_pb * pb_int_lung
    d/dt(bis_int_lung) <- psp_lung * (bis_v_lung / vv_lung - bis_int_lung / vi_lung) +
      krec * bis_endo_lung + koff * bis_bnd_lung -
      kon / vi_lung * rfree_lung * bis_int_lung +
      lam_gr * pb_int_lung - lam_bi * bis_int_lung
    d/dt(bi_int_lung) <- psf_lung * (bi_v_lung / vv_lung - bi_int_lung / vi_lung) -
      kup_lung * bi_int_lung - lam_bi * bi_int_lung + lam_ex * (pb_int_lung + pb_bnd_lung)
    d/dt(tl_int_lung) <- psf_lung * (tl_v_lung / vv_lung - tl_int_lung / vi_lung) +
      lam_bia * (bis_int_lung + bi_int_lung + bis_bnd_lung) - lam_tl * tl_int_lung
    d/dt(cold_int_lung) <- psp_lung * (cold_v_lung / vv_lung - cold_int_lung / vi_lung) +
      krec * cold_endo_lung + koff * cold_bnd_lung -
      kon / vi_lung * rfree_lung * cold_int_lung +
      lam_ex * pb_int_lung + lam_bi * bis_int_lung
    d/dt(pb_bnd_lung) <- -koff * pb_bnd_lung + kon / vi_lung * rfree_lung * pb_int_lung -
      kint * pb_bnd_lung - lam_pb * pb_bnd_lung
    d/dt(bis_bnd_lung) <- -koff * bis_bnd_lung + kon / vi_lung * rfree_lung * bis_int_lung -
      kint * bis_bnd_lung +
      lam_gr * pb_bnd_lung - lam_bi * bis_bnd_lung
    d/dt(cold_bnd_lung) <- -koff * cold_bnd_lung + kon / vi_lung * rfree_lung * cold_int_lung -
      kint * cold_bnd_lung +
      lam_ex * pb_bnd_lung + lam_bi * bis_bnd_lung
    d/dt(pb_endo_lung) <- kint * pb_bnd_lung - krec * pb_endo_lung - ksort * pb_endo_lung -
      lam_pb * pb_endo_lung
    d/dt(bis_endo_lung) <- kint * bis_bnd_lung - krec * bis_endo_lung - ksort * bis_endo_lung +
      lam_gr * pb_endo_lung - lam_bi * bis_endo_lung
    d/dt(bi_endo_lung) <- kup_lung * bi_int_lung + lam_ex * pb_endo_lung - lam_bi * bi_endo_lung
    d/dt(tl_endo_lung) <- lam_bia * (bis_endo_lung + bi_endo_lung) - lam_tl * tl_endo_lung
    d/dt(cold_endo_lung) <- kint * cold_bnd_lung - krec * cold_endo_lung - ksort * cold_endo_lung +
      lam_ex * pb_endo_lung + lam_bi * bis_endo_lung
    d/dt(pb_sort_lung) <- ksort * pb_endo_lung - krel_lung * pb_sort_lung - lam_pb * pb_sort_lung
    d/dt(bis_sort_lung) <- ksort * bis_endo_lung - krel_lung * bis_sort_lung +
      lam_gr * pb_sort_lung - lam_bi * bis_sort_lung
    d/dt(bi_sort_lung) <- lam_ex * pb_sort_lung - lam_bi * bi_sort_lung
    d/dt(tl_sort_lung) <- lam_bia * (bis_sort_lung + bi_sort_lung) - lam_tl * tl_sort_lung
    d/dt(cold_sort_lung) <- ksort * cold_endo_lung - krel_lung * cold_sort_lung +
      lam_ex * pb_sort_lung + lam_bi * bis_sort_lung
    # ---- kidney: vascular, receptor-bound, endosomal and sorting ----
    d/dt(pb_v_kidney) <- qp_kidney * pb_art / vv_art -
      qp_kidney * pb_v_kidney / vv_kidney +
      krec * pb_endo_kidney + koff * pb_bnd_kidney -
      kon / vv_kidney * rfree_kidney * pb_v_kidney -
      gfr * theta_p * pb_v_kidney / vv_kidney -
      lam_pb * pb_v_kidney
    d/dt(bis_v_kidney) <- qp_kidney * bis_art / vv_art -
      qp_kidney * bis_v_kidney / vv_kidney +
      krec * bis_endo_kidney + koff * bis_bnd_kidney -
      kon / vv_kidney * rfree_kidney * bis_v_kidney -
      gfr * theta_p * bis_v_kidney / vv_kidney +
      lam_gr * pb_v_kidney - lam_bi * bis_v_kidney
    d/dt(bi_v_kidney) <- qp_kidney * bi_art / vv_art -
      qp_kidney * bi_v_kidney / vv_kidney -
      gfr * theta_f * bi_v_kidney / vv_kidney +
      lam_ex * pb_bnd_kidney +
      kdiss_hwpp * bihwpp_v_kidney - kass_hwpp * bi_v_kidney +
      kdiss_rbc * birbc_v_kidney - kass_rbc * bi_v_kidney -
      lam_bi * bi_v_kidney + lam_ex * pb_v_kidney
    d/dt(bihwpp_v_kidney) <- qp_kidney * bihwpp_art / vv_art -
      qp_kidney * bihwpp_v_kidney / vv_kidney +
      kass_hwpp * bi_v_kidney - kdiss_hwpp * bihwpp_v_kidney - lam_bi * bihwpp_v_kidney
    d/dt(birbc_v_kidney) <- qr_kidney * birbc_art / vr_art -
      qr_kidney * birbc_v_kidney / vr_kidney +
      kass_rbc * bi_v_kidney - kdiss_rbc * birbc_v_kidney - lam_bi * birbc_v_kidney
    d/dt(tl_v_kidney) <- qp_kidney * tl_art / vv_art -
      qp_kidney * tl_v_kidney / vv_kidney -
      gfr * theta_f * tl_v_kidney / vv_kidney +
      lam_bia * bis_bnd_kidney +
      lam_bia * (bis_v_kidney + bi_v_kidney + bihwpp_v_kidney + birbc_v_kidney) - lam_tl * tl_v_kidney
    d/dt(cold_v_kidney) <- qp_kidney * cold_art / vv_art -
      qp_kidney * cold_v_kidney / vv_kidney +
      krec * cold_endo_kidney + koff * cold_bnd_kidney -
      kon / vv_kidney * rfree_kidney * cold_v_kidney -
      gfr * theta_p * cold_v_kidney / vv_kidney +
      lam_ex * pb_v_kidney + lam_bi * bis_v_kidney
    d/dt(pb_bnd_kidney) <- -koff * pb_bnd_kidney + kon / vv_kidney * rfree_kidney * pb_v_kidney -
      kint * pb_bnd_kidney - lam_pb * pb_bnd_kidney
    d/dt(bis_bnd_kidney) <- -koff * bis_bnd_kidney + kon / vv_kidney * rfree_kidney * bis_v_kidney -
      kint * bis_bnd_kidney +
      lam_gr * pb_bnd_kidney - lam_bi * bis_bnd_kidney
    d/dt(cold_bnd_kidney) <- -koff * cold_bnd_kidney + kon / vv_kidney * rfree_kidney * cold_v_kidney -
      kint * cold_bnd_kidney +
      lam_ex * pb_bnd_kidney + lam_bi * bis_bnd_kidney
    d/dt(pb_endo_kidney) <- kint * pb_bnd_kidney - krec * pb_endo_kidney - ksort * pb_endo_kidney -
      lam_pb * pb_endo_kidney
    d/dt(bis_endo_kidney) <- kint * bis_bnd_kidney - krec * bis_endo_kidney - ksort * bis_endo_kidney +
      lam_gr * pb_endo_kidney - lam_bi * bis_endo_kidney
    d/dt(bi_endo_kidney) <- lam_ex * pb_endo_kidney - lam_bi * bi_endo_kidney
    d/dt(tl_endo_kidney) <- lam_bia * (bis_endo_kidney + bi_endo_kidney) - lam_tl * tl_endo_kidney
    d/dt(cold_endo_kidney) <- kint * cold_bnd_kidney - krec * cold_endo_kidney - ksort * cold_endo_kidney +
      lam_ex * pb_endo_kidney + lam_bi * bis_endo_kidney
    d/dt(pb_sort_kidney) <- ksort * pb_endo_kidney - krel_kidney * pb_sort_kidney - lam_pb * pb_sort_kidney
    d/dt(bis_sort_kidney) <- ksort * bis_endo_kidney - krel_kidney * bis_sort_kidney +
      lam_gr * pb_sort_kidney - lam_bi * bis_sort_kidney
    d/dt(bi_sort_kidney) <- lam_ex * pb_sort_kidney - lam_bi * bi_sort_kidney
    d/dt(tl_sort_kidney) <- lam_bia * (bis_sort_kidney + bi_sort_kidney) - lam_tl * tl_sort_kidney
    d/dt(cold_sort_kidney) <- ksort * cold_endo_kidney - krel_kidney * cold_sort_kidney +
      lam_ex * pb_sort_kidney + lam_bi * bis_sort_kidney
    # Proximal tubule (Eqs 51-55): filtered peptide and free nuclides;
    # a fraction f_exc leaves to the urine, the rest is taken up
    # non-specifically (peptide) or by the tubular cells (free 212Bi).
    d/dt(pb_tub_kidney) <- gfr * theta_p * pb_v_kidney / vv_kidney -
      gfr * theta_p * fexc_p * pb_tub_kidney / v_tub - kendo_kidney * pb_tub_kidney -
      lam_pb * pb_tub_kidney
    d/dt(bis_tub_kidney) <- gfr * theta_p * bis_v_kidney / vv_kidney -
      gfr * theta_p * fexc_p * bis_tub_kidney / v_tub - kendo_kidney * bis_tub_kidney +
      lam_gr * pb_tub_kidney - lam_bi * bis_tub_kidney
    d/dt(bi_tub_kidney) <- gfr * theta_f * bi_v_kidney / vv_kidney -
      gfr * theta_f * fexc_f * bi_tub_kidney / v_tub - kup_kidney * bi_tub_kidney +
      lam_ex * pb_tub_kidney - lam_bi * bi_tub_kidney
    d/dt(tl_tub_kidney) <- gfr * theta_f * tl_v_kidney / vv_kidney -
      gfr * theta_f * fexc_f * tl_tub_kidney / v_tub +
      lam_bia * (bis_tub_kidney + bi_tub_kidney) - lam_tl * tl_tub_kidney
    d/dt(cold_tub_kidney) <- gfr * theta_p * cold_v_kidney / vv_kidney -
      gfr * theta_p * fexc_p * cold_tub_kidney / v_tub - kendo_kidney * cold_tub_kidney +
      lam_ex * pb_tub_kidney + lam_bi * bis_tub_kidney
    # Renal non-specific uptake compartment (Eqs 56-59)
    d/dt(pb_nsp_kidney) <- kendo_kidney * pb_tub_kidney - lam_pb * pb_nsp_kidney
    d/dt(bis_nsp_kidney) <- kendo_kidney * bis_tub_kidney + lam_gr * pb_nsp_kidney -
      lam_bi * bis_nsp_kidney
    d/dt(bi_nsp_kidney) <- kup_kidney * bi_tub_kidney + lam_ex * pb_nsp_kidney -
      lam_bi * bi_nsp_kidney
    d/dt(tl_nsp_kidney) <- lam_bia * (bis_nsp_kidney + bi_nsp_kidney) - lam_tl * tl_nsp_kidney
    d/dt(cold_nsp_kidney) <- kendo_kidney * cold_tub_kidney + lam_ex * pb_nsp_kidney +
      lam_bi * bis_nsp_kidney
    # Urinary excretion (Eqs 60-64); cumulative amounts, decaying in place
    d/dt(pb_exc_kidney) <- gfr * theta_p * fexc_p * pb_tub_kidney / v_tub - lam_pb * pb_exc_kidney
    d/dt(bis_exc_kidney) <- gfr * theta_p * fexc_p * bis_tub_kidney / v_tub +
      lam_gr * pb_exc_kidney - lam_bi * bis_exc_kidney
    d/dt(bi_exc_kidney) <- gfr * theta_f * fexc_f * bi_tub_kidney / v_tub +
      lam_ex * pb_exc_kidney - lam_bi * bi_exc_kidney
    d/dt(tl_exc_kidney) <- gfr * theta_f * fexc_f * tl_tub_kidney / v_tub +
      lam_bia * (bis_exc_kidney + bi_exc_kidney) - lam_tl * tl_exc_kidney
    d/dt(cold_exc_kidney) <- gfr * theta_p * fexc_p * cold_tub_kidney / v_tub +
      lam_ex * pb_exc_kidney + lam_bi * bis_exc_kidney
    # ---- muscle (no sst2) ----
    d/dt(pb_v_muscle) <- qp_muscle * pb_art / vv_art -
      qp_muscle * pb_v_muscle / vv_muscle +
      psp_muscle * (pb_int_muscle / vi_muscle - pb_v_muscle / vv_muscle) -
      lam_pb * pb_v_muscle
    d/dt(bis_v_muscle) <- qp_muscle * bis_art / vv_art -
      qp_muscle * bis_v_muscle / vv_muscle +
      psp_muscle * (bis_int_muscle / vi_muscle - bis_v_muscle / vv_muscle) +
      lam_gr * pb_v_muscle - lam_bi * bis_v_muscle
    d/dt(bi_v_muscle) <- qp_muscle * bi_art / vv_art -
      qp_muscle * bi_v_muscle / vv_muscle +
      psf_muscle * (bi_int_muscle / vi_muscle - bi_v_muscle / vv_muscle) +
      kdiss_hwpp * bihwpp_v_muscle - kass_hwpp * bi_v_muscle +
      kdiss_rbc * birbc_v_muscle - kass_rbc * bi_v_muscle -
      lam_bi * bi_v_muscle + lam_ex * pb_v_muscle
    d/dt(bihwpp_v_muscle) <- qp_muscle * bihwpp_art / vv_art -
      qp_muscle * bihwpp_v_muscle / vv_muscle +
      kass_hwpp * bi_v_muscle - kdiss_hwpp * bihwpp_v_muscle - lam_bi * bihwpp_v_muscle
    d/dt(birbc_v_muscle) <- qr_muscle * birbc_art / vr_art -
      qr_muscle * birbc_v_muscle / vr_muscle +
      kass_rbc * bi_v_muscle - kdiss_rbc * birbc_v_muscle - lam_bi * birbc_v_muscle
    d/dt(tl_v_muscle) <- qp_muscle * tl_art / vv_art -
      qp_muscle * tl_v_muscle / vv_muscle +
      psf_muscle * (tl_int_muscle / vi_muscle - tl_v_muscle / vv_muscle) +
      lam_bia * (bis_v_muscle + bi_v_muscle + bihwpp_v_muscle + birbc_v_muscle) - lam_tl * tl_v_muscle
    d/dt(cold_v_muscle) <- qp_muscle * cold_art / vv_art -
      qp_muscle * cold_v_muscle / vv_muscle +
      psp_muscle * (cold_int_muscle / vi_muscle - cold_v_muscle / vv_muscle) +
      lam_ex * pb_v_muscle + lam_bi * bis_v_muscle
    d/dt(pb_int_muscle) <- psp_muscle * (pb_v_muscle / vv_muscle - pb_int_muscle / vi_muscle) -
      lam_pb * pb_int_muscle
    d/dt(bis_int_muscle) <- psp_muscle * (bis_v_muscle / vv_muscle - bis_int_muscle / vi_muscle) +
      lam_gr * pb_int_muscle - lam_bi * bis_int_muscle
    d/dt(bi_int_muscle) <- psf_muscle * (bi_v_muscle / vv_muscle - bi_int_muscle / vi_muscle) -
      kup_muscle * bi_int_muscle - lam_bi * bi_int_muscle + lam_ex * pb_int_muscle
    d/dt(tl_int_muscle) <- psf_muscle * (tl_v_muscle / vv_muscle - tl_int_muscle / vi_muscle) +
      lam_bia * (bis_int_muscle + bi_int_muscle) - lam_tl * tl_int_muscle
    d/dt(cold_int_muscle) <- psp_muscle * (cold_v_muscle / vv_muscle - cold_int_muscle / vi_muscle) +
      lam_ex * pb_int_muscle + lam_bi * bis_int_muscle
    d/dt(bi_cell_muscle) <- kup_muscle * bi_int_muscle - lam_bi * bi_cell_muscle
    d/dt(tl_cell_muscle) <- lam_bia * bi_cell_muscle - lam_tl * tl_cell_muscle
    # ---- skin (no sst2) ----
    d/dt(pb_v_skin) <- qp_skin * pb_art / vv_art -
      qp_skin * pb_v_skin / vv_skin +
      psp_skin * (pb_int_skin / vi_skin - pb_v_skin / vv_skin) -
      lam_pb * pb_v_skin
    d/dt(bis_v_skin) <- qp_skin * bis_art / vv_art -
      qp_skin * bis_v_skin / vv_skin +
      psp_skin * (bis_int_skin / vi_skin - bis_v_skin / vv_skin) +
      lam_gr * pb_v_skin - lam_bi * bis_v_skin
    d/dt(bi_v_skin) <- qp_skin * bi_art / vv_art -
      qp_skin * bi_v_skin / vv_skin +
      psf_skin * (bi_int_skin / vi_skin - bi_v_skin / vv_skin) +
      kdiss_hwpp * bihwpp_v_skin - kass_hwpp * bi_v_skin +
      kdiss_rbc * birbc_v_skin - kass_rbc * bi_v_skin -
      lam_bi * bi_v_skin + lam_ex * pb_v_skin
    d/dt(bihwpp_v_skin) <- qp_skin * bihwpp_art / vv_art -
      qp_skin * bihwpp_v_skin / vv_skin +
      kass_hwpp * bi_v_skin - kdiss_hwpp * bihwpp_v_skin - lam_bi * bihwpp_v_skin
    d/dt(birbc_v_skin) <- qr_skin * birbc_art / vr_art -
      qr_skin * birbc_v_skin / vr_skin +
      kass_rbc * bi_v_skin - kdiss_rbc * birbc_v_skin - lam_bi * birbc_v_skin
    d/dt(tl_v_skin) <- qp_skin * tl_art / vv_art -
      qp_skin * tl_v_skin / vv_skin +
      psf_skin * (tl_int_skin / vi_skin - tl_v_skin / vv_skin) +
      lam_bia * (bis_v_skin + bi_v_skin + bihwpp_v_skin + birbc_v_skin) - lam_tl * tl_v_skin
    d/dt(cold_v_skin) <- qp_skin * cold_art / vv_art -
      qp_skin * cold_v_skin / vv_skin +
      psp_skin * (cold_int_skin / vi_skin - cold_v_skin / vv_skin) +
      lam_ex * pb_v_skin + lam_bi * bis_v_skin
    d/dt(pb_int_skin) <- psp_skin * (pb_v_skin / vv_skin - pb_int_skin / vi_skin) -
      lam_pb * pb_int_skin
    d/dt(bis_int_skin) <- psp_skin * (bis_v_skin / vv_skin - bis_int_skin / vi_skin) +
      lam_gr * pb_int_skin - lam_bi * bis_int_skin
    d/dt(bi_int_skin) <- psf_skin * (bi_v_skin / vv_skin - bi_int_skin / vi_skin) -
      kup_skin * bi_int_skin - lam_bi * bi_int_skin + lam_ex * pb_int_skin
    d/dt(tl_int_skin) <- psf_skin * (tl_v_skin / vv_skin - tl_int_skin / vi_skin) +
      lam_bia * (bis_int_skin + bi_int_skin) - lam_tl * tl_int_skin
    d/dt(cold_int_skin) <- psp_skin * (cold_v_skin / vv_skin - cold_int_skin / vi_skin) +
      lam_ex * pb_int_skin + lam_bi * bis_int_skin
    d/dt(bi_cell_skin) <- kup_skin * bi_int_skin - lam_bi * bi_cell_skin
    d/dt(tl_cell_skin) <- lam_bia * bi_cell_skin - lam_tl * tl_cell_skin
    # ---- heart (no sst2) ----
    d/dt(pb_v_heart) <- qp_heart * pb_art / vv_art -
      qp_heart * pb_v_heart / vv_heart +
      psp_heart * (pb_int_heart / vi_heart - pb_v_heart / vv_heart) -
      lam_pb * pb_v_heart
    d/dt(bis_v_heart) <- qp_heart * bis_art / vv_art -
      qp_heart * bis_v_heart / vv_heart +
      psp_heart * (bis_int_heart / vi_heart - bis_v_heart / vv_heart) +
      lam_gr * pb_v_heart - lam_bi * bis_v_heart
    d/dt(bi_v_heart) <- qp_heart * bi_art / vv_art -
      qp_heart * bi_v_heart / vv_heart +
      psf_heart * (bi_int_heart / vi_heart - bi_v_heart / vv_heart) +
      kdiss_hwpp * bihwpp_v_heart - kass_hwpp * bi_v_heart +
      kdiss_rbc * birbc_v_heart - kass_rbc * bi_v_heart -
      lam_bi * bi_v_heart + lam_ex * pb_v_heart
    d/dt(bihwpp_v_heart) <- qp_heart * bihwpp_art / vv_art -
      qp_heart * bihwpp_v_heart / vv_heart +
      kass_hwpp * bi_v_heart - kdiss_hwpp * bihwpp_v_heart - lam_bi * bihwpp_v_heart
    d/dt(birbc_v_heart) <- qr_heart * birbc_art / vr_art -
      qr_heart * birbc_v_heart / vr_heart +
      kass_rbc * bi_v_heart - kdiss_rbc * birbc_v_heart - lam_bi * birbc_v_heart
    d/dt(tl_v_heart) <- qp_heart * tl_art / vv_art -
      qp_heart * tl_v_heart / vv_heart +
      psf_heart * (tl_int_heart / vi_heart - tl_v_heart / vv_heart) +
      lam_bia * (bis_v_heart + bi_v_heart + bihwpp_v_heart + birbc_v_heart) - lam_tl * tl_v_heart
    d/dt(cold_v_heart) <- qp_heart * cold_art / vv_art -
      qp_heart * cold_v_heart / vv_heart +
      psp_heart * (cold_int_heart / vi_heart - cold_v_heart / vv_heart) +
      lam_ex * pb_v_heart + lam_bi * bis_v_heart
    d/dt(pb_int_heart) <- psp_heart * (pb_v_heart / vv_heart - pb_int_heart / vi_heart) -
      lam_pb * pb_int_heart
    d/dt(bis_int_heart) <- psp_heart * (bis_v_heart / vv_heart - bis_int_heart / vi_heart) +
      lam_gr * pb_int_heart - lam_bi * bis_int_heart
    d/dt(bi_int_heart) <- psf_heart * (bi_v_heart / vv_heart - bi_int_heart / vi_heart) -
      kup_heart * bi_int_heart - lam_bi * bi_int_heart + lam_ex * pb_int_heart
    d/dt(tl_int_heart) <- psf_heart * (tl_v_heart / vv_heart - tl_int_heart / vi_heart) +
      lam_bia * (bis_int_heart + bi_int_heart) - lam_tl * tl_int_heart
    d/dt(cold_int_heart) <- psp_heart * (cold_v_heart / vv_heart - cold_int_heart / vi_heart) +
      lam_ex * pb_int_heart + lam_bi * bis_int_heart
    d/dt(bi_cell_heart) <- kup_heart * bi_int_heart - lam_bi * bi_cell_heart
    d/dt(tl_cell_heart) <- lam_bia * bi_cell_heart - lam_tl * tl_cell_heart
    # ---- bone (no sst2) ----
    d/dt(pb_v_bone) <- qp_bone * pb_art / vv_art -
      qp_bone * pb_v_bone / vv_bone +
      psp_bone * (pb_int_bone / vi_bone - pb_v_bone / vv_bone) -
      lam_pb * pb_v_bone
    d/dt(bis_v_bone) <- qp_bone * bis_art / vv_art -
      qp_bone * bis_v_bone / vv_bone +
      psp_bone * (bis_int_bone / vi_bone - bis_v_bone / vv_bone) +
      lam_gr * pb_v_bone - lam_bi * bis_v_bone
    d/dt(bi_v_bone) <- qp_bone * bi_art / vv_art -
      qp_bone * bi_v_bone / vv_bone +
      psf_bone * (bi_int_bone / vi_bone - bi_v_bone / vv_bone) +
      kdiss_hwpp * bihwpp_v_bone - kass_hwpp * bi_v_bone +
      kdiss_rbc * birbc_v_bone - kass_rbc * bi_v_bone -
      lam_bi * bi_v_bone + lam_ex * pb_v_bone
    d/dt(bihwpp_v_bone) <- qp_bone * bihwpp_art / vv_art -
      qp_bone * bihwpp_v_bone / vv_bone +
      kass_hwpp * bi_v_bone - kdiss_hwpp * bihwpp_v_bone - lam_bi * bihwpp_v_bone
    d/dt(birbc_v_bone) <- qr_bone * birbc_art / vr_art -
      qr_bone * birbc_v_bone / vr_bone +
      kass_rbc * bi_v_bone - kdiss_rbc * birbc_v_bone - lam_bi * birbc_v_bone
    d/dt(tl_v_bone) <- qp_bone * tl_art / vv_art -
      qp_bone * tl_v_bone / vv_bone +
      psf_bone * (tl_int_bone / vi_bone - tl_v_bone / vv_bone) +
      lam_bia * (bis_v_bone + bi_v_bone + bihwpp_v_bone + birbc_v_bone) - lam_tl * tl_v_bone
    d/dt(cold_v_bone) <- qp_bone * cold_art / vv_art -
      qp_bone * cold_v_bone / vv_bone +
      psp_bone * (cold_int_bone / vi_bone - cold_v_bone / vv_bone) +
      lam_ex * pb_v_bone + lam_bi * bis_v_bone
    d/dt(pb_int_bone) <- psp_bone * (pb_v_bone / vv_bone - pb_int_bone / vi_bone) -
      lam_pb * pb_int_bone
    d/dt(bis_int_bone) <- psp_bone * (bis_v_bone / vv_bone - bis_int_bone / vi_bone) +
      lam_gr * pb_int_bone - lam_bi * bis_int_bone
    d/dt(bi_int_bone) <- psf_bone * (bi_v_bone / vv_bone - bi_int_bone / vi_bone) -
      kup_bone * bi_int_bone - lam_bi * bi_int_bone + lam_ex * pb_int_bone
    d/dt(tl_int_bone) <- psf_bone * (tl_v_bone / vv_bone - tl_int_bone / vi_bone) +
      lam_bia * (bis_int_bone + bi_int_bone) - lam_tl * tl_int_bone
    d/dt(cold_int_bone) <- psp_bone * (cold_v_bone / vv_bone - cold_int_bone / vi_bone) +
      lam_ex * pb_int_bone + lam_bi * bis_int_bone
    d/dt(bi_cell_bone) <- kup_bone * bi_int_bone - lam_bi * bi_cell_bone
    d/dt(tl_cell_bone) <- lam_bia * bi_cell_bone - lam_tl * tl_cell_bone
    # ---- fat (no sst2) ----
    d/dt(pb_v_fat) <- qp_fat * pb_art / vv_art -
      qp_fat * pb_v_fat / vv_fat +
      psp_fat * (pb_int_fat / vi_fat - pb_v_fat / vv_fat) -
      lam_pb * pb_v_fat
    d/dt(bis_v_fat) <- qp_fat * bis_art / vv_art -
      qp_fat * bis_v_fat / vv_fat +
      psp_fat * (bis_int_fat / vi_fat - bis_v_fat / vv_fat) +
      lam_gr * pb_v_fat - lam_bi * bis_v_fat
    d/dt(bi_v_fat) <- qp_fat * bi_art / vv_art -
      qp_fat * bi_v_fat / vv_fat +
      psf_fat * (bi_int_fat / vi_fat - bi_v_fat / vv_fat) +
      kdiss_hwpp * bihwpp_v_fat - kass_hwpp * bi_v_fat +
      kdiss_rbc * birbc_v_fat - kass_rbc * bi_v_fat -
      lam_bi * bi_v_fat + lam_ex * pb_v_fat
    d/dt(bihwpp_v_fat) <- qp_fat * bihwpp_art / vv_art -
      qp_fat * bihwpp_v_fat / vv_fat +
      kass_hwpp * bi_v_fat - kdiss_hwpp * bihwpp_v_fat - lam_bi * bihwpp_v_fat
    d/dt(birbc_v_fat) <- qr_fat * birbc_art / vr_art -
      qr_fat * birbc_v_fat / vr_fat +
      kass_rbc * bi_v_fat - kdiss_rbc * birbc_v_fat - lam_bi * birbc_v_fat
    d/dt(tl_v_fat) <- qp_fat * tl_art / vv_art -
      qp_fat * tl_v_fat / vv_fat +
      psf_fat * (tl_int_fat / vi_fat - tl_v_fat / vv_fat) +
      lam_bia * (bis_v_fat + bi_v_fat + bihwpp_v_fat + birbc_v_fat) - lam_tl * tl_v_fat
    d/dt(cold_v_fat) <- qp_fat * cold_art / vv_art -
      qp_fat * cold_v_fat / vv_fat +
      psp_fat * (cold_int_fat / vi_fat - cold_v_fat / vv_fat) +
      lam_ex * pb_v_fat + lam_bi * bis_v_fat
    d/dt(pb_int_fat) <- psp_fat * (pb_v_fat / vv_fat - pb_int_fat / vi_fat) -
      lam_pb * pb_int_fat
    d/dt(bis_int_fat) <- psp_fat * (bis_v_fat / vv_fat - bis_int_fat / vi_fat) +
      lam_gr * pb_int_fat - lam_bi * bis_int_fat
    d/dt(bi_int_fat) <- psf_fat * (bi_v_fat / vv_fat - bi_int_fat / vi_fat) -
      kup_fat * bi_int_fat - lam_bi * bi_int_fat + lam_ex * pb_int_fat
    d/dt(tl_int_fat) <- psf_fat * (tl_v_fat / vv_fat - tl_int_fat / vi_fat) +
      lam_bia * (bis_int_fat + bi_int_fat) - lam_tl * tl_int_fat
    d/dt(cold_int_fat) <- psp_fat * (cold_v_fat / vv_fat - cold_int_fat / vi_fat) +
      lam_ex * pb_int_fat + lam_bi * bis_int_fat
    d/dt(bi_cell_fat) <- kup_fat * bi_int_fat - lam_bi * bi_cell_fat
    d/dt(tl_cell_fat) <- lam_bia * bi_cell_fat - lam_tl * tl_cell_fat
    # ---- rb (no sst2) ----
    d/dt(pb_v_rb) <- qp_rb * pb_art / vv_art -
      qp_rb * pb_v_rb / vv_rb +
      psp_rb * (pb_int_rb / vi_rb - pb_v_rb / vv_rb) -
      lam_pb * pb_v_rb
    d/dt(bis_v_rb) <- qp_rb * bis_art / vv_art -
      qp_rb * bis_v_rb / vv_rb +
      psp_rb * (bis_int_rb / vi_rb - bis_v_rb / vv_rb) +
      lam_gr * pb_v_rb - lam_bi * bis_v_rb
    d/dt(bi_v_rb) <- qp_rb * bi_art / vv_art -
      qp_rb * bi_v_rb / vv_rb +
      psf_rb * (bi_int_rb / vi_rb - bi_v_rb / vv_rb) +
      kdiss_hwpp * bihwpp_v_rb - kass_hwpp * bi_v_rb +
      kdiss_rbc * birbc_v_rb - kass_rbc * bi_v_rb -
      lam_bi * bi_v_rb + lam_ex * pb_v_rb
    d/dt(bihwpp_v_rb) <- qp_rb * bihwpp_art / vv_art -
      qp_rb * bihwpp_v_rb / vv_rb +
      kass_hwpp * bi_v_rb - kdiss_hwpp * bihwpp_v_rb - lam_bi * bihwpp_v_rb
    d/dt(birbc_v_rb) <- qr_rb * birbc_art / vr_art -
      qr_rb * birbc_v_rb / vr_rb +
      kass_rbc * bi_v_rb - kdiss_rbc * birbc_v_rb - lam_bi * birbc_v_rb
    d/dt(tl_v_rb) <- qp_rb * tl_art / vv_art -
      qp_rb * tl_v_rb / vv_rb +
      psf_rb * (tl_int_rb / vi_rb - tl_v_rb / vv_rb) +
      lam_bia * (bis_v_rb + bi_v_rb + bihwpp_v_rb + birbc_v_rb) - lam_tl * tl_v_rb
    d/dt(cold_v_rb) <- qp_rb * cold_art / vv_art -
      qp_rb * cold_v_rb / vv_rb +
      psp_rb * (cold_int_rb / vi_rb - cold_v_rb / vv_rb) +
      lam_ex * pb_v_rb + lam_bi * bis_v_rb
    d/dt(pb_int_rb) <- psp_rb * (pb_v_rb / vv_rb - pb_int_rb / vi_rb) -
      lam_pb * pb_int_rb
    d/dt(bis_int_rb) <- psp_rb * (bis_v_rb / vv_rb - bis_int_rb / vi_rb) +
      lam_gr * pb_int_rb - lam_bi * bis_int_rb
    d/dt(bi_int_rb) <- psf_rb * (bi_v_rb / vv_rb - bi_int_rb / vi_rb) -
      kup_rb * bi_int_rb - lam_bi * bi_int_rb + lam_ex * pb_int_rb
    d/dt(tl_int_rb) <- psf_rb * (tl_v_rb / vv_rb - tl_int_rb / vi_rb) +
      lam_bia * (bis_int_rb + bi_int_rb) - lam_tl * tl_int_rb
    d/dt(cold_int_rb) <- psp_rb * (cold_v_rb / vv_rb - cold_int_rb / vi_rb) +
      lam_ex * pb_int_rb + lam_bi * bis_int_rb
    d/dt(bi_cell_rb) <- kup_rb * bi_int_rb - lam_bi * bi_cell_rb
    d/dt(tl_cell_rb) <- lam_bia * bi_cell_rb - lam_tl * tl_cell_rb

    # =====================================================================
    # 4. Outputs. <x>tot_<tissue> are tissue amounts (nmol) of 212Pb-peptide
    #    (pbtot), 212Bi-peptide (bistot), free 212Bi in any form (bitot) and
    #    208Tl (tltot), summed over every sub-compartment of the tissue
    #    (kidney excludes the urine). cpb_<tissue> = pbtot / tissue mass
    #    (nmol/g): divided by the injected 212Pb-peptide amount and times
    #    100 it is the non-decay-corrected %ID/g of 212Pb that Table S2
    #    reports. The brain is vascular-only and at arterial concentration.
    # =====================================================================
    pbtot_tumor <- pb_v_tumor + pb_int_tumor + pb_bnd_tumor + pb_endo_tumor +
      pb_sort_tumor
    bistot_tumor <- bis_v_tumor + bis_int_tumor + bis_bnd_tumor + bis_endo_tumor +
      bis_sort_tumor
    bitot_tumor <- bi_v_tumor + bihwpp_v_tumor + birbc_v_tumor + bi_int_tumor +
      bi_endo_tumor + bi_sort_tumor
    tltot_tumor <- tl_v_tumor + tl_int_tumor + tl_endo_tumor + tl_sort_tumor
    cpb_tumor <- pbtot_tumor / (vtot_tumor * 1000)
    pbtot_kidney <- pb_v_kidney + pb_bnd_kidney + pb_endo_kidney + pb_sort_kidney +
      pb_tub_kidney + pb_nsp_kidney
    bistot_kidney <- bis_v_kidney + bis_bnd_kidney + bis_endo_kidney + bis_sort_kidney +
      bis_tub_kidney + bis_nsp_kidney
    bitot_kidney <- bi_v_kidney + bihwpp_v_kidney + birbc_v_kidney + bi_endo_kidney +
      bi_sort_kidney + bi_tub_kidney + bi_nsp_kidney
    tltot_kidney <- tl_v_kidney + tl_endo_kidney + tl_sort_kidney + tl_tub_kidney +
      tl_nsp_kidney
    cpb_kidney <- pbtot_kidney / (vtot_kidney * 1000)
    pbtot_liver <- pb_v_liver + pb_int_liver + pb_bnd_liver + pb_endo_liver +
      pb_sort_liver
    bistot_liver <- bis_v_liver + bis_int_liver + bis_bnd_liver + bis_endo_liver +
      bis_sort_liver
    bitot_liver <- bi_v_liver + bihwpp_v_liver + birbc_v_liver + bi_int_liver +
      bi_endo_liver + bi_sort_liver
    tltot_liver <- tl_v_liver + tl_int_liver + tl_endo_liver + tl_sort_liver
    cpb_liver <- pbtot_liver / (vtot_liver * 1000)
    pbtot_spleen <- pb_v_spleen + pb_int_spleen + pb_bnd_spleen + pb_endo_spleen +
      pb_sort_spleen
    bistot_spleen <- bis_v_spleen + bis_int_spleen + bis_bnd_spleen + bis_endo_spleen +
      bis_sort_spleen
    bitot_spleen <- bi_v_spleen + bihwpp_v_spleen + birbc_v_spleen + bi_int_spleen +
      bi_endo_spleen + bi_sort_spleen
    tltot_spleen <- tl_v_spleen + tl_int_spleen + tl_endo_spleen + tl_sort_spleen
    cpb_spleen <- pbtot_spleen / (vtot_spleen * 1000)
    pbtot_pancreas <- pb_v_pancreas + pb_int_pancreas + pb_bnd_pancreas +
      pb_endo_pancreas + pb_sort_pancreas
    bistot_pancreas <- bis_v_pancreas + bis_int_pancreas + bis_bnd_pancreas +
      bis_endo_pancreas + bis_sort_pancreas
    bitot_pancreas <- bi_v_pancreas + bihwpp_v_pancreas + birbc_v_pancreas +
      bi_int_pancreas + bi_endo_pancreas + bi_sort_pancreas
    tltot_pancreas <- tl_v_pancreas + tl_int_pancreas + tl_endo_pancreas +
      tl_sort_pancreas
    cpb_pancreas <- pbtot_pancreas / (vtot_pancreas * 1000)
    pbtot_gi <- pb_v_gi + pb_int_gi + pb_bnd_gi + pb_endo_gi + pb_sort_gi
    bistot_gi <- bis_v_gi + bis_int_gi + bis_bnd_gi + bis_endo_gi + bis_sort_gi
    bitot_gi <- bi_v_gi + bihwpp_v_gi + birbc_v_gi + bi_int_gi + bi_endo_gi + bi_sort_gi
    tltot_gi <- tl_v_gi + tl_int_gi + tl_endo_gi + tl_sort_gi
    cpb_gi <- pbtot_gi / (vtot_gi * 1000)
    pbtot_adrenal <- pb_v_adrenal + pb_int_adrenal + pb_bnd_adrenal + pb_endo_adrenal +
      pb_sort_adrenal
    bistot_adrenal <- bis_v_adrenal + bis_int_adrenal + bis_bnd_adrenal +
      bis_endo_adrenal + bis_sort_adrenal
    bitot_adrenal <- bi_v_adrenal + bihwpp_v_adrenal + birbc_v_adrenal + bi_int_adrenal +
      bi_endo_adrenal + bi_sort_adrenal
    tltot_adrenal <- tl_v_adrenal + tl_int_adrenal + tl_endo_adrenal + tl_sort_adrenal
    cpb_adrenal <- pbtot_adrenal / (vtot_adrenal * 1000)
    pbtot_lung <- pb_v_lung + pb_int_lung + pb_bnd_lung + pb_endo_lung + pb_sort_lung
    bistot_lung <- bis_v_lung + bis_int_lung + bis_bnd_lung + bis_endo_lung +
      bis_sort_lung
    bitot_lung <- bi_v_lung + bihwpp_v_lung + birbc_v_lung + bi_int_lung + bi_endo_lung +
      bi_sort_lung
    tltot_lung <- tl_v_lung + tl_int_lung + tl_endo_lung + tl_sort_lung
    cpb_lung <- pbtot_lung / (vtot_lung * 1000)
    pbtot_muscle <- pb_v_muscle + pb_int_muscle
    bistot_muscle <- bis_v_muscle + bis_int_muscle
    bitot_muscle <- bi_v_muscle + bihwpp_v_muscle + birbc_v_muscle + bi_int_muscle +
      bi_cell_muscle
    tltot_muscle <- tl_v_muscle + tl_int_muscle + tl_cell_muscle
    cpb_muscle <- pbtot_muscle / (vtot_muscle * 1000)
    pbtot_skin <- pb_v_skin + pb_int_skin
    bistot_skin <- bis_v_skin + bis_int_skin
    bitot_skin <- bi_v_skin + bihwpp_v_skin + birbc_v_skin + bi_int_skin + bi_cell_skin
    tltot_skin <- tl_v_skin + tl_int_skin + tl_cell_skin
    cpb_skin <- pbtot_skin / (vtot_skin * 1000)
    pbtot_heart <- pb_v_heart + pb_int_heart
    bistot_heart <- bis_v_heart + bis_int_heart
    bitot_heart <- bi_v_heart + bihwpp_v_heart + birbc_v_heart + bi_int_heart +
      bi_cell_heart
    tltot_heart <- tl_v_heart + tl_int_heart + tl_cell_heart
    cpb_heart <- pbtot_heart / (vtot_heart * 1000)
    pbtot_bone <- pb_v_bone + pb_int_bone
    bistot_bone <- bis_v_bone + bis_int_bone
    bitot_bone <- bi_v_bone + bihwpp_v_bone + birbc_v_bone + bi_int_bone + bi_cell_bone
    tltot_bone <- tl_v_bone + tl_int_bone + tl_cell_bone
    cpb_bone <- pbtot_bone / (vtot_bone * 1000)
    pbtot_fat <- pb_v_fat + pb_int_fat
    bistot_fat <- bis_v_fat + bis_int_fat
    bitot_fat <- bi_v_fat + bihwpp_v_fat + birbc_v_fat + bi_int_fat + bi_cell_fat
    tltot_fat <- tl_v_fat + tl_int_fat + tl_cell_fat
    cpb_fat <- pbtot_fat / (vtot_fat * 1000)
    pbtot_rb <- pb_v_rb + pb_int_rb
    bistot_rb <- bis_v_rb + bis_int_rb
    bitot_rb <- bi_v_rb + bihwpp_v_rb + birbc_v_rb + bi_int_rb + bi_cell_rb
    tltot_rb <- tl_v_rb + tl_int_rb + tl_cell_rb
    cpb_rb <- pbtot_rb / (vtot_rb * 1000)
    pbtot_brain <- pb_art / vv_art * vv_brain
    bistot_brain <- bis_art / vv_art * vv_brain
    bitot_brain <- (bi_art + bihwpp_art) / vv_art * vv_brain + birbc_art / vr_art * vr_brain
    tltot_brain <- tl_art / vv_art * vv_brain
    cpb_brain <- pbtot_brain / (vtot_brain * 1000)

    # Total peptide, labelled plus unlabelled, in every state. Decay only
    # moves peptide between the labelled and unlabelled systems, so this
    # equals the injected amount whenever the release rates are zero.
    peptide_total <- pb_art + bis_art + cold_art + pb_ven + bis_ven + cold_ven +
      pb_v_tumor + bis_v_tumor + cold_v_tumor + pb_int_tumor + bis_int_tumor +
      cold_int_tumor + pb_bnd_tumor + bis_bnd_tumor + cold_bnd_tumor + pb_endo_tumor +
      bis_endo_tumor + cold_endo_tumor + pb_sort_tumor + bis_sort_tumor +
      cold_sort_tumor + pb_v_liver + bis_v_liver + cold_v_liver + pb_int_liver +
      bis_int_liver + cold_int_liver + pb_bnd_liver + bis_bnd_liver + cold_bnd_liver +
      pb_endo_liver + bis_endo_liver + cold_endo_liver + pb_sort_liver + bis_sort_liver +
      cold_sort_liver + pb_v_spleen + bis_v_spleen + cold_v_spleen + pb_int_spleen +
      bis_int_spleen + cold_int_spleen + pb_bnd_spleen + bis_bnd_spleen +
      cold_bnd_spleen + pb_endo_spleen + bis_endo_spleen + cold_endo_spleen +
      pb_sort_spleen + bis_sort_spleen + cold_sort_spleen + pb_v_pancreas +
      bis_v_pancreas + cold_v_pancreas + pb_int_pancreas + bis_int_pancreas +
      cold_int_pancreas + pb_bnd_pancreas + bis_bnd_pancreas + cold_bnd_pancreas +
      pb_endo_pancreas + bis_endo_pancreas + cold_endo_pancreas + pb_sort_pancreas +
      bis_sort_pancreas + cold_sort_pancreas + pb_v_gi + bis_v_gi + cold_v_gi +
      pb_int_gi + bis_int_gi + cold_int_gi + pb_bnd_gi + bis_bnd_gi + cold_bnd_gi +
      pb_endo_gi + bis_endo_gi + cold_endo_gi + pb_sort_gi + bis_sort_gi + cold_sort_gi +
      pb_v_adrenal + bis_v_adrenal + cold_v_adrenal + pb_int_adrenal + bis_int_adrenal +
      cold_int_adrenal + pb_bnd_adrenal + bis_bnd_adrenal + cold_bnd_adrenal +
      pb_endo_adrenal + bis_endo_adrenal + cold_endo_adrenal + pb_sort_adrenal +
      bis_sort_adrenal + cold_sort_adrenal + pb_v_lung + bis_v_lung + cold_v_lung +
      pb_int_lung + bis_int_lung + cold_int_lung + pb_bnd_lung + bis_bnd_lung +
      cold_bnd_lung + pb_endo_lung + bis_endo_lung + cold_endo_lung + pb_sort_lung +
      bis_sort_lung + cold_sort_lung + pb_v_kidney + bis_v_kidney + cold_v_kidney +
      pb_bnd_kidney + bis_bnd_kidney + cold_bnd_kidney + pb_endo_kidney +
      bis_endo_kidney + cold_endo_kidney + pb_sort_kidney + bis_sort_kidney +
      cold_sort_kidney + pb_tub_kidney + bis_tub_kidney + cold_tub_kidney +
      pb_nsp_kidney + bis_nsp_kidney + cold_nsp_kidney + pb_exc_kidney + bis_exc_kidney +
      cold_exc_kidney + pb_v_muscle + bis_v_muscle + cold_v_muscle + pb_int_muscle +
      bis_int_muscle + cold_int_muscle + pb_v_skin + bis_v_skin + cold_v_skin +
      pb_int_skin + bis_int_skin + cold_int_skin + pb_v_heart + bis_v_heart +
      cold_v_heart + pb_int_heart + bis_int_heart + cold_int_heart + pb_v_bone +
      bis_v_bone + cold_v_bone + pb_int_bone + bis_int_bone + cold_int_bone + pb_v_fat +
      bis_v_fat + cold_v_fat + pb_int_fat + bis_int_fat + cold_int_fat + pb_v_rb +
      bis_v_rb + cold_v_rb + pb_int_rb + bis_int_rb + cold_int_rb
  })
}
