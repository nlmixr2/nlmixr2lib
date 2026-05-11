Parhiz_2024_mRNALNP <- function() {
  description <- "Preclinical (mouse, C57BL/6, ~25 g). Whole-body PBPK / luciferase-expression model for systemically administered firefly-luciferase mRNA delivered in lipid nanoparticles (bare, untargeted parameterization). Six anatomical regions (blood, lung, heart, kidney, spleen, liver) plus portal-organ and carcass remainder vasculature; each of the five major tissues additionally tracks an intracellular LNP pool, a translatable mRNA pool, and a luciferase signal. The luciferase observation uses the bare-LNP homogenate-assay parameter set (LU/mg protein); the vignette documents how to switch to the bare-LNP BLI, IgG-coated, and PECAM-targeted parameterizations from paper Tables 1-3."
  reference <- "Parhiz H, Shuvaev VV, Li Q, et al. Physiologically based modeling of LNP-mediated delivery of mRNA in the vascular system. Mol Ther Nucleic Acids. 2024;35(2):102175. doi:10.1016/j.omtn.2024.102175. ADAPT 5 control stream archived with the article (Data S1, mmc2.zip)."
  vignette <- "Parhiz_2024_mRNALNP"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (C57BL/6, male, 6-8 weeks old, ~25 g body weight)",
    n_subjects     = 3,
    n_studies      = 1,
    age_range      = "6-8 weeks",
    weight_range   = "~25 g (reference body weight for BioDMET physiology)",
    sex_female_pct = 0,
    disease_state  = "Healthy adult mice (no induced pathology).",
    dose_range     = "8 ug mRNA IV (retro-orbital) single dose; 100 uL bolus per mouse.",
    regions        = "Single-center preclinical (University of Pennsylvania).",
    notes          = "Sample size is n = 3 mice per time point destructively sampled across blood, lung, heart, kidney, spleen, and liver. Physiological parameters (tissue volumes, blood flows) come from the BioDMET database and are fixed; only the 12 LNP / luciferase parameters listed in ini() were estimated.",
    scope_note     = "Mechanistic platform model: no IIV, no residual error - intended for typical-value simulation. See vignette for blood-PK and tissue-biodistribution replication checks against paper Figure 2 (bare LNP, homogenate assay)."
  )

  ini({
    # === Estimated PK parameters (Parhiz 2024 Table 1, Bare LNPs column) ===
    # mL/h tissue-specific non-specific uptake clearances and blood elimination
    lcl_lu  <- log(0.134);  label("Lung non-specific uptake clearance CL_lung (mL/h)")    # Table 1, Bare LNPs
    lcl_he  <- log(0.0999); label("Heart non-specific uptake clearance CL_heart (mL/h)")  # Table 1, Bare LNPs
    lcl_ki  <- log(0.428);  label("Kidney non-specific uptake clearance CL_kidney (mL/h)")# Table 1, Bare LNPs
    lcl_sp  <- log(0.656);  label("Spleen non-specific uptake clearance CL_spleen (mL/h)")# Table 1, Bare LNPs
    lcl_li  <- log(16.3);   label("Liver non-specific uptake clearance CL_liver (mL/h)")  # Table 1, Bare LNPs
    lcl_bl  <- log(1.96);   label("Blood elimination clearance CL_blood (mL/h)")          # Table 1, Bare LNPs
    lklnp   <- log(1.00);   label("Intracellular LNP degradation rate k_deg (1/h)")       # Table 1, Bare LNPs

    # === Estimated PD (luciferase expression) parameters ===
    # Bare LNP, homogenate-assay column of Table 2.
    lslu    <- log(2.57e5); label("Lung luciferase production rate S_lung (LU/mg protein per ug RNA / h)")    # Table 2, Bare LNPs homogenate
    lshe    <- log(6.70e4); label("Heart luciferase production rate S_heart (LU/mg protein per ug RNA / h)")  # Table 2, Bare LNPs homogenate
    lski    <- log(1.12e5); label("Kidney luciferase production rate S_kidney (LU/mg protein per ug RNA / h)")# Table 2, Bare LNPs homogenate
    lssp    <- log(3.82e5); label("Spleen luciferase production rate S_spleen (LU/mg protein per ug RNA / h)")# Table 2, Bare LNPs homogenate
    lsli    <- log(8.64e5); label("Liver luciferase production rate S_liver (LU/mg protein per ug RNA / h)")  # Table 2, Bare LNPs homogenate
    lkluc   <- log(0.0940); label("Luciferase signal decay rate k_luc, non-liver tissues (1/h)")             # Table 2, Bare LNPs homogenate
    lkliv   <- log(0.247);  label("Luciferase signal decay rate k_luc,liver, liver only (1/h)")              # Table 2, Bare LNPs homogenate

    # === Fixed parameter ===
    # Modified-mRNA degradation rate carried over from in vitro work
    # (Anderson 2011, Nucleic Acids Res 39:9329-9338); held fixed during the fit
    # per Methods "Luciferase expression model" subsection.
    lkrna   <- fixed(log(0.114)); label("Pseudouridine-modified mRNA degradation rate k_RNA (1/h)")  # Methods, Anderson 2011 ref 36
  })

  model({
    # === Individual (typical-value) parameter assignments ===
    cl_lu <- exp(lcl_lu)
    cl_he <- exp(lcl_he)
    cl_ki <- exp(lcl_ki)
    cl_sp <- exp(lcl_sp)
    cl_li <- exp(lcl_li)
    cl_bl <- exp(lcl_bl)
    klnp  <- exp(lklnp)
    slu   <- exp(lslu)
    she   <- exp(lshe)
    ski   <- exp(lski)
    ssp   <- exp(lssp)
    sli   <- exp(lsli)
    kluc  <- exp(lkluc)
    kliv  <- exp(lkliv)
    krna  <- exp(lkrna)

    # === Physiological parameters (BioDMET, 25-g mouse) ===
    # Volumes in mL, flows in mL/h. Mirrors the ADAPT 5 control stream
    # (Data S1, mmc2.zip "LNP Model.for"); flow values are mL/min in
    # Table S3 and converted to mL/h here (60 x).
    v_bl   <- 1.533       # Blood total volume (mL)                                  # Table S3, BioDMET
    q_co   <- 605.4       # Cardiac output (mL/h)                                    # Table S3 (10.1 mL/min x 60)
    v_lu   <- 0.182       # Lung total volume (mL)                                   # Table S3
    vv_lu  <- 0.0479      # Lung vascular volume (mL)                                # Table S3
    q_he   <- 59.3        # Heart blood flow (mL/h)                                  # Table S3 (0.988 mL/min x 60)
    v_he   <- 0.136       # Heart total volume (mL)                                  # Table S3
    vv_he  <- 0.0095      # Heart vascular volume (mL)                               # Table S3
    q_ki   <- 111.25      # Kidney blood flow (mL/h)                                 # Table S3 (1.85 mL/min x 60)
    v_ki   <- 0.469       # Kidney total volume (mL)                                 # Table S3
    vv_ki  <- 0.0558      # Kidney vascular volume (mL)                              # Table S3
    q_sp   <- 13.3        # Spleen blood flow (mL/h)                                 # Table S3 (0.222 mL/min x 60)
    v_sp   <- 0.113       # Spleen total volume (mL)                                 # Table S3
    vv_sp  <- 0.025       # Spleen vascular volume (mL)                              # Table S3
    q_li   <- 16.7        # Liver hepatic-artery flow (mL/h)                         # Table S3 (0.278 mL/min x 60, hepatic artery only)
    v_li   <- 1.72        # Liver total volume (mL)                                  # Table S3
    vv_li  <- 0.266       # Liver vascular volume (mL)                               # Table S3
    q_po   <- 132.4       # Portal organs blood flow (mL/h)                          # Table S3 (2.21 mL/min x 60)
    vv_po  <- 0.0313      # Portal organs vascular volume (mL)                       # Table S3
    q_re   <- 272.86      # Remainder blood flow (mL/h)                              # Table S3 (4.55 mL/min x 60)
    vv_re  <- 0.886       # Remainder vascular volume (mL)                           # Table S3

    # === Concentration aliases (vascular spaces) ===
    # For the bare-LNP parameterization there is no target binding,
    # so the "free" vascular concentration equals the total vascular
    # concentration in every organ (supplement equation S5 with R_or = 0
    # collapses to C_free = C_total).
    cv_bl <- blood   / v_bl
    cv_lu <- vp_lu   / vv_lu
    cv_he <- vp_he   / vv_he
    cv_ki <- vp_ki   / vv_ki
    cv_sp <- vp_sp   / vv_sp
    cv_li <- vp_li   / vv_li
    cv_po <- vp_po   / vv_po
    cv_re <- vp_re   / vv_re

    # === Vascular ODEs (amount, ug) ===
    # Venous blood -- supplement "Pharmacokinetic Model / Venous Blood"
    # and ADAPT 5 XP(1).
    d/dt(blood) <- q_he * cv_he +
                    q_ki * cv_ki +
                    q_re * cv_re +
                    (q_li + q_sp + q_po) * cv_li -
                    q_co * cv_bl -
                    cl_bl * cv_bl

    d/dt(bldeg) <- cl_bl * cv_bl - klnp * bldeg

    # Lung vasculature -- ADAPT 5 XP(3). Pre-pulmonary mixing chamber
    # for the entire cardiac output.
    d/dt(vp_lu) <- q_co * cv_bl - q_co * cv_lu - cl_lu * cv_lu

    # Heart, kidney, spleen, portal organs, remainder -- ADAPT 5 XP(5),
    # XP(7), XP(9), XP(2), XP(13). All receive arterial inflow at the
    # lung-vasculature concentration cv_lu.
    d/dt(vp_he) <- q_he * cv_lu - q_he * cv_he - cl_he * cv_he
    d/dt(vp_ki) <- q_ki * cv_lu - q_ki * cv_ki - cl_ki * cv_ki
    d/dt(vp_sp) <- q_sp * cv_lu - q_sp * cv_sp - cl_sp * cv_sp
    d/dt(vp_po) <- q_po * cv_lu - q_po * cv_po
    d/dt(vp_re) <- q_re * cv_lu - q_re * cv_re

    # Liver vasculature -- ADAPT 5 XP(11). Receives hepatic-artery
    # inflow plus portal returns from spleen and portal organs.
    d/dt(vp_li) <- q_li * cv_lu + q_sp * cv_sp + q_po * cv_po -
                    (q_li + q_sp + q_po) * cv_li -
                    cl_li * cv_li

    # === Intracellular LNP, mRNA, luciferase ODEs ===
    # Per tissue: LNP enters intracellular via CL_org * C_free,vasc, is
    # degraded at klnp; the released mRNA pool is mRNA-degraded at krna;
    # mRNA drives a tissue-specific luciferase production rate S_org;
    # luciferase signal decays at kluc (kliv in liver). ADAPT 5 XP(4),
    # XP(15), XP(20) for lung and analogous indices for the other tissues.
    d/dt(int_lu)  <- cl_lu * cv_lu - klnp * int_lu
    d/dt(mrna_lu) <- klnp * int_lu - krna * mrna_lu
    d/dt(luc_lu)  <- slu * mrna_lu - kluc * luc_lu

    d/dt(int_he)  <- cl_he * cv_he - klnp * int_he
    d/dt(mrna_he) <- klnp * int_he - krna * mrna_he
    d/dt(luc_he)  <- she * mrna_he - kluc * luc_he

    d/dt(int_ki)  <- cl_ki * cv_ki - klnp * int_ki
    d/dt(mrna_ki) <- klnp * int_ki - krna * mrna_ki
    d/dt(luc_ki)  <- ski * mrna_ki - kluc * luc_ki

    d/dt(int_sp)  <- cl_sp * cv_sp - klnp * int_sp
    d/dt(mrna_sp) <- klnp * int_sp - krna * mrna_sp
    d/dt(luc_sp)  <- ssp * mrna_sp - kluc * luc_sp

    # Liver uses kliv (separate luciferase decay constant per Table 2).
    d/dt(int_li)  <- cl_li * cv_li - klnp * int_li
    d/dt(mrna_li) <- klnp * int_li - krna * mrna_li
    d/dt(luc_li)  <- sli * mrna_li - kliv * luc_li

    # === Derived observation aliases ===
    # Blood concentration (ug/mL) -- canonical Cc observation.
    Cc <- cv_bl

    # Percent-injected-dose-per-gram tissue (paper Figure 2 axis).
    # ADAPT 5 OUTPUT block Y(1)-Y(6); the divisor 8 ug = injected mRNA
    # dose / mouse appears explicitly so that re-scaling to a different
    # dose preserves the %ID/g interpretation.
    pid_g_blood  <- 100 * (blood + bldeg)   / (8 * v_bl)
    pid_g_lung   <- 100 * (vp_lu + int_lu)  / (8 * v_lu)
    pid_g_heart  <- 100 * (vp_he + int_he)  / (8 * v_he)
    pid_g_kidney <- 100 * (vp_ki + int_ki)  / (8 * v_ki)
    pid_g_spleen <- 100 * (vp_sp + int_sp)  / (8 * v_sp)
    pid_g_liver  <- 100 * (vp_li + int_li)  / (8 * v_li)

    # Luciferase signal per tissue -- ADAPT 5 OUTPUT block Y(7)-Y(11).
    luc_lung   <- luc_lu
    luc_heart  <- luc_he
    luc_kidney <- luc_ki
    luc_spleen <- luc_sp
    luc_liver  <- luc_li
  })
}
