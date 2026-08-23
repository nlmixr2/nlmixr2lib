Wang_2024_ionizableLipid_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat, Sprague-Dawley, 225-250 g).",
    "PBPK (semi-mechanistic six-compartment, MATLAB SimBiology 2022a) model for",
    "the disposition of the ionizable lipid component of an intravenous",
    "mRNA-lipid-nanoparticle (LNP). Vein, artery, lung blood vessel, liver,",
    "spleen and a lumped 'other organs' compartment; liver and spleen are each",
    "resolved into vascular, interstitial and intracellular spaces. LNP crosses",
    "the sinusoidal endothelium by passive permeation (P x S), is internalised",
    "by receptor-mediated uptake whose rate is proportional to an LDL-receptor",
    "concentration, then disassembles (kdis) to release free ionizable lipid",
    "which is hydrolysed by esterase (kel). Three ionizable lipids are covered:",
    "DLin-MC3-DMA (MC3, the reference) plus SM-102 and Lipid 5, selected with",
    "FORM_LNP_SM102 / FORM_LNP_LIPID5; the uptake, disassembly and metabolism",
    "rates of the two comparators are expressed as scaling factors on the MC3",
    "values exactly as the authors parameterised them.",
    sep = " "
  )
  reference <- paste(
    "Wang W, Deng S, Lin J, Ouyang D (2024).",
    "Modeling on in vivo disposition and cellular transportation of RNA lipid",
    "nanoparticles via quantum mechanics/physiologically-based pharmacokinetic",
    "approaches. Acta Pharm Sin B 14(10):4591-4607.",
    "doi:10.1016/j.apsb.2024.06.011.",
    "Model equations Eqs. (1)-(9) and the full SimBiology ODE export in the",
    "Supporting Information (mmc1.pdf) section 1; physiology from Table 1",
    "('Rat' column); fitted parameters from Figure 4A.",
    sep = " "
  )
  vignette <- "Wang_2024_rnaLipidNanoparticle"
  units <- list(time = "h", dosing = "mg", concentration = "mg/mL")

  covariateData <- list(
    FORM_LNP_SM102 = list(
      description        = "1 = LNP formulated with the ionizable lipid SM-102; 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (DLin-MC3-DMA / MC3 reference formulation)",
      notes              = paste(
        "Selects the SM-102 column of Wang 2024 Figure 4A. SM-102 shares the",
        "MC3 liver permeability (the paper states 'the permeability rate of MC3",
        "applies to SM-102') but has its own spleen and 'other' permeabilities,",
        "its own 'other'-organ metabolism rate, and its own scaling factors on",
        "kin / kdis / kel. Set FORM_LNP_SM102 = FORM_LNP_LIPID5 = 0 for MC3;",
        "exactly one of the two may be 1."
      ),
      source_name        = "SM-102"
    ),
    FORM_LNP_LIPID5 = list(
      description        = "1 = LNP formulated with the ionizable lipid Lipid 5; 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (DLin-MC3-DMA / MC3 reference formulation)",
      notes              = paste(
        "Selects the Lipid 5 column of Wang 2024 Figure 4A. Unlike SM-102,",
        "Lipid 5 required a down-regulated liver permeability (paper section",
        "3.1.2: 'that of Lipid 5 had to be down-regulated to get a satisfactory",
        "fitting'). Mutually exclusive with FORM_LNP_SM102."
      ),
      source_name        = "Lipid 5"
    )
  )

  # Every state below is the ionizable-lipid MASS (mg) in one anatomical
  # sub-space. The `_np` suffix marks lipid still travelling inside an intact
  # LNP; the two bare `int_<organ>` states hold free ionizable lipid released
  # by LNP disassembly, per the registered `np` (nanoparticle-conjugated
  # species) suffix convention.
  compartmentData <- list(
    venous_np      = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_np    = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_lung_np     = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_liver_np    = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_liver_np    = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_np   = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver      = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_spleen_np   = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_spleen_np   = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_np  = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen     = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_other_np    = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    int_other_np   = list(analyte = "ionizable lipid (in LNP)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species      = "rat (Sprague-Dawley)",
    n_subjects   = NA_integer_,
    n_studies    = 2L,
    weight_range = "225-250 g (0.2375 kg reference body weight, Table 1)",
    dose_range   = paste(
      "0.2 mg/kg mRNA as a single intravenous mRNA-LNP dose, equivalent to",
      "about 2.23 mg/kg MC3 and 2.46 mg/kg SM-102 or Lipid 5 (Figure 4A legend)"
    ),
    disease_state = "healthy",
    notes = paste(
      "Digitised mean ionizable-lipid concentrations in liver and spleen from",
      "Sabnis 2018 (Mol Ther 26:1509-1519) and the Benenato US9868691B2 patent,",
      "which Wang 2024 cites as references 13 and 25. LNP composition was",
      "ionizable lipid : DSPC : cholesterol : PEG-lipid at 50:10:38.5:1.5 molar",
      "ratio with an estimated N/P ratio of 5.67. The model was fitted to pooled",
      "mean profiles, not to individual animals, so no subject count is defined",
      "and the fit carries no between-subject variability."
    )
  )

  ini({
    # ================================================================
    # Endothelial permeability P (Figure 4A, "P to <organ>" rows).
    # Reported by the authors in cm/s; the value is used numerically as
    # the mL/h-scale coefficient of the P x S permeation term, with S the
    # endothelium surface area in cm^2 from Table 1 (see the note on
    # ps_liver in model()).
    # MC3 and SM-102 share one liver permeability (section 3.1.2:
    # "The permeability rate of MC3 applies to SM-102, but that of Lipid 5
    # had to be down-regulated to get a satisfactory fitting").
    # ================================================================
    lp_liver_mc3      <- log(2.6112) ; label("Liver endothelial permeability, MC3 and SM-102 (cm/s)")   # Fig 4A, "P to liver", MC3 = SM-102 column
    lp_liver_lipid5   <- log(0.5222) ; label("Liver endothelial permeability, Lipid 5 (cm/s)")          # Fig 4A, "P to liver", Lipid 5 column
    lp_spleen_mc3     <- log(0.0026) ; label("Spleen endothelial permeability, MC3 (cm/s)")             # Fig 4A, "P to spleen", MC3 column
    lp_spleen_sm102   <- log(0.0020) ; label("Spleen endothelial permeability, SM-102 (cm/s)")          # Fig 4A, "P to spleen", SM-102 column
    lp_spleen_lipid5  <- log(5.0e-4) ; label("Spleen endothelial permeability, Lipid 5 (cm/s)")         # Fig 4A, "P to spleen", Lipid 5 column
    lp_other_mc3      <- log(0.0135) ; label("'Other organs' endothelial permeability, MC3 (cm/s)")     # Fig 4A, "P to other", MC3 column
    lp_other_sm102    <- log(2.1e-5) ; label("'Other organs' endothelial permeability, SM-102 (cm/s)")  # Fig 4A, "P to other", SM-102 column
    lp_other_lipid5   <- log(0.0020) ; label("'Other organs' endothelial permeability, Lipid 5 (cm/s)") # Fig 4A, "P to other", Lipid 5 column

    # ================================================================
    # Metabolism (esterase hydrolysis) of free ionizable lipid, Eq. (6).
    # The liver and spleen rates are MC3 values scaled by scale_kel; the
    # 'other organs' rate is fitted independently per lipid and carries no
    # scaling factor (SimBiology export: kel_other has no kel_scale term).
    # ================================================================
    lkel_liver        <- log(3.7943) ; label("Free-lipid metabolism rate in liver, MC3 reference (1/h)")   # Fig 4A, "kel in liver", MC3 column
    lkel_spleen       <- log(0.0627) ; label("Free-lipid metabolism rate in spleen, MC3 reference (1/h)")  # Fig 4A, "kel in spleen", MC3 column
    lkel_other_mc3    <- log(0.0143) ; label("Lipid elimination rate in 'other organs', MC3 (1/h)")        # Fig 4A, "kel in other", MC3 column
    lkel_other_sm102  <- log(0.0200) ; label("Lipid elimination rate in 'other organs', SM-102 (1/h)")     # Fig 4A, "kel in other", SM-102 column
    lkel_other_lipid5 <- log(0.0495) ; label("Lipid elimination rate in 'other organs', Lipid 5 (1/h)")    # Fig 4A, "kel in other", Lipid 5 column

    # ================================================================
    # Receptor-mediated cellular uptake, Eq. (4), and LNP disassembly,
    # Eq. (5). Both are MC3 reference values; SM-102 and Lipid 5 enter
    # through the Eq. (8) scaling factors X = X_MC3 * Scale_X.
    # ================================================================
    lkin              <- log(20.0013) ; label("Receptor-mediated LNP uptake rate, MC3 reference (mL/mmol/h)") # Fig 4A, "kin", MC3 column
    lkdis             <- log(0.0193)  ; label("LNP disassembly rate, MC3 reference (1/h)")                    # Fig 4A, "kdis", MC3 column

    # Eq. (8) scaling factors. The MC3 scales are 1 by construction and are
    # therefore not carried as parameters.
    lscale_kin_sm102    <- log(0.2961)   ; label("Scale of kin for SM-102 relative to MC3 (unitless)")   # Fig 4A, "Scale of kin", SM-102 column
    lscale_kin_lipid5   <- log(0.1043)   ; label("Scale of kin for Lipid 5 relative to MC3 (unitless)")  # Fig 4A, "Scale of kin", Lipid 5 column
    lscale_kdis_sm102   <- log(85.3820)  ; label("Scale of kdis for SM-102 relative to MC3 (unitless)")  # Fig 4A, "Scale of kdis", SM-102 column
    lscale_kdis_lipid5  <- log(137.2342) ; label("Scale of kdis for Lipid 5 relative to MC3 (unitless)") # Fig 4A, "Scale of kdis", Lipid 5 column
    lscale_kel_sm102    <- log(1.1206)   ; label("Scale of kel for SM-102 relative to MC3 (unitless)")   # Fig 4A, "Scale of kel", SM-102 column
    lscale_kel_lipid5   <- log(2.9830)   ; label("Scale of kel for Lipid 5 relative to MC3 (unitless)")  # Fig 4A, "Scale of kel", Lipid 5 column

    # ================================================================
    # Uptake-receptor concentration, Eq. (9). The liver concentration is an
    # assumed anchor (Table 1 footnote b: "The 'receptor concentration in the
    # liver' is assumed to be 1"); only the liver-to-spleen ratio is fitted.
    # ================================================================
    creceptor_liver   <- fixed(1.000) ; label("LDL-receptor concentration in liver, assumed anchor (mmol/mL)") # Table 1, "Receptor concentration in liver" + footnote b
    lscale_receptor   <- log(0.3729)  ; label("Receptor concentration ratio, spleen to liver (unitless)")      # Table 1, "Receptor scaling from the liver to spleen" (fitted, footnote b)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology of the 0.2375 kg Sprague-Dawley rat.
    #    Volumes in mL, blood flows in mL/h, endothelium surface areas in
    #    cm^2, all from Wang 2024 Table 1, "Rat" column. Sources cited by
    #    the authors: Charles River growth curves and the PK-Sim 11 (Open
    #    Systems Pharmacology) physiology database.
    # ------------------------------------------------------------------
    v_venous       <- 7.830    # Vein (mL)                          # Table 1
    v_arterial     <- 3.840    # Artery (mL)                        # Table 1
    v_lung_vas     <- 0.6552   # Lung blood vessel (mL)             # Table 1
    v_liver_vas    <- 1.289    # Liver blood vessel (mL)            # Table 1
    v_liver_inter  <- 1.718    # Liver interstitium (mL)            # Table 1
    v_liver_cell   <- 7.733    # Liver cell (mL)                    # Table 1
    v_spleen_vas   <- 0.1753   # Spleen blood vessel (mL)           # Table 1
    v_spleen_inter <- 0.09390  # Spleen interstitium (mL)           # Table 1
    v_spleen_cell  <- 0.3568   # Spleen cell (mL)                   # Table 1
    v_other_vas    <- 7.043    # Other organs' blood vessel (mL)    # Table 1
    v_other_cell   <- 209.5    # Other organs' cell (mL)            # Table 1

    q_lung         <- 2400     # Lung blood flow (mL/h)             # Table 1
    q_hepart       <- 697.7    # Hepatic artery blood flow (mL/h)   # Table 1 + footnote a (hepatic artery + intestines + pancreas)
    q_liver        <- 737.3    # Hepatic vein blood flow (mL/h)     # Table 1
    q_spleen       <- 39.60    # Spleen blood flow (mL/h)           # Table 1
    q_other        <- 1663     # Other organs' blood flow (mL/h)    # Table 1

    sa_liver       <- 1173     # Liver endothelium surface area (cm^2)          # Table 1
    sa_spleen      <- 167.6    # Spleen endothelium surface area (cm^2)         # Table 1
    sa_other       <- 6123     # Other organs' endothelium surface area (cm^2)  # Table 1
    hct            <- 0.4500   # Hematocrit (unitless)                          # Table 1

    # ------------------------------------------------------------------
    # 2. Formulation selection. Both indicators zero selects MC3.
    # ------------------------------------------------------------------
    is_mc3 <- (1 - FORM_LNP_SM102) * (1 - FORM_LNP_LIPID5)

    # MC3 and SM-102 share the liver permeability; only Lipid 5 differs.
    p_liver  <- exp(lp_liver_mc3) * (1 - FORM_LNP_LIPID5) +
      exp(lp_liver_lipid5) * FORM_LNP_LIPID5
    p_spleen <- exp(lp_spleen_mc3) * is_mc3 +
      exp(lp_spleen_sm102) * FORM_LNP_SM102 +
      exp(lp_spleen_lipid5) * FORM_LNP_LIPID5
    p_other  <- exp(lp_other_mc3) * is_mc3 +
      exp(lp_other_sm102) * FORM_LNP_SM102 +
      exp(lp_other_lipid5) * FORM_LNP_LIPID5
    kel_other <- exp(lkel_other_mc3) * is_mc3 +
      exp(lkel_other_sm102) * FORM_LNP_SM102 +
      exp(lkel_other_lipid5) * FORM_LNP_LIPID5

    # Eq. (8): X = X_MC3 * Scale_X for X in {kin, kdis, kel}.
    scale_kin  <- is_mc3 +
      exp(lscale_kin_sm102) * FORM_LNP_SM102 +
      exp(lscale_kin_lipid5) * FORM_LNP_LIPID5
    scale_kdis <- is_mc3 +
      exp(lscale_kdis_sm102) * FORM_LNP_SM102 +
      exp(lscale_kdis_lipid5) * FORM_LNP_LIPID5
    scale_kel  <- is_mc3 +
      exp(lscale_kel_sm102) * FORM_LNP_SM102 +
      exp(lscale_kel_lipid5) * FORM_LNP_LIPID5

    kin        <- exp(lkin) * scale_kin
    kdis       <- exp(lkdis) * scale_kdis
    kel_liver  <- exp(lkel_liver) * scale_kel
    kel_spleen <- exp(lkel_spleen) * scale_kel

    # Eq. (9): C_receptor,spleen = C_receptor,liver * Scale_receptor.
    creceptor_spleen <- creceptor_liver * exp(lscale_receptor)

    # ------------------------------------------------------------------
    # 3. Permeation coefficients, Eq. (3): dM/dt = P * S * (Cblood - Cinter).
    #    P is tabulated in nominal cm/s and S in cm^2; the product is used
    #    directly as a mL/h exchange coefficient rather than being converted
    #    from cm^3/s (which would multiply it by 3600). The paper's own
    #    parameter-sensitivity analysis settles this: Figure S4 shows that
    #    permeability "shows a minimal influence on the liver PK profile, but
    #    the spleen PK profile is very sensitive to the change of
    #    permeability". Without the 3600 factor P*S is 3062 mL/h for liver
    #    (4x the 737 mL/h hepatic flow, hence perfusion-limited and
    #    insensitive) and 0.44 mL/h for spleen (1/90 of the 39.6 mL/h splenic
    #    flow, hence permeability-limited and highly sensitive), reproducing
    #    both statements. With the 3600 factor both organs would be strongly
    #    perfusion-limited and neither would be sensitive to P.
    # ------------------------------------------------------------------
    ps_liver  <- p_liver  * sa_liver
    ps_spleen <- p_spleen * sa_spleen
    ps_other  <- p_other  * sa_other

    # ------------------------------------------------------------------
    # 4. Sub-space concentrations (mg/mL) driving the mass-transfer terms.
    # ------------------------------------------------------------------
    c_venous       <- venous_np     / v_venous
    c_arterial     <- arterial_np   / v_arterial
    c_lung_vas     <- vp_lung_np    / v_lung_vas
    c_liver_vas    <- vp_liver_np   / v_liver_vas
    c_liver_inter  <- is_liver_np   / v_liver_inter
    c_liver_cell   <- int_liver_np  / v_liver_cell
    c_liver_free   <- int_liver     / v_liver_cell
    c_spleen_vas   <- vp_spleen_np  / v_spleen_vas
    c_spleen_inter <- is_spleen_np  / v_spleen_inter
    c_spleen_cell  <- int_spleen_np / v_spleen_cell
    c_spleen_free  <- int_spleen    / v_spleen_cell
    c_other_vas    <- vp_other_np   / v_other_vas
    c_other_cell   <- int_other_np  / v_other_cell

    # ------------------------------------------------------------------
    # 5. ODE system. States are lipid MASS (mg); the SimBiology export
    #    (Supporting Information section 1.2) writes the same balances as
    #    concentrations by dividing each right-hand side by the compartment
    #    volume, which is undone here.
    #
    #    The lung is a pure vascular pass-through: Figure 2 draws it as a
    #    "Lung blood vessel" box with no route into lung tissue, Table 1
    #    reports no lung endothelium surface area, and the export's lung
    #    permeation coefficient PerLung = ConRatLung2Oth * PerOther has no
    #    reported value for ConRatLung2Oth. Encoded as zero lung uptake.
    # ------------------------------------------------------------------
    # Blood pool
    d/dt(venous_np)   <- -q_lung * c_venous + q_other * c_other_vas + q_liver * c_liver_vas
    d/dt(arterial_np) <-  q_lung * c_lung_vas - q_hepart * c_arterial -
      q_other * c_arterial - q_spleen * c_arterial
    d/dt(vp_lung_np)  <-  q_lung * c_venous - q_lung * c_lung_vas

    # Liver: Eqs. (1)-(7). Splenic outflow drains into the liver vasculature.
    d/dt(vp_liver_np)  <-  q_spleen * c_spleen_vas + q_hepart * c_arterial -
      q_liver * c_liver_vas - ps_liver * (c_liver_vas - c_liver_inter)
    d/dt(is_liver_np)  <-  ps_liver * (c_liver_vas - c_liver_inter) -
      kin * c_liver_inter * creceptor_liver * v_liver_inter
    d/dt(int_liver_np) <-  kin * c_liver_inter * creceptor_liver * v_liver_inter -
      kdis * c_liver_cell * v_liver_cell
    d/dt(int_liver)    <-  kdis * c_liver_cell * v_liver_cell -
      kel_liver * c_liver_free * v_liver_cell

    # Spleen: same structure as the liver.
    d/dt(vp_spleen_np)  <-  q_spleen * c_arterial - q_spleen * c_spleen_vas -
      ps_spleen * (c_spleen_vas - c_spleen_inter)
    d/dt(is_spleen_np)  <-  ps_spleen * (c_spleen_vas - c_spleen_inter) -
      kin * c_spleen_inter * creceptor_spleen * v_spleen_inter
    d/dt(int_spleen_np) <-  kin * c_spleen_inter * creceptor_spleen * v_spleen_inter -
      kdis * c_spleen_cell * v_spleen_cell
    d/dt(int_spleen)    <-  kdis * c_spleen_cell * v_spleen_cell -
      kel_spleen * c_spleen_free * v_spleen_cell

    # 'Other organs': a mass-balancing sink with a single deep-tissue pool
    # and direct first-order elimination (Figure 2, "Other organs" panel).
    d/dt(vp_other_np)  <- q_other * c_arterial - q_other * c_other_vas -
      ps_other * (c_other_vas - c_other_cell)
    d/dt(int_other_np) <- ps_other * (c_other_vas - c_other_cell) -
      kel_other * c_other_cell * v_other_cell

    # ------------------------------------------------------------------
    # 6. Observations, from the SimBiology "repeated assignment" block.
    #
    #    Cliver / Cspleen are the volume-weighted whole-organ lipid
    #    concentrations plotted in Figure 4A. The liver average includes its
    #    vascular space; the spleen average excludes it, because the export
    #    multiplies both spleen vascular terms by zero
    #    (Spleen_vas.Drug*Spleen_vas*vas_drug_control*0 and
    #    Spleen_vas*vas_org_control*0).
    # ------------------------------------------------------------------
    Cliver  <- (vp_liver_np + is_liver_np + int_liver_np + int_liver) /
      (v_liver_vas + v_liver_inter + v_liver_cell)
    Cspleen <- (is_spleen_np + int_spleen_np + int_spleen) /
      (v_spleen_inter + v_spleen_cell)
    Cblood  <- (venous_np + arterial_np + vp_lung_np + vp_liver_np +
                  vp_spleen_np + vp_other_np) /
      (v_venous + v_arterial + v_lung_vas + v_liver_vas +
         v_spleen_vas + v_other_vas)
    Cplasma <- Cblood / (1 - hct)
  })
}
