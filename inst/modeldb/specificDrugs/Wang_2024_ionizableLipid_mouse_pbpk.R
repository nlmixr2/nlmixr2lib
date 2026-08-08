Wang_2024_ionizableLipid_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse, CD-1 and C57BL/6, 6-8 weeks, 27 g).",
    "PBPK (semi-mechanistic six-compartment, MATLAB SimBiology 2022a) model for",
    "the disposition of the lipid component of an intravenous siRNA-lipid",
    "nanoparticle (LNP), used to compare particle sizes. Same six-compartment",
    "skeleton as the rat model (vein, artery, lung blood vessel, liver, spleen,",
    "'other organs') but with LNP disassembly and ionizable-lipid hydrolysis",
    "switched off in the liver and spleen, because the mouse experiments only",
    "resolve processes up to organ distribution; elimination is retained in the",
    "'other organs' compartment to trap superfluous lipid. A first-order",
    "dissociation of ionizable lipid from the LNP (kf) acts in every blood and",
    "interstitial space and drains into free-lipid accumulator states. Three",
    "formulations are covered: MC3 at 80 nm (the reference) plus DMAP-BLP at",
    "78 nm and 45 nm, selected with FORM_LNP_DMAPBLP78 / FORM_LNP_DMAPBLP45.",
    sep = " "
  )
  reference <- paste(
    "Wang W, Deng S, Lin J, Ouyang D (2024).",
    "Modeling on in vivo disposition and cellular transportation of RNA lipid",
    "nanoparticles via quantum mechanics/physiologically-based pharmacokinetic",
    "approaches. Acta Pharm Sin B 14(10):4591-4607.",
    "doi:10.1016/j.apsb.2024.06.011.",
    "Model equations Eqs. (1)-(4) and the full SimBiology ODE export in the",
    "Supporting Information (mmc1.pdf) section 2 ('The in vivo model (with LNP",
    "dissociation)'); physiology from Table 1 ('Mice' column); fitted",
    "parameters from Figures 5A-5C and 6A.",
    sep = " "
  )
  vignette <- "Wang_2024_rnaLipidNanoparticle"
  units <- list(time = "h", dosing = "mg", concentration = "mg/mL")

  covariateData <- list(
    FORM_LNP_DMAPBLP78 = list(
      description        = "1 = LNP formulated with the ionizable lipid DMAP-BLP at ~78 nm diameter; 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MC3 LNP at ~80 nm, Wang 2024 Figure 5A)",
      notes              = paste(
        "Selects the Figure 5B parameter column. Only the spleen and 'other'",
        "permeabilities were re-fitted for this arm; the liver permeability,",
        "'other'-organ elimination rate and uptake rate carry the '#' footnote",
        "'the value was assumed to be the same as the 80 nm LNP'. The plasma",
        "dissociation rate kf = 0.0414 1/h comes from Figure 6A. Wang 2024",
        "spells this lipid DMAP-BLP in Methods section 2.1.1 and in the Figure",
        "5B/5C parameter tables but DMAP-DLP in Results section 3.1.3 and",
        "Figure S6; the Methods / figure-table spelling is used here.",
        "Mutually exclusive with FORM_LNP_DMAPBLP45."
      ),
      source_name        = "DMAP-BLP, 78 nm"
    ),
    FORM_LNP_DMAPBLP45 = list(
      description        = "1 = LNP formulated with the ionizable lipid DMAP-BLP at ~45 nm diameter; 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MC3 LNP at ~80 nm, Wang 2024 Figure 5A)",
      notes              = paste(
        "Selects the Figure 5C parameter column: every parameter was re-fitted",
        "for the 45 nm arm, giving a much higher uptake rate (125.06 vs 5.24",
        "mL/mmol/h) and a lower liver permeability than the larger particles.",
        "The plasma dissociation rate kf = 0.1089 1/h comes from Figure 6A.",
        "Mutually exclusive with FORM_LNP_DMAPBLP78."
      ),
      source_name        = "DMAP-BLP, 45 nm"
    )
  )

  # `_np` states hold lipid still associated with an intact LNP; the bare
  # states are the free lipid released by the kf dissociation process in the
  # blood and interstitial spaces (SimBiology `.FreeLipid` species). The bare
  # states are pure accumulators: dissociated lipid is not re-taken up.
  compartmentData <- list(
    venous_np     = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    venous        = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_np   = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial      = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_lung_np    = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_lung       = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_liver_np   = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_liver      = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_liver_np   = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "tissue", verified = TRUE),
    is_liver      = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_np  = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_spleen_np  = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_spleen     = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_spleen_np  = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "tissue", verified = TRUE),
    is_spleen     = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen_np = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_other_np   = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_other      = list(analyte = "ionizable lipid (free)", units = "mg", specimen = "whole blood", verified = TRUE),
    int_other_np  = list(analyte = "LNP lipid (3H-CHE traced)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species      = "mouse (CD-1 for the MC3 arm; C57BL/6 for the DMAP-BLP arms)",
    n_subjects   = NA_integer_,
    n_studies    = 2L,
    age_range    = "6-8 weeks",
    weight_range = "27 g (0.0270 kg reference body weight, Table 1)",
    dose_range   = paste(
      "MC3 arm: 11.1 mg/kg total lipid intravenously (equivalent to 1 mg/kg",
      "siRNA); DMAP-BLP arms: 0.3 mg/kg FVII siRNA, equivalent to about",
      "3.42 mg/kg ionizable lipid"
    ),
    disease_state = "healthy",
    notes = paste(
      "Digitised mean organ lipid concentrations from the two mouse studies",
      "Wang 2024 cites as references 14 and 16. Both traced particle",
      "distribution with 3H-CHE, a non-exchangeable and non-metabolisable lipid",
      "label, so the measured analyte is LNP lipid rather than free ionizable",
      "lipid; this is why LNP disassembly and hydrolysis are switched off in",
      "the liver and spleen for the mouse fits. The MC3 study showed a linear",
      "dose-exposure relationship over 0.33-11.1 mg/kg lipid, so only the",
      "11.1 mg/kg arm was fitted. Fitted to pooled mean profiles, so the model",
      "carries no between-subject variability."
    )
  )

  ini({
    # ================================================================
    # Endothelial permeability P (Figures 5A-5C, "P to <organ>" rows),
    # reported in nominal cm/s and used numerically as the mL/h-scale
    # coefficient of the P x S permeation term (see ps_liver in model()).
    # The 78 nm arm reuses the 80 nm liver permeability, uptake rate and
    # 'other'-organ elimination rate under the Figure 5B '#' footnote.
    # ================================================================
    lp_liver_mc3        <- log(5.7086) ; label("Liver endothelial permeability, MC3 80 nm and DMAP-BLP 78 nm (cm/s)") # Fig 5A "P to liver"; Fig 5B carries the same value with footnote #
    lp_liver_dmap45     <- log(0.4650) ; label("Liver endothelial permeability, DMAP-BLP 45 nm (cm/s)")               # Fig 5C, "P to liver"
    lp_spleen_mc3       <- log(5.0e-4) ; label("Spleen endothelial permeability, MC3 80 nm (cm/s)")                   # Fig 5A, "P to spleen"
    lp_spleen_dmap78    <- log(0.0012) ; label("Spleen endothelial permeability, DMAP-BLP 78 nm (cm/s)")              # Fig 5B, "P to spleen"
    lp_spleen_dmap45    <- log(2.3e-4) ; label("Spleen endothelial permeability, DMAP-BLP 45 nm (cm/s)")              # Fig 5C, "P to spleen"
    lp_other_mc3        <- log(6.6e-4) ; label("'Other organs' endothelial permeability, MC3 80 nm (cm/s)")           # Fig 5A, "P to other"
    lp_other_dmap78     <- log(5.4e-4) ; label("'Other organs' endothelial permeability, DMAP-BLP 78 nm (cm/s)")      # Fig 5B, "P to other"
    lp_other_dmap45     <- log(3.6e-4) ; label("'Other organs' endothelial permeability, DMAP-BLP 45 nm (cm/s)")      # Fig 5C, "P to other"

    # 'Other organs' first-order elimination, the only elimination route
    # retained in the mouse model.
    lkel_other_mc3      <- log(0.6317)   ; label("Lipid elimination rate in 'other organs', MC3 80 nm and DMAP-BLP 78 nm (1/h)") # Fig 5A "Kel in other"; Fig 5B same value with footnote #
    lkel_other_dmap45   <- log(2.3447)   ; label("Lipid elimination rate in 'other organs', DMAP-BLP 45 nm (1/h)")               # Fig 5C, "Kel in other"

    # Receptor-mediated cellular uptake, Eq. (4).
    lkin_mc3            <- log(5.2395)   ; label("Receptor-mediated LNP uptake rate, MC3 80 nm and DMAP-BLP 78 nm (mL/mmol/h)")  # Fig 5A "Kin"; Fig 5B same value with footnote #
    lkin_dmap45         <- log(125.0626) ; label("Receptor-mediated LNP uptake rate, DMAP-BLP 45 nm (mL/mmol/h)")                # Fig 5C, "Kin"

    # ================================================================
    # First-order dissociation of ionizable lipid from the LNP in plasma
    # and interstitium, fitted to in vitro mouse-plasma incubation data with
    # a one-compartment model. Not measured for the MC3 80 nm arm, which is
    # therefore encoded with no dissociation.
    # ================================================================
    lkf_dmap78          <- log(0.0414)   ; label("LNP dissociation rate in blood and interstitium, DMAP-BLP 78 nm (1/h)")        # Fig 6A, "kdis_plasma" 78 nm curve
    lkf_dmap45          <- log(0.1089)   ; label("LNP dissociation rate in blood and interstitium, DMAP-BLP 45 nm (1/h)")        # Fig 6A, "kdis_plasma" 45 nm curve

    # Uptake-receptor concentration, Eq. (9).
    creceptor_liver     <- fixed(1.000)  ; label("LDL-receptor concentration in liver, assumed anchor (mmol/mL)")                # Table 1, "Receptor concentration in liver" + footnote b
    lscale_receptor     <- log(0.1217)   ; label("Receptor concentration ratio, spleen to liver (unitless)")                     # Table 1, "Receptor scaling from the liver to spleen" (fitted, footnote b)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology of the 27 g mouse (Wang 2024 Table 1, "Mice" column).
    #    Volumes mL, blood flows mL/h, endothelium surface areas cm^2.
    # ------------------------------------------------------------------
    v_venous       <- 0.5700   # Vein (mL)                          # Table 1
    v_arterial     <- 0.2900   # Artery (mL)                        # Table 1
    v_lung_vas     <- 0.05670  # Lung blood vessel (mL)             # Table 1
    v_liver_vas    <- 0.1476   # Liver blood vessel (mL)            # Table 1
    v_liver_inter  <- 0.1968   # Liver interstitium (mL)            # Table 1
    v_liver_cell   <- 0.8856   # Liver cell (mL)                    # Table 1
    v_spleen_vas   <- 0.02520  # Spleen blood vessel (mL)           # Table 1
    v_spleen_inter <- 0.01350  # Spleen interstitium (mL)           # Table 1
    v_spleen_cell  <- 0.05130  # Spleen cell (mL)                   # Table 1
    v_other_vas    <- 0.6855   # Other organs' blood vessel (mL)    # Table 1
    v_other_cell   <- 18.22    # Other organs' cell (mL)            # Table 1

    q_lung         <- 310.8    # Lung blood flow (mL/h)             # Table 1
    q_hepart       <- 114.3    # Hepatic artery blood flow (mL/h)   # Table 1 + footnote a
    q_liver        <- 119.4    # Hepatic vein blood flow (mL/h)     # Table 1
    q_spleen       <- 5.100    # Spleen blood flow (mL/h)           # Table 1
    q_other        <- 191.4    # Other organs' blood flow (mL/h)    # Table 1

    sa_liver       <- 134.5    # Liver endothelium surface area (cm^2)          # Table 1
    sa_spleen      <- 25.37    # Spleen endothelium surface area (cm^2)         # Table 1
    sa_other       <- 589.9    # Other organs' endothelium surface area (cm^2)  # Table 1
    hct            <- 0.4500   # Hematocrit (unitless)                          # Table 1

    # ------------------------------------------------------------------
    # 2. Formulation selection. Both indicators zero selects MC3 at 80 nm.
    # ------------------------------------------------------------------
    is_mc3 <- (1 - FORM_LNP_DMAPBLP78) * (1 - FORM_LNP_DMAPBLP45)

    # The 78 nm arm inherits the 80 nm liver permeability, uptake rate and
    # 'other'-organ elimination rate (Figure 5B footnote #).
    p_liver   <- exp(lp_liver_mc3) * (1 - FORM_LNP_DMAPBLP45) +
      exp(lp_liver_dmap45) * FORM_LNP_DMAPBLP45
    kel_other <- exp(lkel_other_mc3) * (1 - FORM_LNP_DMAPBLP45) +
      exp(lkel_other_dmap45) * FORM_LNP_DMAPBLP45
    kin       <- exp(lkin_mc3) * (1 - FORM_LNP_DMAPBLP45) +
      exp(lkin_dmap45) * FORM_LNP_DMAPBLP45

    p_spleen  <- exp(lp_spleen_mc3) * is_mc3 +
      exp(lp_spleen_dmap78) * FORM_LNP_DMAPBLP78 +
      exp(lp_spleen_dmap45) * FORM_LNP_DMAPBLP45
    p_other   <- exp(lp_other_mc3) * is_mc3 +
      exp(lp_other_dmap78) * FORM_LNP_DMAPBLP78 +
      exp(lp_other_dmap45) * FORM_LNP_DMAPBLP45

    # No plasma dissociation was measured for the MC3 80 nm arm, so kf is
    # zero there and the Figure 5A simulation is recovered exactly.
    kf <- exp(lkf_dmap78) * FORM_LNP_DMAPBLP78 +
      exp(lkf_dmap45) * FORM_LNP_DMAPBLP45

    creceptor_spleen <- creceptor_liver * exp(lscale_receptor)

    # ------------------------------------------------------------------
    # 3. Permeation coefficients, Eq. (3). See the rat model file for the
    #    evidence that P x S is used directly on the mL/h scale rather than
    #    converted from cm^3/s.
    # ------------------------------------------------------------------
    ps_liver  <- p_liver  * sa_liver
    ps_spleen <- p_spleen * sa_spleen
    ps_other  <- p_other  * sa_other

    # ------------------------------------------------------------------
    # 4. Sub-space concentrations (mg/mL) of LNP-associated lipid.
    # ------------------------------------------------------------------
    c_venous       <- venous_np     / v_venous
    c_arterial     <- arterial_np   / v_arterial
    c_lung_vas     <- vp_lung_np    / v_lung_vas
    c_liver_vas    <- vp_liver_np   / v_liver_vas
    c_liver_inter  <- is_liver_np   / v_liver_inter
    c_liver_cell   <- int_liver_np  / v_liver_cell
    c_spleen_vas   <- vp_spleen_np  / v_spleen_vas
    c_spleen_inter <- is_spleen_np  / v_spleen_inter
    c_spleen_cell  <- int_spleen_np / v_spleen_cell
    c_other_vas    <- vp_other_np   / v_other_vas
    c_other_cell   <- int_other_np  / v_other_cell

    # ------------------------------------------------------------------
    # 5. ODE system, Supporting Information section 2.2. Relative to the rat
    #    model: the liver and spleen intracellular pools have no disassembly
    #    and no hydrolysis (mouse data cannot resolve them), and every blood
    #    and interstitial space loses LNP-associated lipid at rate kf into a
    #    free-lipid accumulator that is never re-taken up.
    #    The lung remains a pure vascular pass-through (Figure 2).
    # ------------------------------------------------------------------
    # Blood pool
    d/dt(venous_np)   <- -q_lung * c_venous + q_other * c_other_vas +
      q_liver * c_liver_vas - kf * venous_np
    d/dt(venous)      <-  kf * venous_np
    d/dt(arterial_np) <-  q_lung * c_lung_vas - q_hepart * c_arterial -
      q_other * c_arterial - q_spleen * c_arterial - kf * arterial_np
    d/dt(arterial)    <-  kf * arterial_np
    d/dt(vp_lung_np)  <-  q_lung * c_venous - q_lung * c_lung_vas - kf * vp_lung_np
    d/dt(vp_lung)     <-  kf * vp_lung_np

    # Liver
    d/dt(vp_liver_np)  <-  q_spleen * c_spleen_vas + q_hepart * c_arterial -
      q_liver * c_liver_vas - ps_liver * (c_liver_vas - c_liver_inter) -
      kf * vp_liver_np
    d/dt(vp_liver)     <-  kf * vp_liver_np
    d/dt(is_liver_np)  <-  ps_liver * (c_liver_vas - c_liver_inter) -
      kin * c_liver_inter * creceptor_liver * v_liver_inter - kf * is_liver_np
    d/dt(is_liver)     <-  kf * is_liver_np
    d/dt(int_liver_np) <-  kin * c_liver_inter * creceptor_liver * v_liver_inter

    # Spleen
    d/dt(vp_spleen_np)  <-  q_spleen * c_arterial - q_spleen * c_spleen_vas -
      ps_spleen * (c_spleen_vas - c_spleen_inter) - kf * vp_spleen_np
    d/dt(vp_spleen)     <-  kf * vp_spleen_np
    d/dt(is_spleen_np)  <-  ps_spleen * (c_spleen_vas - c_spleen_inter) -
      kin * c_spleen_inter * creceptor_spleen * v_spleen_inter - kf * is_spleen_np
    d/dt(is_spleen)     <-  kf * is_spleen_np
    d/dt(int_spleen_np) <-  kin * c_spleen_inter * creceptor_spleen * v_spleen_inter

    # 'Other organs'
    d/dt(vp_other_np)  <- q_other * c_arterial - q_other * c_other_vas -
      ps_other * (c_other_vas - c_other_cell) - kf * vp_other_np
    d/dt(vp_other)     <- kf * vp_other_np
    d/dt(int_other_np) <- ps_other * (c_other_vas - c_other_cell) -
      kel_other * c_other_cell * v_other_cell

    # ------------------------------------------------------------------
    # 6. Observations. Cliver / Cspleen / Cblood are the quantities plotted
    #    in Figures 5A-5C and 6B; as in the rat model the liver average
    #    includes its vascular space and the spleen average excludes it.
    #    Only LNP-associated lipid contributes, matching the 3H-CHE signal.
    # ------------------------------------------------------------------
    Cliver  <- (vp_liver_np + is_liver_np + int_liver_np) /
      (v_liver_vas + v_liver_inter + v_liver_cell)
    Cspleen <- (is_spleen_np + int_spleen_np) /
      (v_spleen_inter + v_spleen_cell)
    Cblood  <- (venous_np + arterial_np + vp_lung_np + vp_liver_np +
                  vp_spleen_np + vp_other_np) /
      (v_venous + v_arterial + v_lung_vas + v_liver_vas +
         v_spleen_vas + v_other_vas)
    Cplasma <- Cblood / (1 - hct)

    # Free (dissociated) ionizable lipid accumulated across all blood and
    # interstitial spaces; zero for the MC3 80 nm arm.
    freeLipid <- venous + arterial + vp_lung + vp_liver + is_liver +
      vp_spleen + is_spleen + vp_other
  })
}
