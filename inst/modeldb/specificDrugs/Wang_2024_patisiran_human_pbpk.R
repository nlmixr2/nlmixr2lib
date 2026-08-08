Wang_2024_patisiran_human_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic six-compartment, MATLAB SimBiology 2022a) model for",
    "the plasma disposition of the ionizable lipid DLin-MC3-DMA (MC3) after",
    "intravenous patisiran (Onpattro) in healthy volunteers. The rat / mouse",
    "six-compartment skeleton (vein, artery, lung blood vessel, liver, spleen,",
    "'other organs') further reduced for solvability: LNP entry into the spleen",
    "is switched off, so the spleen is a perfused pass-through, and the",
    "elimination route in the 'other organs' compartment is switched off, so",
    "that compartment acts as a slowly equilibrating deep depot that generates",
    "the long terminal phase. As in the mouse model there is no LNP",
    "disassembly and no ionizable-lipid hydrolysis. Fitted to the 0.5 mg/kg",
    "arm of the ALN-TTR02-001 and ALN-TTR02-005 phase 1 trials and then used,",
    "with those parameters held, to predict the 0.01-0.3 mg/kg arms.",
    sep = " "
  )
  reference <- paste(
    "Wang W, Deng S, Lin J, Ouyang D (2024).",
    "Modeling on in vivo disposition and cellular transportation of RNA lipid",
    "nanoparticles via quantum mechanics/physiologically-based pharmacokinetic",
    "approaches. Acta Pharm Sin B 14(10):4591-4607.",
    "doi:10.1016/j.apsb.2024.06.011.",
    "Model equations Eqs. (1)-(4) and the SimBiology ODE export in the",
    "Supporting Information (mmc1.pdf) section 1; physiology from Table 1",
    "('Human' column); fitted parameters from Figure 7. The observed plasma",
    "MC3 data are from the patisiran phase 1 review document Wang 2024 cites",
    "as reference 26.",
    sep = " "
  )
  vignette <- "Wang_2024_rnaLipidNanoparticle"
  units <- list(time = "h", dosing = "mg", concentration = "mg/mL")

  compartmentData <- list(
    venous_np    = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_np  = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_lung_np   = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    vp_liver_np  = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_liver_np  = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver_np = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_spleen_np = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    is_spleen_np = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_other_np  = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "whole blood", verified = TRUE),
    int_other_np = list(analyte = "DLin-MC3-DMA (in LNP)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = 2L,
    weight_range  = "69.70 kg reference body weight (Table 1)",
    disease_state = "healthy volunteers",
    dose_range    = "0.01, 0.05, 0.15, 0.3 and 0.5 mg/kg patisiran (siRNA) as a single intravenous dose",
    regions       = "not reported",
    notes         = paste(
      "Mean plasma MC3 concentrations digitised from the review document for",
      "the phase 1 trials ALN-TTR02-001 and ALN-TTR02-005 of patisiran",
      "(Onpattro), the first approved siRNA-LNP drug. The 0.5 mg/kg arm was",
      "fitted; the remaining four dose levels were then predicted with the",
      "fitted parameters held fixed and reproduced the observations for at",
      "least 1000 h. Fitted to pooled mean profiles, so the model carries no",
      "between-subject variability."
    )
  )

  ini({
    # ================================================================
    # The only three parameters fitted to the human data (Figure 7),
    # reported in nominal cm/s (permeabilities) and mL/mmol/h (uptake).
    # The permeability x surface-area product is used directly on the mL/h
    # scale; see the rat model file for the sensitivity-analysis evidence.
    # ================================================================
    lp_liver        <- log(5.7e-4)   ; label("Liver endothelial permeability (cm/s)")                        # Fig 7, "P to liver"
    lp_other        <- log(3.6e-5)   ; label("'Other organs' endothelial permeability (cm/s)")               # Fig 7, "P to other"
    lkin            <- log(10.8241)  ; label("Receptor-mediated LNP uptake rate (mL/mmol/h)")                # Fig 7, "Kin"

    # ================================================================
    # Switched-off routes. Section 3.1.4: "the model structure was further
    # reduced by shutting off the permeation to the spleen and elimination
    # in the 'other' organ." Both are encoded as exact zeros rather than
    # deleted so that the human file stays structurally comparable to its
    # rat and mouse siblings.
    # ================================================================
    p_spleen        <- fixed(0)      ; label("Spleen endothelial permeability, route switched off (cm/s)")   # Section 3.1.4, "shutting off the LNP entry into the spleen"
    kel_other       <- fixed(0)      ; label("Lipid elimination rate in 'other organs', route switched off (1/h)") # Section 3.1.4, "the elimination route in the 'other' compartment"

    # Uptake-receptor anchor. Table 1 marks both receptor rows "not
    # applicable" for the human column because the spleen route - and hence
    # the liver-to-spleen receptor ratio - is switched off. The liver anchor
    # itself keeps the assumed value of 1 mmol/mL used for rat and mouse; it
    # is an arbitrary scale that is absorbed into the fitted kin.
    creceptor_liver <- fixed(1.000)  ; label("LDL-receptor concentration in liver, assumed anchor (mmol/mL)") # Table 1 footnote b (rat / mouse columns); "not applicable" in the human column
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology of the 69.70 kg adult (Wang 2024 Table 1, "Human"
    #    column). Volumes mL, blood flows mL/h, surface areas cm^2.
    # ------------------------------------------------------------------
    v_venous       <- 1480     # Vein (mL)                          # Table 1
    v_arterial     <- 940.0    # Artery (mL)                        # Table 1
    v_lung_vas     <- 696.0    # Lung blood vessel (mL)             # Table 1
    v_liver_vas    <- 399.5    # Liver blood vessel (mL)            # Table 1
    v_liver_inter  <- 376.0    # Liver interstitium (mL)            # Table 1
    v_liver_cell   <- 1574     # Liver cell (mL)                    # Table 1
    v_spleen_vas   <- 69.30    # Spleen blood vessel (mL)           # Table 1
    v_spleen_inter <- 31.50    # Spleen interstitium (mL)           # Table 1
    v_other_vas    <- 2055     # Other organs' blood vessel (mL)    # Table 1
    v_other_cell   <- 61710    # Other organs' cell (mL)            # Table 1

    q_lung         <- 360000   # Lung blood flow (mL/h)             # Table 1
    q_hepart       <- 88200    # Hepatic artery blood flow (mL/h)   # Table 1 + footnote a
    q_liver        <- 98400    # Hepatic vein blood flow (mL/h)     # Table 1
    q_spleen       <- 10200    # Spleen blood flow (mL/h)           # Table 1
    q_other        <- 261600   # Other organs' blood flow (mL/h)    # Table 1

    sa_liver       <- 378900   # Liver endothelium surface area (cm^2)          # Table 1
    sa_spleen      <- 64650    # Spleen endothelium surface area (cm^2)         # Table 1
    sa_other       <- 1805000  # Other organs' endothelium surface area (cm^2)  # Table 1
    hct            <- 0.4700   # Hematocrit (unitless)                          # Table 1

    p_liver <- exp(lp_liver)
    p_other <- exp(lp_other)
    kin     <- exp(lkin)

    ps_liver  <- p_liver  * sa_liver
    ps_spleen <- p_spleen * sa_spleen
    ps_other  <- p_other  * sa_other

    # ------------------------------------------------------------------
    # 2. Sub-space concentrations (mg/mL).
    # ------------------------------------------------------------------
    c_venous      <- venous_np    / v_venous
    c_arterial    <- arterial_np  / v_arterial
    c_lung_vas    <- vp_lung_np   / v_lung_vas
    c_liver_vas   <- vp_liver_np  / v_liver_vas
    c_liver_inter <- is_liver_np  / v_liver_inter
    c_spleen_vas   <- vp_spleen_np / v_spleen_vas
    c_spleen_inter <- is_spleen_np / v_spleen_inter
    c_other_vas   <- vp_other_np  / v_other_vas
    c_other_cell  <- int_other_np / v_other_cell

    # ------------------------------------------------------------------
    # 3. ODE system. Structurally this is the rat system with p_spleen,
    #    kel_other, the LNP disassembly rate kdis and the free-lipid
    #    hydrolysis rates kel all set to zero, which is exactly the
    #    reduction the paper describes. With p_spleen = 0 the spleen keeps
    #    its blood flow (artery to spleen vasculature to liver vasculature)
    #    but no lipid crosses into spleen tissue, so is_spleen_np stays at
    #    zero for the whole simulation; it is retained so the human file
    #    remains structurally comparable to its rat and mouse siblings.
    #    With kel_other = 0 the 'other organs' cell pool is a reversible
    #    deep depot rather than an elimination route, and it is that depot
    #    that produces the multi-thousand-hour terminal phase in Figure 7.
    #    The lung is a pure vascular pass-through (Figure 2).
    # ------------------------------------------------------------------
    d/dt(venous_np)   <- -q_lung * c_venous + q_other * c_other_vas + q_liver * c_liver_vas
    d/dt(arterial_np) <-  q_lung * c_lung_vas - q_hepart * c_arterial -
      q_other * c_arterial - q_spleen * c_arterial
    d/dt(vp_lung_np)  <-  q_lung * c_venous - q_lung * c_lung_vas

    d/dt(vp_liver_np)  <-  q_spleen * c_spleen_vas + q_hepart * c_arterial -
      q_liver * c_liver_vas - ps_liver * (c_liver_vas - c_liver_inter)
    d/dt(is_liver_np)  <-  ps_liver * (c_liver_vas - c_liver_inter) -
      kin * c_liver_inter * creceptor_liver * v_liver_inter
    d/dt(int_liver_np) <-  kin * c_liver_inter * creceptor_liver * v_liver_inter

    d/dt(vp_spleen_np) <-  q_spleen * c_arterial - q_spleen * c_spleen_vas -
      ps_spleen * (c_spleen_vas - c_spleen_inter)
    d/dt(is_spleen_np) <-  ps_spleen * (c_spleen_vas - c_spleen_inter)

    d/dt(vp_other_np)  <-  q_other * c_arterial - q_other * c_other_vas -
      ps_other * (c_other_vas - c_other_cell)
    d/dt(int_other_np) <-  ps_other * (c_other_vas - c_other_cell) -
      kel_other * c_other_cell * v_other_cell

    # ------------------------------------------------------------------
    # 4. Observations. Figure 7 plots the plasma MC3 concentration; its
    #    y-axis is labelled mg/mL but the plotted values are ng/mL (see the
    #    vignette Errata), so Cplasma here is in mg/mL and must be
    #    multiplied by 1e6 to overlay Figure 7 directly. The SimBiology
    #    export makes this explicit with
    #    plasma_drug_ngperml = Plasma.Drug_mg_mL * 1e6.
    # ------------------------------------------------------------------
    Cblood  <- (venous_np + arterial_np + vp_lung_np + vp_liver_np +
                  vp_spleen_np + vp_other_np) /
      (v_venous + v_arterial + v_lung_vas + v_liver_vas +
         v_spleen_vas + v_other_vas)
    Cplasma <- Cblood / (1 - hct)
    Cliver  <- (vp_liver_np + is_liver_np + int_liver_np) /
      (v_liver_vas + v_liver_inter + v_liver_cell)
  })
}
