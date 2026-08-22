Zhu_2024_borneol_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat, Sprague-Dawley). PBPK (whole-body, Phoenix WinNonlin 8.3).",
    "SPT-07A (D-borneol), a bicyclic monoterpene under development in China for acute",
    "ischaemic stroke, after intravenous bolus dosing. Fourteen mass-balance compartments",
    "(lung, heart, brain, muscle, skin, adipose, kidney, spleen, stomach, liver, intestine,",
    "rest of body, plus arterial and venous blood) built entirely from in vitro and in silico",
    "data: every tissue except adipose is perfusion-rate limited, and adipose is",
    "permeability-limited (a vascular sub-compartment `vp_adipose` exchanging with the",
    "extravascular tissue `adipose` through a permeability-surface product PS), because rat",
    "tissue-distribution studies showed high adipose accumulation. Elimination occurs in",
    "three organs -- liver, kidney and intestine -- as unbound whole-organ intrinsic",
    "clearances measured in rat liver, kidney and intestinal microsomes and scaled by",
    "MPPGL / MPPGK / MPPGI; SPT-07A is cleared almost entirely by UDP-glucuronosyltransferase",
    "glucuronidation (hepatic UGT CLint,u 2060 vs CYP 5.14 uL/min/mg protein), with the liver,",
    "kidney and intestine contributing 62.2%, 32.6% and 5.2% of systemic clearance.",
    "Tissue:plasma partition coefficients were computed by the Schmitt tissue-composition",
    "method, not fitted. This is the species the model was BUILT and validated on before being",
    "scaled to dog and human (see modellib('Zhu_2024_borneol_dog_pbpk') and",
    "modellib('Zhu_2024_borneol_human_pbpk')); PS is the single parameter estimated from in",
    "vivo data (fitted to rat adipose concentrations) and the two other species inherit it by",
    "body-weight scaling. Deterministic typical-value simulator: the paper's 5th-95th",
    "percentile bands come from 100 virtual subjects whose variability magnitude is never",
    "reported, so no IIV and no residual error are encoded (see the vignette Errata)."
  )
  reference <- paste(
    "Zhu X, Kong W, Wang Z, Liu X, Liu L. (2024).",
    "Prediction of SPT-07A Pharmacokinetics in Rats, Dogs, and Humans Using a",
    "Physiologically-Based Pharmacokinetic Model and In Vitro Data.",
    "Pharmaceutics 16(12):1596. doi:10.3390/pharmaceutics16121596. PMCID PMC11676658.",
    "Model equations: Equations (7)-(16), pages 7-8.",
    "Rat physiology (volumes, blood flows, tissue:plasma partition coefficients, Rb, fu,p):",
    "Table 1, 'Rat (0.25 kg)' columns.",
    "Unbound whole-organ intrinsic clearances: Table 2, rat rows (RLMs / RKMs / RIMs).",
    "The adipose permeability-surface product PS is NOT reported anywhere in the paper;",
    "it is recovered by digitising the paper's own rat adipose simulation (Figure 4J) --",
    "see the inline note on PS and the vignette Errata.",
    "No erratum or corrigendum was located for this article.",
    sep = " "
  )
  vignette <- "Zhu_2024_borneol_pbpk"
  dosing   <- "venous"   # IV into the venous blood pool; else buildModelDb() mislabels it

  units <- list(
    time          = "min",
    dosing        = "mg",
    concentration = "mg/L",
    amount        = "mg",
    weight        = "kg"
  )

  # Fourteen organs; adipose is split into a vascular sub-compartment and the
  # extravascular tissue because it is the one permeability-limited organ.
  # Compartment-name mapping to the paper: `vp_adipose` is the paper's C1
  # (Equation 11, "vascular"), `adipose` is the paper's C2 (Equation 12,
  # "extravascular"). `is_adipose` was deliberately NOT used for C2: the
  # registered `is_` prefix denotes interstitial FLUID, whereas for a small
  # lipophilic monoterpene with Kadipose:plasma = 1.92 the extravascular space is
  # predominantly intracellular lipid. `other` is the paper's "rest of body" (ROB).
  compartmentData <- list(
    lung        = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    heart       = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    brain       = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    muscle      = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    skin        = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    other       = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    spleen      = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    stomach     = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    vp_adipose  = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose     = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    kidney      = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    intestine   = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    liver       = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "tissue", verified = TRUE),
    arterial    = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "whole blood", verified = TRUE),
    venous      = list(analyte = "SPT-07A (D-borneol)", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  # No covariates: a fixed 0.25 kg reference rat. BW enters only the scaling of
  # the per-kilogram intrinsic clearances; organ volumes and blood flows are
  # absolute values for a 0.25 kg rat and are NOT weight-scaled by the model.
  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 258L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "0.25 kg reference rat (Table 1)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy rats (kinetic and tissue-distribution studies; no disease model)",
    dose_range     = "0.5, 1 and 2 mg/kg single IV bolus; 1 mg/kg qd IV for 7 days (multiple dose)",
    regions        = "China (China Pharmaceutical University)",
    notes          = paste(
      "Methods 2.3.1: 180 SD rats received single doses of 0.5, 1 or 2 mg/kg and a further",
      "60 rats received 1 mg/kg qd for 7 days; Methods 2.3.2: 18 additional rats received",
      "2 mg/kg for the tissue-distribution study (plasma, heart, liver, spleen, stomach,",
      "brain, intestine, kidney, muscle, skin, adipose and lung sampled at 5, 30 and 90 min,",
      "six rats per time point). 180 + 60 + 18 = 258 animals. Plasma sampling ran to 240 min",
      "with the FIRST sample at 5 min -- relevant to the Table 3 AUC discrepancy documented",
      "in the vignette Errata. Age and sex distribution are not reported; the microsomal",
      "work used mixed-gender rat liver microsomes. Animal protocol CPU-PCPK-15030016."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Body weight. Enters ONLY the per-kilogram intrinsic-clearance terms
    # (Equations 14-16 multiply CLint,in vivo,u by W). Organ volumes and
    # blood flows below are absolute for this reference animal.
    # ------------------------------------------------------------------
    BW <- fixed(0.25) ; label("Body weight (kg)")                                    # Table 1 column header 'Rat (0.25 kg)'

    # ------------------------------------------------------------------
    # Organ volumes (mL) -- Table 1, rat Volume column, transcribed verbatim.
    # Adipose 19.00 mL total is split by the permeability-limited model into
    # V1 = 2.66 mL (vascular) + V2 = 16.34 mL (extravascular); 2.66 + 16.34 =
    # 19.00 exactly. V1 / V2 are from PK-Sim 11.2 (paper text below Eq 12).
    # ------------------------------------------------------------------
    Vlung      <- fixed(1.25)   ; label("Lung volume (mL)")                          # Table 1 rat Lung Volume
    Vheart     <- fixed(0.83)   ; label("Heart volume (mL)")                         # Table 1 rat Heart Volume
    Vbrain     <- fixed(1.43)   ; label("Brain volume (mL)")                         # Table 1 rat Brain Volume
    Vmuscle    <- fixed(117.50) ; label("Muscle volume (mL)")                        # Table 1 rat Muscle Volume
    Vskin      <- fixed(47.50)  ; label("Skin volume (mL)")                          # Table 1 rat Skin Volume
    Vadipose1  <- fixed(2.66)   ; label("Adipose vascular sub-compartment volume V1 (mL)")      # paper text below Eq 12, 'for rats, they are 2.66 mL and 16.34 mL' (PK-Sim 11.2)
    Vadipose2  <- fixed(16.34)  ; label("Adipose extravascular sub-compartment volume V2 (mL)") # paper text below Eq 12; V1 + V2 = 19.00 = Table 1 rat Adipose Volume
    Vkidney    <- fixed(1.83)   ; label("Kidney volume (mL)")                        # Table 1 rat Kidney Volume
    Vspleen    <- fixed(0.50)   ; label("Spleen volume (mL)")                        # Table 1 rat Spleen Volume
    Vstomach   <- fixed(1.10)   ; label("Stomach volume (mL)")                       # Table 1 rat Stomach Volume
    Vliver     <- fixed(9.15)   ; label("Liver volume (mL)")                         # Table 1 rat Liver Volume
    Vvenous    <- fixed(13.60)  ; label("Venous blood volume (mL)")                  # Table 1 rat Vein Volume
    Varterial  <- fixed(6.80)   ; label("Arterial blood volume (mL)")                # Table 1 rat Artery Volume
    Vintestine <- fixed(10.01)  ; label("Intestine volume (mL)")                     # Table 1 rat intestine Volume
    Vother     <- fixed(19.50)  ; label("Rest-of-body volume (mL)")                  # Table 1 rat 'Rest of body' Volume

    # ------------------------------------------------------------------
    # Organ blood flows (mL/min) -- Table 1, rat Blood Flow column.
    # Qtotal is the cardiac output, printed on the Lung row. The eleven
    # tissue flows sum to 83.90 = Qtotal exactly (mass balance verified;
    # see the vignette source-trace table). Qhepart is the hepatic ARTERY
    # flow (Table 1 'Liver-art' row); the Table 1 'Liver' flow of 12.30 is
    # the TOTAL hepatic inflow = Qhepart + Qspleen + Qstomach + Qintestine
    # = 1.99 + 1.66 + 1.13 + 7.52 = 12.30, and is reconstructed in model().
    # ------------------------------------------------------------------
    Qtotal     <- fixed(83.90) ; label("Cardiac output (mL/min)")                    # Table 1 rat Lung Blood Flow (= cardiac output, Eq 8-10 Qtotal)
    Qheart     <- fixed(4.07)  ; label("Heart blood flow (mL/min)")                  # Table 1 rat Heart Blood Flow
    Qbrain     <- fixed(1.66)  ; label("Brain blood flow (mL/min)")                  # Table 1 rat Brain Blood Flow
    Qmuscle    <- fixed(8.23)  ; label("Muscle blood flow (mL/min)")                 # Table 1 rat Muscle Blood Flow
    Qadipose   <- fixed(5.82)  ; label("Adipose blood flow (mL/min)")                # Table 1 rat Adipose Blood Flow
    Qskin      <- fixed(4.82)  ; label("Skin blood flow (mL/min)")                   # Table 1 rat Skin Blood Flow
    Qkidney    <- fixed(11.71) ; label("Kidney blood flow (mL/min)")                 # Table 1 rat Kidney Blood Flow
    Qspleen    <- fixed(1.66)  ; label("Spleen blood flow (mL/min)")                 # Table 1 rat Spleen Blood Flow
    Qstomach   <- fixed(1.13)  ; label("Stomach blood flow (mL/min)")                # Table 1 rat Stomach Blood Flow
    Qhepart    <- fixed(1.99)  ; label("Hepatic artery blood flow (mL/min)")         # Table 1 rat 'Liver-art' Blood Flow (Qhep in Eq 14)
    Qintestine <- fixed(7.52)  ; label("Intestine blood flow (mL/min)")              # Table 1 rat intestine Blood Flow
    Qother     <- fixed(35.29) ; label("Rest-of-body blood flow (mL/min)")           # Table 1 rat 'Rest of body' Blood Flow

    # ------------------------------------------------------------------
    # Tissue:plasma partition coefficients (unitless) -- Table 1, rat Kt:pl
    # column. CALCULATED by the Schmitt tissue-composition method from
    # species tissue composition and fu,p (Methods 2.5.1 final paragraph),
    # not fitted -> fixed().
    # ------------------------------------------------------------------
    Kplung      <- fixed(0.59) ; label("Lung:plasma partition coefficient (unitless)")        # Table 1 rat Lung Kt:pl
    Kpheart     <- fixed(1.13) ; label("Heart:plasma partition coefficient (unitless)")       # Table 1 rat Heart Kt:pl
    Kpbrain     <- fixed(1.46) ; label("Brain:plasma partition coefficient (unitless)")       # Table 1 rat Brain Kt:pl
    Kpmuscle    <- fixed(0.90) ; label("Muscle:plasma partition coefficient (unitless)")      # Table 1 rat Muscle Kt:pl
    Kpadipose   <- fixed(1.92) ; label("Adipose:plasma partition coefficient (unitless)")     # Table 1 rat Adipose Kt:pl
    Kpskin      <- fixed(7.20) ; label("Skin:plasma partition coefficient (unitless)")        # Table 1 rat Skin Kt:pl
    Kpkidney    <- fixed(1.51) ; label("Kidney:plasma partition coefficient (unitless)")      # Table 1 rat Kidney Kt:pl
    Kpspleen    <- fixed(0.98) ; label("Spleen:plasma partition coefficient (unitless)")      # Table 1 rat Spleen Kt:pl
    Kpstomach   <- fixed(1.00) ; label("Stomach:plasma partition coefficient (unitless)")     # Table 1 rat Stomach Kt:pl
    Kpliver     <- fixed(1.53) ; label("Liver:plasma partition coefficient (unitless)")       # Table 1 rat Liver Kt:pl
    Kpintestine <- fixed(1.00) ; label("Intestine:plasma partition coefficient (unitless)")   # Table 1 rat intestine Kt:pl
    Kpother     <- fixed(1.00) ; label("Rest-of-body:plasma partition coefficient (unitless)") # Table 1 rat 'Rest of body' Kt:pl

    # ------------------------------------------------------------------
    # Blood binding and blood:plasma ratio -- Table 1 rat R_b and f_u,p rows.
    # Both are MEASURED in vitro (Methods 2.2.1 / 2.2.2; Results 3.1.1).
    # ------------------------------------------------------------------
    Rb  <- fixed(1.10)   ; label("Blood:plasma concentration ratio (unitless)")       # Table 1 rat R_b; Results 3.1.1 average of 1.11/1.09/1.11 at 100/300/900 ng/mL
    fup <- fixed(0.2870) ; label("Fraction unbound in plasma (unitless)")             # Table 1 rat f_u,p 28.70%; Results 3.1.1 (mean protein binding 71.30%)

    # ------------------------------------------------------------------
    # Unbound whole-organ in vivo intrinsic clearances (mL/min/kg body
    # weight) -- Table 2, rat rows, 'CL int, in vivo,u' column. The liver
    # value is the TOTAL column (UGT 3380 + CYP 8.43 = 3388 -> 3390 as
    # printed); kidney and intestine are UGT-only (Results 3.1.2 reports no
    # CYP-mediated metabolism in renal or intestinal microsomes).
    # ------------------------------------------------------------------
    CLintLiver     <- fixed(3390) ; label("Unbound hepatic whole-organ intrinsic clearance (mL/min/kg)")     # Table 2 rat RLMs 'Total CL int,in vivo,u' 3390 (UGT 3380 + CYP 8.43)
    CLintKidney    <- fixed(196)  ; label("Unbound renal whole-organ intrinsic clearance (mL/min/kg)")       # Table 2 rat RKMs CL int,in vivo,u 196 (UGT only)
    CLintIntestine <- fixed(17.2) ; label("Unbound intestinal whole-organ intrinsic clearance (mL/min/kg)")  # Table 2 rat RIMs CL int,in vivo,u 17.2 (UGT only)

    # ------------------------------------------------------------------
    # Adipose permeability-surface product (mL/min).
    #
    # NOT REPORTED IN THE PAPER. Methods 2.5.1 states only that "PS is the
    # permeability-surface product, which was estimated by fitting this value
    # to the adipose concentration data in rats" -- the fitted number appears
    # in no table, figure caption or supplement, and this article has no
    # supplementary material. Recovered here by digitising the paper's OWN
    # rat adipose simulation (Figure 4J, 2 mg/kg IV bolus, median line) and
    # solving for the PS that reproduces it: the optimiser returns 0.499
    # mL/min, i.e. the authors' round value 0.5. At PS = 0.5 the model
    # reproduces the digitised median to within 5% on average
    # (Tmax 12 vs ~10 min; peak 452 vs ~435; 30 min 414 vs ~410;
    # 90 min 215 vs ~205 ng/g). PS is the ONLY in-vivo-fitted parameter in
    # the whole model and the only one that is not from Table 1 or Table 2.
    # The dog and human siblings inherit it by Equation (13),
    # PS_i = PS_rat * (W_i / W_rat), a LINEAR body-weight scaling (Equation 13
    # as printed carries no exponent). See the vignette Errata.
    # ------------------------------------------------------------------
    PS <- fixed(0.5) ; label("Adipose permeability-surface product (mL/min) -- figure-derived")  # back-solved from Figure 4J (digitised); not reported in the paper

    # ------------------------------------------------------------------
    # Residual error. The paper performs SIMULATION only and never reports a
    # residual-error model, so this is fixed at zero rather than invented.
    # The 5th-95th percentile bands in Figures 3-5 come from 100 virtual
    # subjects whose parameter-variability magnitude is likewise never
    # reported, so no IIV is encoded either. See the vignette Errata.
    # ------------------------------------------------------------------
    propSd <- fixed(0) ; label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # ================================================================
    # 1. Unit handling. Table 1 volumes and flows are transcribed above in
    #    mL and mL/min exactly as printed; convert to L and L/min so that
    #    with dosing in mg the concentrations come out in mg/L.
    # ================================================================
    v_lung      <- Vlung      * 1e-3
    v_heart     <- Vheart     * 1e-3
    v_brain     <- Vbrain     * 1e-3
    v_muscle    <- Vmuscle    * 1e-3
    v_skin      <- Vskin      * 1e-3
    v_adipose1  <- Vadipose1  * 1e-3
    v_adipose2  <- Vadipose2  * 1e-3
    v_kidney    <- Vkidney    * 1e-3
    v_spleen    <- Vspleen    * 1e-3
    v_stomach   <- Vstomach   * 1e-3
    v_liver     <- Vliver     * 1e-3
    v_venous    <- Vvenous    * 1e-3
    v_arterial  <- Varterial  * 1e-3
    v_intestine <- Vintestine * 1e-3
    v_other     <- Vother     * 1e-3

    q_total     <- Qtotal     * 1e-3
    q_heart     <- Qheart     * 1e-3
    q_brain     <- Qbrain     * 1e-3
    q_muscle    <- Qmuscle    * 1e-3
    q_adipose   <- Qadipose   * 1e-3
    q_skin      <- Qskin      * 1e-3
    q_kidney    <- Qkidney    * 1e-3
    q_spleen    <- Qspleen    * 1e-3
    q_stomach   <- Qstomach   * 1e-3
    q_hepart    <- Qhepart    * 1e-3
    q_intestine <- Qintestine * 1e-3
    q_other     <- Qother     * 1e-3
    ps_adipose  <- PS         * 1e-3

    # Total hepatic inflow: hepatic artery plus the three splanchnic organs
    # that drain into the liver (Equation 14 outflow term).
    q_liver <- q_hepart + q_stomach + q_intestine + q_spleen

    # ================================================================
    # 2. Tissue concentrations (mg/L). States hold AMOUNTS (mg); the printed
    #    equations are in the form V_t * dC_t/dt = flux, and since
    #    A_t = V_t * C_t these are identical to dA_t/dt = flux.
    # ================================================================
    c_lung      <- lung       / v_lung
    c_heart     <- heart      / v_heart
    c_brain     <- brain      / v_brain
    c_muscle    <- muscle     / v_muscle
    c_skin      <- skin       / v_skin
    c_other     <- other      / v_other
    c_spleen    <- spleen     / v_spleen
    c_stomach   <- stomach    / v_stomach
    c_adipose1  <- vp_adipose / v_adipose1
    c_adipose2  <- adipose    / v_adipose2
    c_kidney    <- kidney     / v_kidney
    c_intestine <- intestine  / v_intestine
    c_liver     <- liver      / v_liver
    c_arterial  <- arterial   / v_arterial
    c_venous    <- venous     / v_venous

    # ================================================================
    # 3. Effluent BLOOD concentrations, C_t / (K_t:pl / Rb). Because K_t:pl
    #    is a tissue-to-PLASMA ratio, C_t / K_t:pl is the plasma-equivalent
    #    concentration and multiplying by Rb converts it to blood. Every
    #    blood-side state (arterial, venous, and the adipose vascular
    #    sub-compartment) therefore holds a BLOOD concentration.
    # ================================================================
    e_lung      <- c_lung      / (Kplung      / Rb)
    e_heart     <- c_heart     / (Kpheart     / Rb)
    e_brain     <- c_brain     / (Kpbrain     / Rb)
    e_muscle    <- c_muscle    / (Kpmuscle    / Rb)
    e_skin      <- c_skin      / (Kpskin      / Rb)
    e_other     <- c_other     / (Kpother     / Rb)
    e_spleen    <- c_spleen    / (Kpspleen    / Rb)
    e_stomach   <- c_stomach   / (Kpstomach   / Rb)
    e_adipose   <- c_adipose2  / (Kpadipose   / Rb)
    e_kidney    <- c_kidney    / (Kpkidney    / Rb)
    e_intestine <- c_intestine / (Kpintestine / Rb)
    e_liver     <- c_liver     / (Kpliver     / Rb)

    # ================================================================
    # 4. Elimination (Equations 14-16). The eliminating term is
    #    CLint,in vivo,u * W * fu,p * C_t / K_t:pl -- note it uses the
    #    PLASMA-equivalent concentration C_t / K_t:pl (NOT divided by Rb),
    #    so that fu,p converts it to the unbound plasma concentration
    #    driving metabolism. Encoding this exactly reproduces the paper's
    #    own well-stirred organ clearances (Results 3.1.3: liver 46.6,
    #    kidney 24.4, intestine 3.90 mL/min/kg for the rat).
    # ================================================================
    elim_liver     <- CLintLiver     * 1e-3 * BW * fup * c_liver     / Kpliver
    elim_kidney    <- CLintKidney    * 1e-3 * BW * fup * c_kidney    / Kpkidney
    elim_intestine <- CLintIntestine * 1e-3 * BW * fup * c_intestine / Kpintestine

    # ================================================================
    # 5. Perfusion-rate-limited, non-eliminating tissues (Equation 7).
    # ================================================================
    d/dt(heart)   <- q_heart   * (c_arterial - e_heart)
    d/dt(brain)   <- q_brain   * (c_arterial - e_brain)
    d/dt(muscle)  <- q_muscle  * (c_arterial - e_muscle)
    d/dt(skin)    <- q_skin    * (c_arterial - e_skin)
    d/dt(other)   <- q_other   * (c_arterial - e_other)
    d/dt(spleen)  <- q_spleen  * (c_arterial - e_spleen)
    d/dt(stomach) <- q_stomach * (c_arterial - e_stomach)

    # ================================================================
    # 6. Permeability-limited adipose (Equations 11 and 12). The vascular
    #    sub-compartment exchanges with arterial blood at q_adipose and with
    #    the extravascular tissue through the permeability-surface product;
    #    its effluent to the venous pool is c_adipose1 itself (Equation 11
    #    writes the perfusion term as q_adipose * (c_arterial - c_adipose1),
    #    i.e. the vascular sub-compartment already holds a blood
    #    concentration and carries no partition coefficient).
    # ================================================================
    d/dt(vp_adipose) <- q_adipose * (c_arterial - c_adipose1) -
                        ps_adipose * (c_adipose1 - e_adipose)
    d/dt(adipose)    <- ps_adipose * (c_adipose1 - e_adipose)

    # ================================================================
    # 7. Eliminating organs (Equations 14-16). The kidney and intestine
    #    receive arterial blood; the liver receives hepatic-artery blood
    #    plus the stomach, intestinal and splenic effluents, and drains the
    #    whole of that flow to the venous pool.
    # ================================================================
    d/dt(kidney)    <- q_kidney    * (c_arterial - e_kidney)    - elim_kidney
    d/dt(intestine) <- q_intestine * (c_arterial - e_intestine) - elim_intestine
    d/dt(liver)     <- q_hepart * c_arterial +
                       q_stomach   * e_stomach +
                       q_intestine * e_intestine +
                       q_spleen    * e_spleen -
                       q_liver     * e_liver -
                       elim_liver

    # ================================================================
    # 8. Lung, arterial and venous blood (Equations 10, 8 and 9). The lung
    #    receives the entire venous return and feeds the arterial pool. The
    #    venous sum runs over the eight compartments that drain DIRECTLY to
    #    the vein (Figure 1): the spleen, stomach and intestine drain into
    #    the liver instead and so are represented there by q_liver * e_liver.
    #    The eleven tissue flows sum to the cardiac output exactly.
    # ================================================================
    d/dt(lung)     <- q_total * (c_venous - e_lung)
    d/dt(arterial) <- q_total * (e_lung - c_arterial)
    d/dt(venous)   <- q_heart    * e_heart +
                      q_brain    * e_brain +
                      q_muscle   * e_muscle +
                      q_skin     * e_skin +
                      q_other    * e_other +
                      q_kidney   * e_kidney +
                      q_adipose  * c_adipose1 +
                      q_liver    * e_liver -
                      q_total    * c_venous

    # ================================================================
    # 9. Observation. The blood-side states hold BLOOD concentrations, so
    #    the plasma concentration reported by the paper is obtained by
    #    dividing the venous blood concentration by Rb. Dose IV into
    #    cmt = "venous".
    # ================================================================
    Cc <- c_venous / Rb
    Cc ~ prop(propSd)
  })
}
