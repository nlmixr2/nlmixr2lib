Zhu_2024_borneol_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, Phoenix WinNonlin 8.3). SPT-07A (D-borneol), a bicyclic monoterpene",
    "under development in China for acute ischaemic stroke, in healthy adults after",
    "intravenous infusion. Structurally IDENTICAL to the rat model the framework was built on",
    "(see modellib('Zhu_2024_borneol_rat_pbpk')): fourteen mass-balance compartments,",
    "perfusion-rate-limited distribution everywhere except a permeability-limited adipose,",
    "and elimination in liver, kidney and intestine as unbound whole-organ intrinsic",
    "clearances. What changes is (a) human physiology (Table 1 human column) and (b)",
    "human-specific in vitro intrinsic clearances (Table 2 HLMs / HKMs / HIMs). This is a",
    "pure BOTTOM-UP PREDICTION: no human in vivo data were used to build or fit it -- the",
    "clinical data it is compared against are cited from a separate first-in-human study --",
    "and the one in-vivo-fitted parameter, the adipose permeability-surface product PS, is",
    "inherited from the rat by linear body-weight scaling (Equation 13). SPT-07A is cleared",
    "predominantly by UGT glucuronidation (hepatic UGT CLint,u 745 vs CYP 3.35 uL/min/mg",
    "protein), with liver, kidney and intestine contributing 76.5%, 23.1% and 0.4% of",
    "systemic clearance; UGT phenotyping identified UGT1A1 and UGT2B7 as the responsible",
    "isoforms, though a relative-activity-factor reconstruction recovers only ~36% of the",
    "measured hepatic activity, implying additional uncharacterised UGTs. The notably high",
    "renal contribution is consistent with the paper's opening observation that the reported",
    "human clearance (1942 mL/min) exceeds hepatic blood flow (1450 mL/min). Deterministic",
    "typical-value simulator: no IIV and no residual error are encoded (see the vignette",
    "Errata)."
  )
  reference <- paste(
    "Zhu X, Kong W, Wang Z, Liu X, Liu L. (2024).",
    "Prediction of SPT-07A Pharmacokinetics in Rats, Dogs, and Humans Using a",
    "Physiologically-Based Pharmacokinetic Model and In Vitro Data.",
    "Pharmaceutics 16(12):1596. doi:10.3390/pharmaceutics16121596. PMCID PMC11676658.",
    "Model equations: Equations (7)-(16), pages 7-8.",
    "Human physiology (volumes, blood flows, tissue:plasma partition coefficients, Rb,",
    "fu,p): Table 1, 'Human (70 kg)' columns.",
    "Unbound whole-organ intrinsic clearances: Table 2, human rows (HLMs / HKMs / HIMs).",
    "The clinical observations used for validation are NOT original to this paper: Methods",
    "2.3.4 states human pharmacokinetic data were cited from the paper's reference [7].",
    "The adipose permeability-surface product PS is NOT reported anywhere in the paper;",
    "it is recovered from the rat (Figure 4J, digitised) and scaled to human by",
    "Equation (13) -- see the inline note on PS and the vignette Errata.",
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

  # See the rat sibling for the compartment-name mapping to the paper's
  # notation: `vp_adipose` is the paper's C1 (Equation 11, "vascular"),
  # `adipose` is C2 (Equation 12, "extravascular"), `other` is "rest of body".
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

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "70 kg reference adult (Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (study conducted in China)",
    disease_state  = "Healthy volunteers",
    dose_range     = "10, 20 and 40 mg as a 1 h IV infusion; single dose on day 0, then q12h for 7 days from 48 h",
    regions        = "China",
    notes          = paste(
      "Methods 2.3.4: the clinical data are NOT original to this paper -- they are cited",
      "from the paper's reference [7]. Thirty-six healthy volunteers were divided evenly",
      "into three dose cohorts (10, 20 and 40 mg, i.e. 12 per cohort). SPT-07A was given as",
      "a 1 h intravenous infusion: a single dose on day 0, then after 48 h a multiple-dose",
      "phase q12h for 7 days with the last dose on the morning of day 9. Because dosing is",
      "an infusion rather than a bolus there is no early distribution spike, which is why",
      "the human predictions reproduce Table 3 closely while the rat and dog rows do not",
      "(see the vignette Errata). Age, sex and body-weight distributions are not reported in",
      "this paper; the 70 kg reference adult is the Table 1 physiology, not a study median."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Body weight. Enters ONLY the per-kilogram intrinsic-clearance terms
    # (Equations 14-16 multiply CLint,in vivo,u by W).
    # ------------------------------------------------------------------
    BW <- fixed(70) ; label("Body weight (kg)")                                      # Table 1 column header 'Human (70 kg)'

    # ------------------------------------------------------------------
    # Organ volumes (mL) -- Table 1, human Volume column, transcribed
    # verbatim. Adipose 10,000 mL total splits into V1 = 1800 (vascular) +
    # V2 = 8200 (extravascular) = 10,000 exactly; V1 / V2 from PK-Sim 11.2.
    # ------------------------------------------------------------------
    Vlung      <- fixed(1170)  ; label("Lung volume (mL)")                           # Table 1 human Lung Volume
    Vheart     <- fixed(310)   ; label("Heart volume (mL)")                          # Table 1 human Heart Volume
    Vbrain     <- fixed(1450)  ; label("Brain volume (mL)")                          # Table 1 human Brain Volume
    Vmuscle    <- fixed(35000) ; label("Muscle volume (mL)")                         # Table 1 human Muscle Volume
    Vskin      <- fixed(7800)  ; label("Skin volume (mL)")                           # Table 1 human Skin Volume
    Vadipose1  <- fixed(1800)  ; label("Adipose vascular sub-compartment volume V1 (mL)")      # paper text below Eq 12, 'for humans, 1800 mL and 8200 mL' (PK-Sim 11.2)
    Vadipose2  <- fixed(8200)  ; label("Adipose extravascular sub-compartment volume V2 (mL)") # paper text below Eq 12; V1 + V2 = 10,000 = Table 1 human Adipose Volume
    Vkidney    <- fixed(280)   ; label("Kidney volume (mL)")                         # Table 1 human Kidney Volume
    Vspleen    <- fixed(190)   ; label("Spleen volume (mL)")                         # Table 1 human Spleen Volume
    Vstomach   <- fixed(160)   ; label("Stomach volume (mL)")                        # Table 1 human Stomach Volume
    Vliver     <- fixed(1690)  ; label("Liver volume (mL)")                          # Table 1 human Liver Volume
    Vvenous    <- fixed(3470)  ; label("Venous blood volume (mL)")                   # Table 1 human Vein Volume
    Varterial  <- fixed(1730)  ; label("Arterial blood volume (mL)")                 # Table 1 human Artery Volume
    Vintestine <- fixed(1650)  ; label("Intestine volume (mL)")                      # Table 1 human intestine Volume
    Vother     <- fixed(5100)  ; label("Rest-of-body volume (mL)")                   # Table 1 human 'Rest of body' Volume

    # ------------------------------------------------------------------
    # Organ blood flows (mL/min) -- Table 1, human Blood Flow column.
    # The Table 1 'Liver' flow of 1518.33 is the TOTAL hepatic inflow and is
    # reconstructed in model() as Qhepart + Qspleen + Qstomach + Qintestine
    # = 300 + 80 + 38.33 + 1100 = 1518.33 (exact). The eleven tissue flows
    # sum to 5600.33 against a printed cardiac output of 5600 -- a 0.006%
    # rounding difference in the source table, carried through unchanged.
    # ------------------------------------------------------------------
    Qtotal     <- fixed(5600)  ; label("Cardiac output (mL/min)")                    # Table 1 human Lung Blood Flow (= cardiac output, Eq 8-10 Qtotal)
    Qheart     <- fixed(240)   ; label("Heart blood flow (mL/min)")                  # Table 1 human Heart Blood Flow
    Qbrain     <- fixed(700)   ; label("Brain blood flow (mL/min)")                  # Table 1 human Brain Blood Flow
    Qmuscle    <- fixed(750)   ; label("Muscle blood flow (mL/min)")                 # Table 1 human Muscle Blood Flow
    Qadipose   <- fixed(260)   ; label("Adipose blood flow (mL/min)")                # Table 1 human Adipose Blood Flow
    Qskin      <- fixed(300)   ; label("Skin blood flow (mL/min)")                   # Table 1 human Skin Blood Flow
    Qkidney    <- fixed(1240)  ; label("Kidney blood flow (mL/min)")                 # Table 1 human Kidney Blood Flow
    Qspleen    <- fixed(80)    ; label("Spleen blood flow (mL/min)")                 # Table 1 human Spleen Blood Flow
    Qstomach   <- fixed(38.33) ; label("Stomach blood flow (mL/min)")                # Table 1 human Stomach Blood Flow
    Qhepart    <- fixed(300)   ; label("Hepatic artery blood flow (mL/min)")         # Table 1 human 'Liver-art' Blood Flow (Qhep in Eq 14)
    Qintestine <- fixed(1100)  ; label("Intestine blood flow (mL/min)")              # Table 1 human intestine Blood Flow
    Qother     <- fixed(592)   ; label("Rest-of-body blood flow (mL/min)")           # Table 1 human 'Rest of body' Blood Flow

    # ------------------------------------------------------------------
    # Tissue:plasma partition coefficients (unitless) -- Table 1, human
    # Kt:pl column. Calculated by the Schmitt method, not fitted -> fixed().
    # ------------------------------------------------------------------
    Kplung      <- fixed(0.76) ; label("Lung:plasma partition coefficient (unitless)")        # Table 1 human Lung Kt:pl
    Kpheart     <- fixed(3.62) ; label("Heart:plasma partition coefficient (unitless)")       # Table 1 human Heart Kt:pl
    Kpbrain     <- fixed(4.16) ; label("Brain:plasma partition coefficient (unitless)")       # Table 1 human Brain Kt:pl
    Kpmuscle    <- fixed(1.19) ; label("Muscle:plasma partition coefficient (unitless)")      # Table 1 human Muscle Kt:pl
    Kpadipose   <- fixed(2.84) ; label("Adipose:plasma partition coefficient (unitless)")     # Table 1 human Adipose Kt:pl
    Kpskin      <- fixed(8.18) ; label("Skin:plasma partition coefficient (unitless)")        # Table 1 human Skin Kt:pl
    Kpkidney    <- fixed(2.22) ; label("Kidney:plasma partition coefficient (unitless)")      # Table 1 human Kidney Kt:pl
    Kpspleen    <- fixed(1.27) ; label("Spleen:plasma partition coefficient (unitless)")      # Table 1 human Spleen Kt:pl
    Kpstomach   <- fixed(1.00) ; label("Stomach:plasma partition coefficient (unitless)")     # Table 1 human Stomach Kt:pl
    Kpliver     <- fixed(2.06) ; label("Liver:plasma partition coefficient (unitless)")       # Table 1 human Liver Kt:pl
    Kpintestine <- fixed(1.00) ; label("Intestine:plasma partition coefficient (unitless)")   # Table 1 human intestine Kt:pl
    Kpother     <- fixed(1.00) ; label("Rest-of-body:plasma partition coefficient (unitless)") # Table 1 human 'Rest of body' Kt:pl

    # ------------------------------------------------------------------
    # Blood binding and blood:plasma ratio -- Table 1 human R_b / f_u,p rows.
    # ------------------------------------------------------------------
    Rb  <- fixed(0.92)   ; label("Blood:plasma concentration ratio (unitless)")       # Table 1 human R_b; Results 3.1.1 average of 0.90/0.93/0.92
    fup <- fixed(0.3264) ; label("Fraction unbound in plasma (unitless)")             # Table 1 human f_u,p 32.64%; Results 3.1.1 (mean protein binding 67.36%)

    # ------------------------------------------------------------------
    # Unbound whole-organ in vivo intrinsic clearances (mL/min/kg body
    # weight) -- Table 2, human rows. The liver Total column is 938
    # (UGT 934 + CYP 4.20 = 938.2). Kidney and intestine are UGT-only.
    # ------------------------------------------------------------------
    CLintLiver     <- fixed(938)   ; label("Unbound hepatic whole-organ intrinsic clearance (mL/min/kg)")     # Table 2 human HLMs 'Total CL int,in vivo,u' 938 (UGT 934 + CYP 4.20)
    CLintKidney    <- fixed(26.5)  ; label("Unbound renal whole-organ intrinsic clearance (mL/min/kg)")       # Table 2 human HKMs CL int,in vivo,u 26.5 (UGT only)
    CLintIntestine <- fixed(0.309) ; label("Unbound intestinal whole-organ intrinsic clearance (mL/min/kg)")  # Table 2 human HIMs CL int,in vivo,u 0.309 (UGT only)

    # ------------------------------------------------------------------
    # Adipose permeability-surface product (mL/min).
    #
    # NOT REPORTED IN THE PAPER for any species. PS was fitted by the authors
    # to RAT adipose data only (Methods 2.5.1) and other species are obtained
    # from Equation (13), PS_i = PS_rat * (W_i / W_rat) -- a LINEAR
    # body-weight scaling; Equation 13 as printed carries no exponent.
    # The rat value 0.5 mL/min was recovered by digitising the paper's own
    # rat adipose simulation (Figure 4J); see the rat sibling for the fit
    # diagnostics. Scaled here: 0.5 * (70 / 0.25) = 140 mL/min.
    # See the vignette Errata.
    # ------------------------------------------------------------------
    PS <- fixed(140) ; label("Adipose permeability-surface product (mL/min) -- figure-derived")  # Eq 13 from the rat value back-solved from Figure 4J: 0.5 * (70/0.25)

    # ------------------------------------------------------------------
    # Residual error -- simulation-only source; see the rat sibling.
    # ------------------------------------------------------------------
    propSd <- fixed(0) ; label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # ================================================================
    # 1. Unit handling: Table 1 mL and mL/min -> L and L/min, so that with
    #    dosing in mg the concentrations come out in mg/L.
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

    q_liver <- q_hepart + q_stomach + q_intestine + q_spleen

    # ================================================================
    # 2. Tissue concentrations (mg/L). States hold AMOUNTS (mg); the printed
    #    V_t * dC_t/dt = flux is identical to dA_t/dt = flux.
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
    # 3. Effluent BLOOD concentrations, C_t / (K_t:pl / Rb).
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
    # 4. Elimination (Equations 14-16), driven by the PLASMA-equivalent
    #    concentration C_t / K_t:pl. Reproduces the paper's own well-stirred
    #    organ clearances (Results 3.1.3: human liver 20.4, kidney 6.14,
    #    intestine 0.109 mL/min/kg).
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
    # 6. Permeability-limited adipose (Equations 11 and 12).
    # ================================================================
    d/dt(vp_adipose) <- q_adipose * (c_arterial - c_adipose1) -
                        ps_adipose * (c_adipose1 - e_adipose)
    d/dt(adipose)    <- ps_adipose * (c_adipose1 - e_adipose)

    # ================================================================
    # 7. Eliminating organs (Equations 14-16).
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
    # 8. Lung, arterial and venous blood (Equations 10, 8 and 9).
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
    # 9. Observation: venous PLASMA concentration. Dose IV into
    #    cmt = "venous" (1 h infusion via rate = amt / 60).
    # ================================================================
    Cc <- c_venous / Rb
    Cc ~ prop(propSd)
  })
}
