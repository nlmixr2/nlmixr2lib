Zhu_2024_borneol_dog_pbpk <- function() {
  description <- paste(
    "Preclinical (beagle dog). PBPK (whole-body, Phoenix WinNonlin 8.3).",
    "SPT-07A (D-borneol), a bicyclic monoterpene under development in China for acute",
    "ischaemic stroke, after intravenous bolus dosing. Structurally IDENTICAL to the rat",
    "model that was built first (see modellib('Zhu_2024_borneol_rat_pbpk')): fourteen",
    "mass-balance compartments, perfusion-rate-limited distribution everywhere except a",
    "permeability-limited adipose, and elimination in liver, kidney and intestine as unbound",
    "whole-organ intrinsic clearances. What changes is (a) beagle physiology (Table 1 dog",
    "column) and (b) dog-specific in vitro intrinsic clearances (Table 2 DLMs / DKMs).",
    "This is a pure FORWARD SCALE-UP, not a refit: nothing was estimated from dog in vivo",
    "data, and the one in-vivo-fitted parameter in the whole model -- the adipose",
    "permeability-surface product PS -- is inherited from the rat by linear body-weight",
    "scaling (Equation 13), which is what makes the dog plasma predictions a genuine",
    "cross-species test. NO intestinal elimination: glucuronidation of SPT-07A was not",
    "detectable in dog intestinal microsomes (Table 2 DIMs row is '-'), so CLintIntestine is",
    "fixed at zero and the dog clears the drug through liver (87.3%) and kidney (12.7%) only.",
    "Dog hepatic glucuronidation is by far the fastest of the three species (UGT CLint,u",
    "12,200 uL/min/mg protein vs 2060 in rat and 745 in human). Deterministic typical-value",
    "simulator: no IIV and no residual error are encoded (see the vignette Errata)."
  )
  reference <- paste(
    "Zhu X, Kong W, Wang Z, Liu X, Liu L. (2024).",
    "Prediction of SPT-07A Pharmacokinetics in Rats, Dogs, and Humans Using a",
    "Physiologically-Based Pharmacokinetic Model and In Vitro Data.",
    "Pharmaceutics 16(12):1596. doi:10.3390/pharmaceutics16121596. PMCID PMC11676658.",
    "Model equations: Equations (7)-(16), pages 7-8.",
    "Dog physiology (volumes, blood flows, tissue:plasma partition coefficients, Rb, fu,p):",
    "Table 1, 'Dog (8.5 kg)' columns.",
    "Unbound whole-organ intrinsic clearances: Table 2, dog rows (DLMs / DKMs / DIMs).",
    "The adipose permeability-surface product PS is NOT reported anywhere in the paper;",
    "it is recovered from the rat (Figure 4J, digitised) and scaled to the dog by",
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
    species        = "beagle dog",
    n_subjects     = 6L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "8.5 kg reference beagle (Table 1)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy beagle dogs (kinetic study; no disease model)",
    dose_range     = "0.25, 0.5 and 1 mg/kg single IV bolus; 0.5 mg/kg qd IV for 7 days (multiple dose)",
    regions        = "China (China Pharmaceutical University)",
    notes          = paste(
      "Methods 2.3.3: six beagle dogs received single IV bolus doses of 0.25, 0.5 and",
      "1 mg/kg in a three-period cross-over design with a 3-5 day washout, then after a",
      "further 5-day washout received 1 mg/kg qd for 7 days. NOTE the internal",
      "inconsistency in the source: Methods 2.3.3 states the multiple-dose regimen was",
      "1 mg/kg qd, but Table 3 labels the dog multiple-dose row 0.5 mg/kg and reports the",
      "same predicted AUC (7.19) as the 0.5 mg/kg single dose; the Table 3 label is the one",
      "used for validation here because it is self-consistent with the reported exposure.",
      "Blood was sampled from the forearm vein to 300 min, with the FIRST sample at 5 min --",
      "relevant to the Table 3 AUC discrepancy documented in the vignette Errata.",
      "Age and sex distribution are not reported. Animal protocol CPU-PCPK-15030016."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Body weight. Enters ONLY the per-kilogram intrinsic-clearance terms
    # (Equations 14-16 multiply CLint,in vivo,u by W).
    # ------------------------------------------------------------------
    BW <- fixed(8.5) ; label("Body weight (kg)")                                     # Table 1 column header 'Dog (8.5 kg)'

    # ------------------------------------------------------------------
    # Organ volumes (mL) -- Table 1, dog Volume column, transcribed verbatim.
    # Adipose 1500 mL total splits into V1 = 210 (vascular) + V2 = 1290
    # (extravascular) = 1500 exactly; V1 / V2 from PK-Sim 11.2.
    # ------------------------------------------------------------------
    Vlung      <- fixed(85)   ; label("Lung volume (mL)")                            # Table 1 dog Lung Volume
    Vheart     <- fixed(43)   ; label("Heart volume (mL)")                           # Table 1 dog Heart Volume
    Vbrain     <- fixed(50)   ; label("Brain volume (mL)")                           # Table 1 dog Brain Volume
    Vmuscle    <- fixed(4250) ; label("Muscle volume (mL)")                          # Table 1 dog Muscle Volume
    Vskin      <- fixed(774)  ; label("Skin volume (mL)")                            # Table 1 dog Skin Volume (footnote c: from reference [38])
    Vadipose1  <- fixed(210)  ; label("Adipose vascular sub-compartment volume V1 (mL)")      # paper text below Eq 12, 'for dogs, 210 mL and 1290 mL' (PK-Sim 11.2)
    Vadipose2  <- fixed(1290) ; label("Adipose extravascular sub-compartment volume V2 (mL)") # paper text below Eq 12; V1 + V2 = 1500 = Table 1 dog Adipose Volume
    Vkidney    <- fixed(40)   ; label("Kidney volume (mL)")                          # Table 1 dog Kidney Volume (footnote b: PK-Sim 11.2)
    Vspleen    <- fixed(22)   ; label("Spleen volume (mL)")                          # Table 1 dog Spleen Volume
    Vstomach   <- fixed(24)   ; label("Stomach volume (mL)")                         # Table 1 dog Stomach Volume
    Vliver     <- fixed(213)  ; label("Liver volume (mL)")                           # Table 1 dog Liver Volume
    Vvenous    <- fixed(284)  ; label("Venous blood volume (mL)")                    # Table 1 dog Vein Volume
    Varterial  <- fixed(141)  ; label("Arterial blood volume (mL)")                  # Table 1 dog Artery Volume
    Vintestine <- fixed(203)  ; label("Intestine volume (mL)")                       # Table 1 dog intestine Volume
    Vother     <- fixed(871)  ; label("Rest-of-body volume (mL)")                    # Table 1 dog 'Rest of body' Volume

    # ------------------------------------------------------------------
    # Organ blood flows (mL/min) -- Table 1, dog Blood Flow column.
    # The Table 1 'Liver' flow of 323.33 is the TOTAL hepatic inflow and is
    # reconstructed in model() as Qhepart + Qspleen + Qstomach + Qintestine
    # = 45 + 13.33 + 10 + 255 = 323.33 (exact).
    # The eleven tissue flows sum to 1121.46 against a printed cardiac output
    # of 1120 -- a 0.13% rounding mismatch in the source table, carried
    # through unchanged rather than silently reconciled (see vignette Errata).
    # ------------------------------------------------------------------
    Qtotal     <- fixed(1120)  ; label("Cardiac output (mL/min)")                    # Table 1 dog Lung Blood Flow (= cardiac output, Eq 8-10 Qtotal)
    Qheart     <- fixed(43.3)  ; label("Heart blood flow (mL/min)")                  # Table 1 dog Heart Blood Flow
    Qbrain     <- fixed(145)   ; label("Brain blood flow (mL/min)")                  # Table 1 dog Brain Blood Flow
    Qmuscle    <- fixed(270)   ; label("Muscle blood flow (mL/min)")                 # Table 1 dog Muscle Blood Flow (footnote a: mean of PK-Sim 11.2 and [37])
    Qadipose   <- fixed(50)    ; label("Adipose blood flow (mL/min)")                # Table 1 dog Adipose Blood Flow
    Qskin      <- fixed(71.5)  ; label("Skin blood flow (mL/min)")                   # Table 1 dog Skin Blood Flow (footnote a)
    Qkidney    <- fixed(170)   ; label("Kidney blood flow (mL/min)")                 # Table 1 dog Kidney Blood Flow
    Qspleen    <- fixed(13.33) ; label("Spleen blood flow (mL/min)")                 # Table 1 dog Spleen Blood Flow
    Qstomach   <- fixed(10)    ; label("Stomach blood flow (mL/min)")                # Table 1 dog Stomach Blood Flow
    Qhepart    <- fixed(45)    ; label("Hepatic artery blood flow (mL/min)")         # Table 1 dog 'Liver-art' Blood Flow (Qhep in Eq 14)
    Qintestine <- fixed(255)   ; label("Intestine blood flow (mL/min)")              # Table 1 dog intestine Blood Flow
    Qother     <- fixed(48.33) ; label("Rest-of-body blood flow (mL/min)")           # Table 1 dog 'Rest of body' Blood Flow

    # ------------------------------------------------------------------
    # Tissue:plasma partition coefficients (unitless) -- Table 1, dog Kt:pl
    # column. Calculated by the Schmitt method, not fitted -> fixed().
    # ------------------------------------------------------------------
    Kplung      <- fixed(0.62) ; label("Lung:plasma partition coefficient (unitless)")        # Table 1 dog Lung Kt:pl
    Kpheart     <- fixed(0.76) ; label("Heart:plasma partition coefficient (unitless)")       # Table 1 dog Heart Kt:pl
    Kpbrain     <- fixed(1.44) ; label("Brain:plasma partition coefficient (unitless)")       # Table 1 dog Brain Kt:pl
    Kpmuscle    <- fixed(0.81) ; label("Muscle:plasma partition coefficient (unitless)")      # Table 1 dog Muscle Kt:pl
    Kpadipose   <- fixed(1.72) ; label("Adipose:plasma partition coefficient (unitless)")     # Table 1 dog Adipose Kt:pl
    Kpskin      <- fixed(6.46) ; label("Skin:plasma partition coefficient (unitless)")        # Table 1 dog Skin Kt:pl
    Kpkidney    <- fixed(0.79) ; label("Kidney:plasma partition coefficient (unitless)")      # Table 1 dog Kidney Kt:pl
    Kpspleen    <- fixed(0.49) ; label("Spleen:plasma partition coefficient (unitless)")      # Table 1 dog Spleen Kt:pl
    Kpstomach   <- fixed(1.00) ; label("Stomach:plasma partition coefficient (unitless)")     # Table 1 dog Stomach Kt:pl
    Kpliver     <- fixed(0.81) ; label("Liver:plasma partition coefficient (unitless)")       # Table 1 dog Liver Kt:pl
    Kpintestine <- fixed(1.00) ; label("Intestine:plasma partition coefficient (unitless)")   # Table 1 dog intestine Kt:pl
    Kpother     <- fixed(1.00) ; label("Rest-of-body:plasma partition coefficient (unitless)") # Table 1 dog 'Rest of body' Kt:pl

    # ------------------------------------------------------------------
    # Blood binding and blood:plasma ratio -- Table 1 dog R_b and f_u,p rows.
    # ------------------------------------------------------------------
    Rb  <- fixed(0.92)   ; label("Blood:plasma concentration ratio (unitless)")       # Table 1 dog R_b; Results 3.1.1 average of 0.92/0.92/0.91
    fup <- fixed(0.2577) ; label("Fraction unbound in plasma (unitless)")             # Table 1 dog f_u,p 25.77%; Results 3.1.1 (mean protein binding 74.23%)

    # ------------------------------------------------------------------
    # Unbound whole-organ in vivo intrinsic clearances (mL/min/kg body
    # weight) -- Table 2, dog rows. The liver Total column is 25,500
    # (UGT 25,500 + CYP 28.4, printed to three significant figures).
    # The INTESTINAL value is fixed at ZERO: Results 3.1.3 states "No
    # glucuronidation of SPT-07A was observed in DIMs" and the Table 2 DIMs
    # row is '-' throughout. This is a genuine species difference, not a
    # reporting gap.
    # ------------------------------------------------------------------
    CLintLiver     <- fixed(25500) ; label("Unbound hepatic whole-organ intrinsic clearance (mL/min/kg)")     # Table 2 dog DLMs 'Total CL int,in vivo,u' 25,500 (UGT 25,500 + CYP 28.4)
    CLintKidney    <- fixed(27.1)  ; label("Unbound renal whole-organ intrinsic clearance (mL/min/kg)")       # Table 2 dog DKMs CL int,in vivo,u 27.1 (UGT only)
    CLintIntestine <- fixed(0)     ; label("Unbound intestinal whole-organ intrinsic clearance (mL/min/kg)")  # Table 2 dog DIMs row is '-'; Results 3.1.3 'No glucuronidation of SPT-07A was observed in DIMs'

    # ------------------------------------------------------------------
    # Adipose permeability-surface product (mL/min).
    #
    # NOT REPORTED IN THE PAPER for any species. PS was fitted by the authors
    # to RAT adipose data only (Methods 2.5.1) and other species are obtained
    # from Equation (13), PS_i = PS_rat * (W_i / W_rat) -- a LINEAR
    # body-weight scaling; Equation 13 as printed carries no exponent.
    # The rat value 0.5 mL/min was recovered by digitising the paper's own
    # rat adipose simulation (Figure 4J); see the rat sibling for the fit
    # diagnostics. Scaled here: 0.5 * (8.5 / 0.25) = 17.0 mL/min.
    # See the vignette Errata.
    # ------------------------------------------------------------------
    PS <- fixed(17.0) ; label("Adipose permeability-surface product (mL/min) -- figure-derived")  # Eq 13 from the rat value back-solved from Figure 4J: 0.5 * (8.5/0.25)

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
    #    organ clearances (Results 3.1.3: dog liver 37.8, kidney 5.50
    #    mL/min/kg). The intestinal term is identically zero for the dog.
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
    #    cmt = "venous".
    # ================================================================
    Cc <- c_venous / Rb
    Cc ~ prop(propSd)
  })
}
