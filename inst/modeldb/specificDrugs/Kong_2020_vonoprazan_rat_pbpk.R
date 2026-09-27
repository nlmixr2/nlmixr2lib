Kong_2020_vonoprazan_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat). PBPK-PD (whole-body, perfusion-limited, WinNonlin",
    "8.1). Vonoprazan disposition and gastric-acid antisecretory effect in a",
    "typical 0.25-kg Sprague-Dawley rat after intravenous or oral dosing",
    "(Kong et al. 2020). Tissue-to-plasma partition coefficients were",
    "measured in rats; hepatic Michaelis-Menten metabolism of unbound drug",
    "is scaled from rat liver microsomes, and a linear 'other' clearance",
    "(extrahepatic metabolism) acts on venous blood. The stomach wall (the",
    "target organ) is permeability-limited with a vascular and an",
    "extravascular space, its permeability-surface product fitted to rat",
    "stomach concentrations. Oral doses enter the stomach lumen and transit",
    "a five-segment gut lumen into per-segment gut-wall compartments. The",
    "PD layer drives an H+/K+-ATPase inhibition state from the unbound",
    "extravascular stomach concentration and reports gastric perfusate pH",
    "as basal pH + I. Deterministic typical-value model with no IIV.",
    "Companion models: modellib('Kong_2020_vonoprazan_dog_pbpk') and",
    "modellib('Kong_2020_vonoprazan_human_pbpk')."
  )
  reference <- paste(
    "Kong WM, Sun BB, Wang ZJ, Zheng XK, Zhao KJ, Chen Y, Zhang JX, Liu PH,",
    "Zhu L, Xu RJ, Li P, Liu L, Liu XD. Physiologically based",
    "pharmacokinetic-pharmacodynamic modeling for prediction of vonoprazan",
    "pharmacokinetics and its inhibition on gastric acid secretion following",
    "intravenous/oral administration to rats, dogs and humans. Acta",
    "Pharmacol Sin. 2020;41(6):852-865. doi:10.1038/s41401-019-0353-2.",
    "Equations 2-22 of the paper; physiology from Table 1; drug-specific",
    "values from the Methods and Results text; ODE routing from the",
    "WinNonlin 'Mode code of human following oral administration' listing",
    "in the Supplementary Information."
  )
  vignette <- "Kong_2020_vonoprazan_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c(
    "ev_stomach",
    "gw_duodenum",
    "gw_jejunum",
    "gw_ileum",
    "gw_cecum",
    "gw_colon",
    "inhibition"
  )

  compartmentData <- list(
    stomach = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    ileum = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    cecum = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    colon = list(analyte = "vonoprazan", units = "mg", specimen = "administration site", verified = TRUE),
    gw_duodenum = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    gw_jejunum = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    gw_ileum = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    gw_cecum = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    gw_colon = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    vp_stomach = list(analyte = "vonoprazan", units = "mg", specimen = "whole blood", verified = TRUE),
    ev_stomach = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "vonoprazan", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "vonoprazan", units = "mg", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "vonoprazan", units = "mg", specimen = "whole blood", verified = TRUE),
    inhibition = list(analyte = "H+/K+-ATPase inhibition", units = "pH units", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 40L,
    n_studies = 1L,
    weight_range = "200-250 g; 0.25 kg reference rat used for the physiology (Table 1)",
    sex_female_pct = 100,
    disease_state = "healthy animals",
    dose_range = paste(
      "Intravenous vonoprazan acetate 0.5, 1 and 2 mg/kg (base) single dose",
      "and 1 mg/kg once daily for 7 days; oral vonoprazan fumarate 2 mg/kg",
      "(base); tissue distribution after 1 mg/kg intravenous"
    ),
    regions = "China",
    notes = paste(
      "Twenty-five rats in five pharmacokinetic groups (n = 5) and 15 rats",
      "in a three-time-point tissue-distribution study. Gastric pH",
      "observations (histamine-stimulated anaesthetised rats, 0.5-1.0 mg/kg",
      "intravenous) were taken from the literature (paper ref 12) for",
      "comparison only. Bottom-up prediction model: only the stomach",
      "permeability-surface product was fitted to rat data."
    )
  )

  ini({
    # ---- Physiology: 0.25-kg rat (Table 1, 'Rat (0.25 kg)' columns) -------
    # Volumes in mL, blood flows in mL/min. The supplement WinNonlin listing
    # carries the same values (fixef tvV*, tvQ*).
    v_lung    <- fixed(1.25);    label("Lung volume (mL)")          # Table 1 Lungs
    v_heart   <- fixed(0.83);     label("Heart volume (mL)")         # Table 1 Heart
    v_brain   <- fixed(1.43);    label("Brain volume (mL)")         # Table 1 Brain
    v_muscle  <- fixed(101.0);   label("Muscle volume (mL)")        # Table 1 Muscle
    v_adipose <- fixed(19.0);   label("Adipose volume (mL)")       # Table 1 Adipose
    v_skin    <- fixed(47.5);    label("Skin volume (mL)")          # Table 1 Skin
    v_kidney  <- fixed(1.83);     label("Kidney volume (mL)")        # Table 1 Kidneys
    v_spleen  <- fixed(0.5);     label("Spleen volume (mL)")        # Table 1 Spleen
    v_liver   <- fixed(9.15);    label("Liver volume (mL)")         # Table 1 Liver
    v_venous  <- fixed(13.6);    label("Venous blood volume (mL)")  # Table 1 Vein
    v_arterial <- fixed(6.8);   label("Arterial blood volume (mL)") # Table 1 Artery
    v_other   <- fixed(36.03);    label("Rest-of-body volume (mL)")  # Table 1 Rest of body

    # Stomach wall: vascular (V1) and extravascular (V2) spaces, Eqs 12-13.
    v_vp_stomach <- fixed(0.2);  label("Stomach vascular volume V1 (mL)")      # Methods 'PBPK model development': rat V1 0.2 mL
    v_ev_stomach <- fixed(0.9); label("Stomach extravascular volume V2 (mL)") # Methods: rat V2 0.9 mL (V1 + V2 = Table 1 Stomach volume)

    # Gut-wall (enterocyte) segments, Eq 20.
    v_gw_duodenum <- fixed(0.48);   label("Duodenum wall volume (mL)")  # Table 1 Duodenum
    v_gw_jejunum  <- fixed(5.36);  label("Jejunum wall volume (mL)")   # Table 1 Jejunum
    v_gw_ileum    <- fixed(0.12);  label("Ileum wall volume (mL)")     # Table 1 Ileum
    v_gw_cecum    <- fixed(1.55);  label("Cecum wall volume (mL)")     # Table 1 Cecum
    v_gw_colon    <- fixed(2.5); label("Colon wall volume (mL)")     # Table 1 Colon

    q_co      <- fixed(83.9);    label("Cardiac output = lung blood flow (mL/min)") # Table 1 Lungs blood flow; supplement tvQtotal
    q_heart   <- fixed(4.07);     label("Heart blood flow (mL/min)")    # Table 1 Heart
    q_brain   <- fixed(1.66);     label("Brain blood flow (mL/min)")    # Table 1 Brain
    q_muscle  <- fixed(23.1);     label("Muscle blood flow (mL/min)")   # Table 1 Muscle
    q_adipose <- fixed(5.82);     label("Adipose blood flow (mL/min)")  # Table 1 Adipose
    q_skin    <- fixed(4.82);     label("Skin blood flow (mL/min)")     # Table 1 Skin
    q_kidney  <- fixed(11.71);    label("Kidney blood flow (mL/min)")   # Table 1 Kidneys
    q_spleen  <- fixed(1.66);      label("Spleen blood flow (mL/min)")   # Table 1 Spleen
    q_stomach <- fixed(1.13);   label("Stomach blood flow (mL/min)")  # Table 1 Stomach
    q_liver   <- fixed(12.3); label("Total liver (portal + hepatic-artery) blood flow (mL/min)") # Table 1 Liver; supplement liver outflow (Qliver + Qstomach + Qg1..Qg5 + Qspleen)
    q_other   <- fixed(20.42);     label("Rest-of-body blood flow (mL/min)") # Table 1 Rest of body
    q_gw_duodenum <- fixed(0.36); label("Duodenum wall blood flow (mL/min)") # Table 1 Duodenum
    q_gw_jejunum  <- fixed(4.03); label("Jejunum wall blood flow (mL/min)")  # Table 1 Jejunum
    q_gw_ileum    <- fixed(0.09); label("Ileum wall blood flow (mL/min)")    # Table 1 Ileum
    q_gw_cecum    <- fixed(1.16);  label("Cecum wall blood flow (mL/min)")    # Table 1 Cecum
    q_gw_colon    <- fixed(1.88); label("Colon wall blood flow (mL/min)")    # Table 1 Colon

    # ---- Tissue-to-plasma partition coefficients (Table 1 'Kt:p') ---------
    # Rat values scaled by fu/fu_rat (Methods). The liver value already
    # carries the Eq 11 extraction-ratio correction.
    kp_lung    <- fixed(103.76); label("Lung:plasma partition coefficient (unitless)")    # Table 1 Lungs
    kp_heart   <- fixed(6.78);  label("Heart:plasma partition coefficient (unitless)")   # Table 1 Heart
    kp_brain   <- fixed(0.85);  label("Brain:plasma partition coefficient (unitless)")   # Table 1 Brain
    kp_muscle  <- fixed(5.58);  label("Muscle:plasma partition coefficient (unitless)")  # Table 1 Muscle
    kp_adipose <- fixed(0.27);  label("Adipose:plasma partition coefficient (unitless)") # Table 1 Adipose
    kp_skin    <- fixed(2.28);  label("Skin:plasma partition coefficient (unitless)")    # Table 1 Skin
    kp_kidney  <- fixed(24.63); label("Kidney:plasma partition coefficient (unitless)")  # Table 1 Kidneys
    kp_spleen  <- fixed(22.93); label("Spleen:plasma partition coefficient (unitless)")  # Table 1 Spleen
    kp_stomach <- fixed(128.67); label("Stomach:plasma partition coefficient (unitless)") # Table 1 Stomach
    kp_liver   <- fixed(12.07);  label("Liver:plasma partition coefficient, Eq 11 corrected (unitless)") # Table 1 Liver
    kp_other   <- fixed(0.01);  label("Rest-of-body:plasma partition coefficient (unitless)") # Table 1 Rest of body; footnote a 'assumed as 0.01'
    kp_gw      <- fixed(5.19);  label("Gut-wall:plasma partition coefficient, all five segments (unitless)") # Table 1 Duodenum..Colon; supplement tvKpintestines

    # ---- Drug-specific disposition parameters ----------------------------
    bp     <- fixed(0.91);   label("Blood-to-plasma concentration ratio Rb (unitless)")   # Methods 'Rb was assumed to be 0.91'
    fu     <- fixed(0.32);   label("Fraction unbound in plasma (unitless)")               # Methods 'fu were 0.32 (rat), 0.17 (dog)'
    vmax   <- fixed(1.5);   label("Microsomal Vmax of M-I formation (nmol/min/mg protein)") # Results 'Vmax,M-I values ... of rats, dogs ... were 1.50, 0.16' nmol/min/mg protein
    km     <- fixed(12.94);  label("Microsomal Km of M-I formation (uM)")                 # Results 'Km,M-I values were 12.94 (rat), 90.97 (dog)' uM
    pbsf   <- fixed(409.92);  label("Physiologically based scaling factor, microsomal protein per body (mg)") # Methods 'PBSF values were 409.92 (rat), 16592.7 (dog)' mg protein/body
    cl_other <- fixed(25);  label("Other (extrahepatic + non-M-I) clearance from venous blood (mL/min)") # Methods '25, 150, and 616 mL/min in rats, dogs and humans'; Results 'clearance of the extrahepatic metabolism was 100 mL/min/kg' (Eq 4)
    mw     <- fixed(345.08); label("Vonoprazan molecular weight used for uM <-> ng/mL conversion (g/mol)") # supplement tvKm 4693.14 ng/mL / Results Km 13.60 uM
    ps_stomach <- 0.6; label("Stomach permeability-surface area product PS (mL/min)") # Results 'PS value of vonoprazan in the rat stomach was estimated to be 0.6 mL/min' (fitted to rat stomach data)

    # ---- Oral absorption and gut-lumen transit (Table 1, Eqs 15-19) ------
    peff        <- fixed(0.002); label("Effective intestinal permeability Peff (cm/min)") # Methods 'values of Peff were 0.002 (rat), 0.008 (dog)' cm/min (Eq 19; dog assumed equal to human)
    r_duodenum  <- fixed(0.2);  label("Duodenum radius r1 (cm)")  # Table 1 r1
    r_jejunum   <- fixed(0.2);  label("Jejunum radius r2 (cm)")   # Table 1 r2
    r_ileum     <- fixed(0.2);  label("Ileum radius r3 (cm)")     # Table 1 r3
    ktr_stomach  <- fixed(0.034);  label("Gastric emptying rate constant k0 (1/min)")  # Table 1 k0
    ktr_duodenum <- fixed(0.48);  label("Duodenum transit rate constant k1 (1/min)")  # Table 1 k1
    ktr_jejunum  <- fixed(0.301);  label("Jejunum transit rate constant k2 (1/min)")   # Table 1 k2
    ktr_ileum    <- fixed(0.02);  label("Ileum transit rate constant k3 (1/min)")     # Table 1 k3
    ktr_cecum    <- fixed(0); label("Cecum transit rate constant k4 (1/min)")     # Table 1 k4 'ND' (not detected); no transit out of this segment
    ktr_colon    <- fixed(0); label("Colon transit rate constant k5 (1/min)")     # Table 1 k5 'ND' (not detected); no transit out of this segment

    # ---- PD: H+/K+-ATPase inhibition (Eqs 21-22) --------------------------
    kd       <- fixed(0.00246); label("Dissociation rate constant of vonoprazan from H+/K+-ATPase (1/min)") # Methods 'kd was estimated to be 0.00246 min-1' (t1/2 4.7 h)
    ki       <- fixed(0.035);   label("Unbound concentration at 50% maximal inhibition KI (uM)")           # Methods 'KI value was set to 35 nM'
    imax     <- fixed(4.0);     label("Maximal increase in gastric perfusate pH Imax (pH units)")              # Results rat 'basal pH and I max were set to 2.0 and 4.0'
    ph_basal <- fixed(2.0);     label("Basal (predose) gastric perfusate pH (pH units)")                       # Results rat 'basal pH and I max were set to 2.0 and 4.0'

    # ---- Residual error --------------------------------------------------
    # Deterministic prediction model; the paper reports no residual-error
    # estimate, so it is fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported)") # no residual-error model reported for the animal predictions
  })

  model({
    # ================= Derived quantities ================================
    # Hepatic-artery flow closes the Table 1 liver flow (supplement code
    # 'Qliver' = 300 mL/min for the human).
    q_ha <- q_liver - q_stomach - q_spleen - q_gw_duodenum - q_gw_jejunum -
      q_gw_ileum - q_gw_cecum - q_gw_colon

    # Segmental absorption rate constants, Eq 17 (ka,i = 2 Peff / r_i).
    ka_duodenum <- 2 * peff / r_duodenum
    ka_jejunum  <- 2 * peff / r_jejunum
    ka_ileum    <- 2 * peff / r_ileum

    # Binding rate constant k = kd / KI (Supplementary Eq 5).
    kon <- kd / ki

    # ================= Concentrations (mg/mL; blood for the pools) =======
    c_venous   <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_vp_stomach <- vp_stomach / v_vp_stomach
    c_ev_stomach <- ev_stomach / v_ev_stomach

    # Emergent venous-blood concentrations of each perfusion-limited tissue
    # (C_t / (Kt:p / Rb), Eq 2).
    cv_lung    <- lung    / v_lung    / (kp_lung    / bp)
    cv_heart   <- heart   / v_heart   / (kp_heart   / bp)
    cv_brain   <- brain   / v_brain   / (kp_brain   / bp)
    cv_muscle  <- muscle  / v_muscle  / (kp_muscle  / bp)
    cv_adipose <- adipose / v_adipose / (kp_adipose / bp)
    cv_skin    <- skin    / v_skin    / (kp_skin    / bp)
    cv_kidney  <- kidney  / v_kidney  / (kp_kidney  / bp)
    cv_spleen  <- spleen  / v_spleen  / (kp_spleen  / bp)
    cv_liver   <- liver   / v_liver   / (kp_liver   / bp)
    cv_other   <- other   / v_other   / (kp_other   / bp)
    cv_gw_duodenum <- gw_duodenum / v_gw_duodenum / (kp_gw / bp)
    cv_gw_jejunum  <- gw_jejunum  / v_gw_jejunum  / (kp_gw / bp)
    cv_gw_ileum    <- gw_ileum    / v_gw_ileum    / (kp_gw / bp)
    cv_gw_cecum    <- gw_cecum    / v_gw_cecum    / (kp_gw / bp)
    cv_gw_colon    <- gw_colon    / v_gw_colon    / (kp_gw / bp)

    # ================= Hepatic metabolism (Eq 9) =========================
    # Unbound liver concentration as coded in the supplement
    # (Cliver / (Kpliver / Rb) * fub), in mg/mL. Km and Vmax are converted
    # from uM and nmol/min/mg to mg/mL and mg/min/mg with mw.
    cu_liver <- fu * cv_liver
    km_mass  <- km * mw * 1e-6
    met_liver <- pbsf * vmax * mw * 1e-6 * cu_liver / (km_mass + cu_liver)

    # ================= Gut lumen (Eqs 15-16), amounts in mg ==============
    d/dt(stomach)  <- -ktr_stomach * stomach
    d/dt(duodenum) <- ktr_stomach * stomach - ktr_duodenum * duodenum - ka_duodenum * duodenum
    d/dt(jejunum)  <- ktr_duodenum * duodenum - ktr_jejunum * jejunum - ka_jejunum * jejunum
    d/dt(ileum)    <- ktr_jejunum * jejunum - ktr_ileum * ileum - ka_ileum * ileum
    d/dt(cecum)    <- ktr_ileum * ileum - ktr_cecum * cecum
    d/dt(colon)    <- ktr_cecum * cecum - ktr_colon * colon

    # ================= Gut wall (Eq 20) ==================================
    d/dt(gw_duodenum) <- q_gw_duodenum * c_arterial + ka_duodenum * duodenum - q_gw_duodenum * cv_gw_duodenum
    d/dt(gw_jejunum)  <- q_gw_jejunum  * c_arterial + ka_jejunum  * jejunum  - q_gw_jejunum  * cv_gw_jejunum
    d/dt(gw_ileum)    <- q_gw_ileum    * c_arterial + ka_ileum    * ileum    - q_gw_ileum    * cv_gw_ileum
    d/dt(gw_cecum)    <- q_gw_cecum    * c_arterial - q_gw_cecum * cv_gw_cecum
    d/dt(gw_colon)    <- q_gw_colon    * c_arterial - q_gw_colon * cv_gw_colon

    # ================= Stomach wall, permeability-limited (Eqs 12-13) ====
    d/dt(vp_stomach) <- q_stomach * (c_arterial - c_vp_stomach) -
      ps_stomach * (c_vp_stomach - c_ev_stomach / (kp_stomach / bp))
    d/dt(ev_stomach) <- ps_stomach * (c_vp_stomach - c_ev_stomach / (kp_stomach / bp))

    # ================= Perfusion-limited tissues (Eq 2) ==================
    d/dt(heart)   <- q_heart   * (c_arterial - cv_heart)
    d/dt(brain)   <- q_brain   * (c_arterial - cv_brain)
    d/dt(muscle)  <- q_muscle  * (c_arterial - cv_muscle)
    d/dt(adipose) <- q_adipose * (c_arterial - cv_adipose)
    d/dt(skin)    <- q_skin    * (c_arterial - cv_skin)
    d/dt(kidney)  <- q_kidney  * (c_arterial - cv_kidney)
    d/dt(spleen)  <- q_spleen  * (c_arterial - cv_spleen)
    d/dt(other)   <- q_other   * (c_arterial - cv_other)

    # ================= Liver (Eq 8, as coded in the supplement) ==========
    # Portal inflow follows the supplement listing verbatim: the stomach
    # vascular outflow enters as Qstomach * C1 / (Kp,stomach / Rb), and only
    # the three absorbing gut-wall segments drain to the liver. The cecum
    # and colon wall outflows are not routed anywhere. Both features make
    # these organs act as extra elimination; they reproduce the paper's
    # published oral predictions exactly (see the vignette Errata).
    d/dt(liver) <- q_ha * c_arterial +
      q_stomach * c_vp_stomach / (kp_stomach / bp) +
      q_gw_duodenum * cv_gw_duodenum + q_gw_jejunum * cv_gw_jejunum + q_gw_ileum * cv_gw_ileum +
      q_spleen * cv_spleen -
      q_liver * cv_liver - met_liver

    # ================= Lung and blood pools (Eqs 3, 6, 7) ================
    d/dt(lung) <- q_co * (c_venous - cv_lung)
    d/dt(arterial) <- q_co * (cv_lung - c_arterial)
    d/dt(venous) <- q_heart * cv_heart + q_brain * cv_brain + q_muscle * cv_muscle +
      q_adipose * cv_adipose + q_skin * cv_skin + q_kidney * cv_kidney +
      q_liver * cv_liver + q_other * cv_other -
      q_co * c_venous - cl_other * c_venous

    # ================= PD: H+/K+-ATPase inhibition (Eq 22) ===============
    # Unbound extravascular stomach concentration in uM drives binding.
    cu_stomach_um <- fu * c_ev_stomach * 1e6 / mw
    d/dt(inhibition) <- kon * cu_stomach_um * (imax - inhibition) - kd * inhibition

    # ================= Outputs ===========================================
    Cc <- 1e6 * c_venous / bp            # venous plasma, ng/mL (supplement Cvenous / Rb)
    Cstomach <- 1e6 * c_ev_stomach       # extravascular stomach wall, ng/g (Fig 2e)
    gastric_ph <- ph_basal + inhibition  # gastric perfusate pH (Fig 2f)

    Cc ~ prop(propSd)
  })
}
