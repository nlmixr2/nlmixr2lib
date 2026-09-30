Kong_2020_vonoprazan_human_pbpk <- function() {
  description <- paste(
    "PBPK-PD (whole-body, perfusion-limited, WinNonlin 8.1).",
    "Vonoprazan disposition and gastric-acid antisecretory effect in a",
    "typical 70-kg adult after oral dosing, extrapolated by Kong et al.",
    "(2020) from in vitro metabolism (hepatic microsomes), permeability",
    "(Caco-2) and rat tissue-distribution data. Thirteen perfusion-limited",
    "tissues plus venous and arterial blood pools; the stomach wall (the",
    "target organ) is permeability-limited with a vascular and an",
    "extravascular space. Oral doses enter the stomach lumen and transit a",
    "five-segment gut lumen (duodenum, jejunum, ileum absorb; cecum and",
    "colon do not) into per-segment gut-wall compartments. Hepatic",
    "Michaelis-Menten metabolism of unbound drug is scaled from microsomes",
    "by a physiologically based scaling factor, and an allometrically",
    "scaled linear 'other' clearance acts on venous blood. The PD layer",
    "drives an H+/K+-ATPase inhibition state from the unbound",
    "extravascular stomach concentration (dI/dt = k * fu * C2 * (Imax - I)",
    "- kd * I) and reports intragastric pH as basal pH + I. Deterministic",
    "typical-value model with no IIV. Companion rat and dog models:",
    "modellib('Kong_2020_vonoprazan_rat_pbpk') and",
    "modellib('Kong_2020_vonoprazan_dog_pbpk')."
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
    inhibition = list(
      analyte = "H+/K+-ATPase inhibition",
      units = "pH units",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 6L,
    age_range = "adults (healthy volunteers in the cited clinical studies)",
    weight_range = "70 kg reference adult used for the physiology (Table 1)",
    disease_state = "healthy volunteers",
    dose_range = "Oral vonoprazan 5-40 mg single dose and 10-40 mg once daily for 7 days",
    regions = "Japan and the United Kingdom (literature data)",
    notes = paste(
      "Bottom-up prediction model, not a population fit. Human",
      "pharmacokinetic and intragastric-pH observations were taken from",
      "published healthy-volunteer studies (paper refs 5, 14, 15, 28, 29,",
      "32) and used only for comparison; no human data informed the",
      "structural parameters. The paper's visual predictive check",
      "estimated variances on CLliver, fu and k from the literature data",
      "but did not report them, so the model is typical-value only."
    )
  )

  ini({
    # ---- Physiology: 70-kg human (Table 1, 'Human (70 kg)' columns) -------
    # Volumes in mL, blood flows in mL/min. The supplement WinNonlin listing
    # carries the same values (fixef tvV*, tvQ*).
    v_lung    <- fixed(1170);    label("Lung volume (mL)")          # Table 1 Lungs
    v_heart   <- fixed(310);     label("Heart volume (mL)")         # Table 1 Heart
    v_brain   <- fixed(1450);    label("Brain volume (mL)")         # Table 1 Brain
    v_muscle  <- fixed(35000);   label("Muscle volume (mL)")        # Table 1 Muscle
    v_adipose <- fixed(10000);   label("Adipose volume (mL)")       # Table 1 Adipose
    v_skin    <- fixed(7800);    label("Skin volume (mL)")          # Table 1 Skin
    v_kidney  <- fixed(280);     label("Kidney volume (mL)")        # Table 1 Kidneys
    v_spleen  <- fixed(190);     label("Spleen volume (mL)")        # Table 1 Spleen
    v_liver   <- fixed(1690);    label("Liver volume (mL)")         # Table 1 Liver
    v_venous  <- fixed(3470);    label("Venous blood volume (mL)")  # Table 1 Vein
    v_arterial <- fixed(1730);   label("Arterial blood volume (mL)") # Table 1 Artery
    v_other   <- fixed(5100);    label("Rest-of-body volume (mL)")  # Table 1 Rest of body

    # Stomach wall: vascular (V1) and extravascular (V2) spaces, Eqs 12-13.
    v_vp_stomach <- fixed(25.6);  label("Stomach vascular volume V1 (mL)")      # Methods 'PBPK model development': human V1 25.6 mL
    v_ev_stomach <- fixed(134.4); label("Stomach extravascular volume V2 (mL)") # Methods: human V2 134.4 mL (V1 + V2 = Table 1 Stomach 160 mL)

    # Gut-wall (enterocyte) segments, Eq 20.
    v_gw_duodenum <- fixed(70);   label("Duodenum wall volume (mL)")  # Table 1 Duodenum
    v_gw_jejunum  <- fixed(209);  label("Jejunum wall volume (mL)")   # Table 1 Jejunum
    v_gw_ileum    <- fixed(139);  label("Ileum wall volume (mL)")     # Table 1 Ileum
    v_gw_cecum    <- fixed(116);  label("Cecum wall volume (mL)")     # Table 1 Cecum
    v_gw_colon    <- fixed(1116); label("Colon wall volume (mL)")     # Table 1 Colon

    q_co      <- fixed(5600);    label("Cardiac output = lung blood flow (mL/min)") # Table 1 Lungs blood flow; supplement tvQtotal
    q_heart   <- fixed(240);     label("Heart blood flow (mL/min)")    # Table 1 Heart
    q_brain   <- fixed(700);     label("Brain blood flow (mL/min)")    # Table 1 Brain
    q_muscle  <- fixed(750);     label("Muscle blood flow (mL/min)")   # Table 1 Muscle
    q_adipose <- fixed(260);     label("Adipose blood flow (mL/min)")  # Table 1 Adipose
    q_skin    <- fixed(300);     label("Skin blood flow (mL/min)")     # Table 1 Skin
    q_kidney  <- fixed(1240);    label("Kidney blood flow (mL/min)")   # Table 1 Kidneys
    q_spleen  <- fixed(80);      label("Spleen blood flow (mL/min)")   # Table 1 Spleen
    q_stomach <- fixed(38.33);   label("Stomach blood flow (mL/min)")  # Table 1 Stomach
    q_liver   <- fixed(1518.33); label("Total liver (portal + hepatic-artery) blood flow (mL/min)") # Table 1 Liver; supplement liver outflow (Qliver + Qstomach + Qg1..Qg5 + Qspleen)
    q_other   <- fixed(592);     label("Rest-of-body blood flow (mL/min)") # Table 1 Rest of body
    q_gw_duodenum <- fixed(118); label("Duodenum wall blood flow (mL/min)") # Table 1 Duodenum
    q_gw_jejunum  <- fixed(413); label("Jejunum wall blood flow (mL/min)")  # Table 1 Jejunum
    q_gw_ileum    <- fixed(244); label("Ileum wall blood flow (mL/min)")    # Table 1 Ileum
    q_gw_cecum    <- fixed(44);  label("Cecum wall blood flow (mL/min)")    # Table 1 Cecum
    q_gw_colon    <- fixed(281); label("Colon wall blood flow (mL/min)")    # Table 1 Colon

    # ---- Tissue-to-plasma partition coefficients (Table 1 'Kt:p') ---------
    # Rat values scaled by fu/fu_rat (Methods). The liver value already
    # carries the Eq 11 extraction-ratio correction.
    kp_lung    <- fixed(48.64); label("Lung:plasma partition coefficient (unitless)")    # Table 1 Lungs
    kp_heart   <- fixed(3.18);  label("Heart:plasma partition coefficient (unitless)")   # Table 1 Heart
    kp_brain   <- fixed(0.40);  label("Brain:plasma partition coefficient (unitless)")   # Table 1 Brain
    kp_muscle  <- fixed(2.62);  label("Muscle:plasma partition coefficient (unitless)")  # Table 1 Muscle
    kp_adipose <- fixed(0.13);  label("Adipose:plasma partition coefficient (unitless)") # Table 1 Adipose
    kp_skin    <- fixed(1.07);  label("Skin:plasma partition coefficient (unitless)")    # Table 1 Skin
    kp_kidney  <- fixed(11.55); label("Kidney:plasma partition coefficient (unitless)")  # Table 1 Kidneys
    kp_spleen  <- fixed(10.75); label("Spleen:plasma partition coefficient (unitless)")  # Table 1 Spleen
    kp_stomach <- fixed(60.31); label("Stomach:plasma partition coefficient (unitless)") # Table 1 Stomach
    kp_liver   <- fixed(2.89);  label("Liver:plasma partition coefficient, Eq 11 corrected (unitless)") # Table 1 Liver
    kp_other   <- fixed(0.01);  label("Rest-of-body:plasma partition coefficient (unitless)") # Table 1 Rest of body; footnote a 'assumed as 0.01'
    kp_gw      <- fixed(2.43);  label("Gut-wall:plasma partition coefficient, all five segments (unitless)") # Table 1 Duodenum..Colon; supplement tvKpintestines

    # ---- Drug-specific disposition parameters ----------------------------
    bp     <- fixed(0.91);   label("Blood-to-plasma concentration ratio Rb (unitless)")   # Methods 'Rb was assumed to be 0.91'
    fu     <- fixed(0.15);   label("Fraction unbound in plasma (unitless)")               # Methods 'fu were ... 0.15 (human)'
    vmax   <- fixed(0.24);   label("Microsomal Vmax of M-I formation (nmol/min/mg protein)") # Results 'Vmax,M-I ... 0.24 nmol/min/mg protein' (human)
    km     <- fixed(13.60);  label("Microsomal Km of M-I formation (uM)")                 # Results 'Km,M-I ... 13.60 (human) uM'
    pbsf   <- fixed(82472);  label("Physiologically based scaling factor, microsomal protein per body (mg)") # Methods 'PBSF ... 82,472 (human) mg protein/body'
    cl_other <- fixed(616);  label("Other (extrahepatic + non-M-I) clearance from venous blood (mL/min)") # Results 'CLother, was estimated to be 616 mL/min'; Eq 5 from dog 150 mL/min
    mw     <- fixed(345.08); label("Vonoprazan molecular weight used for uM <-> ng/mL conversion (g/mol)") # supplement tvKm 4693.14 ng/mL / Results Km 13.60 uM
    ps_stomach <- fixed(26.17); label("Stomach permeability-surface area product PS (mL/min)") # Results 'PS values in the stomach of ... humans ... 26.17 mL/min' (Eq 14)

    # ---- Oral absorption and gut-lumen transit (Table 1, Eqs 15-19) ------
    peff        <- fixed(0.008); label("Effective intestinal permeability Peff (cm/min)") # Methods 'values of Peff were ... 0.008 (human) cm/min' (Eq 18)
    r_duodenum  <- fixed(2.00);  label("Duodenum radius r1 (cm)")  # Table 1 r1
    r_jejunum   <- fixed(1.63);  label("Jejunum radius r2 (cm)")   # Table 1 r2
    r_ileum     <- fixed(1.45);  label("Ileum radius r3 (cm)")     # Table 1 r3
    ktr_stomach  <- fixed(0.08);  label("Gastric emptying rate constant k0 (1/min)")  # Table 1 k0
    ktr_duodenum <- fixed(0.07);  label("Duodenum transit rate constant k1 (1/min)")  # Table 1 k1
    ktr_jejunum  <- fixed(0.03);  label("Jejunum transit rate constant k2 (1/min)")   # Table 1 k2
    ktr_ileum    <- fixed(0.04);  label("Ileum transit rate constant k3 (1/min)")     # Table 1 k3
    ktr_cecum    <- fixed(0.003); label("Cecum transit rate constant k4 (1/min)")     # Table 1 k4
    ktr_colon    <- fixed(0.001); label("Colon transit rate constant k5 (1/min)")     # Table 1 k5

    # ---- PD: H+/K+-ATPase inhibition (Eqs 21-22) --------------------------
    kd       <- fixed(0.00246); label("Dissociation rate constant of vonoprazan from H+/K+-ATPase (1/min)") # Methods 'kd was estimated to be 0.00246 min-1' (t1/2 4.7 h)
    ki       <- fixed(0.035);   label("Unbound concentration at 50% maximal inhibition KI (uM)")           # Methods 'KI value was set to 35 nM'
    imax     <- fixed(5.0);     label("Maximal increase in intragastric pH Imax (pH units)")              # Results 'I max was set to 5.0' (human)
    ph_basal <- fixed(2.0);     label("Basal (predose) intragastric pH (pH units)")                       # Results 'basal pH (predose) was assumed to be 2.0'

    # ---- Residual error --------------------------------------------------
    # Deterministic prediction model; the paper reports no residual-error
    # estimate, so it is fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported)") # not reported (VPC variances unpublished)
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
    Cstomach <- 1e6 * c_ev_stomach       # extravascular stomach wall, ng/g (Fig 5d)
    gastric_ph <- ph_basal + inhibition  # intragastric pH (Fig 5e)

    Cc ~ prop(propSd)
  })
}
