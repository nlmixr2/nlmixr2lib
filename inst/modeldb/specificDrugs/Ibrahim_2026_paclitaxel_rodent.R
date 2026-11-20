Ibrahim_2026_paclitaxel_rodent <- function() {
  description <- paste(
    "Preclinical (rat and mouse). Semi-mechanistic PK-PD model of intravenous",
    "Cremophor EL (CrEL)-formulated paclitaxel and serum neurofilament light chain",
    "(sNfL) in rodents (Ibrahim 2026). Plasma: two-compartment paclitaxel model",
    "(paclitaxel not bound to CrEL) with linear plus saturable (Tmax/Kmt)",
    "distribution to the peripheral compartment, coupled to an allometrically",
    "scaled three-compartment human CrEL model with Michaelis-Menten elimination;",
    "total plasma paclitaxel = unbound-to-CrEL paclitaxel x (1 + Bslope x CrEL",
    "concentration). Peripheral-neurotoxicity sites (dorsal root ganglia, peripheral",
    "nerves, brain): unbound paclitaxel crosses each blood-tissue barrier without",
    "depleting plasma (CLin = Kp,uu x CLout), the intracellular unbound",
    "concentration is Kp,uu,cell x the extracellular one, and paclitaxel-tubulin",
    "binding is capacity-limited with a binding capacity that grows with the",
    "cumulative total plasma AUC (Hill function in DRG and brain, exponential in",
    "nerves); the brain extracellular space has a second exchanging compartment.",
    "NfL: brain and peripheral-nerve NfL pools leak into CSF and serum at a rate",
    "proportional to the amount of tubulin-bound paclitaxel; CSF and serum NfL",
    "follow turnover models at baseline steady state. Typical-value model (no",
    "between-animal variability was estimated); log-scale residual error on each",
    "of the six observed outputs. Paclitaxel and CrEL are dosed into separate",
    "compartments (mg and mL).",
    sep = " "
  )
  reference <- paste(
    "Ibrahim EIK, Girdenyte M, Hu Y, Di Cesare Mannelli L, Balayssac D,",
    "Busserolles J, Theil D, Roussignol G, Perrault O, Le Berre N, Chanut F,",
    "Slaoui M, Loryan I, Friberg LE. Evaluation of sNfL as a Biomarker for",
    "Paclitaxel-Induced Peripheral Neurotoxicity Through an Integrated PKPD Model.",
    "Pharm Res. 2026. doi:10.1007/s11095-026-04053-z.",
    sep = " "
  )
  vignette <- "Ibrahim_2026_paclitaxel_rodent"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling: paclitaxel clearances (exponent 0.75) and volumes",
        "(exponent 1) are referenced to a 0.300 kg rat; the CrEL distribution",
        "parameters keep their human 70 kg reference (Henningsson 2005 model),",
        "while the re-estimated CrEL Vmax is referenced to 0.300 kg; the NfL",
        "serum elimination rate is scaled from the human 10.1 1/day with exponent",
        "-0.25 and a 70 kg reference. DRG and peripheral-nerve weights scale",
        "linearly with WT (40 mg and 2781 mg at 0.300 kg), as does rat brain",
        "weight (1.8 g at 0.300 kg). Supplement 1 control stream, $PK.",
        sep = " "
      ),
      source_name = "WT"
    ),
    SPECIES_MOUSE = list(
      description = "Mouse species indicator (1 = mouse, 0 = rat)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = rat",
      notes = paste(
        "Only changes the brain weight: 1.8 g x WT/0.300 for rats and a fixed",
        "0.4 g for mice (control stream 'IF(SRC.EQ.4) / IF(SRC.EQ.6) WT_BRN =",
        "(0.4*1E3)/1E6', sources 4 and 6 being the CD1 and BALB/c mouse",
        "studies). All other species differences are carried by WT. The NfL",
        "sub-model was fitted to rat data only.",
        sep = " "
      ),
      source_name = "SRC (data source 4 or 6 = mouse)"
    )
  )

  # States that are not in the canonical compartment register: the three CrEL
  # vehicle compartments, the extracellular paclitaxel pools of DRG and
  # peripheral nerve, the second brain compartment, and the four NfL pools.
  paper_specific_compartments <- c(
    "crel_central",
    "crel_peripheral1",
    "crel_peripheral2",
    "ecf_drg",
    "ecf_pn",
    "brain_ecf2",
    "nfl_serum",
    "nfl_csf",
    "nfl_brain",
    "nfl_pns"
  )

  compartmentData <- list(
    central = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    crel_central = list(analyte = "Cremophor EL", units = "mL", specimen = "plasma", verified = TRUE),
    crel_peripheral1 = list(analyte = "Cremophor EL", units = "mL", specimen = "tissue", verified = TRUE),
    crel_peripheral2 = list(analyte = "Cremophor EL", units = "mL", specimen = "tissue", verified = TRUE),
    auc_central = list(analyte = "paclitaxel", units = "mg*h/L", specimen = "not applicable", verified = TRUE),
    ecf_drg = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    ecf_pn = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "paclitaxel", units = "mg", specimen = "brain ISF", verified = TRUE),
    brain_ecf2 = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    nfl_serum = list(analyte = "neurofilament light chain", units = "pg", specimen = "serum", verified = TRUE),
    nfl_csf = list(analyte = "neurofilament light chain", units = "pg", specimen = "CSF", verified = TRUE),
    nfl_brain = list(analyte = "neurofilament light chain", units = "pg", specimen = "tissue", verified = TRUE),
    nfl_pns = list(analyte = "neurofilament light chain", units = "pg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Wistar Hanover, Sprague-Dawley) + mouse (CD1, BALB/c)",
    n_subjects = 73,
    n_studies = 7,
    sex_female_pct = NA,
    weight_range = "about 0.3 kg (rat reference weight); mice about 0.025-0.03 kg",
    disease_state = "healthy (paclitaxel-induced peripheral neurotoxicity models)",
    dose_range = paste(
      "Paclitaxel (CrEL formulation) IV: 2.6-30 mg/kg single doses, 5 mg/kg as a",
      "4-h infusion, 5 mg/kg weekly for 2 or 4 weeks, 7.5 mg/kg weekly for 6",
      "weeks, 5/10/15 mg/kg three times in week 1 and four times in week 2, and",
      "30 mg/kg every 2 days x 3 for 2 weeks",
      sep = " "
    ),
    regions = "Sweden, Italy, France, plus digitized literature data (Nakamura 2016, Wozniak 2017)",
    notes = paste(
      "Ibrahim 2026 Table I: six data sources, 57 animals with individual data",
      "(Uppsala University, University of Florence, University Clermont",
      "Auvergne, Sanofi) plus digitized mean data from Nakamura 2016 (Sprague",
      "Dawley rats) and Wozniak 2017 (BALB/c mice). The NfL kinetics were",
      "built on a de novo longitudinal study in 16 male Sprague-Dawley rats",
      "(8-12 weeks, pilot included; 10 mg/kg CrEL-paclitaxel IV bolus, N = 7,",
      "or saline) and externally validated against single-time-point sNfL",
      "from sources 2-4. n_subjects = 57 + 16.",
      sep = " "
    )
  )

  ini({
    # ---- Paclitaxel plasma PK (not bound to CrEL); reference rat 0.300 kg ----
    lcl <- log(1.07323); label("Clearance of paclitaxel not bound to CrEL, CL_PTX,c (L/h)") # Table III CL_PTX,c = 1.07; Supplement 1 $THETA(1) = 1.07323
    lvc <- log(0.558984); label("Central volume of paclitaxel not bound to CrEL, V_PTX,c (L)") # Table III V_PTX,c = 0.559; $THETA(2) = 0.558984
    lq <- log(16.8238); label("Linear intercompartmental clearance of paclitaxel, Q_PTX,p (L/h)") # Table III Q_PTX,p = 16.8; $THETA(3) = 16.8238
    lvp <- log(7.49171); label("Peripheral volume of paclitaxel not bound to CrEL, V_PTX,p (L)") # Table III V_PTX,p = 7.49; $THETA(4) = 7.49171
    ltmax <- log(0.547305); label("Maximum transport capacity to the peripheral compartment, Tmax (mg/h)") # Table III Tmax = 0.547; $THETA(5) = 0.547305
    lkmt <- log(0.100924); label("Concentration of half-maximal saturable transport, Kmt (ug/L)") # Table III Kmt = 0.101 ug/L; $THETA(6) = 0.100924 (divided by 1000 to mg/L in $PK)

    # ---- CrEL plasma PK (human 70 kg reference except Vmax) ----
    lvc_crel <- fixed(log(4.54)); label("CrEL central volume at 70 kg, V_CrEL,c (L)") # Table III V_CrEL,c = 4.54 FIX (Henningsson 2005)
    lvp_crel <- fixed(log(1.32)); label("CrEL first peripheral volume at 70 kg, V_CrEL,p1 (L)") # Table III V_CrEL,p1 = 1.32 FIX
    lvp2_crel <- fixed(log(3.53)); label("CrEL second peripheral volume at 70 kg, V_CrEL,p2 (L)") # Table III V_CrEL,p2 = 3.53 FIX
    lq_crel <- fixed(log(1.17)); label("CrEL first intercompartmental clearance at 70 kg, Q_CrEL,p1 (L/h)") # Table III Q_CrEL,p1 = 1.17 FIX
    lq2_crel <- fixed(log(0.479)); label("CrEL second intercompartmental clearance at 70 kg, Q_CrEL,p2 (L/h)") # Table III Q_CrEL,p2 = 0.479 FIX
    lvmax_crel <- log(0.175047); label("Maximum CrEL elimination capacity at 0.300 kg, Vmax (mL/h)") # Table III Vmax = 0.175; $THETA(8) = 0.175047
    lkm_crel <- fixed(log(2.57)); label("CrEL concentration of half-maximal elimination, Km_CrEL (mL/L)") # Table III Km_CrEL = 2.57 FIX; units mL/L per Supplement 1 $PK
    bslope_crel <- fixed(4.46); label("Slope of paclitaxel binding to CrEL concentration, Bslope,CrEL (L/mL)") # Table III Bslope,CrEL = 4.46 FIX

    # ---- Blood-to-tissue transfer (all sites) ----
    fu <- fixed(0.069); label("Unbound fraction of paclitaxel in plasma without CrEL, fu (fraction)") # Table III fu = 0.069 FIX
    fcrel <- 0.525151; label("Fraction of CrEL-bound paclitaxel crossing the blood-tissue barriers, fCrEL (fraction)") # Table III fCrEL = 0.525; $THETA(12) = 0.525151

    # ---- Dorsal root ganglia (DRG) ----
    kpuu_drg <- fixed(4.25); label("Unbound DRG-to-plasma partition coefficient, Kp,uu,DRG (ratio)") # Table III Kp,uu,DRG = 4.25 FIX
    kpuucell_drg <- fixed(1.45); label("Unbound ICF-to-ECF concentration ratio in DRG, Kp,uu,cell,DRG (ratio)") # Table III Kp,uu,cell,DRG = 1.45 FIX
    lclout_drg <- log(0.57698); label("Efflux clearance across the blood-DRG barrier, CLout,BDB (uL/h)") # Table III CLout,BDB = 0.577 uL/h; $THETA(9) = 0.57698
    lbmax_drg <- log(0.286703); label("Initial tubulin binding capacity in DRG cells, Bmax,DRG,initial (mg/L)") # Table III Bmax,DRG,initial = 0.287; $THETA(10) = 0.286703
    lkd <- log(0.352706); label("Paclitaxel-tubulin dissociation constant in DRG and nerve cells, Kd (ug/L)") # Table III Kd_DRG = Kd_PN = 0.353 ug/L; $THETA(11) = 0.352706
    lbmax_max <- log(14.9417); label("Maximum relative increase in DRG and brain Bmax, MAX (ratio)") # Table III MAX_DRG = MAX_BRN = 14.9; $THETA(20) = 14.9417
    lauc50 <- log(97.9133); label("Cumulative plasma AUC at half-maximal Bmax increase, AUC50 (mg*h/L)") # Table III AUC50,DRG = AUC50,BRN = 97.9; $THETA(22) = 97.9133
    lhill <- log(2.75525); label("Hill coefficient of the Bmax increase in DRG and brain, h (unitless)") # Table III h_DRG = h_BRN = 2.76; $THETA(21) = 2.75525
    wt_drg <- fixed(40); label("DRG weight at 0.300 kg body weight, WT_DRG (mg)") # Table III WT_DRG = 40 mg FIX
    fecf_drg <- fixed(0.2); label("Extracellular fluid volume of DRG, V_ECF,DRG (mL/g tissue)") # Table III V_ECF,DRG = 0.2 FIX
    ficf_drg <- fixed(0.8); label("Intracellular fluid volume of DRG, V_ICF,DRG (mL/g tissue)") # Table III V_ICF,DRG = 0.8 FIX

    # ---- Peripheral nerves (PN, sciatic nerve) ----
    kpuu_pn <- fixed(4.48); label("Unbound nerve-to-plasma partition coefficient, Kp,uu,PN (ratio)") # Table III Kp,uu,PN = 4.48 FIX
    kpuucell_pn <- fixed(0.55); label("Unbound ICF-to-ECF concentration ratio in nerve, Kp,uu,cell,PN (ratio)") # Table III Kp,uu,cell,PN = 0.55 FIX
    lclout_pn <- log(772.494); label("Efflux clearance across the blood-nerve barrier, CLout,BNB (uL/h)") # Table III CLout,BNB = 772.5 uL/h; $THETA(13) = 772.494
    lbmax_pn <- log(0.144341); label("Initial tubulin binding capacity in nerve cells, Bmax,PN,initial (mg/L)") # Table III Bmax,PN,initial = 0.144; $THETA(14) = 0.144341
    slp_bmax_pn <- 0.0193666; label("Exponential coefficient of nerve Bmax on cumulative plasma AUC, SLP_PN (L/(mg*h))") # Table III SLP_PN = 0.0194; $THETA(26) = 0.0193666
    wt_pn <- fixed(2781); label("Peripheral nerve weight at 0.300 kg body weight, WT_PN (mg)") # Table III WT_PN = 2781 mg FIX
    fecf_pn <- fixed(0.5); label("Extracellular fluid volume of nerve, V_ECF,PN (mL/g tissue)") # Table III V_ECF,PN = 0.5 FIX
    ficf_pn <- fixed(0.5); label("Intracellular fluid volume of nerve, V_ICF,PN (mL/g tissue)") # Table III V_ICF,PN = 0.5 FIX

    # ---- Brain ----
    kpuu_brain <- fixed(0.032); label("Unbound brain-to-plasma partition coefficient, Kp,uu,BRN (ratio)") # Table III Kp,uu,BRN = 0.032 FIX
    kpuucell_brain <- fixed(11.2); label("Unbound ICF-to-ECF concentration ratio in brain, Kp,uu,cell,BRN (ratio)") # Table III Kp,uu,cell,BRN = 11.2 FIX
    lclout_brain <- log(23.6918); label("Efflux clearance across the blood-brain barrier, CLout,BBB (uL/h)") # Table III CLout,BBB = 23.7 uL/h; $THETA(16) = 23.6918
    lq_brain <- log(6.46732); label("Intercompartmental clearance between brain ECF compartments, Qend,BRN (uL/h)") # Table III Qend,BRN = 6.47 uL/h; $THETA(17) = 6.46732
    lvp_brain <- log(3.22582); label("Volume of the second brain compartment, Vend,BRN (mL)") # Table III Vend,BRN = 3.23 mL; $THETA(18) = 3.22582
    lbmax_brain <- log(0.025915); label("Initial tubulin binding capacity in brain cells, Bmax,BRN,initial (mg/L)") # Table III Bmax,BRN,initial = 0.0259; $THETA(19) = 0.025915
    kd_brain <- fixed(4.208248); label("Paclitaxel-tubulin dissociation constant in brain cells, Kd_BRN (ug/L)") # Table III Kd_BRN = 4.208 ug/L FIX (Kuh 2000); Supplement 1 KD_BRN = 0.004208248 mg/L
    wt_brain_rat <- fixed(1.8); label("Rat brain weight at 0.300 kg body weight, WT_BRN (g)") # Table III WT_BRN = 1.8 g (rat) FIX
    wt_brain_mouse <- fixed(0.4); label("Mouse brain weight, WT_BRN (g)") # Table III WT_BRN = 0.4 g (mice) FIX
    fecf_brain <- fixed(0.2); label("Extracellular fluid volume of brain, V_ECF,BRN (mL/g tissue)") # Table III V_ECF,BRN = 0.2 FIX
    ficf_brain <- fixed(0.8); label("Intracellular fluid volume of brain, V_ICF,BRN (mL/g tissue)") # Table III V_ICF,BRN = 0.8 FIX

    # ---- NfL kinetics (rat) ----
    lbl_nfl_serum <- log(17.6446); label("Baseline serum NfL concentration, NfL_baseline,serum (pg/mL)") # Table IV NfL_baseline,serum = 17.64; $THETA(27) = 17.6446
    lbl_nfl_csf <- log(222.501); label("Baseline CSF NfL concentration, NfL_baseline,csf (pg/mL)") # Table IV NfL_baseline,csf = 222.5; $THETA(28) = 222.501
    bl_nfl_brain <- fixed(0.43e6); label("Baseline brain NfL concentration, NfL_baseline,brain (pg/mL)") # Table IV NfL_baseline,brain = 0.43 x 10^6 FIX
    bl_nfl_pns <- fixed(538.74e6); label("Baseline DRG plus nerve NfL concentration, NfL_baseline,PN (pg/mL)") # Table IV NfL_baseline,PN = 538.74 x 10^6 FIX
    lslp_nfl <- log(0.029005); label("Slope of NfL leakage per amount of tubulin-bound paclitaxel, slp (10^9 pg/(h*mg))") # Table IV slp_NfL,brain,PTX = 0.029 x 10^9; $THETA(32) = 0.029005; same slope for brain and DRG + nerve (Table IV footnote)
    kout_nfl_serum <- fixed(10.1); label("Serum NfL elimination rate constant at 70 kg, kout,serum (1/day)") # Table IV kout,serum = 10.1 1/day FIX (Ref 27; allometric exponent -0.25 at 70 kg)
    v_serum <- fixed(5); label("Serum volume, V_serum (mL)") # Table IV V_serum = 5 mL FIX
    v_csf <- fixed(0.25); label("CSF volume, V_csf (mL)") # Table IV V_csf = 0.25 mL FIX

    # ---- Residual error: additive on log-transformed data (SD on log scale) ----
    expSd <- 0.356814; label("Residual SD of log total plasma paclitaxel (log scale)") # Table III RUV_plasma = 36%; $THETA(7) = 0.356814
    expSd_Cdrg <- 0.493365; label("Residual SD of log total DRG paclitaxel (log scale)") # Table III RUV_DRG = 49%; $THETA(24) = 0.493365
    expSd_Cpn <- 0.69105; label("Residual SD of log total nerve paclitaxel (log scale)") # Table III RUV_PN = 69%; $THETA(23) = 0.69105
    expSd_Cbrain <- 0.417205; label("Residual SD of log total brain paclitaxel (log scale)") # Table III RUV_BRN = 42%; $THETA(25) = 0.417205
    expSd_Cnfl_serum <- 0.465805; label("Residual SD of log serum NfL (log scale)") # Table IV RUV_serum = 46%; $THETA(34) = 0.465805
    expSd_Cnfl_csf <- 0.371217; label("Residual SD of log CSF NfL (log scale)") # Table IV RUV_csf = 37%; $THETA(35) = 0.371217
  })

  model({
    # ---- Allometric scaling (Supplement 1 $PK) ----
    wt_rat <- WT / 0.300 # paclitaxel parameters referenced to a 0.300 kg rat
    wt_hum <- WT / 70 # CrEL distribution parameters and NfL kout from the human model

    cl <- exp(lcl) * wt_rat^0.75
    vc <- exp(lvc) * wt_rat
    q <- exp(lq) * wt_rat^0.75
    vp <- exp(lvp) * wt_rat
    tmax <- exp(ltmax)
    kmt <- exp(lkmt) / 1000 # ug/L -> mg/L

    vc_crel <- exp(lvc_crel) * wt_hum
    vp_crel <- exp(lvp_crel) * wt_hum
    vp2_crel <- exp(lvp2_crel) * wt_hum
    q_crel <- exp(lq_crel) * wt_hum^0.75
    q2_crel <- exp(lq2_crel) * wt_hum^0.75
    vmax_crel <- exp(lvmax_crel)
    km_crel <- exp(lkm_crel)

    # Tissue weights (L, density 1 g/mL) and fluid volumes (L)
    wtl_drg <- wt_drg * wt_rat / 1e6
    wtl_pn <- wt_pn * wt_rat / 1e6
    wtl_brain <- (1 - SPECIES_MOUSE) * wt_brain_rat * wt_rat / 1e3 + SPECIES_MOUSE * wt_brain_mouse / 1e3
    vecf_drg <- fecf_drg * wtl_drg
    vicf_drg <- ficf_drg * wtl_drg
    vecf_pn <- fecf_pn * wtl_pn
    vicf_pn <- ficf_pn * wtl_pn
    vecf_brain <- fecf_brain * wtl_brain
    vicf_brain <- ficf_brain * wtl_brain

    # Blood-tissue barrier clearances (uL/h -> L/h); CLin = Kp,uu x CLout
    clout_drg <- exp(lclout_drg) / 1e6 * wt_rat^0.75
    clout_pn <- exp(lclout_pn) / 1e6 * wt_rat^0.75
    clout_brain <- exp(lclout_brain) / 1e6 * wt_rat^0.75
    clin_drg <- kpuu_drg * clout_drg
    clin_pn <- kpuu_pn * clout_pn
    clin_brain <- kpuu_brain * clout_brain
    q_brain <- exp(lq_brain) / 1e6 * wt_rat^0.75
    vp_brain <- exp(lvp_brain) / 1e3 * wt_rat

    kd <- exp(lkd) / 1000 # ug/L -> mg/L
    kd_br <- kd_brain / 1000 # ug/L -> mg/L
    bmax_drg <- exp(lbmax_drg)
    bmax_pn <- exp(lbmax_pn)
    bmax_brain <- exp(lbmax_brain)
    bmax_max <- exp(lbmax_max)
    auc50 <- exp(lauc50)
    hill <- exp(lhill)

    # ---- NfL baseline steady state (Eqs. 13-14, written in amounts) ----
    bl_nfl_serum <- exp(lbl_nfl_serum)
    bl_nfl_csf <- exp(lbl_nfl_csf)
    amt_bl_serum <- bl_nfl_serum * v_serum # pg
    amt_bl_csf <- bl_nfl_csf * v_csf # pg
    amt_bl_brain <- bl_nfl_brain * wtl_brain * 1e3 # pg
    amt_bl_pns <- bl_nfl_pns * (wtl_pn + wtl_drg) * 1e3 # pg
    kout_serum <- kout_nfl_serum / 24 * wt_hum^(-0.25) # 1/day -> 1/h
    kout_csf <- kout_serum * amt_bl_serum / amt_bl_csf
    kin_csf <- amt_bl_csf * kout_csf # pg/h
    slp_nfl <- exp(lslp_nfl) * 1e9 # pg/(h*mg)

    # ---- Plasma paclitaxel and CrEL ----
    cu_ptx <- central / vc # mg/L, paclitaxel not bound to CrEL
    c_crel <- crel_central / vc_crel # mL/L, CrEL not bound to paclitaxel
    q_nl <- tmax / (kmt + cu_ptx) * wt_rat^0.75 # L/h, saturable distribution
    cl_crel <- vmax_crel / (km_crel + c_crel) * wt_rat^0.75 # L/h, saturable CrEL elimination
    ctot_ptx <- cu_ptx * (1 + bslope_crel * c_crel) # mg/L, total plasma paclitaxel (Eq. 5)

    d/dt(central) <- -(cl / vc) * central - (q_nl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q_nl / vc) * central + (q / vc) * central - (q / vp) * peripheral1
    d/dt(crel_central) <- -(cl_crel / vc_crel) * crel_central - (q_crel / vc_crel) * crel_central - (q2_crel / vc_crel) * crel_central + (q_crel / vp_crel) * crel_peripheral1 + (q2_crel / vp2_crel) * crel_peripheral2
    d/dt(crel_peripheral1) <- (q_crel / vc_crel) * crel_central - (q_crel / vp_crel) * crel_peripheral1
    d/dt(crel_peripheral2) <- (q2_crel / vc_crel) * crel_central - (q2_crel / vp2_crel) * crel_peripheral2
    d/dt(auc_central) <- ctot_ptx # mg*h/L, cumulative total plasma AUC

    # Paclitaxel available to cross each barrier: unbound fraction of the
    # CrEL-free amount plus fCrEL of the CrEL-bound amount (Eq. 6). Transfer
    # does not deplete plasma (negligible mass transfer).
    a_avail <- central * fu + bslope_crel * c_crel * cu_ptx * vc * fcrel

    d/dt(ecf_drg) <- clin_drg / vc * a_avail - clout_drg / vecf_drg * ecf_drg
    d/dt(ecf_pn) <- clin_pn / vc * a_avail - clout_pn / vecf_pn * ecf_pn
    d/dt(brain_ecf) <- clin_brain / vc * a_avail - clout_brain / vecf_brain * brain_ecf - q_brain / vecf_brain * brain_ecf + q_brain / vp_brain * brain_ecf2
    d/dt(brain_ecf2) <- q_brain / vecf_brain * brain_ecf - q_brain / vp_brain * brain_ecf2

    # ---- Intracellular paclitaxel and tubulin binding (Eqs. 2-3, 7-8) ----
    auc_cum <- auc_central
    cu_icf_drg <- ecf_drg / vecf_drg * kpuucell_drg
    cu_icf_pn <- ecf_pn / vecf_pn * kpuucell_pn
    cu_icf_brain <- brain_ecf / vecf_brain * kpuucell_brain
    bmax_drg_t <- bmax_drg * (1 + bmax_max * auc_cum^hill / (auc50^hill + auc_cum^hill))
    bmax_brain_t <- bmax_brain * (1 + bmax_max * auc_cum^hill / (auc50^hill + auc_cum^hill))
    bmax_pn_t <- bmax_pn * exp(slp_bmax_pn * auc_cum)
    cb_icf_drg <- bmax_drg_t * cu_icf_drg / (kd + cu_icf_drg)
    cb_icf_pn <- bmax_pn_t * cu_icf_pn / (kd + cu_icf_pn)
    cb_icf_brain <- bmax_brain_t * cu_icf_brain / (kd_br + cu_icf_brain)
    bound_drg <- cb_icf_drg * vicf_drg # mg tubulin-bound paclitaxel
    bound_pn <- cb_icf_pn * vicf_pn
    bound_brain <- cb_icf_brain * vicf_brain

    # ---- NfL kinetics (Eqs. 9-12) ----
    nfl_leak_brain <- slp_nfl * bound_brain # pg/h
    nfl_leak_pns <- slp_nfl * (bound_drg + bound_pn) # pg/h
    d/dt(nfl_serum) <- kout_csf * nfl_csf + nfl_leak_pns - kout_serum * nfl_serum
    d/dt(nfl_csf) <- nfl_leak_brain + kin_csf - kout_csf * nfl_csf
    d/dt(nfl_brain) <- -nfl_leak_brain
    d/dt(nfl_pns) <- -nfl_leak_pns
    nfl_serum(0) <- amt_bl_serum
    nfl_csf(0) <- amt_bl_csf
    nfl_brain(0) <- amt_bl_brain
    nfl_pns(0) <- amt_bl_pns

    # ---- Outputs ----
    tot_drg <- ecf_drg + (cu_icf_drg + cb_icf_drg) * vicf_drg # mg
    tot_pn <- ecf_pn + (cu_icf_pn + cb_icf_pn) * vicf_pn # mg
    tot_brain <- brain_ecf + brain_ecf2 + (cu_icf_brain + cb_icf_brain) * vicf_brain # mg

    Cc <- ctot_ptx * 1000 # ng/mL total plasma paclitaxel
    Cdrg <- tot_drg / wtl_drg * 1000 # ng/mL (ng/g) total DRG paclitaxel
    Cpn <- tot_pn / wtl_pn * 1000 # ng/mL (ng/g) total nerve paclitaxel
    Cbrain <- tot_brain / wtl_brain * 1000 # ng/mL (ng/g) total brain paclitaxel
    Cnfl_serum <- nfl_serum / v_serum # pg/mL
    Cnfl_csf <- nfl_csf / v_csf # pg/mL
    Cnfl_brain <- nfl_brain / (wtl_brain * 1e3) # pg/mL
    Cnfl_pns <- nfl_pns / ((wtl_pn + wtl_drg) * 1e3) # pg/mL

    Cc ~ lnorm(expSd)
    Cdrg ~ lnorm(expSd_Cdrg)
    Cpn ~ lnorm(expSd_Cpn)
    Cbrain ~ lnorm(expSd_Cbrain)
    Cnfl_serum ~ lnorm(expSd_Cnfl_serum)
    Cnfl_csf ~ lnorm(expSd_Cnfl_csf)
  })
}
