Kutumova_2024_albuminNanoparticles_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse, BALB/c, male, 6-8 weeks, 0.024 kg).",
    "PBPK (whole-body, BioUML) model for the biodistribution of intravenous",
    "albumin nanoparticles (ANP, ~120 nm, zeta ~ -30 mV) during induction of and",
    "recovery from lipopolysaccharide-induced acute lung injury.",
    "Seven compartments: venous and arterial plasma plus four mononuclear",
    "phagocyte system organs (lungs, spleen, liver, kidneys) and a rest-of-body",
    "compartment. Each organ carries a membrane-limited triple of capillary",
    "blood, tissue interstitium, and an internalised (phagocytosed) pool, with",
    "time-dependent Hill-function endocytic uptake and first-order exocytic",
    "release; the liver additionally phagocytoses directly from capillary blood",
    "(Kupffer cells) and excretes into bile, and the kidneys excrete into urine.",
    "The permeability (PAC) and distribution (P) coefficients are cohort-specific:",
    "STUDY_LPS30M / STUDY_LPS6H / STUDY_LPS24H select the ANP-after-LPS arm, with",
    "all three zero giving the LPS-naive control. The observed quantity is organ",
    "total radiant efficiency per luminous area (TRE), obtained from the",
    "simulated tissue concentration through a single fitted scale factor k.",
    "Physiologic blood flows, organ volumes, and capillary blood fractions are",
    "inherited from the Cheng 2020 nanoparticle PBPK framework, which the",
    "authors cite as their source for these values.",
    sep = " "
  )
  reference <- paste(
    "Kutumova EO, Akberdin IR, Egorova VS, Kolesova EP, Parodi A,",
    "Pokrovsky VS, Zamyatnin AA Jr, Kolpakov FA (2024).",
    "Physiologically based pharmacokinetic model for predicting the",
    "biodistribution of albumin nanoparticles after induction and recovery",
    "from acute lung injury. Heliyon 10(10):e30962.",
    "doi:10.1016/j.heliyon.2024.e30962.",
    "Physiologic parameters (fractional blood flows, fractional organ volumes,",
    "capillary blood volume fractions) and the IV infusion duration are",
    "inherited from the framework model of Cheng Y-H, He C, Riviere JE,",
    "Monteiro-Riviere NA, Lin Z (2020) ACS Nano 14:3075-3095,",
    "doi:10.1021/acsnano.9b08142, which Kutumova 2024 cites (reference 69) as",
    "the source of those values.",
    sep = " "
  )
  vignette <- "Kutumova_2024_albuminNanoparticles"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # `bile` is the cumulative amount excreted into bile; it is a bookkeeping
  # accumulator rather than a distribution compartment (same role as the
  # canonical `urine` state). Precedent: Mi_2023_cefquinome_pbpk.R.
  paper_specific_compartments <- c("bile")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales cardiac output allometrically (Q = QC * QCC * BW^0.75), all",
        "organ and plasma volumes linearly (V = BW * VC), and the biliary and",
        "urinary excretion clearances allometrically. Kutumova 2024 section 3.3",
        "uses BW = 0.024 kg and states in section 3.5 (Limitations) that",
        "LPS-induced weight loss was deliberately NOT modelled, so WT is",
        "time-invariant in this model."
      ),
      source_name        = "BW"
    ),
    STUDY_LPS30M = list(
      description        = "1 = ANP administered 30 min after intraperitoneal LPS (exp 2); 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LPS-naive control cohort, exp 1)",
      notes              = paste(
        "Selects the LPS-30-min column of Kutumova 2024 Table 1 for every",
        "organ's permeability (PAC) and distribution (P) coefficient."
      ),
      source_name        = "LPS 30 min"
    ),
    STUDY_LPS6H = list(
      description        = "1 = ANP administered 6 h after intraperitoneal LPS (exp 3); 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LPS-naive control cohort, exp 1)",
      notes              = paste(
        "Selects the LPS-6-h column of Kutumova 2024 Table 1. This is the",
        "cohort with peak pulmonary ANP accumulation."
      ),
      source_name        = "LPS 6 h"
    ),
    STUDY_LPS24H = list(
      description        = "1 = ANP administered 24 h after intraperitoneal LPS (exp 4); 0 otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LPS-naive control cohort, exp 1)",
      notes              = paste(
        "Selects the LPS-24-h column of Kutumova 2024 Table 1. Exactly one of",
        "STUDY_LPS30M / STUDY_LPS6H / STUDY_LPS24H may be 1; all three zero",
        "gives the control cohort."
      ),
      source_name        = "LPS 24 h"
    )
  )

  compartmentData <- list(
    venous     = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    arterial   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    vp_lung    = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    is_lung    = list(analyte = "albumin nanoparticles", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "endosome", verified = TRUE),
    vp_spleen  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    is_spleen  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "tissue", verified = TRUE),
    int_spleen = list(analyte = "albumin nanoparticles", units = "mg", specimen = "endosome", verified = TRUE),
    vp_liver   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    is_liver   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "tissue", verified = TRUE),
    int_liver  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "endosome", verified = TRUE),
    vp_kidney  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    is_kidney  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "tissue", verified = TRUE),
    int_kidney = list(analyte = "albumin nanoparticles", units = "mg", specimen = "endosome", verified = TRUE),
    vp_other   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "plasma", verified = TRUE),
    is_other   = list(analyte = "albumin nanoparticles", units = "mg", specimen = "tissue", verified = TRUE),
    int_other  = list(analyte = "albumin nanoparticles", units = "mg", specimen = "endosome", verified = TRUE),
    urine      = list(analyte = "albumin nanoparticles", units = "mg", specimen = "urine", verified = TRUE),
    bile       = list(analyte = "albumin nanoparticles", units = "mg", specimen = "bile", verified = TRUE)
  )

  population <- list(
    species        = "mouse (BALB/c, male, 6-8 weeks old)",
    n_subjects     = 60,
    n_studies      = 1,
    age_range      = "6-8 weeks",
    weight_median  = "0.024 kg",
    sex_female_pct = 0,
    disease_state  = paste(
      "Non-lethal lipopolysaccharide-induced acute lung injury (LPS 6 mg/kg",
      "intraperitoneally, 200 uL per mouse) and LPS-naive controls."
    ),
    dose_range     = paste(
      "Single intravenous ANP dose, 0.5 mg/mouse (about 20.8 mg/kg) in 100 uL",
      "per the experimental protocol (section 2.4); the model-calibration",
      "section 3.3 states PDOSEiv = 208 mg/kg (5 mg / 0.024 kg). See the",
      "vignette Assumptions and deviations section -- this 10-fold discrepancy",
      "is internal to the paper and the parameter set (in particular the fitted",
      "scale factor k) is coherent only with the value the authors simulated."
    ),
    regions        = "Russia (N.N. Blokhin National Medical Research Center of Oncology, Moscow)",
    notes          = paste(
      "n = 3 mice per biodistribution time point per ANP-administration time",
      "point; 60 mice total (section 2.4). ANP given intravenously at 0.5, 6,",
      "or 24 h after LPS, or without LPS (control). Organs (lungs, liver,",
      "spleen, kidneys) were harvested at 10, 180, 360, 1440, and 2880 min",
      "after ANP administration and read by ex-vivo IVIS as total radiant",
      "efficiency normalised to organ area (section 2.4, Fig. 3). The model was",
      "calibrated simultaneously against these four datasets (section 3.3)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Physiologic parameters -- inherited from the Cheng 2020 framework
    # ------------------------------------------------------------------
    # Kutumova 2024 section 3.3: "The remaining parameters (fractional blood
    # flow rates, compartment volumes, and compartment capillary blood volumes)
    # were taken from the model by Cheng et al. [69]." All values below are read
    # from the Cheng 2020 model source code (supporting information
    # nn9b08142_si_001.pdf, "{Physiological parameters}" block). None was
    # estimated here, hence fixed().
    #
    # Cross-check on the one physiologic value Kutumova states directly:
    # Cheng's QCC = 16.5 L/h/kg^0.75; 16.5 / 60 = 0.275 L/min/kg^0.75, which is
    # exactly the "QCC = 0.275 L/min" quoted in Kutumova Appendix A1. The two
    # sources therefore agree and the inheritance is confirmed, not assumed.
    qcc <- fixed(0.275)   ; label("Cardiac output coefficient (L/min/kg^0.75)")            # Kutumova 2024 Appendix A1 (prints "L/min"; the BW^0.75 factor in Q = QC*QCC*BW^0.75 makes the effective unit L/min/kg^0.75). Cheng 2020 SI code QCC = 16.5 L/h/kg^0.75 = 0.275 L/min/kg^0.75

    # Fractional blood flows (fraction of cardiac output, unitless).
    # The lungs sit in series and receive the entire cardiac output; the liver,
    # spleen, kidneys, and rest of body sit in parallel off the arterial side.
    # Kutumova drops Cheng's brain, muscle, and tumour compartments and its
    # portal circulation, so QRC = 1 - (QLC + QSC + QKC) per Appendix A10.
    qc_lung   <- fixed(1)      ; label("Fraction of cardiac output to lungs (unitless)")    # Cheng 2020 SI code QLuC = 1 (Brown 1997)
    qc_liver  <- fixed(0.02)   ; label("Fraction of cardiac output to liver (unitless)")    # Cheng 2020 SI code QLC = 0.02 (Brown 1997 Table 23)
    qc_spleen <- fixed(0.011)  ; label("Fraction of cardiac output to spleen (unitless)")   # Cheng 2020 SI code QSC = 0.011 (Lin 2008; Davies & Morris 1993)
    qc_kidney <- fixed(0.091)  ; label("Fraction of cardiac output to kidneys (unitless)")  # Cheng 2020 SI code QKC = 0.091 (Brown 1997 Table 23)

    # Fractional organ volumes (fraction of body weight, unitless).
    # VRC = 1 - (VLC + VSC + VKC + VLuC + VPlC) per Kutumova Appendix A10.
    vc_lung   <- fixed(0.007)  ; label("Lung volume as a fraction of body weight (unitless)")    # Cheng 2020 SI code VLuC = 0.007 (Brown 1997 Table 21)
    vc_liver  <- fixed(0.055)  ; label("Liver volume as a fraction of body weight (unitless)")   # Cheng 2020 SI code VLC = 0.055 (Brown 1997 Table 21)
    vc_spleen <- fixed(0.005)  ; label("Spleen volume as a fraction of body weight (unitless)")  # Cheng 2020 SI code VSC = 0.005 (Lin 2008; Davies & Morris 1993)
    vc_kidney <- fixed(0.017)  ; label("Kidney volume as a fraction of body weight (unitless)")  # Cheng 2020 SI code VKC = 0.017 (Brown 1997 Table 21)
    vc_plasma <- fixed(0.0355) ; label("Total plasma volume as a fraction of body weight (unitless)")  # Cheng 2020 SI code VPlasmaC = 0.0355 (Davies & Morris 1993, haematocrit 0.41); Kutumova calls this VPlC in Appendix A8

    # Capillary (vascular) blood volume fractions (fraction of organ volume).
    bv_lung   <- fixed(0.5)    ; label("Capillary blood fraction of lung volume (unitless)")    # Cheng 2020 SI code BVLu = 0.5 (Brown 1997 Table 30)
    bv_liver  <- fixed(0.31)   ; label("Capillary blood fraction of liver volume (unitless)")   # Cheng 2020 SI code BVL = 0.31 (Brown 1997 Table 30)
    bv_spleen <- fixed(0.17)   ; label("Capillary blood fraction of spleen volume (unitless)")  # Cheng 2020 SI code BVS = 0.17 (Brown 1997 Table 30)
    bv_kidney <- fixed(0.24)   ; label("Capillary blood fraction of kidney volume (unitless)")  # Cheng 2020 SI code BVK = 0.24 (Brown 1997 Table 30)
    bv_other  <- fixed(0.04)   ; label("Capillary blood fraction of rest-of-body volume (unitless)")  # Cheng 2020 SI code BVR = 0.04 (Brown 1997 Table 30, assumed equal to muscle)

    # ------------------------------------------------------------------
    # Distribution coefficients P -- Kutumova 2024 Table 1 (fitted, cohort-specific)
    # ------------------------------------------------------------------
    # P is the tissue-to-plasma distribution coefficient in r2out
    # (Appendix A1: v2out = PAC * Q * At / (Vt * P)). Higher P means slower
    # return of nanoparticles from interstitium to blood, i.e. more tissue
    # accumulation. Fitted separately for each of the four experiments.
    kp_lung_ctrl     <- 0.32912 ; label("Lung distribution coefficient P, control (unitless)")           # Table 1, P, Lungs, Control
    kp_lung_lps30m   <- 0.53913 ; label("Lung distribution coefficient P, LPS 30 min (unitless)")        # Table 1, P, Lungs, LPS 30 min
    kp_lung_lps6h    <- 0.89661 ; label("Lung distribution coefficient P, LPS 6 h (unitless)")           # Table 1, P, Lungs, LPS 6 h
    kp_lung_lps24h   <- 0.68027 ; label("Lung distribution coefficient P, LPS 24 h (unitless)")          # Table 1, P, Lungs, LPS 24 h

    kp_liver_ctrl    <- 0.33850 ; label("Liver distribution coefficient P, control (unitless)")          # Table 1, P, Liver, Control
    kp_liver_lps30m  <- 0.15396 ; label("Liver distribution coefficient P, LPS 30 min (unitless)")       # Table 1, P, Liver, LPS 30 min
    kp_liver_lps6h   <- 0.18851 ; label("Liver distribution coefficient P, LPS 6 h (unitless)")          # Table 1, P, Liver, LPS 6 h
    kp_liver_lps24h  <- 0.12654 ; label("Liver distribution coefficient P, LPS 24 h (unitless)")         # Table 1, P, Liver, LPS 24 h

    kp_spleen_ctrl   <- 0.77655 ; label("Spleen distribution coefficient P, control (unitless)")         # Table 1, P, Spleen, Control
    kp_spleen_lps30m <- 0.54406 ; label("Spleen distribution coefficient P, LPS 30 min (unitless)")      # Table 1, P, Spleen, LPS 30 min
    kp_spleen_lps6h  <- 0.48332 ; label("Spleen distribution coefficient P, LPS 6 h (unitless)")         # Table 1, P, Spleen, LPS 6 h
    kp_spleen_lps24h <- 0.49543 ; label("Spleen distribution coefficient P, LPS 24 h (unitless)")        # Table 1, P, Spleen, LPS 24 h

    kp_kidney_ctrl   <- 0.26674 ; label("Kidney distribution coefficient P, control (unitless)")         # Table 1, P, Kidneys, Control
    kp_kidney_lps30m <- 0.18484 ; label("Kidney distribution coefficient P, LPS 30 min (unitless)")      # Table 1, P, Kidneys, LPS 30 min
    kp_kidney_lps6h  <- 0.22915 ; label("Kidney distribution coefficient P, LPS 6 h (unitless)")         # Table 1, P, Kidneys, LPS 6 h
    kp_kidney_lps24h <- 0.23102 ; label("Kidney distribution coefficient P, LPS 24 h (unitless)")        # Table 1, P, Kidneys, LPS 24 h

    kp_other_ctrl    <- 1.27424 ; label("Rest-of-body distribution coefficient P, control (unitless)")     # Table 1, P, Rest of body, Control
    kp_other_lps30m  <- 1.63449 ; label("Rest-of-body distribution coefficient P, LPS 30 min (unitless)")  # Table 1, P, Rest of body, LPS 30 min
    kp_other_lps6h   <- 2.61706 ; label("Rest-of-body distribution coefficient P, LPS 6 h (unitless)")     # Table 1, P, Rest of body, LPS 6 h
    kp_other_lps24h  <- 8.19636 ; label("Rest-of-body distribution coefficient P, LPS 24 h (unitless)")    # Table 1, P, Rest of body, LPS 24 h

    # ------------------------------------------------------------------
    # Permeability coefficients PAC -- Kutumova 2024 Table 1 (fitted, cohort-specific)
    # ------------------------------------------------------------------
    # PAC gates transcapillary diffusion in both directions
    # (Appendix A1: v2in = PAC * Q * Ab / Vb, v2out = PAC * Q * At / (Vt * P)).
    # The optimisation constrained PAC to be non-decreasing from control through
    # LPS 24 h for the lungs, liver, spleen, and kidneys (section 2.8,
    # equation 9); the rest-of-body PAC was unconstrained, which is why its
    # column is non-monotone.
    pa_lung_ctrl     <- 0.00010 ; label("Lung permeability coefficient PAC, control (unitless)")         # Table 1, PAC, Lungs, Control
    pa_lung_lps30m   <- 0.00400 ; label("Lung permeability coefficient PAC, LPS 30 min (unitless)")      # Table 1, PAC, Lungs, LPS 30 min
    pa_lung_lps6h    <- 0.00444 ; label("Lung permeability coefficient PAC, LPS 6 h (unitless)")         # Table 1, PAC, Lungs, LPS 6 h
    pa_lung_lps24h   <- 0.04834 ; label("Lung permeability coefficient PAC, LPS 24 h (unitless)")        # Table 1, PAC, Lungs, LPS 24 h

    pa_liver_ctrl    <- 0.00082 ; label("Liver permeability coefficient PAC, control (unitless)")        # Table 1, PAC, Liver, Control
    pa_liver_lps30m  <- 0.00083 ; label("Liver permeability coefficient PAC, LPS 30 min (unitless)")     # Table 1, PAC, Liver, LPS 30 min
    pa_liver_lps6h   <- 0.00090 ; label("Liver permeability coefficient PAC, LPS 6 h (unitless)")        # Table 1, PAC, Liver, LPS 6 h
    pa_liver_lps24h  <- 0.00099 ; label("Liver permeability coefficient PAC, LPS 24 h (unitless)")       # Table 1, PAC, Liver, LPS 24 h

    pa_spleen_ctrl   <- 0.03112 ; label("Spleen permeability coefficient PAC, control (unitless)")       # Table 1, PAC, Spleen, Control
    pa_spleen_lps30m <- 0.13851 ; label("Spleen permeability coefficient PAC, LPS 30 min (unitless)")    # Table 1, PAC, Spleen, LPS 30 min
    pa_spleen_lps6h  <- 0.31223 ; label("Spleen permeability coefficient PAC, LPS 6 h (unitless)")       # Table 1, PAC, Spleen, LPS 6 h
    pa_spleen_lps24h <- 0.58405 ; label("Spleen permeability coefficient PAC, LPS 24 h (unitless)")      # Table 1, PAC, Spleen, LPS 24 h

    pa_kidney_ctrl   <- 0.00178 ; label("Kidney permeability coefficient PAC, control (unitless)")       # Table 1, PAC, Kidneys, Control
    pa_kidney_lps30m <- 0.03171 ; label("Kidney permeability coefficient PAC, LPS 30 min (unitless)")    # Table 1, PAC, Kidneys, LPS 30 min
    pa_kidney_lps6h  <- 0.03176 ; label("Kidney permeability coefficient PAC, LPS 6 h (unitless)")       # Table 1, PAC, Kidneys, LPS 6 h
    pa_kidney_lps24h <- 0.06059 ; label("Kidney permeability coefficient PAC, LPS 24 h (unitless)")      # Table 1, PAC, Kidneys, LPS 24 h

    pa_other_ctrl    <- 0.65139 ; label("Rest-of-body permeability coefficient PAC, control (unitless)")     # Table 1, PAC, Rest of body, Control
    pa_other_lps30m  <- 0.17960 ; label("Rest-of-body permeability coefficient PAC, LPS 30 min (unitless)")  # Table 1, PAC, Rest of body, LPS 30 min
    pa_other_lps6h   <- 1.29770 ; label("Rest-of-body permeability coefficient PAC, LPS 6 h (unitless)")     # Table 1, PAC, Rest of body, LPS 6 h
    pa_other_lps24h  <- 0.00246 ; label("Rest-of-body permeability coefficient PAC, LPS 24 h (unitless)")    # Table 1, PAC, Rest of body, LPS 24 h

    # ------------------------------------------------------------------
    # Endocytic uptake / exocytic release -- Kutumova 2024 Table 2
    # (fitted, shared by all four cohorts)
    # ------------------------------------------------------------------
    # Uptake is a Hill function of ABSOLUTE TIME SINCE THE ANP DOSE, not of
    # concentration (Appendix A1 r3in: KRESUP = KRESmax * t^KRESn /
    # (KRES50^KRESn + t^KRESn), v3in = KRESUP * At). The model therefore
    # requires the ANP dose to be administered at t = 0.
    kres_max_lung   <- 0.09782  ; label("Lung maximum endocytic uptake rate KRESmax (1/min)")         # Table 2, KRESmax, Lungs
    kres_50_lung    <- 12.5763  ; label("Lung half-maximum uptake time KRES50 (min)")                 # Table 2, KRES 50, Lungs
    kres_n_lung     <- 1.69178  ; label("Lung uptake Hill coefficient KRESn (unitless)")              # Table 2, KRESn, Lungs
    kres_rel_lung   <- 0.00195  ; label("Lung exocytic release rate KRESrelease (1/min)")             # Table 2, KRESrelease, Lungs

    kres_max_liver  <- 4.17318  ; label("Liver maximum endocytic uptake rate KRESmax (1/min)")        # Table 2, KRESmax, Liver
    kres_50_liver   <- 162.343  ; label("Liver half-maximum uptake time KRES50 (min)")                # Table 2, KRES 50, Liver
    kres_n_liver    <- 1.15243  ; label("Liver uptake Hill coefficient KRESn (unitless)")             # Table 2, KRESn, Liver
    kres_rel_liver  <- 0.02254  ; label("Liver exocytic release rate KRESrelease (1/min)")            # Table 2, KRESrelease, Liver

    kres_max_spleen <- 20.8142  ; label("Spleen maximum endocytic uptake rate KRESmax (1/min)")       # Table 2, KRESmax, Spleen
    kres_50_spleen  <- 138.387  ; label("Spleen half-maximum uptake time KRES50 (min)")               # Table 2, KRES 50, Spleen
    kres_n_spleen   <- 1.97529  ; label("Spleen uptake Hill coefficient KRESn (unitless)")            # Table 2, KRESn, Spleen
    kres_rel_spleen <- 0.09512  ; label("Spleen exocytic release rate KRESrelease (1/min)")           # Table 2, KRESrelease, Spleen

    kres_max_kidney <- 0.77376  ; label("Kidney maximum endocytic uptake rate KRESmax (1/min)")       # Table 2, KRESmax, Kidneys
    kres_50_kidney  <- 345.945  ; label("Kidney half-maximum uptake time KRES50 (min)")               # Table 2, KRES 50, Kidneys
    kres_n_kidney   <- 0.65725  ; label("Kidney uptake Hill coefficient KRESn (unitless)")            # Table 2, KRESn, Kidneys
    kres_rel_kidney <- 0.00302  ; label("Kidney exocytic release rate KRESrelease (1/min)")           # Table 2, KRESrelease, Kidneys

    kres_max_other  <- 0.81059  ; label("Rest-of-body maximum endocytic uptake rate KRESmax (1/min)") # Table 2, KRESmax, Rest of body
    kres_50_other   <- 53.6192  ; label("Rest-of-body half-maximum uptake time KRES50 (min)")         # Table 2, KRES 50, Rest of body
    kres_n_other    <- 25.4498  ; label("Rest-of-body uptake Hill coefficient KRESn (unitless)")      # Table 2, KRESn, Rest of body
    kres_rel_other  <- 0.16775  ; label("Rest-of-body exocytic release rate KRESrelease (1/min)")     # Table 2, KRESrelease, Rest of body

    # ------------------------------------------------------------------
    # Excretion -- Kutumova 2024 Table 2 / Appendix A5, A6
    # ------------------------------------------------------------------
    # Table 2 prints the unit as "L/min", but the rate laws multiply by BW^0.75
    # (vb = KbileC * BW^0.75 * ALt/VLt; vu = KurineC * BW^0.75 * AKb/VKb), so
    # the effective unit is L/min/kg^0.75 -- the same loose convention the paper
    # uses for QCC. Cheng 2020 writes this explicitly as
    # `Kbile = KbileC*BW**0.75 ; L/h`, confirming the reading.
    kbile_c  <- 3.4e-06 ; label("Biliary excretion coefficient KbileC (L/min/kg^0.75)")   # Table 2, KbileC, Liver (printed "L/min"; see note above)
    kurine_c <- 2.7e-06 ; label("Urinary excretion coefficient KurineC (L/min/kg^0.75)")  # Table 2, KurineC, Kidneys (printed "L/min"; see note above)

    # ------------------------------------------------------------------
    # Observation scaling and dosing
    # ------------------------------------------------------------------
    # Appendix A10: TRELu = k * CLuT, TRES = k * CST, TREK = k * CKT,
    # TREL = k * CLT. A single fitted constant converts every organ's
    # simulated tissue concentration to the measured average total radiant
    # efficiency per luminous area.
    k_tre <- 1358630.5 ; label("Tissue concentration to total radiant efficiency scale factor k (unitless)")  # Table 2, k

    # Kutumova 2024 gives the intravenous input as a zero-order infusion over
    # Timeiv (Appendix A8) but never states Timeiv. It is one of the framework
    # values inherited from Cheng 2020, whose model code sets
    # `Timeiv = 0.005 ; IV infusion time (h), set, approximately 15-20 seconds,
    # on average 18 sec` -- 0.005 h = 0.3 min. See vignette Assumptions.
    dur_iv <- fixed(0.3) ; label("Intravenous infusion duration Timeiv (min)")  # Cheng 2020 SI code Timeiv = 0.005 h = 0.3 min; not reported by Kutumova 2024
  })

  model({
    # ================================================================
    # 1. Cohort selector -- Kutumova 2024 Table 1
    # ================================================================
    # Exactly one of the three LPS indicators is 1, or all are 0 for the
    # LPS-naive control arm.
    f_ctrl <- 1 - STUDY_LPS30M - STUDY_LPS6H - STUDY_LPS24H

    kp_lung <- kp_lung_ctrl * f_ctrl + kp_lung_lps30m * STUDY_LPS30M +
      kp_lung_lps6h * STUDY_LPS6H + kp_lung_lps24h * STUDY_LPS24H
    kp_liver <- kp_liver_ctrl * f_ctrl + kp_liver_lps30m * STUDY_LPS30M +
      kp_liver_lps6h * STUDY_LPS6H + kp_liver_lps24h * STUDY_LPS24H
    kp_spleen <- kp_spleen_ctrl * f_ctrl + kp_spleen_lps30m * STUDY_LPS30M +
      kp_spleen_lps6h * STUDY_LPS6H + kp_spleen_lps24h * STUDY_LPS24H
    kp_kidney <- kp_kidney_ctrl * f_ctrl + kp_kidney_lps30m * STUDY_LPS30M +
      kp_kidney_lps6h * STUDY_LPS6H + kp_kidney_lps24h * STUDY_LPS24H
    kp_other <- kp_other_ctrl * f_ctrl + kp_other_lps30m * STUDY_LPS30M +
      kp_other_lps6h * STUDY_LPS6H + kp_other_lps24h * STUDY_LPS24H

    pa_lung <- pa_lung_ctrl * f_ctrl + pa_lung_lps30m * STUDY_LPS30M +
      pa_lung_lps6h * STUDY_LPS6H + pa_lung_lps24h * STUDY_LPS24H
    pa_liver <- pa_liver_ctrl * f_ctrl + pa_liver_lps30m * STUDY_LPS30M +
      pa_liver_lps6h * STUDY_LPS6H + pa_liver_lps24h * STUDY_LPS24H
    pa_spleen <- pa_spleen_ctrl * f_ctrl + pa_spleen_lps30m * STUDY_LPS30M +
      pa_spleen_lps6h * STUDY_LPS6H + pa_spleen_lps24h * STUDY_LPS24H
    pa_kidney <- pa_kidney_ctrl * f_ctrl + pa_kidney_lps30m * STUDY_LPS30M +
      pa_kidney_lps6h * STUDY_LPS6H + pa_kidney_lps24h * STUDY_LPS24H
    pa_other <- pa_other_ctrl * f_ctrl + pa_other_lps30m * STUDY_LPS30M +
      pa_other_lps6h * STUDY_LPS6H + pa_other_lps24h * STUDY_LPS24H

    # ================================================================
    # 2. Physiology -- Kutumova 2024 Appendix A1, A8, A9, A10
    # ================================================================
    # Q = QC * QCC * BW^0.75 (Appendix A1 r1)
    q_co     <- qcc * WT^0.75
    q_lung   <- qc_lung   * q_co
    q_liver  <- qc_liver  * q_co
    q_spleen <- qc_spleen * q_co
    q_kidney <- qc_kidney * q_co
    # QRC = 1 - (QLC + QSC + QKC) (Appendix A10)
    q_other  <- (1 - qc_liver - qc_spleen - qc_kidney) * q_co

    # V = BW * VC; Vb = V * BV; Vt = V - Vb (Appendix A1 r2in, r2out)
    v_lung   <- WT * vc_lung
    v_liver  <- WT * vc_liver
    v_spleen <- WT * vc_spleen
    v_kidney <- WT * vc_kidney
    # VRC = 1 - (VLC + VSC + VKC + VLuC + VPlC) (Appendix A10)
    v_other  <- WT * (1 - vc_liver - vc_spleen - vc_kidney - vc_lung - vc_plasma)

    vb_lung   <- v_lung   * bv_lung
    vb_liver  <- v_liver  * bv_liver
    vb_spleen <- v_spleen * bv_spleen
    vb_kidney <- v_kidney * bv_kidney
    vb_other  <- v_other  * bv_other

    vt_lung   <- v_lung   - vb_lung
    vt_liver  <- v_liver  - vb_liver
    vt_spleen <- v_spleen - vb_spleen
    vt_kidney <- v_kidney - vb_kidney
    vt_other  <- v_other  - vb_other

    # VPv = BW * VPlC * 0.8 (Appendix A8); APv = 0.25 * VPv (Appendix A9)
    v_venous   <- WT * vc_plasma * 0.8
    v_arterial <- 0.25 * v_venous

    c_venous   <- venous   / v_venous     # CV, Appendix A8
    c_arterial <- arterial / v_arterial   # CA, Appendix A9

    # ================================================================
    # 3. Time-dependent endocytic uptake rate constants
    # ================================================================
    # KRESUP = KRESmax * t^KRESn / (KRES50^KRESn + t^KRESn), Appendix A1 r3in.
    # `t` is time since the ANP dose, so the dose must be at t = 0.
    kup_lung   <- kres_max_lung   * t^kres_n_lung /
      (kres_50_lung^kres_n_lung + t^kres_n_lung)
    kup_liver  <- kres_max_liver  * t^kres_n_liver /
      (kres_50_liver^kres_n_liver + t^kres_n_liver)
    kup_spleen <- kres_max_spleen * t^kres_n_spleen /
      (kres_50_spleen^kres_n_spleen + t^kres_n_spleen)
    kup_kidney <- kres_max_kidney * t^kres_n_kidney /
      (kres_50_kidney^kres_n_kidney + t^kres_n_kidney)
    kup_other  <- kres_max_other  * t^kres_n_other /
      (kres_50_other^kres_n_other + t^kres_n_other)

    # ================================================================
    # 4. Reaction rates -- Kutumova 2024 Appendix A1
    # ================================================================
    # Lungs sit in series: input is venous plasma, output is arterial plasma.
    r1_lung    <- q_lung * c_venous                        # r1
    r2in_lung  <- pa_lung * q_lung * vp_lung / vb_lung      # r2in
    r2out_lung <- pa_lung * q_lung * is_lung / (vt_lung * kp_lung)  # r2out
    r3in_lung  <- kup_lung * is_lung                        # r3in
    r3out_lung <- kres_rel_lung * int_lung                  # r3out
    r4_lung    <- q_lung * vp_lung / vb_lung                # r4

    # Spleen, kidneys, and rest of body sit in parallel: input is arterial
    # plasma, output is venous plasma.
    r1_spleen    <- q_spleen * c_arterial
    r2in_spleen  <- pa_spleen * q_spleen * vp_spleen / vb_spleen
    r2out_spleen <- pa_spleen * q_spleen * is_spleen / (vt_spleen * kp_spleen)
    r3in_spleen  <- kup_spleen * is_spleen
    r3out_spleen <- kres_rel_spleen * int_spleen
    r4_spleen    <- q_spleen * vp_spleen / vb_spleen

    r1_kidney    <- q_kidney * c_arterial
    r2in_kidney  <- pa_kidney * q_kidney * vp_kidney / vb_kidney
    r2out_kidney <- pa_kidney * q_kidney * is_kidney / (vt_kidney * kp_kidney)
    r3in_kidney  <- kup_kidney * is_kidney
    r3out_kidney <- kres_rel_kidney * int_kidney
    r4_kidney    <- q_kidney * vp_kidney / vb_kidney
    # vu = KurineC * BW^0.75 * AKb / VKb (Appendix A6) -- from capillary blood
    r_urine      <- kurine_c * WT^0.75 * vp_kidney / vb_kidney

    r1_other    <- q_other * c_arterial
    r2in_other  <- pa_other * q_other * vp_other / vb_other
    r2out_other <- pa_other * q_other * is_other / (vt_other * kp_other)
    r3in_other  <- kup_other * is_other
    r3out_other <- kres_rel_other * int_other
    r4_other    <- q_other * vp_other / vb_other

    # The liver differs (Appendix A5, Fig. 1c): Kupffer cells phagocytose
    # directly from CAPILLARY BLOOD rather than from the interstitium, so
    # r3in / r3out act on vp_liver, and biliary excretion drains the
    # interstitium.
    r1_liver    <- q_liver * c_arterial
    r2in_liver  <- pa_liver * q_liver * vp_liver / vb_liver
    r2out_liver <- pa_liver * q_liver * is_liver / (vt_liver * kp_liver)
    r3in_liver  <- kup_liver * vp_liver
    r3out_liver <- kres_rel_liver * int_liver
    r4_liver    <- q_liver * vp_liver / vb_liver
    # vb = KbileC * BW^0.75 * ALt / VLt (Appendix A5)
    r_bile      <- kbile_c * WT^0.75 * is_liver / vt_liver

    # ================================================================
    # 5. ODE system
    # ================================================================
    # Venous plasma collects the venous outflow of every parallel organ and
    # supplies the lungs; the intravenous dose lands here (Appendix A8).
    d/dt(venous) <- r4_liver + r4_spleen + r4_kidney + r4_other - r1_lung
    # Arterial plasma is filled by the lungs and drains to the parallel organs
    # (Appendix A9).
    d/dt(arterial) <- r4_lung - r1_liver - r1_spleen - r1_kidney - r1_other

    # Common organ module (Appendix A2)
    d/dt(vp_lung)  <- r1_lung - r2in_lung + r2out_lung - r4_lung
    d/dt(is_lung)  <- r2in_lung - r2out_lung - r3in_lung + r3out_lung
    d/dt(int_lung) <- r3in_lung - r3out_lung

    d/dt(vp_spleen)  <- r1_spleen - r2in_spleen + r2out_spleen - r4_spleen
    d/dt(is_spleen)  <- r2in_spleen - r2out_spleen - r3in_spleen + r3out_spleen
    d/dt(int_spleen) <- r3in_spleen - r3out_spleen

    # Liver module (Appendix A5)
    d/dt(vp_liver)  <- r1_liver - r2in_liver + r2out_liver -
      r3in_liver + r3out_liver - r4_liver
    d/dt(is_liver)  <- r2in_liver - r2out_liver - r_bile
    d/dt(int_liver) <- r3in_liver - r3out_liver

    # Kidney module (Appendix A6)
    d/dt(vp_kidney)  <- r1_kidney - r2in_kidney + r2out_kidney -
      r4_kidney - r_urine
    d/dt(is_kidney)  <- r2in_kidney - r2out_kidney - r3in_kidney + r3out_kidney
    d/dt(int_kidney) <- r3in_kidney - r3out_kidney

    # Rest of body (Appendix A7)
    d/dt(vp_other)  <- r1_other - r2in_other + r2out_other - r4_other
    d/dt(is_other)  <- r2in_other - r2out_other - r3in_other + r3out_other
    d/dt(int_other) <- r3in_other - r3out_other

    # Excretion accumulators (for mass balance; not in the paper's state list)
    d/dt(urine) <- r_urine
    d/dt(bile)  <- r_bile

    # Zero-order intravenous input over Timeiv (Appendix A8). Requires
    # rate = -2 on the dose row so rxode2 takes the duration from the model.
    dur(venous) <- dur_iv

    # ================================================================
    # 6. Observations
    # ================================================================
    # Tissue concentration Ctissue = (At + ARES) / Vt (Appendix A2). For the
    # liver, ALivert = ALt + ALRES and CLivert = ALivert / VLt (Appendix A5) --
    # the internalised pool counts toward the tissue signal even though the
    # liver takes it up from blood.
    Clung   <- (is_lung   + int_lung)   / vt_lung
    Cspleen <- (is_spleen + int_spleen) / vt_spleen
    Ckidney <- (is_kidney + int_kidney) / vt_kidney
    Cliver  <- (is_liver  + int_liver)  / vt_liver

    # Total organ concentration Ctotal = (Ab + At + ARES) / V (Appendix A2)
    Clung_total   <- (vp_lung   + is_lung   + int_lung)   / v_lung
    Cspleen_total <- (vp_spleen + is_spleen + int_spleen) / v_spleen
    Ckidney_total <- (vp_kidney + is_kidney + int_kidney) / v_kidney
    Cliver_total  <- (vp_liver  + is_liver  + int_liver)  / v_liver

    # Measured observable: average total radiant efficiency per luminous area
    # (Appendix A10). Reproduces Supplementary Figure 1.
    TRElung   <- k_tre * Clung
    TREspleen <- k_tre * Cspleen
    TREkidney <- k_tre * Ckidney
    TREliver  <- k_tre * Cliver

    # Liver-normalised signal -- the y axis of Figure 4. The scale factor k
    # cancels, so these ratios are independent of k.
    TREratio_lung   <- Clung   / Cliver
    TREratio_spleen <- Cspleen / Cliver
    TREratio_kidney <- Ckidney / Cliver

    # Venous plasma concentration. Not measured in this study (no plasma
    # sampling); provided as the canonical PK observation.
    Cc <- c_venous

    # Total nanoparticle mass in the system, for the mass-balance check in the
    # vignette (zero the excretion coefficients and this must stay constant
    # after the infusion ends).
    amt_total <- venous + arterial +
      vp_lung + is_lung + int_lung +
      vp_spleen + is_spleen + int_spleen +
      vp_liver + is_liver + int_liver +
      vp_kidney + is_kidney + int_kidney +
      vp_other + is_other + int_other +
      urine + bile
  })
}
