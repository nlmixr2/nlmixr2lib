Elmokadem_2023_mavoglurant_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 16 compartments; Bayesian, fit with R/Stan/Torsten",
    "and Julia/SciML/Turing.jl). Mavoglurant disposition after a single",
    "short intravenous infusion in 20 healthy adults (study A2121).",
    "Perfusion-limited (well-stirred) organ compartments - lung, heart,",
    "brain, muscle, adipose, skin, spleen, pancreas, liver, stomach, gut,",
    "bone, kidney and a lumped rest-of-body - plus arterial and venous",
    "blood. The lung sits in series between venous and arterial blood so",
    "the whole cardiac output passes through it; spleen, pancreas,",
    "stomach and gut drain into the liver, which is the only eliminating",
    "organ (unbound intrinsic clearance on the well-stirred outflow",
    "concentration). Every physiological parameter is fixed and driven by",
    "body weight: cardiac output scales as 187 * WT^0.81 mL/min, each",
    "organ blood flow is a fixed fraction of cardiac output, and each",
    "organ volume is a fixed fraction of body weight divided by that",
    "organ's density. Six drug-specific parameters were estimated by full",
    "Bayesian inference - the intrinsic clearance and the brain, muscle,",
    "adipose, bone and rest-of-body tissue:plasma partition coefficients -",
    "with lognormal between-subject variability on intrinsic clearance",
    "only and a lognormal residual error. The remaining nine partition",
    "coefficients were held fixed. Values here are the posterior medians",
    "of the Stan/Torsten general-ODE fit; the linear-ODE and Turing.jl",
    "fits of the same model agree to within about 1%."
  )
  reference <- paste(
    "Elmokadem A, Zhang Y, Knab T, Jordie E, Gillespie WR. Bayesian PBPK",
    "modeling using R/Stan/Torsten and Julia/SciML/Turing.jl.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(3):300-310.",
    "doi:10.1002/psp4.12926.",
    "Structural and fixed physiological parameters are taken from the",
    "paper's designated online supplement, the model file",
    "model/mavoPBPKGenODE.stan in",
    "https://github.com/metrumresearchgroup/BayesPBPK-tutorial,",
    "which the paper names as the reference for those values.",
    "The PBPK structure follows Wendling T, Dumitras S, Ogungbenro K,",
    "Aarons L. J Pharmacokinet Pharmacodyn. 2015;42(6):639-657.",
    "doi:10.1007/s10928-015-9430-4.",
    sep = " "
  )
  vignette <- "Elmokadem_2023_mavoglurant_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate. Drives the entire physiology: cardiac output",
        "is 187 * WT^0.81 mL/min and every organ volume is a fixed",
        "fraction of WT divided by that organ's density. Observed range",
        "61.9-103.5 kg over the 20 analysed subjects."
      ),
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    lung     = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    heart    = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    brain    = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    adipose  = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    skin     = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    spleen   = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    stomach  = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    gut      = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    bone     = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    other    = list(analyte = "mavoglurant", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "mavoglurant", units = "mg", specimen = "whole blood", verified = TRUE),
    venous   = list(analyte = "mavoglurant", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    weight_range   = "61.9-103.5 kg",
    weight_median  = "78.9 kg",
    age_range      = "19-50 years",
    age_median     = "30 years",
    disease_state  = "healthy volunteers",
    dose_range     = paste(
      "single i.v. infusion of 25 mg (infusion rate 75 or 150 mg/h) or",
      "37.5 mg (infusion rate 225 mg/h); nominal 10 min infusion"
    ),
    notes          = paste(
      "DATA section of the paper: the mavoglurant PK data are from study",
      "A2121, a healthy-volunteer study of a 10 min i.v. infusion, shared",
      "publicly as a csv by the nlmixr team. The dataset holds 120",
      "subjects; the first 20 were analysed here, contributing 268",
      "concentration observations over 24-48 h. Body weight is the only",
      "covariate used by the model. The paper itself reports no other",
      "demographics; the age range above is computed from the AGE column",
      "of the analysis dataset in the paper's online supplement",
      "(data/Mavoglurant_A2121_nmpk.csv). That dataset also carries SEX",
      "(16 subjects coded 1, 4 coded 2) and RACE (12 coded 1, 7 coded 2,",
      "1 coded 88) columns, but neither coding is documented in the paper",
      "or the supplement, so no sex or race breakdown is asserted here."
    )
  )

  ini({
    # ---- Estimated drug-specific parameters -------------------------------
    # Table 1, "Stan/Torsten (general ODE solver)" column: posterior medians.
    lclint      <- log(1395)  ; label("Unbound intrinsic hepatic clearance (L/h)")                 # Table 1 CLintHat 1395 (90% CI 1239, 1604)
    lkp_brain   <- log(4.83)  ; label("Brain:plasma partition coefficient (unitless)")             # Table 1 KbBR 4.83 (3.64, 6.36)
    lkp_muscle  <- log(1.75)  ; label("Muscle:plasma partition coefficient (unitless)")            # Table 1 KbMU 1.75 (1.49, 2.04)
    lkp_adipose <- log(10.3)  ; label("Adipose:plasma partition coefficient (unitless)")           # Table 1 KbAD 10.3 (9.31, 11.51)
    lkp_bone    <- log(1.03)  ; label("Bone:plasma partition coefficient (unitless)")              # Table 1 KbBO 1.03 (0.74, 1.47)
    lkp_other   <- log(1.77)  ; label("Rest-of-body:plasma partition coefficient (unitless)")      # Table 1 KbRB 1.77 (1.26, 2.46)

    # ---- Fixed partition coefficients -------------------------------------
    # PBPK MODEL section: "several tissue:plasma partition coefficients (Kb)"
    # were fixed. The supplement writes each as exp(x), so x is the value
    # transcribed here on the log scale.
    lkp_lung     <- fixed(0.8334)  ; label("Lung:plasma partition coefficient, log scale")         # supplement mavoPBPKGenODE.stan: KbLU = exp(0.8334)
    lkp_heart    <- fixed(1.1205)  ; label("Heart:plasma partition coefficient, log scale")        # supplement mavoPBPKGenODE.stan: KbHT = exp(1.1205)
    lkp_skin     <- fixed(-0.5238) ; label("Skin:plasma partition coefficient, log scale")         # supplement mavoPBPKGenODE.stan: KbSK = exp(-0.5238)
    lkp_spleen   <- fixed(0.3224)  ; label("Spleen:plasma partition coefficient, log scale")       # supplement mavoPBPKGenODE.stan: KbSP = exp(0.3224)
    lkp_pancreas <- fixed(0.3224)  ; label("Pancreas:plasma partition coefficient, log scale")     # supplement mavoPBPKGenODE.stan: KbPA = exp(0.3224)
    lkp_liver    <- fixed(1.7604)  ; label("Liver:plasma partition coefficient, log scale")        # supplement mavoPBPKGenODE.stan: KbLI = exp(1.7604)
    lkp_stomach  <- fixed(0.3224)  ; label("Stomach:plasma partition coefficient, log scale")      # supplement mavoPBPKGenODE.stan: KbST = exp(0.3224)
    lkp_gut      <- fixed(1.2026)  ; label("Gut:plasma partition coefficient, log scale")          # supplement mavoPBPKGenODE.stan: KbGU = exp(1.2026)
    lkp_kidney   <- fixed(1.3171)  ; label("Kidney:plasma partition coefficient, log scale")       # supplement mavoPBPKGenODE.stan: KbKI = exp(1.3171)

    # ---- Blood binding ----------------------------------------------------
    bp <- fixed(0.61)  ; label("Blood:plasma concentration ratio (unitless)")                      # supplement mavoPBPKGenODE.stan: BP = 0.61
    fu <- fixed(0.028) ; label("Fraction unbound in plasma (unitless)")                            # supplement mavoPBPKGenODE.stan: fup = 0.028

    # ---- Cardiac output ---------------------------------------------------
    qc_coef <- fixed(187)  ; label("Cardiac output coefficient (mL/min at 1 kg)")                  # supplement mavoPBPKGenODE.stan: CO = (187.00*WT^0.81)*60/1000, from White et al (1968)
    e_wt_qc <- fixed(0.81) ; label("Body-weight exponent on cardiac output (unitless)")            # supplement mavoPBPKGenODE.stan: CO = (187.00*WT^0.81)*60/1000

    # ---- Fractional organ blood flows (fraction of cardiac output) --------
    # The supplement writes each as <pct>*CO/100; the fractions are recorded here.
    fq_heart    <- fixed(0.040) ; label("Fractional blood flow, heart (fraction of cardiac output)")     # supplement: QHT = 4.0*CO/100
    fq_brain    <- fixed(0.120) ; label("Fractional blood flow, brain (fraction of cardiac output)")     # supplement: QBR = 12.0*CO/100
    fq_muscle   <- fixed(0.170) ; label("Fractional blood flow, muscle (fraction of cardiac output)")    # supplement: QMU = 17.0*CO/100
    fq_adipose  <- fixed(0.050) ; label("Fractional blood flow, adipose (fraction of cardiac output)")   # supplement: QAD = 5.0*CO/100
    fq_skin     <- fixed(0.050) ; label("Fractional blood flow, skin (fraction of cardiac output)")      # supplement: QSK = 5.0*CO/100
    fq_spleen   <- fixed(0.030) ; label("Fractional blood flow, spleen (fraction of cardiac output)")    # supplement: QSP = 3.0*CO/100
    fq_pancreas <- fixed(0.010) ; label("Fractional blood flow, pancreas (fraction of cardiac output)")  # supplement: QPA = 1.0*CO/100
    fq_liver    <- fixed(0.255) ; label("Fractional blood flow, total liver (fraction of cardiac output)") # supplement: QLI = 25.5*CO/100
    fq_stomach  <- fixed(0.010) ; label("Fractional blood flow, stomach (fraction of cardiac output)")   # supplement: QST = 1.0*CO/100
    fq_gut      <- fixed(0.140) ; label("Fractional blood flow, gut (fraction of cardiac output)")       # supplement: QGU = 14.0*CO/100
    fq_bone     <- fixed(0.050) ; label("Fractional blood flow, bone (fraction of cardiac output)")      # supplement: QBO = 5.0*CO/100
    fq_kidney   <- fixed(0.190) ; label("Fractional blood flow, kidney (fraction of cardiac output)")    # supplement: QKI = 19.0*CO/100

    # ---- Fractional organ volumes (fraction of body weight) ---------------
    # The supplement writes each as (<pct>*WT/100)/<density>; the volume
    # fraction and the density are recorded separately so both are auditable.
    fvol_lung     <- fixed(0.0076) ; label("Fractional volume, lung (kg per kg body weight)")            # supplement: VLU = (0.76*WT/100)/1.051
    fvol_heart    <- fixed(0.0047) ; label("Fractional volume, heart (kg per kg body weight)")           # supplement: VHT = (0.47*WT/100)/1.030
    fvol_brain    <- fixed(0.0200) ; label("Fractional volume, brain (kg per kg body weight)")           # supplement: VBR = (2.00*WT/100)/1.036
    fvol_muscle   <- fixed(0.4000) ; label("Fractional volume, muscle (kg per kg body weight)")          # supplement: VMU = (40.00*WT/100)/1.041
    fvol_adipose  <- fixed(0.2142) ; label("Fractional volume, adipose (kg per kg body weight)")         # supplement: VAD = (21.42*WT/100)/0.916
    fvol_skin     <- fixed(0.0371) ; label("Fractional volume, skin (kg per kg body weight)")            # supplement: VSK = (3.71*WT/100)/1.116
    fvol_spleen   <- fixed(0.0026) ; label("Fractional volume, spleen (kg per kg body weight)")          # supplement: VSP = (0.26*WT/100)/1.054
    fvol_pancreas <- fixed(0.0014) ; label("Fractional volume, pancreas (kg per kg body weight)")        # supplement: VPA = (0.14*WT/100)/1.045
    fvol_liver    <- fixed(0.0257) ; label("Fractional volume, liver (kg per kg body weight)")           # supplement: VLI = (2.57*WT/100)/1.040
    fvol_stomach  <- fixed(0.0021) ; label("Fractional volume, stomach (kg per kg body weight)")         # supplement: VST = (0.21*WT/100)/1.050
    fvol_gut      <- fixed(0.0144) ; label("Fractional volume, gut (kg per kg body weight)")             # supplement: VGU = (1.44*WT/100)/1.043
    fvol_bone     <- fixed(0.1429) ; label("Fractional volume, bone (kg per kg body weight)")            # supplement: VBO = (14.29*WT/100)/1.990
    fvol_kidney   <- fixed(0.0044) ; label("Fractional volume, kidney (kg per kg body weight)")          # supplement: VKI = (0.44*WT/100)/1.050
    fvol_arterial <- fixed(0.0281) ; label("Fractional volume, arterial blood (kg per kg body weight)")  # supplement: VAB = (2.81*WT/100)/1.040
    fvol_venous   <- fixed(0.0562) ; label("Fractional volume, venous blood (kg per kg body weight)")    # supplement: VVB = (5.62*WT/100)/1.040
    fvol_other    <- fixed(0.0386) ; label("Fractional volume, rest of body (kg per kg body weight)")    # supplement: VRB = (3.86*WT/100)/1.040

    # ---- Organ densities --------------------------------------------------
    density_lung     <- fixed(1.051) ; label("Tissue density, lung (kg/L)")            # supplement: VLU = (0.76*WT/100)/1.051
    density_heart    <- fixed(1.030) ; label("Tissue density, heart (kg/L)")           # supplement: VHT = (0.47*WT/100)/1.030
    density_brain    <- fixed(1.036) ; label("Tissue density, brain (kg/L)")           # supplement: VBR = (2.00*WT/100)/1.036
    density_muscle   <- fixed(1.041) ; label("Tissue density, muscle (kg/L)")          # supplement: VMU = (40.00*WT/100)/1.041
    density_adipose  <- fixed(0.916) ; label("Tissue density, adipose (kg/L)")         # supplement: VAD = (21.42*WT/100)/0.916
    density_skin     <- fixed(1.116) ; label("Tissue density, skin (kg/L)")            # supplement: VSK = (3.71*WT/100)/1.116
    density_spleen   <- fixed(1.054) ; label("Tissue density, spleen (kg/L)")          # supplement: VSP = (0.26*WT/100)/1.054
    density_pancreas <- fixed(1.045) ; label("Tissue density, pancreas (kg/L)")        # supplement: VPA = (0.14*WT/100)/1.045
    density_liver    <- fixed(1.040) ; label("Tissue density, liver (kg/L)")           # supplement: VLI = (2.57*WT/100)/1.040
    density_stomach  <- fixed(1.050) ; label("Tissue density, stomach (kg/L)")         # supplement: VST = (0.21*WT/100)/1.050
    density_gut      <- fixed(1.043) ; label("Tissue density, gut (kg/L)")             # supplement: VGU = (1.44*WT/100)/1.043
    density_bone     <- fixed(1.990) ; label("Tissue density, bone (kg/L)")            # supplement: VBO = (14.29*WT/100)/1.990
    density_kidney   <- fixed(1.050) ; label("Tissue density, kidney (kg/L)")          # supplement: VKI = (0.44*WT/100)/1.050
    density_arterial <- fixed(1.040) ; label("Blood density, arterial (kg/L)")         # supplement: VAB = (2.81*WT/100)/1.040
    density_venous   <- fixed(1.040) ; label("Blood density, venous (kg/L)")           # supplement: VVB = (5.62*WT/100)/1.040
    density_other    <- fixed(1.040) ; label("Tissue density, rest of body (kg/L)")    # supplement: VRB = (3.86*WT/100)/1.040

    # ---- Between-subject variability --------------------------------------
    # Intrinsic clearance is the only parameter carrying IIV (STATISTICAL
    # MODEL section). Table 1 reports omega[1] = 0.358 as the standard
    # deviation of the log-transformed CLint, so the variance is 0.358^2.
    etalclint ~ 0.128164  # Table 1 omega[1] 0.358 (0.263, 0.499), SD of log CLint; 0.358^2 = 0.128164

    # ---- Residual error ---------------------------------------------------
    # STATISTICAL MODEL section: observed concentrations are lognormally
    # distributed about the prediction with variance sigma^2 on the log scale
    # (supplement: logCObs ~ normal(log(cHatObs), sigma)).
    expSd <- 0.32 ; label("Lognormal residual error, SD on the log scale")  # Table 1 sigma 0.32 (0.298, 0.345)
  })

  model({
    # ---- 1. Physiology derived from body weight ---------------------------
    # Cardiac output, mL/min converted to L/h.
    qc <- qc_coef * WT^e_wt_qc * 60 / 1000

    q_heart    <- fq_heart    * qc
    q_brain    <- fq_brain    * qc
    q_muscle   <- fq_muscle   * qc
    q_adipose  <- fq_adipose  * qc
    q_skin     <- fq_skin     * qc
    q_spleen   <- fq_spleen   * qc
    q_pancreas <- fq_pancreas * qc
    q_liver    <- fq_liver    * qc
    q_stomach  <- fq_stomach  * qc
    q_gut      <- fq_gut      * qc
    q_bone     <- fq_bone     * qc
    q_kidney   <- fq_kidney   * qc

    # Hepatic artery flow is the part of total liver flow not supplied by the
    # portal drainage of spleen, pancreas, stomach and gut.
    q_liver_arterial <- q_liver - (q_spleen + q_pancreas + q_stomach + q_gut)

    # Rest of body takes the cardiac output not claimed by the organs that
    # drain directly to venous blood.
    q_other <- qc - (q_heart + q_brain + q_muscle + q_adipose + q_skin +
                       q_liver + q_bone + q_kidney)

    # The lung is in series: the whole cardiac output passes through it.
    q_lung <- q_heart + q_brain + q_muscle + q_adipose + q_skin + q_liver +
      q_bone + q_kidney + q_other

    v_lung     <- fvol_lung     * WT / density_lung
    v_heart    <- fvol_heart    * WT / density_heart
    v_brain    <- fvol_brain    * WT / density_brain
    v_muscle   <- fvol_muscle   * WT / density_muscle
    v_adipose  <- fvol_adipose  * WT / density_adipose
    v_skin     <- fvol_skin     * WT / density_skin
    v_spleen   <- fvol_spleen   * WT / density_spleen
    v_pancreas <- fvol_pancreas * WT / density_pancreas
    v_liver    <- fvol_liver    * WT / density_liver
    v_stomach  <- fvol_stomach  * WT / density_stomach
    v_gut      <- fvol_gut      * WT / density_gut
    v_bone     <- fvol_bone     * WT / density_bone
    v_kidney   <- fvol_kidney   * WT / density_kidney
    v_arterial <- fvol_arterial * WT / density_arterial
    v_venous   <- fvol_venous   * WT / density_venous
    v_other    <- fvol_other    * WT / density_other

    # ---- 2. Individual drug-specific parameters ---------------------------
    clint <- exp(lclint + etalclint)

    kp_lung     <- exp(lkp_lung)
    kp_heart    <- exp(lkp_heart)
    kp_brain    <- exp(lkp_brain)
    kp_muscle   <- exp(lkp_muscle)
    kp_adipose  <- exp(lkp_adipose)
    kp_skin     <- exp(lkp_skin)
    kp_spleen   <- exp(lkp_spleen)
    kp_pancreas <- exp(lkp_pancreas)
    kp_liver    <- exp(lkp_liver)
    kp_stomach  <- exp(lkp_stomach)
    kp_gut      <- exp(lkp_gut)
    kp_bone     <- exp(lkp_bone)
    kp_kidney   <- exp(lkp_kidney)
    kp_other    <- exp(lkp_other)

    # Fraction unbound in blood; hepatic elimination is defined on unbound
    # drug in blood because the partition coefficients are tissue:plasma.
    fub <- fu / bp

    # ---- 3. Tissue concentrations -----------------------------------------
    # c_<organ> is the total tissue concentration; the well-stirred venous
    # outflow concentration is c_<organ> / kp_<organ>. The blood compartments
    # carry no partition coefficient.
    c_arterial <- arterial / v_arterial
    c_venous   <- venous   / v_venous
    c_lung     <- lung     / v_lung
    c_heart    <- heart    / v_heart
    c_brain    <- brain    / v_brain
    c_muscle   <- muscle   / v_muscle
    c_adipose  <- adipose  / v_adipose
    c_skin     <- skin     / v_skin
    c_spleen   <- spleen   / v_spleen
    c_pancreas <- pancreas / v_pancreas
    c_liver    <- liver    / v_liver
    c_stomach  <- stomach  / v_stomach
    c_gut      <- gut      / v_gut
    c_bone     <- bone     / v_bone
    c_kidney   <- kidney   / v_kidney
    c_other    <- other    / v_other

    # ---- 4. ODE system ----------------------------------------------------
    # Lung is fed by venous blood and drains to arterial blood.
    d/dt(lung)     <- q_lung     * (c_venous   - c_lung     / kp_lung)

    # Organs draining directly to venous blood.
    d/dt(heart)    <- q_heart    * (c_arterial - c_heart    / kp_heart)
    d/dt(brain)    <- q_brain    * (c_arterial - c_brain    / kp_brain)
    d/dt(muscle)   <- q_muscle   * (c_arterial - c_muscle   / kp_muscle)
    d/dt(adipose)  <- q_adipose  * (c_arterial - c_adipose  / kp_adipose)
    d/dt(skin)     <- q_skin     * (c_arterial - c_skin     / kp_skin)
    d/dt(bone)     <- q_bone     * (c_arterial - c_bone     / kp_bone)
    d/dt(kidney)   <- q_kidney   * (c_arterial - c_kidney   / kp_kidney)
    d/dt(other)    <- q_other    * (c_arterial - c_other    / kp_other)

    # Splanchnic organs draining into the liver.
    d/dt(spleen)   <- q_spleen   * (c_arterial - c_spleen   / kp_spleen)
    d/dt(pancreas) <- q_pancreas * (c_arterial - c_pancreas / kp_pancreas)
    d/dt(stomach)  <- q_stomach  * (c_arterial - c_stomach  / kp_stomach)
    d/dt(gut)      <- q_gut      * (c_arterial - c_gut      / kp_gut)

    # Liver: hepatic artery inflow plus portal inflow, less hepatic
    # elimination and venous outflow. The only eliminating organ.
    d/dt(liver) <- q_liver_arterial * c_arterial +
      q_spleen   * c_spleen   / kp_spleen +
      q_pancreas * c_pancreas / kp_pancreas +
      q_stomach  * c_stomach  / kp_stomach +
      q_gut      * c_gut      / kp_gut -
      clint * fub * c_liver / kp_liver -
      q_liver * c_liver / kp_liver

    # Arterial blood is fed only by the lung.
    d/dt(arterial) <- q_lung * (c_lung / kp_lung - c_arterial)

    # Venous blood collects every organ that drains systemically.
    d/dt(venous) <-
      q_heart   * c_heart   / kp_heart +
      q_brain   * c_brain   / kp_brain +
      q_muscle  * c_muscle  / kp_muscle +
      q_adipose * c_adipose / kp_adipose +
      q_skin    * c_skin    / kp_skin +
      q_liver   * c_liver   / kp_liver +
      q_bone    * c_bone    / kp_bone +
      q_kidney  * c_kidney  / kp_kidney +
      q_other   * c_other   / kp_other -
      q_lung    * c_venous

    # ---- 5. Observation ---------------------------------------------------
    # Venous PLASMA concentration in ng/mL: the venous state holds blood, so
    # divide by the blood:plasma ratio, and by 1000 to turn mg/L into ng/mL.
    Cc <- venous / (v_venous * bp / 1000)
    Cc ~ lnorm(expSd)
  })
}
