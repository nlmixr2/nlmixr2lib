Chen_2015_zincNitrate_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse, ICR, male, 6 weeks, 0.031-0.032 kg).",
    "PBPK (whole-body, flow-limited, Berkeley Madonna 8.3.9) model for the",
    "biodistribution of a single intravenous dose of radiolabelled zinc nitrate",
    "(65Zn(NO3)2, i.e. soluble zinc ion) in mice. Nine flow-limited",
    "compartments (blood, lung, gut, spleen, liver, heart, brain, kidney and a",
    "carcass remainder mapped to `other`), with biliary transfer from liver to",
    "gut and faecal and urinary excretion from gut and kidney. Heart, brain and",
    "carcass partitioning rises with time as a four-parameter Hill function and",
    "the renal elimination rate constant falls with time as a Hill function;",
    "the remaining partition coefficients and the biliary and faecal rate",
    "constants are time-invariant. Unlike the two 65ZnO nanoparticle models",
    "from the same paper, this parameter set fitted the whole 28-day",
    "observation period without a late-phase recalibration -- it is in fact the",
    "source of the post-day-7 parameters those models switch to. States hold",
    "CONCENTRATION (ug/mL in blood, ug/g in tissue), as printed in the source,",
    "so the intravenous dose is converted to a blood concentration by",
    "`f(blood) <- 1 / v_blood`. The paper reports no between-subject",
    "variability and no residual error model -- it was fitted by minimising",
    "mean absolute percentage error -- so the model is for typical-value",
    "simulation only. One printed mass-balance anomaly in the liver-to-gut",
    "biliary term is reproduced verbatim; see the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Chen W-Y, Cheng Y-H, Hsieh N-H, Wu B-C, Chou W-C, Ho C-C, Chen J-K,",
    "Liao C-M, Lin P (2015). Physiologically based pharmacokinetic modeling",
    "of zinc oxide nanoparticles and zinc nitrate in mice.",
    "International Journal of Nanomedicine 10:6277-6292.",
    "doi:10.2147/IJN.S86785.",
    "The biodistribution data the model was fitted to were published in",
    "Chen J-K, Ho C-C, Chang H, Lin J-F, Yang CS, Tsai M-H, Tsai H-T, Lin P",
    "(2015) Nanotoxicology 9:43-53 (reference 7 of the present paper);",
    "only the parameter values reproduced in Chen 2015 Tables 2-3 and S1-S4",
    "are used here.",
    sep = " "
  )
  vignette <- "Chen_2015_zincOxideNanoparticles"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")
  dosing <- "blood"

  # `other` carries the paper's "carcass" (muscle and bone) compartment. It is
  # the rest-of-body remainder in the strict sense: the eight named organ
  # volume fractions of Table 2 sum to 0.294 of body weight and the carcass
  # fraction is exactly 1 - 0.294 = 0.706.
  compartmentData <- list(
    blood  = list(analyte = "zinc-65", units = "ug/mL", specimen = "whole blood", verified = TRUE),
    lung   = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    gut    = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    liver  = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    heart  = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    brain  = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    other  = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales cardiac output allometrically (QBl = QC * BW^0.75, Equation 1)",
        "and every compartment volume linearly (Vi = Vi,f * BW, Table 2).",
        "Chen 2015 Table 2 uses BW = 0.032 kg; the dosed animals weighed",
        "0.031-0.032 kg. Tissue density is taken as 1 g/mL so that the Table 2",
        "volume fractions (fractions of body weight) give volumes in mL."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species      = "mouse (ICR, male, 6 weeks old)",
    n_subjects   = NA_integer_,
    n_studies    = 1L,
    age_range    = "6 weeks",
    weight_range = "0.031-0.032 kg",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range   = paste(
      "Single intravenous (tail vein) dose of 120 ug of 65Zn(NO3)2 dissolved",
      "in 400 uL of distilled water."
    ),
    regions      = "Taiwan (National Health Research Institutes, Zhunan)",
    notes        = paste(
      "Biodistribution was measured in blood, liver, lung, kidney, spleen,",
      "brain, heart, gastrointestinal tract and carcass (muscle and bone) at",
      "1, 2, 4 and 7 hours and 1, 2, 3, 7 and 28 days after injection",
      "(Study data section). The number of animals per time point is not",
      "stated in Chen 2015; it is reported in the upstream biodistribution",
      "paper (reference 7), which is not on disk here.",
      "The model was fitted by minimising the mean absolute percentage error",
      "(Equation 6) rather than by maximum likelihood, so no between-subject",
      "variance or residual error is reported. The 95% confidence intervals in",
      "Tables S1-S3 are Monte Carlo uncertainty on the AUC-ratio partition",
      "coefficients, not fitted variance components."
    )
  )

  ini({
    # ==================================================================
    # Physiology -- Chen 2015 Table 2. Nothing here was estimated.
    # ==================================================================
    # Equation 1: QBl = QC * BW^0.75. QC is quoted as L/h for a 1 kg mouse;
    # written here as mL/h so that dividing an amount in ug by a volume in mL
    # gives ug/mL, the unit of Table 1.
    qcc     <- fixed(9025)   ; label("Cardiac output constant QC (mL/h/kg^0.75)")   # Table 2, QC = 9.025 L/h = 9025 mL/h
    e_wt_qco <- fixed(0.75)  ; label("Allometric exponent on cardiac output (unitless)")  # Equation 1, BW^0.75

    # Organ volumes as a fraction of body weight (Table 2). Tissue density is
    # taken as 1 g/mL, so a fraction of BW in kg times 1000 gives mL.
    vc_blood  <- fixed(0.060) ; label("Blood volume as a fraction of body weight (unitless)")     # Table 2, VBl
    vc_lung   <- fixed(0.007) ; label("Lung volume as a fraction of body weight (unitless)")      # Table 2, VLu
    vc_gut    <- fixed(0.127) ; label("GI tract volume as a fraction of body weight (unitless)")  # Table 2, VGI
    vc_spleen <- fixed(0.004) ; label("Spleen volume as a fraction of body weight (unitless)")    # Table 2, VSp
    vc_liver  <- fixed(0.059) ; label("Liver volume as a fraction of body weight (unitless)")     # Table 2, VLi
    vc_heart  <- fixed(0.007) ; label("Heart volume as a fraction of body weight (unitless)")     # Table 2, VHe
    vc_brain  <- fixed(0.014) ; label("Brain volume as a fraction of body weight (unitless)")     # Table 2, VBr
    vc_kidney <- fixed(0.016) ; label("Kidney volume as a fraction of body weight (unitless)")    # Table 2, VKi
    vc_other  <- fixed(0.706) ; label("Carcass volume as a fraction of body weight (unitless)")   # Table 2, VCa

    # Blood flow to organ as a fraction of cardiac output (Table 2). The lung
    # takes the whole cardiac output and, in the Table 1 blood equation, sits in
    # PARALLEL with the other organs rather than in series; the eight fractions
    # therefore sum to 1.995 rather than to 1. This is the authors' structure
    # and is reproduced as printed.
    qc_lung   <- fixed(1)     ; label("Fraction of cardiac output to lung (unitless)")     # Table 2, QLu = 1 (Brown et al)
    qc_gut    <- fixed(0.188) ; label("Fraction of cardiac output to GI tract (unitless)") # Table 2, QGI (Davies and Morris)
    qc_spleen <- fixed(0.011) ; label("Fraction of cardiac output to spleen (unitless)")   # Table 2, QSp (Davies and Morris)
    qc_liver  <- fixed(0.161) ; label("Fraction of cardiac output to liver (unitless)")    # Table 2, QLi (Davies and Morris)
    qc_heart  <- fixed(0.060) ; label("Fraction of cardiac output to heart (unitless)")    # Table 2, QHe (Brown et al)
    qc_brain  <- fixed(0.030) ; label("Fraction of cardiac output to brain (unitless)")    # Table 2, QBr (Brown et al)
    qc_kidney <- fixed(0.091) ; label("Fraction of cardiac output to kidney (unitless)")   # Table 2, QKi (Brown et al)
    qc_other  <- fixed(0.454) ; label("Fraction of cardiac output to carcass (unitless)")  # Table 2, QCa

    # ==================================================================
    # Physicochemical parameters -- Chen 2015 Table 3, 65Zn(NO3)2 column
    # ==================================================================
    # Blank Xmax / TX50 / Xn cells in Table 3 mean the parameter is a CONSTANT
    # equal to the Xini column. Table 3 footnote a (early-phase-only fitting)
    # applies to the 65ZnO NP columns; the 65Zn(NO3)2 simulation "was well
    # fitted with the experimental data" over the whole period and needed no
    # late-phase recalibration (Results, "Simulation and model validation").
    # The time-invariant partition coefficients carry the canonical
    # `lkp_<tissue>` log-scale names and are back-transformed in model().
    lkp_liver  <- fixed(log(2.1))  ; label("Liver:blood partition coefficient (unitless)")     # Table 3, p li, 65Zn(NO3)2, Xini = 2.1 (constant)
    lkp_kidney <- fixed(log(2.1))  ; label("Kidney:blood partition coefficient (unitless)")    # Table 3, p Ki, 65Zn(NO3)2, Xini = 2.1 (constant)
    lkp_spleen <- fixed(log(1.3))  ; label("Spleen:blood partition coefficient (unitless)")    # Table 3, p sp, 65Zn(NO3)2, Xini = 1.3 (constant)
    lkp_lung   <- fixed(log(0.95)) ; label("Lung:blood partition coefficient (unitless)")      # Table 3, p lu, 65Zn(NO3)2, Xini = 0.95 (constant)
    lkp_gut    <- fixed(log(1.7))  ; label("GI tract:blood partition coefficient (unitless)")  # Table 3, p gI, 65Zn(NO3)2, Xini = 1.7 (constant)

    # Heart, brain and carcass partitioning rises sigmoidally with time and is
    # fitted with the four-parameter Hill function of Equation 2:
    #   p(t) = p_ini + (p_max - p_ini) * t^n / (T50^n + t^n)
    # The heart is time-dependent for 65Zn(NO3)2 but not for either 65ZnO NP
    # size (Table 3 Notes).
    kp_heart_ini  <- fixed(0.24) ; label("Heart:blood partition coefficient at time zero (unitless)")            # Table 3, p he, 65Zn(NO3)2, Xini
    kp_heart_max  <- fixed(0.93) ; label("Maximum heart:blood partition coefficient (unitless)")                 # Table 3, p he, 65Zn(NO3)2, Xmax
    kp_heart_t50  <- fixed(3.31) ; label("Time to half-maximum heart:blood partition coefficient (h)")           # Table 3, p he, 65Zn(NO3)2, TX50
    kp_heart_hill <- fixed(2.96) ; label("Hill coefficient of the heart partition-coefficient rise (unitless)")  # Table 3, p he, 65Zn(NO3)2, Xn (r2 = 0.99)

    kp_brain_ini  <- fixed(0.024)  ; label("Brain:blood partition coefficient at time zero (unitless)")            # Table 3, p Br, 65Zn(NO3)2, Xini
    kp_brain_max  <- fixed(0.748)  ; label("Maximum brain:blood partition coefficient (unitless)")                 # Table 3, p Br, 65Zn(NO3)2, Xmax
    kp_brain_t50  <- fixed(55.614) ; label("Time to half-maximum brain:blood partition coefficient (h)")           # Table 3, p Br, 65Zn(NO3)2, TX50
    kp_brain_hill <- fixed(0.929)  ; label("Hill coefficient of the brain partition-coefficient rise (unitless)")  # Table 3, p Br, 65Zn(NO3)2, Xn (r2 = 0.99)

    kp_other_ini  <- fixed(0.085)  ; label("Carcass:blood partition coefficient at time zero (unitless)")            # Table 3, p ca, 65Zn(NO3)2, Xini
    kp_other_max  <- fixed(1.157)  ; label("Maximum carcass:blood partition coefficient (unitless)")                 # Table 3, p ca, 65Zn(NO3)2, Xmax
    kp_other_t50  <- fixed(22.680) ; label("Time to half-maximum carcass:blood partition coefficient (h)")           # Table 3, p ca, 65Zn(NO3)2, TX50
    kp_other_hill <- fixed(0.567)  ; label("Hill coefficient of the carcass partition-coefficient rise (unitless)")  # Table 3, p ca, 65Zn(NO3)2, Xn (r2 = 0.99)

    # Biliary and faecal excretion are time-invariant; renal elimination FALLS
    # with time and uses the DIFFERENT Hill form of Equation 4 (note: not
    # Equation 2):
    #   k(t) = k_min + (k_max - k_min) / (1 + (t / T50)^n)
    # so k(0) = k_max and k(inf) = k_min.
    lkbile <- fixed(log(0.04)) ; label("Biliary excretion rate constant (1/h)")    # Table 3, k li, 65Zn(NO3)2, Xini = 0.04 (constant)
    lkfec  <- fixed(log(0.07)) ; label("GI tract excretion rate constant (1/h)")   # Table 3, k gI, 65Zn(NO3)2, Xini = 0.07 (constant)

    kurine_min  <- fixed(0.02) ; label("Minimum renal elimination rate constant (1/h)")                      # Table 3, k Ki, 65Zn(NO3)2, Xini
    kurine_max  <- fixed(0.07) ; label("Maximum renal elimination rate constant (1/h)")                       # Table 3, k Ki, 65Zn(NO3)2, Xmax
    kurine_t50  <- fixed(16.4) ; label("Time to half-maximum renal elimination rate constant (h)")            # Table 3, k Ki, 65Zn(NO3)2, TX50
    kurine_hill <- fixed(1.90) ; label("Hill coefficient of the renal elimination-rate decline (unitless)")   # Table 3, k Ki, 65Zn(NO3)2, Xn (r2 = 0.99)
  })

  model({
    # ================================================================
    # 1. Physiology -- Chen 2015 Equation 1 and Table 2
    # ================================================================
    # Volumes in mL (density 1 g/mL), flows in mL/h.
    q_co <- qcc * WT^e_wt_qco

    v_blood  <- 1000 * vc_blood  * WT
    v_lung   <- 1000 * vc_lung   * WT
    v_gut    <- 1000 * vc_gut    * WT
    v_spleen <- 1000 * vc_spleen * WT
    v_liver  <- 1000 * vc_liver  * WT
    v_heart  <- 1000 * vc_heart  * WT
    v_brain  <- 1000 * vc_brain  * WT
    v_kidney <- 1000 * vc_kidney * WT
    v_other  <- 1000 * vc_other  * WT

    q_lung   <- qc_lung   * q_co
    q_gut    <- qc_gut    * q_co
    q_spleen <- qc_spleen * q_co
    q_liver  <- qc_liver  * q_co
    q_heart  <- qc_heart  * q_co
    q_brain  <- qc_brain  * q_co
    q_kidney <- qc_kidney * q_co
    q_other  <- qc_other  * q_co

    # ================================================================
    # 2. Partition coefficients and elimination rates
    # ================================================================
    # Time-invariant partition coefficients and the biliary rate constant,
    # back-transformed from their canonical log-scale ini() names.
    kp_liver  <- exp(lkp_liver)
    kp_kidney <- exp(lkp_kidney)
    kp_spleen <- exp(lkp_spleen)
    kp_lung   <- exp(lkp_lung)
    kp_gut    <- exp(lkp_gut)
    kbile     <- exp(lkbile)
    kfeces    <- exp(lkfec)

    # Equation 2 (four-parameter Hill, RISING) for heart, brain and carcass.
    kp_heart <- kp_heart_ini + (kp_heart_max - kp_heart_ini) *
      t^kp_heart_hill / (kp_heart_t50^kp_heart_hill + t^kp_heart_hill)
    kp_brain <- kp_brain_ini + (kp_brain_max - kp_brain_ini) *
      t^kp_brain_hill / (kp_brain_t50^kp_brain_hill + t^kp_brain_hill)
    kp_other <- kp_other_ini + (kp_other_max - kp_other_ini) *
      t^kp_other_hill / (kp_other_t50^kp_other_hill + t^kp_other_hill)

    # Equation 4 (Hill, FALLING -- a different algebraic form from Equation 2).
    kurine <- kurine_min +
      (kurine_max - kurine_min) / (1 + (t / kurine_t50)^kurine_hill)

    # ================================================================
    # 3. ODE system -- Chen 2015 Table 1, equations (A) to (I)
    # ================================================================
    # The states hold CONCENTRATION exactly as Table 1 prints them, so each
    # line below is a term-by-term transcription of its printed equation.
    #
    # Note on equation (A): the lung sits in parallel with the other organs,
    # so QLu appears in both the influx and the efflux sum and the bracketed
    # flow total is 1.995 * QBl. Reproduced as printed.
    d/dt(blood) <-
      (q_lung   * lung   / kp_lung +
       q_spleen * spleen / kp_spleen +
       (q_liver + q_gut) * liver / kp_liver +
       q_heart  * heart  / kp_heart +
       q_brain  * brain  / kp_brain +
       q_kidney * kidney / kp_kidney +
       q_other  * other  / kp_other) / v_blood -
      blood / v_blood *
      (q_lung + q_spleen + q_liver + q_heart +
       q_brain + q_kidney + q_other + q_gut)

    d/dt(lung) <- q_lung / v_lung * (blood - lung / kp_lung)

    # Equation (C). The liver-to-gut biliary term `Cli * kli / pli` is added to
    # the gut CONCENTRATION derivative without rescaling by the volume ratio,
    # while equation (E) removes the same expression from the liver
    # concentration derivative. Because VGI (0.127 BW) differs from VLi
    # (0.059 BW) the transfer is not mass-conserving as printed. Reproduced
    # verbatim; see the vignette Errata.
    d/dt(gut) <- q_gut / v_gut * (blood - gut / kp_gut) +
      liver * kbile / kp_liver -
      gut * kfeces / kp_gut

    d/dt(spleen) <- q_spleen / v_spleen * (blood - spleen / kp_spleen)

    # Equation (E). The excretion term sits OUTSIDE the 1/VLi bracket as
    # printed, so it acts on the liver concentration rather than on its amount.
    d/dt(liver) <-
      (q_liver * blood +
       q_gut * gut / kp_gut -
       (q_liver + q_gut) * liver / kp_liver) / v_liver -
      liver * kbile / kp_liver

    d/dt(heart) <- q_heart / v_heart * (blood - heart / kp_heart)

    d/dt(brain) <- q_brain / v_brain * (blood - brain / kp_brain)

    # Equation (H). As in (E), the elimination term sits outside the QKi/VKi
    # bracket as printed.
    d/dt(kidney) <- q_kidney / v_kidney * (blood - kidney / kp_kidney) -
      kidney * kurine / kp_kidney

    d/dt(other) <- q_other / v_other * (blood - other / kp_other)

    # ================================================================
    # 4. Dosing and observation
    # ================================================================
    # The intravenous dose is given in ug; `blood` holds ug/mL, so the dose is
    # divided by the blood volume to give the initial blood concentration.
    f(blood) <- 1 / v_blood

    # Blood is the sampled matrix in Figure 4; the tissue states are themselves
    # the concentrations plotted in Figure 3 and need no further scaling.
    Cc <- blood
  })
}
