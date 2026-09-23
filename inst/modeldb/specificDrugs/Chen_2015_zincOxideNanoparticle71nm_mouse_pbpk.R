Chen_2015_zincOxideNanoparticle71nm_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse, ICR, male, 6 weeks, 0.031-0.032 kg).",
    "PBPK (whole-body, flow-limited, Berkeley Madonna 8.3.9) model for the",
    "biodistribution of a single intravenous dose of 71 nm radiolabelled zinc",
    "oxide nanoparticles (65ZnO NPs) in mice. Nine flow-limited compartments",
    "(blood, lung, gut, spleen, liver, heart, brain, kidney and a carcass",
    "remainder mapped to `other`), with biliary transfer from liver to gut and",
    "faecal and urinary excretion from gut and kidney. Brain and carcass",
    "partitioning rises with time as a four-parameter Hill function; the spleen",
    "partition coefficient is piecewise constant (it drops sharply over days",
    "3-7, which is what the authors optimised it to capture); the three",
    "excretion/elimination rate constants are constant. The authors' calibrated",
    "model additionally assumes the nanoparticles have largely decomposed to",
    "zinc ion by day 7 and switches every partition coefficient and excretion",
    "rate to a second, constant parameter set at that point (Table 4); the",
    "switch time is exposed as `tdecomp` so that setting it beyond the",
    "simulation horizon reproduces the paper's uncalibrated simulation. States",
    "hold CONCENTRATION (ug/mL in blood, ug/g in tissue), as printed in the",
    "source, so the intravenous dose is converted to a blood concentration by",
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
    "only the parameter values reproduced in Chen 2015 Tables 2-4 and S1-S4",
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
    blood = list(analyte = "zinc-65", units = "ug/mL", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE),
    other = list(analyte = "zinc-65", units = "ug/g", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales cardiac output allometrically (QBl = QC * BW^0.75, Equation 1)",
        "and every compartment volume linearly (Vi = Vi,f * BW, Table 2).",
        "Chen 2015 Table 2 uses BW = 0.032 kg; the dosed animals weighed",
        "0.031-0.032 kg. Tissue density is taken as 1 g/mL so that the Table 2",
        "volume fractions (fractions of body weight) give volumes in mL."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "mouse (ICR, male, 6 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "6 weeks",
    weight_range = "0.031-0.032 kg",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = paste(
      "Single intravenous (tail vein) dose of 120 ug of suspended 65ZnO",
      "nanoparticles in 400 uL of distilled water."
    ),
    regions = "Taiwan (National Health Research Institutes, Zhunan)",
    notes = paste(
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
    # Early phase (up to day 7) -- Chen 2015 Table 3, 71 nm 65ZnO NP column
    # ==================================================================
    # Table 3 footnote a: the 65ZnO NP physicochemical parameters are fitted for
    # use on the early phase (first 7 days) only. Blank Xmax / TX50 / Xn cells
    # in Table 3 mean the parameter is a CONSTANT equal to the Xini column.
    kp_liver_early  <- fixed(6.28) ; label("Liver:blood partition coefficient, days 0-7 (unitless)")     # Table 3, p li, 71 nm, Xini = 6.28 (constant)
    kp_kidney_early <- fixed(3.52) ; label("Kidney:blood partition coefficient, days 0-7 (unitless)")    # Table 3, p Ki, 71 nm, Xini = 3.52 (constant)
    kp_lung_early   <- fixed(3.78) ; label("Lung:blood partition coefficient, days 0-7 (unitless)")      # Table 3, p lu, 71 nm, Xini = 3.78 (constant)
    kp_heart_early  <- fixed(1.22) ; label("Heart:blood partition coefficient, days 0-7 (unitless)")     # Table 3, p he, 71 nm, Xini = 1.22 (constant)
    kp_gut_early    <- fixed(2.70) ; label("GI tract:blood partition coefficient, days 0-7 (unitless)")  # Table 3, p gI, 71 nm, Xini = 2.70 (constant)

    # The spleen partition coefficient is piecewise constant within the early
    # phase (Table 3 footnotes c and d). The authors state in Results
    # ("Simulation and model validation") that these are OPTIMIZED values,
    # chosen because the observed 71 nm spleen concentration "decreased
    # substantially to almost zero from days 3 to 7".
    kp_spleen_early <- fixed(3.7)  ; label("Spleen:blood partition coefficient, days 0-3 (unitless)")  # Table 3, p sp, 71 nm, first value 3.7 (footnote c, first 3 days)
    kp_spleen_mid   <- fixed(1.44) ; label("Spleen:blood partition coefficient, days 3-7 (unitless)")  # Table 3, p sp, 71 nm, second value 1.44 (footnote d, days 3-7)

    # Brain and carcass partitioning rises sigmoidally with time and is fitted
    # with the four-parameter Hill function of Equation 2:
    #   p(t) = p_ini + (p_max - p_ini) * t^n / (T50^n + t^n)
    kp_brain_ini  <- fixed(0.002)  ; label("Brain:blood partition coefficient at time zero (unitless)")            # Table 3, p Br, 71 nm, Xini
    kp_brain_max  <- fixed(1.969)  ; label("Maximum brain:blood partition coefficient (unitless)")                 # Table 3, p Br, 71 nm, Xmax
    kp_brain_t50  <- fixed(61.970) ; label("Time to half-maximum brain:blood partition coefficient (h)")           # Table 3, p Br, 71 nm, TX50
    kp_brain_hill <- fixed(0.711)  ; label("Hill coefficient of the brain partition-coefficient rise (unitless)")  # Table 3, p Br, 71 nm, Xn (r2 = 0.98)

    kp_other_ini  <- fixed(0.0002)  ; label("Carcass:blood partition coefficient at time zero (unitless)")           # Table 3, p ca, 71 nm, Xini
    kp_other_max  <- fixed(5.463)   ; label("Maximum carcass:blood partition coefficient (unitless)")                # Table 3, p ca, 71 nm, Xmax
    kp_other_t50  <- fixed(103.700) ; label("Time to half-maximum carcass:blood partition coefficient (h)")          # Table 3, p ca, 71 nm, TX50
    kp_other_hill <- fixed(0.483)   ; label("Hill coefficient of the carcass partition-coefficient rise (unitless)") # Table 3, p ca, 71 nm, Xn (r2 = 0.98)

    # For 71 nm particles all three excretion/elimination rate constants are
    # time-INDEPENDENT (blank Xmax / TX50 / Xn cells in Table 3).
    kbile_early  <- fixed(0.023) ; label("Biliary excretion rate constant, days 0-7 (1/h)")    # Table 3, k li, 71 nm, Xini = 0.023 (constant)
    kfeces_early <- fixed(0.079) ; label("GI tract excretion rate constant, days 0-7 (1/h)")   # Table 3, k gI, 71 nm, Xini = 0.079 (constant)
    kurine_early <- fixed(0.128) ; label("Renal elimination rate constant, days 0-7 (1/h)")    # Table 3, k Ki, 71 nm, Xini = 0.128 (constant)

    # ==================================================================
    # Late phase (after day 7) -- Chen 2015 Table 4, 71 nm column
    # ==================================================================
    # Results, "Simulation and model validation": the uncalibrated simulation
    # barely fit the days 7-28 data, so the authors replaced every partition
    # coefficient and excretion rate with a constant set derived from the
    # 65Zn(NO3)2 data, on the hypothesis that the nanoparticles have decomposed
    # to zinc ion by day 7. This is the paper's CALIBRATED (final) model.
    tdecomp <- fixed(168) ; label("Time at which the nanoparticles are assumed decomposed to zinc ion (h)")  # Results / Discussion: "after day 7" = 168 h

    kp_liver_post  <- fixed(3.53) ; label("Liver:blood partition coefficient after day 7 (unitless)")     # Table 4, p li, 71 nm
    kp_kidney_post <- fixed(0.91) ; label("Kidney:blood partition coefficient after day 7 (unitless)")    # Table 4, p Ki, 71 nm
    kp_spleen_post <- fixed(0.26) ; label("Spleen:blood partition coefficient after day 7 (unitless)")    # Table 4, p sp, 71 nm
    kp_lung_post   <- fixed(0.34) ; label("Lung:blood partition coefficient after day 7 (unitless)")      # Table 4, p lu, 71 nm
    kp_heart_post  <- fixed(0.20) ; label("Heart:blood partition coefficient after day 7 (unitless)")     # Table 4, p he, 71 nm
    kp_gut_post    <- fixed(1.63) ; label("GI tract:blood partition coefficient after day 7 (unitless)")  # Table 4, p gI, 71 nm
    kp_brain_post  <- fixed(0.93) ; label("Brain:blood partition coefficient after day 7 (unitless)")     # Table 4, p Br, 71 nm
    kp_other_post  <- fixed(4.45) ; label("Carcass:blood partition coefficient after day 7 (unitless)")   # Table 4, p ca, 71 nm

    kbile_post  <- fixed(0.0229) ; label("Biliary excretion rate constant after day 7 (1/h)")    # Table 4, k li, 71 nm
    kfeces_post <- fixed(0.0788) ; label("GI tract excretion rate constant after day 7 (1/h)")   # Table 4, k gI, 71 nm
    kurine_post <- fixed(0.1277) ; label("Renal elimination rate constant after day 7 (1/h)")    # Table 4, k Ki, 71 nm
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
    # 2. Time-dependent partition coefficients
    # ================================================================
    # Phase indicator: 1 before the assumed decomposition time, 0 after.
    early <- (t <= tdecomp)

    # Equation 2 (four-parameter Hill, RISING) for brain and carcass.
    kp_brain_t <- kp_brain_ini + (kp_brain_max - kp_brain_ini) *
      t^kp_brain_hill / (kp_brain_t50^kp_brain_hill + t^kp_brain_hill)
    kp_other_t <- kp_other_ini + (kp_other_max - kp_other_ini) *
      t^kp_other_hill / (kp_other_t50^kp_other_hill + t^kp_other_hill)

    # Table 3 footnotes c and d: the spleen is piecewise constant within the
    # early phase, switching at day 3 (72 h).
    kp_spleen_t <- kp_spleen_early * (t <= 72) + kp_spleen_mid * (t > 72)

    # Switch to the Table 4 post-decomposition constants after day 7.
    kp_liver  <- early * kp_liver_early  + (1 - early) * kp_liver_post
    kp_kidney <- early * kp_kidney_early + (1 - early) * kp_kidney_post
    kp_lung   <- early * kp_lung_early   + (1 - early) * kp_lung_post
    kp_heart  <- early * kp_heart_early  + (1 - early) * kp_heart_post
    kp_gut    <- early * kp_gut_early    + (1 - early) * kp_gut_post
    kp_spleen <- early * kp_spleen_t     + (1 - early) * kp_spleen_post
    kp_brain  <- early * kp_brain_t      + (1 - early) * kp_brain_post
    kp_other  <- early * kp_other_t      + (1 - early) * kp_other_post

    kbile  <- early * kbile_early  + (1 - early) * kbile_post
    kfeces <- early * kfeces_early + (1 - early) * kfeces_post
    kurine <- early * kurine_early + (1 - early) * kurine_post

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
