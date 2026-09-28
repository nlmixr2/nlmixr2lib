Zhou_2021_enrofloxacin_pig_pbpk <- function() {
  description <- paste(
    "Veterinary (pig). PBPK (whole-body, flow-limited, acslXtreme 3.0) for",
    "oral enrofloxacin granules and its main metabolite ciprofloxacin in",
    "three-way hybrid pigs, built to predict edible-tissue withdrawal",
    "intervals and a liver-toxicity dose (Zhou et al. 2021, Antibiotics",
    "10:955). The parent (enrofloxacin) sub-model has a stomach and an",
    "intestinal lumen (first-order gastric emptying, first-order",
    "absorption into the liver, first-order faecal loss), venous and",
    "arterial blood pools, an in-series lung that the whole cardiac output",
    "passes through, and perfusion-limited, well-stirred liver, kidney,",
    "muscle, fat and lumped rest-of-body compartments; only unbound",
    "arterial drug perfuses the tissues. Enrofloxacin is cleared by",
    "first-order hepatic metabolism (a fraction of it forming",
    "ciprofloxacin in the liver) and by renal excretion from the kidney",
    "venous concentration. The ciprofloxacin sub-model has the same",
    "blood-lung-tissue structure and is cleared renally. States hold",
    "amounts in umol; the concentration outputs are in mg/L (ug/g). The",
    "model is deterministic; the fixed etas are the Monte Carlo",
    "coefficients of variation the authors assumed for the withdrawal-time",
    "analysis (Table S4), not estimated between-animal variances."
  )
  reference <- paste(
    "Zhou K, Liu A, Ma W, Sun L, Mi K, Xu X, Algharib SA, Xie S, Huang L.",
    "Apply a Physiologically Based Pharmacokinetic Model to Promote the",
    "Development of Enrofloxacin Granules: Predict Withdrawal Interval and",
    "Toxicity Dose. Antibiotics (Basel). 2021;10(8):955.",
    "doi:10.3390/antibiotics10080955.",
    "Model equations transcribed from Supplementary File S2 (acslX code),",
    "which the authors adapted from Lin et al. 2016; parameter values from",
    "Supplementary Tables S7 and S8, Monte Carlo distributions from",
    "Supplementary Table S4.",
    sep = " "
  )
  vignette <- "Zhou_2021_enrofloxacin_pig_pbpk"

  # Dose in mg; f(stomach) converts it to umol, because the acslX code
  # integrates every amount in umol (`DOSEoral = PDOSEoral*BW*MWmol`).
  # Concentrations are reported in mg/L, which the paper treats as ug/g of
  # tissue and ug/mL of plasma.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  # The oral dose goes to the `stomach` state, not `depot` or `central`, so
  # state the route explicitly or the registry records dosing = NA.
  dosing <- "stomach"

  # What each ODE state holds, in what amount units, in what matrix.
  compartmentData <- list(
    stomach = list(analyte = "enrofloxacin", units = "umol", specimen = "administration site", verified = TRUE),
    intestine = list(analyte = "enrofloxacin", units = "umol", specimen = "administration site", verified = TRUE),
    venous = list(analyte = "enrofloxacin", units = "umol", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "enrofloxacin", units = "umol", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    other = list(analyte = "enrofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    venous_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "plasma", verified = TRUE),
    arterial_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "plasma", verified = TRUE),
    lung_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    liver_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    kidney_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    muscle_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    adipose_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    other_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "tissue", verified = TRUE),
    # Cumulative process accumulators, carried so that the two mass-balance
    # identities of the acslX code (`Bal`, `Bal1`) can be checked.
    a_oral_absorbed = list(analyte = "enrofloxacin", units = "umol", specimen = "not applicable", verified = TRUE),
    a_metabolized = list(analyte = "enrofloxacin", units = "umol", specimen = "not applicable", verified = TRUE),
    urine = list(analyte = "enrofloxacin", units = "umol", specimen = "urine", verified = TRUE),
    a_feces = list(analyte = "enrofloxacin", units = "umol", specimen = "faeces", verified = TRUE),
    urine_cipro = list(analyte = "ciprofloxacin", units = "umol", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2021 Table S7 and the File S2 code fix BW = 55 kg, the mean",
        "weight of the study pigs (55 +/- 10 kg, Section 4.2). Organ volumes",
        "and cardiac output are fractions of BW and the renal clearances",
        "are per-kg, so those scale linearly; the hepatic metabolic rate",
        "constant is ALSO multiplied by BW (`Km = KmC*BW`, a 1/h rate",
        "constant), which makes tissue concentrations fall with body weight",
        "at a fixed mg/kg dose. That is why Table S3 reports BW as the most",
        "sensitive parameter (normalised sensitivity coefficients -1.8 to",
        "-2.2). The Monte Carlo analysis drew BW from a normal distribution",
        "with mean 55 kg and SD 10 kg, bounded to 45-65 kg (Table S4)."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "pig (three-way hybrid)",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "55 +/- 10 kg (mean +/- SD; Zhou 2021 Section 4.2)",
    sex_female_pct = NA_real_,
    disease_state = "Clinically healthy pigs",
    dose_range = "Enrofloxacin granules (5% content) 5 mg/kg twice daily in feed for 5 days (10 doses); oral doses up to 600 mg/kg were simulated for the liver-toxicity analysis (Zhou 2021 Figure 5b)",
    regions = "China",
    notes = paste(
      "Residue depletion study (Zhou 2021 Section 4.3): 18 treated pigs and",
      "2 controls; three treated pigs were slaughtered at each of 0.042,",
      "0.5, 1, 2, 3 and 5 days after the last dose and muscle, fat, liver",
      "and kidney were assayed for enrofloxacin and ciprofloxacin by",
      "fluorescence HPLC (LOD = LOQ = 0.02 ug/mL; Table S2). The model was",
      "not fitted statistically. Physiological parameters are pig",
      "literature values (Table S7); the chemical-specific parameters come",
      "from Lin et al. 2016, with the absorption rate constant Ka, the",
      "ciprofloxacin urinary rate constant Kurine1C and the partition",
      "coefficients adjusted by hand to the residue data (Section 4.5).",
      "Plasma predictions were compared against an earlier plasma study by",
      "the same group (Figure S1, ref. 3). The Monte Carlo withdrawal-time",
      "analysis used 1000 virtual pigs (Section 4.7)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Oral absorption (Zhou 2021 Table S8; File S2 code). Kst and Kfeces
    # are inherited from Lin et al. 2016 and held fixed; Ka is one of the
    # three quantities the authors adjusted by hand to the residue data
    # (Section 4.5: 'the visually reasonable values for KurineC1, Ka, and
    # the partition coefficients (PCs) were obtained by an iterative
    # manual adjustment approach').
    # ---------------------------------------------------------------
    lkst <- fixed(log(2))
    label("Gastric emptying rate constant Kst, stomach to intestine (1/h)")                           # Table S8 (Kst = 2.0 /h); File S2 `Kst = 2`
    lka <- log(0.55)
    label("Intestinal absorption rate constant Ka, intestine to liver (1/h)")                         # Table S8 (Ka = 0.55 /h); File S2 `Ka = 0.55`
    lkfec <- fixed(log(0.01))
    label("Faecal elimination rate constant Kfeces, from the intestine (1/h)")                        # Table S8 (Kfeces = 0.01 /h); File S2 `Kfeces = 0.01`

    # ---------------------------------------------------------------
    # Tissue:plasma partition coefficients for enrofloxacin (Table S8;
    # File S2). Log-transformed typical values; the Monte Carlo etas are
    # nevertheless added on the LINEAR scale in model(), because the
    # analysis sampled these parameters from NORMAL distributions
    # (Section 4.7, Table S4).
    # ---------------------------------------------------------------
    lkp_liver <- log(3.2)
    label("Liver:plasma partition coefficient of enrofloxacin PL (unitless)")                         # Table S8 (Pl = 3.2); File S2 `PL = 3.2`
    lkp_kidney <- log(5.5)
    label("Kidney:plasma partition coefficient of enrofloxacin PK (unitless)")                        # Table S8 (Pk = 5.5); File S2 `PK = 5.5`
    lkp_muscle <- log(2.5)
    label("Muscle:plasma partition coefficient of enrofloxacin PM (unitless)")                        # Table S8 (Pm = 2.5); File S2 `PM = 2.5`
    lkp_adipose <- log(0.6)
    label("Fat:plasma partition coefficient of enrofloxacin PF (unitless)")                           # Table S8 (Pf = 0.6); File S2 `PF = 0.6`
    lkp_lung <- log(4.3)
    label("Lung:plasma partition coefficient PLu, shared by enrofloxacin and ciprofloxacin (unitless)")            # Table S8 (Plu = 4.3); File S2 `PLu = 4.3`, which the code also uses in the ciprofloxacin lung equation `CVLu1 = CLu1/PLu`
    lkp_other <- log(8)
    label("Rest-of-body:plasma partition coefficient of enrofloxacin Prest (unitless)")               # Table S8 (Prest = 8.0); File S2 `Prest = 8`

    # Partition coefficients for ciprofloxacin (Table S8 'main metabolite';
    # File S2). Table S8 also lists Plu1 = 4.3, and File S2 declares
    # `PLu1 = 4.3`, but no equation uses it -- the ciprofloxacin lung
    # takes the enrofloxacin PLu (same value). It is therefore not carried.
    lkp_liver_cipro <- log(4.3)
    label("Liver:plasma partition coefficient of ciprofloxacin PL1 (unitless)")                       # Table S8 (PL1 = 4.3); File S2 `PL1 = 4.3`
    lkp_kidney_cipro <- log(5.5)
    label("Kidney:plasma partition coefficient of ciprofloxacin PK1 (unitless)")                      # Table S8 (PK1 = 5.5); File S2 `PK1 = 5.5`
    lkp_muscle_cipro <- log(1.5)
    label("Muscle:plasma partition coefficient of ciprofloxacin PM1 (unitless)")                      # Table S8 (PM1 = 1.5); File S2 `PM1 = 1.5`
    lkp_adipose_cipro <- log(0.53)
    label("Fat:plasma partition coefficient of ciprofloxacin PF1 (unitless)")                         # Table S8 (PF1 = 0.53); File S2 `PF1 = 0.53`
    lkp_other_cipro <- log(8)
    label("Rest-of-body:plasma partition coefficient of ciprofloxacin Prest1 (unitless)")             # Table S8 (Prest1 = 8.0); File S2 `Prest1 = 8`

    # ---------------------------------------------------------------
    # Metabolism, protein binding and excretion (Table S8; File S2).
    # ---------------------------------------------------------------
    kmet <- fixed(0.035)
    label("Hepatic metabolic rate constant KmC per kg body weight (1/(h*kg))")                        # Table S8 (KmC = 0.035 /(h*kg)); File S2 `KmC = 0.035`, used as `Km = KmC*BW`
    fm <- fixed(0.35)
    label("Fraction of metabolised enrofloxacin forming ciprofloxacin Frac (unitless)")               # Table S8 (Frac = 0.35); File S2 `Frac = 0.35`
    # Table S8 reports the BOUND fractions (PB = 0.46, PB1 = 0.19) and the
    # code uses them as `CAfree = CA*(1-PB)`. Encoded on the canonical
    # unbound scale.
    fu <- fixed(0.54)
    label("Fraction of enrofloxacin unbound in plasma (unitless) = 1 - PB")                           # Table S8 (PB = 0.46 bound); File S2 `PB = 0.46`; fu = 1 - 0.46
    fu_cipro <- fixed(0.81)
    label("Fraction of ciprofloxacin unbound in plasma (unitless) = 1 - PB1")                         # Table S8 (PB1 = 0.19 bound); File S2 `PB1 = 0.19`; fu = 1 - 0.19
    lcl_renal <- fixed(log(0.12))
    label("Urinary elimination rate constant of enrofloxacin KurineC (L/h/kg)")                       # Table S8 (KurineC = 0.12 L/h/kg); File S2 `KurineC = 0.12`
    lcl_renal_cipro <- log(0.35)
    label("Urinary elimination rate constant of ciprofloxacin Kurine1C (L/h/kg)")                     # Table S8 (Kurine1C = 0.35 L/h/kg); File S2 `Kurine1C = 0.35`

    # ---------------------------------------------------------------
    # Monte Carlo variability (Zhou 2021 Table S4; Section 4.7). These are
    # NOT estimated variances. The authors sampled 12 chemical-specific
    # parameters (plus BW) from NORMAL distributions with the Table S4
    # mean and SD, bounded to mean - SD and mean + SD, giving 20% CV for
    # partition coefficients and 30% otherwise. Each eta is ADDITIVE on
    # the linear-scale parameter with variance SD^2 from the Table S4 SD
    # column. The names follow the package's eta + transformed-parameter
    # convention (etalkp_liver for lkp_liver) but the eta is NOT on the
    # log scale: 0.4096 is a variance in (unitless)^2, i.e. SD 0.64 = 20%
    # of 3.2. The +/- 1 SD truncation cannot be expressed in an eta; the
    # vignette draws truncated etas explicitly to reproduce the analysis.
    # Table S4 has several internal inconsistencies (see the vignette
    # Errata); the SD column is followed throughout.
    # ---------------------------------------------------------------
    etalkp_liver ~ fixed(0.4096)
    # Table S4 Pl: mean 3.2, SD 0.64 -> 0.64^2, additive on the linear scale
    etalkp_kidney ~ fixed(1.21)
    # Table S4 Pk: mean 5.5, SD 1.1 -> 1.1^2, additive on the linear scale
    etalkp_muscle ~ fixed(0.25)
    # Table S4 Pm: mean 2.5, SD 0.5 -> 0.5^2, additive on the linear scale
    etalkp_adipose ~ fixed(0.0144)
    # Table S4 Pf: mean 0.6, SD 0.12 -> 0.12^2, additive on the linear scale
    etalkp_liver_cipro ~ fixed(0.7396)
    # Table S4 Pl1: mean 4.3, SD 0.86 -> 0.86^2, additive on the linear scale
    etalkp_kidney_cipro ~ fixed(1.21)
    # Table S4 Pk1: mean 5.5, SD 1.1 -> 1.1^2, additive on the linear scale
    etalkp_muscle_cipro ~ fixed(0.25)
    # Table S4 Pm1: mean 1.5, SD 0.5 -> 0.5^2, additive on the linear scale (the printed CV 0.2 would give SD 0.3; SD and bounds 1.0-2.0 both say 0.5)
    etalkp_adipose_cipro ~ fixed(0.011236)
    # Table S4 Pf1: mean 0.53, SD 0.106 -> 0.106^2, additive on the linear scale
    etakmet ~ fixed(0.00011025)
    # Table S4 Kmc: mean 0.035, SD 0.0105 -> 0.0105^2, additive on the linear scale
    etafm ~ fixed(0.011025)
    # Table S4 Frac: mean 0.35, SD 0.105 -> 0.105^2, additive on the linear scale
    etalcl_renal ~ fixed(0.001296)
    # Table S4 Kurinec: SD 0.036 -> 0.036^2, additive on the linear scale (= 0.3 x 0.12, the Table S8 value; Table S4 prints the mean as 0.2)
    etalcl_renal_cipro ~ fixed(0.011025)
    # Table S4 Kurinelc: mean 0.35, SD 0.105 -> 0.105^2, additive on the linear scale

    # ---------------------------------------------------------------
    # Zhou 2021 compared predictions with residue data by linear
    # regression (Figure S2) and reports no residual-error model. Fixed
    # to zero so that simulations return model predictions; not a
    # paper-derived value.
    # ---------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error on plasma enrofloxacin, not reported (fraction)")              # not reported in Zhou 2021
  })

  model({
    # =================================================================
    # Pig physiology (Zhou 2021 Table S7, from Lin et al. 2016; File S2
    # INITIAL block). Literature constants with no Monte Carlo
    # variability, so they are carried as literals.
    # =================================================================
    q_co <- 5 * WT                     # L/h; Table S7 QCC = 5 L/h/kg; File S2 `QC = QCC*BW`
    q_liver <- 0.2725 * q_co           # Table S7 QLC = 0.2725
    q_kidney <- 0.12 * q_co            # Table S7 QKC = 0.12
    q_muscle <- 0.251 * q_co           # Table S7 QMC = 0.251
    q_adipose <- 0.1275 * q_co         # Table S7 QFC = 0.1275
    q_other <- q_co - q_liver - q_kidney - q_muscle - q_adipose
    # File S2 `Qrest = QC-QL-QK-QM-QF`, i.e. 0.229 * QC (Table S7 QrestC = 0.229)

    v_liver <- 0.0247 * WT             # L; Table S7 VLC = 0.0247
    v_kidney <- 0.004 * WT             # Table S7 VKC = 0.004
    v_muscle <- 0.4 * WT               # Table S7 VMC = 0.4
    v_adipose <- 0.32 * WT             # Table S7 VFC = 0.32
    v_lung <- 0.01 * WT                # Table S7 VLuC = 0.01
    v_blood <- 0.06 * WT               # Table S7 VBloodC = 0.06
    v_venous <- 0.74 * v_blood         # File S2 `Vven = VBlood*0.74`
    v_arterial <- 0.26 * v_blood       # File S2 `Vart = VBlood*0.26`
    v_other <- WT - v_liver - v_kidney - v_muscle - v_adipose - v_lung - v_blood
    # File S2 `Vrest = BW-VL-VK-VM-VF-VLu-VBlood`, i.e. 0.1813 * BW (Table S7 VrestC = 0.1813)

    # Molar mass conversions (File S2). The dose is converted mg -> umol
    # with MWmol = 2.78 umol/mg (f(stomach) below) and enrofloxacin
    # amounts back to mg with MWmg = 0.36 mg/umol. File S2 converts
    # ciprofloxacin with `MW1mg` but never assigns it; MW1mg = MW1 / 1000
    # = 331.34 / 1000 mg/umol is the reading that reproduces the paper's
    # 51.9 ug/mL liver ciprofloxacin at the 130 mg/kg Cmax (Section 2.4).
    mw_mg <- 0.36                      # mg/umol; File S2 `MWmg = 0.36`
    mw_mg_cipro <- 331.34 / 1000       # mg/umol; File S2 `MW1 = 331.34` g/mol

    # =================================================================
    # Individual chemical-specific parameters. Monte Carlo etas are
    # additive on the LINEAR scale, exp(typical) + eta, so that each
    # parameter is normally distributed with the Table S4 mean and SD.
    # =================================================================
    kst <- exp(lkst)
    ka <- exp(lka)
    kfec <- exp(lkfec)

    pc_liver <- exp(lkp_liver) + etalkp_liver
    pc_kidney <- exp(lkp_kidney) + etalkp_kidney
    pc_muscle <- exp(lkp_muscle) + etalkp_muscle
    pc_adipose <- exp(lkp_adipose) + etalkp_adipose
    pc_lung <- exp(lkp_lung)
    pc_other <- exp(lkp_other)
    pc_liver_cipro <- exp(lkp_liver_cipro) + etalkp_liver_cipro
    pc_kidney_cipro <- exp(lkp_kidney_cipro) + etalkp_kidney_cipro
    pc_muscle_cipro <- exp(lkp_muscle_cipro) + etalkp_muscle_cipro
    pc_adipose_cipro <- exp(lkp_adipose_cipro) + etalkp_adipose_cipro
    pc_other_cipro <- exp(lkp_other_cipro)

    k_hep <- (kmet + etakmet) * WT                        # 1/h; File S2 `Km = KmC*BW`
    fcip <- fm + etafm
    cl_urine <- (exp(lcl_renal) + etalcl_renal) * WT             # L/h; File S2 `Kurine = KurineC*BW`
    cl_urine_cipro <- (exp(lcl_renal_cipro) + etalcl_renal_cipro) * WT
    # File S2 `Kurine1 = Kurine1C*BW`

    # =================================================================
    # Enrofloxacin concentrations (umol/L). File S2: venous `CV = AV/Vven`;
    # arterial `CAfree = CA*(1-PB)`; each organ `C = A/V`, `CV = C/P`.
    # =================================================================
    c_venous <- venous / v_venous
    c_arterial_free <- fu * arterial / v_arterial
    c_lung <- lung / v_lung
    cv_lung <- c_lung / pc_lung
    c_liver <- liver / v_liver
    cv_liver <- c_liver / pc_liver
    c_kidney <- kidney / v_kidney
    cv_kidney <- c_kidney / pc_kidney
    c_muscle <- muscle / v_muscle
    cv_muscle <- c_muscle / pc_muscle
    c_adipose <- adipose / v_adipose
    cv_adipose <- c_adipose / pc_adipose
    c_other <- other / v_other
    cv_other <- c_other / pc_other

    r_absorb <- ka * intestine        # File S2 `RAO = Ka*AI`
    r_feces <- kfec * intestine       # File S2 `Rfeces = Kfeces*AI`
    r_met <- k_hep * liver            # File S2 `Rmet = Km*CL*VL`
    r_urine <- cl_urine * cv_kidney   # File S2 `Rurine = Kurine*CVK`

    # =================================================================
    # Ciprofloxacin concentrations (umol/L), same structure (File S2
    # ciprofloxacin sub-model). The lung uses the enrofloxacin PLu, as
    # coded (`CVLu1 = CLu1/PLu`).
    # =================================================================
    c_venous_cipro <- venous_cipro / v_venous
    c_arterial_free_cipro <- fu_cipro * arterial_cipro / v_arterial
    c_lung_cipro <- lung_cipro / v_lung
    cv_lung_cipro <- c_lung_cipro / pc_lung
    c_liver_cipro <- liver_cipro / v_liver
    cv_liver_cipro <- c_liver_cipro / pc_liver_cipro
    c_kidney_cipro <- kidney_cipro / v_kidney
    cv_kidney_cipro <- c_kidney_cipro / pc_kidney_cipro
    c_muscle_cipro <- muscle_cipro / v_muscle
    cv_muscle_cipro <- c_muscle_cipro / pc_muscle_cipro
    c_adipose_cipro <- adipose_cipro / v_adipose
    cv_adipose_cipro <- c_adipose_cipro / pc_adipose_cipro
    c_other_cipro <- other_cipro / v_other
    cv_other_cipro <- c_other_cipro / pc_other_cipro

    r_urine_cipro <- cl_urine_cipro * cv_kidney_cipro   # File S2 `Rurine1 = Kurine1*CVK1`

    # =================================================================
    # Enrofloxacin ODEs (File S2 DERIVATIVE block). The absorbed drug
    # enters the liver directly (`RL = QL*(CAfree-CVL)+RAO-Rmet`).
    # =================================================================
    d/dt(stomach) <- -kst * stomach                                   # File S2 `RAST = RDOSEoral-Kst*AST`
    d/dt(intestine) <- kst * stomach - r_absorb - r_feces             # File S2 `RAI = Kst*AST-Ka*AI-Kfeces*AI`
    d/dt(venous) <- q_liver * cv_liver + q_kidney * cv_kidney + q_muscle * cv_muscle +
      q_adipose * cv_adipose + q_other * cv_other - q_co * c_venous   # File S2 `RV` (IV/IM/SC inputs are zero in this study)
    d/dt(arterial) <- q_co * cv_lung - q_co * c_arterial_free        # File S2 `RA = QC*CVLu-QC*CAfree`
    d/dt(lung) <- q_co * (c_venous - cv_lung)                         # File S2 `RALu = QC*(CV-CVLu)`
    d/dt(liver) <- q_liver * (c_arterial_free - cv_liver) + r_absorb - r_met   # File S2 `RL`
    d/dt(kidney) <- q_kidney * (c_arterial_free - cv_kidney) - r_urine          # File S2 `RK`
    d/dt(muscle) <- q_muscle * (c_arterial_free - cv_muscle)                    # File S2 `RM`
    d/dt(adipose) <- q_adipose * (c_arterial_free - cv_adipose)                 # File S2 `RF`
    d/dt(other) <- q_other * (c_arterial_free - cv_other)                       # File S2 `Rrest`

    # =================================================================
    # Ciprofloxacin ODEs (File S2). Formed in the liver at
    # `Rmet1 = Rmet*Frac` (molar).
    # =================================================================
    d/dt(venous_cipro) <- q_liver * cv_liver_cipro + q_kidney * cv_kidney_cipro +
      q_muscle * cv_muscle_cipro + q_adipose * cv_adipose_cipro +
      q_other * cv_other_cipro - q_co * c_venous_cipro                           # File S2 `RV1`
    d/dt(arterial_cipro) <- q_co * cv_lung_cipro - q_co * c_arterial_free_cipro  # File S2 `RA1`
    d/dt(lung_cipro) <- q_co * (c_venous_cipro - cv_lung_cipro)                  # File S2 `RALu1`
    d/dt(liver_cipro) <- q_liver * (c_arterial_free_cipro - cv_liver_cipro) + fcip * r_met   # File S2 `RL1`
    d/dt(kidney_cipro) <- q_kidney * (c_arterial_free_cipro - cv_kidney_cipro) - r_urine_cipro   # File S2 `RK1`
    d/dt(muscle_cipro) <- q_muscle * (c_arterial_free_cipro - cv_muscle_cipro)   # File S2 `RM1`
    d/dt(adipose_cipro) <- q_adipose * (c_arterial_free_cipro - cv_adipose_cipro)   # File S2 `RF1`
    d/dt(other_cipro) <- q_other * (c_arterial_free_cipro - cv_other_cipro)      # File S2 `Rrest1`

    # Cumulative accumulators for the File S2 mass balances
    # `Bal = AAO - Tmass` and `Bal1 = Amet1 - Tmass1`.
    d/dt(a_oral_absorbed) <- r_absorb      # File S2 `AAO = Integ(RAO,0.0)`
    d/dt(a_metabolized) <- r_met           # File S2 `Amet = Integ(Rmet,0.0)`
    d/dt(urine) <- r_urine                 # File S2 `Aurine = Integ(Rurine,0.0)`
    d/dt(a_feces) <- r_feces               # File S2 `Afeces = Integ(Rfeces,0.0)`
    d/dt(urine_cipro) <- r_urine_cipro     # File S2 `Aurine1 = Integ(Rurine1,0.0)`

    # Oral dose in mg -> umol (File S2 `DOSEoral = PDOSEoral*BW*MWmol`,
    # `MWmol = 2.78` umol/mg). No bioavailability term: incomplete
    # absorption arises from the competing faecal loss.
    f(stomach) <- 2.78

    # =================================================================
    # Outputs (mg/L = ug/mL = ug/g). File S2 `CVmg`, `CLmg`, `CKmg`,
    # `CMmg`, `CFmg`; the `_cipro` outputs are `CL1mg`, ... and the
    # `_total` outputs are the enrofloxacin + ciprofloxacin residue
    # marker (`CLtotalmg`, ...), which is compared with the maximum
    # residue limits for the withdrawal time.
    # =================================================================
    Cc <- c_venous * mw_mg
    Cliver <- c_liver * mw_mg
    Ckidney <- c_kidney * mw_mg
    Cmuscle <- c_muscle * mw_mg
    Cadipose <- c_adipose * mw_mg
    Clung <- c_lung * mw_mg
    Cc_cipro <- c_venous_cipro * mw_mg_cipro
    Cliver_cipro <- c_liver_cipro * mw_mg_cipro
    Ckidney_cipro <- c_kidney_cipro * mw_mg_cipro
    Cmuscle_cipro <- c_muscle_cipro * mw_mg_cipro
    Cadipose_cipro <- c_adipose_cipro * mw_mg_cipro
    Cliver_total <- Cliver + Cliver_cipro
    Ckidney_total <- Ckidney + Ckidney_cipro
    Cmuscle_total <- Cmuscle + Cmuscle_cipro
    Cadipose_total <- Cadipose + Cadipose_cipro

    Cc ~ prop(propSd)
  })
}
