Sharma_2023_nitrofurantoin_rat_pbpk <- function() {
  description <- "Preclinical (rat). PBPK (whole-body, GNU MCSim 6.1, Bayesian MCMC). Nitrofurantoin (NFT) disposition in rats -- Sharma 2023 'Model V5a', the paper's final structure. Structurally IDENTICAL to the rabbit sibling (see modellib('Sharma_2023_nitrofurantoin_rabbit_pbpk')): thirteen mass-balance amount states with perfusion-rate-limited distribution, non-linear renal handling (glomerular filtration + saturable active tubular secretion - first-order tubular reabsorption), Michaelis-Menten hepatic metabolism, and saturable hepatobiliary efflux feeding an enterohepatic recirculation loop. What differs is (a) rat physiology and (b) the rat arm of the hierarchical enterohepatic-recirculation fit. The renal and metabolic parameters were NOT refitted to rat data -- they are the rabbit model V4 estimates carried over unchanged and rescaled by rat body weight, which is what makes the rat plasma and urine predictions a genuine cross-species test rather than a fit (Results 3.3: 'the fact that the model was extrapolated from rabbits by adapting to known rat physiology and that biochemical parameters were scaled just based on rat body weight provides strong confidence'). The rat is the species that supplied the biliary-excretion data, so its hepatobiliary efflux capacity is ~24-fold higher than the rabbit's and enterohepatic recirculation is a real feature of the rat plasma profile (visible late-phase biphasic kinetics). The rat EHR estimates are also the ones the human model inherits by allometric scaling. Deterministic typical-value simulator: the MCMC posterior describes PARAMETER UNCERTAINTY, not between-subject variability, so no IIV is encoded (see the vignette Errata)."
  reference <- paste(
    "Sharma RP, Burgers EJ, Beltman JB. (2023).",
    "Development of a Physiologically Based Pharmacokinetic Model for",
    "Nitrofurantoin in Rabbits, Rats, and Humans.",
    "Pharmaceutics 15(9):2199. doi:10.3390/pharmaceutics15092199. PMCID PMC10535763.",
    "Model equations: Supplementary Materials 'Standard ordinary differential equations",
    "used in PBPK model for NFT'. Biochemical parameters: Supplementary Table S4.",
    "Rat physiology: Supplementary Table S1 is internally corrupt (its organ-volume block",
    "is byte-identical to the human Table S3 and the fractional gut volume is missing",
    "entirely), so rat physiology is taken from the author-deposited model file, which is",
    "the artifact that generated every published rat figure -- see the vignette Errata.",
    "Author-deposited GNU MCSim source and MCMC posterior chains:",
    "doi:10.5281/zenodo.8276305 (file ExtrapolatedRatehrRabbitRenal_mixed.model.R,",
    "identified by the deposit README as ModelV5a).",
    "The article carries a publisher Correction Statement (republished with a minor change",
    "not affecting scientific content); no parameter value is revised by it.",
    sep = " "
  )
  vignette <- "Sharma_2023_nitrofurantoin"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L",
    amount        = "mg",
    weight        = "kg"
  )

  compartmentData <- list(
    a_gut_lumen      = list(analyte = "nitrofurantoin", units = "mg", specimen = "administration site", verified = TRUE),
    a_gut            = list(analyte = "nitrofurantoin", units = "mg", specimen = "tissue",              verified = TRUE),
    a_liver          = list(analyte = "nitrofurantoin", units = "mg", specimen = "tissue",              verified = TRUE),
    a_hepatic        = list(analyte = "nitrofurantoin metabolites (lumped)", units = "mg", specimen = "tissue", verified = TRUE),
    a_bile           = list(analyte = "nitrofurantoin", units = "mg", specimen = "bile",                verified = TRUE),
    a_feces          = list(analyte = "nitrofurantoin", units = "mg", specimen = "faeces",              verified = TRUE),
    a_kidney         = list(analyte = "nitrofurantoin", units = "mg", specimen = "tissue",              verified = TRUE),
    a_filtrate       = list(analyte = "nitrofurantoin", units = "mg", specimen = "urine",               verified = TRUE),
    a_urine_storage  = list(analyte = "nitrofurantoin", units = "mg", specimen = "urine",               verified = TRUE),
    a_fat            = list(analyte = "nitrofurantoin", units = "mg", specimen = "tissue",              verified = TRUE),
    a_rest_of_body   = list(analyte = "nitrofurantoin", units = "mg", specimen = "tissue",              verified = TRUE),
    a_plasma         = list(analyte = "nitrofurantoin", units = "mg", specimen = "plasma",              verified = TRUE),
    a_urine          = list(analyte = "nitrofurantoin", units = "mg", specimen = "urine",               verified = TRUE)
  )

  # No covariates -- fixed 0.25 kg reference rat. Override BW via
  # rxSolve(params = c(BW = ...)) to rescale volumes, flows and BW-scaled
  # rate parameters consistently.
  covariateData <- list()

  population <- list(
    species        = "rat",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = "0.25 kg reference rat (author-deposited rat model file)",
    weight_median  = "0.25 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy rats (kinetic studies, no disease model)",
    dose_range     = "3 to 25 mg/kg single IV (Figure S2); oral dosing also simulated (Figure S6)",
    notes          = paste(
      "Calibration and test data are DIGITISED literature values, not an individual-level",
      "dataset: Sharma 2023 Methods 2.3 states rat kinetics data upon IV or oral dosing were",
      "taken from references 23 and 28 using WebPlotDigitizer, and the fitted points are the",
      "MEAN of the experimental data at each time point. No subject count, age range, strain",
      "or sex split is recoverable from the source paper, and the strain is not stated --",
      "hence species is recorded as plain 'rat' rather than a named strain.",
      "Reference 23 is the study that reported biliary NFT excretion and identified BCRP /",
      "ABCG2 involvement; it supplied the single biliary-excretion dose that anchors the",
      "enterohepatic-recirculation parameters (Results 3.3 notes this limited biliary data",
      "availability -- only one dose -- as a source of uncertainty)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Rat physiology -- author-deposited ExtrapolatedRatehrRabbitRenal_mixed.model.R
    # (ModelV5a). Supplementary Table S1 is NOT used for the organ-volume block:
    # Table S1's Fliver 0.026 / Fkidney 0.0073 / Ffiltrate 0.00073 / Ffat 0.187 /
    # Fplasma 0.0428 are byte-identical to the HUMAN Table S3 values, and Table S1
    # has no fractional-gut-volume row at all even though the rat gut ODE requires
    # one (c_gut = a_gut / v_gut). Table S1 agrees with the deposited file on
    # FQliver, FQkidney, FQfat and Qfiltrate. See vignette Errata.
    # ------------------------------------------------------------------
    BW         <- fixed(0.25)     ; label("Body weight (kg)")                                          # deposited rat file BW = 0.25; Table S1 'Body weight' 0.25 (agrees)
    QCC        <- fixed(15.5)     ; label("Cardiac output coefficient (L/h/kg^0.75)")                   # deposited rat file QCC = 15.5 (Table S1 prints 15.7)
    FQliver    <- fixed(0.174)    ; label("Fraction of cardiac output to liver (unitless)")             # deposited rat file; Table S1 0.174 (agrees; Campbell 2016)
    FQkidney   <- fixed(0.141)    ; label("Fraction of cardiac output to kidney (unitless)")            # deposited rat file; Table S1 0.141 (agrees; Brown 1997)
    FQfat      <- fixed(0.07)     ; label("Fraction of cardiac output to fat (unitless)")               # deposited rat file; Table S1 0.07 (agrees; Campbell 2016)
    FQgut      <- fixed(0.075)    ; label("Fraction of cardiac output to gut (unitless)")               # deposited rat file FQgut = 0.075 (Table S1 prints 0.021 -- not used, see Errata)
    Fliver     <- fixed(0.036)    ; label("Liver volume as a fraction of body weight (unitless)")       # deposited rat file (Table S1 prints the human 0.026)
    Fkidney    <- fixed(0.0073)   ; label("Kidney volume as a fraction of body weight (unitless)")      # deposited rat file; Table S1 0.0073 (agrees)
    Ffiltrate  <- fixed(0.00073)  ; label("Tubular filtrate volume as a fraction of body weight (unitless)")  # deposited rat file; Table S1 0.00073, 10% of kidney volume (agrees)
    Ffat       <- fixed(0.07)     ; label("Fat volume as a fraction of body weight (unitless)")         # deposited rat file (Table S1 prints the human 0.187)
    Fgut       <- fixed(0.027)    ; label("Gut volume as a fraction of body weight (unitless)")         # deposited rat file (Table S1 has no fractional-gut-volume row)
    Fplasma    <- fixed(0.074)    ; label("Plasma volume as a fraction of body weight (unitless)")      # deposited rat file (Table S1 prints the human 0.0428)
    Qfiltrate  <- fixed(0.129)    ; label("Glomerular filtration rate (L/h)")                           # deposited rat file; Table S1 0.129 (agrees; Katayama 2010)

    # ------------------------------------------------------------------
    # Tissue:plasma partition coefficients -- Table S4, shared across all three
    # species. Calculated (Berezhkovskiy) -> fixed(); the rest-of-body
    # coefficient is the one fitted partition coefficient.
    # ------------------------------------------------------------------
    Kgut_plasma      <- fixed(0.622) ; label("Gut:plasma partition coefficient (unitless)")             # Table S4 Kgut_plasma, Calculated
    Kliver_plasma    <- fixed(0.651) ; label("Liver:plasma partition coefficient (unitless)")           # Table S4 Kliver_plasma, Calculated
    Kkidney_plasma   <- fixed(0.671) ; label("Kidney:plasma partition coefficient (unitless)")          # Table S4 Kkidney_plasma, Calculated
    Kfat_plasma      <- fixed(0.159) ; label("Fat:plasma partition coefficient (unitless)")             # Table S4 Kfat_plasma, Calculated
    Krestbody_plasma <- 0.423        ; label("Rest-of-body:plasma partition coefficient (unitless)")    # Table S4, Fitted, 0.423 (0.39-0.45); V4 posterior mean 0.4213 (deposited rat file rounds to 0.421)

    fu  <- fixed(0.4149) ; label("Fraction unbound in plasma (unitless)")                               # deposited rat file fu = 0.4149; Table S4 rounds to 0.42 (Watari 1985)
    Kbp <- fixed(0.76)   ; label("Blood:plasma partition coefficient (unitless)")                       # Table S4 Kbp = 0.76 (Zhang 2022)

    # ------------------------------------------------------------------
    # Renal and hepatic-metabolism parameters -- Table S4, Fitted on RABBIT IV
    # data (model V4) and carried over to the rat UNCHANGED, rescaled only by
    # rat body weight (Results 3.3). See the rabbit sibling for the QurineC
    # 11.45-vs-13.45 Table S4 typo resolution; the same value applies here.
    # ------------------------------------------------------------------
    QurineC <- 13.45 ; label("Urine excretion rate coefficient (1/h/kg)")                               # deposited files + V4 posterior mean 13.456; Table S4 typo 11.45
    Trc     <- 1.33  ; label("Tubular reabsorption rate coefficient (L/h/kg)")                          # Table S4 krc, Fitted, 1.33 (1.11-1.56); V4 posterior mean 1.3349
    Tmc     <- 8.02  ; label("Maximum active tubular secretion rate coefficient (mg/h/kg)")             # Table S4 Tmc, Fitted, 8.02 (6.78-9.30); V4 posterior mean 8.0241
    Kt      <- 0.059 ; label("Concentration at half-maximal active tubular secretion (mg/L)")           # Table S4 Kt, Fitted, 0.059 (0.043-0.079); V4 posterior mean 0.0593
    VmaxC   <- 0.47  ; label("Maximum hepatic metabolism rate coefficient (mg/h/g liver)")              # Table S4 VmaxC, Fitted, 0.47 (0.42-0.53); V4 posterior mean 0.4729
    km      <- 5.83  ; label("Concentration at half-maximal hepatic metabolism (mg/L)")                 # Table S4 Km, Fitted, 5.83 (4.96-6.82); V4 posterior mean 5.8324

    # ------------------------------------------------------------------
    # Enterohepatic recirculation -- Table S4 RAT column, the rat arm of the
    # hierarchical rat+rabbit model V5a fit (Results 3.3). These are the values
    # the HUMAN model inherits by allometric scaling from rat body weight
    # (Table S4 note: parameters marked ** use rat body weight as the basis).
    # ------------------------------------------------------------------
    Vehrc   <- 0.52   ; label("Maximum hepatobiliary excretion rate coefficient (mg/h/kg)")             # Table S4 Vehrc rat column, Fitted, 0.52 (0.43-0.66)
    Kehr    <- 0.017  ; label("Concentration at half-maximal hepatobiliary excretion (mg/L)")           # Table S4 Kehr rat column, Fitted, 0.017 (0.0014-0.063)
    kbile   <- 0.256  ; label("Bile-to-gut-lumen transfer rate coefficient (1/h/kg)")                   # Table S4 kbile rat column, Fitted, 0.256 (0.007-0.83)
    kgutabs <- 2.11   ; label("Gut-lumen absorption rate coefficient (1/h/kg)")                         # Table S4 kgutabs rat column, Fitted, 2.11 (1.80-2.48)
    kfeces  <- 0.0187 ; label("Faecal excretion rate constant from gut lumen (1/h)")                    # Table S4 kfeces rat column, Fitted, 0.0187 (0.00063-0.0649)

    # ------------------------------------------------------------------
    # Residual error -- Methods 2.3, assumed 10% CV -> fixed(). No IIV; see the
    # rabbit sibling and the vignette Errata.
    # ------------------------------------------------------------------
    propSd <- fixed(0.10) ; label("Proportional residual error (fraction)")                             # Methods 2.3 'coefficient of variation of 10%'
  })

  model({
    # ================================================================
    # 1. Cardiac output and organ plasma flows (L/h). The haematocrit
    #    correction is commented out in every deposited model file
    #    (QCplasma = QCblood), so HCT is not used and is not declared.
    # ================================================================
    qc_plasma <- QCC * BW^0.75
    q_liver   <- FQliver  * qc_plasma
    q_gut     <- FQgut    * qc_plasma
    q_kidney  <- FQkidney * qc_plasma
    q_fat     <- FQfat    * qc_plasma
    q_rest    <- qc_plasma - (q_liver + q_kidney + q_fat + q_gut)

    # ================================================================
    # 2. Organ volumes (L); rest-of-body closes total body volume (0.84 * BW)
    #    by subtraction.
    # ================================================================
    v_gut      <- Fgut      * BW
    v_liver    <- Fliver    * BW
    v_kidney   <- Fkidney   * BW
    v_filtrate <- Ffiltrate * BW
    v_fat      <- Ffat      * BW
    v_plasma   <- Fplasma   * BW
    v_rest     <- 0.84 * BW - v_liver - v_kidney - v_fat - v_plasma - v_gut - v_filtrate

    # ================================================================
    # 3. Body-weight scaling of the fitted rate parameters.
    # ================================================================
    qurine    <- QurineC * BW
    tr        <- Trc     * BW
    tm        <- Tmc     * BW
    vmax      <- VmaxC   * 1000 * v_liver
    vehr      <- Vehrc   * BW
    kbile_s   <- kbile   * BW
    kgutabs_s <- kgutabs * BW

    # ================================================================
    # 4. Tissue concentrations (mg/L)
    # ================================================================
    c_gut      <- a_gut          / v_gut
    c_liver    <- a_liver        / v_liver
    c_kidney   <- a_kidney       / v_kidney
    c_filtrate <- a_filtrate     / v_filtrate
    c_fat      <- a_fat          / v_fat
    c_rest     <- a_rest_of_body / v_rest
    Cc         <- a_plasma       / v_plasma

    # ================================================================
    # 5. Unbound, partition-corrected driving concentrations
    # ================================================================
    fu_bp     <- fu / Kbp
    cu_plasma <- Cc       * fu_bp
    cu_gut    <- c_gut    * fu_bp / Kgut_plasma
    cu_liver  <- c_liver  * fu_bp / Kliver_plasma
    cu_kidney <- c_kidney * fu_bp / Kkidney_plasma
    cu_fat    <- c_fat    * fu_bp / Kfat_plasma
    cu_rest   <- c_rest   * fu_bp / Krestbody_plasma

    # ================================================================
    # 6. Renal processes. Secretion is driven by the unbound KIDNEY-TISSUE
    #    concentration, per the deposited code; see the rabbit sibling's
    #    comment and the vignette Errata for the printed-equation conflict.
    # ================================================================
    filtration   <- Qfiltrate * cu_kidney
    reabsorption <- tr * c_filtrate
    secretion    <- tm * cu_kidney / (Kt + cu_kidney)

    # ================================================================
    # 7. Hepatic processes
    # ================================================================
    metabolism     <- vmax * c_liver * fu / (km   + c_liver * fu)
    biliary_efflux <- vehr * c_liver * fu / (Kehr + c_liver * fu)

    # ================================================================
    # 8. ODE system -- 13 mass-balance amount states (mg).
    #    Dose oral into a_gut_lumen, IV into a_plasma.
    #    The faecal gain uses a_gut_lumen (the mass-conserving human-file form),
    #    not the deposited rat file's a_gut; see the rabbit sibling's comment
    #    and the vignette Errata.
    # ================================================================
    d/dt(a_gut_lumen)     <- -kgutabs_s * a_gut_lumen - kfeces * a_gut_lumen + kbile_s * a_bile
    d/dt(a_gut)           <-  kgutabs_s * a_gut_lumen + q_gut * (cu_plasma - cu_gut)
    d/dt(a_liver)         <-  q_liver * cu_plasma + q_gut * cu_gut -
                              (q_liver + q_gut) * cu_liver - metabolism - biliary_efflux
    d/dt(a_hepatic)       <-  metabolism
    d/dt(a_bile)          <-  biliary_efflux - kbile_s * a_bile
    d/dt(a_feces)         <-  kfeces * a_gut_lumen
    d/dt(a_kidney)        <-  q_kidney * cu_plasma - q_kidney * cu_kidney -
                              filtration + reabsorption - secretion
    d/dt(a_filtrate)      <-  filtration - Qfiltrate * c_filtrate - reabsorption + secretion
    d/dt(a_urine_storage) <-  Qfiltrate * c_filtrate - qurine * a_urine_storage
    d/dt(a_fat)           <-  q_fat  * cu_plasma - q_fat  * cu_fat
    d/dt(a_rest_of_body)  <-  q_rest * cu_plasma - q_rest * cu_rest
    d/dt(a_plasma)        <-  (q_liver + q_gut) * cu_liver + q_kidney * cu_kidney +
                              q_fat * cu_fat + q_rest * cu_rest - qc_plasma * cu_plasma
    d/dt(a_urine)         <-  qurine * a_urine_storage

    # ================================================================
    # 9. Diagnostics from the deposited CalcOutputs{} block
    # ================================================================
    mass_balance   <- a_gut_lumen + a_gut + a_liver + a_hepatic + a_bile + a_feces +
                      a_kidney + a_filtrate + a_urine_storage + a_fat +
                      a_rest_of_body + a_plasma + a_urine
    total_renal_cl <- filtration + secretion - reabsorption

    # ================================================================
    # 10. Observation: total plasma concentration (mg/L)
    # ================================================================
    Cc ~ prop(propSd)
  })
}
