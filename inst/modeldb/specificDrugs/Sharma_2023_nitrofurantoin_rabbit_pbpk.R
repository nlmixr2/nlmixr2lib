Sharma_2023_nitrofurantoin_rabbit_pbpk <- function() {
  description <- "Preclinical (rabbit). PBPK (whole-body, GNU MCSim 6.1, Bayesian MCMC). Nitrofurantoin (NFT) disposition in rabbits -- Sharma 2023 'Model V5a', the paper's final structure. Thirteen mass-balance amount states: gut lumen, gut tissue, liver, cumulative hepatic metabolism, bile, cumulative faeces, kidney, renal tubular filtrate, pre-void urine storage, fat, lumped rest of body, plasma, cumulative urine. Perfusion-rate-limited tissue distribution with a single combined plasma pool (no arterial/venous split); only unbound drug (fu/Kbp, divided by each tissue:plasma partition coefficient) distributes. Renal handling is deliberately non-linear: passive glomerular filtration plus SATURABLE active tubular secretion out of the kidney into the tubules, minus first-order tubular reabsorption from the tubules back into the kidney -- the combination the paper shows is required to reproduce the observed fall in fractional urinary recovery from ~60% to ~10% over a 30-fold IV dose increase. Hepatic elimination is Michaelis-Menten metabolism plus a second saturable hepatobiliary efflux into bile, which returns to the gut lumen at a first-order rate (enterohepatic recirculation, EHR). Renal parameters were calibrated on rabbit IV plasma and urine data (model V4); the EHR parameters are the rabbit arm of a hierarchical rat+rabbit fit (model V5a), and are small in rabbits (biliary excretion < 1% of dose), so EHR mainly adds late-phase biphasic plasma behaviour rather than changing urinary recovery. Deterministic typical-value simulator: the MCMC posterior describes PARAMETER UNCERTAINTY, not between-subject variability, so no IIV is encoded (see the vignette Errata)."
  reference <- paste(
    "Sharma RP, Burgers EJ, Beltman JB. (2023).",
    "Development of a Physiologically Based Pharmacokinetic Model for",
    "Nitrofurantoin in Rabbits, Rats, and Humans.",
    "Pharmaceutics 15(9):2199. doi:10.3390/pharmaceutics15092199. PMCID PMC10535763.",
    "Model equations and parameters: Supplementary Materials Tables S2 and S4 and the",
    "'Standard ordinary differential equations used in PBPK model for NFT' section.",
    "Author-deposited GNU MCSim source and MCMC posterior chains:",
    "doi:10.5281/zenodo.8276305 (file Extrapolated_ratEHR_rabbit.model.R = rabbit model V5a).",
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

  # Every state is an AMOUNT in mg (the MCSim `States` block). Concentrations are
  # derived algebraically in model() as amount / organ volume.
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

  # No covariates. Rabbit physiology is a fixed 2.5 kg reference animal
  # (Supplementary Table S2); the deposited rabbit model file hard-codes BW and
  # every organ fraction. Users who want a different body weight override BW via
  # rxSolve(params = c(BW = ...)), which rescales every volume, flow and
  # BW-scaled rate parameter consistently.
  covariateData <- list()

  population <- list(
    species        = "rabbit",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "2.5 kg reference rabbit (Supplementary Table S2)",
    weight_median  = "2.5 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy rabbits (kinetic study, no disease model)",
    dose_range     = "0.25, 1.25, 2.5, 5, 10 and 15 mg/kg single IV; oral dosing also simulated",
    notes          = paste(
      "Calibration data are DIGITISED literature values, not an individual-level dataset:",
      "Sharma 2023 Methods 2.3 states plasma and urine kinetics after single oral or IV",
      "doses of 0.25, 1.25, 2.5, 5, 10 and 15 mg/kg to rabbits were extracted from",
      "Watari N, Aizawa K, Kaneniwa N (1985) J Pharm Sci 74:165-170 using WebPlotDigitizer.",
      "Fitted points are the MEAN of the experimental data at each time point, so no",
      "subject count, age range or sex split is recoverable from the source.",
      "Figure 4 of the paper reports IV doses as 0.5, 1.25, 2.5, 5, 10 and 15 mg/kg",
      "whereas Methods 2.3 lists the lowest dose as 0.25 mg/kg -- see vignette Errata."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Rabbit physiology -- Supplementary Table S2, cross-checked against the
    # deposited Extrapolated_ratEHR_rabbit.model.R. Where the two differ the
    # deposited file simply carries an extra significant figure (FQfat 0.0604
    # vs 0.06; FQgut 0.2094 vs 0.209; Fplasma 0.08594 vs 0.085; Qfiltrate
    # 0.631 vs 0.62), so the deposited value is used. All literature-sourced
    # or study-design constants -> fixed().
    # ------------------------------------------------------------------
    BW         <- fixed(2.5)      ; label("Body weight (kg)")                                          # Table S2 'Body weight'; deposited file BW = 2.5
    QCC        <- fixed(15.96)    ; label("Cardiac output coefficient (L/h/kg^0.75)")                   # Table S2 'Cardiac blood output' (Brown 1997; Davies and Morris 1993)
    FQliver    <- fixed(0.1245)   ; label("Fraction of cardiac output to liver (unitless)")             # Table S2 'Fractional liver blood flow'
    FQkidney   <- fixed(0.151)    ; label("Fraction of cardiac output to kidney (unitless)")            # Table S2 'Fractional kidney blood flow'
    FQfat      <- fixed(0.0604)   ; label("Fraction of cardiac output to fat (unitless)")               # deposited file FQfat = 0.0604 (Table S2 rounds to 0.06)
    FQgut      <- fixed(0.2094)   ; label("Fraction of cardiac output to gut (unitless)")               # deposited file FQgut = 0.2094 (Table S2 rounds to 0.209)
    Fliver     <- fixed(0.04)     ; label("Liver volume as a fraction of body weight (unitless)")       # Table S2 'Fractional liver volume'
    Fkidney    <- fixed(0.006)    ; label("Kidney volume as a fraction of body weight (unitless)")      # Table S2 'Fractional kidney volume'
    Ffiltrate  <- fixed(0.0006)   ; label("Tubular filtrate volume as a fraction of body weight (unitless)")  # Table S2 'Fractional filtrate volume' (10% of kidney volume)
    Ffat       <- fixed(0.048)    ; label("Fat volume as a fraction of body weight (unitless)")         # Table S2 'Fractional fat volume'
    Fgut       <- fixed(0.048)    ; label("Gut volume as a fraction of body weight (unitless)")         # Table S2 'Fractional gut volume'
    Fplasma    <- fixed(0.08594)  ; label("Plasma volume as a fraction of body weight (unitless)")      # deposited file Fplasma = 0.08594 (Table S2 rounds to 0.085)
    Qfiltrate  <- fixed(0.631)    ; label("Glomerular filtration rate (L/h)")                           # deposited file Qfilterate = 0.631 (Table S2 rounds to 0.62; Michigoshi 2012)

    # ------------------------------------------------------------------
    # Tissue:plasma partition coefficients. Table S4 Type = "Calculated"
    # (Berezhkovskiy method, Methods 2.2) -> fixed(). Krestbody_plasma is the
    # one Type = "Fitted" partition coefficient (Methods 2.2: "except for the
    # Restbody:plasma partition coefficient, which we calibrated"), so it is
    # left estimable.
    # ------------------------------------------------------------------
    Kgut_plasma      <- fixed(0.622) ; label("Gut:plasma partition coefficient (unitless)")             # Table S4 Kgut_plasma, Calculated
    Kliver_plasma    <- fixed(0.651) ; label("Liver:plasma partition coefficient (unitless)")           # Table S4 Kliver_plasma, Calculated
    Kkidney_plasma   <- fixed(0.671) ; label("Kidney:plasma partition coefficient (unitless)")          # Table S4 Kkidney_plasma, Calculated
    Kfat_plasma      <- fixed(0.159) ; label("Fat:plasma partition coefficient (unitless)")             # Table S4 Kfat_plasma, Calculated
    Krestbody_plasma <- 0.423        ; label("Rest-of-body:plasma partition coefficient (unitless)")    # Table S4, Fitted, 0.423 (0.39-0.45); V4 posterior mean 0.4213

    fu  <- fixed(0.4149) ; label("Fraction unbound in plasma (unitless)")                               # deposited file fu = 0.4149; Table S4 rounds to 0.42 (Watari 1985)
    Kbp <- fixed(0.76)   ; label("Blood:plasma partition coefficient (unitless)")                       # Table S4 Kbp = 0.76 (Zhang 2022)

    # ------------------------------------------------------------------
    # Renal parameters -- Table S4, Type = "Fitted" (MCMC on rabbit IV plasma
    # and urine data, model V4). Estimable (not fixed()).
    # QurineC: Table S4 PRINTS 11.45 (7.92-20.24), but all five deposited
    # MCSim model files and the human Monte-Carlo input file use 13.45, and the
    # V4 posterior chains give mean 13.456 with 2.5-97.5 percentiles
    # 7.929-20.249 -- i.e. the printed interval matches the chains exactly
    # while the printed point estimate does not. 11.45 is a digit
    # transposition of 13.45; the fitted value is used. See vignette Errata.
    # ------------------------------------------------------------------
    QurineC <- 13.45 ; label("Urine excretion rate coefficient (1/h/kg)")                               # deposited files + V4 posterior mean 13.456; Table S4 typo 11.45
    Trc     <- 1.33  ; label("Tubular reabsorption rate coefficient (L/h/kg)")                          # Table S4 krc, Fitted, 1.33 (1.11-1.56); V4 posterior mean 1.3349
    Tmc     <- 8.02  ; label("Maximum active tubular secretion rate coefficient (mg/h/kg)")             # Table S4 Tmc, Fitted, 8.02 (6.78-9.30); V4 posterior mean 8.0241
    Kt      <- 0.059 ; label("Concentration at half-maximal active tubular secretion (mg/L)")           # Table S4 Kt, Fitted, 0.059 (0.043-0.079); V4 posterior mean 0.0593

    # ------------------------------------------------------------------
    # Hepatic metabolism -- Table S4, Fitted (model V4).
    # ------------------------------------------------------------------
    VmaxC <- 0.47 ; label("Maximum hepatic metabolism rate coefficient (mg/h/g liver)")                 # Table S4 VmaxC, Fitted, 0.47 (0.42-0.53); V4 posterior mean 0.4729
    Km    <- 5.83 ; label("Concentration at half-maximal hepatic metabolism (mg/L)")                    # Table S4 Km, Fitted, 5.83 (4.96-6.82); V4 posterior mean 5.8324

    # ------------------------------------------------------------------
    # Enterohepatic recirculation -- Table S4 RABBIT column, the rabbit arm of
    # the hierarchical rat+rabbit model V5a fit (Results 3.3). Rabbit biliary
    # efflux capacity is ~24-fold lower than the rat's (Vehrc 0.022 vs 0.52),
    # which is how the model reproduces "in rabbits this excretion appears to
    # be low (less than 1%)" (Discussion).
    # ------------------------------------------------------------------
    Vehrc   <- 0.022 ; label("Maximum hepatobiliary excretion rate coefficient (mg/h/kg)")              # Table S4 Vehrc rabbit column, Fitted, 0.022 (0.0012-0.104)
    Kehr    <- 3.69  ; label("Concentration at half-maximal hepatobiliary excretion (mg/L)")            # Table S4 Kehr rabbit column, Fitted, 3.69 (0.74-7.81)
    kbile   <- 3.36  ; label("Bile-to-gut-lumen transfer rate coefficient (1/h/kg)")                    # Table S4 kbile rabbit column, Fitted, 3.36 (0.2-8.98)
    kgutabs <- 0.30  ; label("Gut-lumen absorption rate coefficient (1/h/kg)")                          # Table S4 kgutabs rabbit column, Fitted, 0.30 (0.07-0.66)
    kfeces  <- 3.34  ; label("Faecal excretion rate constant from gut lumen (1/h)")                     # Table S4 kfeces rabbit column, Fitted, 3.34 (0.23-7.46)

    # ------------------------------------------------------------------
    # Residual error. Methods 2.3: "The likelihood of the data was considered
    # to follow a normal distribution with a coefficient of variation of 10%."
    # That is a 10% proportional residual error, and it was an assumption of
    # the fit rather than an estimated quantity -> fixed().
    # No IIV: the MCMC posterior quantifies PARAMETER UNCERTAINTY across a
    # digitised mean-data fit, not between-animal variability.
    # ------------------------------------------------------------------
    propSd <- fixed(0.10) ; label("Proportional residual error (fraction)")                             # Methods 2.3 'coefficient of variation of 10%'
  })

  model({
    # ================================================================
    # 1. Cardiac output and organ plasma flows (L/h)
    #    Deposited Initialize{}: QCblood = QCC * BW^0.75 and
    #    QCplasma = QCblood -- the haematocrit correction
    #    (QCplasma = QCblood * (1 - HCT)) is commented out in every deposited
    #    model file, so HCT is NOT used and is not declared in ini().
    # ================================================================
    qc_plasma <- QCC * BW^0.75
    q_liver   <- FQliver  * qc_plasma
    q_gut     <- FQgut    * qc_plasma
    q_kidney  <- FQkidney * qc_plasma
    q_fat     <- FQfat    * qc_plasma
    # Rest-of-body flow closes the circulation by subtraction (Methods 2.1).
    q_rest    <- qc_plasma - (q_liver + q_kidney + q_fat + q_gut)

    # ================================================================
    # 2. Organ volumes (L). Rest-of-body volume closes total body VOLUME by
    #    subtraction; the deposited file uses 0.84 * BW as total body volume
    #    (i.e. it excludes the ~16% of body mass -- bone, and other
    #    non-perfused mass -- that the model does not distribute drug into).
    # ================================================================
    v_gut      <- Fgut      * BW
    v_liver    <- Fliver    * BW
    v_kidney   <- Fkidney   * BW
    v_filtrate <- Ffiltrate * BW
    v_fat      <- Ffat      * BW
    v_plasma   <- Fplasma   * BW
    v_rest     <- 0.84 * BW - v_liver - v_kidney - v_fat - v_plasma - v_gut - v_filtrate

    # ================================================================
    # 3. Body-weight scaling of the fitted rate parameters (deposited
    #    Initialize{} block). Vmax additionally converts mg/h/g liver to
    #    mg/h for the whole liver: * 1000 g/kg * v_liver (L ~ kg).
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
    c_gut      <- a_gut           / v_gut
    c_liver    <- a_liver         / v_liver
    c_kidney   <- a_kidney        / v_kidney
    c_filtrate <- a_filtrate      / v_filtrate
    c_fat      <- a_fat           / v_fat
    c_rest     <- a_rest_of_body  / v_rest
    Cc         <- a_plasma        / v_plasma

    # ================================================================
    # 5. Unbound, partition-corrected driving concentrations.
    #    Every flow term in the supplement's ODEs carries fu/Kbp on the
    #    plasma side and (fu/Kbp)/K<tissue>_plasma on the tissue side.
    # ================================================================
    fu_bp     <- fu / Kbp
    cu_plasma <- Cc      * fu_bp
    cu_gut    <- c_gut   * fu_bp / Kgut_plasma
    cu_liver  <- c_liver * fu_bp / Kliver_plasma
    cu_kidney <- c_kidney * fu_bp / Kkidney_plasma
    cu_fat    <- c_fat   * fu_bp / Kfat_plasma
    cu_rest   <- c_rest  * fu_bp / Krestbody_plasma

    # ================================================================
    # 6. Renal processes. Secretion is driven by the UNBOUND KIDNEY-TISSUE
    #    concentration (active transport out of the kidney into the tubular
    #    lumen), per every deposited model file. The supplement's printed
    #    Akidney / Afiltrate equations render the secretion term with
    #    c_tubules in place of the kidney concentration; that form would make
    #    lumen-driven transport INTO the lumen a positive feedback and is
    #    inconsistent with the paper's own "active secretion saturates at low
    #    PLASMA concentrations" analysis (Results 3.2, Figure 4). The code
    #    form is used. See vignette Errata.
    # ================================================================
    filtration   <- Qfiltrate * cu_kidney
    reabsorption <- tr * c_filtrate
    secretion    <- tm * cu_kidney / (Kt + cu_kidney)

    # ================================================================
    # 7. Hepatic processes: Michaelis-Menten metabolism plus saturable
    #    hepatobiliary efflux. Both are driven by the unbound TOTAL liver
    #    concentration (c_liver * fu) -- the supplement's ALiver and ABile
    #    equations do NOT divide by Kliver_plasma in these two terms.
    # ================================================================
    metabolism     <- vmax * c_liver * fu / (Km   + c_liver * fu)
    biliary_efflux <- vehr * c_liver * fu / (Kehr + c_liver * fu)

    # ================================================================
    # 8. ODE system -- 13 mass-balance amount states (mg).
    #    Dose oral into a_gut_lumen, IV into a_plasma.
    #
    #    NOTE on the faecal term: the loss from the gut LUMEN is
    #    -kfeces * a_gut_lumen in the supplement's printed Agutlumen equation
    #    and in every deposited model file. The deposited RABBIT and RAT files
    #    then write the matching gain as `dt(Afeces) = kfeces * Agut` (gut
    #    TISSUE), which does not conserve mass; the deposited HUMAN file
    #    writes `kfeces * Agutlumen`, which does. The mass-conserving human
    #    form is used in all three species files -- it is the only reading
    #    under which the paper's own Massbalance output is constant. See
    #    vignette Errata.
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
    # 9. Diagnostics reproduced from the deposited CalcOutputs{} block.
    #    mass_balance is the paper's Massbalance output and must equal the
    #    administered dose at all times -- the primary validation gate for
    #    this model (see the vignette).
    # ================================================================
    mass_balance   <- a_gut_lumen + a_gut + a_liver + a_hepatic + a_bile + a_feces +
                      a_kidney + a_filtrate + a_urine_storage + a_fat +
                      a_rest_of_body + a_plasma + a_urine
    total_renal_cl <- filtration + secretion - reabsorption

    # ================================================================
    # 10. Observation: total plasma concentration (mg/L).
    # ================================================================
    Cc ~ prop(propSd)
  })
}
