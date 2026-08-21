Boone_2025_vinylchloride_pbpk <- function() {
  description <- "PBPK (whole-body, flow-limited, seven tissue compartments; recoded from Berkeley Madonna 10.4.2). Generic ATSDR volatile-organic-compound human PBPK model parameterised for vinyl chloride (VC) and applied by Boone et al. 2025 to reconstruct community exposure after the 2012 Paulsboro, New Jersey train derailment. Venous blood, rapidly perfused, slowly perfused, fat, liver, kidney and skin are flow-limited well-mixed compartments; arterial blood is an algebraic steady-state lung equation coupling cardiac output to alveolar ventilation through the blood:air partition coefficient, so inhalation exposure enters as an air concentration rather than as a dose record. Saturable Michaelis-Menten hepatic metabolism. Three exposure routes are implemented simultaneously and driven by exposure parameters exactly as in the published Berkeley Madonna listing: inhalation of a square-wave air concentration in ppm, repeated oral ingestion of contaminated drinking water, and dermal transfer from water across the skin. Deterministic typical-value model (no IIV, no residual error): the authors recoded the published Clewell VC model onto the ATSDR toolkit platform and evaluated it against published human kinetic data rather than fitting it. Outputs are arterial blood, venous blood and exhaled (alveolar) breath concentrations plus cumulative metabolised, inhaled, exhaled, ingested and dermally transferred amounts for mass-balance checking."
  reference <- paste(
    "Boone S, Sun W, Gonnabathula P, Wu J, Orr MF, Mumtaz MM, Ruiz P. (2025).",
    "Assessing the Application of Physiologically Based Pharmacokinetic Models",
    "in Acute Chemical Incidents.",
    "Journal of Xenobiotics 15(2):42.",
    "doi:10.3390/jox15020042. PMCID PMC11932312.",
    "Model code: Supplementary Materials, 'Berkeley Madonna code for PBPK model of Vinyl Chloride'.",
    sep = " "
  )
  vignette <- "Boone_2025_vinylchloride"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L",
    amount        = "mg",
    weight        = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "The Berkeley Madonna listing hard-codes BW = 70.0 kg (cited to Clewell et al. 2000)",
        "and every published simulation in the paper uses that reference adult. Carried here as",
        "a covariate so the model can be re-scaled: all seven tissue volumes are linear in body",
        "weight, and cardiac output, alveolar ventilation and the metabolic Vmax are allometric",
        "at exponent 0.75. The dermal surface area SA and the skin permeability Kp are NOT",
        "body-weight scaled in the published code and are left as absolute constants.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    a_venous             = list(analyte = "Vinyl chloride", units = "mg", specimen = "whole blood", verified = TRUE),
    a_rapidly_perfused   = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_slowly_perfused    = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_fat                = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_liver              = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_kidney             = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_skin               = list(analyte = "Vinyl chloride", units = "mg", specimen = "tissue",       verified = TRUE),
    a_metabolized        = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",     verified = TRUE),
    a_inhaled            = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",         verified = TRUE),
    a_exhaled            = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",         verified = TRUE),
    a_oral               = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",        verified = TRUE),
    a_dermal_absorbed    = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",   verified = TRUE),
    a_dermal_eliminated  = list(analyte = "Vinyl chloride", units = "mg", specimen = "not applicable",  verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 0L,
    n_studies      = 0L,
    age_range      = "adult (single 70 kg reference adult; no age term in the model)",
    weight_range   = "70 kg (Berkeley Madonna listing, BW = 70.0, cited to Clewell et al. 2000)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "healthy reference adult",
    dose_range     = paste(
      "Inhalation only in the published simulations. Acute Exposure Guideline Level (AEGL)",
      "air concentrations of 70-12,000 ppm over 10 min to 8 h (Table 4 / Table S5), and",
      "measured Paulsboro air concentrations of 0.19-1649.2 ppm simulated as 24 h exposures",
      "(Table 3 / Tables S3-S4). The oral and dermal routes are implemented in the code but",
      "are switched off in every published simulation (pdose = 0, Cliq = 0, skin_time = 0).",
      sep = " "
    ),
    regions        = "United States (ATSDR / CDC; case study in Paulsboro, New Jersey)",
    notes          = paste(
      "This is a forward EXPOSURE-RECONSTRUCTION simulation, not a fit, so n_subjects = 0.",
      "ATSDR recoded the published Clewell et al. vinyl-chloride PBPK model into a generic",
      "volatile-organic-compound structure on the Berkeley Madonna platform; the paper states",
      "the recoded model's performance 'was found adequate based on a comparison with the",
      "published human kinetic data for VC'. Physiological parameters are attributed in the",
      "code listing to Brown et al. 1997, Clewell et al. 2000/2001/2005, Covington et al. 2007",
      "and Fisher, Mahle and Abbas 1998; the chemical-specific metabolic constants to Reitz 1996;",
      "the dermal constants to Poet et al. 2000. Metabolite kinetics were deliberately excluded",
      "from the generic model (Methods 2.2). The exposed community comprised adults and children",
      "and 250 hospital visits followed the release, but no individual-level kinetic data were",
      "collected, so the model carries no between-subject variability. The paper instead applies",
      "an assessment factor of 3.16 (the square root of 10) outside the model to account for",
      "interindividual toxicokinetic variability in the AEGL blood-level ranges.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Physiological parameters -- Supplementary Berkeley Madonna listing.
    # Every value in this model is a recoded literature constant, so every
    # parameter is fixed(); nothing was estimated from data.
    # =====================================================================
    qp_per_kg           <- fixed(24)     ; label("Alveolar ventilation rate (L/h/kg^0.75)")             # BM listing, QPC=24, T. R. Covington et al., 2007
    qc_per_kg           <- fixed(16.5)   ; label("Cardiac output (L/h/kg^0.75)")                        # BM listing, QCC=16.5, Clewell et al., 2000
    e_wt_flow           <- fixed(0.75)   ; label("Allometric exponent on cardiac output, ventilation and Vmax (unitless)")  # BM listing, QC=QCC*BW**.75, QP=QPC*BW**.75, VMAX=Vmaxc*BW**.75

    # Tissue volumes, fraction of body weight. NB the slowly-perfused and
    # rapidly-perfused fractions are TOP-LEVEL splits from which the named
    # organs are carved out (VS = VSC*BW - VF - VSk; VR = VRC*BW - VL - VK),
    # so they are not the volumes of the lumped remainders themselves.
    fv_fat              <- fixed(0.214)  ; label("Fat volume, fraction of body weight")                 # BM listing, VFC=.214, Clewell et al., 2000
    fv_liver            <- fixed(0.026)  ; label("Liver volume, fraction of body weight")               # BM listing, VLC=.026, Clewell et al., 2000
    fv_blood            <- fixed(0.079)  ; label("Venous blood volume, fraction of body weight")        # BM listing, VBloodC=.079
    fv_rapidly_perfused <- fixed(0.09)   ; label("Rapidly-perfused group volume (liver + kidney + remainder), fraction of body weight")  # BM listing, VRC=.09
    fv_slowly_perfused  <- fixed(0.82)   ; label("Slowly-perfused group volume (fat + skin + remainder), fraction of body weight")       # BM listing, VSC=.82
    fv_skin             <- fixed(0.051)  ; label("Skin volume, fraction of body weight")                # BM listing, VSkC=.051
    fv_kidney           <- fixed(0.044)  ; label("Kidney volume, fraction of body weight")              # BM listing, VKC=.044 (see vignette Errata: an order of magnitude above the Brown et al. 1997 human kidney fraction)

    # Blood flows, fraction of cardiac output. As with the volumes, the
    # rapidly- and slowly-perfused fractions are top-level splits: they sum
    # to exactly 1.0 and the named organs are carved out of them.
    fq_fat              <- fixed(0.052)  ; label("Fat blood flow, fraction of cardiac output")          # BM listing, QFC=.052, Clewell et al., 2000
    fq_liver            <- fixed(0.24)   ; label("Liver blood flow, fraction of cardiac output")        # BM listing, QLC=.24, Fisher Mahle and Abbas, 1998
    fq_rapidly_perfused <- fixed(0.7)    ; label("Rapidly-perfused group blood flow (liver + kidney + remainder), fraction of cardiac output")  # BM listing, QRC=0.7, Clewell, 2005
    fq_slowly_perfused  <- fixed(0.30)   ; label("Slowly-perfused group blood flow (fat + skin + remainder), fraction of cardiac output")       # BM listing, QSC=.30, Clewell, 2005
    fq_kidney           <- fixed(0.197)  ; label("Kidney blood flow, fraction of cardiac output")       # BM listing, QKC=.197
    fq_skin             <- fixed(0.05)   ; label("Skin blood flow, fraction of cardiac output")         # BM listing, QSkC=.05

    sa_skin             <- fixed(19975)  ; label("Exposed skin surface area, body minus head (cm^2)")   # BM listing, SA=19975

    # =====================================================================
    # Chemical-specific parameters for vinyl chloride -- BM listing block
    # "{Chemical-Specific Parameters for VC}".
    # =====================================================================
    pc_blood_air        <- fixed(1.16)   ; label("Blood:air partition coefficient (unitless)")          # BM listing, PB=1.16, Clewell, 2001
    pc_liver            <- fixed(1.45)   ; label("Liver:blood partition coefficient (unitless)")        # BM listing, PL=1.45, Clewell, 2001
    pc_fat              <- fixed(20.7)   ; label("Fat:blood partition coefficient (unitless)")          # BM listing, PF=20.7, Clewell, 2001
    pc_rapidly_perfused <- fixed(1.45)   ; label("Rapidly-perfused:blood partition coefficient (unitless)")  # BM listing, PR=1.45, Clewell, 2001
    pc_slowly_perfused  <- fixed(0.83)   ; label("Slowly-perfused:blood partition coefficient (unitless)")   # BM listing, PS=.83, Clewell, 2001
    pc_kidney           <- fixed(1.45)   ; label("Kidney:blood partition coefficient (unitless)")       # BM listing, PK=1.45 (see liver)
    pc_skin             <- fixed(1.45)   ; label("Skin:blood partition coefficient (unitless)")         # BM listing, PSk=1.45, Poet et al., 2000
    mw                  <- fixed(62.5)   ; label("Molecular weight of vinyl chloride (g/mol)")          # BM listing, MW=62.5, EPA 2000
    vmax_per_kg         <- fixed(3.97)   ; label("Maximum rate of hepatic metabolism (mg/h/kg^0.75)")   # BM listing, Vmaxc=3.97, Reitz, 1996
    km                  <- fixed(0.04)   ; label("Michaelis-Menten constant for hepatic metabolism (mg/L)")  # BM listing, Km=.04, Reitz, 1996

    # Molar volume used by the listing to convert ppm to mg/L and back
    # (CIX=CONC*MW/24450; CXppm=(.7*CX+.3*CI)*24450/MW). 24450 mL/mol is the
    # ideal-gas molar volume at 25 degrees C and 1 atm.
    molar_vol           <- fixed(24450)  ; label("Ideal-gas molar volume used for the ppm <-> mg/L conversion (mL/mol)")  # BM listing, CIX=CONC*MW/24450
    f_alveolar          <- fixed(0.7)    ; label("Alveolar fraction of an exhaled breath sample; the remainder is unabsorbed inhaled air (unitless)")  # BM listing, CXppm=(.7*CX+.3*CI)*24450/MW

    # =====================================================================
    # Dermal-route constants -- BM listing block ";Dermal".
    # =====================================================================
    perm_skin             <- fixed(0.015)  ; label("Skin permeability constant (cm/h)")                   # BM listing, Kp=.015, Poet et al., 2000
    pc_skin_water       <- fixed(53)     ; label("Skin:water partition coefficient (unitless)")          # BM listing, PSkliq=53, Poet et al., 2000

    # =====================================================================
    # Exposure parameters -- BM listing block "{Exposure Parameters}".
    # These are the sliders a Berkeley Madonna user sets to define a
    # scenario; the published defaults switch every route off. They are
    # fixed() at the listing's defaults and overridden per scenario through
    # rxSolve(params = ...). The paper's simulations set conc_air_ppm and
    # inhale_time and leave the oral and dermal routes at zero.
    # =====================================================================
    conc_air_ppm        <- fixed(0)      ; label("Inhaled vinyl-chloride air concentration during an exposure window (ppm)")  # BM listing, CONC=0, default set=0
    inhale_time         <- fixed(5)      ; label("Duration of each inhalation exposure window (h)")     # BM listing, inhale_time=5
    inhale_interval     <- fixed(100000) ; label("Period between repeated inhalation exposures (h); a large value gives a single window")  # BM listing, inhale_interval=100000
    dose_oral           <- fixed(0)      ; label("Oral vinyl-chloride dose from drinking water (mg/kg/day)")  # BM listing, pdose=0, default=0
    drink_time          <- fixed(0.25)   ; label("Duration of each drinking-water intake (h)")          # BM listing, drink_time=.25
    drink_interval      <- fixed(6)      ; label("Period between drinks (h)")                           # BM listing, drink_interval=6
    conc_water          <- fixed(0)      ; label("Vinyl-chloride concentration in water contacting the skin (mg/L)")  # BM listing, Cliq=0, Dr. Fisher, default=0
    skin_time           <- fixed(0)      ; label("Duration of each dermal exposure window (h); 0 switches the dermal route off")  # BM listing, skin_time=0, SET =to 0 when no dermal exposure
    skin_interval       <- fixed(100000) ; label("Period between repeated dermal exposures (h); a large value gives a single window")  # BM listing, skin_interval=100000

    # =====================================================================
    # No IIV and no residual error. The authors recoded a published
    # deterministic PBPK model and evaluated it against digitised human
    # kinetic data; no individual-level data were fitted and the paper
    # explicitly defers inter-individual variability to future work
    # (Discussion), handling it outside the model with a factor of 3.16.
    # =====================================================================
  })

  model({
    # -------------------------------------------------------------------
    # 1. Calculated volumes (L) and blood flows (L/h).
    #    BM listing block "{Calculated Parameters}". Note the carve-out
    #    ordering: VS and VR (and QS and QR) are the group totals MINUS the
    #    organs modelled explicitly inside each group.
    # -------------------------------------------------------------------
    v_fat               <- fv_fat    * WT
    v_skin              <- fv_skin   * WT
    v_kidney            <- fv_kidney * WT
    v_liver             <- fv_liver  * WT
    v_blood             <- fv_blood  * WT
    v_slowly_perfused   <- fv_slowly_perfused  * WT - v_fat   - v_skin
    v_rapidly_perfused  <- fv_rapidly_perfused * WT - v_liver - v_kidney

    qc                  <- qc_per_kg * WT^e_wt_flow
    qp                  <- qp_per_kg * WT^e_wt_flow

    q_fat               <- fq_fat    * qc
    q_liver             <- fq_liver  * qc
    q_kidney            <- fq_kidney * qc
    q_skin              <- fq_skin   * qc
    q_rapidly_perfused  <- fq_rapidly_perfused * qc - q_liver - q_kidney
    q_slowly_perfused   <- fq_slowly_perfused  * qc - q_fat   - q_skin

    vmax                <- vmax_per_kg * WT^e_wt_flow

    # -------------------------------------------------------------------
    # 2. Exposure schedules. Berkeley Madonna writes each route as a
    #    repeated square wave, e.g.
    #      AIR = IF MOD(TIME,inhale_interval)>=inhale_time THEN 0 ELSE 1
    #    MOD(TIME, P) is rendered here as time - P*floor(time/P), and the
    #    IF/THEN/ELSE as the complement of the comparison, which rxode2
    #    evaluates to 1 or 0. With the listing's default interval of 100000 h
    #    the wave is a single window starting at time 0.
    # -------------------------------------------------------------------
    conc_air_mgl        <- conc_air_ppm * mw / molar_vol            # BM listing, CIX=CONC*MW/24450
    time_mod_inhale     <- time - inhale_interval * floor(time / inhale_interval)
    air_on              <- 1 - (time_mod_inhale >= inhale_time)
    ci                  <- conc_air_mgl * air_on                    # BM listing, CI=CIX*AIR

    # Repeated drinking water. daily_drink_interval caps the interval at 24 h
    # so that a single daily drink still delivers the full daily dose
    # (BM: daily_drink_interval = IF drink_interval>24 THEN 24 ELSE drink_interval).
    daily_drink_interval <- min(drink_interval, 24)
    number_drink        <- 24 / daily_drink_interval
    dose_per_drink      <- dose_oral * WT / number_drink
    doser               <- dose_per_drink / drink_time
    time_mod_drink      <- time - drink_interval * floor(time / drink_interval)
    oral_on             <- 1 - (time_mod_drink >= drink_time)
    r_oral              <- doser * oral_on                          # BM listing, Roral=doser*ORALH20

    time_mod_skin       <- time - skin_interval * floor(time / skin_interval)
    skin_on             <- 1 - (time_mod_skin >= skin_time)

    # -------------------------------------------------------------------
    # 3. Tissue concentrations (mg/L) and the capillary (venous-effluent)
    #    concentrations that leave each flow-limited compartment.
    # -------------------------------------------------------------------
    c_rapidly_perfused  <- a_rapidly_perfused / v_rapidly_perfused
    cv_rapidly_perfused <- c_rapidly_perfused / pc_rapidly_perfused
    c_slowly_perfused   <- a_slowly_perfused / v_slowly_perfused
    cv_slowly_perfused  <- c_slowly_perfused / pc_slowly_perfused
    c_fat               <- a_fat    / v_fat
    cv_fat              <- c_fat    / pc_fat
    c_liver             <- a_liver  / v_liver
    cv_liver            <- c_liver  / pc_liver
    c_kidney            <- a_kidney / v_kidney
    cv_kidney           <- c_kidney / pc_kidney
    c_skin              <- a_skin   / v_skin
    cv_skin             <- c_skin   / pc_skin

    c_venous            <- a_venous / v_blood                       # BM listing, CV=AVBlood/VBlood

    # Arterial blood is algebraic, not a state: the lung is assumed to be at
    # steady state, so arterial blood equilibrates instantaneously between
    # the returning venous blood and the inhaled air across PB.
    c_arterial          <- (qc * c_venous + qp * ci) / (qc + qp / pc_blood_air)   # BM listing, CA=(QC*CV+QP*CI)/(QC+(QP/PB))
    c_alveolar          <- c_arterial / pc_blood_air                # BM listing, CX=CA/PB

    # Saturable hepatic metabolism on the liver capillary concentration.
    r_metab             <- vmax * cv_liver / (km + cv_liver)        # BM listing, RAM=Vmax*CVL/(Km+CVL)

    # Dermal flux into and out of the skin compartment.
    r_dermal            <- (perm_skin * sa_skin / 1000) * (conc_water - c_skin / pc_skin_water) * skin_on   # BM listing, RASkin
    r_dermal_out        <- (perm_skin * sa_skin / 1000) * (c_skin / pc_skin_water) * skin_on                # BM listing, ASkinOut'

    # -------------------------------------------------------------------
    # 4. Mass-balance ODE system, BM listing block "{Model Equations}".
    # -------------------------------------------------------------------
    d/dt(a_venous)            <- (q_fat * cv_fat + q_liver * cv_liver +
                                  q_slowly_perfused * cv_slowly_perfused +
                                  q_rapidly_perfused * cv_rapidly_perfused +
                                  q_kidney * cv_kidney + q_skin * cv_skin) -
                                 qc * c_venous
    d/dt(a_rapidly_perfused)  <- q_rapidly_perfused * (c_arterial - cv_rapidly_perfused)
    d/dt(a_slowly_perfused)   <- q_slowly_perfused  * (c_arterial - cv_slowly_perfused)
    d/dt(a_fat)               <- q_fat    * (c_arterial - cv_fat)
    d/dt(a_liver)             <- q_liver  * (c_arterial - cv_liver) - r_metab + r_oral
    d/dt(a_kidney)            <- q_kidney * (c_arterial - cv_kidney)
    d/dt(a_skin)              <- q_skin   * (c_arterial - cv_skin) + r_dermal

    # Cumulative-amount accumulators, carried by the listing so that the
    # inhalation, oral and dermal mass balances can be checked.
    d/dt(a_metabolized)       <- r_metab                            # BM listing, AM'=RAM
    d/dt(a_inhaled)           <- qp * ci                            # BM listing, AINH'=QP*CI
    d/dt(a_exhaled)           <- qp * c_alveolar                    # BM listing, AX'=QP*CX
    d/dt(a_oral)              <- r_oral                             # BM listing, AORAL'=Roral
    d/dt(a_dermal_absorbed)   <- r_dermal                           # BM listing, ASkin'=RASkin
    d/dt(a_dermal_eliminated) <- r_dermal_out                       # BM listing, ASkinOut'

    # -------------------------------------------------------------------
    # 5. Outputs. Cc is the arterial blood concentration, the quantity the
    #    paper reports as CA throughout Tables 3 and 4 and the one that
    #    would be compared against a biomonitoring blood sample.
    # -------------------------------------------------------------------
    Cc                  <- c_arterial
    Cvenous             <- c_venous
    Cexhaled            <- c_alveolar
    # Exhaled-breath sample as it would be measured at the mouth: 70 percent
    # alveolar air plus 30 percent unabsorbed inhaled air, expressed in ppm.
    Cexhaled_ppm        <- (f_alveolar * c_alveolar + (1 - f_alveolar) * ci) * molar_vol / mw
    Cliver              <- c_liver
    Cfat                <- c_fat
    Ckidney             <- c_kidney
    Cskin               <- c_skin
    Crapid              <- c_rapidly_perfused
    Cslow               <- c_slowly_perfused
    # Net amount absorbed by inhalation (BM listing InhDOSE' = QP*(CI-CX));
    # identical to the difference of the two accumulators above, so it is
    # derived here rather than carried as a redundant state.
    Ainhaled_net        <- a_inhaled - a_exhaled
    # Total body burden: chemical still in tissues plus chemical already
    # metabolised (BM listing, Mass).
    Amass               <- a_fat + a_rapidly_perfused + a_slowly_perfused + a_liver +
                           a_metabolized + a_kidney + a_skin + a_venous
  })
}
