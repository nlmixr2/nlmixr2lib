Sharma_2023_nitrofurantoin_human_pbpk <- function() {
  description <- "PBPK (whole-body, GNU MCSim 6.1, age-dependent physiology). Nitrofurantoin (NFT) disposition in human adults -- Sharma 2023 'Model V5a' extrapolated across species WITHOUT further calibration, the paper's final human model. Same thirteen mass-balance amount states and same structure as the rabbit and rat siblings (perfusion-rate-limited distribution, glomerular filtration + saturable active tubular secretion - first-order tubular reabsorption, Michaelis-Menten hepatic metabolism, saturable hepatobiliary efflux with enterohepatic recirculation). Three things are human-specific. (1) Physiology is not a fixed reference adult: height, body weight, cardiac output, liver / kidney / fat / plasma volumes and liver / kidney / fat blood flows are all AGE- and SEX-dependent polynomials, so AGE and SEXF are genuine covariates that reshape the whole body. (2) Renal and metabolic parameters are allometrically scaled from the RABBIT fit and the enterohepatic-recirculation and gut-absorption parameters from the RAT fit, each as parameter * (BW_species / BW_human)^0.25. (3) Renal function is an explicit input: CRCL carries the subject's absolute glomerular filtration rate in mL/min, which lets the model reproduce the paper's renal-insufficiency analysis (Figure 7) showing that a fall from normal to severely compromised GFR raises plasma Cmax ~1.3-fold, trough ~2-fold and AUC ~1.3-fold while substantially cutting urinary NFT delivery -- simultaneously raising hepatic exposure (a DILI risk) and undercutting efficacy at the bladder. Between-subject variability is a lognormal geometric SD of 1.17 (~16% CV) applied to all 18 drug-specific parameters, as in the paper's Monte-Carlo simulations; the model has no human-fitted parameters at all, so every human prediction is a genuine extrapolation."
  reference <- paste(
    "Sharma RP, Burgers EJ, Beltman JB. (2023).",
    "Development of a Physiologically Based Pharmacokinetic Model for",
    "Nitrofurantoin in Rabbits, Rats, and Humans.",
    "Pharmaceutics 15(9):2199. doi:10.3390/pharmaceutics15092199. PMCID PMC10535763.",
    "Model equations: Supplementary Materials 'Standard ordinary differential equations",
    "used in PBPK model for NFT'. Biochemical parameters and the allometric-scaling note:",
    "Supplementary Table S4. Static human reference physiology: Supplementary Table S3",
    "(superseded in this file by the age-dependent equations actually implemented -- see",
    "the vignette Errata).",
    "Author-deposited GNU MCSim source and MCMC posterior chains:",
    "doi:10.5281/zenodo.8276305 (model file Agedynamics_NFThuman.model.R; parameter means",
    "and the between-subject lognormal geometric SD from Agedynamic_montecarlo.in.R;",
    "renal-insufficiency scenarios from GFRdiseasedAgedynamic_montecarlo.in.R).",
    "Age-dependent physiology equations are attributed by Sharma 2023 Methods 2.2 to its",
    "reference 26; that upstream paper is not on disk, so the equations are transcribed",
    "verbatim from the deposited model file rather than re-derived.",
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

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT a covariate on a clearance term -- age drives the whole body. It enters the",
        "sex-specific 6th-order height and body-weight polynomials, and from those the",
        "Du Bois body surface area, cardiac output, and every organ volume and blood flow.",
        "Sharma 2023 sampled age over 25-40 years for its human simulations, and the",
        "deposited polynomials are only validated over that adult range: they are 6th-order",
        "fits and will diverge outside it. The deposited harness also applies mg/kg rather",
        "than absolute-mg dosing when age <= 12 years, which is a dosing convention of the",
        "authors' simulation script rather than model structure and is therefore not encoded",
        "here -- supply the absolute dose amount in the event table.",
        "The paper's Results 3.4 states age was varied uniformly over 25-40, whereas the",
        "deposited Monte-Carlo input file draws a truncated normal (mean 35, SD 10, bounds",
        "25-40); see the vignette Errata."
      ),
      source_name        = "age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Selects between two complete sets of physiology equations (height, body weight,",
        "liver / kidney / fat / plasma volume, cardiac output, liver / kidney / fat blood",
        "flow, fractional gut blood flow), not a single coefficient.",
        "Orientation was resolved from the deposited code, NOT from its comments: the",
        "deposited file's header says 'model 1 is for boys, and model 2 is for girls', but",
        "its model-1 branch consumes the Fgut_f / FQgut_f constants and its model-2 branch",
        "the Fgut_m / FQgut_m constants, and evaluating the two anthropometry branches at",
        "age 30 gives 169 cm / 66 kg for branch 1 versus 182 cm / 81 kg for branch 2.",
        "Both lines of evidence agree that branch 1 is FEMALE and branch 2 is MALE, so the",
        "header comment is wrong. Mapped here to the canonical SEXF orientation.",
        "Note that the two published human figures used different sexes: the plasma",
        "validation (Figure 6) ran with the deposited default sex = 1, which selects the",
        "MALE branch, while the renal-insufficiency analysis (Figure 7) set sex = 2, which",
        "falls through the deposited ternary to the FEMALE branch. See the vignette Errata."
      ),
      source_name        = "sex"
    ),
    CRCL = list(
      description        = "Absolute glomerular filtration rate (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Absolute GFR in mL/min, converted to filtrate plasma flow in L/h inside model() by",
        "the factor 0.06. This is the model's only renal-function input and it drives",
        "glomerular filtration directly (it is not a scalar on a clearance parameter).",
        "The source paper's renal-insufficiency scenarios are CRCL = 70 (labelled",
        "'moderate'), 45 ('mild') and 20 ('severe') mL/min, versus normal (> 90 mL/min).",
        "NOTE the source's own mild/moderate labels are inverted relative to standard",
        "CKD staging, where 45 mL/min is more impaired than 70 mL/min; the NUMBERS are",
        "unambiguous and are what this model consumes. See the vignette Errata.",
        "For a normal-renal-function subject, set CRCL to the model's own age- and",
        "sex-predicted normal GFR, which is returned on every solve as the output column",
        "`gfr_normal` (mL/min) -- solve once to read it, then re-solve with",
        "CRCL = gfr_normal. That two-pass idiom reproduces the deposited code's CKD == 0",
        "branch exactly without duplicating the body-surface-area polynomial outside the",
        "model."
      ),
      source_name        = "GFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "25-40 years (simulated)",
    age_median     = "35 years",
    weight_range   = "Not an input -- body weight is PREDICTED from age and sex (about 63-68 kg female, 77-83 kg male over the simulated 25-40 year range)",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Healthy adults for the plasma validation; renal insufficiency explored in silico as",
      "absolute GFR 70, 45 and 20 mL/min versus normal. No human data were used to fit any",
      "parameter -- the human model is a pure cross-species extrapolation from the rabbit",
      "(renal, metabolic) and rat (enterohepatic recirculation, gut absorption) fits."
    ),
    dose_range     = "50, 100 and 200 mg single oral (Figure 6); 50 mg orally four times daily for five days (Figure 7)",
    renal_function = "Normal (age/BSA-predicted, about 120-150 mL/min) for the validation; 70, 45 and 20 mL/min for the renal-insufficiency scenarios",
    notes          = paste(
      "The human plasma comparison data are DIGITISED literature values, not an",
      "individual-level dataset: Sharma 2023 Methods 2.3 states human time-course plasma",
      "data after single oral doses of 50, 100 and 200 mg were extracted from its reference",
      "4 using WebPlotDigitizer. No subject count, age range or sex split is reported.",
      "The deposited digitisation is filed as NFT_female_PK.csv, indicating a female source",
      "cohort, even though the deposited Figure-6 simulation ran with male physiology.",
      "Between-subject variability was NOT estimated from data: Results 3.4 states the",
      "authors drew individual parameter values around each fitted mean with a lognormal",
      "geometric SD of 1.17 (~16% CV), a value 'informed by the estimated range of the",
      "posterior parameters' rather than by inter-individual measurements, and explicitly",
      "note that 'detailed NFT measurements amongst multiple individuals would allow for a",
      "more data-informed choice of the CoV across individuals'."
    )
  )

  ini({
    # ==================================================================
    # Human physiology. Supplementary Table S3 tabulates a STATIC 70 kg adult
    # (QCC 4.8 L/h/kg, FQliver 0.25, Qfiltrate 7.2 L/h, Fliver 0.026, ...), but
    # the deposited human model file does not use those values: it implements
    # age- and sex-dependent equations for height, weight, cardiac output, and
    # every organ volume and flow (Methods 2.2: "We used published equations to
    # include the dependence of such physiological parameters ... on age").
    # Those equations are what generated Figures 6 and 7, so they are what this
    # file encodes; the equation COEFFICIENTS are hard-coded in model() rather
    # than declared here because they are polynomial coefficients of a cited
    # upstream body-composition model, not parameters of the NFT model.
    # Only the three physiology constants the deposited file carries as named
    # parameters are declared here.
    # ==================================================================
    Ffiltrate <- fixed(0.0004) ; label("Tubular filtrate volume as a fraction of body weight (unitless)")  # deposited human file Ffilterate = 0.0004 (Table S3 prints 0.00073)
    Fgut      <- fixed(0.013)  ; label("Gut volume as a fraction of body weight (unitless)")               # deposited human file Fgut_f = Fgut_m = 0.013 (Table S3 prints 0.016)
    FQgut_f   <- fixed(0.14)   ; label("Fraction of cardiac output to gut, female (unitless)")             # deposited human file FQgut_f = 0.14
    FQgut_m   <- fixed(0.18)   ; label("Fraction of cardiac output to gut, male (unitless)")               # deposited human file FQgut_m = 0.18

    # ==================================================================
    # Allometric cross-species scaling. Table S4 note: "To scale the value to
    # human we used allometric scaling ... Parameters with a single star (*) use
    # rabbit body weight and with double stars (**) use rat bodyweight as a
    # basis". The deposited human file implements this as
    # parameter * (BW_reference_species / BW_human)^scaling with scaling = 0.25,
    # i.e. rate parameters fall as BW^-0.25 with increasing body size.
    # ==================================================================
    BWrabbit <- fixed(2.5)  ; label("Rabbit reference body weight for allometric scaling (kg)")   # deposited human file BWrabbit = 2.5
    BWrat    <- fixed(0.25) ; label("Rat reference body weight for allometric scaling (kg)")      # deposited human file BWrat = 0.25
    scaling  <- fixed(0.25) ; label("Allometric scaling exponent for rate parameters (unitless)") # deposited human file scaling = 0.25

    # ==================================================================
    # Tissue:plasma partition coefficients -- Table S4, shared across species.
    # Calculated (Berezhkovskiy) -> fixed(); rest-of-body is the fitted one.
    # ==================================================================
    Kgut_plasma      <- fixed(0.622) ; label("Gut:plasma partition coefficient (unitless)")            # Table S4 Kgut_plasma, Calculated
    Kliver_plasma    <- fixed(0.651) ; label("Liver:plasma partition coefficient (unitless)")          # Table S4 Kliver_plasma, Calculated
    Kkidney_plasma   <- fixed(0.671) ; label("Kidney:plasma partition coefficient (unitless)")         # Table S4 Kkidney_plasma, Calculated
    Kfat_plasma      <- fixed(0.159) ; label("Fat:plasma partition coefficient (unitless)")            # Table S4 Kfat_plasma, Calculated
    Krestbody_plasma <- 0.42         ; label("Rest-of-body:plasma partition coefficient (unitless)")   # deposited human Monte-Carlo input 0.42 (Table S4 prints 0.423; V4 posterior mean 0.4213)

    fu  <- fixed(0.41) ; label("Fraction unbound in plasma (unitless)")                                # deposited human Monte-Carlo input 0.41 (Table S4 prints 0.42; Watari 1985)
    Kbp <- fixed(0.76) ; label("Blood:plasma partition coefficient (unitless)")                        # Table S4 Kbp = 0.76 (Zhang 2022)

    # ==================================================================
    # Renal parameters -- fitted on RABBIT IV data (model V4), scaled by rabbit
    # body weight. See the rabbit sibling for the QurineC Table S4 typo
    # resolution (printed 11.45; fitted and deposited 13.45).
    # ==================================================================
    QurineC <- 13.45 ; label("Urine excretion rate coefficient (1/h/kg), rabbit-scaled")              # deposited human file + V4 posterior mean 13.456; Table S4 typo 11.45
    Trc     <- 1.334 ; label("Tubular reabsorption rate coefficient (L/h/kg), rabbit-scaled")          # deposited human file 1.334; Table S4 krc 1.33; V4 posterior mean 1.3349
    Tmc     <- 8.024 ; label("Maximum active tubular secretion rate coefficient (mg/h/kg), rabbit-scaled")  # deposited human file 8.024; Table S4 8.02; V4 posterior mean 8.0241
    Kt      <- 0.059 ; label("Concentration at half-maximal active tubular secretion (mg/L)")          # Table S4 Kt 0.059 (0.043-0.079); V4 posterior mean 0.0593

    # ==================================================================
    # Hepatic metabolism -- fitted on rabbit data, scaled by rabbit body weight.
    # ==================================================================
    VmaxC <- 0.472 ; label("Maximum hepatic metabolism rate coefficient (mg/h/g liver), rabbit-scaled") # deposited human Monte-Carlo input 0.472; Table S4 0.47; V4 posterior mean 0.4729
    km    <- 5.83  ; label("Concentration at half-maximal hepatic metabolism (mg/L)")                   # Table S4 Km 5.83 (4.96-6.82); V4 posterior mean 5.8324

    # ==================================================================
    # Enterohepatic recirculation and gut absorption -- Table S4 RAT column
    # (hierarchical model V5a), scaled by rat body weight.
    # ==================================================================
    Vehrc   <- 0.527  ; label("Maximum hepatobiliary excretion rate coefficient (mg/h/kg), rat-scaled") # deposited human file 0.527; Table S4 rat column 0.52 (0.43-0.66)
    Kehr    <- 0.017  ; label("Concentration at half-maximal hepatobiliary excretion (mg/L)")           # Table S4 rat column 0.017 (0.0014-0.063)
    kbile   <- 0.25   ; label("Bile-to-gut-lumen transfer rate coefficient (1/h/kg), rat-scaled")       # deposited human file 0.25; Table S4 rat column 0.256 (0.007-0.83)
    kgutabs <- 2.11   ; label("Gut-lumen absorption rate coefficient (1/h/kg), rat-scaled")             # Table S4 rat column 2.11 (1.80-2.48)
    kfeces  <- 0.0187 ; label("Faecal excretion rate constant from gut lumen (1/h)")                    # Table S4 rat column 0.0187 (0.00063-0.0649)

    # ==================================================================
    # Between-subject variability. Results 3.4: "we ran 2000 simulations using
    # individual parameter values drawn randomly around their fitted mean value,
    # with a standard deviation of 1.17 (on a log scale), equivalent to a
    # coefficient of variation (CoV) of 16%". The deposited Monte-Carlo input
    # file writes each as Distrib(<par>, LogNormal, <mean>, 1.17); in GNU MCSim
    # the second LogNormal argument is the GEOMETRIC standard deviation, so the
    # log-scale variance is log(1.17)^2 = 0.02465 and the CV is
    # sqrt(exp(0.02465) - 1) = 15.8% -- which is what reconciles "1.17" with the
    # paper's stated "16%". Applied multiplicatively as parameter * exp(eta).
    # All 18 variances are FIXED at the same assumed value (they were not
    # estimated from inter-individual data), so each is written ~ fixed(...).
    # ==================================================================
    etaKgut_plasma      ~ fixed(0.02465) # Agedynamic_montecarlo.in.R Distrib(Kgut_plasma, LogNormal, 0.622, 1.17); log(1.17)^2
    etaKliver_plasma    ~ fixed(0.02465) # Distrib(Kliver_plasma, LogNormal, 0.651, 1.17)
    etaKkidney_plasma   ~ fixed(0.02465) # Distrib(Kkidney_plasma, LogNormal, 0.671, 1.17)
    etaKfat_plasma      ~ fixed(0.02465) # Distrib(Kfat_plasma, LogNormal, 0.159, 1.17)
    etaKrestbody_plasma ~ fixed(0.02465) # Distrib(Krestbody_plasma, LogNormal, 0.42, 1.17)
    etafu               ~ fixed(0.02465) # Distrib(fu, LogNormal, 0.41, 1.17)
    etaKbp              ~ fixed(0.02465) # Distrib(Kbp, LogNormal, 0.76, 1.17)
    etaQurineC          ~ fixed(0.02465) # Distrib(QurineC, LogNormal, 13.45, 1.17)
    etaTrc              ~ fixed(0.02465) # Distrib(Trc, LogNormal, 1.334, 1.17)
    etaTmc              ~ fixed(0.02465) # Distrib(Tmc, LogNormal, 8.024, 1.17)
    etaKt               ~ fixed(0.02465) # Distrib(Kt, LogNormal, 0.059, 1.17)
    etaVmaxC            ~ fixed(0.02465) # Distrib(VmaxC, LogNormal, 0.472, 1.17)
    etakm               ~ fixed(0.02465) # Distrib(Km, LogNormal, 5.83, 1.17)
    etaVehrc            ~ fixed(0.02465) # Distrib(Vehrc, LogNormal, 0.527, 1.17)
    etaKehr             ~ fixed(0.02465) # Distrib(Kehr, LogNormal, 0.017, 1.17)
    etakbile            ~ fixed(0.02465) # Distrib(kbile, LogNormal, 0.25, 1.17)
    etakgutabs          ~ fixed(0.02465) # Distrib(kgutabs, LogNormal, 2.11, 1.17)
    etakfeces           ~ fixed(0.02465) # Distrib(kfeces, LogNormal, 0.0187, 1.17)

    # ==================================================================
    # Residual error -- Methods 2.3, assumed 10% CV -> fixed().
    # ==================================================================
    propSd <- fixed(0.10) ; label("Proportional residual error (fraction)")                             # Methods 2.3 'coefficient of variation of 10%'
  })

  model({
    # ================================================================
    # 1. Individual drug-specific parameters: mean * exp(eta), lognormal
    #    between-subject variability with geometric SD 1.17.
    # ================================================================
    kgut_p   <- Kgut_plasma      * exp(etaKgut_plasma)
    kliver_p <- Kliver_plasma    * exp(etaKliver_plasma)
    kkid_p   <- Kkidney_plasma   * exp(etaKkidney_plasma)
    kfat_p   <- Kfat_plasma      * exp(etaKfat_plasma)
    krest_p  <- Krestbody_plasma * exp(etaKrestbody_plasma)
    fu_i     <- fu               * exp(etafu)
    kbp_i    <- Kbp              * exp(etaKbp)
    qurinec  <- QurineC          * exp(etaQurineC)
    trc_i    <- Trc              * exp(etaTrc)
    tmc_i    <- Tmc              * exp(etaTmc)
    kt_i     <- Kt               * exp(etaKt)
    vmaxc_i  <- VmaxC            * exp(etaVmaxC)
    km_i     <- km               * exp(etakm)
    vehrc_i  <- Vehrc            * exp(etaVehrc)
    kehr_i   <- Kehr             * exp(etaKehr)
    kbile_i  <- kbile            * exp(etakbile)
    kgutabs_i <- kgutabs         * exp(etakgutabs)
    kfeces_i  <- kfeces          * exp(etakfeces)

    # ================================================================
    # 2. Age- and sex-dependent anthropometry.
    #    Deposited Agedynamics_NFThuman.model.R Dynamics{} block. The
    #    deposited code selects branches through
    #      model = ((sex <= 0) ? 1 : (sex <= 1) ? 2 : 0)
    #    and then tests `model <= 1` / `model <= 2`; branch 1 consumes the
    #    _f constants and branch 2 the _m constants, so branch 1 is FEMALE and
    #    branch 2 is MALE (the deposited header comment states the opposite and
    #    is wrong -- see covariateData$SEXF$notes and the vignette Errata).
    #    HT is in cm and BW in kg; both are 6th-order polynomials in age and
    #    are only validated over the simulated 25-40 year adult range.
    # ================================================================
    ht <- ifelse(
      SEXF > 0.5,
      5.373e+01 + 1.296e+01 * AGE - 5.506e-01 * AGE^2 + 1.113e-02 * AGE^3 -
        1.106e-04 * AGE^4 + 4.697e-07 * AGE^5 - 4.416e-10 * AGE^6,
      5.869e+01 + 1.265e+01 * AGE - 4.665e-01 * AGE^2 + 7.198e-03 * AGE^3 -
        3.224e-05 * AGE^4 - 2.512e-07 * AGE^5 + 2.071e-09 * AGE^6
    )
    bw <- ifelse(
      SEXF > 0.5,
      2.354e+00 + 4.050e+00 * AGE - 3.240e-02 * AGE^2 - 3.057e-03 * AGE^3 +
        9.353e-05 * AGE^4 - 1.022e-06 * AGE^5 + 3.918e-09 * AGE^6,
      3.382e+00 + 2.866e+00 * AGE + 1.694e-01 * AGE^2 - 1.169e-02 * AGE^3 +
        2.577e-04 * AGE^4 - 2.484e-06 * AGE^5 + 8.891e-09 * AGE^6
    )
    # Du Bois body surface area (m^2) from weight in kg and height in cm.
    sa <- 0.007184 * bw^0.425 * ht^0.725

    # ================================================================
    # 3. Age- and sex-dependent organ volumes (L). Rest-of-body closes TOTAL
    #    body volume by subtraction -- note the human file subtracts from BW
    #    itself, not from 0.84 * BW as the rabbit and rat files do.
    #    The deposited file also computes skin volume, bone-marrow volume and a
    #    haematocrit trajectory, none of which enter any ODE or the rest-of-body
    #    closure; they are omitted here.
    # ================================================================
    v_liver <- ifelse(
      SEXF > 0.5,
      0.0017717 - 0.0030113 * AGE + 0.0253455 * bw,
      -0.0143744 - 0.0044728 * AGE + 0.0264591 * bw
    )
    v_kidney <- ifelse(
      SEXF > 0.5,
      0.0458676 - 0.0003957 * AGE + 0.0035115 * bw,
      5.668e-02 - 4.962e-04 * AGE + 3.501e-03 * bw
    )
    v_fat <- ifelse(
      SEXF > 0.5,
      6.132e-01 + 8.475e-02 * AGE + 8.151e-05 * AGE^2 + 1.341e-01 * bw + 2.297e-03 * bw^2,
      1.3054356 + 0.3622685 * AGE - 0.0025165 * AGE^2 + 0.0906119 * bw + 0.0001731 * bw^2
    )
    v_plasma <- ifelse(
      SEXF > 0.5,
      (1423 * sa - 194) / 1000,
      (1587 * sa - 304) / 1000
    )
    v_gut      <- Fgut      * bw
    v_filtrate <- Ffiltrate * bw
    v_rest     <- bw - (v_liver + v_kidney + v_filtrate + v_fat + v_gut + v_plasma)

    # ================================================================
    # 4. Age- and sex-dependent cardiac output and organ plasma flows (L/h).
    #    As in the animal files, plasma flow is set equal to blood flow
    #    (QCplasma = QCC; no haematocrit correction).
    # ================================================================
    qc_plasma <- ifelse(
      SEXF > 0.5,
      5.528076 - 2.834486 * AGE + 0.012591 * AGE^2 + 204.262351 * sa + 19.274290 * sa^2,
      6.48370 - 1.59948 * AGE + 214.68572 * sa
    )
    q_liver <- ifelse(
      SEXF > 0.5,
      (2.590e-01 - 1.042e-03 * AGE + 4.265e-04 * bw) * qc_plasma,
      (2.502e-01 - 1.062e-03 * AGE + 3.439e-04 * bw) * qc_plasma
    )
    q_kidney <- ifelse(
      SEXF > 0.5,
      (1.819e-01 - 1.160e-03 * AGE + 4.873e-04 * bw) * qc_plasma,
      (1.975e-01 - 1.182e-03 * AGE + 3.937e-04 * bw) * qc_plasma
    )
    q_fat <- ifelse(
      SEXF > 0.5,
      (8.298e-02 + 6.850e-04 * AGE - 2.804e-04 * bw) * qc_plasma,
      (5.075e-02 + 4.327e-04 * AGE - 1.401e-04 * bw) * qc_plasma
    )
    q_gut  <- ifelse(SEXF > 0.5, FQgut_f, FQgut_m) * qc_plasma
    q_rest <- qc_plasma - (q_liver + q_kidney + q_fat + q_gut)

    # ================================================================
    # 5. Renal function. `gfr_normal` is the deposited file's CKD == 0 branch:
    #    the age- and sex-predicted NORMAL glomerular filtration rate in mL/min
    #    for this subject's body surface area. It is returned as an output
    #    column so a normal-renal-function cohort can be built by solving once,
    #    reading gfr_normal, and re-solving with CRCL = gfr_normal.
    #    The ODEs are always driven by the supplied CRCL, matching the
    #    deposited file's CKD != 0 branch (Qfilterate = GFR * 0.06).
    #    Factor 0.06 converts mL/min to L/h.
    # ================================================================
    gfr_normal  <- -6.616 * sa^2 + 99.054 * sa - 17.74
    q_filtrate  <- CRCL * 0.06

    # ================================================================
    # 6. Allometric scaling of the fitted animal parameters to this subject's
    #    body weight, then conversion to whole-body / whole-organ terms.
    # ================================================================
    rabbit_scalar <- (BWrabbit / bw)^scaling
    rat_scalar    <- (BWrat    / bw)^scaling

    qurine    <- qurinec  * rabbit_scalar
    tr        <- trc_i    * rabbit_scalar
    tm        <- tmc_i    * rabbit_scalar
    vmax      <- vmaxc_i  * rabbit_scalar * 1000 * v_liver
    vehr      <- vehrc_i  * rat_scalar
    kbile_s   <- kbile_i  * rat_scalar
    kgutabs_s <- kgutabs_i * rat_scalar

    # ================================================================
    # 7. Tissue concentrations (mg/L)
    # ================================================================
    c_gut      <- a_gut          / v_gut
    c_liver    <- a_liver        / v_liver
    c_kidney   <- a_kidney       / v_kidney
    c_filtrate <- a_filtrate     / v_filtrate
    c_fat      <- a_fat          / v_fat
    c_rest     <- a_rest_of_body / v_rest
    Cc         <- a_plasma       / v_plasma

    # ================================================================
    # 8. Unbound, partition-corrected driving concentrations
    # ================================================================
    fu_bp     <- fu_i / kbp_i
    cu_plasma <- Cc       * fu_bp
    cu_gut    <- c_gut    * fu_bp / kgut_p
    cu_liver  <- c_liver  * fu_bp / kliver_p
    cu_kidney <- c_kidney * fu_bp / kkid_p
    cu_fat    <- c_fat    * fu_bp / kfat_p
    cu_rest   <- c_rest   * fu_bp / krest_p

    # ================================================================
    # 9. Renal processes. Secretion is driven by the unbound KIDNEY-TISSUE
    #    concentration, per the deposited code; see the rabbit sibling and the
    #    vignette Errata for the printed-equation conflict.
    # ================================================================
    filtration   <- q_filtrate * cu_kidney
    reabsorption <- tr * c_filtrate
    secretion    <- tm * cu_kidney / (kt_i + cu_kidney)

    # ================================================================
    # 10. Hepatic processes
    # ================================================================
    metabolism     <- vmax * c_liver * fu_i / (km_i   + c_liver * fu_i)
    biliary_efflux <- vehr * c_liver * fu_i / (kehr_i + c_liver * fu_i)

    # ================================================================
    # 11. ODE system -- 13 mass-balance amount states (mg).
    #     Dose oral into a_gut_lumen, IV into a_plasma.
    # ================================================================
    d/dt(a_gut_lumen)     <- -kgutabs_s * a_gut_lumen - kfeces_i * a_gut_lumen + kbile_s * a_bile
    d/dt(a_gut)           <-  kgutabs_s * a_gut_lumen + q_gut * (cu_plasma - cu_gut)
    d/dt(a_liver)         <-  q_liver * cu_plasma + q_gut * cu_gut -
                              (q_liver + q_gut) * cu_liver - metabolism - biliary_efflux
    d/dt(a_hepatic)       <-  metabolism
    d/dt(a_bile)          <-  biliary_efflux - kbile_s * a_bile
    d/dt(a_feces)         <-  kfeces_i * a_gut_lumen
    d/dt(a_kidney)        <-  q_kidney * cu_plasma - q_kidney * cu_kidney -
                              filtration + reabsorption - secretion
    d/dt(a_filtrate)      <-  filtration - q_filtrate * c_filtrate - reabsorption + secretion
    d/dt(a_urine_storage) <-  q_filtrate * c_filtrate - qurine * a_urine_storage
    d/dt(a_fat)           <-  q_fat  * cu_plasma - q_fat  * cu_fat
    d/dt(a_rest_of_body)  <-  q_rest * cu_plasma - q_rest * cu_rest
    d/dt(a_plasma)        <-  (q_liver + q_gut) * cu_liver + q_kidney * cu_kidney +
                              q_fat * cu_fat + q_rest * cu_rest - qc_plasma * cu_plasma
    d/dt(a_urine)         <-  qurine * a_urine_storage

    # ================================================================
    # 12. Diagnostics from the deposited CalcOutputs{} block
    # ================================================================
    mass_balance   <- a_gut_lumen + a_gut + a_liver + a_hepatic + a_bile + a_feces +
                      a_kidney + a_filtrate + a_urine_storage + a_fat +
                      a_rest_of_body + a_plasma + a_urine
    total_renal_cl <- filtration + secretion - reabsorption

    # ================================================================
    # 13. Observation: total plasma concentration (mg/L)
    # ================================================================
    Cc ~ prop(propSd)
  })
}
