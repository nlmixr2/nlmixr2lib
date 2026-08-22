Chen_2025_iohexol_creatinine <- function() {
  description <- paste(
    "Joint population pharmacokinetic model of intravenous iohexol and creatinine in 14 healthy adults",
    "(Chen 2025), fit simultaneously to dense plasma and urine data for both analytes. Iohexol follows",
    "three-compartment linear disposition and its clearance IS the glomerular filtration rate (GFR).",
    "Creatinine follows one-compartment disposition driven by two inputs: a zero-order endogenous",
    "generation rate (CGR, replaced during model development by the Cockcroft-Gault expression in age,",
    "total body weight and sex) and first-order absorption of a cooked-beef creatinine load with a lag",
    "time and an estimated bioavailability. Creatinine clearance is the sum of GFR and a net tubular",
    "secretion arm (nCTS), the OCT2/MATE-mediated secretory flux net of tubular reabsorption, which",
    "accounted for 31% of total creatinine clearance in this cohort. Both analytes are assumed to be",
    "solely renally eliminated, so each carries a cumulative urinary-excretion state and a urine output",
    "alongside its plasma output. A piecewise-sine circadian rhythm (14 h daytime rise, 10 h nocturnal",
    "fall) multiplies GFR and nCTS with a single shared pair of amplitudes. All clearances and volumes",
    "are allometrically scaled to a 70 kg reference weight with fixed exponents of 0.75 and 1.",
    sep = " "
  )
  reference <- paste(
    "Chen Z, Dong Q, Dokos C, Boland J, Fuhr U, Taubert M.",
    "A Joint Pharmacometric Model of Iohexol and Creatinine Administered through a Meat Meal",
    "to Assess GFR and Renal OCT2/MATE Activity.",
    "Clin Pharmacol Ther. 2025;118(2):510-520.",
    "doi:10.1002/cpt.3612.",
    "Parameter values are taken from Table 3 and from the deposited NONMEM control stream",
    "in the Supporting Information (CPT-118-510-s001.docx), which is the authoritative source",
    "wherever the two disagree; see the validation vignette Errata.",
    sep = " "
  )
  vignette <- "Chen_2025_iohexol_creatinine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Declared explicitly: buildModelDb()'s auto-detection only recognises the
  # literal names "depot" and "central", so it would miss the creatinine meal
  # depot. Iohexol is dosed intravenously into `central`; the cooked-beef
  # creatinine load enters `depot_creatinine`.
  dosing <- c("central", "depot_creatinine")

  compartmentData <- list(
    central             = list(analyte = "iohexol",    units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1         = list(analyte = "iohexol",    units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral2         = list(analyte = "iohexol",    units = "mg", specimen = "plasma",              verified = TRUE),
    urine               = list(analyte = "iohexol",    units = "mg", specimen = "urine",               verified = TRUE),
    depot_creatinine    = list(analyte = "creatinine", units = "mg", specimen = "administration site", verified = TRUE),
    central_creatinine  = list(analyte = "creatinine", units = "mg", specimen = "plasma",              verified = TRUE),
    urine_creatinine    = list(analyte = "creatinine", units = "mg", specimen = "urine",               verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed in the source study. Used for standard allometric scaling to a 70 kg reference",
        "weight with exponents fixed at 0.75 for GFR, nCTS, Qp1 and Qp2 and at 1 for iohexol Vc, Vp1,",
        "Vp2 and creatinine Vd (Chen 2025 Methods, 'Covariate model'; control-stream TBWonCL / TBWonV).",
        "Also enters the Cockcroft-Gault creatinine generation rate. Cohort mean 78.5 kg",
        "(range 59.1-95.8 kg; Table 2)."
      ),
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the model only through the Cockcroft-Gault creatinine generation rate",
        "CGR = (140 - AGE) * WT / 72 * 0.85^SEXF * 60/100 (Chen 2025 Table 3, CGR row).",
        "Age was screened as a covariate on the PK parameters themselves and was NOT retained",
        "(Discussion: 'Age ... was not identified as significant in this study'). Cohort mean 33 years",
        "(range 23-48; Table 2)."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Two distinct uses: (a) the 0.85 female factor of the Cockcroft-Gault",
        "creatinine generation rate, and (b) a multiplicative effect of 0.628 on nCTS, i.e. net tubular",
        "secretion in females is 62.8% of the male value (Chen 2025 Table 3 row 'SEX on CTS';",
        "Discussion attributes it to higher transporter abundance in males).",
        "POLARITY NOTE: the source data column is inverted relative to this canonical column. Chen 2025",
        "Table 3 footnote a states 'SEX is a categorical covariate of 0 for male and 1 for female', but",
        "the deposited control stream codes both effects as theta^(1 - SEX) -- CG = ... * (0.85**(1-SEX))",
        "and SEXonTS = THETA(15)**(1-SEX). Since Cockcroft-Gault applies the 0.85 factor to females, the",
        "dataset column must have been 1 = male / 0 = female, contradicting the footnote. Table 2's",
        "demographics settle it independently: the tabulated estimated CGR is 78.4 mg/h for males",
        "(= (140-31)*86.2/72*0.6, no 0.85 factor) and 48.0 mg/h for females",
        "(= (140-37)*64.5/72*0.6*0.85 = 47.1). Encoded here on the canonical SEXF polarity as",
        "0.85^SEXF and 0.628^SEXF, which reproduces the paper's numbers; see vignette Errata."
      ),
      source_name        = "SEX"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  # Documented so the paper's covariate screen is preserved without implying
  # these columns are required to run the model.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a body-size covariate (Chen 2025 Methods, 'Covariate model'); not retained. Cohort mean 178 cm (Table 2)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a body-size covariate; not retained. Cohort mean 24.7 kg/m^2 (Table 2)."
    ),
    LBM = list(
      description = "Lean body mass",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a body-size covariate; not retained. Table S2 records the head-to-head test:",
        "replacing total body weight with fat-free mass (models 2 and 3) WORSENED the objective function",
        "by 7.47 and 0.82 points, and estimating the scaling exponents instead of fixing them improved it",
        "by only 0.85, so standard allometry on total body weight was kept. Cohort mean 58.6 kg (Table 2)."
      )
    ),
    ALB = list(
      description = "Plasma albumin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate on nCTS; not retained. The Discussion notes that lower serum albumin has",
        "been associated with higher net tubular secretion in previous studies but that no effect was",
        "detectable here, plausibly because only healthy participants were enrolled.",
        "Cohort mean 45.6 g/L (Table 2). Estimated total body water (Watson equation) was screened as well",
        "and has no canonical covariate column; it is recorded here only in prose."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 14L,
    n_studies      = 2L,
    n_observations = 2475L,
    age_range      = "23-48 years",
    age_median     = "33 years (mean)",
    weight_range   = "59.1-95.8 kg",
    weight_median  = "78.5 kg (mean)",
    sex_female_pct = 35.7,
    race_ethnicity = c(White = 100),
    disease_state  = "Healthy volunteers with normal renal function (mean estimated GFR 102 mL/min/1.73 m^2, range 80-116)",
    dose_range     = paste(
      "Iohexol 259 mg or 3,235 mg as a single intravenous dose;",
      "creatinine administered as 250 g cooked beef (mean creatinine content 401 mg) eaten 25 minutes",
      "after the iohexol dose"
    ),
    regions        = "Germany (single centre, Cologne)",
    renal_function = "Normal; healthy volunteers only. The model was not fit to any renally impaired data, so extrapolation to reduced GFR is untested.",
    notes          = paste(
      "Demographics from Chen 2025 Table 2 (n = 14: 9 male, 5 female; 2 in a pilot study, 12 in the main",
      "study). Mean height 178 cm (163-196), BMI 24.7 kg/m^2 (21.2-28.9), body surface area 1.96 m^2",
      "(1.63-2.22), plasma albumin 45.6 g/L (40.0-50.0), screening plasma creatinine 0.90 mg/dL",
      "(0.60-1.17). The dataset comprised 771 iohexol and 826 creatinine plasma concentrations plus 439",
      "urine measurements for each analyte (2,475 records total); observations with missing data (< 5%)",
      "were discarded. Crossover design with a fasting reference period, a fasting low-dose test period",
      "and a meat period; washout at least 7 days. Participants were healthy Caucasian volunteers and",
      "received approximately 240 mL of water at every urine collection interval, which the Discussion",
      "notes may have kept them in a rehydrated state and so contributed to the relatively high nCTS",
      "fraction. Fit in NONMEM 7.4.0 with FOCE-I; 1,000-sample non-parametric bootstrap."
    )
  )

  ini({
    # ---- Iohexol: three-compartment disposition -----------------------------
    # Typical values at the 70 kg allometric reference weight. Table 3 reports
    # GFR in mL/min; converted to L/h as value * 60 / 1000. The deposited
    # control stream $THETA carries more significant digits than Table 3 and is
    # used here (Table 3's rounded values are given in each comment).
    lcl  <- log(5.22317)   ; label("Iohexol clearance = glomerular filtration rate GFR (L/h)")          # $THETA 7 CL_IOX = 5.22317 L/h; Table 3 'GFR (mL/min)' = 87.0 (RSE 3.4%, 95% CI 81.0-92.7) = 5.220 L/h
    lvc  <- log(8.69091)   ; label("Iohexol central volume Vc (L)")                                     # $THETA 8 V1_IOX = 8.69091; Table 3 'Vc (L)' = 8.69 (RSE 4.7%, 95% CI 7.91-9.54)
    lq   <- log(0.130821)  ; label("Iohexol first inter-compartmental clearance Qp1 (L/h)")             # $THETA 9 Q_IOX = 0.130821; Table 3 'Qp1 (L/h)' = 0.131 (RSE 9.2%, 95% CI 0.107-0.163)
    lvp  <- log(1.15193)   ; label("Iohexol first peripheral volume Vp1 (L)")                           # $THETA 10 V2_IOX = 1.15193; Table 3 'Vp1 (L)' = 1.15 (RSE 5.2%, 95% CI 1.05-1.33)
    lq2  <- log(4.00936)   ; label("Iohexol second inter-compartmental clearance Qp2 (L/h)")            # $THETA 11 Q2_IOX = 4.00936; Table 3 'Qp2 (L/h)' = 4.01 (RSE 8.0%, 95% CI 3.37-4.78)
    lvp2 <- log(4.21713)   ; label("Iohexol second peripheral volume Vp2 (L)")                          # $THETA 12 V3_IOX = 4.21713; Table 3 'Vp2 (L)' = 4.22 (RSE 3.1%, 95% CI 3.93-4.53)

    # ---- Creatinine: one-compartment disposition ---------------------------
    lka_creatinine      <- log(1.70996)  ; label("Creatinine absorption rate constant from the beef meal Ka (1/h)")            # $THETA 1 KA_CRE = 1.70996; Table 3 'Ka (1/h)' = 1.71 (RSE 13.3%, 95% CI 1.35-2.19)
    lcl_tsnet_creatinine <- log(2.38177) ; label("Net creatinine tubular secretion clearance nCTS, male reference (L/h)")      # $THETA 2 CL_SEC = 2.38177 L/h; Table 3 'nCTS (mL/min)' = 39.7 (RSE 10.2%, 95% CI 31.2-46.8) = 2.382 L/h
    lvc_creatinine      <- log(28.943)   ; label("Creatinine volume of distribution Vd (L)")                                   # $THETA 3 V_CRE = 28.943; Table 3 'Vd (L)' = 28.9 (RSE 6.5%, 95% CI 25.6-33.4); Results: 41.3% of total body weight
    lfdepot_creatinine  <- log(0.523045) ; label("Bioavailability of creatinine from the cooked-beef meal F1 (fraction)")      # $THETA 4 F1_beef = 0.523045; Table 3 'F1 (%)' = 52.3 (RSE 5.2%, 95% CI 47.7-58.4)
    ltlag_creatinine    <- log(0.291197) ; label("Lag time before creatinine absorption from the meal (h)")                    # $THETA 6 LAG_CRE = 0.291197; Table 3 'Lag time (h)' = 0.291 (RSE 3.0%, 95% CI 0.277-0.308)

    # Endogenous creatinine generation rate. THETA(5) is a FIXED multiplier of 1
    # on the Cockcroft-Gault expression, which is evaluated in model() from AGE,
    # WT and SEXF; only the between-subject variability around it was estimated.
    lksyn_creatinine <- fixed(log(1))    ; label("Multiplier on the Cockcroft-Gault creatinine generation rate CGR (unitless)")  # $THETA 5 CGR = '1 FIX'; Table 3 CGR row prints the formula rather than an estimate

    # ---- Allometric scaling, reference weight 70 kg -------------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of total body weight on all clearance terms (unitless)")   # Table 3 'TBW on GFR, Qp1, Qp2, and nCTS' = 0.75 FIX; control-stream TBWonCL = (TBW/70)**0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of total body weight on all volume terms (unitless)")      # Methods: exponent 'of 1.0 for iohexol central compartment volume (Vc), creatinine Vd, and iohexol peripheral compartment volumes (Vp1 and Vp2)'; control-stream TBWonV = (TBW/70)**1. Table 3's row label 'TBW on iohexol Vc and creatinine Vd' omits Vp1/Vp2 -- see vignette Errata.

    # ---- Covariate effect ---------------------------------------------------
    e_sexf_cl_tsnet_creatinine <- 0.628122 ; label("Multiplicative effect of female sex on nCTS (unitless)")       # $THETA 15 SEXonTS = 0.628122, entering as SEXonTS = THETA(15)**(1-SEX); Table 3 'SEX on CTS' = 0.627 (RSE 15.0%, 95% CI 0.455-0.857)

    # ---- Circadian rhythm on GFR and nCTS -----------------------------------
    # Signed fractional amplitudes; stored as the positive magnitudes Table 3
    # prints, with the nocturnal minus sign carried in model(). A single pair is
    # shared by GFR and nCTS ("Circadian rhythms of GFR and nCTS were merged").
    cl_circ_famp_day   <- 0.0370398 ; label("Daytime circadian fractional amplitude on GFR and nCTS (unitless)")     # $THETA 13 CIR_DAY = 0.0370398; Table 3 'Circadian rhythm during daytime (%)' = 3.70 (RSE 21.0%, 95% CI 2.12-5.72)
    cl_circ_famp_night <- 0.0842088 ; label("Nighttime circadian fractional amplitude on GFR and nCTS (unitless)")   # $THETA 14 CIR_NIGHT = 0.0842088; Table 3 'Circadian rhythm during nighttime (%)' = 8.42 (RSE 17.7%, 95% CI 5.10-11.1)

    # ---- Inter-individual variability --------------------------------------
    # Exponential IIV, theta_i = theta * exp(eta_i). Variances are the deposited
    # $OMEGA diagonal. Table 3's "Estimate" column disagrees with $OMEGA for
    # three parameters (GFR, creatinine Vd, F1); in every case Table 3's own
    # CV(%) column back-transforms to the $OMEGA value via
    # CV = sqrt(exp(omega^2) - 1), so the $OMEGA values are used. See vignette
    # Errata. Ka, lag time, Qp1 and Qp2 carry '0 FIX' omegas and so have no eta.
    etalcl  ~ 0.0140315   # $OMEGA 3 CL_IOX; CV = sqrt(exp(0.0140315)-1) = 11.9%, matching Table 3's CV column and the Results narrative that GFR IIV fell to 11.8% once allometry was added. Table 3's Estimate cell reads 0.0226, which back-transforms to 15.1% and would be an INCREASE over the pre-covariate 14.6%.
    etalvc  ~ 0.0211092   # $OMEGA 5 V1_IOX; CV = 14.6%, matching Table 3 (Estimate 0.0211, CV 14.6%)
    etalvp  ~ 0.012076    # $OMEGA 10 V2_IOX; CV = 11.0%, matching Table 3 (Estimate 0.0121, CV 11.0%)
    etalvp2 ~ 0.0113318   # $OMEGA 12 V3_IOX; CV = 10.7%, matching Table 3 (Estimate 0.0113, CV 10.7%)

    etalcl_tsnet_creatinine ~ 0.0519899   # $OMEGA 2 CL_SEC; CV = 23.1%, matching Table 3's CV column (Estimate 0.0506)
    etalvc_creatinine       ~ 0.0226214   # $OMEGA 4 V_CRE;  CV = 15.1%, matching Table 3's CV column. Table 3's Estimate cell reads 0.0211, an exact duplicate of the Vc row above, which back-transforms to 14.6% instead.
    etalksyn_creatinine     ~ 0.0160459   # $OMEGA 6 CGR;    CV = 12.7%, matching Table 3 (Estimate 0.0163). Results: replacing the estimated CGR with Cockcroft-Gault cut CGR IIV from 30.4% to 12.8%.
    etalfdepot_creatinine   ~ 0.0107643   # $OMEGA 8 F1_beef; CV = 10.4%, matching Table 3's CV column. Table 3's Estimate cell reads 0.00810, which back-transforms to 9.0%.

    # ---- Residual error -----------------------------------------------------
    # Four independent proportional models. The deposited $SIGMA pairs each
    # proportional term with an additive term fixed at 0, so all four are purely
    # proportional. propSd = sqrt(sigma^2).
    propSd                   <- 0.130781  ; label("Proportional residual error, iohexol plasma concentration (fraction)")   # $SIGMA prop_IOX_plasma = 0.0171037; sqrt = 0.130781. Table 3 RV 'Iohexol / Plasma concentration' Estimate 0.0171 (RSE 18.4%)
    propSd_Aurine            <- 0.248338  ; label("Proportional residual error, iohexol amount excreted in urine (fraction)")   # $SIGMA prop_IOX_urine = 0.061672; sqrt = 0.248338. Table 3 RV 'Iohexol / Excreted amount in urine' Estimate 0.0617 (RSE 27.2%)
    propSd_creatinine     <- 0.0505066 ; label("Proportional residual error, creatinine plasma concentration (fraction)")    # $SIGMA prop_CRE_plasma = 0.00255091; sqrt = 0.0505066. Table 3 RV 'Creatinine / Plasma concentration' Estimate 0.00255 (RSE 7.3%)
    propSd_Aurine_creatinine <- 0.203219  ; label("Proportional residual error, creatinine amount excreted in urine (fraction)")  # $SIGMA prop_CRE_urine = 0.041298; sqrt = 0.203219. Table 3 RV 'Creatinine / Excreted amount in urine' Estimate 0.0413 (RSE 41.4%)
  })

  model({
    # ---- 1. Derived covariate terms ---------------------------------------
    # Standard allometric scaling to a 70 kg reference weight.
    wt_cl <- (WT / 70)^e_wt_cl
    wt_v  <- (WT / 70)^e_wt_vc

    # Circadian rhythm shared by GFR and nCTS. The deposited control stream
    # splits the day into a 14 h daytime window in which the rhythm rises and a
    # 10 h night window in which it falls:
    #   day   (0 <= tod < 14):  CIR = 1 + famp_day   * sin(tod / 14 * pi)
    #   night (14 <= tod < 24): CIR = 1 - famp_night * sin((tod - 14) / 10 * pi)
    # Both branches equal 1 at tod = 0, 14 and 24 because sin(0) = sin(pi) = 0,
    # so the function is continuous and introduces no solver discontinuity.
    # `tod` is time of day, i.e. time reset to 0 every 24 h; this generalises
    # the control stream's DAY = 1 / DAY = 2 latch beyond its 38 h horizon and
    # is numerically identical to it over the study window (see vignette Errata).
    tod     <- t - 24 * floor(t / 24)
    isnight <- (tod >= 14)
    cir     <- 1 +
      (1 - isnight) * cl_circ_famp_day   * sin(tod / 14 * 3.141592653589793) -
      isnight       * cl_circ_famp_night * sin((tod - 14) / 10 * 3.141592653589793)

    # Cockcroft-Gault creatinine generation rate, in mg/h. This is the
    # Cockcroft-Gault creatinine clearance (mL/min) multiplied by the serum
    # creatinine it would be divided by, so serum creatinine cancels; the
    # trailing 60/100 converts mL/min * mg/dL to mg/h. Encoded here on the
    # canonical SEXF polarity (1 = female), which reproduces the 0.85 female
    # factor of Cockcroft-Gault -- see the SEXF covariateData note.
    cg <- (140 - AGE) * WT / 72 * 0.85^SEXF * 60 / 100

    # ---- 2. Individual parameters ------------------------------------------
    # Iohexol clearance IS the glomerular filtration rate.
    cl  <- exp(lcl  + etalcl)  * wt_cl * cir
    vc  <- exp(lvc  + etalvc)  * wt_v
    q   <- exp(lq)             * wt_cl
    vp  <- exp(lvp  + etalvp)  * wt_v
    q2  <- exp(lq2)            * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v

    ka_creatinine <- exp(lka_creatinine)
    vc_creatinine <- exp(lvc_creatinine + etalvc_creatinine) * wt_v
    cl_tsnet_creatinine <- exp(lcl_tsnet_creatinine + etalcl_tsnet_creatinine) *
      wt_cl * cir * e_sexf_cl_tsnet_creatinine^SEXF
    ksyn_creatinine <- exp(lksyn_creatinine + etalksyn_creatinine) * cg

    # Total creatinine clearance = glomerular filtration + net tubular secretion.
    cl_creatinine <- cl + cl_tsnet_creatinine

    # ---- 3. Micro-constants -------------------------------------------------
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2
    kel_creatinine <- cl_creatinine / vc_creatinine

    # ---- 4. ODE system ------------------------------------------------------
    # Iohexol: three-compartment linear disposition with the whole elimination
    # flux accumulating in the urine state (renal elimination assumed complete).
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(urine)       <- kel * central

    # Creatinine: one compartment with a first-order meal input, a zero-order
    # endogenous generation term, and renal elimination into its own urine state.
    d/dt(depot_creatinine)   <- -ka_creatinine * depot_creatinine
    d/dt(central_creatinine) <- ka_creatinine * depot_creatinine +
      ksyn_creatinine - kel_creatinine * central_creatinine
    d/dt(urine_creatinine)   <- kel_creatinine * central_creatinine

    # ---- 5. Dose modifiers and initial conditions ---------------------------
    f(depot_creatinine)    <- exp(lfdepot_creatinine + etalfdepot_creatinine)
    alag(depot_creatinine) <- exp(ltlag_creatinine)

    # Creatinine starts at its endogenous steady state, C0 = CGR / CrCL. The
    # circadian factor is exactly 1 at t = 0, so the steady-state clearance
    # below is the t = 0 value of cl_creatinine (control stream:
    # C0 = CGR/CL_CRE; A_0(2) = C0*V_CRE).
    cl_creatinine_t0 <- exp(lcl + etalcl) * wt_cl +
      exp(lcl_tsnet_creatinine + etalcl_tsnet_creatinine) * wt_cl *
      e_sexf_cl_tsnet_creatinine^SEXF
    central_creatinine(0) <- ksyn_creatinine / cl_creatinine_t0 * vc_creatinine

    # ---- 6. Outputs and residual error --------------------------------------
    Cc                <- central / vc
    Aurine            <- urine
    Cc_creatinine     <- central_creatinine / vc_creatinine
    Aurine_creatinine <- urine_creatinine

    Cc                ~ prop(propSd)
    Aurine            ~ prop(propSd_Aurine)
    Cc_creatinine     ~ prop(propSd_creatinine)
    Aurine_creatinine ~ prop(propSd_Aurine_creatinine)
  })
}
