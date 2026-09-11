Lv_2024_delafloxacin <- function() {
  description <- "Three-compartment population PK model for intravenous delafloxacin in 58 healthy Chinese adults (phase I trial CTR20213308; single ascending doses of 150, 300 and 600 mg and multiple doses of 300 mg q12h, each given as a 1-hour infusion). Elimination from the central compartment is the sum of a linear clearance CL and a parallel saturable Michaelis-Menten pathway whose maximum rate VM and half-saturation constant KM were fixed at 40 mg/h and 5 ug/mL. Body weight enters CL and the central volume V1 as power functions centred on the cohort median 61.9 kg; no other covariate was retained. Observation is plasma delafloxacin concentration (ug/mL) with a combined additive plus proportional residual error."
  reference <- paste(
    "Lv JX, Huang YH, Kafauit F, Wang YH, Su C, Ma JH, Xu Y, Huang CC, Zhang Q, Su YW. (2024).",
    "Pharmacokinetics and pharmacodynamics of intravenous delafloxacin in healthy subjects:",
    "model-based dose optimization.",
    "Antimicrobial Agents and Chemotherapy 68(7):e00428-24.",
    "doi:10.1128/aac.00428-24.",
    sep = " "
  )
  vignette <- "Lv_2024_delafloxacin"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Lv 2024 Fig. 2 (schematic of the final
  # PopPK model) and the $MODEL / $DES blocks of the final NONMEM control
  # stream reproduced in the supplement ("The code of final delafloxacin
  # model"), which declare COMP=(CENTRAL), COMP=(PERIPH1), COMP=(PERIPH2) with
  # S1 = V1, S2 = V2, S3 = V3 so every state is an amount in mg.
  compartmentData <- list(
    central     = list(analyte = "delafloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "delafloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "delafloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained by the stepwise covariate model (Lv 2024 Results",
        "'Delafloxacin PopPK modeling'). The Methods 'Covariate analysis' quote forward",
        "P < 0.05 and backward P < 0.01; the SCM configuration in the supplement uses",
        "p_forward = 0.05 but p_backward = 0.001. Either way WT was the sole survivor of",
        "the candidate set CL ~ WT/eGFR/LDH/NEUT/PTA/SEX, V1 ~ WT/eGFR/SEX,",
        "V2 ~ HDL/ALT/NEUT/SEX.",
        "Entered as a power function of weight normalised to the cohort median, per",
        "Eq. 3 of the paper. The normaliser is 61.9 kg, read verbatim from the final",
        "NONMEM control stream in the supplement",
        "(CLWT = ((WT/61.9)**THETA(11)); V1WT = ((WT/61.9)**THETA(12)));",
        "this is the overall median weight of the 58-subject PK concentration set in",
        "Table 1, and is the same reference the paper rounds to '62 kg' in the Fig. 4",
        "forest plot. Applied to CL and to the central volume V1 only -- V2, V3, Q2, Q3,",
        "VM and KM carry no covariate. Cohort range 45.0-81.2 kg.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 58L,
    n_studies      = 1L,
    age_range      = "18-43 years",
    age_median     = "28 years",
    weight_range   = "45.0-81.2 kg",
    weight_median  = "61.9 kg",
    sex_female_pct = 46.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy volunteers. Inclusion required men weighing >= 50 kg and women >= 45 kg, aged 18-75 years, in good overall health by medical history, physical examination, vital signs, 12-lead ECG and laboratory tests (Lv 2024 Methods 'Subjects'; exclusion criteria in the supplement). No subject had renal or hepatic impairment: median renal creatinine clearance 126 mL/min (range 85.0-173, Cockcroft-Gault) and median eGFR 135 mL/min/1.73 m^2 (range 85.7-215, modified MDRD for Chinese).",
    renal_function = "Normal; the cohort contained no subjects with renal insufficiency, which is why renal function was not tested as a covariate on CL (Lv 2024 Discussion).",
    dose_range     = "Single ascending intravenous doses of 150 mg (N = 12), 300 mg (N = 23) and 600 mg (N = 11), plus a multiple-dose arm of 300 mg q12h (N = 12) dosed on the evening of day 1, morning and evening of days 2-5 and the morning of day 6. Every dose was a 1-hour intravenous infusion. The 300 mg single-dose group additionally served a 2-period cross-over bioequivalence comparison against Baxdela.",
    regions        = "China (Sir Run Run Hospital, Nanjing Medical University); trial registration CTR20213308.",
    notes          = "Demographics from Lv 2024 Table 1 (main paper) and Table S1 (supplement). 60 subjects were enrolled; one 300 mg single-dose subject withdrew for concomitant medication and one 600 mg subject withdrew for a vasovagal reaction, leaving 58 subjects in the PK concentration set that the model was fit to. Below-quantification-limit records were fewer than 10% of observations and were discarded (M1 method); the assay LLOQ was 0.04 ug/mL."
  )

  ini({
    # ---- Structural PK parameters ----
    # Lv 2024 Table 2, "Delafloxacin final model" / "Final estimation" column. The
    # supplement's control stream lists slightly different $THETA values (4.54924,
    # 7.3707, 16.9656, 17.8772, 25.7329, 0.944349) -- those are the initial estimates
    # of the final run, not its results: the control stream's V2 of 16.97 sits outside
    # the Table 2 bootstrap 95% CI for V2 (14.27-15.66). Table 2 is used throughout.
    # NONMEM V1/V2/V3 map to the canonical vc/vp/vp2; Q2/Q3 map to q/q2.
    lcl  <- log(4.54); label("CL: linear clearance from the central compartment (L/h)")                    # Table 2: CL = 4.54 L/h (RSE 3.6%; bootstrap 4.56 [4.24-4.83])
    lvc  <- log(7.36); label("V1: central volume of distribution (L)")                                     # Table 2: V1 = 7.36 L (RSE 3.3%; bootstrap 7.34 [6.94-7.78])
    lvp  <- log(15.0); label("V2: first peripheral volume of distribution (L)")                            # Table 2: V2 = 15.0 L (RSE 2.7%; bootstrap 14.97 [14.27-15.66])
    lvp2 <- log(18.1); label("V3: second peripheral volume of distribution (L)")                           # Table 2: V3 = 18.1 L (RSE 6.5%; bootstrap 17.80 [12.39-23.75])
    lq   <- log(25.8); label("Q2: intercompartmental clearance central<->peripheral1 (L/h)")               # Table 2: Q2 = 25.8 L/h (RSE 3.1%; bootstrap 25.87 [24.40-27.20])
    lq2  <- log(0.96); label("Q3: intercompartmental clearance central<->peripheral2 (L/h)")               # Table 2: Q3 = 0.96 L/h (RSE 6.1%; bootstrap 0.96 [0.84-1.07])

    # Parallel saturable (Michaelis-Menten) elimination from the central compartment.
    # Both constants were held fixed during estimation (Table 2 records "FIX" and no
    # bootstrap interval for each). VM is a maximum elimination RATE: the $DES block of
    # the supplement's control stream removes VM*A(1)/(KM*V1 + A(1)) from the central
    # compartment, i.e. VM * C1 / (KM + C1) with A(1) in mg, so VM carries mg/h. The
    # "L/h" printed against VM in Table 2 is a table artefact -- the main-text Eq. 5
    # (CLN = Vmax/(KM + C1)) requires an amount-per-time numerator for CLN to come out
    # in L/h, and the Fig. 2 caption's alternative "Vmax*KM/(KM + C1)" contradicts both
    # Eq. 5 and the control stream (see the vignette Assumptions and deviations).
    lvmax <- fixed(log(40)); label("VM: maximum rate of the saturable elimination pathway (mg/h)")         # Table 2: VM = 40, FIX
    lkm   <- fixed(log(5));  label("KM: central concentration at half-maximal saturable elimination (ug/mL)") # Table 2: KM = 5 ug/mL, FIX

    # ---- Covariate effects ----
    # Power model on weight normalised to the cohort median 61.9 kg (Eq. 3; normaliser
    # from the supplement control stream). Both exponents were estimated.
    e_wt_cl <- 1.13; label("Power exponent of (WT / 61.9 kg) on CL (unitless)")                            # Table 2: "The effect of weight on CL" = 1.13 (RSE 19.8%; bootstrap 1.13 [0.76-1.50])
    e_wt_vc <- 1.38; label("Power exponent of (WT / 61.9 kg) on V1 (unitless)")                            # Table 2: "The effect of weight on V1" = 1.38 (RSE 15.0%; bootstrap 1.38 [1.01-1.74])

    # ---- Inter-individual variability ----
    # Exponential IIV (Eq. 1: Pi = theta * exp(eta_i), eta ~ N(0, omega^2)). The values
    # below are the NONMEM $OMEGA diagonal, i.e. VARIANCES on the log scale -- the
    # "(%CV)" in the Table 2 sub-header is a labelling artefact, as the supplement's
    # control stream places the same quantities in an unadorned $OMEGA block. The
    # corresponding coefficients of variation are sqrt(exp(omega^2) - 1):
    #   CL 0.066 -> 26.1%   V1 0.041 -> 20.5%   V2 0.027 -> 16.6%   Q3 0.1 -> 32.4%
    # The control stream fixes IIV on V3, Q2, VM and KM to zero, so no eta is carried
    # for those; the omega matrix is diagonal (no $OMEGA BLOCK).
    etalcl ~ 0.066            # Table 2: IIV_CL = 0.066 (RSE 17%, shrinkage 2.4%; bootstrap 0.063 [0.048-0.084])
    etalvc ~ 0.041            # Table 2: IIV_V1 = 0.041 (RSE 20.9%, shrinkage 9.2%; bootstrap 0.039 [0.026-0.084])
    etalvp ~ 0.027            # Table 2: IIV_V2 = 0.027 (RSE 19.6%, shrinkage 8%; bootstrap 0.027 [0.018-0.036])
    etalq2 ~ fixed(0.1)       # Table 2: IIV_Q3 = 0.1, FIX (matches '0.1 FIX ; IIV_Q3' in the control stream)

    # ---- Residual error ----
    # Eq. 2 combined error. The control stream's $ERROR block is
    #   W = SQRT(THETA(9)**2 * IPRED**2 + THETA(10)**2); Y = IPRED + W*EPS(1); $SIGMA 1 FIX
    # i.e. SD = sqrt(addSd^2 + (propSd * f)^2), which is exactly nlmixr2's default
    # combined2 form for add() + prop(). Both quantities are standard deviations, not
    # variances, and the "(CV%)" in the Table 2 sub-header is again a labelling artefact
    # -- an additive term cannot be a CV.
    propSd <- 0.073; label("Proportional residual error (fraction of the predicted concentration)")        # Table 2: prop.err = 0.073 (RSE 9.8%; bootstrap 0.073 [0.06-0.08])
    addSd  <- 0.1;   label("Additive residual error (ug/mL)")                                              # Table 2: add.err = 0.1 (RSE 13.9%; bootstrap 0.1 [0.07-0.13])
  })

  model({
    # ---- Individual PK parameters ----
    # Weight enters CL and V1 only, as a power of WT/61.9 (Eq. 3 and the supplement's
    # CLCOV / V1COV assignments).
    cl   <- exp(lcl + etalcl) * (WT / 61.9)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 61.9)^e_wt_vc
    vp   <- exp(lvp + etalvp)
    vp2  <- exp(lvp2)
    q    <- exp(lq)
    q2   <- exp(lq2 + etalq2)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # ---- Micro-constants (control stream: K10 = CL/V1, K12 = Q2/V1, K21 = Q2/V2,
    #      K13 = Q3/V1, K31 = Q3/V3) ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Central-compartment concentration driving the saturable pathway (ug/mL = mg/L).
    Cc <- central / vc

    # ---- ODE system (supplement $DES) ----
    # DADT(1) = K21*A(2) + K31*A(3) - (K12 + K13 + K10)*A(1) - VM*A(1)/(KM*V1 + A(1))
    # The saturable term is written here in its equivalent concentration form
    # VM * Cc / (KM + Cc), which is what Eq. 5 (CLN = VM/(KM + C1)) multiplies by Cc.
    d/dt(central)     <- k21 * peripheral1 + k31 * peripheral2 -
      (k12 + k13 + kel) * central - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ---- Observation and error model ----
    Cc ~ add(addSd) + prop(propSd)
  })
}
