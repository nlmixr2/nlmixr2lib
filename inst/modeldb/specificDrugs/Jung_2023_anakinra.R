Jung_2023_anakinra <- function() {
  description <- paste(
    "One-compartment quasi-steady-state target-mediated drug disposition (QSS-TMDD) model for",
    "subcutaneous anakinra, with a fractal (Kopelman) absorption rate. This is Model Case 5 of",
    "Jung 2023, whose contribution is to replace the constant first-order absorption rate",
    "constant by the time-dependent instantaneous rate coefficient rate = q / time^h (Eq. 1),",
    "where h is a heterogeneity exponent. Drug leaves the injection-site depot at the fractal",
    "rate kf = ka / time^h into a central compartment holding TOTAL (free + IL1R-bound)",
    "anakinra; free concentration is recovered from the total by the QSS quadratic with",
    "dissociation constant kss and total IL1R concentration rtot. Free drug is cleared linearly",
    "(CL) and the drug-IL1R complex is internalised at kint. The observed quantity is the FREE",
    "anakinra concentration, not the total. All amounts and concentrations are nmol and nmol/L."
  )
  reference <- paste(
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304. Model Case 5 (fractal model):",
    "Supplementary Table S5 (estimates) and Code S11 (NONMEM control stream).",
    "The base (non-fractal) anakinra model and the clinical study are reported in",
    "Ngo L, Oh J, Kim A, Back H-m, Kang W-h, Chae J-w, Yun H-y, Lee H. Development of a",
    "Pharmacokinetic Model Describing Neonatal Fc Receptor-Mediated Recycling of HL2351,",
    "a Novel Hybrid Fc-Fused Interleukin-1 Receptor Antagonist, to Optimize Dosage Regimen.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(10):584-595. doi:10.1002/psp4.12552;",
    "see modellib('Ngo_2020_HL2351') for the companion HL2351 model."
  )
  vignette <- "Jung_2023_fractal_kinetics"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot   = list(analyte = "anakinra", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "anakinra", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 8,
    n_studies      = 1,
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult Korean male volunteers",
    dose_range     = "Single subcutaneous anakinra 100 mg",
    regions        = "Republic of Korea",
    n_observations = 93,
    notes          = paste(
      "Jung 2023 Table 2 reports Model Case 5 as 8 subjects and 93 observations; Jung 2023",
      "Section 2.3 states the concentration was observed to a maximum of 48 h and that one",
      "subcutaneous amount was dosed. The anakinra cohort is the comparator arm (100 mg SC,",
      "n = 8) of the HL2351 first-in-human study NCT02175056 reported by Ngo 2020, in healthy",
      "adult Korean men. Jung 2023 does not restate the demographics, so age and weight ranges",
      "are not recorded here."
    )
  )

  ini({
    # Structural parameters. Code S11 $THETA carries the converged fractal-model
    # values; they reproduce Supplementary Table S5 to three significant figures.
    lcl   <- log(9.02825)  ; label("Apparent clearance of free anakinra CL/F (L/h)")             # Code S11 $THETA 1 (0,9.02825); Table S5 fractal CL = 9.03
    lvc   <- log(53.6953)  ; label("Apparent central volume of distribution Vc/F (L)")           # Code S11 $THETA 2 (0,53.6953); Table S5 fractal Vc = 53.7
    lka   <- log(0.469732) ; label("Absorption rate coefficient Ka from the injection site (1/h)") # Code S11 $THETA 3 (0,0.469732); Table S5 fractal Ka = 0.47
    lrtot <- log(1.67483)  ; label("Total IL1R concentration Rtot (nmol/L)")                     # Code S11 $THETA 6 (0,1.67483); Table S5 fractal Rtot = 1.67
    lkss  <- log(1.391)    ; label("Quasi-steady-state constant for IL1R-anakinra binding Kss (nmol/L)") # Code S11 $THETA 7 (0,1.391); Table S5 fractal Kss = 1.39

    # Held fixed by the authors (carried from the Ngo 2020 base model).
    lkint <- fixed(log(0.206)) ; label("Internalisation rate constant of the IL1R-anakinra complex Kint (1/h)") # Code S11 $THETA 5 "0.206 FIX"; Table S5 Kint = 0.206 (fixed)

    # Fractal-kinetics heterogeneity exponent.
    h_abs <- 0.139228      ; label("Heterogeneity exponent h of the fractal absorption rate (unitless)") # Code S11 $THETA 4 (0,0.139228); Table S5 fractal h = 0.139

    # IIV. Code S11 $OMEGA holds variances; Table S5 CV% = sqrt(exp(omega) - 1).
    # Note h carries a log-normal eta in this model: Code S11 H = THETA(4)*EXP(ETA(4)).
    etalcl   ~ 0.0122006   # Code S11 $OMEGA 1; Table S5 fractal CL IIV 11.1% CV [Shr 0.00%]
    etalvc   ~ 0.0396066   # Code S11 $OMEGA 2; Table S5 fractal Vc IIV 20.1% CV [Shr 0.00%]
    etalka   ~ 0.0705846   # Code S11 $OMEGA 3; Table S5 fractal Ka IIV 27.0% CV [Shr 24.59%]
    etah_abs ~ 1.23662     # Code S11 $OMEGA 4; Table S5 fractal h IIV 157% CV [Shr 9.87%]

    # Residual error. Code S11 $ERROR: W = SQRT(THETA(8)^2 + THETA(9)^2 * IPRED^2)
    # applied to the FREE concentration, with $SIGMA 1 FIX.
    addSd  <- 0.0293577    ; label("Additive residual error (nmol/L)")                            # Code S11 $THETA 8 (0,0.0293577); Table S5 fractal Add-error = 0.0294
    propSd <- 0.112586     ; label("Proportional residual error (fraction)")                      # Code S11 $THETA 9 (0,0.112586); Table S5 fractal Prop-error = 0.113
  })

  model({
    # 1. Individual parameters
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    rtot <- exp(lrtot)
    kss  <- exp(lkss)
    kint <- exp(lkint)

    # Code S11: H = THETA(4) * EXP(ETA(4)) -- a multiplicative log-normal eta on a
    # bare (untransformed) exponent, so h is reproduced as h_abs * exp(eta).
    h <- h_abs * exp(etah_abs)

    # 2. Micro-constants
    kel <- cl / vc

    # 3. Fractal absorption rate. Jung 2023 Eq. (1): rate = q_rate / time^h.
    #    Code S11 $PK: KF = EXP(LOG(KaA) - H*LOG(TIME)). Case 5 is single-dose, so
    #    absolute solver time equals time after dose.
    kf <- ka / time^h

    # 4. QSS approximation. `central` holds TOTAL anakinra amount (free + bound);
    #    the free concentration is the positive root of the QSS quadratic.
    #    Code S11 $DES: Ct = A(2)/VP; D = Ct - Rtot - KSSA;
    #    CP = 0.5*(D + SQRT(D^2 + 4*KSSA*Ct)).
    ctot  <- central / vc
    dqss  <- ctot - rtot - kss
    cfree <- 0.5 * (dqss + sqrt(dqss * dqss + 4 * kss * ctot))

    # 5. ODE system. Code S11 $MODEL / $DES: 1 DEPOT (DEFDOSE), 2 CENTRAL.
    #    The third NONMEM compartment is an AUC accumulator (DADT(3) = CP) that
    #    does not feed back into the system; it is omitted here.
    d/dt(depot)   <- -kf * depot
    d/dt(central) <-  kf * depot - kel * cfree * vc -
      kint * rtot * cfree * vc / (kss + cfree)

    # 6. Observation: FREE anakinra concentration (Code S11 $ERROR: IPRED = Cfree).
    Cc <- cfree
    Cc ~ add(addSd) + prop(propSd)
  })
}
