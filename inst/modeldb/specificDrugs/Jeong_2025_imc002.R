Jeong_2025_imc002 <- function() {
  description <- paste(
    "Semi-mechanistic target-mediated drug disposition (TMDD) model for intravenous IMC-002,",
    "a fully human anti-CD47 IgG4 monoclonal antibody, in patients with advanced solid tumours.",
    "Free drug in the central compartment binds CD47 by full (non-QSS) second-order kinetics to",
    "form a complex that is catabolised, and is also taken up into a peripheral space where it",
    "equilibrates with FcRn by a quasi-steady-state quadratic and is recycled back to plasma.",
    "CD47 turns over by zero-order synthesis and first-order degradation. Every amount and",
    "concentration in the model is molar (umol and umol/L); the observed quantity is the FREE",
    "IMC-002 concentration."
  )
  reference <- paste(
    "Jeong S, Lee SY, Kim SH, Kim HT, Yun H-y, Chae J-w, Lee S.",
    "Model-Informed Optimal Dosing of Anti-CD47 Antibody Using Target-Mediated Drug Disposition Model.",
    "Clin Transl Sci. 2025;18(8):e70321. doi:10.1111/cts.70321.",
    "Structure and fixed constants from Data S1 (NONMEM control stream);",
    "final parameter values from Table 2."
  )
  vignette <- "Jeong_2025_imc002"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    central     = list(analyte = "IMC-002", units = "umol", specimen = "plasma", verified = TRUE),
    target      = list(analyte = "CD47", units = "umol/L", specimen = "plasma", verified = TRUE),
    complex     = list(analyte = "IMC-002-CD47 complex", units = "umol/L", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "IMC-002", units = "umol", specimen = "plasma", verified = TRUE)
  )

  # Jeong 2025 tested weight, sex and age by stepwise covariate modelling
  # (forward p = 0.05, backward p = 0.01) and retained NONE of them, so the final
  # model carries no covariate effects. They are recorded here for provenance.
  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by stepwise covariate modelling on the PK parameters; not retained in the final model (Jeong 2025 sec. 2.2 and 3.2)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by stepwise covariate modelling on the PK parameters; not retained in the final model (Jeong 2025 sec. 2.2 and 3.2)."
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened by stepwise covariate modelling on the PK parameters; not retained in the final model (Jeong 2025 sec. 2.2 and 3.2). Reported in Table 1 as male/female counts (8/4)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    study_id       = "NCT05276310 (phase Ia open-label dose escalation, single centre)",
    age_mean       = "58.8 years (SD 8.0)",
    weight_mean    = "64.5 kg (SD 10.2)",
    sex_female_pct = 33.3,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Adults with advanced solid tumours who had failed standard therapy: hepatocellular carcinoma (9), breast cancer (2), gallbladder cancer (1)",
    dose_range     = "5, 10, 20 or 30 mg/kg IMC-002 as a 3 h intravenous infusion every 2 weeks; 3 patients per dose level",
    regions        = "Republic of Korea",
    n_observations = 213,
    notes          = paste(
      "Demographics from Jeong 2025 Table 1; the per-cohort sampling schedule is Table S1.",
      "Dosing is prescribed in mg/kg but the model is parameterised in molar units, so a",
      "molecular weight is needed to convert. IMC-002 is a fully human IgG4 monoclonal",
      "antibody; Jeong 2025 does not report its molecular weight. The validation vignette",
      "uses 146000 g/mol, a representative human IgG4 value, and flags it as an assumption."
    )
  )

  ini({
    # ---- Structural parameters. Final estimates are Jeong 2025 Table 2 ----------------
    # ("Population estimate" column). Data S1 $THETA carries the same values to the
    # precision it prints (Vc 4.19, CL 0.0115, Krec 0.0142, Kel_dR 0.0188 identical;
    # Ksyn 0.00273 vs 0.0027, Vp 85 vs 85.1, Kup 0.00968 vs 0.0097). Where the two
    # differ, the published table wins and the control stream supplies the structure.
    lvc   <- log(4.19)   ; label("Central volume of distribution (Vc, L)")                                          # Table 2 row "Vc (L)" 4.19 (RSE 4%); Data S1 $THETA 1
    lcl   <- log(0.0115) ; label("Linear clearance of free IMC-002 from the central compartment (CL, L/h)")         # Table 2 row "CL (L/h)" 0.0115 (RSE 12%); Data S1 $THETA 6
    lvp   <- log(85.1)   ; label("Peripheral (FcRn recycling space) volume of distribution (Vp, L)")                # Table 2 row "Vp (L)" 85.1 (RSE 18%); Data S1 $THETA 7
    lkup  <- log(0.0097) ; label("First-order uptake rate constant from central into the peripheral space (Kup, 1/h)") # Table 2 row "Kup" 0.0097 (RSE 5%); Data S1 $THETA 8

    # ---- CD47 target turnover --------------------------------------------------------
    lksyn <- log(0.0027) ; label("Zero-order synthesis rate constant of the CD47 receptor (Ksyn, umol/L/h)")        # Table 2 row "Ksyn (uM/h)" 0.0027 (RSE 8%); Data S1 $THETA 2
    lkdeg <- fixed(log(0.0213)) ; label("First-order degradation rate constant of the CD47 receptor (Kdeg, 1/h), from the CD47 protein half-life reported by Du 2023") # Table 2 row "Kdeg (1/h)" 0.0213 (fixed); sec. 3.2 cites ref. 14; Data S1 $THETA 3 FIX

    # ---- IMC-002 / CD47 binding (full second-order TMDD) -----------------------------
    # The control stream fixes KD and Koff from in-house experiments and derives the
    # association rate constant as Kon = Koff / KD,CD47 inside $PK.
    lkd   <- fixed(log(0.046)) ; label("Equilibrium dissociation constant of IMC-002 for CD47 (KD,CD47, umol/L), from in-house binding experiments") # Table 2 row "KD,CD47 (uM)" 0.046 (fixed); sec. 3.2; Data S1 $THETA 4 FIX
    lk2   <- fixed(log(90.4))  ; label("Dissociation rate constant of the IMC-002-CD47 complex (Koff, 1/h), from in-house binding experiments")      # Table 2 row "Koff (1/h)" 90.4 (fixed); sec. 3.2; Data S1 $THETA 5 FIX
    lkint <- log(0.0188) ; label("Elimination (catabolism) rate constant of the IMC-002-CD47 complex (Kel,D-R complex, 1/h)")                        # Table 2 row "Kel,D-R complex (1/h)" 0.0188 (RSE 13%); Data S1 $THETA 12

    # ---- IMC-002 / FcRn binding in the peripheral space (quasi-steady state) ---------
    lkss     <- fixed(log(0.117)) ; label("Quasi-steady-state dissociation constant of IMC-002 for FcRn (KD,FcRn, umol/L)")   # Table 2 row "KD,FcRn (uM)" 0.117 (fixed); Data S1 $THETA 9 FIX ("KSS1")
    lc_fcrn_t <- fixed(log(0.291)) ; label("Total FcRn concentration in the peripheral space (Rtot,FcRn, umol/L), literature value from Li 2018 human in-vitro data") # Table 2 row "FcRn (uM)" 0.291 (fixed); sec. 3.2 cites ref. 15; Data S1 $THETA 10 FIX
    lkrec    <- log(0.0142) ; label("Recycling rate constant of the IMC-002-FcRn complex back into central (Krec, 1/h)")      # Table 2 row "Krec (1/h)" 0.0142 (RSE 11%); Data S1 $THETA 11

    # ---- IIV. Log-normal on Vc and CL, correlated. -----------------------------------
    # Table 2 reports the two variances (omega-Vc 0.037 RSE 35%, omega-CL 0.319 RSE 45%)
    # but no covariance. Data S1 $OMEGA BLOCK(2) is the only source for the off-diagonal
    # (0.0199, implying a correlation of 0.183); its diagonals reproduce Table 2 exactly.
    etalvc + etalcl ~ c(0.037,
                        0.0199, 0.319)  # Table 2 "Interindividual variability" rows; off-diagonal from Data S1 $OMEGA BLOCK(2)

    # ---- Residual error --------------------------------------------------------------
    # Data S1 $ERROR is additive on the LOG scale -- IPRED = LOG(A(1)/VC),
    # Y = IPRED + W*EPS(1) with W = THETA(13) and $SIGMA 1 FIX -- which is a log-normal
    # residual, encoded here as lnorm(). Table 2 labels the row "Proportional error";
    # at this magnitude the two forms are numerically close (11.9 percent), but the
    # control stream is what was fitted. The Table 2 "Residual variability 1*" row is
    # the fixed $SIGMA scaffolding, not a separate parameter.
    expSd <- 0.119 ; label("Log-scale residual standard deviation for free IMC-002 (log units)")  # Table 2 row "Proportional error" 0.119 (RSE 9%); Data S1 $THETA 13
  })

  model({
    # 1. Individual parameters
    vc       <- exp(lvc + etalvc)
    cl       <- exp(lcl + etalcl)
    vp       <- exp(lvp)
    kup      <- exp(lkup)
    ksyn     <- exp(lksyn)
    kdeg     <- exp(lkdeg)
    kd       <- exp(lkd)
    k2       <- exp(lk2)
    kint     <- exp(lkint)
    kss      <- exp(lkss)
    c_fcrn_t <- exp(lc_fcrn_t)
    krec     <- exp(lkrec)

    # 2. Micro-constants. Data S1 $PK: Kel = CL/VC and Kon = Koff/KCD47.
    kel   <- cl / vc
    k1    <- k2 / kd
    rbase <- ksyn / kdeg   # Data S1 $PK: BASE = Ksyn/Kdeg, the baseline CD47 concentration

    # 3. Quasi-steady-state IMC-002 / FcRn binding in the peripheral space.
    #    Data S1 $DES:
    #      DAA   = A(4)/VP - FcRn - KSS1
    #      Dfree = 0.5*(DAA + SQRT(DAA**2 + 4*KSS1*A(4)/VP))
    #      Dcpx  = FcRn*Dfree/(KSS1 + Dfree)
    #    peripheral1 holds the TOTAL (free + FcRn-bound) amount; ctot_p is its
    #    concentration, cfree_p the free concentration and ccpx_p the bound
    #    concentration, which is capped by the total FcRn pool c_fcrn_t.
    ctot_p  <- peripheral1 / vp
    daa     <- ctot_p - c_fcrn_t - kss
    cfree_p <- 0.5 * (daa + sqrt(daa * daa + 4 * kss * ctot_p))
    ccpx_p  <- c_fcrn_t * cfree_p / (kss + cfree_p)

    # 4. ODE system. Data S1 $DES / paper sec. 3.2, reproduced term for term.
    #    NOTE ON STATE UNITS: `central` and `peripheral1` carry AMOUNTS (umol) while
    #    `target` and `complex` carry CONCENTRATIONS (umol/L) referred to Vc. That is
    #    the published parameterisation, and it is what makes the volume factors below
    #    asymmetric between the drug and receptor equations.
    d/dt(central) <- -kup * central - kel * central -
      k1 * central * target + k2 * complex * vc +
      krec * ccpx_p * vp
    d/dt(target) <- ksyn - kdeg * target -
      k1 * (central / vc) * target + k2 * complex
    #    The complex elimination term is written as Kel,D-R complex * A(3)/Vc in BOTH
    #    the printed equation and Data S1 $DES, i.e. the published estimate 0.0188
    #    acts on the complex concentration divided by Vc. Reproduced verbatim: the
    #    estimate was obtained under this form, so rescaling it would change the model.
    d/dt(complex) <- k1 * (central / vc) * target - k2 * complex - kint * complex / vc
    d/dt(peripheral1) <- kup * central - krec * ccpx_p * vp

    # 5. Initial conditions. Data S1 $PK: A_0(2) = BASE, A_0(3) = 0, A_0(4) = 0;
    #    paper sec. 3.2: AR_CD47(t=0) = BASE, ADR_CD47(t=0) = 0, Ap(t=0) = 0.
    target(0) <- rbase

    # 6. Observation: FREE IMC-002 concentration (Data S1 $ERROR: Cf = LOG(A(1)/VC)).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
