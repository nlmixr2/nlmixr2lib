Jung_2023_HL2351 <- function() {
  description <- paste(
    "Two-compartment quasi-steady-state target-mediated drug disposition (QSS-TMDD) model for",
    "subcutaneous HL2351 (hIL-1Ra-hyFc, ~97 kDa) coupled with FcRn-mediated recycling, with a",
    "fractal (Kopelman) absorption rate. This is Model Case 4 of Jung 2023, whose contribution",
    "is to replace the constant first-order rate carrying free drug from the absorption site to",
    "plasma by the time-dependent instantaneous rate coefficient rate = q / time^h (Eq. 1),",
    "where h is a heterogeneity exponent bounded in (0, 1). Drug leaves the injection-site",
    "depot at Kinj into an absorption site where free drug equilibrates with FcRn by the QSS",
    "quadratic (dissociation constant a_kss1, total FcRn amount a_fcrn_t), is degraded at Kdeg,",
    "and reaches plasma either directly at the fractal rate or by FcRn-mediated recycling of",
    "the complex. In the central compartment free drug equilibrates with IL1R, is internalised",
    "as the complex, exchanged with one peripheral compartment, taken back up to the absorption",
    "site at Kup, and eliminated linearly. The observed quantity is the FREE concentration.",
    "All amounts and concentrations are nmol and nmol/L."
  )
  reference <- paste(
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304. Model Case 4 (fractal model):",
    "Supplementary Table S4 (estimates) and Code S9 (NONMEM control stream).",
    "The base (non-fractal) model and the clinical study are reported in",
    "Ngo L, Lee J, Lim L, Lim H, Bae KS, Hong T, Bae S, Hong Y. Development of a",
    "Pharmacokinetic Model Describing Neonatal Fc Receptor-Mediated Recycling of HL2351,",
    "a Novel Hybrid Fc-Fused Interleukin-1 Receptor Antagonist, to Optimize Dosage Regimen.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(10):584-595. doi:10.1002/psp4.12552;",
    "see modellib('Ngo_2020_HL2351') for the upstream model."
  )
  vignette <- "Jung_2023_fractal_kinetics"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot       = list(analyte = "HL2351", units = "nmol", specimen = "administration site", verified = TRUE),
    abs_site    = list(analyte = "HL2351", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "HL2351", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "HL2351", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 40,
    n_studies      = 1,
    study_id       = "NCT02175056 (Phase I single-ascending-dose, Seoul National University Hospital)",
    age_range      = "20-45 years",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult Korean men",
    dose_range     = "Single SC dose: HL2351 1, 2, 4, 8, or 12 mg/kg (n = 8 per dose group)",
    regions        = "Republic of Korea",
    n_observations = 472,
    notes          = paste(
      "Jung 2023 Table 2 reports Model Case 4 as 40 subjects and 472 observations, fitted over",
      "672 h. The cohort is the five HL2351 arms of the first-in-human single-ascending-dose",
      "study NCT02175056 reported by Ngo 2020 (the sixth arm, anakinra 100 mg, is Model Case 5).",
      "Jung 2023 does not restate the demographics, so those above are taken from Ngo 2020.",
      "Convert mg dosing with molecular weight ~97 kDa (1 mg HL2351 = approximately 10.31 nmol)."
    )
  )

  ini({
    # Structural parameters. Values are the FINAL fractal-model estimates: Code S9
    # $THETA carries the converged values (they reproduce Supplementary Table S4 to
    # three significant figures), so these are final estimates and not initial ones.
    #
    # PARAMETER-NAME MAPPING. Jung 2023 renumbers the two rate constants that carry
    # drug from the absorption site to plasma relative to Ngo 2020, the upstream
    # paper for this same model and dataset. The names below follow the ROLE-based
    # scheme already used by modellib('Ngo_2020_HL2351') so the two files can be
    # compared line by line; each source-trace comment names the Jung 2023 label, so
    # nothing is lost. The correspondence is:
    #
    #   this file   Jung 2023 (Table S4 / Code S9)   Ngo 2020 (Table 1)   role
    #   ---------   -----------------------------   ------------------   ------------------------------
    #   ka1         Kinj                            Ka1                  depot -> absorption site
    #   ka2         Ka1  (made fractal)             Ka2                  free drug, abs site -> plasma
    #   krec        Ka2                             Krec                 FcRn-recycled complex -> plasma
    #   kdeg        Kdeg                            Kdeg1                degradation at absorption site
    #   kup         Kq                              Kup                  plasma -> absorption site
    #   kint        Kint                            Kdeg2                loss of the IL1R-drug complex
    #
    # That the two papers' estimates agree parameter-for-parameter and in order of
    # magnitude for every one of the fourteen shared quantities is what pins this
    # mapping (e.g. ka2 0.0246 vs Ngo 0.0171; krec 0.0472 vs Ngo 0.0338).
    lka1      <- log(0.72295)   ; label("Absorption rate from injection-site depot to absorption site (1/h)")   # Code S9 $THETA 1 (0,0.72295) "Kinj"; Table S4 fractal Kinj = 0.723 (RSE 33.1%)
    lka2      <- log(0.024555)  ; label("Direct rate of free drug from absorption site to central (1/h)")       # Code S9 $THETA 4 (0,0.024555) "Ka1"; Table S4 fractal Ka1 = 0.0246 (RSE 78.4%)
    lkrec     <- log(0.0472069) ; label("FcRn-mediated recycling rate of the FcRn-drug complex into central (1/h)") # Code S9 $THETA 5 (0,0.0472069) "Ka2"; Table S4 fractal Ka2 = 0.0472 (RSE 49.6%)
    lkdeg     <- log(0.022881)  ; label("Degradation rate of free drug at the absorption site (1/h)")           # Code S9 $THETA 6 (0,0.022881) "Kdeg"; Table S4 fractal Kdeg = 0.0229 (RSE 56%)
    lcl       <- log(0.22751)   ; label("Apparent clearance of free drug from central (CL/F, L/h)")             # Code S9 $THETA 7 (0,0.22751); Table S4 fractal CL = 0.228 (RSE 39.6%)
    lvc       <- log(11.2913)   ; label("Apparent central volume of distribution (Vc/F, L)")                    # Code S9 $THETA 8 (0,11.2913) "V3"; Table S4 fractal Vc = 11.3 (RSE 27.7%)
    lq        <- log(0.0304478) ; label("Apparent inter-compartmental clearance (Q/F, L/h)")                    # Code S9 $THETA 9 (0,0.0304478); Table S4 fractal Q = 0.0304 (RSE 38.1%)
    la_kss1   <- log(426.482)   ; label("Drug-FcRn QSS dissociation constant at the absorption site (nmol)")    # Code S9 $THETA 2 (0,426.482) "KSS1"; Table S4 fractal Kss1 = 426 (RSE 129%)
    la_fcrn_t <- log(685.457)   ; label("Total active amount of FcRn at the absorption site (nmol)")            # Code S9 $THETA 3 (0,685.457) "FcRn"; Table S4 fractal FcRn = 685 (RSE 72.7%)
    lkss2     <- log(11.595)    ; label("Drug-IL1R QSS dissociation constant in central (nmol/L)")              # Code S9 $THETA 13 (0,11.595) "KSS2"; Table S4 fractal Kss2 = 11.6 (RSE 63.9%)
    lc_il1r_t <- log(1.58768)   ; label("Total active concentration of IL1R in central (nmol/L)")               # Code S9 $THETA 12 (0,1.58768) "Rtot"; Table S4 fractal Rtot = 1.59 (RSE 87.5%)
    ltlag     <- log(0.286036)  ; label("Subcutaneous absorption lag time (h)")                                 # Code S9 $THETA 17 (0,0.286036) "Alag"; Table S4 fractal Alag = 0.286 (RSE 32.3%)

    # Held fixed by the authors, carried over from the Ngo 2020 base model. Code S9
    # writes these three with a bare "FIX" flag and no bounds.
    lvp   <- fixed(log(5.06))     ; label("Apparent peripheral volume of distribution (Vp/F, L)")               # Code S9 $THETA 10 "5.06 FIX" "V4"; Table S4 Vp = 5.06 (fixed)
    lkint <- fixed(log(0.206))    ; label("Internalisation rate constant of the IL1R-drug complex (1/h)")       # Code S9 $THETA 11 "0.206 FIX" "Kint"; Table S4 Kint = 0.206 (fixed)
    lkup  <- fixed(log(0.00952))  ; label("Uptake rate of free drug from central back to the absorption site (1/h)") # Code S9 $THETA 14 "0.00952 FIX" "Kup"; Table S4 Kup = 0.00952 (fixed)

    # Fractal-kinetics heterogeneity exponent. Bounded (0, 1) in the control
    # stream; h = 0 recovers the constant-rate (Fick) base model of Code S8.
    h_abs <- 0.276798           ; label("Heterogeneity exponent h of the fractal absorption rate (unitless)")   # Code S9 $THETA 18 (0,0.276798,1); Table S4 fractal h = 0.277 (RSE 60.8%)

    # IIV. Code S9 $OMEGA holds variances of log-normal etas; the CV% printed in
    # Table S4 is recovered as sqrt(exp(omega) - 1). Kup is a FIXED typical value
    # that nonetheless carries an estimated IIV, exactly as in Ngo 2020.
    etalka1  ~ 0.721661         # Code S9 $OMEGA 1 (ETA(1) on Kinj); Table S4 fractal Kinj IIV 103% CV [Shr 11.01%]
    etalka2  ~ 0.642134         # Code S9 $OMEGA 2 (ETA(2) on Ka1);  Table S4 fractal Ka1 IIV 94.9% CV [Shr 15.06%]
    etalkrec ~ 0.124918         # Code S9 $OMEGA 3 (ETA(3) on Ka2);  Table S4 fractal Ka2 IIV 36.5% CV [Shr 33.02%]
    etalkdeg ~ 0.0427415        # Code S9 $OMEGA 4 (ETA(4) on Kdeg); Table S4 fractal Kdeg IIV 20.9% CV [Shr 45.58%]
    etalkup  ~ 1.32123          # Code S9 $OMEGA 5 (ETA(5) on Kq);   Table S4 fractal Kup IIV 166% CV [Shr 30.22%]
    etaltlag ~ 0.329387         # Code S9 $OMEGA 6 (ETA(6) on Alag); Table S4 fractal Alag IIV 62.4% CV [Shr 32.61%]
    etalcl   ~ 0.053764         # Code S9 $OMEGA 7 (ETA(7) on CL);   Table S4 fractal CL IIV 23.5% CV [Shr 16.03%]

    # Residual error. Code S9 $ERROR: W = SQRT(THETA(15)^2 + THETA(16)^2 * IPRED^2)
    # applied to the FREE concentration, with $SIGMA 1 FIX.
    addSd  <- 0.198217          ; label("Additive residual error (nmol/L)")                                     # Code S9 $THETA 15 (0,0.198217) "add"; Table S4 fractal Add-error = 0.198 (RSE 33.7%)
    propSd <- 0.107837          ; label("Proportional residual error (fraction)")                               # Code S9 $THETA 16 (0,0.107837) "pro"; Table S4 fractal Prop-error = 0.108 (RSE 4.81%)
  })

  model({
    # 1. Individual parameters
    ka1      <- exp(lka1  + etalka1)
    ka2      <- exp(lka2  + etalka2)
    krec     <- exp(lkrec + etalkrec)
    kdeg     <- exp(lkdeg + etalkdeg)
    kup      <- exp(lkup  + etalkup)
    cl       <- exp(lcl   + etalcl)
    tlag     <- exp(ltlag + etaltlag)
    vc       <- exp(lvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    kint     <- exp(lkint)
    a_kss1   <- exp(la_kss1)
    a_fcrn_t <- exp(la_fcrn_t)
    kss2     <- exp(lkss2)
    c_il1r_t <- exp(lc_il1r_t)

    # 2. Micro-constants. Code S9 $PK: Kel = CL/V3, Kpt = Q/V3, Ktp = Q/V4.
    kel <- cl / vc
    kcp <- q  / vc
    kpc <- q  / vp

    # 3. Fractal absorption rate. Jung 2023 Eq. (1): rate = q_rate / time^h.
    #    Code S9 $PK: KF = EXP(LOG(Ka1) - H*LOG(TIME)). Case 4 is a
    #    single-ascending-dose study, so absolute solver time equals time after
    #    dose. Unlike Case 2 the authors specify no onset-time regulariser here;
    #    h < 1 makes the singularity at time = 0 integrable, and the absorption
    #    lag keeps the absorption site empty until time > tlag in any case.
    kf <- ka2 / time^h_abs

    # 4. QSS approximations. Code S9 $DES.
    #    (a) At the absorption site the state holds the TOTAL amount (free +
    #        FcRn-bound); free amount is the positive root of the QSS quadratic.
    #        DAA = A(2) - FcRn - KSS1;
    #        Afree = 0.5*(DAA + SQRT(DAA^2 + 4*KSS1*A(2))).
    daa      <- abs_site - a_fcrn_t - a_kss1
    a_free   <- 0.5 * (daa + sqrt(daa * daa + 4 * a_kss1 * abs_site))
    a_fcrn_d <- abs_site - a_free

    #    (b) In central the state holds the TOTAL amount; free concentration is
    #        the positive root of the second QSS quadratic.
    #        Ct = A(3)/V3; D = Ct - Rtot - KSS2;
    #        CP = 0.5*(D + SQRT(D^2 + 4*KSS2*Ct)).
    c_tot  <- central / vc
    dcen   <- c_tot - c_il1r_t - kss2
    c_free <- 0.5 * (dcen + sqrt(dcen * dcen + 4 * kss2 * c_tot))

    # 5. ODE system. Code S9 $MODEL declares the compartments in the order
    #    1 DEPOT (DEFDOSE), 2 ABS, 3 CENTRAL, 4 PERIPH, reproduced here.
    #
    #    The recycling flux is written two equivalent ways in Code S9: DADT(2)
    #    removes Ka2*FcRn*Afree/(KSS1+Afree) while DADT(3) adds Ka2*(A(2)-Afree).
    #    Under the QSS these are the same quantity -- the FcRn-bound amount -- so
    #    the single form krec * a_fcrn_d is used on both sides, which makes the
    #    mass balance exact rather than merely algebraically equal.
    d/dt(depot)       <- -ka1 * depot
    d/dt(abs_site)    <-  ka1 * depot - (kdeg + kf) * a_free - krec * a_fcrn_d +
      kup * c_free * vc
    d/dt(central)     <-  kf * a_free + krec * a_fcrn_d + kpc * peripheral1 -
      (kel + kcp) * c_free * vc -
      kint * c_il1r_t * c_free * vc / (kss2 + c_free) -
      kup * c_free * vc
    d/dt(peripheral1) <-  kcp * c_free * vc - kpc * peripheral1

    # 6. Subcutaneous absorption lag time. Code S9 $PK: Alag1 = THETA(17)*EXP(ETA(6)),
    #    i.e. NONMEM's reserved ALAG1 acting on compartment 1, the dosing depot.
    alag(depot) <- tlag

    # 7. Observation: FREE HL2351 concentration (Code S9 $ERROR: IPRED = Cfree).
    Cc <- c_free
    Cc ~ add(addSd) + prop(propSd)
  })
}
