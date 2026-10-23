Byun_2024_HL2351_tmdd <- function() {
  description <- paste(
    "Full target-mediated drug disposition (TMDD) population PK model for subcutaneous HL2351",
    "(hIL-1Ra-hyFc, ~97 kDa) in healthy adult Korean men, with FcRn-mediated recycling; case 1",
    "'TMDD' model of Byun 2024, which refits the Ngo 2020 HL2351 structure four ways to compare",
    "TMDD approximations. Drug moves from the injection-site depot (Ka1, lag time) into an",
    "absorption/distribution space where free drug binds FcRn by explicit association /",
    "dissociation (kon1, koff1) into a drug-FcRn complex that recycles to plasma (Krec); free",
    "drug there is degraded (Kdeg) or moves directly to plasma (Ka2). In central, free drug binds",
    "IL1R by explicit kinetics (kon2, koff2) into a complex internalised at Kint, distributes to",
    "one peripheral compartment, is taken back up to the absorption space (Kup), and is",
    "eliminated linearly (CL). The observed quantity is the FREE concentration. Amounts nmol,",
    "concentrations nmol/L."
  )
  reference <- paste(
    "Byun JH, Jeon HS, Yun H-y, Kim JK. Validity conditions of approximations for a",
    "target-mediated drug disposition model: A novel first-order approximation and its",
    "comparison to other approximations. PLoS Comput Biol. 2024;20(4):e1012066.",
    "doi:10.1371/journal.pcbi.1012066. Case 1 TMDD model: Table 4 (estimates) and",
    "S1 Text (a) (NONMEM control stream). The clinical study and base structure are from",
    "Ngo L et al. CPT Pharmacometrics Syst Pharmacol. 2020;9(10):584-595.",
    "doi:10.1002/psp4.12555; see modellib('Ngo_2020_HL2351')."
  )
  vignette <- "Byun_2024_tmdd_approximations"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(analyte = "HL2351", units = "nmol", specimen = "administration site", verified = TRUE),
    abs_site = list(analyte = "HL2351 (free)", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "HL2351 (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "HL2351", units = "nmol", specimen = "tissue", verified = TRUE),
    complex_fcrn = list(
      analyte = "HL2351-FcRn complex",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    complex = list(analyte = "HL2351-IL1R complex", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 40,
    n_studies = 1,
    study_id = "NCT02175056 (Phase I single-ascending-dose, Seoul National University Hospital)",
    age_range = "20-45 years",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy adult Korean men",
    dose_range = "Single SC dose: HL2351 1, 2, 4, 8, or 12 mg/kg (n = 8 per dose group)",
    regions = "Republic of Korea",
    n_observations = 472,
    notes = paste(
      "Byun 2024 case 1 refits the five HL2351 cohorts of the first-in-human study reported by",
      "Ngo 2020; Byun 2024 states 640 data records, which matches the S1 Table dataset (40",
      "subjects; 40 dose records, 128 MDV = 1 sample rows and 472 observations over 672 h,",
      "the same 472 observations as Jung 2023 Model Case 4). Byun 2024 does",
      "not restate the demographics, so those above are taken from Ngo 2020. Doses in the",
      "dataset are in nmol (molecular weight ~97 kDa; 1 mg = approximately 10.31 nmol)."
    )
  )

  ini({
    # Final estimates are Byun 2024 Table 4 (case 1, TMDD column). Parameters that
    # Table 4 does not list (Ka2, kon1, koff1, kon2, koff2, Alag1, and the IIVs
    # on Ka2 and Alag1) are taken from the S1 Text (a) $THETA / $OMEGA, which
    # carry the final estimates re-used as starting values (they agree with
    # Table 4 to within rounding wherever both print a value).
    lka1 <- log(1.22); label("Absorption rate from injection-site depot to the absorption space (1/h)") # Table 4 case 1 TMDD Ka = 1.22; S1 Text (a) THETA(1) 1.22
    lka2 <- log(0.0166); label("Direct rate of free drug from the absorption space to central (1/h)") # S1 Text (a) THETA(4) Ka2 = 0.0166 (not in Table 4)
    lkrec <- log(0.0346); label("Recycling rate of the drug-FcRn complex into central (1/h)") # Table 4 case 1 TMDD Krecycle = 0.0346; S1 Text (a) THETA(5) 0.0345
    lcl <- log(0.210); label("Apparent clearance of free drug from central (L/h)") # Table 4 case 1 TMDD CL = 0.210; S1 Text (a) THETA(7) 0.208
    lvc <- log(11.4); label("Apparent central volume of distribution (L)") # Table 4 case 1 TMDD Vc = 11.4; S1 Text (a) THETA(8) V3
    lq <- log(0.0285); label("Apparent inter-compartmental clearance (L/h)") # Table 4 case 1 TMDD Q = 0.0285; S1 Text (a) THETA(9)
    lv_abs_site <- log(117); label("Apparent volume of the absorption space for FcRn binding (L)") # Table 4 case 1 TMDD Vd = 117; S1 Text (a) THETA(20) V2
    lkon1 <- log(0.249); label("Drug-FcRn association rate constant (L/nmol/h)") # S1 Text (a) THETA(2) KON1 = 0.249 (not in Table 4)
    lkoff1 <- log(0.466); label("Drug-FcRn dissociation rate constant (1/h)") # S1 Text (a) THETA(18) KOFF1 = 0.466 (not in Table 4)
    lkon2 <- log(6.47); label("Drug-IL1R association rate constant (L/nmol/h)") # S1 Text (a) THETA(13) KON2 = 6.47 (not in Table 4)
    lkoff2 <- log(95.4); label("Drug-IL1R dissociation rate constant (1/h)") # S1 Text (a) THETA(19) KOFF2 = 95.4 (not in Table 4)
    ltlag <- log(0.241); label("Subcutaneous absorption lag time (h)") # S1 Text (a) THETA(15) Alag1 = 0.241 (not in Table 4)

    # Held constant at literature values (Table 4 asterisk; S1 Text (a) FIX).
    lkdeg <- fixed(log(0.0264)); label("Degradation rate of free drug in the absorption space (1/h)") # S1 Text (a) THETA(6) 0.0264 FIX; Table 4 prints 0.264 (decimal-point slip, see vignette)
    lkint <- fixed(log(0.206)); label("Internalisation rate of the IL1R-drug complex (1/h)") # Table 4 case 1 Kint = 0.206 (literature); S1 Text (a) THETA(11) Kdeg2
    lkup <- fixed(log(0.00952)); label("Uptake rate of free drug from central into the absorption space (1/h)") # Table 4 case 1 Kuptake = 0.00952 (literature); S1 Text (a) THETA(14)
    lc_il1r_t <- fixed(log(2.23)); label("Total IL1R concentration in central (nmol/L)") # Table 4 case 1 Rtot = 2.23 (literature); S1 Text (a) THETA(12)
    la_fcrn_t <- fixed(log(749)); label("Total FcRn amount in the absorption space (nmol)") # S1 Text (a) THETA(3) FcRn = 749 FIX
    lvp <- fixed(log(5.06)); label("Apparent peripheral volume of distribution (L)") # S1 Text (a) THETA(10) V4 = 5.06 FIX

    # IIV: log-normal; Table 4 prints variances (omega^2).
    etalka1 ~ 0.88 # Table 4 case 1 TMDD IIV Ka = 0.88; S1 Text (a) OMEGA 1 = 0.883
    etalka2 ~ 0.549 # S1 Text (a) OMEGA 2 on Ka2 = 0.549 (not in Table 4)
    etalkrec ~ 0.098 # Table 4 case 1 TMDD IIV Krecycle = 0.098; S1 Text (a) OMEGA 3 = 0.0972
    etalkdeg ~ 0.0851 # Table 4 case 1 TMDD IIV Kdeg = 0.0851; S1 Text (a) OMEGA 4 = 0.0867
    etalcl ~ 0.0497 # Table 4 case 1 TMDD IIV CL = 0.0497; S1 Text (a) OMEGA 5 = 0.051
    etalkup ~ 0.886 # Table 4 case 1 TMDD IIV Kuptake = 0.886; S1 Text (a) OMEGA 6 = 0.874
    etaltlag ~ 0.403 # S1 Text (a) OMEGA 7 on Alag1 = 0.403 (not in Table 4)

    # Residual error: W = SQRT(THETA(16)^2 + THETA(17)^2 * IPRED^2), SIGMA 1 FIX.
    addSd <- 0.188; label("Additive residual error (nmol/L)") # Table 4 case 1 TMDD additive = 0.188; S1 Text (a) THETA(16) 0.187
    propSd <- 0.115; label("Proportional residual error (fraction)") # Table 4 case 1 TMDD proportional = 0.115; S1 Text (a) THETA(17)
  })

  model({
    ka1 <- exp(lka1 + etalka1)
    ka2 <- exp(lka2 + etalka2)
    krec <- exp(lkrec + etalkrec)
    kdeg <- exp(lkdeg + etalkdeg)
    cl <- exp(lcl + etalcl)
    kup <- exp(lkup + etalkup)
    tlag <- exp(ltlag + etaltlag)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    v_abs_site <- exp(lv_abs_site)
    kon1 <- exp(lkon1)
    koff1 <- exp(lkoff1)
    kon2 <- exp(lkon2)
    koff2 <- exp(lkoff2)
    kint <- exp(lkint)
    c_il1r_t <- exp(lc_il1r_t)
    a_fcrn_t <- exp(la_fcrn_t)

    # S1 Text (a) $PK: Kel = CL/V3, Kpt = Q/V3, Ktp = Q/V4.
    kel <- cl / vc
    kcp <- q / vc
    kpc <- q / vp

    # S1 Text (a) $DES: Afree = A(2)/V2, Cfree = A(3)/V3. abs_site and central
    # hold FREE drug; the two complexes are separate states. The binding fluxes
    # are kon * C * (R_total - complex), written as in the control stream.
    c_free_abs <- abs_site / v_abs_site
    c_free <- central / vc

    # States in the S1 Text (a) $MODEL order: DEPOT, ABS, CENTRAL, PERIPH, CPX1, CPX2.
    d/dt(depot) <- -ka1 * depot
    d/dt(abs_site) <- ka1 * depot - (kdeg + ka2) * abs_site - kon1 * c_free_abs * a_fcrn_t +
      (kon1 * c_free_abs + koff1) * complex_fcrn + kup * central
    d/dt(central) <- krec * complex_fcrn + ka2 * abs_site + kpc * peripheral1 - (kel + kcp) * central -
      kon2 * central * c_il1r_t + (kon2 * c_free + koff2) * complex - kup * central
    d/dt(peripheral1) <- kcp * central - kpc * peripheral1
    d/dt(complex_fcrn) <- kon1 * c_free_abs * a_fcrn_t - (kon1 * c_free_abs + koff1 + krec) * complex_fcrn
    d/dt(complex) <- kon2 * central * c_il1r_t - (kon2 * c_free + koff2 + kint) * complex

    # ALAG1 acts on compartment 1, the dosing depot.
    alag(depot) <- tlag

    # S1 Text (a) $ERROR: IPRED = A(3)/V3, the FREE concentration.
    Cc <- c_free
    Cc ~ add(addSd) + prop(propSd)
  })
}
