Byun_2024_HL2351_mtmdd <- function() {
  description <- paste(
    "Michaelis-Menten (standard quasi-steady-state, mTMDD) approximation of the target-mediated",
    "drug disposition population PK model for subcutaneous HL2351 (hIL-1Ra-hyFc, ~97 kDa) in",
    "healthy adult Korean men, with FcRn-mediated recycling; case 1 'mTMDD' model of Byun 2024.",
    "Same structure as modellib('Byun_2024_HL2351_tmdd') but both bindings are replaced by",
    "Michaelis-Menten (sQSSA) terms on the TOTAL amounts: the FcRn-bound amount in the",
    "absorption space is FcRn_tot * A / (KSS1 + A) and the IL1R-bound concentration in central",
    "is Rtot * C / (KSS2 + C), with C the total central concentration treated as free. Linear",
    "clearance, uptake, and peripheral distribution act on the total central amount. The",
    "observed quantity is the total central concentration. Amounts nmol, concentrations nmol/L."
  )
  reference <- paste(
    "Byun JH, Jeon HS, Yun H-y, Kim JK. Validity conditions of approximations for a",
    "target-mediated drug disposition model: A novel first-order approximation and its",
    "comparison to other approximations. PLoS Comput Biol. 2024;20(4):e1012066.",
    "doi:10.1371/journal.pcbi.1012066. Case 1 mTMDD model: Table 4 (estimates) and",
    "S1 Text (b) (NONMEM control stream). The clinical study and base structure are from",
    "Ngo L et al. CPT Pharmacometrics Syst Pharmacol. 2020;9(10):584-595.",
    "doi:10.1002/psp4.12555; see modellib('Ngo_2020_HL2351')."
  )
  vignette <- "Byun_2024_tmdd_approximations"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(analyte = "HL2351", units = "nmol", specimen = "administration site", verified = TRUE),
    abs_site = list(
      analyte = "HL2351 (total: free + FcRn-bound)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "HL2351 (total: free + IL1R-bound)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "HL2351", units = "nmol", specimen = "tissue", verified = TRUE)
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
      "subjects; 40 dose records, 128 MDV = 1 sample rows and 472 observations over 672 h).",
      "Byun 2024 does not restate the demographics, so those above are taken from Ngo 2020.",
      "Doses in the dataset are in nmol (molecular weight ~97 kDa; 1 mg = approximately",
      "10.31 nmol)."
    )
  )

  ini({
    # Final estimates are Byun 2024 Table 4 (case 1, mTMDD column). Parameters that
    # Table 4 does not list (Ka2, Alag1, and the IIVs on Ka2 and Alag1) are taken
    # from the S1 Text (b) $THETA / $OMEGA, which carry the final estimates re-used
    # as starting values.
    lka1 <- log(1.17); label("Absorption rate from injection-site depot to the absorption space (1/h)") # Table 4 case 1 mTMDD Ka = 1.17; S1 Text (b) THETA(1)
    lka2 <- log(0.018); label("Direct rate of drug from the absorption space to central (1/h)") # S1 Text (b) THETA(4) Ka2 = 0.018 (not in Table 4)
    lkrec <- log(0.0104); label("Recycling rate of the drug-FcRn complex into central (1/h)") # Table 4 case 1 mTMDD Krecycle = 0.0104; S1 Text (b) THETA(5)
    lcl <- log(0.209); label("Apparent clearance from central (L/h)") # Table 4 case 1 mTMDD CL = 0.209; S1 Text (b) THETA(7) 0.21
    lvc <- log(11.6); label("Apparent central volume of distribution (L)") # Table 4 case 1 mTMDD Vc = 11.6; S1 Text (b) THETA(8) 11.7
    lq <- log(0.0145); label("Apparent inter-compartmental clearance (L/h)") # Table 4 case 1 mTMDD Q = 0.0145; S1 Text (b) THETA(9)
    la_kss1 <- log(5.97); label("Drug-FcRn Michaelis-Menten constant in the absorption space (nmol)") # Table 4 case 1 mTMDD KSS1 = 5.97; S1 Text (b) THETA(2) 5.99
    lkss2 <- log(67.8); label("Drug-IL1R Michaelis-Menten constant in central (nmol/L)") # Table 4 case 1 mTMDD KSS2 = 67.8; S1 Text (b) THETA(13)
    ltlag <- log(0.407); label("Subcutaneous absorption lag time (h)") # S1 Text (b) THETA(15) Alag1 = 0.407 (not in Table 4)

    # Held constant at literature values (Table 4 asterisk; S1 Text (b) FIX).
    lkdeg <- fixed(log(0.0264)); label("Degradation rate of drug in the absorption space (1/h)") # S1 Text (b) THETA(6) 0.0264 FIX; Table 4 prints 0.264 (decimal-point slip, see vignette)
    lkint <- fixed(log(0.206)); label("Internalisation rate of the IL1R-drug complex (1/h)") # Table 4 case 1 Kint = 0.206 (literature); S1 Text (b) THETA(11) Kdeg2
    lkup <- fixed(log(0.00952)); label("Uptake rate from central into the absorption space (1/h)") # Table 4 case 1 Kuptake = 0.00952 (literature); S1 Text (b) THETA(14)
    lc_il1r_t <- fixed(log(2.23)); label("Total IL1R concentration in central (nmol/L)") # Table 4 case 1 Rtot = 2.23 (literature); S1 Text (b) THETA(12)
    la_fcrn_t <- fixed(log(749)); label("Total FcRn amount in the absorption space (nmol)") # S1 Text (b) THETA(3) FcRn = 749 FIX
    lvp <- fixed(log(5.06)); label("Apparent peripheral volume of distribution (L)") # S1 Text (b) THETA(10) V4 = 5.06 FIX

    # IIV: log-normal; Table 4 prints variances (omega^2).
    etalka1 ~ 1.21 # Table 4 case 1 mTMDD IIV Ka = 1.21; S1 Text (b) OMEGA 1
    etalka2 ~ 0.36 # S1 Text (b) OMEGA 2 on Ka2 = 0.36 (not in Table 4)
    etalkrec ~ 0.398 # Table 4 case 1 mTMDD IIV Krecycle = 0.398; S1 Text (b) OMEGA 3
    etalkdeg ~ 0.0798 # Table 4 case 1 mTMDD IIV Kdeg = 0.0798; S1 Text (b) OMEGA 4 = 0.0789
    etalcl ~ 0.0389 # Table 4 case 1 mTMDD IIV CL = 0.0389; S1 Text (b) OMEGA 5 = 0.039
    etalkup ~ 0.927 # Table 4 case 1 mTMDD IIV Kuptake = 0.927; S1 Text (b) OMEGA 6
    etaltlag ~ 0.323 # S1 Text (b) OMEGA 7 on Alag1 = 0.323 (not in Table 4)

    # Residual error: W = SQRT(THETA(16)^2 + THETA(17)^2 * IPRED^2), SIGMA 1 FIX.
    addSd <- 0.237; label("Additive residual error (nmol/L)") # Table 4 case 1 mTMDD additive = 0.237; S1 Text (b) THETA(16)
    propSd <- 0.0984; label("Proportional residual error (fraction)") # Table 4 case 1 mTMDD proportional = 0.0984; S1 Text (b) THETA(17) 0.0987
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
    a_kss1 <- exp(la_kss1)
    kss2 <- exp(lkss2)
    kint <- exp(lkint)
    c_il1r_t <- exp(lc_il1r_t)
    a_fcrn_t <- exp(la_fcrn_t)

    kel <- cl / vc
    kcp <- q / vc
    kpc <- q / vp

    # S1 Text (b) $DES, standard QSSA (Michaelis-Menten) on the totals:
    # Acpx = FcRn*A(2)/(KSS1+A(2)); Cf = A(3)/V3; Ccpx = Rtot*Cf/(KSS2+Cf).
    a_fcrn_d <- a_fcrn_t * abs_site / (a_kss1 + abs_site)
    c_tot <- central / vc
    c_cpx <- c_il1r_t * c_tot / (kss2 + c_tot)

    # States in the S1 Text (b) $MODEL order: DEPOT, ABS, CENTRAL, PERIPH.
    d/dt(depot) <- -ka1 * depot
    d/dt(abs_site) <- ka1 * depot - (kdeg + ka2) * abs_site - krec * a_fcrn_d + kup * central
    d/dt(central) <- krec * a_fcrn_d + ka2 * abs_site + kpc * peripheral1 - (kel + kcp) * central -
      kint * c_cpx * vc - kup * central
    d/dt(peripheral1) <- kcp * central - kpc * peripheral1

    alag(depot) <- tlag

    # S1 Text (b) $ERROR: IPRED = A(3)/V3, the TOTAL central concentration.
    Cc <- c_tot
    Cc ~ add(addSd) + prop(propSd)
  })
}
