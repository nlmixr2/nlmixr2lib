Byun_2024_HL2351_qtmdd <- function() {
  description <- paste(
    "Quasi-steady-state (total QSSA, qTMDD) approximation of the target-mediated drug",
    "disposition population PK model for subcutaneous HL2351 (hIL-1Ra-hyFc, ~97 kDa) in healthy",
    "adult Korean men, with FcRn-mediated recycling; case 1 'qTMDD' model of Byun 2024. Same",
    "structure as modellib('Byun_2024_HL2351_tmdd') but the states hold TOTAL drug and the free",
    "drug is the positive root of the QSS quadratic, both for FcRn in the absorption space",
    "(amount basis, KSS1 in nmol) and for IL1R in central (concentration basis, KSS2 in nmol/L).",
    "Linear clearance, uptake, and peripheral distribution act on the free central drug. The",
    "observed quantity is the FREE concentration. Amounts nmol, concentrations nmol/L."
  )
  reference <- paste(
    "Byun JH, Jeon HS, Yun H-y, Kim JK. Validity conditions of approximations for a",
    "target-mediated drug disposition model: A novel first-order approximation and its",
    "comparison to other approximations. PLoS Comput Biol. 2024;20(4):e1012066.",
    "doi:10.1371/journal.pcbi.1012066. Case 1 qTMDD model: Table 4 (estimates) and",
    "S1 Text (c) (NONMEM control stream). The clinical study and base structure are from",
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
    # Final estimates are Byun 2024 Table 4 (case 1, qTMDD column). Parameters that
    # Table 4 does not list (Ka2, Alag1, and the IIVs on Ka2 and Alag1) are taken
    # from the S1 Text (c) $THETA / $OMEGA, which carry the final estimates re-used
    # as starting values.
    lka1 <- log(1.21); label("Absorption rate from injection-site depot to the absorption space (1/h)") # Table 4 case 1 qTMDD Ka = 1.21; S1 Text (c) THETA(1) 1.19
    lka2 <- log(0.017); label("Direct rate of free drug from the absorption space to central (1/h)") # S1 Text (c) THETA(4) Ka2 = 0.017 (not in Table 4)
    lkrec <- log(0.0332); label("Recycling rate of the drug-FcRn complex into central (1/h)") # Table 4 case 1 qTMDD Krecycle = 0.0332; S1 Text (c) THETA(5) 0.0331
    lcl <- log(0.207); label("Apparent clearance of free drug from central (L/h)") # Table 4 case 1 qTMDD CL = 0.207; S1 Text (c) THETA(7) 0.208
    lvc <- log(11.2); label("Apparent central volume of distribution (L)") # Table 4 case 1 qTMDD Vc = 11.2; S1 Text (c) THETA(8) 11.1
    lq <- log(0.0292); label("Apparent inter-compartmental clearance (L/h)") # Table 4 case 1 qTMDD Q = 0.0292; S1 Text (c) THETA(9)
    la_kss1 <- log(233); label("Drug-FcRn QSS constant in the absorption space (nmol)") # Table 4 case 1 qTMDD KSS1 = 233; S1 Text (c) THETA(2) 227
    lkss2 <- log(14.1); label("Drug-IL1R QSS constant in central (nmol/L)") # Table 4 case 1 qTMDD KSS2 = 14.1; S1 Text (c) THETA(13) 14
    ltlag <- log(0.304); label("Subcutaneous absorption lag time (h)") # S1 Text (c) THETA(15) Alag1 = 0.304 (not in Table 4)

    # Held constant at literature values (Table 4 asterisk; S1 Text (c) FIX).
    lkdeg <- fixed(log(0.0264)); label("Degradation rate of free drug in the absorption space (1/h)") # S1 Text (c) THETA(6) 0.0264 FIX; Table 4 prints 0.264 (decimal-point slip, see vignette)
    lkint <- fixed(log(0.206)); label("Internalisation rate of the IL1R-drug complex (1/h)") # Table 4 case 1 Kint = 0.206 (literature); S1 Text (c) THETA(11) Kdeg2
    lkup <- fixed(log(0.00952)); label("Uptake rate of free drug from central into the absorption space (1/h)") # Table 4 case 1 Kuptake = 0.00952 (literature); S1 Text (c) THETA(14)
    lc_il1r_t <- fixed(log(2.23)); label("Total IL1R concentration in central (nmol/L)") # Table 4 case 1 Rtot = 2.23 (literature); S1 Text (c) THETA(12)
    la_fcrn_t <- fixed(log(749)); label("Total FcRn amount in the absorption space (nmol)") # S1 Text (c) THETA(3) FcRn = 749 FIX
    lvp <- fixed(log(5.06)); label("Apparent peripheral volume of distribution (L)") # S1 Text (c) THETA(10) V4 = 5.06 FIX

    # IIV: log-normal; Table 4 prints variances (omega^2).
    etalka1 ~ 0.94 # Table 4 case 1 qTMDD IIV Ka = 0.94; S1 Text (c) OMEGA 1 = 0.927
    etalka2 ~ 0.524 # S1 Text (c) OMEGA 2 on Ka2 = 0.524 (not in Table 4)
    etalkrec ~ 0.110 # Table 4 case 1 qTMDD IIV Krecycle = 0.110; S1 Text (c) OMEGA 3 = 0.109
    etalkdeg ~ 0.0751 # Table 4 case 1 qTMDD IIV Kdeg = 0.0751; S1 Text (c) OMEGA 4 = 0.0752
    etalcl ~ 0.0532 # Table 4 case 1 qTMDD IIV CL = 0.0532; S1 Text (c) OMEGA 5 = 0.0535
    etalkup ~ 0.813 # Table 4 case 1 qTMDD IIV Kuptake = 0.813; S1 Text (c) OMEGA 6 = 0.801
    etaltlag ~ 0.319 # S1 Text (c) OMEGA 7 on Alag1 = 0.319 (not in Table 4)

    # Residual error: W = SQRT(THETA(16)^2 + THETA(17)^2 * IPRED^2), SIGMA 1 FIX.
    addSd <- 0.184; label("Additive residual error (nmol/L)") # Table 4 case 1 qTMDD additive = 0.184; S1 Text (c) THETA(16)
    propSd <- 0.115; label("Proportional residual error (fraction)") # Table 4 case 1 qTMDD proportional = 0.115; S1 Text (c) THETA(17)
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

    # S1 Text (c) $DES, total QSSA. Absorption space (amounts):
    # DAA = A(2)-FcRn-KSS1; Afree = 0.5*(DAA+SQRT(DAA^2+4*KSS1*A(2)));
    # Acpx = FcRn*Afree/(KSS1+Afree).
    daa <- abs_site - a_fcrn_t - a_kss1
    a_free <- 0.5 * (daa + sqrt(daa * daa + 4 * a_kss1 * abs_site))
    a_fcrn_d <- a_fcrn_t * a_free / (a_kss1 + a_free)
    # Central (concentrations): D = Ct-Rtot-KSS2; CP = 0.5*(D+SQRT(D^2+4*KSS2*Ct));
    # Ccpx = Rtot*CP/(KSS2+CP).
    c_tot <- central / vc
    dcen <- c_tot - c_il1r_t - kss2
    c_free <- 0.5 * (dcen + sqrt(dcen * dcen + 4 * kss2 * c_tot))
    c_cpx <- c_il1r_t * c_free / (kss2 + c_free)

    # States in the S1 Text (c) $MODEL order: DEPOT, ABS, CENTRAL, PERIPH.
    d/dt(depot) <- -ka1 * depot
    d/dt(abs_site) <- ka1 * depot - (kdeg + ka2) * a_free - krec * a_fcrn_d + kup * c_free * vc
    d/dt(central) <- krec * a_fcrn_d + ka2 * a_free + kpc * peripheral1 - (kel + kcp) * c_free * vc -
      kint * c_cpx * vc - kup * c_free * vc
    d/dt(peripheral1) <- kcp * c_free * vc - kpc * peripheral1

    alag(depot) <- tlag

    # S1 Text (c) $ERROR: IPRED = Cfree from the same quadratic, the FREE concentration.
    Cc <- c_free
    Cc ~ add(addSd) + prop(propSd)
  })
}
