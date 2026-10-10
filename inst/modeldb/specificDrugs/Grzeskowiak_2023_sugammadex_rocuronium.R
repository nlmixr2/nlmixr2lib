Grzeskowiak_2023_sugammadex_rocuronium <- function() {
  description <- paste(
    "Bayesian population PK-PD model of rocuronium-induced neuromuscular",
    "blockade and its reversal by sugammadex in ASA I-II children aged",
    "2.5-17 years (Grzeskowiak 2023). The structural model of Kleijn 2011 was",
    "refitted with informative priors: two-compartment PK for free",
    "rocuronium, free sugammadex and the rocuronium-sugammadex complex (complex",
    "PK set equal to sugammadex), dynamic complexation in the central",
    "compartments (kd fixed at 0.0559 uM, dissociation rate k2 estimated), a",
    "rocuronium effect compartment (ke0) driving a sigmoid Emax model of the",
    "train-of-four (TOF) ratio with E0 = Emax = 100 fixed, and a",
    "sugammadex-driven second-order elimination of rocuronium from the effect",
    "compartment (ks). Body weight is the only covariate: allometric exponents",
    "1 on volumes, 0.75 on clearances and inter-compartmental clearances, and",
    "-0.25 on ke0 and ks, standardised to 70 kg. Posterior-mean parameter",
    "values; doses in umol, concentrations in uM, time in min."
  )
  reference <- paste(
    "Grzeskowiak M, Bienert A, Wiczling P, Malec M, Grzelak J, Jarosz K, Ber J,",
    "Ksiazkiewicz M, Rosada-Kurasinska J, Grzeskowiak E, Bartkowska-Sniatkowska",
    "A. Population Pharmacokinetic-Pharmacodynamic Modeling and Probability of",
    "Target Attainment Analysis of Rocuronium and Sugammadex in Children",
    "Undergoing Surgery. Eur J Drug Metab Pharmacokinet. 2023;48:101-114.",
    "doi:10.1007/s13318-022-00809-1. Parameter values from Table S1 and the",
    "NONMEM control stream (S2) of Supplementary Material 2."
  )
  vignette <- "Grzeskowiak_2023_sugammadex_rocuronium"
  units <- list(
    time = "min",
    dosing = "umol",
    concentration = "umol/L"
  )

  compartmentData <- list(
    central = list(analyte = "sugammadex", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sugammadex", units = "umol", specimen = "plasma", verified = TRUE),
    central_roc = list(analyte = "rocuronium", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_roc = list(analyte = "rocuronium", units = "umol", specimen = "plasma", verified = TRUE),
    central_complex = list(
      analyte = "sugammadex-rocuronium complex",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_complex = list(
      analyte = "sugammadex-rocuronium complex",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    effect_roc = list(
      analyte = "rocuronium (biophase concentration)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling standardised to 70 kg (control stream: LNWT - 4.2485",
        "with log(70) = 4.2485). Exponents fixed at 1 for VR1, VR2, VS1, VS2;",
        "0.75 for CLR1, CLR2, CLS1, CLS2 (CLR2 / CLS2 are the inter-compartmental",
        "clearances); -0.25 for ke0 and ks; 0 for EC50, k2 and the Hill",
        "coefficient (Methods Section 2.3 and S2 $PK). Unlike Kleijn 2011, the",
        "sugammadex clearance is allometrically scaled because the children had",
        "normal renal function and body weight was assumed to be the only",
        "size factor."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "In the dataset ($INPUT AGE) and inspected graphically against the",
        "posterior-mean etas (Supplementary Material 2, S5); no age effect was",
        "included in the final model."
      )
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Dataset column SEX (1 = male, 2 = female; 22 / 8 subjects, matching",
        "Results Section 3.1). Inspected graphically against the posterior-mean",
        "etas (S5); no sex effect was included in the final model."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "2.5-17 years",
    age_median = "8.5 years",
    weight_range = "14-85 kg",
    weight_median = "36.5 kg",
    sex_female_pct = round(8 / 30 * 100, 1),
    disease_state = paste(
      "ASA physical status I (n = 19) or II (n = 11) children undergoing",
      "surgery under propofol total intravenous anaesthesia that required",
      "muscle relaxation for more than 20 min"
    ),
    dose_range = paste(
      "Rocuronium 0.6 mg/kg IV bolus at induction with maintenance boluses of",
      "half the initial dose when the TOF ratio exceeded 40%; sugammadex 0.5,",
      "1.0 or 2.0 mg/kg IV bolus at the end of surgery"
    ),
    regions = "Poland (single centre, Poznan)",
    notes = paste(
      "Prospective observational cohort, NCT04851574 (Methods 2.1, Results 3.1).",
      "223 rocuronium and 130 sugammadex plasma concentrations and 992 TOF-ratio",
      "measurements (Results 3.2). Fitted by full Bayesian MCMC (NONMEM",
      "METHOD=BAYES with $PRIOR NWPRI) using the Kleijn 2011 adult/paediatric",
      "estimates as informative priors."
    )
  )

  ini({
    # Typical values are the posterior means (Table S1 part A, column
    # 'Posterior thetaP') at a 70 kg reference subject.

    # Sugammadex PK (S2: VS1, CLS1, VS2, CLS2)
    lcl <- log(0.093); label("Sugammadex clearance CLS1 at 70 kg (L/min)") # Table S1A CLS1 posterior thetaP 0.093 (log -2.37), RSE 1.7%
    lvc <- log(4.31); label("Sugammadex central volume VS1 at 70 kg (L)") # Table S1A VS1 posterior thetaP 4.31 (log 1.46), RSE 3%
    lq <- log(0.206); label("Sugammadex inter-compartmental clearance CLS2 at 70 kg (L/min)") # Table S1A CLS2 posterior thetaP 0.206 (log -1.58), RSE 2.9%
    lvp <- log(6.49); label("Sugammadex peripheral volume VS2 at 70 kg (L)") # Table S1A VS2 posterior thetaP 6.49 (log 1.87), RSE 2.4%

    # Rocuronium PK (S2: VR1, CLR1, VR2, CLR2)
    lcl_roc <- log(0.278); label("Rocuronium clearance CLR1 at 70 kg (L/min)") # Table S1A CLR1 posterior thetaP 0.278 (log -1.28), RSE 3.5%
    lvc_roc <- log(5.00); label("Rocuronium central volume VR1 at 70 kg (L)") # Table S1A VR1 posterior thetaP 5.00 (log 1.61), RSE 2.8%
    lq_roc <- log(0.284); label("Rocuronium inter-compartmental clearance CLR2 at 70 kg (L/min)") # Table S1A CLR2 posterior thetaP 0.284 (log -1.26), RSE 3.9%
    lvp_roc <- log(7.17); label("Rocuronium peripheral volume VR2 at 70 kg (L)") # Table S1A VR2 posterior thetaP 7.17 (log 1.97), RSE 2.5%

    # Complexation (S2: KD = 0.0559; KON = KOFF/KD)
    kd <- fixed(0.0559); label("Equilibrium dissociation constant of the rocuronium-sugammadex complex (uM)") # Methods 2.3 'Kd included as a fixed parameter (0.0559 mM)'; S2 $PK KD = 0.0559 on the uM scale of $DES
    lk2 <- log(0.022); label("Complex dissociation rate constant kOFF (1/min)") # Table S1A kOFF posterior thetaP 0.022 (log -3.83), RSE 6.6%

    # Rocuronium PD (S2: EC50R, GAM, KEO; E0 = Emax = 100 in $ERROR)
    lke0 <- log(0.185); label("Rocuronium effect-compartment equilibration rate kEO at 70 kg (1/min)") # Table S1A kEO posterior thetaP 0.185 (log -1.69), RSE 6.2%
    lec50 <- log(1.43); label("Biophase rocuronium concentration giving 50% TOF-ratio reduction EC50R (uM)") # Table S1A EC50R posterior thetaP 1.43 (log 0.36), RSE 23.3%
    lhill <- log(6.96); label("Hill coefficient GAM of the TOF-ratio Emax model (unitless)") # Table S1A GAM posterior thetaP 6.96 (log 1.94), RSE 7.9%
    e0 <- fixed(100); label("Baseline TOF ratio E0 (%)") # Methods 2.3 'baseline TOF ratio (E0: fixed to 100)'; S2 $ERROR IPRE = 100*(...)
    emax <- fixed(100); label("Maximal reduction in TOF ratio Emax (%)") # Methods 2.3 'maximal reduction in TOF ratio (Emax: fixed to 100)'

    # Sugammadex-mediated elimination of rocuronium from the biophase (S2: KS)
    lks <- log(0.080); label("Second-order rate constant kS of sugammadex-driven biophase rocuronium elimination at 70 kg (1/(uM*min))") # Table S1A kS posterior thetaP 0.080 (log -2.53), RSE 6.7%

    # IIV: Table S1 part B posterior Omega in % is 100*sqrt(omega^2) (the
    # prior Omega0 = 30% corresponds to $OMEGAP diagonal 0.09), so the
    # variance is (Omega/100)^2. The posterior off-diagonals of the
    # estimated $OMEGA BLOCK(8) / BLOCK(5) are not reported; diagonal only.
    etalcl ~ 0.153664 # Table S1B CLS1 posterior Omega 39.2% -> 0.392^2
    etalvc ~ 0.172225 # Table S1B VS1 posterior Omega 41.5% -> 0.415^2
    etalq ~ 0.259081 # Table S1B CLS2 posterior Omega 50.9% -> 0.509^2
    etalvp ~ 0.203401 # Table S1B VS2 posterior Omega 45.1% -> 0.451^2
    etalcl_roc ~ 0.148225 # Table S1B CLR1 posterior Omega 38.5% -> 0.385^2
    etalvc_roc ~ 0.1764 # Table S1B VR1 posterior Omega 42.0% -> 0.420^2
    etalq_roc ~ 0.2704 # Table S1B CLR2 posterior Omega 52.0% -> 0.520^2
    etalvp_roc ~ 0.2809 # Table S1B VR2 posterior Omega 53.0% -> 0.530^2
    etalk2 ~ 0.177241 # Table S1B kOFF posterior Omega 42.1% -> 0.421^2
    etalke0 ~ 0.207936 # Table S1B kEO posterior Omega 45.6% -> 0.456^2
    etalec50 ~ 0.142884 # Table S1B EC50R posterior Omega 37.8% -> 0.378^2
    etalhill ~ 0.772641 # Table S1B GAM posterior Omega 87.9% -> 0.879^2
    etalks ~ 0.683929 # Table S1B kS posterior Omega 82.7% -> 0.827^2

    # Residual error (Table S1 part C posterior sigma; S2 $ERROR forms)
    propSd <- 0.269; label("Proportional residual SD, total sugammadex plasma concentration (fraction)") # Table S1C SUG prop posterior sigma 26.9%, RSE 9.5%
    propSd_roc <- 0.252; label("Proportional residual SD, total rocuronium plasma concentration (fraction)") # Table S1C ROC prop posterior sigma 25.2%, RSE 6.8%
    addSd_NMB <- 11.6; label("Additive residual SD, TOF ratio (%)") # Table S1C TOF add posterior sigma 11.6, RSE 2.5%
  })

  model({
    # Allometric scaling to 70 kg (S2 $PK: exponents 1 / 0.75 / -0.25 / 0)
    wt_ratio <- WT / 70

    vc <- exp(lvc + etalvc) * wt_ratio
    cl <- exp(lcl + etalcl) * wt_ratio^0.75
    vp <- exp(lvp + etalvp) * wt_ratio
    q <- exp(lq + etalq) * wt_ratio^0.75

    vc_roc <- exp(lvc_roc + etalvc_roc) * wt_ratio
    cl_roc <- exp(lcl_roc + etalcl_roc) * wt_ratio^0.75
    vp_roc <- exp(lvp_roc + etalvp_roc) * wt_ratio
    q_roc <- exp(lq_roc + etalq_roc) * wt_ratio^0.75

    ec50 <- exp(lec50 + etalec50)
    k2 <- exp(lk2 + etalk2)
    hill <- exp(lhill + etalhill)
    ke0 <- exp(lke0 + etalke0) * wt_ratio^(-0.25)
    ks <- exp(lks + etalks) * wt_ratio^(-0.25)

    # Association rate from the fixed equilibrium constant (S2: KON = KOFF/KD)
    k1 <- k2 / kd

    # Concentrations (uM); the complex shares the sugammadex volumes
    cs1 <- central / vc
    cs2 <- peripheral1 / vp
    cr1 <- central_roc / vc_roc
    cr2 <- peripheral1_roc / vp_roc
    crs1 <- central_complex / vc
    crs2 <- peripheral1_complex / vp

    # S2 $DES, transcribed term for term. The complexation flux is scaled by
    # the rocuronium central volume in the rocuronium equation and by the
    # sugammadex central volume in the sugammadex and complex equations.
    d/dt(central) <- -q * (cs1 - cs2) - cl * cs1 - k1 * cr1 * cs1 * vc + k2 * crs1 * vc
    d/dt(peripheral1) <- q * (cs1 - cs2)
    d/dt(central_roc) <- -q_roc * (cr1 - cr2) - cl_roc * cr1 - k1 * cr1 * cs1 * vc_roc + k2 * crs1 * vc_roc
    d/dt(peripheral1_roc) <- q_roc * (cr1 - cr2)
    d/dt(central_complex) <- k1 * cr1 * cs1 * vc - k2 * crs1 * vc - cl * crs1 - q * (crs1 - crs2)
    d/dt(peripheral1_complex) <- q * (crs1 - crs2)
    # Biophase rocuronium concentration (uM); sugammadex removes it at ks*cs1
    d/dt(effect_roc) <- ke0 * (cr1 - effect_roc) - ks * cs1 * effect_roc

    # Assays measure total drug: free + complex (S2 $ERROR IPRE)
    Cc <- cs1 + crs1
    Cc_roc <- cr1 + crs1

    # TOF ratio (%), S2 $ERROR: 100*(1 - (Ce/EC50)^GAM / (1 + (Ce/EC50)^GAM))
    NMB <- e0 - emax * effect_roc^hill / (ec50^hill + effect_roc^hill)

    Cc ~ prop(propSd)
    Cc_roc ~ prop(propSd_roc)
    NMB ~ add(addSd_NMB)
  })
}
