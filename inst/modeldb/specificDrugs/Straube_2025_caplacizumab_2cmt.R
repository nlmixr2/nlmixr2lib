Straube_2025_caplacizumab_2cmt <- function() {
  description <- "Species not stated in the source. Two-compartment Mager-Jusko TMDD model with explicit drug-target binding for intravenous caplacizumab (ALX-0081) and its target von Willebrand factor (Straube 2025 supplement Table S1); the improved-fit refit of the one-compartment model in Table 2"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in supplement Table S1 were estimated by Straube by fitting the",
    "two-compartment TMDD model in Equations (S18)-(S19) to ALX-0081 total-drug",
    "time courses digitised from Glassman PM, Muzykantov VR.",
    "Target-Mediated Exposure Enhancement: A Previously Unexplored Limit of TMDD. J Pharmacokinet Pharmacodyn. 2020;47(5):411-420; doi:10.1007/s10928-020-09693-1,",
    "which is also the source of the fixed koff, Kd, Rb and Tacc values.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source; specimen is not, hence
  # verified = FALSE. peripheral1 holds FREE drug only -- Equation (S18) gives
  # the peripheral compartment no target and no binding.
  compartmentData <- list(
    central     = list(analyte = "caplacizumab",                              units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "caplacizumab",                              units = "nmol", specimen = "tissue", verified = FALSE),
    target      = list(analyte = "von Willebrand factor",                     units = "nmol", specimen = "plasma", verified = FALSE),
    complex     = list(analyte = "caplacizumab-von Willebrand factor complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "not stated in the source",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    disease_state = "ALX-0081 is a single-domain antibody against von Willebrand factor (vWF) developed to treat acquired thrombotic thrombocytopenic purpura (Straube 2025 section 3.2).",
    dose_range    = "Single intravenous doses of 0.02, 0.4 and 8 mg/kg (Figure S3b).",
    regions       = NA_character_,
    notes         = paste(
      "SPECIES AND BODY WEIGHT ARE NOT STATED; see Straube_2025_caplacizumab_1cmt for the full note.",
      "Same digitised data as the one-compartment fit, refitted with the two-compartment TMDD model of",
      "Equations (S18)-(S19). Table S1 footnote #: koff, Kd, Rb and Tacc were FIXED at values reported in",
      "Glassman and Muzykantov; Vc, Vp, CL, Q and keR were estimated. Straube reports no IIV and no",
      "residual error, so this is a deterministic typical-value fit."
    )
  )

  ini({
    # Drug disposition. Straube 2025 supplement Table S1, ALX-0081 column.
    lvc    <- log(0.042);          label("Central volume of distribution (L)")                    # Table S1, ALX-0081: Vc = 0.042 L. Table S1 also prints the derived keD = CL/Vc = 65.6 1/day, which implies Vc = 0.0415 L; the Vc column is rounded to 3 decimals, so cl/vc here gives 64.8 1/day (1.2% low). See vignette Errata.
    lvp    <- log(0.246);          label("Peripheral volume of distribution (L)")                 # Table S1, ALX-0081: Vp = 0.246 L
    lcl    <- log(2.723);          label("Systemic clearance (L/day)")                            # Table S1, ALX-0081: CL = 2.723 L/day
    lq     <- log(0.871);          label("Intercompartmental clearance (L/day)")                   # Table S1, ALX-0081: Q = 0.871 L/day

    # Drug-target binding: koff and Kd fixed, kon = koff/Kd derived (Equation 7).
    lk2    <- fixed(log(83.76));   label("Dissociation (off) rate constant of drug-target binding (1/day)")  # Table S1 footnote #: koff fixed at 83.76 1/day from Glassman and Muzykantov
    lkd    <- fixed(log(0.0036));  label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table S1 footnote #: Kd fixed at 0.0036 nM from Glassman and Muzykantov. Table S1 reports the derived kon = koff/Kd = 2.3e4 1/(nM*day)

    # Target turnover (Figure 1): ksyn = keR*Rb, baseline Rb = ksyn/keR.
    lrbase <- fixed(log(32.8));    label("Baseline (basal) free target concentration (nM)")       # Table S1 footnote #: Rb fixed at 32.8 nM from Glassman and Muzykantov. Table S1 reports the derived ksyn = keR*Rb = 10.24 nM/day
    lkdeg  <- log(0.312);          label("Free target elimination rate constant (1/day)")        # Table S1, ALX-0081: keR = 0.312 1/day (estimated)
    lkint  <- log(0.935);          label("Drug-target complex elimination rate constant (1/day)")  # Table S1, ALX-0081: keDR = keR/Tacc = 0.312/0.334 = 0.935 1/day. NOT itself fixed: Tacc = 0.334 was fixed (footnote #) but keR was estimated. Tacc (Equation 14) is recoverable as kdeg/kint and so is not carried separately.
  })

  model({
    # Deterministic typical-value fit: no IIV and no residual error reported.
    vc    <- exp(lvc)
    vp    <- exp(lvp)
    cl    <- exp(lcl)
    q     <- exp(lq)
    k2    <- exp(lk2)
    kd    <- exp(lkd)
    rbase <- exp(lrbase)
    kdeg  <- exp(lkdeg)
    kint  <- exp(lkint)

    kel  <- cl / vc      # keD = CL/Vc
    k1   <- k2 / kd      # kon = koff/Kd (Equation 7)
    ksyn <- kdeg * rbase # ksyn = keR*Rb (Figure 1)

    # Two-compartment TMDD system, Straube 2025 Equations (S18)-(S19),
    # intravenous dosing into central. The paper writes these in
    # concentrations; the states below hold AMOUNTS in nmol, central within vc
    # and peripheral1 within vp, so the distribution terms are q*(central/vc)
    # and q*(peripheral1/vp). Only free drug distributes: Equation (S18) gives
    # the periphery no target and no binding.
    d/dt(central)     <- -k1 * central * target / vc + k2 * complex - kel * central +
                          q * (peripheral1 / vp - central / vc)
    d/dt(peripheral1) <-  q * (central / vc - peripheral1 / vp)
    d/dt(target)      <- -k1 * central * target / vc + k2 * complex + ksyn * vc - kdeg * target
    d/dt(complex)     <-  k1 * central * target / vc - k2 * complex - kint * complex

    # Equation (S19): the target starts at its drug-free steady state Rb.
    target(0) <- rbase * vc

    # Observations. Figure S3b plots total drug (DT) and the FTBR (Equation 28).
    # Total drug is the central (plasma) species only.
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
