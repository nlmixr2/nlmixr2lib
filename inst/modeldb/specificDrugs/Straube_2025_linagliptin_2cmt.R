Straube_2025_linagliptin_2cmt <- function() {
  description <- "Two-compartment Mager-Jusko TMDD model with explicit drug-target binding for intravenous linagliptin and its target dipeptidyl peptidase-4 (Straube 2025 supplement Table S1); the improved-fit refit of the one-compartment model in Table 2"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in supplement Table S1 were estimated by Straube by fitting the",
    "two-compartment TMDD model in Equations (S18)-(S19) to linagliptin",
    "total-drug time courses digitised from Glassman PM, Muzykantov VR.",
    "Target-Mediated Exposure Enhancement: A Previously Unexplored Limit of TMDD. J Pharmacokinet Pharmacodyn. 2020;47(5):411-420; doi:10.1007/s10928-020-09693-1.",
    "koff and Kd were fixed at values reported in Wu N, An G. AAPS J. 2020;22:125;",
    "doi:10.1208/s12248-020-00514-4.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source; specimen is not, hence
  # verified = FALSE. peripheral1 holds FREE drug only -- Equation (S18) gives
  # the peripheral compartment no target and no binding.
  compartmentData <- list(
    central     = list(analyte = "linagliptin",              units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "linagliptin",              units = "nmol", specimen = "tissue", verified = FALSE),
    target      = list(analyte = "dipeptidyl peptidase-4",   units = "nmol", specimen = "plasma", verified = FALSE),
    complex     = list(analyte = "linagliptin-DPP-4 complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    disease_state = "Linagliptin is a small molecule inhibitor of dipeptidyl peptidase-4 (DPP-4) used to treat type 2 diabetes (Straube 2025 section 3.2).",
    dose_range    = "Single intravenous doses of 0.5, 2.5 and 10 mg (Figure S3c).",
    regions       = NA_character_,
    notes         = paste(
      "Same digitised data as Straube_2025_linagliptin_1cmt, refitted with the two-compartment TMDD model",
      "of Equations (S18)-(S19). Species is not stated explicitly; see the one-compartment file's note.",
      "Table S1 footnote +: only koff and Kd were fixed (at Wu and An 2020 values); Vc, Vp, CL, Q, Rb, keR",
      "and Tacc were estimated. The two-compartment fit gives a markedly slower target turnover than the",
      "one-compartment fit (keR 0.550 vs 9.988 1/day, Tacc 5.399 vs 86.804), and Straube notes that with",
      "these estimates 'we see marked deviations of the approximations ... especially for the FTBR'",
      "(Figure S3 caption). Straube reports no IIV and no residual error, so this is a deterministic",
      "typical-value fit."
    )
  )

  ini({
    # Drug disposition. Straube 2025 supplement Table S1, Linagliptin column.
    lvc    <- log(56.25);          label("Central volume of distribution (L)")                    # Table S1, Linagliptin: Vc = 56.25 L
    lvp    <- log(117.865);        label("Peripheral volume of distribution (L)")                 # Table S1, Linagliptin: Vp = 117.865 L
    lcl    <- log(639.614);        label("Systemic clearance (L/day)")                            # Table S1, Linagliptin: CL = 639.614 L/day; Table S1 reports the derived keD = CL/Vc = 11.37 1/day (reproduced exactly)
    lq     <- log(495.636);        label("Intercompartmental clearance (L/day)")                   # Table S1, Linagliptin: Q = 495.636 L/day

    # Drug-target binding: koff and Kd fixed, kon = koff/Kd derived (Equation 7).
    lk2    <- fixed(log(1.675));   label("Dissociation (off) rate constant of drug-target binding (1/day)")  # Table S1 footnote +: koff fixed at 1.675 1/day from Wu and An
    lkd    <- fixed(log(0.074));   label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table S1 footnote +: Kd fixed at 0.074 nM from Wu and An. Table S1 prints kon = 22.73 1/(nM*day), which implies Kd = 0.0737; koff/Kd from the rounded printed values gives 22.64. See vignette Errata.

    # Target turnover (Figure 1): ksyn = keR*Rb, baseline Rb = ksyn/keR.
    lrbase <- log(2.861);          label("Baseline (basal) free target concentration (nM)")       # Table S1, Linagliptin: Rb = 2.861 nM (estimated); Table S1 reports the derived ksyn = keR*Rb = 1.573 nM/day
    lkdeg  <- log(0.550);          label("Free target elimination rate constant (1/day)")        # Table S1, Linagliptin: keR = 0.550 1/day (estimated)
    lkint  <- log(0.102);          label("Drug-target complex elimination rate constant (1/day)")  # Table S1, Linagliptin: keDR = keR/Tacc = 0.550/5.399 = 0.102 1/day. Tacc = 5.399 (Equation 14) was estimated here, not fixed; it is recoverable as kdeg/kint and so is not carried separately.
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

    # Observations. Figure S3c plots total drug (DT) and the FTBR (Equation 28).
    # Total drug is the central (plasma) species only.
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
