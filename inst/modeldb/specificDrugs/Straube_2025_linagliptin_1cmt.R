Straube_2025_linagliptin_1cmt <- function() {
  description <- "One-compartment Mager-Jusko TMDD model with explicit drug-target binding for intravenous linagliptin and its target dipeptidyl peptidase-4 (Straube 2025 Table 2); a high-affinity small-molecule example showing target-mediated exposure enhancement (TMEE)"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in Table 2 were estimated by Straube by fitting the",
    "one-compartment TMDD model in Equation (2) to linagliptin total-drug time",
    "courses digitised from Glassman PM, Muzykantov VR.",
    "Target-Mediated Exposure Enhancement: A Previously Unexplored Limit of TMDD. J Pharmacokinet Pharmacodyn. 2020;47(5):411-420; doi:10.1007/s10928-020-09693-1.",
    "koff and Kd were fixed at values reported in Wu N, An G. AAPS J. 2020;22:125;",
    "doi:10.1208/s12248-020-00514-4.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source (linagliptin is a small
  # molecule inhibitor of dipeptidyl peptidase-4; Straube 2025 section 3.2).
  # specimen is not stated, hence verified = FALSE.
  compartmentData <- list(
    central = list(analyte = "linagliptin",              units = "nmol", specimen = "plasma", verified = FALSE),
    target  = list(analyte = "dipeptidyl peptidase-4",   units = "nmol", specimen = "plasma", verified = FALSE),
    complex = list(analyte = "linagliptin-DPP-4 complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    disease_state = "Linagliptin is a small molecule inhibitor of dipeptidyl peptidase-4 (DPP-4) used to treat type 2 diabetes (Straube 2025 section 3.2).",
    dose_range    = "Single intravenous doses of 0.5, 2.5 and 10 mg (Figure 5c).",
    regions       = NA_character_,
    notes         = paste(
      "Straube 2025 does not state the species explicitly, but the absolute (not per-kg) mg dose levels and",
      "the fitted Vc of 104.854 L are consistent with adult human intravenous data; nlmixr2lib already",
      "carries human linagliptin popPK as Retlich_2015_linagliptin. Data were digitised by Straube from",
      "Glassman and Muzykantov (2020), which is not open access and could not be consulted.",
      "Table 2 footnote b: only koff and Kd were fixed (at Wu and An 2020 values); Vc, CL, Rb, keR and Tacc",
      "were estimated. Straube reports no IIV and no residual error, so this is a deterministic",
      "typical-value fit. A two-compartment refit of the same data is Straube_2025_linagliptin_2cmt, and it",
      "gives a markedly different target turnover (keR 0.550 vs 9.988 1/day, Tacc 5.4 vs 86.8)."
    )
  )

  ini({
    # Drug disposition. Straube 2025 Table 2, Linagliptin column.
    lvc    <- log(104.854);        label("Central volume of distribution (L)")                    # Table 2, Linagliptin: Vc = 104.854 L
    lcl    <- log(937.378);        label("Systemic clearance (L/day)")                            # Table 2, Linagliptin: CL = 937.378 L/day; Table 2 reports the derived keD = CL/Vc = 8.94 1/day (reproduced exactly)

    # Drug-target binding. Both koff and Kd were fixed from Wu and An, so
    # kon = koff/Kd is derived in model() (Equation 7).
    lk2    <- fixed(log(1.675));   label("Dissociation (off) rate constant of drug-target binding (1/day)")  # Table 2 footnote b: koff fixed at 1.675 1/day from Wu and An
    lkd    <- fixed(log(0.074));   label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table 2 footnote b: Kd fixed at 0.074 nM from Wu and An. Table 2 prints kon = 22.73 1/(nM*day), which implies Kd = 0.0737; koff/Kd from the rounded printed values gives 22.64 (0.4% low). See vignette Errata.

    # Target turnover (Figure 1): ksyn = keR*Rb, baseline Rb = ksyn/keR.
    lrbase <- log(1.855);          label("Baseline (basal) free target concentration (nM)")       # Table 2, Linagliptin: Rb = 1.855 nM (estimated); Table 2 reports the derived ksyn = keR*Rb = 18.529 nM/day
    lkdeg  <- log(9.988);          label("Free target elimination rate constant (1/day)")        # Table 2, Linagliptin: keR = 9.988 1/day (estimated)
    lkint  <- log(0.115);          label("Drug-target complex elimination rate constant (1/day)")  # Table 2, Linagliptin: keDR = keR/Tacc = 9.988/86.804 = 0.115 1/day. Tacc = 86.804 (Equation 14) was estimated here, not fixed; it is recoverable as kdeg/kint and so is not carried separately.
  })

  model({
    # Deterministic typical-value fit: no IIV and no residual error reported.
    vc    <- exp(lvc)
    cl    <- exp(lcl)
    k2    <- exp(lk2)
    kd    <- exp(lkd)
    rbase <- exp(lrbase)
    kdeg  <- exp(lkdeg)
    kint  <- exp(lkint)

    kel  <- cl / vc      # keD = CL/Vc
    k1   <- k2 / kd      # kon = koff/Kd (Equation 7)
    ksyn <- kdeg * rbase # ksyn = keR*Rb (Figure 1)

    # Mager-Jusko TMDD system, Straube 2025 Equation (2), intravenous dosing
    # into central (Equation 3). The paper writes Equation (2) in
    # concentrations (D, R, DR in nM); the states below hold the corresponding
    # AMOUNTS in nmol within the central volume vc, so every concentration is
    # state/vc and each bimolecular term carries one 1/vc.
    d/dt(central) <- -k1 * central * target / vc + k2 * complex - kel * central
    d/dt(target)  <- -k1 * central * target / vc + k2 * complex + ksyn * vc - kdeg * target
    d/dt(complex) <-  k1 * central * target / vc - k2 * complex - kint * complex

    # Equation (3): the target starts at its drug-free steady state Rb.
    target(0) <- rbase * vc

    # Observations. Figure 5c plots total drug (DT) and the FTBR (Equation 28).
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
