Straube_2025_caplacizumab_1cmt <- function() {
  description <- "Species not stated in the source. One-compartment Mager-Jusko TMDD model with explicit drug-target binding for intravenous caplacizumab (ALX-0081) and its target von Willebrand factor (Straube 2025 Table 2); a high-affinity example showing target-mediated exposure enhancement (TMEE)"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in Table 2 were estimated by Straube by fitting the",
    "one-compartment TMDD model in Equation (2) to ALX-0081 total-drug time",
    "courses digitised from Glassman PM, Muzykantov VR.",
    "J Pharmacokinet Pharmacodyn. 2020;47:573-591; doi:10.1007/s10928-020-09700-5,",
    "which is also the source of the fixed koff, Kd, Rb and Tacc values.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source (ALX-0081 is a single-domain
  # antibody against von Willebrand factor; Straube 2025 section 3.2).
  # specimen is not stated, hence verified = FALSE.
  compartmentData <- list(
    central = list(analyte = "caplacizumab",                            units = "nmol", specimen = "plasma", verified = FALSE),
    target  = list(analyte = "von Willebrand factor",                   units = "nmol", specimen = "plasma", verified = FALSE),
    complex = list(analyte = "caplacizumab-von Willebrand factor complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "not stated in the source",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    disease_state = "ALX-0081 is a single-domain antibody against von Willebrand factor (vWF) developed to treat acquired thrombotic thrombocytopenic purpura (Straube 2025 section 3.2).",
    dose_range    = "Single intravenous doses of 0.02, 0.4 and 8 mg/kg (Figure 5b).",
    regions       = NA_character_,
    notes         = paste(
      "SPECIES AND BODY WEIGHT ARE NOT STATED. Straube 2025 gives only the mg/kg dose levels and reports",
      "Vc = 0.046 L, which is not consistent with an adult human; the upstream data source",
      "(Glassman and Muzykantov 2020) is not open access and could not be consulted, so the species and the",
      "body weight behind the mg/kg doses could not be resolved. This matters only for converting the mg/kg",
      "doses into the nmol dose the model takes; it does not affect any parameter. See the vignette for how",
      "the nmol doses were derived from the paper's own Vc and Figure 5b.",
      "Table 2 footnote a: koff, Kd, Rb and Tacc were FIXED at values reported in Glassman and Muzykantov;",
      "Vc, CL and keR were estimated. Straube reports no IIV and no residual error, so this is a",
      "deterministic typical-value fit. A two-compartment refit of the same data is",
      "Straube_2025_caplacizumab_2cmt."
    )
  )

  ini({
    # Drug disposition. Straube 2025 Table 2, ALX-0081 column.
    lvc    <- log(0.046);          label("Central volume of distribution (L)")                    # Table 2, ALX-0081: Vc = 0.046 L. Table 2 also prints the derived keD = CL/Vc = 58.80 1/day, which implies Vc = 0.0464 L; the Vc column is rounded to 3 decimals, so cl/vc here gives 59.3 1/day (0.8% high). See vignette Errata.
    lcl    <- log(2.728);          label("Systemic clearance (L/day)")                            # Table 2, ALX-0081: CL = 2.728 L/day

    # Drug-target binding. Both koff and Kd were fixed from Glassman and
    # Muzykantov, so kon = koff/Kd is derived in model() (Equation 7).
    lk2    <- fixed(log(83.76));   label("Dissociation (off) rate constant of drug-target binding (1/day)")  # Table 2 footnote a: koff fixed at 83.76 1/day from Glassman and Muzykantov
    lkd    <- fixed(log(0.0036));  label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table 2 footnote a: Kd fixed at 0.0036 nM from Glassman and Muzykantov. Table 2 reports the derived kon = koff/Kd = 2.3e4 1/(nM*day)

    # Target turnover (Figure 1): ksyn = keR*Rb, baseline Rb = ksyn/keR.
    lrbase <- fixed(log(32.8));    label("Baseline (basal) free target concentration (nM)")       # Table 2 footnote a: Rb fixed at 32.8 nM from Glassman and Muzykantov. Table 2 reports the derived ksyn = keR*Rb = 10.394 nM/day
    lkdeg  <- log(0.317);          label("Free target elimination rate constant (1/day)")        # Table 2, ALX-0081: keR = 0.317 1/day (estimated)
    lkint  <- log(0.948);          label("Drug-target complex elimination rate constant (1/day)")  # Table 2, ALX-0081: keDR = keR/Tacc = 0.317/0.334 = 0.948 1/day. NOT itself fixed: Tacc = 0.334 was fixed (footnote a) but keR was estimated, so keDR inherits keR's uncertainty. Tacc (Equation 14) is recoverable as kdeg/kint and so is not carried separately.
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

    # Observations. Figure 5b plots total drug (DT) and the FTBR (Equation 28).
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
