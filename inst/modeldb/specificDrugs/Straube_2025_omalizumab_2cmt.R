Straube_2025_omalizumab_2cmt <- function() {
  description <- "Two-compartment Mager-Jusko TMDD model with explicit drug-target binding for subcutaneous omalizumab and its target IgE (Straube 2025 supplement Table S1, Figure S3a patient); the improved-fit refit of the one-compartment model in Table 2"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in supplement Table S1 were estimated by Straube by fitting the",
    "two-compartment TMDD model in Equations (S18)-(S21) to omalizumab",
    "total-drug, total-target and free-target time courses digitised from",
    "Meno-Tetang GML, Lowe PJ. Basic Clin Pharmacol Toxicol. 2005;96:182-192.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source; specimen is not, so the
  # circulating states are verified = FALSE. peripheral1 holds FREE drug only
  # -- Equation (S18) gives the peripheral compartment no target and no
  # binding, so no complex accumulates there.
  compartmentData <- list(
    depot       = list(analyte = "omalizumab",             units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "omalizumab",             units = "nmol", specimen = "plasma",              verified = FALSE),
    peripheral1 = list(analyte = "omalizumab",             units = "nmol", specimen = "tissue",              verified = FALSE),
    target      = list(analyte = "IgE",                    units = "nmol", specimen = "plasma",              verified = FALSE),
    complex     = list(analyte = "omalizumab-IgE complex", units = "nmol", specimen = "plasma",              verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 2L,
    n_studies     = 1L,
    disease_state = "Atopic disease; omalizumab is a monoclonal antibody against IgE used to treat atopic diseases (Straube 2025 section 3.2).",
    dose_range    = "Single subcutaneous dose. 90 mg for the Figure S3a patient encoded here; the second patient (Table S1 'Fig. S2b' column) received 270 mg.",
    regions       = NA_character_,
    notes         = paste(
      "Same two digitised phase I patients as Straube_2025_omalizumab_1cmt, refitted with the",
      "two-compartment TMDD model of Equations (S18)-(S21). Straube states that the fits to the observed",
      "data 'can be improved with a two-compartment model' (section 3.2). Table S1 footnote: Vc, Vp, CL and",
      "Rb were estimated patient-specific and the remaining parameters are pooled, so this file encodes an",
      "INDIVIDUAL fit with no IIV and no residual error. The second patient's values are Vc = 1.34 L,",
      "Vp = 4.41 L, CL = 0.13 L/day and Rb = 1.1 nM (Table S1, 'Fig. S2b' column).",
      "Unlike the one-compartment fit in Table 2, Kd was ESTIMATED here (1.96 nM): the Table S1 footnote",
      "marks only F as fixed for omalizumab."
    )
  )

  ini({
    # Drug disposition. Straube 2025 supplement Table S1, Omalizumab / "Fig. S3a" column.
    lvc     <- log(0.68);         label("Central volume of distribution (Vc, L)")                    # Table S1, Fig. S3a column: Vc = 0.68 L (estimated patient-specific)
    lvp     <- log(2.81);         label("Peripheral volume of distribution (Vp, L)")                 # Table S1, Fig. S3a column: Vp = 2.81 L (estimated patient-specific)
    lcl     <- log(0.09);         label("Systemic clearance (CL, L/day)")                            # Table S1, Fig. S3a column: CL = 0.09 L/day (estimated patient-specific); Table S1 reports the derived keD = CL/Vc = 0.13 1/day
    lq      <- log(4.3);          label("Intercompartmental clearance (Q, L/day)")                   # Table S1, Omalizumab: Q = 4.3 L/day (pooled across the two patients)
    lka     <- log(0.31);         label("First-order absorption rate constant (ka, 1/day)")          # Table S1, Omalizumab: ka = 0.31 1/day (pooled)
    lfdepot <- fixed(log(0.42));  label("Subcutaneous bioavailability (F, fraction)")                # Table S1 footnote *: "F fixed at value reported in Stein & Ramakrishna" -> F = 0.42

    # Drug-target binding. For the two-compartment refit BOTH koff and Kd were
    # estimated (only F carries the fixed-value footnote), so kon = koff/Kd is
    # still derived in model() (Equation 7).
    lk2     <- log(1.98);         label("Dissociation (off) rate constant of drug-target binding (koff, 1/day)")  # Table S1, Omalizumab: koff = 1.98 1/day (pooled)
    lkd     <- log(1.96);         label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table S1, Omalizumab: Kd = 1.96 nM (estimated; NOT fixed here, unlike Table 2). Table S1 reports the derived kon = koff/Kd = 1.01 1/(nM*day)

    # Target turnover (Figure 1): ksyn = keR*Rb, baseline Rb = ksyn/keR.
    lrbase  <- log(1.28);         label("Baseline (basal) free target concentration (Rb, nM)")       # Table S1, Fig. S3a column: Rb = 1.28 nM (estimated patient-specific); Table S1 reports the derived ksyn = keR*Rb = 1.15 nM/day
    lkdeg   <- log(0.90);         label("Free target elimination rate constant (keR, 1/day)")        # Table S1, Omalizumab: keR = 0.90 1/day (pooled)
    lkint   <- log(0.16);         label("Drug-target complex elimination rate constant (keDR, 1/day)")  # Table S1, Omalizumab: keDR = keR/Tacc = 0.90/5.7 = 0.16 1/day. Tacc = 5.7 (Equation 14) is recoverable as kdeg/kint and so is not carried separately
  })

  model({
    # Deterministic individual fit: no IIV and no residual error are reported.
    vc    <- exp(lvc)
    vp    <- exp(lvp)
    cl    <- exp(lcl)
    q     <- exp(lq)
    ka    <- exp(lka)
    k2    <- exp(lk2)
    kd    <- exp(lkd)
    rbase <- exp(lrbase)
    kdeg  <- exp(lkdeg)
    kint  <- exp(lkint)

    kel  <- cl / vc      # keD = CL/Vc
    k1   <- k2 / kd      # kon = koff/Kd (Equation 7)
    ksyn <- kdeg * rbase # ksyn = keR*Rb (Figure 1)

    # Two-compartment TMDD system, Straube 2025 Equations (S18)-(S21). The
    # paper writes these in concentrations; the states below hold AMOUNTS in
    # nmol, central within vc and peripheral1 within vp, so the distribution
    # terms are q*(central/vc) and q*(peripheral1/vp). Only free drug
    # distributes: Equation (S18) gives the periphery no target and no binding.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - k1 * central * target / vc + k2 * complex - kel * central +
                          q * (peripheral1 / vp - central / vc)
    d/dt(peripheral1) <-  q * (central / vc - peripheral1 / vp)
    d/dt(target)      <- -k1 * central * target / vc + k2 * complex + ksyn * vc - kdeg * target
    d/dt(complex)     <-  k1 * central * target / vc - k2 * complex - kint * complex

    # Equation (S19): the target starts at its drug-free steady state Rb.
    target(0) <- rbase * vc

    # Bioavailability applied to the dose; see the vignette Errata for the
    # numerical check that selects this form over the Equation (S20) depot
    # loss term.
    f(depot) <- exp(lfdepot)

    # Observations. Figure S3a plots total drug (DT), total target (RT) and
    # free target (R); the FTBR is Equation (28). Total drug is the central
    # (plasma) species only -- peripheral1 is a separate tissue compartment.
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
