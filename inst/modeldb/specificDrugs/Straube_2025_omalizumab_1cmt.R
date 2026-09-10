Straube_2025_omalizumab_1cmt <- function() {
  description <- "One-compartment Mager-Jusko TMDD model with explicit drug-target binding for subcutaneous omalizumab and its target IgE (Straube 2025 Table 2, Figure 5a patient); the low-affinity (Michaelis-Menten) example of the paper's high-vs-low affinity classification"
  reference <- paste(
    "Straube R. Target-Mediated Drug Disposition (TMDD) Revisited: High Versus",
    "Low-Affinity Approximations of the TMDD Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1262-1272.",
    "doi:10.1002/psp4.70048.",
    "Parameters in Table 2 were estimated by Straube by fitting the",
    "one-compartment TMDD model in Equation (2) to omalizumab total-drug,",
    "total-target and free-target time courses digitised from",
    "Meno-Tetang GML, Lowe PJ. Basic Clin Pharmacol Toxicol. 2005;96:182-192.",
    sep = " "
  )
  vignette <- "Straube_2025_TMDD_high_vs_low_affinity"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482. analyte is stated by the source (omalizumab is an anti-IgE
  # monoclonal antibody; Straube 2025 section 3.2). specimen is NOT stated --
  # Straube reports only "measured for two patients" -- so verified = FALSE
  # for the circulating states. The depot is verified because the source
  # states the route explicitly ("single dose of subcutaneous injection").
  compartmentData <- list(
    depot   = list(analyte = "omalizumab",             units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "omalizumab",             units = "nmol", specimen = "plasma",              verified = FALSE),
    target  = list(analyte = "IgE",                    units = "nmol", specimen = "plasma",              verified = FALSE),
    complex = list(analyte = "omalizumab-IgE complex", units = "nmol", specimen = "plasma",              verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 2L,
    n_studies     = 1L,
    disease_state = "Atopic disease; omalizumab is a monoclonal antibody against IgE used to treat atopic diseases (Straube 2025 section 3.2).",
    dose_range    = "Single subcutaneous dose. 90 mg for the Figure 5a patient encoded here; the second patient (Table 2 'Figure S2a' column) received 270 mg.",
    regions       = NA_character_,
    notes         = paste(
      "Phase I data for two patients, digitised by Straube from Meno-Tetang and Lowe (2005) and",
      "re-fitted with the one-compartment TMDD model of Equation (2). Table 2 footnote: Vc, CL and Rb",
      "were estimated patient-specific and the remaining parameters are pooled estimates across the two",
      "patients, so this file encodes an INDIVIDUAL fit, not a population model -- there is no IIV and no",
      "residual error to encode. The second patient's values are Vc = 6.566 L, CL = 0.051 L/day and",
      "Rb = 1.139 nM (Table 2, 'Figure S2a' column); all other parameters are shared with this file.",
      "A two-compartment refit of the same data is Straube_2025_omalizumab_2cmt."
    )
  )

  ini({
    # Drug disposition. Straube 2025 Table 2, Omalizumab / "Figure 5a" column.
    lvc     <- log(3.925);        label("Central volume of distribution (L)")                    # Table 2, Figure 5a column: Vc = 3.925 L (estimated patient-specific)
    lcl     <- log(0.01);         label("Systemic clearance (L/day)")                            # Table 2, Figure 5a column: CL = 0.01 L/day (estimated patient-specific); Table 2 reports the derived keD = CL/Vc = 0.0025 1/day
    lka     <- log(1.03);         label("First-order absorption rate constant (1/day)")          # Table 2, Omalizumab: ka = 1.03 1/day (pooled across the two patients)
    lfdepot <- fixed(log(0.42));  label("Subcutaneous bioavailability (fraction)")                # Table 2 footnote *: "F and Kd fixed at values reported in Stein and Ramakrishna" -> F = 0.42

    # Drug-target binding. Kd is fixed from an external source and koff is
    # estimated, so kon = koff / Kd is derived in model() (Equation 7).
    lk2     <- log(3.288);        label("Dissociation (off) rate constant of drug-target binding (1/day)")  # Table 2, Omalizumab: koff = 3.288 1/day (pooled)
    lkd     <- fixed(log(2.3));   label("Equilibrium dissociation constant (Kd = koff/kon, nM)")     # Table 2 footnote *: Kd fixed at 2.3 nM from Stein and Ramakrishna. Table 2 reports the derived kon = koff/Kd = 1.430 1/(nM*day)

    # Target turnover. Figure 1: free target is synthesised at ksyn = keR*Rb
    # and eliminated with rate constant keR, so the baseline is Rb = ksyn/keR.
    lrbase  <- log(1.342);        label("Baseline (basal) free target concentration (nM)")       # Table 2, Figure 5a column: Rb = 1.342 nM (estimated patient-specific); Table 2 reports the derived ksyn = keR*Rb = 1.130 nM/day
    lkdeg   <- log(0.842);        label("Free target elimination rate constant (1/day)")        # Table 2, Omalizumab: keR = 0.842 1/day (pooled)
    lkint   <- log(0.167);        label("Drug-target complex elimination rate constant (1/day)")  # Table 2, Omalizumab: keDR = keR/Tacc = 0.842/5.058 = 0.167 1/day. Tacc = 5.058 is the target accumulation ratio keR/keDR (Equation 14); it is recoverable from this file as kdeg/kint and so is not carried as a separate parameter
  })

  model({
    # Individual (typical-value) parameters. This is a deterministic
    # individual fit: Straube 2025 reports no IIV and no residual error.
    vc    <- exp(lvc)
    cl    <- exp(lcl)
    ka    <- exp(lka)
    k2    <- exp(lk2)
    kd    <- exp(lkd)
    rbase <- exp(lrbase)
    kdeg  <- exp(lkdeg)
    kint  <- exp(lkint)

    # Derived constants, exactly as Table 2 and Figure 1 define them.
    kel  <- cl / vc      # keD = CL/Vc
    k1   <- k2 / kd      # kon = koff/Kd (Equation 7)
    ksyn <- kdeg * rbase # ksyn = keR*Rb (Figure 1)

    # Mager-Jusko TMDD system, Straube 2025 Equation (2), with the
    # subcutaneous depot of Equation (S20). The paper writes Equation (2) in
    # concentrations (D, R, DR in nM); the states below hold the corresponding
    # AMOUNTS in nmol within the central volume vc, so every concentration is
    # state/vc and each bimolecular term carries one 1/vc.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - k1 * central * target / vc + k2 * complex - kel * central
    d/dt(target)  <- -k1 * central * target / vc + k2 * complex + ksyn * vc - kdeg * target
    d/dt(complex) <-  k1 * central * target / vc - k2 * complex - kint * complex

    # Equation (3): the target starts at its drug-free steady state Rb.
    target(0) <- rbase * vc

    # Straube 2025 encodes bioavailability by scaling the dose, not by the
    # Equation (S20) depot-loss term; see the vignette Errata for the
    # numerical check that selects this form.
    f(depot) <- exp(lfdepot)

    # Observations. Figure 5a plots total drug (DT), total target (RT) and
    # free target (R); the FTBR is Equation (28).
    Cc          <- (central + complex) / vc
    totalTarget <- (target + complex) / vc
    freeTarget  <- target / vc
    ftbr        <- freeTarget / rbase
  })
}
