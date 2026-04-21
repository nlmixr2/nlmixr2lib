PK_1cmt_tmdd_qss <- function() {
  description <- "One-compartment TMDD archetype, quasi-steady-state (QSS) approximation (Gibiansky et al. 2008)"
  reference <- "Gibiansky L, Gibiansky E, Kakkar T, Ma P. Approximations of the target-mediated drug disposition model and identifiability of model parameters. J Pharmacokinet Pharmacodyn. 2008;35(5):573-591. doi:10.1007/s10928-008-9102-8"
  vignette <- "tmdd_archetypes"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects = NA_integer_,
    notes = paste(
      "Generic TMDD QSS archetype; initial values are plausible mAb-scale defaults, not fit to a specific dataset.",
      "Drug and target are carried as concentrations in mg/L so the QSS binding relation (Kss in mg/L) is",
      "unit-consistent without a molecular-weight conversion. The central-compartment state holds the total",
      "drug amount (Atot = free + bound); the observed Cc is total drug concentration."
    )
  )

  ini({
    # Drug disposition
    lka     <- log(0.3);  label("Absorption rate (Ka, 1/day)")                                # Generic mAb-scale default (cf. Davda 2014 Table 3)
    lfdepot <- log(0.7);  label("Extravascular bioavailability (F, fraction)")                # Generic mAb-scale default
    lcl     <- log(0.2);  label("Linear (non-specific) clearance (CL, L/day)")                # Gibiansky 2008 Eq 8 (k_el * V); generic mAb-scale default
    lvc     <- log(3);    label("Central volume of distribution (Vc, L)")                     # Gibiansky 2008 Eq 8 (V); generic mAb-scale default

    # Target turnover
    lT0     <- log(0.1);  label("Baseline total target concentration (T0, mg/L)")             # Gibiansky 2008 Eq 9 (R_tot(0) = ksyn/kdeg); generic TMDD default
    lkdeg   <- log(0.1);  label("Free target first-order degradation rate (kdeg, 1/day)")     # Gibiansky 2008 Eq 9 (k_deg)
    lkint   <- log(1);    label("Drug-target complex internalization rate (kint, 1/day)")     # Gibiansky 2008 Eq 9 (k_int)

    # Binding (QSS)
    lKss    <- log(1.1);  label("Steady-state binding constant (Kss = (koff + kint)/kon, mg/L)")  # Gibiansky 2008 Eq 7 (Kss)

    # IIV (generic archetype; single eta on each of CL, Vc, Ka)
    etalcl ~ 0.09
    etalvc ~ 0.09
    etalka ~ 0.25

    propSd <- 0.15; label("Proportional residual error (fraction)")
  })

  model({
    # Individual drug-disposition parameters
    ka   <- exp(lka + etalka)
    fdep <- exp(lfdepot)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    kel  <- cl / vc

    # Target-turnover and binding parameters
    t0   <- exp(lT0)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    kss  <- exp(lKss)
    ksyn <- kdeg * t0  # steady-state relation: total_target(0) = T0

    # Initial conditions
    total_target(0) <- t0

    # QSS relation (Gibiansky 2008 Eq 7): solve Ctot = Cfree + complex, where
    #   complex = total_target * Cfree / (Kss + Cfree).
    # central is the total-drug amount in the central compartment (Atot),
    # so ctot = Atot / Vc.
    ctot    <- central / vc
    disc    <- ctot - total_target - kss
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # Reported concentration is total drug (common mAb assay).
    Cc <- ctot

    # Gibiansky 2008 Eq 8-9 (QSS state equations)
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * cfree * vc - kint * complex * vc
    d/dt(total_target) <-  ksyn - kdeg * (total_target - complex) - kint * complex

    f(depot) <- fdep

    Cc ~ prop(propSd)
  })
}
