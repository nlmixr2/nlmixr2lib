PK_2cmt_tmdd_qss <- function() {
  description <- "Two-compartment TMDD archetype, quasi-steady-state (QSS) approximation (Gibiansky et al. 2008)"
  reference <- "Gibiansky L, Gibiansky E, Kakkar T, Ma P. Approximations of the target-mediated drug disposition model and identifiability of model parameters. J Pharmacokinet Pharmacodyn. 2008;35(5):573-591. doi:10.1007/s10928-008-9102-8"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects = NA_integer_,
    notes = paste(
      "Generic two-compartment TMDD QSS archetype; initial values are plausible mAb-scale defaults, not fit",
      "to a specific dataset. Drug distributes between central and peripheral compartments; target turnover",
      "and binding occur in the central compartment only. The central state holds the total drug amount",
      "(free + bound) and Cc is the reported total drug concentration."
    )
  )

  ini({
    # Drug disposition
    lka     <- log(0.3);   label("Absorption rate (Ka, 1/day)")                                # Generic mAb-scale default (cf. Davda 2014 Table 3)
    lfdepot <- log(0.7);   label("Extravascular bioavailability (F, fraction)")                # Generic mAb-scale default
    lcl     <- log(0.2);   label("Linear (non-specific) clearance (CL, L/day)")                # Gibiansky 2008 Eq 8 (k_el * V); generic mAb-scale default
    lvc     <- log(3);     label("Central volume of distribution (Vc, L)")                     # Gibiansky 2008 Eq 8 (V); generic mAb-scale default
    lvp     <- log(3);     label("Peripheral volume of distribution (Vp, L)")                  # Gibiansky 2008 two-compartment extension; generic mAb-scale default
    lq      <- log(0.5);   label("Intercompartmental clearance (Q, L/day)")                    # Gibiansky 2008 two-compartment extension; generic mAb-scale default

    # Target turnover
    lT0     <- log(0.1);   label("Baseline total target concentration (T0, mg/L)")             # Gibiansky 2008 Eq 9 (R_tot(0) = ksyn/kdeg); generic TMDD default
    lkdeg   <- log(0.1);   label("Free target first-order degradation rate (kdeg, 1/day)")     # Gibiansky 2008 Eq 9 (k_deg)
    lkint   <- log(1);     label("Drug-target complex internalization rate (kint, 1/day)")     # Gibiansky 2008 Eq 9 (k_int)

    # Binding (QSS)
    lKss    <- log(1.1);   label("Steady-state binding constant (Kss = (koff + kint)/kon, mg/L)")  # Gibiansky 2008 Eq 7 (Kss)

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
    vp   <- exp(lvp)
    q    <- exp(lq)
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp

    # Target-turnover and binding parameters
    t0   <- exp(lT0)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    kss  <- exp(lKss)
    ksyn <- kdeg * t0

    total_target(0) <- t0

    # QSS relation in the central compartment.
    ctot    <- central / vc
    disc    <- ctot - total_target - kss
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # Reported concentration is total drug.
    Cc <- ctot

    # Gibiansky 2008 Eq 8-9 extended with linear peripheral distribution on free drug.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * cfree * vc - kint * complex * vc - k12 * cfree * vc + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * cfree * vc - k21 * peripheral1
    d/dt(total_target) <-  ksyn - kdeg * (total_target - complex) - kint * complex

    f(depot) <- fdep

    Cc ~ prop(propSd)
  })
}
