PK_2cmt_tmdd_mm <- function() {
  description <- "Two-compartment TMDD archetype, Michaelis-Menten (MM) approximation (Gibiansky et al. 2008)"
  reference <- "Gibiansky L, Gibiansky E, Kakkar T, Ma P. Approximations of the target-mediated drug disposition model and identifiability of model parameters. J Pharmacokinet Pharmacodyn. 2008;35(5):573-591. doi:10.1007/s10928-008-9102-8"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects = NA_integer_,
    notes = paste(
      "Generic two-compartment TMDD Michaelis-Menten archetype; initial values are plausible mAb-scale",
      "defaults, not fit to a specific dataset. Drug distributes between central and peripheral compartments",
      "with saturable target-mediated elimination from central (Vm, Km) plus linear clearance. The MM form",
      "is valid when total target changes slowly relative to drug disposition."
    )
  )

  ini({
    # Drug disposition
    lka     <- log(0.3);  label("Absorption rate (Ka, 1/day)")                              # Generic mAb-scale default (cf. Davda 2014 Table 3)
    lfdepot <- log(0.7);  label("Extravascular bioavailability (F, fraction)")              # Generic mAb-scale default
    lcl     <- log(0.2);  label("Linear (non-specific) clearance (CL, L/day)")              # Gibiansky 2008 Eq 10 (k_el * V); generic mAb-scale default
    lvc     <- log(3);    label("Central volume of distribution (Vc, L)")                   # Gibiansky 2008 Eq 10 (V); generic mAb-scale default
    lvp     <- log(3);    label("Peripheral volume of distribution (Vp, L)")                # Gibiansky 2008 two-compartment extension; generic mAb-scale default
    lq      <- log(0.5);  label("Intercompartmental clearance (Q, L/day)")                  # Gibiansky 2008 two-compartment extension; generic mAb-scale default

    # Saturable (target-mediated) elimination
    lVm     <- log(0.1);  label("Maximum target-mediated elimination rate (Vm, mg/L/day)")  # Gibiansky 2008 Eq 10 (V_m = k_int * R0); generic TMDD default
    lKm     <- log(1.1);  label("MM constant for target-mediated elimination (Km, mg/L)")   # Gibiansky 2008 Eq 10 (K_m approx Kss)

    # IIV (generic archetype; single eta on each of CL, Vc, Ka)
    etalcl ~ 0.09
    etalvc ~ 0.09
    etalka ~ 0.25

    propSd <- 0.15; label("Proportional residual error (fraction)")
  })

  model({
    # Individual parameters
    ka   <- exp(lka + etalka)
    fdep <- exp(lfdepot)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    vm   <- exp(lVm)
    km   <- exp(lKm)

    # central holds total drug amount; Cc is concentration in mg/L.
    Cc <- central / vc

    # Gibiansky 2008 Eq 10 with a linear peripheral distribution compartment.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - vm * Cc * vc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdep

    Cc ~ prop(propSd)
  })
}
