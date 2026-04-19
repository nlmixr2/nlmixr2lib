PK_1cmt_tmdd_full <- function() {
  description <- "One-compartment TMDD archetype with explicit drug-target binding (Mager & Jusko 2001 full model)"
  reference <- "Mager DE, Jusko WJ. General pharmacokinetic model for drugs exhibiting target-mediated drug disposition. J Pharmacokinet Pharmacodyn. 2001;28(6):507-532. doi:10.1023/A:1014414520282"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects = NA_integer_,
    notes = paste(
      "Generic TMDD archetype; initial values are plausible mAb-scale defaults, not fit to a specific dataset.",
      "Free drug, free target, and drug-target complex are all carried as concentrations in mg/L so that kon",
      "has consistent units (L/(mg*day)) without a molecular-weight conversion. Users who prefer molar units",
      "(nM) can re-interpret all concentration-bearing parameters (lT0, lkon, propSd) in that unit system."
    )
  )

  ini({
    # Drug disposition
    lka     <- log(0.3);  label("Absorption rate (Ka, 1/day)")                                # Mager & Jusko 2001 Fig 2 (generic mAb-scale default; cf. Davda 2014 Table 3 Ka = 0.282)
    lfdepot <- log(0.7);  label("Extravascular bioavailability (F, fraction)")                # Generic mAb-scale default (cf. Davda 2014 Table 3 F = 0.744)
    lcl     <- log(0.2);  label("Linear (non-specific) clearance (CL, L/day)")                # Mager & Jusko 2001 Eq 1 (k_el * V); generic mAb-scale default
    lvc     <- log(3);    label("Central volume of distribution (Vc, L)")                     # Mager & Jusko 2001 Eq 1 (V); generic mAb-scale default

    # Target turnover
    lT0     <- log(0.1);  label("Baseline free target concentration (T0, mg/L)")              # Mager & Jusko 2001 Eq 2 (R0 = ksyn/kdeg); generic TMDD default
    lkdeg   <- log(0.1);  label("Free target first-order degradation rate (kdeg, 1/day)")     # Mager & Jusko 2001 Eq 2 (k_deg)
    lkint   <- log(1);    label("Drug-target complex internalization rate (kint, 1/day)")     # Mager & Jusko 2001 Eq 3 (k_int)

    # Binding
    lkon    <- log(1);    label("Association rate constant (kon, L/(mg*day))")                # Mager & Jusko 2001 Eq 1 (k_on)
    lkoff   <- log(0.1);  label("Dissociation rate constant (koff, 1/day)")                   # Mager & Jusko 2001 Eq 1 (k_off)

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
    kon  <- exp(lkon)
    koff <- exp(lkoff)
    ksyn <- kdeg * t0  # steady-state relation: target(0) = ksyn/kdeg = T0

    # Initial conditions: target at baseline, no complex
    target(0)  <- t0
    complex(0) <- 0

    # Free drug concentration (central is amount in mg; vc in L; Cc in mg/L)
    Cc <- central / vc

    # Mager & Jusko 2001 Eq 1-3; drug/complex binding written in concentration
    # units then scaled by vc so that the central-compartment ODE is in mg/day.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central - kon * Cc * target * vc + koff * complex * vc
    d/dt(target)  <-  ksyn - kdeg * target - kon * Cc * target + koff * complex
    d/dt(complex) <-  kon * Cc * target - (koff + kint) * complex

    f(depot) <- fdep

    Cc ~ prop(propSd)
  })
}
