Mohamed_2012_gentamicin <- function() {
  description <- "In vitro (Escherichia coli ATCC 25922). Semi-mechanistic PKPD model of gentamicin bactericidal activity with adaptive resistance: drug-susceptible growing bacteria (bact_growing) plus insusceptible resting bacteria (bact_resting), with a binding model (ar_off / ar_on) by which gentamicin reduces its own Emax. Fit jointly to static and dynamic in-vitro time-kill curves."
  reference <- "Mohamed AF, Nielsen EI, Cars O, Friberg LE. Pharmacokinetic-pharmacodynamic model for gentamicin and its adaptive resistance with predictions of dosing schedules in newborn infants. Antimicrob Agents Chemother. 2012 Jan;56(1):179-188. doi:10.1128/AAC.00694-11. Model differential equations (Eqs 1-7), Figure 1 schematic, and final-model parameter estimates (Table 1) are in the main text."
  vignette <- "Mohamed_2012_gentamicin"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  # paper-mechanistic states (bacterial populations, adaptive-resistance
  # binding states, and the in-vitro flask gentamicin concentration);
  # none map onto canonical PK compartments.
  paper_specific_compartment_pattern <- "^bact_|^ar_|^cgent$"

  # No patient covariates: this is an in-vitro mechanism-based PKPD model.
  # The gentamicin exposure is a state variable (cgent) dosed by the user,
  # not a covariate column.
  covariateData <- list()

  population <- list(
    species          = "in vitro (Escherichia coli ATCC 25922; MIC 2 mg/L by macrodilution)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Gram-negative neonatal infection (in-vitro time-kill experiments; linked downstream to a published 3-compartment neonatal popPK of gentamicin for dosing predictions)",
    model_system     = "In-vitro time-kill curve experiments in cation-adjusted Mueller-Hinton broth at 35 C (48 static experiments + 25 dynamic experiments in a 2-compartment kinetic flask)",
    initial_inoculum = "approximately 5e5 CFU/mL; mean across experiments 4.83e5 CFU/mL (high-inoculum subset pre-grown 12 h to ~1e9 CFU/mL)",
    dose_range       = "Static gentamicin 0.125-16 mg/L; dynamic peak concentrations 2.0, 3.9, 7.8, 16 mg/L (simulating neonatal 1, 2, 4, 8 mg/kg doses) at 6-, 12-, or 24-h intervals",
    n_observations   = 1695L,
    notes            = paste(
      "Semi-mechanistic PKPD model fitted to bacterial counts (natural log of CFU/mL) jointly from static and dynamic in-vitro time-kill experiments using the NONMEM Laplacian method with ADVAN9; bacterial counts below the limit of detection (10 CFU/mL) were handled with the M3 method.",
      "kdeath was fixed to 0.179 /h from the prior Nielsen et al. semi-mechanistic antibiotic model (Nielsen 2007 / ref 36 in the paper); the data did not strongly identify it.",
      "koff was fixed to 0.0139 /h (50-h half-life of return to susceptibility), the lowest value that did not worsen model fit; lower values increased OFV.",
      "Three residual-error magnitudes were estimated: RE_static = 1.69 ln-units, RE_dynamic = 2.80 ln-units, and a replicate-specific RRE = 0.618 ln-units; the model file ships RE_static as the default residual error (see vignette Assumptions for the other two).",
      "kgrowth = 2.00 /h gives a 21-min mean generation time (60*ln(2)/kgrowth) consistent with E. coli in broth.",
      "Inter-experiment variability was not estimated, so the model has no eta/IIV terms."
    )
  )

  ini({
    # --- Bacterial population dynamics (Table 1) -----------------------------
    kgrowth <- 2.00          ; label("Rate constant of bacterial growth (1/h)")                              # Table 1: k_growth = 2.00 (95% CI 1.89-2.27)
    kdeath  <- fixed(0.179)  ; label("Rate constant of natural bacterial death (1/h; FIXED)")               # Table 1: k_death = 0.179 (fix; from Nielsen 2007 / ref 36)
    bp      <- 2.09e6        ; label("Breakpoint total bacteria for k_SR > 0 (CFU/mL)")                     # Table 1: BP = 2.09e6 (95% CI 1.01e6-3.32e6)
    bmax    <- 8.26e8        ; label("Stationary-phase total bacterial count B_max (CFU/mL)")               # Table 1: B_max = 8.26e8 (95% CI 6.18e8-11.10e8)

    # --- Gentamicin Emax model on the susceptible compartment (Table 1) ------
    emax0   <- 51.0          ; label("Maximum gentamicin kill rate constant at zero adaptive resistance (1/h)") # Table 1: E_max(0) = 51.0 (95% CI 44.6-61.6)
    ec50    <- 9.93          ; label("Gentamicin concentration for 50% of E_max(0) (mg/L)")                 # Table 1: EC_50 = 9.93 (95% CI 8.45-12.1)

    # --- Adaptive-resistance binding model (Table 1) -------------------------
    ar50    <- 0.113         ; label("AR_on value at which E_max is reduced by 50% (unitless)")             # Table 1: AR_50 = 0.113 (95% CI 0.0983-0.146)
    kon     <- 0.0426        ; label("Rate constant for development of adaptive resistance (L/(mg*h))")     # Table 1: k_on = 0.0426 (95% CI 0.0376-0.0478)
    koff    <- fixed(0.0139) ; label("Rate constant for return to susceptibility (1/h; FIXED at ~50 h half-life)") # Table 1: k_off = 0.0139 (fix; ln(2)/0.0139 ~= 49.9 h)

    # --- Initial bacterial inoculum (typical value used in the paper's predictions) ---
    s0      <- 4.83e5        ; label("Initial susceptible-bacteria inoculum S(0) (CFU/mL)")                 # Methods: avg of all experimental start inocula 4.83e5; predictions used this value

    # --- In-vitro flask gentamicin elimination (input, not estimated) --------
    # The paper used three values depending on experimental setting:
    #   0       /h for static experiments,
    #   0.33    /h for the first 4 h of dynamic experiments,
    #   0.037   /h thereafter in dynamic experiments,
    # all parameterising the in-vitro kinetic system that mimicked
    # preterm-neonate gentamicin PK. The default below is the dynamic
    # phase-2 (terminal) value; override `lkel` at rxSolve() time for
    # static experiments (e.g. `params = c(lkel = log(1e-12))`) or
    # split the simulation around the 4-h flow-rate switch (see vignette).
    lkel    <- fixed(log(0.037)) ; label("Log first-order gentamicin elimination from in-vitro flask (1/h; FIXED, dynamic phase-2 default)") # Methods: ke = 0 (static); 0.33 then 0.037 /h (dynamic in-vitro kinetic system, switched at 4 h)

    # --- Residual error (additive on natural-log CFU/mL scale) ---------------
    addSd   <- 1.69          ; label("Additive residual SD on ln(CFU/mL) scale (static experiments; see vignette for RE_dynamic and RRE)") # Table 1: RE_static = 1.69 (95% CI 1.52-1.83); RE_dynamic = 2.80 and RRE = 0.618 also reported but not encoded here
  })

  model({
    # 1. Adaptive-resistance attenuation of the maximum killing rate (Eq 4).
    #    Inhibition fraction by AR is ar_on/(ar_50 + ar_on); preserved
    #    fraction is ar_50/(ar_50 + ar_on). Mass balance keeps
    #    ar_off + ar_on = 1, so ar_on saturates at 1 and the maximum
    #    inhibition is 1/(1 + ar_50) = 0.898 (~ 90%; paper Results).
    emax <- emax0 * ar50 / (ar50 + ar_on)

    # 2. Gentamicin killing rate constant on the susceptible compartment
    #    (Eq 3; Hill = 1 -- paper showed no improvement from estimating gamma).
    drug <- emax * cgent / (ec50 + cgent)

    # 3. Transfer rate constant S -> R (Eq 1 narrative and Methods).
    #    k_SR = beta * (S + R) when total bacteria exceed the breakpoint,
    #    else 0. beta is derived from the steady-state identity
    #    B_max = (k_growth - k_death) / beta.
    total_bact <- bact_growing + bact_resting
    beta <- (kgrowth - kdeath) / bmax
    ksr  <- beta * total_bact * (total_bact > bp)

    # 4. Bacterial population ODEs (Eqs 1 + 5: drug acts only on S; no
    #    transfer R -> S, paper fixed k_RS = 0).
    d/dt(bact_growing) <- (kgrowth - kdeath - ksr - drug) * bact_growing
    d/dt(bact_resting) <- ksr * bact_growing - kdeath * bact_resting

    # 5. Adaptive-resistance binding model (Eqs 6-7). Gentamicin drives
    #    transfer from ar_off to ar_on at rate kon * cgent * ar_off; the
    #    reverse first-order rate koff drives return to susceptibility.
    d/dt(ar_off) <- -kon * cgent * ar_off + koff * ar_on
    d/dt(ar_on)  <- +kon * cgent * ar_off - koff * ar_on

    # 6. Gentamicin concentration in the in-vitro flask (paper Methods).
    #    The user adds gentamicin to cgent via dosing events; kel is set
    #    by the experiment type (0 for static; 0.33 then 0.037 /h for
    #    the dynamic in-vitro kinetic system at the 4-h phase switch).
    kel <- exp(lkel)
    d/dt(cgent) <- -kel * cgent

    # 7. Initial conditions: all bacteria start in the susceptible/growing
    #    compartment at the typical inoculum; adaptive resistance starts
    #    in the off state (paper Methods: ar_off(0) = 1, ar_on(0) = 0).
    bact_growing(0) <- s0
    bact_resting(0) <- 0
    ar_off(0)       <- 1
    ar_on(0)        <- 0

    # 8. Observation: natural log of total viable count (the scale on which
    #    the paper transformed data and fitted residual error).
    Cc <- log(bact_growing + bact_resting)
    Cc ~ add(addSd)
  })
}
