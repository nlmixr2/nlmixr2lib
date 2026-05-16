Kimura_2023_glucarpidase <- function() {
  description <- "Modified Michaelis-Menten PK/PD simulation model for glucarpidase (CPG2) rescue after high-dose methotrexate (Kimura 2023). MTX disposition is 2-compartment IV with renal-only first-order elimination (Kr fixed at ~10% of literature total MTX CL from Fukahara 2008); the remaining elimination is captured by a saturable hydrolysis term coupled to a 1-compartment IV CPG2 disposition. All structural parameters are literature-sourced point values (no estimation in the source paper)."
  reference <- "Kimura T, Fukaya Y, Hamada Y, Yoshimura K, Kawamoto H. Pharmacokinetics and Pharmacodynamics of Glucarpidase Rescue Treatment After High-dose Methotrexate Therapy Based on Modeling and Simulation. Anticancer Res. 2023;43(5):1919-1924. doi:10.21873/anticanres.16351"
  vignette <- "Kimura_2023_glucarpidase"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 500L,
    n_studies      = 0L,
    age_range      = "adult model parameters; planned phase II included pediatric patients",
    weight_range   = "60 kg (assumed for all virtual subjects)",
    sex_female_pct = NA_real_,
    disease_state  = "Patients receiving high-dose methotrexate (1 g/m^2 over 4 h) with subsequent glucarpidase rescue.",
    dose_range     = "MTX 1 g/m^2 IV over 4 h; CPG2 10-80 U/kg IV bolus (50 U/kg recommended).",
    regions        = "Japan (planned phase II study; JMA-IIA00097).",
    notes          = "Simulation-only study: PK parameters were taken from prior popPK / NCA publications (Fukahara 2008 for MTX, Table I ref 3; Phillips 2008 for CPG2, Table I ref 16) and the EMA Voraxaze assessment 2021 for Km and alpha (Table I ref 18). Monte-Carlo simulations of 500 virtual patients per arm sampled each PK parameter independently from the Normal distributions in Table II. Virtual subjects had body weight 60 kg, BSA 1.73 m^2, creatinine clearance 80 ml/min, and the MTX renal-only clearance fraction was 10% of the literature total CL."
  )

  ini({
    # Methotrexate (MTX) disposition (2-compartment, IV). All values from
    # Kimura 2023 Table II (Mean column) -- the parameter point estimates
    # that were used in the published Monte-Carlo runs. The MTX 'Clearance'
    # of Table II corresponds to Kr * Vc (renal-only clearance, ~10% of the
    # literature total CL = 5.57 L/h at CLcr=80 ml/min from Fukahara 2008,
    # Table I ref 3); the remaining elimination is captured by the
    # saturable CPG2-mediated hydrolysis term in model().
    lcl <- fixed(log(0.615));  label("MTX renal clearance (Kr*Vc, L/h)")            # Kimura 2023 Table II
    lvc <- fixed(log(26.683)); label("MTX central volume of distribution (L)")      # Kimura 2023 Table II
    lvp <- fixed(log(2.253));  label("MTX peripheral volume of distribution (L)")   # Kimura 2023 Table II
    lq  <- fixed(log(0.078));  label("MTX intercompartmental clearance (L/h)")      # Kimura 2023 Table II

    # Glucarpidase (CPG2) disposition (1-compartment, IV bolus). Values from
    # Kimura 2023 Table II (Mean column); the body-weight-normalized source
    # values in Table I ref 16 (Phillips 2008) become 0.310 L/h and 4.114 L
    # at the 60 kg virtual-patient weight used in the simulation.
    lcl_cpg2 <- fixed(log(0.310)); label("CPG2 clearance (L/h)")                    # Kimura 2023 Table II
    lvc_cpg2 <- fixed(log(4.114)); label("CPG2 volume of distribution at steady state (L)") # Kimura 2023 Table II

    # Michaelis-Menten coupling (MTX hydrolysis by CPG2) -- Kimura 2023
    # Equations 2-4. Both constants are taken from Table I ref 18 (EMA
    # Voraxaze assessment 2021) and are held fixed across all simulations.
    # alpha is defined by Vmax = alpha * [catalyst]; the paper reports
    # Vmax = 800 mol/s at 1 mol/L of CPG2, giving alpha = 800 L/s =
    # 2.88e6 L/h. With concentrations in umol/L the same numeric value of
    # alpha applies (the L/h units are independent of the mole-prefix).
    km_cpg2 <- fixed(86);      label("Michaelis-Menten constant for MTX hydrolysis by CPG2 (umol/L)") # Kimura 2023 Table I ref 18
    alpha   <- fixed(2.88e6);  label("MM conversion constant: Vmax = alpha * [CPG2] (L/h)")          # Kimura 2023 Table I ref 18 and Equation 3
  })

  model({
    # Back-transform structural parameters (all log-fixed).
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    cl_cpg2 <- exp(lcl_cpg2)
    vc_cpg2 <- exp(lvc_cpg2)

    # Micro-rate constants.
    kr      <- cl / vc          # MTX renal elimination rate constant (1/h)
    k12     <- q  / vc
    k21     <- q  / vp
    ke_cpg2 <- cl_cpg2 / vc_cpg2

    # Concentrations (amounts in umol; volumes in L; concentrations in umol/L).
    Cc      <- central      / vc
    Cc_cpg2 <- central_cpg2 / vc_cpg2

    # Saturable MM degradation of MTX by CPG2 (Kimura 2023 Equation 4).
    # Units: alpha (L/h) * Cc_cpg2 (umol/L) * dimensionless -> umol/h.
    mm_rate <- alpha * Cc_cpg2 * Cc / (km_cpg2 + Cc)

    # ODE system.
    d/dt(central)      <- -(kr + k12) * central + k21 * peripheral1 - mm_rate
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(central_cpg2) <- -ke_cpg2 * central_cpg2
  })
}
