Gibiansky_2010_mab_two_target_tmdd_full <- function() {
  description <- "Methodology reference. Full two-target target-mediated drug disposition (TMDD) model for a hypothetical fully human monoclonal antibody that binds, one target at a time, a soluble target (S) and a membrane-bound target (M), each with its own reversible binding, complex internalisation and zero-order synthesis / first-order degradation (Gibiansky 2010 Eqs 1-3 with N = 2). Drug disposition is two-compartment with linear elimination from the central compartment, i.v. bolus and first-order s.c. absorption. Cc is the free (unbound) drug concentration and Ctotal_target_soluble the total (free + drug-bound) soluble-target concentration, both nmol/L. This is the 'true' model the authors simulated their 224-subject identifiability study from; defaults are the authors' own nominal parameter set (Table 1) with the assumed simulation variability (20% CV inter-individual variability on every parameter except F, kon and koff; 15% / 20% log-normal residual error). No drug is named and no patient data were fitted. Sibling: modellib('Gibiansky_2010_mab_two_target_tmdd_qss'), the two-target QSS (S) + Michaelis-Menten (M) approximation the paper recommends."
  reference <- paste(
    "Gibiansky L, Gibiansky E. Target-mediated drug disposition model for",
    "drugs that bind to more than one target. J Pharmacokinet Pharmacodyn.",
    "2010;37(4):323-346. doi:10.1007/s10928-010-9163-3.",
    "Structural equations are Eqs 1-3 (N = 2 targets); parameter values are",
    "Table 1; the simulation variability is stated in 'Simulation study: methods'.",
    sep = " "
  )
  vignette <- "Gibiansky_2010_mab_two_target_tmdd"

  # Time is days, doses are nmol and every concentration is nmol/L, because
  # the binding terms of Eqs 1-3 multiply drug by target concentration and
  # Table 1 quotes kon in L/nmol/day and ksyn in nmol/L/day.
  units <- list(time = "day", dosing = "nmol", concentration = "nmol/L")

  # The soluble and membrane targets are not standard library compartments:
  # the paper names them only S and M, so the states carry that identity in
  # spelled-out form. The registered generic pair target_a / target_b is
  # reserved for papers whose targets carry no identity at all.
  paper_specific_compartments <- c(
    "target_soluble",
    "complex_soluble",
    "target_membrane",
    "complex_membrane"
  )

  # The paper writes the drug states as Ad (depot amount), C (central
  # CONCENTRATION) and AT (tissue amount). Here `central` holds the AMOUNT
  # C * Vc so that an i.v. bolus enters as a dose amount; Eq 1 is multiplied
  # through by Vc below. The four target states keep the paper's
  # concentration units.
  compartmentData <- list(
    depot = list(
      analyte = "monoclonal antibody (unnamed), s.c. depot",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "monoclonal antibody (unnamed), free",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "monoclonal antibody (unnamed), free",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    target_soluble = list(analyte = "soluble target S, free", units = "nmol/L", specimen = "serum", verified = TRUE),
    complex_soluble = list(
      analyte = "antibody-soluble target S complex",
      units = "nmol/L",
      specimen = "serum",
      verified = TRUE
    ),
    target_membrane = list(
      analyte = "membrane-bound target M, free",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_membrane = list(
      analyte = "antibody-membrane target M complex",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  dosing <- c("central", "depot")

  covariateData <- list()

  population <- list(
    species = "None (methodology paper; hypothetical monoclonal antibody, simulated subjects, no fitted patient data)",
    n_subjects = 224L,
    n_studies = 2L,
    disease_state = "None. The parameters 'were chosen to imitate a monoclonal antibody that binds to the soluble (S) and membrane-bound (M) targets' (Abstract); no drug or indication is named.",
    dose_range = "Study 1: single 100, 300 or 600 nmol i.v. or 1000 nmol s.c. (6 subjects per arm); Study 2: 600 nmol i.v. or 1000 nmol s.c. on days 0, 28 and 56 (100 subjects per arm) (Table 2).",
    notes = paste(
      "All 224 subjects are SIMULATED from this model with the Table 1 parameters; the paper",
      "removed simulated values below 0.1 nmol/L, leaving 3,298 free-drug and 3,257 total",
      "S-target observations. The linear-PK and binding values are described as typical for fully",
      "human therapeutic antibodies and the target turnover values as consistent with the literature",
      "(Table 1, refs 20-27), but the specific numbers are the authors' own nominal choices.",
      "Every value is therefore encoded with fixed()."
    )
  )

  ini({
    # ---- Linear part of the model (Gibiansky 2010 Table 1) ---------------------
    # Table 1 prints the units of CL and Q as 'day^-1'; they are clearances in
    # L/day, as the Table 1 footnote kel = CL/Vc, kpt = Q/Vc, ktp = Q/Vp requires.
    lcl <- fixed(log(0.3)); label("Linear clearance of free drug (L/day)") # Table 1, CL = 0.3
    lvc <- fixed(log(3)); label("Central volume of distribution (L)") # Table 1, Vc = 3.0
    lq <- fixed(log(0.2)); label("Intercompartmental clearance (L/day)") # Table 1, Q = 0.2
    lvp <- fixed(log(3)); label("Peripheral volume of distribution (L)") # Table 1, Vp = 3.0
    lfdepot <- fixed(log(0.7)); label("Bioavailability of the s.c. dose (fraction)") # Table 1, F_SC = 0.7
    lka <- fixed(log(0.5)); label("First-order s.c. absorption rate constant (1/day)") # Table 1, ka = 0.5

    # ---- Soluble target S (Gibiansky 2010 Table 1) -----------------------------
    lkon_s <- fixed(log(10)); label("Association rate constant, drug + soluble target (L/nmol/day)") # Table 1, kon_S = 10
    lkoff_s <- fixed(log(0.1)); label("Dissociation rate constant of the drug-soluble target complex (1/day)") # Table 1, koff_S = 0.1
    lkint_s <- fixed(log(0.05)); label("Elimination (internalisation) rate constant of the drug-soluble target complex (1/day)") # Table 1, kint_S = 0.05
    lksyn_s <- fixed(log(1)); label("Zero-order synthesis rate of the soluble target (nmol/L/day)") # Table 1, ksyn_S = 1
    lkdeg_s <- fixed(log(10)); label("First-order degradation rate constant of the free soluble target (1/day)") # Table 1, kdeg_S = 10

    # ---- Membrane-bound target M (Gibiansky 2010 Table 1) ----------------------
    lkon_m <- fixed(log(5)); label("Association rate constant, drug + membrane target (L/nmol/day)") # Table 1, kon_M = 5
    lkoff_m <- fixed(log(0.25)); label("Dissociation rate constant of the drug-membrane target complex (1/day)") # Table 1, koff_M = 0.25
    lkint_m <- fixed(log(15)); label("Internalisation rate constant of the drug-membrane target complex (1/day)") # Table 1, kint_M = 15
    lksyn_m <- fixed(log(1.5)); label("Zero-order synthesis rate of the membrane target (nmol/L/day)") # Table 1, ksyn_M = 1.5
    lkdeg_m <- fixed(log(15)); label("First-order degradation rate constant of the free membrane target (1/day)") # Table 1, kdeg_M = 15

    # ---- Inter-individual variability (simulation assumption) ------------------
    # 'Moderate inter-individual variability (log-normal distributions with 20%
    # CV) was assumed for all the parameters except SC bioavailability and the
    # binding constants kon and koff' (Simulation study: methods). Variance of a
    # log-normal with 20% CV: omega^2 = log(1 + 0.2^2) = 0.0392207. No
    # correlation is stated.
    etalcl ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalvc ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalq ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalvp ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalka ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkint_s ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalksyn_s ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkdeg_s ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkint_m ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalksyn_m ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkdeg_m ~ fixed(0.0392207) # 20% CV, Simulation study: methods

    # ---- Residual error (simulation assumption) --------------------------------
    # 'lognormal distributions with 15 and 20% CV for the drug and target
    # concentrations, respectively' (Simulation study: methods); log-scale SD =
    # sqrt(log(1 + CV^2)).
    expSd <- fixed(0.149167); label("Log-normal residual error on free drug Cc (log-scale SD; 15% CV)") # Simulation study: methods, 15% CV
    expSd_Ctotal_target_soluble <- fixed(0.198042); label("Log-normal residual error on total soluble target (log-scale SD; 20% CV)") # Simulation study: methods, 20% CV
  })

  model({
    # ---- Individual parameters ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot)

    kon_s <- exp(lkon_s)
    koff_s <- exp(lkoff_s)
    kint_s <- exp(lkint_s + etalkint_s)
    ksyn_s <- exp(lksyn_s + etalksyn_s)
    kdeg_s <- exp(lkdeg_s + etalkdeg_s)

    kon_m <- exp(lkon_m)
    koff_m <- exp(lkoff_m)
    kint_m <- exp(lkint_m + etalkint_m)
    ksyn_m <- exp(lksyn_m + etalksyn_m)
    kdeg_m <- exp(lkdeg_m + etalkdeg_m)

    # ---- Micro-constants (Table 1 footnote b) ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Free drug concentration (the paper's C)
    cfree <- central / vc

    # ---- ODE system (Gibiansky 2010 Eqs 1-3, N = 2) ----
    # Eq 1, multiplied through by Vc: the F_SC * ka * Ad input is carried by
    # f(depot) below, and the binding terms are concentration fluxes scaled
    # to amounts by vc.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - (kel + k12) * central -
      vc * (kon_s * cfree * target_soluble - koff_s * complex_soluble) -
      vc * (kon_m * cfree * target_membrane - koff_m * complex_membrane)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Eq 2 (free targets)
    d/dt(target_soluble) <- ksyn_s - kdeg_s * target_soluble -
      kon_s * cfree * target_soluble + koff_s * complex_soluble
    d/dt(target_membrane) <- ksyn_m - kdeg_m * target_membrane -
      kon_m * cfree * target_membrane + koff_m * complex_membrane

    # Eq 3 (drug-target complexes)
    d/dt(complex_soluble) <- kon_s * cfree * target_soluble - (kint_s + koff_s) * complex_soluble
    d/dt(complex_membrane) <- kon_m * cfree * target_membrane - (kint_m + koff_m) * complex_membrane

    f(depot) <- fdepot

    # Eq 2 initial conditions: R_i(0) = R0_i = ksyn_i / kdeg_i; complexes start at 0.
    target_soluble(0) <- ksyn_s / kdeg_s
    target_membrane(0) <- ksyn_m / kdeg_m

    # ---- Outputs ----
    # The simulated observations are the free drug and the total (free +
    # drug-bound) S-target concentrations; the M-target was simulated but
    # treated as unobservable (Simulation study: methods).
    Cc <- cfree
    Ctotal_target_soluble <- target_soluble + complex_soluble

    Cc ~ lnorm(expSd)
    Ctotal_target_soluble ~ lnorm(expSd_Ctotal_target_soluble)
  })
}
