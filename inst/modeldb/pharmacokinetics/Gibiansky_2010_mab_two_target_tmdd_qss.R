Gibiansky_2010_mab_two_target_tmdd_qss <- function() {
  description <- "Methodology reference. Two-target quasi-steady-state (QSS) approximation of target-mediated drug disposition (TMDD) for a hypothetical fully human monoclonal antibody that binds a soluble target (S) and a membrane-bound target (M) (Gibiansky 2010 Eqs 12-17, model M3). The soluble target is described by the QSS approximation with a dynamic total-target state, while M-target-mediated elimination reduces to a Michaelis-Menten term because the total membrane-target concentration is assumed constant and small relative to the free drug. Drug disposition is two-compartment with linear elimination, i.v. bolus and first-order s.c. absorption; the central state holds the total (free + S-bound) drug and the free drug is recovered algebraically (Eq 12). Cc is the free drug concentration and Ctotal_target_soluble the total soluble-target concentration, both nmol/L. Defaults are the true values of the authors' simulation (Table 1, which the paper also lists as the true values for M3 in Table 5) with the assumed simulation variability; no drug is named and no patient data were fitted. Sibling: modellib('Gibiansky_2010_mab_two_target_tmdd_full'), the full two-target TMDD model this approximates."
  reference <- paste(
    "Gibiansky L, Gibiansky E. Target-mediated drug disposition model for",
    "drugs that bind to more than one target. J Pharmacokinet Pharmacodyn.",
    "2010;37(4):323-346. doi:10.1007/s10928-010-9163-3.",
    "Structural equations are Eqs 12-17 (model M3); parameter values are the",
    "Table 1 true values, repeated as the 'True value' column of Table 5.",
    sep = " "
  )
  vignette <- "Gibiansky_2010_mab_two_target_tmdd"

  units <- list(time = "day", dosing = "nmol", concentration = "nmol/L")

  # Spelled-out soluble-target identity, matching the S / M naming of the
  # sibling full model. The membrane target has no state here: Eqs 12-17
  # replace it with a Michaelis-Menten elimination term.
  paper_specific_compartments <- c("total_target_soluble")

  # `central` holds the AMOUNT Ctot * Vc of total (free + drug-S-target
  # complex) drug; the paper's Eq 14 is written for the concentration Ctot
  # and is multiplied through by Vc below.
  compartmentData <- list(
    depot = list(
      analyte = "monoclonal antibody (unnamed), s.c. depot",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "monoclonal antibody (unnamed), free plus soluble-target-bound",
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
    total_target_soluble = list(
      analyte = "soluble target S, free plus antibody-bound",
      units = "nmol/L",
      specimen = "serum",
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
      "The 224 subjects were SIMULATED from the full two-target TMDD model",
      "(Gibiansky_2010_mab_two_target_tmdd_full) and this QSS approximation (M3) was fitted to them;",
      "its Set 1 estimates (Table 4) were within 16% of the true values except KSS_S (29% low).",
      "This file carries the TRUE values, not the M3 estimates, as the paper's reference parameter set.",
      "Every value is therefore encoded with fixed()."
    )
  )

  ini({
    # ---- Linear part of the model (Gibiansky 2010 Table 1 / Table 5 true values) ----
    # Table 1 prints the units of CL and Q as 'day^-1'; they are clearances in
    # L/day, as the Table 1 footnote kel = CL/Vc, kpt = Q/Vc, ktp = Q/Vp requires.
    lcl <- fixed(log(0.3)); label("Linear clearance of free drug (L/day)") # Table 1 and Table 5, CL = 0.3
    lvc <- fixed(log(3)); label("Central volume of distribution (L)") # Table 1 and Table 5, Vc = 3
    lq <- fixed(log(0.2)); label("Intercompartmental clearance (L/day)") # Table 1 and Table 5, Q = 0.2
    lvp <- fixed(log(3)); label("Peripheral volume of distribution (L)") # Table 1 and Table 5, Vp = 3
    lfdepot <- fixed(log(0.7)); label("Bioavailability of the s.c. dose (fraction)") # Table 1 and Table 5, F_SC = 0.7
    lka <- fixed(log(0.5)); label("First-order s.c. absorption rate constant (1/day)") # Table 1 and Table 5, ka = 0.5

    # ---- Soluble target S, QSS parameterisation of M3 (Table 4 / Table 5 rows) ----
    lkint_s <- fixed(log(0.05)); label("Elimination (internalisation) rate constant of the drug-soluble target complex (1/day)") # Table 1 and Table 5, kint_S = 0.05
    lkdeg_s <- fixed(log(10)); label("First-order degradation rate constant of the free soluble target (1/day)") # Table 1 and Table 5, kdeg_S = 10
    lrbase_s <- fixed(log(0.1)); label("Baseline soluble-target concentration R0_S = ksyn_S/kdeg_S (nmol/L)") # Table 1 and Table 5, R0_S = 0.1 (derived from ksyn_S = 1, kdeg_S = 10)
    lkss_s <- fixed(log(0.015)); label("QSS constant of the soluble target KSS_S = (koff_S + kint_S)/kon_S (nmol/L)") # Table 1 and Table 5, KSS_S = 0.015 (derived)

    # ---- Membrane target M, Michaelis-Menten term (Eq 14, Eq 17) ----
    lvmax_m <- fixed(log(1.5)); label("Maximum elimination rate via the membrane target Vmax_M = R0_M * kint_M (nmol/L/day)") # Table 1 and Table 5, Vmax_M = 1.5 (derived)
    lkss_m <- fixed(log(3.05)); label("QSS constant of the membrane target KSS_M = (koff_M + kint_M)/kon_M (nmol/L)") # Table 1 and Table 5, KSS_M = 3.05 (derived)

    # ---- Inter-individual variability (simulation assumption) ----
    # 20% CV log-normal IIV was assumed in the simulation 'for all the
    # parameters except SC bioavailability and the binding constants kon and
    # koff' (Simulation study: methods): omega^2 = log(1 + 0.2^2) = 0.0392207.
    # It is applied here to the M3 parameters that are themselves simulation
    # parameters (CL, Vc, Q, Vp, ka, kint_S, kdeg_S). R0_S, KSS_S, Vmax_M and
    # KSS_M are composites of simulation parameters and carry no IIV; the
    # paper does not report the M3 random-effect estimates.
    etalcl ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalvc ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalq ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalvp ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalka ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkint_s ~ fixed(0.0392207) # 20% CV, Simulation study: methods
    etalkdeg_s ~ fixed(0.0392207) # 20% CV, Simulation study: methods

    # ---- Residual error (simulation assumption) ----
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

    kint_s <- exp(lkint_s + etalkint_s)
    kdeg_s <- exp(lkdeg_s + etalkdeg_s)
    rbase_s <- exp(lrbase_s)
    kss_s <- exp(lkss_s)
    ksyn_s <- rbase_s * kdeg_s # Eq 17: R0_S = ksyn_S / kdeg_S

    vmax_m <- exp(lvmax_m)
    kss_m <- exp(lkss_m)

    # ---- Micro-constants (Table 1 footnote b) ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Free drug from total drug (Eq 12) ----
    ctot <- central / vc
    disc <- ctot - total_target_soluble - kss_s
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kss_s * ctot))

    # ---- ODE system (Eqs 13-16) ----
    # Eq 14 multiplied through by Vc: linear elimination and distribution act
    # on the FREE drug, the S-target term eliminates complex at kint_S, and
    # the M-target term is the Michaelis-Menten flux.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - (kel + k12) * cfree * vc -
      vc * total_target_soluble * kint_s * cfree / (kss_s + cfree) -
      vc * vmax_m * cfree / (kss_m + cfree)
    d/dt(peripheral1) <- k12 * cfree * vc - k21 * peripheral1
    d/dt(total_target_soluble) <- ksyn_s - kdeg_s * total_target_soluble -
      (kint_s - kdeg_s) * total_target_soluble * cfree / (kss_s + cfree)

    f(depot) <- fdepot

    # Eq 17: Rtot_S(0) = R0_S
    total_target_soluble(0) <- rbase_s

    # ---- Outputs ----
    Cc <- cfree
    Ctotal_target_soluble <- total_target_soluble

    Cc ~ lnorm(expSd)
    Ctotal_target_soluble ~ lnorm(expSd_Ctotal_target_soluble)
  })
}
