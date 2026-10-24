Dua_2015_tmdd_multitarget_qss <- function() {
  description <- "QSP. Two-compartment target-mediated drug disposition (TMDD) model for a drug that binds BOTH a soluble target and a membrane-bound target in the central compartment (Gibiansky & Gibiansky 2010 two-target approximation, as implemented in Dua 2015 Supplementary Model 6). The soluble target is handled by the quasi-steady-state (QSS) approximation with its own zero-order synthesis and first-order degradation; the membrane target is assumed constant in total and enters as a Michaelis-Menten elimination term with Vmax = Rtot_M * kint_M. Only free drug distributes to the tissue compartment and is eliminated linearly. Defaults are the tutorial's illustrative parameter set for a generic monoclonal antibody (i.v. bolus, mg/kg dosing converted to nmol/kg with MW 150 kDa); the model is not fitted to any dataset, so no drug and no real population is attached. Sibling in this library: modellib('Dua_2015_tmdd_fcrn_recycling')."
  reference <- paste(
    "Dua P, Hawkins E, van der Graaf PH. A Tutorial on Target-Mediated Drug",
    "Disposition (TMDD) Models. CPT Pharmacometrics Syst Pharmacol.",
    "2015;4(6):324-337. doi:10.1002/psp4.41.",
    "Equations are main-text Eqs V.1-V.6 ('Multiple targets'); parameter values",
    "are Supplementary Materials 'Berkeley Madonna TMDD model library', Model 6.",
    "Model structure after Gibiansky L, Gibiansky E. Target-mediated drug",
    "disposition model for drugs that bind to more than one target.",
    "J Pharmacokinet Pharmacodyn. 2010;37(4):323-346. doi:10.1007/s10928-010-9163-3.",
    sep = " "
  )
  vignette <- "Dua_2015_tmdd_tutorial"

  # Time is days. Every concentration is nM because the binding constants
  # (Kss_S, Kss_M) and both target levels are printed in nM (Supplement Model 6).
  # The drug states hold AMOUNTS per kg body weight (nmol/kg) so that a dose
  # enters directly; the supplement's tissue variable Lt is a concentration
  # referenced to the CENTRAL volume (Lt' = kpt*Lc - ktp*Lt, no volume ratio),
  # so the tissue amount is peripheral1 = Lt * vc and Lt = peripheral1 / vc.
  # The supplement doses 10 mg/kg and converts with
  # Dose = (Dosemg * 1e-3 / MWlig) * 1e9 nmol/kg, MWlig = 150000 Da, i.e.
  # 1 mg/kg = 6.667 nmol/kg; the conversion is applied in the vignette.
  units <- list(time = "day", dosing = "nmol/kg", concentration = "nM")

  compartmentData <- list(
    central = list(
      analyte = "drug (total: free + soluble-target complex + membrane-target complex)",
      units = "nmol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "drug (free)", units = "nmol/kg", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "soluble target (total: free + drug-bound)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    )
  )

  dosing <- "central"

  covariateData <- list()

  population <- list(
    species = "not applicable (theoretical model)",
    n_subjects = NA_integer_,
    disease_state = "none; a teaching implementation of a generic monoclonal antibody binding a soluble and a membrane-bound target",
    dose_range = "10 mg/kg single i.v. bolus into the central compartment (Supplement Model 6 'Dosing Regime')",
    notes = paste(
      "TUTORIAL model. Dua 2015 is a review/tutorial of TMDD model structures; Supplementary",
      "Model 6 is the authors' runnable Berkeley Madonna implementation of the two-target",
      "(soluble S + membrane M) QSS/MM approximation of Gibiansky & Gibiansky 2010, with an",
      "illustrative parameter set that is not fitted to any data. Every value is printed in the",
      "supplement and is encoded with fixed(). There is no inter-individual variability and no",
      "residual error in the source; propSd is fixed(0) so the model can be simulated as written.",
      "Weight-normalised: the central volume is 0.04 L/kg and doses are per kg."
    )
  )

  ini({
    # ---- Drug disposition (Supplement Model 6, 'Values for Rate Constants') ----
    lvc <- fixed(log(0.04)); label("Central volume of distribution (L/kg)") # Supplement Model 6: Vc = 0.04 L/kg
    lkel <- fixed(log(0.046)); label("Linear elimination rate constant of free drug from central (1/day)") # Supplement Model 6: kel = 0.046 1/day
    lk12 <- fixed(log(0.902)); label("Central-to-tissue transfer rate constant of free drug (1/day)") # Supplement Model 6: kpt = 0.902 1/day
    lk21 <- fixed(log(0.725)); label("Tissue-to-central transfer rate constant of free drug (1/day)") # Supplement Model 6: ktp = 0.725 1/day

    # ---- Soluble target (QSS) -------------------------------------------------
    lrbase_total_target <- fixed(log(0.00657)); label("Baseline soluble target concentration Rs0 (nM)") # Supplement Model 6: Rs0 = 0.00657 nM
    lkdeg_sol <- fixed(log(17.3)); label("First-order degradation rate constant of free soluble target (1/day)") # Supplement Model 6: kouts = 17.3 1/day
    lkint_sol <- fixed(log(0.170)); label("Elimination rate constant of the drug-soluble-target complex (1/day)") # Supplement Model 6: keps = 0.170 1/day
    lkon_sol <- fixed(log(30.2)); label("Association rate constant, drug to soluble target (1/(nM*day))") # Supplement Model 6: kons = 30.2 1/(nM day)
    lkoff_sol <- fixed(log(169)); label("Dissociation rate constant of the drug-soluble-target complex (1/day)") # Supplement Model 6: koffs = 169 1/day

    # ---- Membrane target (constant total; Michaelis-Menten term) -------------
    lrtot_mem <- fixed(log(0.002)); label("Total membrane-target concentration, held constant Rm0 (nM)") # Supplement Model 6: Rm0 = 0.002 nM
    lkint_mem <- fixed(log(0.170)); label("Elimination rate constant of the drug-membrane-target complex (1/day)") # Supplement Model 6: kepm = 0.170 1/day
    lkon_mem <- fixed(log(30.2)); label("Association rate constant, drug to membrane target (1/(nM*day))") # Supplement Model 6: konm = 30.2 1/(nM day)
    lkoff_mem <- fixed(log(169)); label("Dissociation rate constant of the drug-membrane-target complex (1/day)") # Supplement Model 6: koffm = 169 1/day

    # The source is a deterministic teaching implementation with no residual
    # error; fixed(0) keeps the model simulable without inventing a value.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - none in the source)")
  })

  model({
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    rbase <- exp(lrbase_total_target)
    kdeg_sol <- exp(lkdeg_sol)
    kint_sol <- exp(lkint_sol)
    kon_sol <- exp(lkon_sol)
    koff_sol <- exp(lkoff_sol)

    rtot_mem <- exp(lrtot_mem)
    kint_mem <- exp(lkint_mem)
    kon_mem <- exp(lkon_mem)
    koff_mem <- exp(lkoff_mem)

    # Supplement Model 6: Ksss = (koffs + keps)/kons, Kssm = (koffm + kepm)/konm,
    # Vmaxm = Rm0 * kepm, kins = kouts * Rs0 (main text Eq V.6: Rtot_S(0) = R0_S = kin_S/kout_S).
    kss_sol <- (koff_sol + kint_sol) / kon_sol
    kss_mem <- (koff_mem + kint_mem) / kon_mem
    vmax_mem <- rtot_mem * kint_mem
    ksyn_sol <- kdeg_sol * rbase

    total_target(0) <- rbase

    # Concentrations (nM). ctot is the supplement's Lctot; ct is its Lt.
    ctot <- central / vc
    ct <- peripheral1 / vc

    # Eq V.5 / Supplement: free drug from the soluble-target QSS quadratic.
    disc <- ctot - total_target - kss_sol
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kss_sol * ctot))

    # Supplement Model 6 'Finding Amount of Free Receptor': the soluble complex
    # is total drug less free drug less the membrane-bound drug Rm0*Lc/(Kssm+Lc).
    complex_sol <- ctot - cfree - rtot_mem * cfree / (kss_mem + cfree)
    target_sol <- total_target - complex_sol
    # The quadratic above solves Lctot = Lc + Rstot*Lc/(Ksss + Lc), which omits
    # the membrane-bound drug, so the printed complex_sol is a small difference
    # of two large numbers and is biased by up to Rm0 when drug is in excess.
    # The QSS relation itself gives the soluble complex without that bias.
    complex_sol_qss <- total_target * cfree / (kss_sol + cfree)
    target_sol_qss <- total_target - complex_sol_qss

    # Eq V.3 (i.v. form, Supplement Model 6), multiplied by vc for amounts.
    d/dt(central) <- vc * (-(kel + k12) * cfree - total_target * kint_sol * cfree / (kss_sol + cfree) -
      vmax_mem * cfree / (kss_mem + cfree)) + k21 * peripheral1
    # Eq V.2
    d/dt(peripheral1) <- k12 * cfree * vc - k21 * peripheral1
    # Eq V.4
    d/dt(total_target) <- ksyn_sol - kdeg_sol * total_target -
      (kint_sol - kdeg_sol) * total_target * cfree / (kss_sol + cfree)

    # Free drug concentration in the central compartment (supplement Lc).
    Cc <- cfree
    Cc ~ prop(propSd)
  })
}
