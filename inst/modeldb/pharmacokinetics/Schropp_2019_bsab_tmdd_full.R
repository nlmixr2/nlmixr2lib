Schropp_2019_bsab_tmdd_full <- function() {
  description <- "QSP. Full target-mediated drug disposition (TMDD) model for a bispecific antibody (BsAb) that binds two different targets carried on the membranes of two different cells. Free BsAb binds free receptor A or free receptor B to form the two binary complexes, and each binary complex cross-links with the opposite receptor to form the pharmacologically active ternary complex (trimer) that bridges the two cells. All four binding events are written as reversible mass action, both receptors turn over with zero-order synthesis and first-order degradation, and all three complexes internalise. Drug disposition is linear elimination from the central compartment plus a peripheral compartment, with i.v. and s.c. routes. The signature behaviour is a BELL-SHAPED trimer response: raising the BsAb dose eventually DELAYS and then suppresses trimer build-up, because saturating both receptors with binary complexes starves the cross-linking reaction. Defaults are the paper's own generic simulation parameter set (Table 2, left-hand block); this is a theoretical model with no fitted dataset, so no drug and no real population is attached. Siblings in this library: modellib('Schropp_2019_bsab_tmdd_qe') and modellib('Schropp_2019_bsab_tmdd_qeconst')."
  reference <- paste(
    "Schropp J, Khot A, Shah DK, Koch G. Target-Mediated Drug Disposition Model",
    "for Bispecific Antibodies: Properties, Approximation, and Optimal Dosing",
    "Strategy. CPT Pharmacometrics Syst Pharmacol. 2019;8(3):177-187.",
    "doi:10.1002/psp4.12369.",
    "Structural equations are Eqs 1-9; parameter values are Table 2 (left-hand",
    "'Simulation' block). The authors' own NONMEM encoding of this model is",
    "Supplementary Material S6, and the MONOLIX encoding is S4.",
    sep = " "
  )
  vignette <- "Schropp_2019_bsab_tmdd"

  # Time is days and every concentration is nM, because the four binding terms
  # of Eqs 2-7 multiply drug by receptor concentration and the published
  # affinities are quoted in nM (Table 2).
  #
  # The paper deliberately MIXES concentration states and amount states, and
  # this file reproduces that exactly rather than renormalising, because the
  # published equations and the authors' own NONMEM code (S6) are written that
  # way. `central` holds the free BsAb CONCENTRATION C (nM); `peripheral1` and
  # `depot_sc` hold AMOUNTS (nmol). That is why Eq 2 carries `ka * AD / V` and
  # `k21 * AP / V` while Eq 8 carries `k12 * V * C` -- the volume appears
  # explicitly at each amount/concentration boundary. It is also why an i.v.
  # dose enters through `f(central) <- 1 / vc` (NONMEM `F1 = 1/VC` in S6).
  #
  # Doses are nmol. The paper works in mg and states its own conversion factor
  # of 1 mg = 6.7 nmol (Modeling and simulation, "Parameter setting for
  # simulations"); that factor is confirmed internally by the authors' NONMEM
  # dataset comment `dose = 335` for the 50 mg i.v. bolus of Figure 3 (S9/S10)
  # and by the worked optimal dose of 44.8 mg for 300 nmol (Eq 39, Results).
  # It is applied in the vignette, not here, so the model file stays molar.
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Both targets are deliberately GENERIC -- the paper names them only A and B
  # and never commits to an antigen pair, so the states carry the paper's own
  # A/B identity rather than an invented antigen name. `trimer` is the existing
  # canonical for the productive ternary complex of a cell-bridging bispecific.
  # `specimen` is "not applicable" for the receptor and complex states: they are
  # membrane-borne species on two unnamed cell types, not anything drawn into a
  # sampled matrix.
  compartmentData <- list(
    depot_sc = list(
      analyte = "bispecific antibody (s.c. depot)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bispecific antibody (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    target_a = list(analyte = "target receptor A (free)", units = "nM", specimen = "not applicable", verified = TRUE),
    target_b = list(analyte = "target receptor B (free)", units = "nM", specimen = "not applicable", verified = TRUE),
    complex_a = list(
      analyte = "BsAb-receptor A binary complex",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_b = list(
      analyte = "BsAb-receptor B binary complex",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    trimer = list(
      analyte = "receptor A-BsAb-receptor B ternary complex",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  dosing <- c("central", "depot_sc")

  covariateData <- list()

  population <- list(
    species = "not applicable (theoretical model)",
    n_subjects = NA_integer_,
    disease_state = "none; the paper motivates the model with immuno-oncology BsAbs that bridge an effector T cell to a tumor cell, but fits no clinical or preclinical dataset",
    dose_range = "50 mg and 250 mg i.v. bolus (Figure 3), 500 mg i.v. bolus (Figure 4), up to 700 mg s.c. (Modeling and simulation)",
    notes = paste(
      "This is a THEORETICAL / methodological paper. There is no subject-level dataset, no",
      "demographic table and no real drug: the model is a general structure for the BsAb class,",
      "and the defaults in ini() are the authors' own generic simulation parameter set (Table 2,",
      "left-hand 'Simulation' block), which they describe as 'roughly based on literature reported",
      "values' citing refs 19 and 26-28. Every value is therefore encoded with fixed(), matching",
      "the authors' own NONMEM control stream (S6), in which all 21 THETAs carry FIX and every",
      "OMEGA is 0 FIX. The only place the paper attaches inter-individual variability and residual",
      "error to anything is its simulation-estimation study, which used a DIFFERENT parameter set",
      "and fitted the quasi-equilibrium approximation; those values live on the sibling models",
      "Schropp_2019_bsab_tmdd_qe and Schropp_2019_bsab_tmdd_qeconst, not here.",
      "A lower limit of quantification of 0.01 nM was assumed for the paper's simulations."
    )
  )

  ini({
    # ---- Disposition (Table 2, left-hand 'Simulation' block) ------------------
    lvc <- fixed(log(3)); label("Central volume of distribution (L)")                    # Table 2 simulation block, V = 3 L
    lka <- fixed(log(0.2)); label("First-order s.c. absorption rate constant (1/day)")   # Table 2 simulation block, k_a = 0.2 1/day
    lkel <- fixed(log(0.1)); label("Linear elimination rate constant from central (1/day)")  # Table 2 simulation block, k_el = 0.1 1/day
    lk12 <- fixed(log(0.1)); label("Central-to-peripheral distribution rate constant (1/day)")  # Table 2 simulation block, k_12 = 0.1 1/day
    lk21 <- fixed(log(0.03)); label("Peripheral-to-central distribution rate constant (1/day)")  # Table 2 simulation block, k_21 = 0.03 1/day
    lfdepot <- fixed(log(0.75)); label("Bioavailability of the s.c. dose (fraction)")    # Table 2 simulation block, F = 0.75

    # ---- Binding (Eq 1 numbering: Z = 1 drug+A, 2 drug+B, 3 RC_A+B, 4 RC_B+A) -
    # The paper states 'For simplicity, we assume k_on3 = k_on2, k_on4 = k_on1,
    # k_off3 = k_off2, and k_off4 = k_off1' (Table 2 footnote), which is the
    # alpha = 1 case of the microscopic-reversibility constraint K_D1*K_D3 =
    # K_D2*K_D4 (Eq 14). The cross-link rates are therefore DERIVED in model()
    # from the first-binding rates rather than carried as separate parameters,
    # exactly as the authors' NONMEM stream S6 fixes THETA(6)=THETA(4),
    # THETA(7)=THETA(5), THETA(8)=THETA(2), THETA(9)=THETA(3).
    lkon1 <- fixed(log(10)); label("Association rate constant, free BsAb to free receptor A (1/(nM*day))")   # Table 2 simulation block, k_on1 = 10 1/(nM day)
    lkoff1 <- fixed(log(0.01)); label("Dissociation rate constant of the BsAb-receptor A binary complex (1/day)")  # Table 2 simulation block, k_off1 = 0.01 1/day
    lkon2 <- fixed(log(1)); label("Association rate constant, free BsAb to free receptor B (1/(nM*day))")    # Table 2 simulation block, k_on2 = 1 1/(nM day)
    lkoff2 <- fixed(log(0.01)); label("Dissociation rate constant of the BsAb-receptor B binary complex (1/day)")  # Table 2 simulation block, k_off2 = 0.01 1/day

    # ---- Target turnover ------------------------------------------------------
    # R_A(0) = ksynA/kdegA = 10 nM and R_B(0) = ksynB/kdegB = 100 nM (Eq 9),
    # which reproduces the R_A^0 = 10 nM and R_B^0 = 100 nM rows of Table 2.
    lksyna <- fixed(log(1)); label("Zero-order synthesis rate of receptor A (nM/day)")   # Table 2 simulation block, k_synA = 1 nM/day
    lkdega <- fixed(log(0.1)); label("First-order degradation rate constant of free receptor A (1/day)")  # Table 2 simulation block, k_degA = 0.1 1/day
    lksynb <- fixed(log(10)); label("Zero-order synthesis rate of receptor B (nM/day)")  # Table 2 simulation block, k_synB = 10 nM/day
    lkdegb <- fixed(log(0.1)); label("First-order degradation rate constant of free receptor B (1/day)")  # Table 2 simulation block, k_degB = 0.1 1/day

    # ---- Complex internalisation ---------------------------------------------
    lkinta <- fixed(log(0.05)); label("Internalisation rate constant of the BsAb-receptor A binary complex (1/day)")  # Table 2 simulation block, k_intA = 0.05 1/day
    lkintb <- fixed(log(0.05)); label("Internalisation rate constant of the BsAb-receptor B binary complex (1/day)")  # Table 2 simulation block, k_intB = 0.05 1/day
    lkintab <- fixed(log(0.1)); label("Internalisation rate constant of the ternary complex (1/day)")  # Table 2 simulation block, k_intAB = 0.1 1/day

    # The paper reports no residual error for the full model used in this way:
    # its NONMEM stream S6 is a pure simulation control stream with $ESTIMATION
    # commented out, an additive error fixed to 0 and a placeholder proportional
    # error fixed to 1. Encoded as fixed(0) rather than inventing a value; the
    # only residual errors the paper estimates (b1 = b2 = b3 = 0.2) belong to
    # the simulation-estimation study on the QE approximation and are carried by
    # the sibling model files. See the vignette Errata.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported for this model in the source)")
  })

  model({
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot <- exp(lfdepot)

    kon1 <- exp(lkon1)
    koff1 <- exp(lkoff1)
    kon2 <- exp(lkon2)
    koff2 <- exp(lkoff2)
    # Cross-linking rates, per the Table 2 footnote (alpha = 1).
    kon3 <- kon2
    koff3 <- koff2
    kon4 <- kon1
    koff4 <- koff1

    ksyna <- exp(lksyna)
    kdega <- exp(lkdega)
    ksynb <- exp(lksynb)
    kdegb <- exp(lkdegb)
    kinta <- exp(lkinta)
    kintb <- exp(lkintb)
    kintab <- exp(lkintab)

    # Eq 9: the BsAb is an engineered molecule with no endogenous baseline, so
    # every drug and complex state starts at zero, and each free receptor starts
    # at its own drug-free turnover steady state.
    target_a(0) <- ksyna / kdega
    target_b(0) <- ksynb / kdegb

    # Eq 1: s.c. depot.
    d/dt(depot_sc) <- -ka * depot_sc

    # Eq 2: free BsAb concentration. Linear elimination, both first-binding
    # events, distribution to the periphery, and s.c. input.
    d/dt(central) <-
      -kel * central -
      kon1 * central * target_a + koff1 * complex_a -
      kon2 * central * target_b + koff2 * complex_b -
      k12 * central + k21 * peripheral1 / vc +
      ka * depot_sc / vc

    # Eq 3: free receptor A. Turnover, binding to free BsAb, and consumption by
    # the cross-linking reaction that converts the B binary complex to trimer.
    d/dt(target_a) <-
      ksyna - kdega * target_a -
      kon1 * central * target_a + koff1 * complex_a -
      kon4 * target_a * complex_b + koff4 * trimer

    # Eq 4: free receptor B, the structural mirror of Eq 3.
    d/dt(target_b) <-
      ksynb - kdegb * target_b -
      kon2 * central * target_b + koff2 * complex_b -
      kon3 * target_b * complex_a + koff3 * trimer

    # Eq 5: BsAb-receptor A binary complex.
    d/dt(complex_a) <-
      kon1 * central * target_a - (koff1 + kinta) * complex_a -
      kon3 * target_b * complex_a + koff3 * trimer

    # Eq 6: BsAb-receptor B binary complex.
    d/dt(complex_b) <-
      kon2 * central * target_b - (koff2 + kintb) * complex_b -
      kon4 * target_a * complex_b + koff4 * trimer

    # Eq 7: the ternary complex, formed from EITHER binary complex by
    # cross-linking with the opposite free receptor. This is the
    # pharmacologically active species.
    d/dt(trimer) <-
      kon4 * target_a * complex_b + kon3 * target_b * complex_a -
      (koff3 + koff4 + kintab) * trimer

    # Eq 8: peripheral compartment, carrying an AMOUNT.
    d/dt(peripheral1) <- k12 * central * vc - k21 * peripheral1

    # An i.v. dose is supplied in nmol and enters as a central CONCENTRATION
    # (NONMEM `F1 = 1/VC`, S6); an s.c. dose enters the depot in nmol.
    f(central) <- 1 / vc
    f(depot_sc) <- fdepot

    # Eqs 11-13: the total-drug and total-receptor variables the paper's optimal
    # working range Eq 36 is written on.
    Cc <- central
    CcTotal <- central + complex_a + complex_b + trimer
    RtotA <- target_a + complex_a + trimer
    RtotB <- target_b + complex_b + trimer

    Cc ~ prop(propSd)
  })
}
