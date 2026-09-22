Schropp_2019_bsab_tmdd_qeconst <- function() {
  description <- "QSP. Quasi-equilibrium (QE) approximation of the bispecific-antibody (BsAb) TMDD model with CONSTANT total receptors - the most reduced member of the family, and the only one the authors actually fitted to data. Adding the assumption that each receptor degrades at the same rate as its complexes internalise (k_degX = k_intX = k_intAB = k_int) holds each total receptor pool at its baseline, which collapses the model to a SINGLE differential equation for free BsAb: both free receptors and all three complexes then follow algebraically, the free receptor B from the positive root of a quadratic. Parameter count falls from 19 to 10. As in the non-constant sibling, rapid binding puts the i.v. input function inside the equation, so an i.v. dose must be delivered through a dummy input compartment rather than added to the central state. Defaults are the true (data-generating) parameter values of the paper's simulation-estimation study, which is the one place the paper attaches inter-individual variability and residual error to this model. Siblings: modellib('Schropp_2019_bsab_tmdd_full') and modellib('Schropp_2019_bsab_tmdd_qe')."
  reference <- paste(
    "Schropp J, Khot A, Shah DK, Koch G. Target-Mediated Drug Disposition Model",
    "for Bispecific Antibodies: Properties, Approximation, and Optimal Dosing",
    "Strategy. CPT Pharmacometrics Syst Pharmacol. 2019;8(3):177-187.",
    "doi:10.1002/psp4.12369.",
    "Structural equations are Eqs 15, 17 and 23-30, with m_11 and the",
    "determinant from Table 1; parameter values are Table 2 (right-hand",
    "'Simulation-estimation' block). The authors' own NONMEM encoding of this",
    "model is Supplementary Material S10, and the MONOLIX encoding is S8.",
    sep = " "
  )
  vignette <- "Schropp_2019_bsab_tmdd"

  # Units and the amount/concentration split follow the rest of the family; see
  # Schropp_2019_bsab_tmdd_full.R. Doses are nmol (1 mg = 6.7 nmol).
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # No peripheral compartment. This is the authors' own configuration for this
  # model in BOTH places they use it: their NONMEM stream S10 declares NCOMP=1,
  # and the Figure 5a-c optimal-dosing example states 'for simplicity, the
  # peripheral compartment was neglected'. The paper's s.c.-plus-peripheral
  # example (Figure 5d-f) is simulated with the FULL model instead, which is
  # available as Schropp_2019_bsab_tmdd_full.
  compartmentData <- list(
    depot_ivdum = list(
      analyte = "bispecific antibody (i.v. input-function device)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    depot_sc = list(
      analyte = "bispecific antibody (s.c. depot)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "plasma", verified = TRUE)
  )

  # An i.v. dose MUST be given into `depot_ivdum`, never into `central`.
  dosing <- c("depot_ivdum", "depot_sc")

  covariateData <- list()

  population <- list(
    species = "not applicable (simulated data)",
    n_subjects = 30L,
    disease_state = "none; simulated subjects in an identifiability (simulation-estimation) study",
    dose_range = "i.v. bolus; the authors' own S10 example stream doses 1675 nmol = 250 mg",
    notes = paste(
      "The paper's simulation-estimation study simulated 30 individuals with the FULL model and",
      "refitted them with this constant-total-receptor QE approximation. Study 1 used 12 free-BsAb",
      "measurements per subject over 40 days (until the 0.01 nM LLOQ was reached) and could not",
      "identify K_D1 / K_D2, which had to be fixed. Study 2 added 18 measurements each of free",
      "receptor A and free receptor B over 70 days and recovered K_D1 and K_D2 as well, with alpha",
      "still fixed. The ini() defaults here are the TRUE (data-generating) values of Table 2, not",
      "either fitted column, because the true values are the authoritative generating set and the",
      "two fitted columns (MONOLIX and NONMEM, Study 1 and Study 2) disagree with each other in the",
      "third significant figure. All four fitted columns are reproduced in the vignette. Relative",
      "standard errors were all below 5% except for the IIV standard deviations (below 100%)."
    )
  )

  ini({
    # ---- Disposition (Table 2, 'True values' column) --------------------------
    lvc <- fixed(log(3)); label("Central volume of distribution (L)")                    # Table 2 simulation-estimation, true V = 3 L
    lkel <- fixed(log(0.1)); label("Linear elimination rate constant from central (1/day)")  # Table 2 simulation-estimation, true k_el = 0.1 1/day
    lka <- fixed(log(0.2)); label("First-order s.c. absorption rate constant (1/day)")   # Table 2 simulation block, k_a = 0.2 1/day (the simulation-estimation study dosed i.v. only)
    lfdepot <- fixed(log(0.75)); label("Bioavailability of the s.c. dose (fraction)")    # Table 2 simulation block, F = 0.75 (the simulation-estimation study dosed i.v. only)

    # ---- Rapid-binding parameters ---------------------------------------------
    lkd1 <- fixed(log(0.1)); label("Equilibrium dissociation constant of the BsAb for receptor A (nM)")  # Table 2 simulation-estimation, true K_D1 = 0.1 nM
    lkd2 <- fixed(log(1)); label("Equilibrium dissociation constant of the BsAb for receptor B (nM)")    # Table 2 simulation-estimation, true K_D2 = 1 nM
    lalpha <- fixed(log(1)); label("Cross-linking cooperativity factor alpha (unitless)")  # Table 2 simulation-estimation, alpha = 1 and FIXED by the authors in every fit

    # ---- Constant total receptors (Eqs 23-24) ---------------------------------
    lrtota <- fixed(log(10)); label("Constant total concentration of receptor A (nM)")    # Table 2 simulation-estimation, true R_A^0 = 10 nM
    lrtotb <- fixed(log(100)); label("Constant total concentration of receptor B (nM)")   # Table 2 simulation-estimation, true R_B^0 = 100 nM

    # The single internalisation/degradation rate that the constant-total-receptor
    # assumption forces: k_degA = k_degB = k_intA = k_intB = k_intAB = k_int
    # (Eq 23). The authors' NONMEM stream S10 makes exactly this assignment.
    lkint <- fixed(log(0.1)); label("Shared receptor degradation and complex internalisation rate constant (1/day)")  # Table 2 simulation-estimation, true k_int = 0.1 1/day

    # ---- i.v. input-function device ------------------------------------------
    # See Schropp_2019_bsab_tmdd_qe.R for the full reasoning. 1 / TDUR with the
    # authors' own TDUR = 0.0001 day (S10).
    kdum <- fixed(1e4); label("Rate constant of the dummy i.v. input compartment (1/day); numerical device, not a PK parameter")  # 1 / TDUR with TDUR = 0.0001 day (S10)

    # ---- Inter-individual variability (Table 2, 'True values' column) ----------
    # The paper states these are 'standard deviations of the lognormally
    # distributed interindividual variability', so the variances recorded here
    # are omega^2 = 0.05^2 = 0.0025.
    etalkel ~ 0.0025 # Table 2 simulation-estimation, true 'omega kel' = 0.05, reported as an SD
    etalvc ~ 0.0025 # Table 2 simulation-estimation, true 'omega V' = 0.05, reported as an SD

    # ---- Residual error --------------------------------------------------------
    # The paper fits a proportional residual error to each of the three
    # observed profiles: b1 on free BsAb, b2 on free receptor A, b3 on free
    # receptor B (Table 2 footnote). All three true values are 0.2.
    propSd <- fixed(0.2); label("Proportional residual error on free BsAb concentration (fraction)")   # Table 2 simulation-estimation, true b_1 = 0.2
    propSd_RA <- fixed(0.2); label("Proportional residual error on free receptor A (fraction)")           # Table 2 simulation-estimation, true b_2 = 0.2
    propSd_RB <- fixed(0.2); label("Proportional residual error on free receptor B (fraction)")           # Table 2 simulation-estimation, true b_3 = 0.2
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    kd1 <- exp(lkd1)
    kd2 <- exp(lkd2)
    alpha <- exp(lalpha)
    rtota0 <- exp(lrtota)
    rtotb0 <- exp(lrtotb)
    kint <- exp(lkint)

    # Eq 42: dummy i.v. input compartment supplying the paper's In_IV(t).
    d/dt(depot_ivdum) <- -kdum * depot_ivdum
    in_iv <- kdum * depot_ivdum

    # Eq 15: s.c. depot.
    d/dt(depot_sc) <- -ka * depot_sc

    # Eqs 27-30: free receptor B is the positive root of a*R_B^2 + b*R_B + d = 0.
    # The branch reproduces the authors' own `IF (C.GT.0)` guard in S10: at
    # C = 0 the quadratic degenerates (its leading coefficient vanishes) and the
    # free receptor is simply its own constant total, R_totB^0.
    qa <- (1 + central / kd2) * (central / (alpha * kd1 * kd2))
    qb <- (central * (rtota0 - rtotb0)) / (alpha * kd1 * kd2) +
      (1 + central / kd1) * (1 + central / kd2)
    qd <- -rtotb0 * (1 + central / kd1)
    if (central > 0) {
      rb <- (-qb + sqrt(qb * qb - 4 * qa * qd)) / (2 * qa)
    } else {
      rb <- rtotb0
    }

    # Eq 26: free receptor A, given R_B.
    ra <- rtota0 / (1 + central / kd1 + rb * central / (alpha * kd1 * kd2))

    # Table 1: only m_11 and the determinant are needed, because the constant
    # total receptors leave a single differential equation. Transcribed from the
    # authors' NONMEM stream S10.
    det1 <- central * kd2 * ra^2 + central^2 * kd2 * ra +
      central * kd1 * rb^2 + central^2 * kd1 * rb + central^2 * ra * rb
    det2 <- alpha * kd1^2 * kd2^2 + central * kd1 * kd2 * ra +
      central * kd1 * kd2 * rb + kd1 * kd2 * ra * rb
    det3 <- alpha * central * kd1 * kd2^2 + alpha * central * kd1^2 * kd2 +
      alpha * central^2 * kd1 * kd2
    det4 <- alpha * kd1 * kd2^2 * ra + alpha * kd1^2 * kd2 * rb +
      alpha * central * kd1 * kd2 * ra + alpha * central * kd1 * kd2 * rb
    det <- (det1 + det2 + det3 + det4) / (alpha * kd1^2 * kd2^2)

    m111 <- central^2 * kd2 * ra + central^2 * kd1 * rb +
      alpha * kd1^2 * kd2^2 + central * kd1 * kd2 * ra
    m112 <- central * kd1 * kd2 * rb + alpha * central * kd1 * kd2^2 +
      alpha * central * kd1^2 * kd2 + alpha * central^2 * kd1 * kd2
    m11 <- (m111 + m112) / (alpha * kd1^2 * kd2^2)

    # Eq 25: the single free-BsAb differential equation. Every internalisation
    # term shares the one rate k_int, which is what makes the total receptors
    # constant.
    g1 <- in_iv / vc + ka * depot_sc / vc - kel * central -
      kint * (ra * central) / kd1 -
      kint * (rb * central) / kd2 -
      kint * (ra * rb * central) / (alpha * kd1 * kd2)

    d/dt(central) <- (m11 / det) * g1

    f(depot_sc) <- fdepot

    # Eqs 20-22: the complexes, algebraic under rapid binding.
    complexA <- central * ra / kd1
    complexB <- central * rb / kd2
    trimerAB <- central * ra * rb / (alpha * kd1 * kd2)

    # Eq 11. R_totA and R_totB are constant by construction (Eq 23), so they are
    # not recomputed here - they are the parameters rtota0 and rtotb0.
    Cc <- central
    CcTotal <- central + complexA + complexB + trimerAB
    RA <- ra
    RB <- rb

    # The three observed profiles of the simulation-estimation study.
    Cc ~ prop(propSd)
    RA ~ prop(propSd_RA)
    RB ~ prop(propSd_RB)
  })
}
