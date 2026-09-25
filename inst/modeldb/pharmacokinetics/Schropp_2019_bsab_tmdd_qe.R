Schropp_2019_bsab_tmdd_qe <- function() {
  description <- "QSP. Quasi-equilibrium (QE) approximation of the full bispecific-antibody (BsAb) TMDD model, with NON-constant total receptors. Assuming all four binding events are infinitely fast collapses the eight-state full model onto three differential equations (free BsAb, free receptor A, free receptor B) plus the peripheral and s.c.-depot states, and replaces the eight association/dissociation rate constants by just three binding parameters: K_D1, K_D2 and the cross-linking cooperativity alpha (microscopic reversibility forces K_D3 = alpha*K_D2 and K_D4 = alpha*K_D1). The two binary complexes and the ternary complex are then ALGEBRAIC rather than integrated. The price of the approximation is that an i.v. dose can no longer be added straight to the central state: rapid binding instantly partitions the incoming drug between free and bound species, so the input function is multiplied into the 3x3 matrix and must be delivered through a dummy input compartment (the paper's Eq 42). Parameter count falls from 19 to 14. Defaults are the paper's generic simulation set (Table 2, left-hand block) re-expressed in K_D form. Siblings: modellib('Schropp_2019_bsab_tmdd_full') and modellib('Schropp_2019_bsab_tmdd_qeconst')."
  reference <- paste(
    "Schropp J, Khot A, Shah DK, Koch G. Target-Mediated Drug Disposition Model",
    "for Bispecific Antibodies: Properties, Approximation, and Optimal Dosing",
    "Strategy. CPT Pharmacometrics Syst Pharmacol. 2019;8(3):177-187.",
    "doi:10.1002/psp4.12369.",
    "Structural equations are Eqs 14-22 and 24, with the matrix entries m_ij and",
    "the determinant taken from Table 1. The authors' own NONMEM encoding of this",
    "model is Supplementary Material S5, and the MONOLIX encoding is S3.",
    sep = " "
  )
  vignette <- "Schropp_2019_bsab_tmdd"

  # Units and the amount/concentration split are identical to the full model;
  # see Schropp_2019_bsab_tmdd_full.R for the reasoning. `central` holds the
  # free BsAb CONCENTRATION (nM), `peripheral1` / `depot_sc` / `depot_iv`
  # hold AMOUNTS (nmol), and doses are nmol (1 mg = 6.7 nmol per the paper's
  # own conversion factor).
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  compartmentData <- list(
    depot_iv = list(
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
    central = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bispecific antibody (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    target_a = list(analyte = "target receptor A (free)", units = "nM", specimen = "not applicable", verified = TRUE),
    target_b = list(analyte = "target receptor B (free)", units = "nM", specimen = "not applicable", verified = TRUE)
  )

  # An i.v. dose MUST be given into `depot_iv`, never into `central`. See the
  # `kdum` comment in ini() and the vignette for why.
  dosing <- c("depot_iv", "depot_sc")

  covariateData <- list()

  population <- list(
    species = "not applicable (theoretical model)",
    n_subjects = NA_integer_,
    disease_state = "none; a general structural approximation for the BsAb class, fitted to no clinical or preclinical dataset",
    dose_range = "50 mg i.v. bolus (the authors' own S5 example stream doses 335 nmol = 50 mg); up to 700 mg s.c.",
    notes = paste(
      "Theoretical / methodological paper with no subject-level dataset. Defaults are the authors'",
      "generic simulation parameter set (Table 2, left-hand 'Simulation' block), converted to the",
      "QE parameterisation by K_D1 = k_off1 / k_on1 = 0.01 / 10 = 0.001 nM and K_D2 = k_off2 /",
      "k_on2 = 0.01 / 1 = 0.01 nM (Eq 14a). Those are exactly the values the authors hard-code as",
      "THETA(2) = 1E-3 and THETA(3) = 0.01 in their own NONMEM stream S5, which confirms the",
      "conversion. alpha = 1 follows from the paper's stated simplification k_on3 = k_on2,",
      "k_on4 = k_on1, k_off3 = k_off2, k_off4 = k_off1 (Table 2 footnote), which the text states is",
      "equivalent to K_D1 = K_D4 and K_D2 = K_D3, hence alpha = 1. Every value is fixed(), matching",
      "the authors' all-FIX control stream. The paper attaches inter-individual variability and",
      "residual error only in its simulation-estimation study, which fitted the CONSTANT-total-receptor",
      "reduction; those values are carried by Schropp_2019_bsab_tmdd_qeconst."
    )
  )

  ini({
    lvc <- fixed(log(3)); label("Central volume of distribution (L)")                    # Table 2 simulation block, V = 3 L
    lka <- fixed(log(0.2)); label("First-order s.c. absorption rate constant (1/day)")   # Table 2 simulation block, k_a = 0.2 1/day
    lkel <- fixed(log(0.1)); label("Linear elimination rate constant from central (1/day)")  # Table 2 simulation block, k_el = 0.1 1/day
    lk12 <- fixed(log(0.1)); label("Central-to-peripheral distribution rate constant (1/day)")  # Table 2 simulation block, k_12 = 0.1 1/day
    lk21 <- fixed(log(0.03)); label("Peripheral-to-central distribution rate constant (1/day)")  # Table 2 simulation block, k_21 = 0.03 1/day
    lfdepot <- fixed(log(0.75)); label("Bioavailability of the s.c. dose (fraction)")    # Table 2 simulation block, F = 0.75

    # ---- Rapid-binding parameters (Eq 14a, Eq 14b) ----------------------------
    # K_DZ = k_offZ / k_onZ. From the Table 2 simulation block: K_D1 = 0.01 / 10
    # and K_D2 = 0.01 / 1. Confirmed by the authors' NONMEM stream S5, which
    # fixes THETA(2) = 1E-3 and THETA(3) = 0.01.
    lkd1 <- fixed(log(0.001)); label("Equilibrium dissociation constant of the BsAb for receptor A (nM)")   # derived from Table 2, k_off1 / k_on1 = 0.01 / 10; S5 THETA(2) = 1E-3
    lkd2 <- fixed(log(0.01)); label("Equilibrium dissociation constant of the BsAb for receptor B (nM)")    # derived from Table 2, k_off2 / k_on2 = 0.01 / 1; S5 THETA(3) = 0.01
    # alpha is the cross-linking cooperativity: the affinity of C for RC_B
    # relative to R_A (and of C for RC_A relative to R_B). Microscopic
    # reversibility (Eq 14b) gives K_D3 = alpha*K_D2 and K_D4 = alpha*K_D1.
    lalpha <- fixed(log(1)); label("Cross-linking cooperativity factor alpha (unitless)")  # Table 2, alpha = 1 (fixed by the authors); equivalent to the Table 2 footnote simplification

    # ---- Target turnover ------------------------------------------------------
    lksyna <- fixed(log(1)); label("Zero-order synthesis rate of receptor A (nM/day)")   # Table 2 simulation block, k_synA = 1 nM/day
    lkdega <- fixed(log(0.1)); label("First-order degradation rate constant of free receptor A (1/day)")  # Table 2 simulation block, k_degA = 0.1 1/day
    lksynb <- fixed(log(10)); label("Zero-order synthesis rate of receptor B (nM/day)")  # Table 2 simulation block, k_synB = 10 nM/day
    lkdegb <- fixed(log(0.1)); label("First-order degradation rate constant of free receptor B (1/day)")  # Table 2 simulation block, k_degB = 0.1 1/day

    # ---- Complex internalisation ---------------------------------------------
    lkinta <- fixed(log(0.05)); label("Internalisation rate constant of the BsAb-receptor A binary complex (1/day)")  # Table 2 simulation block, k_intA = 0.05 1/day
    lkintb <- fixed(log(0.05)); label("Internalisation rate constant of the BsAb-receptor B binary complex (1/day)")  # Table 2 simulation block, k_intB = 0.05 1/day
    lkintab <- fixed(log(0.1)); label("Internalisation rate constant of the ternary complex (1/day)")  # Table 2 simulation block, k_intAB = 0.1 1/day

    # ---- i.v. input-function device ------------------------------------------
    # NOT a pharmacological parameter. Because rapid binding makes the i.v.
    # input appear INSIDE the matrix product of Eq 16, an i.v. dose cannot be
    # added straight to `central`; the paper devotes a whole subsection to this
    # ('Implementation of the QE approximation with an i.v. administration') and
    # instructs that the bolus be mimicked by a very short infusion. The authors
    # hard-code TDUR = 0.0001 day in their own NONMEM streams (S5, S9, S10). The
    # dummy compartment of the paper's Eq 42 is reproduced here by a first-order
    # drain whose mean residence time is that same 0.0001 day, i.e.
    # kdum = 1 / 0.0001 = 1e4 1/day. It is mass-conserving (the time integral of
    # the input equals the dose exactly) and, being smooth, is kinder to the
    # solver than a discontinuous rate.
    kdum <- fixed(1e4); label("Rate constant of the dummy i.v. input compartment (1/day); numerical device, not a PK parameter")  # 1 / TDUR with TDUR = 0.0001 day, the authors' own short-infusion duration (S5, S9, S10)

    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported for this model in the source)")
  })

  model({
    vc <- exp(lvc)
    ka <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot <- exp(lfdepot)

    kd1 <- exp(lkd1)
    kd2 <- exp(lkd2)
    alpha <- exp(lalpha)

    ksyna <- exp(lksyna)
    kdega <- exp(lkdega)
    ksynb <- exp(lksynb)
    kdegb <- exp(lkdegb)
    kinta <- exp(lkinta)
    kintb <- exp(lkintb)
    kintab <- exp(lkintab)

    # Initial free receptors, unchanged from the full model (Eq 9).
    target_a(0) <- ksyna / kdega
    target_b(0) <- ksynb / kdegb

    # Eq 42: the dummy i.v. input compartment. `in_iv` is the paper's In_IV(t).
    d/dt(depot_iv) <- -kdum * depot_iv
    in_iv <- kdum * depot_iv

    # Eq 15: s.c. depot, which needs no special treatment - the paper notes that
    # 'in case of an absorption compartment, the internal dosing mechanisms from
    # the PK/PD software can be used as usual due to the structure of Eq. 15'.
    d/dt(depot_sc) <- -ka * depot_sc

    # Table 1: entries of the matrix M_BsAb(C, R_A, R_B) and its determinant.
    # Transcribed from the authors' NONMEM stream S5, which is Table 1 written
    # out in FORTRAN and therefore the least ambiguous form of it.
    det1 <- central * kd2 * target_a^2 + central^2 * kd2 * target_a +
      central * kd1 * target_b^2 + central^2 * kd1 * target_b +
      central^2 * target_a * target_b
    det2 <- alpha * kd1^2 * kd2^2 + central * kd1 * kd2 * target_a +
      central * kd1 * kd2 * target_b + kd1 * kd2 * target_a * target_b
    det3 <- alpha * central * kd1 * kd2^2 + alpha * central * kd1^2 * kd2 +
      alpha * central^2 * kd1 * kd2
    det4 <- alpha * kd1 * kd2^2 * target_a + alpha * kd1^2 * kd2 * target_b +
      alpha * central * kd1 * kd2 * target_a + alpha * central * kd1 * kd2 * target_b
    det <- (det1 + det2 + det3 + det4) / (alpha * kd1^2 * kd2^2)

    m111 <- central^2 * kd2 * target_a + central^2 * kd1 * target_b +
      alpha * kd1^2 * kd2^2 + central * kd1 * kd2 * target_a
    m112 <- central * kd1 * kd2 * target_b + alpha * central * kd1 * kd2^2 +
      alpha * central * kd1^2 * kd2 + alpha * central^2 * kd1 * kd2
    m11 <- (m111 + m112) / (alpha * kd1^2 * kd2^2)
    m12 <- -(central^2 * target_a + central * kd1 * target_b +
      alpha * central^2 * kd1 + alpha * central * kd1 * kd2) / (alpha * kd1^2 * kd2)
    m13 <- -(central^2 * target_b + central * kd2 * target_a +
      alpha * central^2 * kd2 + alpha * central * kd1 * kd2) / (alpha * kd1 * kd2^2)

    m21 <- -(central * target_a^2 + kd1 * target_a * target_b +
      alpha * central * kd1 * target_a + alpha * kd1 * kd2 * target_a) / (alpha * kd1^2 * kd2)
    m221 <- central * target_a^2 + central * kd1 * target_a +
      kd1 * target_a * target_b + alpha * central * kd1^2 + alpha * kd1^2 * kd2
    m222 <- alpha * kd1^2 * target_b + alpha * central * kd1 * target_a +
      alpha * kd1 * kd2 * target_a
    m22 <- (m221 + m222) / (alpha * kd1^2 * kd2)
    m23 <- -(central * target_a - alpha * central * target_a) / (alpha * kd1 * kd2)

    m31 <- -(central * target_b^2 + kd2 * target_a * target_b +
      alpha * central * kd2 * target_b + alpha * kd1 * kd2 * target_b) / (alpha * kd1 * kd2^2)
    m32 <- -(central * target_b - alpha * central * target_b) / (alpha * kd1 * kd2)
    m331 <- central * target_b^2 + central * kd2 * target_b +
      kd2 * target_a * target_b + alpha * central * kd2^2 + alpha * kd1 * kd2^2
    m332 <- alpha * kd2^2 * target_a + alpha * central * kd2 * target_b +
      alpha * kd1 * kd2 * target_b
    m33 <- (m331 + m332) / (alpha * kd1 * kd2^2)

    # Eq 19: the vector g_BsAb(AD, C, R_A, R_B, AP).
    g1 <- in_iv / vc + ka * depot_sc / vc - kel * central -
      kinta * (target_a * central) / kd1 -
      kintb * (target_b * central) / kd2 -
      kintab * (target_a * target_b * central) / (alpha * kd1 * kd2) -
      k12 * central + k21 * peripheral1 / vc
    g2 <- ksyna - kdega * target_a -
      kinta * (target_a * central) / kd1 -
      kintab * (target_a * target_b * central) / (alpha * kd1 * kd2)
    g3 <- ksynb - kdegb * target_b -
      kintb * (target_b * central) / kd2 -
      kintab * (target_a * target_b * central) / (alpha * kd1 * kd2)

    # Eq 16: the matrix-vector product that replaces the three free-species ODEs.
    d/dt(central) <- (m11 / det) * g1 + (m12 / det) * g2 + (m13 / det) * g3
    d/dt(target_a) <- (m21 / det) * g1 + (m22 / det) * g2 + (m23 / det) * g3
    d/dt(target_b) <- (m31 / det) * g1 + (m32 / det) * g2 + (m33 / det) * g3

    # Eq 17: peripheral compartment, carrying an AMOUNT.
    d/dt(peripheral1) <- k12 * central * vc - k21 * peripheral1

    f(depot_sc) <- fdepot

    # Eqs 20-22: the complexes are ALGEBRAIC under rapid binding, not states.
    complexA <- central * target_a / kd1
    complexB <- central * target_b / kd2
    trimerAB <- central * target_a * target_b / (alpha * kd1 * kd2)

    # Eqs 11-13.
    Cc <- central
    CcTotal <- central + complexA + complexB + trimerAB
    RtotA <- target_a + complexA + trimerAB
    RtotB <- target_b + complexB + trimerAB

    Cc ~ prop(propSd)
  })
}
