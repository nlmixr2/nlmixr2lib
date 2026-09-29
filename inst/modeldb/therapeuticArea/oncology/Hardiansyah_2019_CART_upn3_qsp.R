# Quantitative systems pharmacology (QSP) model of anti-CD19 chimeric
# antigen receptor T-cell (CART) therapy and the pro-inflammatory
# cytokines associated with cytokine release syndrome (CRS), published
# by Hardiansyah & Ng 2019 (Clin Transl Sci 12:343-349).
#
# This file carries the UPN3 parameter set. Its sibling
# Hardiansyah_2019_CART_upn1_qsp.R carries UPN1. Table S2 of the
# supplement reports a complete, independent set of seven estimated
# parameters per subject, so the two subjects are extracted as two files
# per references/replicate-author-structure.md ("per subpopulation ->
# Author_Year_<drug>_<population>.R"). Both point at one vignette.
#
# UPN3 is the low-burden / low-dose contrast to UPN1: a 4.5-fold lower
# baseline disease burden and a 79-fold lower CART dose, which the paper
# uses to show that cytokine elevation tracks disease burden rather than
# administered dose.
#
# SOURCES. The main article is a Brief Report and contains no equations
# and no parameter values. Every structural equation (Eqs. 1-9) and both
# parameter tables live in the Supplementary File, distributed in the
# EuropePMC open-access package for PMC6662387 as CTS-12-343-s001.docx.
# Table S1 holds the literature-sourced (fixed) parameters, Table S2 the
# estimated ones. The cytokine baselines are in neither table and were
# digitized from Figure 1b; see the vignette Errata. The three
# unit/typographic defects resolved in the model body are documented in
# the sibling file and in the vignette Errata.

Hardiansyah_2019_CART_upn3_qsp <- function() {
  description <- "QSP. Anti-CD19 CAR T-cell (CART) therapy and cytokine release syndrome in advanced chronic lymphocytic leukaemia, subject UPN3 (low disease burden, low dose). Eight ODE states: CLL B cells in peripheral blood, effector and memory CART each distributed between peripheral blood and tissue, and IL-6, IL-10 and IFN-gamma. CART expansion is driven by the B-cell burden (a target-mediated 'living drug' structure in which the cells grow in the presence of target rather than being cleared by it), effector CART kill B cells and drive cytokine secretion, and IL-10 inhibits IFN-gamma secretion through a decreasing Hill function. UPN3 received a 79-fold lower CART dose and carried a 4.5-fold lower baseline disease burden than UPN1. Deterministic: all parameters are fixed at the published point estimates and the paper reports no between-subject or residual variability."
  reference <- paste(
    "Hardiansyah D, Ng CM (2019).",
    "Quantitative systems pharmacology model of chimeric antigen receptor T-cell therapy.",
    "Clinical and Translational Science 12(4):343-349.",
    "doi:10.1111/cts.12636.",
    "Equations 1-9 and Tables S1-S2 are in the Supplementary File",
    "(CTS-12-343-s001.docx in the PMC6662387 open-access package).",
    "Fit to CART and cytokine kinetics digitized from",
    "Kalos et al. 2011 (Sci Transl Med 3:95ra73; doi:10.1126/scitranslmed.3002842).",
    sep = " "
  )
  vignette <- "Hardiansyah_2019_CART_crs_qsp"

  # The CLL B-cell pool and the four CAR T-cell phenotype/site states are
  # paper-mechanistic QSP compartments with no canonical analogue in
  # inst/references/compartment-names.md; see that file's
  # "Paper-specific compartments" section. The three cytokine states are
  # canonical members of the inflammatory-mediator family (`il6` was
  # already registered; `il10` and `ifng` are registered in this PR).
  paper_specific_compartments <- c(
    "b_pb",
    "carte_pb",
    "carte_t",
    "cartm_pb",
    "cartm_t"
  )

  units <- list(
    time = "day",
    dosing = "1e9 cells (anti-CD19 CART; input as amt on the carte_pb compartment)",
    concentration = "pg/mL (IL-6, IL-10, IFN-gamma); cell states are counts in units of 1e9 cells"
  )

  compartmentData <- list(
    b_pb = list(
      analyte = "CLL B cells",
      units = "1e9 cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_pb = list(
      analyte = "activated (effector) anti-CD19 CAR T-cells",
      units = "1e9 cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    carte_t = list(
      analyte = "activated (effector) anti-CD19 CAR T-cells",
      units = "1e9 cells",
      specimen = "tissue",
      verified = TRUE
    ),
    cartm_pb = list(
      analyte = "memory anti-CD19 CAR T-cells",
      units = "1e9 cells",
      specimen = "whole blood",
      verified = TRUE
    ),
    cartm_t = list(
      analyte = "memory anti-CD19 CAR T-cells",
      units = "1e9 cells",
      specimen = "tissue",
      verified = TRUE
    ),
    il6 = list(analyte = "interleukin-6", units = "pg/mL", specimen = "whole blood", verified = TRUE),
    il10 = list(analyte = "interleukin-10", units = "pg/mL", specimen = "whole blood", verified = TRUE),
    ifng = list(analyte = "interferon gamma", units = "pg/mL", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1,
    n_studies = 1,
    disease_state = "advanced, chemotherapy-resistant chronic lymphocytic leukaemia (CLL)",
    dose_range = "1.4e7 anti-CD19 CART cells total, split 10% / 30% / 60% over days 0, 1 and 2",
    disease_burden = "9.3e8 CLL B cells in peripheral blood at baseline",
    notes = paste(
      "Subject UPN3 of the three-patient CLL cohort reported by Kalos et al. 2011.",
      "The QSP model was developed on UPN1 and UPN3 only; UPN2 was excluded by the",
      "authors because concurrent corticosteroid therapy masked the cytokine response",
      "(Hardiansyah 2019 Methods). Baseline peripheral-blood B-cell counts were derived",
      "in the supplement from total white cell counts assuming 80% lymphocytes and 50%",
      "B cells among lymphocytes. This file carries the UPN3 parameter set;",
      "modellib('Hardiansyah_2019_CART_upn1_qsp') carries UPN1.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Table S1 -- parameters taken from the published literature (fixed).
    # Identical across subjects.
    # ------------------------------------------------------------------
    ldil6 <- fixed(log(3.96))
    label("IL-6 natural elimination rate constant (1/day)") # Table S1, d_IL6 = 3.96 1/d
    ldil10 <- fixed(log(4.68))
    label("IL-10 natural elimination rate constant (1/day)") # Table S1, d_IL10 = 4.68 1/d
    ldifng <- fixed(log(1.51))
    label("IFN-gamma natural elimination rate constant (1/day)") # Table S1, d_IFNgamma = 1.51 1/d
    ag <- fixed(0.59)
    label("Asymptotic residual IFN-gamma secretion under maximal IL-10 inhibition (unitless)") # Table S1, a_G = 0.59
    lbg <- fixed(log(650))
    label("IL-10 concentration giving half-maximal inhibition of IFN-gamma secretion (pg/mL)") # Table S1, b_G = 650; see Errata on the printed 'ng' unit
    lrm <- fixed(log(5.09e-2))
    label("Memory CART natural growth rate constant (1/day)") # Table S1, r_M = 5.09e-2 1/d
    ldcartm <- fixed(log(4.08e-1))
    label("Memory CART death rate constant in peripheral blood (1/day)") # Table S1, d_CARTM = 4.08e-1 1/d
    lam <- fixed(log(1))
    label("Activation rate constant, memory CART to effector CART (1/day)") # Table S1, a_M = 1 1/d
    lkin <- fixed(log(76.78))
    label("CART distribution rate constant, peripheral blood to tissue (1/day)") # Table S1, k_in = 76.78 1/d
    lkout <- fixed(log(1.84))
    label("CART distribution rate constant, tissue to peripheral blood (1/day)") # Table S1, k_out = 1.84 1/d
    lh <- fixed(log(1e5))
    label("Half-saturation constant of the B-cell saturation function f(B) (1e9 cells)") # Table S1, h = 1e5; see Errata on the printed 'Cells/ml' unit
    lrb <- fixed(log(9.30e-3))
    label("CLL B-cell natural reproduction rate constant (1/day)") # Table S1, r_B = 9.30e-3 1/d
    ldb <- fixed(log(7.60e-3))
    label("CLL B-cell natural death rate constant (1/day)") # Table S1, d_B = 7.60e-3 1/d

    # ------------------------------------------------------------------
    # Table S2 -- parameters estimated by fitting the QSP model to the
    # observed CART and cytokine kinetics. Values are the UPN3 column;
    # the parenthesised percentages in Table S2 are estimation precision
    # (%RSE), not between-subject variability, so no IIV is encoded.
    # ------------------------------------------------------------------
    lpil6 <- log(0.6e-2)
    label("CART-driven IL-6 production rate constant (pg/uL per 1e15 cell^2 per day)") # Table S2, P_IL6 UPN3 = 0.6e-2 (33% RSE)
    lpil10 <- log(1.7e-2)
    label("CART-driven IL-10 production rate constant (pg/uL per 1e15 cell^2 per day)") # Table S2, P_IL10 UPN3 = 1.7e-2 (27% RSE)
    lpifng <- log(0.6e-2)
    label("CART-driven IFN-gamma production rate constant (pg/uL per 1e15 cell^2 per day)") # Table S2, P_IFNgamma UPN3 = 0.6e-2 (31% RSE)
    lre <- log(14.40)
    label("Effector CART replication rate constant (1 per 1e9 cells per day)") # Table S2, r_E UPN3 = 14.40 (5% RSE)
    ldcarte <- log(1.44)
    label("Effector CART elimination rate constant (1/day)") # Table S2, d_CARTE UPN3 = 1.44 (27% RSE)
    lae <- log(1.2e-1)
    label("Activation rate constant, effector CART to memory CART (1/day)") # Table S2, a_E UPN3 = 1.2e-1 (11% RSE); see Errata on the printed unit
    lkbc <- log(2.21)
    label("Effector-CART-mediated CLL B-cell killing rate constant (1 per 1e9 cells per day)") # Table S2, K_BC UPN3 = 2.21 (21% RSE)

    # ------------------------------------------------------------------
    # Not reported in Table S1, Table S2 or the main text.
    # Digitized from the UPN3 cytokine panels of Figure 1b (the fitted
    # line's post-therapy plateau, which by the supplement's definition
    # equals P_endo / d). Read off a six-decade log axis; treat as
    # approximately +/- 30%. See the vignette Errata.
    # ------------------------------------------------------------------
    bl_il6 <- fixed(1)
    label("Baseline IL-6 concentration (pg/mL) -- figure-derived") # Figure 1b, UPN3 IL-6 panel plateau
    bl_il10 <- fixed(0.5)
    label("Baseline IL-10 concentration (pg/mL) -- figure-derived") # Figure 1b, UPN3 IL-10 panel plateau
    bl_ifng <- fixed(1.7)
    label("Baseline IFN-gamma concentration (pg/mL) -- figure-derived") # Figure 1b, UPN3 IFN-gamma panel plateau

    lb0 <- fixed(log(0.93))
    label("Baseline CLL B-cell disease burden in peripheral blood (1e9 cells)") # Figure 2 caption, 'disease burden baseline in PB: ... 9.3 x 10^8 cells (UPN3)'
  })

  model({
    # ---- back-transform ----
    dil6 <- exp(ldil6)
    dil10 <- exp(ldil10)
    difng <- exp(ldifng)
    bg <- exp(lbg)
    rm_ <- exp(lrm)
    dcartm <- exp(ldcartm)
    am <- exp(lam)
    kin <- exp(lkin)
    kout <- exp(lkout)
    h <- exp(lh)
    rb <- exp(lrb)
    db <- exp(ldb)
    pil6 <- exp(lpil6)
    pil10 <- exp(lpil10)
    pifng <- exp(lpifng)
    re <- exp(lre)
    dcarte <- exp(ldcarte)
    ae <- exp(lae)
    kbc <- exp(lkbc)

    # ---- initial conditions ----
    # CART states start empty (the supplement assumes no memory cells at
    # time zero); the infusion enters carte_pb as a dose record.
    b_pb(0) <- exp(lb0)
    il6(0) <- bl_il6
    il10(0) <- bl_il10
    ifng(0) <- bl_ifng

    # ---- structural functions ----
    # Eq. 9: B-cell saturation function. Governs the direction of the
    # effector/memory interconversion: while burden is high (f(B) -> 1)
    # memory cells are mobilised into effectors; as burden falls
    # (f(B) -> 0) effectors convert to long-lived memory cells.
    fb <- b_pb / (b_pb + h)
    # Eq. 4 bracket: decreasing function for IL-10 inhibition of
    # IFN-gamma secretion. Fd(0) = 1, Fd(inf) = a_G, Fd(b_G) = (a_G+1)/2.
    fd_il10 <- ag + (1 - ag) * bg / (il10 + bg)

    # ---- ODE system (Supplementary File Eqs. 1-8) ----
    # Eq. 1: CLL B cells in peripheral blood.
    d/dt(b_pb) <- rb * b_pb - db * b_pb - kbc * carte_pb * b_pb

    # Eq. 5: activated (effector) CART in peripheral blood. The infusion
    # term D_inj is supplied as a dose record on this compartment.
    d/dt(carte_pb) <- re * carte_pb * b_pb - dcarte * carte_pb -
      kin * carte_pb + kout * carte_t -
      ae * carte_pb * (1 - fb) + am * cartm_pb * fb
    # Eq. 6: activated CART in tissue.
    d/dt(carte_t) <- kin * carte_pb - kout * carte_t

    # Eq. 7: memory CART in peripheral blood. The gain term is
    # + a_E * carte_pb * (1 - f(B)), i.e. the exact flux Eq. 5 loses.
    # Eq. 7 as printed reads + a_E * CARTM_PB * (1 - f(B)); that is a
    # typographical error -- it breaks the transfer pairing that the
    # a_M terms obey, and with CARTM(0) = 0 and no other source it makes
    # cartm_pb identically zero for all time, contradicting both the
    # prose ("Concentration of CARTM increased as the activated CARTE
    # were converted into CARTM via activation rate a_E") and the
    # non-zero CART Memory curve in Figure 1b. Figure 1a labels this
    # arrow a_E*CARTE_PB*(1-f(B)) and is the arbiter. See Errata.
    d/dt(cartm_pb) <- rm_ * cartm_pb - dcartm * cartm_pb -
      kin * cartm_pb + kout * cartm_t +
      ae * carte_pb * (1 - fb) - am * cartm_pb * fb
    # Eq. 8: memory CART in tissue.
    d/dt(cartm_t) <- kin * cartm_pb - kout * cartm_t

    # Eqs. 2-4: pro-inflammatory cytokines. P_endo = d * baseline, per
    # the supplement's definition, so each cytokine starts and returns to
    # its baseline. The 1e6 factor converts the Table S2 production
    # constants (pg/uL per 1e15 cell^2 per day) to pg/mL per day given
    # cell states in 1e9 cells: 1e9 * 1e9 / 1e15 * 1e3 = 1e6.
    d/dt(il6) <- dil6 * bl_il6 + pil6 * 1e6 * b_pb * carte_pb - dil6 * il6
    d/dt(il10) <- dil10 * bl_il10 + pil10 * 1e6 * b_pb * carte_pb - dil10 * il10
    d/dt(ifng) <- difng * bl_ifng + pifng * 1e6 * b_pb * carte_pb * fd_il10 -
      difng * ifng

    # ---- reported outputs ----
    # Total CART in peripheral blood: the quantity plotted as 'Predicted
    # CART' in Figure 1b and compared against the observed symbols.
    cart_pb <- carte_pb + cartm_pb
  })
}
