Almquist_2015_Mig1_exp4 <- function() {
  description <- paste(
    "In vitro (yeast, Saccharomyces cerevisiae).",
    "Phenomenological two-state perfect-adaptation model of nuclear Mig1",
    "dynamics in single yeast cells after an instantaneous downshift in",
    "extracellular glucose, fitted by nonlinear mixed effects (FOCE) to",
    "time-lapse fluorescence-microscopy data from 46 individual cells in",
    "EXPERIMENT 4 (4% -> 0.5% glucose, the lowest post-shift level at which",
    "the transient response is still seen). Nuclear Mig1 (observed as",
    "Mig1-GFP light intensity) is produced at a glucose-proportional rate and",
    "removed at a rate modulated by a lumped, unmeasured adaptation component",
    "X; the circuit adapts perfectly, so nuclear Mig1 returns to its basal",
    "level regardless of the post-shift glucose level. Cell-to-cell",
    "variability is carried by a full 3x3 log-normal random-effect block on",
    "the basal Mig1 level and the two rate constants. One of four",
    "per-experiment fits; see modellib('Almquist_2015_Mig1_exp1') and",
    "modellib('Almquist_2015_Mig1_exp2') (4% -> 1.5%), and",
    "modellib('Almquist_2015_Mig1_exp3') (4% -> 1%)."
  )
  reference <- paste(
    "Almquist J, Bendrioua L, Adiels CB, Goksor M, Hohmann S, Jirstrand M.",
    "(2015). A Nonlinear Mixed Effects Approach for Modeling the",
    "Cell-To-Cell Variability of Mig1 Dynamics in Yeast.",
    "PLoS ONE 10(4):e0124050.",
    "doi:10.1371/journal.pone.0124050."
  )
  vignette <- "Almquist_2015_Mig1"

  # `mig1` (nuclear Mig1-GFP light intensity) and `adapt` (the paper's lumped
  # adaptation component X) are paper-mechanistic states with no canonical
  # analogue; declared here per the compartment-naming reference. Closest
  # in-library precedent: the Bizzotto 2016 glucose-insulin phase states and
  # the Denti 2010 `insulin_action` state.
  paper_specific_compartments <- c("mig1", "adapt")

  units <- list(
    time = "s",
    dosing = "not applicable (no drug is administered; the input is a step change in extracellular glucose)",
    concentration = "fluorescence light intensity (a.u.); `adapt` is dimensionless"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are latent / non-matrix quantities in a
  # single yeast cell, so `specimen` is "not applicable" for both.
  compartmentData <- list(
    mig1  = list(analyte = "nuclear Mig1-GFP", units = "fluorescence light intensity (a.u.)", specimen = "not applicable", verified = TRUE),
    adapt = list(analyte = "lumped adaptation component X (unidentified)", units = "dimensionless", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "in vitro (yeast, Saccharomyces cerevisiae)",
    n_subjects     = 46L,
    n_studies      = 1L,
    disease_state  = paste(
      "Glucose-grown Mig1-GFP expressing S. cerevisiae cells, held in a",
      "microfluidic device by optical tweezers and imaged by fluorescence",
      "microscopy (Almquist 2015 Methods; experimental setup previously",
      "published as reference 39 of the source paper, Bendrioua 2014",
      "J Biol Chem 289:12863)."
    ),
    dose_range     = paste(
      "No drug. Extracellular glucose is stepped instantaneously from 4% to",
      "0.5% at time 0 (Almquist 2015 Table 1, Exp Nr 4). Up to 15",
      "fluorescence observations per cell."
    ),
    notes          = paste(
      "Almquist 2015 Table 1: experiment 4 used 46 cells shifted from 4% to",
      "0.5% glucose. This experiment is kinetically distinct from the other",
      "three: k4 is roughly double and k2 slightly larger (Table 2), which",
      "the authors read as evidence that 0.5% glucose is close to a threshold",
      "in the Snf1-Mig1 system -- the transient behaviour disappears entirely",
      "below 0.5% (Discussion, 'Parameter estimates'). Consistent with this,",
      "when all four data sets were fitted simultaneously the empirical Bayes",
      "estimates from experiment 4 formed a separate cluster (Fig 5), which",
      "violates the assumption that the random effects are identically",
      "distributed across experiments and is the reason the separate",
      "per-experiment analyses are the trustworthy ones."
    )
  )

  ini({
    # =====================================================================
    # Fixed effects -- Almquist 2015 Table 2, column "Exp 4".
    # Relative standard errors from the same table are quoted per line.
    # The paper models Ms, k2 and k4 as log-normal (Ms = Ms_hat * exp(eta1)
    # etc.), which is exactly the log-transformed-theta plus additive-eta
    # parameterization used here.
    # =====================================================================

    lrbase <- log(3.14e3)
    label("Basal nuclear Mig1 level Ms (fluorescence intensity, a.u.)")
    # Almquist 2015 Table 2, Exp 4: Ms = 3.14 x 10^3 (RSE 1%).
    # Ms is both the pre-shift steady-state level of nuclear Mig1 and, by
    # perfect adaptation, the level the system returns to after the shift.

    lkout_mig1 <- log(0.00815)
    label("First-order loss rate constant of nuclear Mig1, k2 (1/s)")
    # Almquist 2015 Table 2, Exp 4: k2 = 0.00815 /s (RSE 7%) -- the largest
    # of the four experiments. Paper text: "k2 will determine the
    # turnover-timescale of Mig1(t)". Enters the Mig1 mass balance as the
    # X-modulated loss term k2 * X(t) * Mig1(t); X is dimensionless so k2
    # has units 1/s.

    lkout_adapt <- log(0.0229)
    label("First-order relaxation rate constant of the adaptation component X, k4 (1/s)")
    # Almquist 2015 Table 2, Exp 4: k4 = 0.0229 /s (RSE 8%) -- roughly
    # double the value in experiments 1 to 3, the signature of the 0.5%
    # threshold behaviour discussed in the paper. k4 sets the relaxation
    # time-scale of X and appears in BOTH the production and the loss term
    # of the X mass balance after the identifiability reduction (see
    # model() for the derivation).

    # =====================================================================
    # Experimental input levels -- Almquist 2015 Table 1 and the definition
    # Glu(t) = 4 - (4 - g) * H(t). These are controlled experimental
    # settings, not estimated quantities, hence fixed().
    # =====================================================================

    glu_pre <- fixed(4)
    label("Pre-shift extracellular glucose, Glu(-30) (% w/v)")
    # Almquist 2015 Table 1, column "From": all four experiments start at 4%.

    glu_post <- fixed(0.5)
    label("Post-shift extracellular glucose, g (% w/v)")
    # Almquist 2015 Table 1, Exp Nr 4, column "To": 0.5%.

    # =====================================================================
    # Random effects -- Almquist 2015 Table 3, Exp Nr 4 (the covariance
    # matrix Omega for eta = (eta1, eta2, eta3), associated with Ms, k2 and
    # k4 respectively per the Table 3 footnote).
    #
    # The values below are reproduced to FULL precision from the Table 2
    # decomposition parameters rather than from the two-significant-figure
    # Table 3 print. Methods, "Parameterization of the random effect
    # covariance matrix": Omega = U %*% t(U) with U upper triangular,
    #   U = [[w11, w12, w13], [0, w22, w23], [0, 0, w33]].
    # Substituting the Table 2 Exp 4 values
    #   w11 = 0.0228, w12 = 0.0691, w13 = -0.0322,
    #   w22 = 0.252,  w23 = 0.281,  w33 = 0.432
    # reproduces every printed entry of Table 3 Exp Nr 4 exactly:
    #   Omega11 = 0.0063315 (Table 3: 0.0063)
    #   Omega21 = 0.0083650 (Table 3: 0.0083)
    #   Omega22 = 0.1424650 (Table 3: 0.14)
    #   Omega31 = -0.0139104 (Table 3: -0.014)
    #   Omega32 = 0.1213920 (Table 3: 0.12)
    #   Omega33 = 0.1866240 (Table 3: 0.19)
    # and the implied correlations 0.28 (eta1,eta2), -0.41 (eta1,eta3) and
    # 0.74 (eta2,eta3) match the Table 3 correlation block. The construction
    # Omega = U %*% t(U) guarantees positive definiteness.
    # =====================================================================

    etalrbase + etalkout_mig1 + etalkout_adapt ~
      c( 0.0063314900,
         0.0083650000, 0.1424650000,
        -0.0139104000, 0.1213920000, 0.1866240000)

    # =====================================================================
    # Residual error -- Almquist 2015 Table 2, row "s", Exp 4.
    # =====================================================================

    addSd_mig1 <- sqrt(24.1e3)
    label("Additive residual SD on nuclear Mig1 intensity (a.u.)")
    # Almquist 2015 Table 2, Exp 4: s = 24.1 x 10^3 (RSE 6%). The paper
    # defines the observation model as y_t = Mig1(t) + e_t with
    # "e_t ~ N(0, s), with s denoting the VARIANCE of the measurement
    # error", so the nlmixr2 additive SD is sqrt(s) = 155.2 a.u.
    # Corroborated by the Methods section "Starting values for the
    # optimization algorithm": "The measurement noise appear to be on the
    # scale of a hundred to a few hundreds and its variance, s, was set to
    # 40 000" -- i.e. sqrt(40000) = 200 a.u., a starting SD in the stated
    # "hundred to a few hundreds" range.
  })

  model({
    # =====================================================================
    # 1. Extracellular glucose input (experimentally controlled).
    #    Almquist 2015: Glu(t) = 4 - (4 - g) * H(t), where H is the
    #    Heaviside step function and the shift happens at time 0. The
    #    integration convention is that the experiment starts at
    #    t = -30 s, where the system is at steady state.
    # =====================================================================
    glu <- glu_pre - (glu_pre - glu_post) * (time >= 0)

    # =====================================================================
    # 2. Individual parameters. Almquist 2015:
    #      Ms = Ms_hat * exp(eta1), k2 = k2_hat * exp(eta2),
    #      k4 = k4_hat * exp(eta3),  eta ~ N(0, Omega).
    # =====================================================================
    rbase      <- exp(lrbase      + etalrbase)
    kout_mig1  <- exp(lkout_mig1  + etalkout_mig1)
    kout_adapt <- exp(lkout_adapt + etalkout_adapt)

    # =====================================================================
    # 3. Zero-order Mig1 production. Almquist 2015 reparameterizes k1 in
    #    terms of the other parameters using the pre-shift steady state
    #    (0 = k1 * Glu(-30) - k2 * Xs * Ms with Xs = 1):
    #      k1 = k2 * Ms / Glu(-30).
    #    This makes the basal level Ms an explicit model parameter.
    # =====================================================================
    kin_mig1 <- kout_mig1 * rbase / glu_pre

    # =====================================================================
    # 4. Mass balances for the reduced (structurally identifiable) model.
    #    Almquist 2015, "the equations defining the model in Fig 2":
    #      dMig1/dt = k1 * Glu(t) - k2 * X(t) * Mig1(t)
    #      dX/dt    = k4 * Glu(t) / Glu(-30) - k4 * X(t)
    #    The original model had a separate k3 in the X production term;
    #    the substitution X~ = alpha * X with alpha = k4 / (k3 * Glu(-30))
    #    removes k3 without changing the observable Mig1(t) (paper's
    #    identifiability reduction), which is why k4 appears in both the
    #    production and the loss term of the X equation.
    #
    #    Dimensional check (state units a.u., time s, glucose %):
    #      kin_mig1 * glu          = [a.u./(s*%)] * [%]  = a.u./s   OK
    #      kout_mig1 * adapt * mig1= [1/s] * [-] * [a.u.]= a.u./s   OK
    #      kout_adapt * glu/glu_pre= [1/s] * [-]         = 1/s      OK
    # =====================================================================
    d/dt(mig1)  <- kin_mig1 * glu - kout_mig1 * adapt * mig1
    d/dt(adapt) <- kout_adapt * (glu / glu_pre) - kout_adapt * adapt

    # =====================================================================
    # 5. Initial conditions at t = -30 s. Almquist 2015: Mig1(-30) = Ms and
    #    X(-30) = Xs, with the steady-state constraint forcing Xs = 1.
    #    Because the pre-shift state IS the steady state, starting the
    #    integration at any time before 0 gives the same trajectory.
    # =====================================================================
    mig1(0)  <- rbase
    adapt(0) <- 1

    # =====================================================================
    # 6. Observation. Almquist 2015: y_t = Mig1(t) + e_t, e_t ~ N(0, s)
    #    with s the variance; addSd_mig1 = sqrt(s). Fluorophore bleaching
    #    was deliberately NOT modelled (paper: the average recovery level
    #    at 20 min was 96% of the pre-shift intensity).
    # =====================================================================
    mig1 ~ add(addSd_mig1)
  })
}
