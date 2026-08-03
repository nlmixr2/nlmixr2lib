Noe_1996_factor_viii <- function() {
  description <- "Deterministic mechanistic model of coagulation factor VIII kinetics with von Willebrand factor (vWF) binding equilibrium in adult humans; captures endogenous synthesis, reversible binding to vWF, and differential elimination of the free and vWF-bound forms of factor VIII"
  reference <- "Noe DA. A mathematical model of coagulation factor VIII kinetics. Haemostasis. 1996;26(6):289-303. doi:10.1159/000217222. PMID: 8979143."
  vignette <- "Noe_1996_factor_viii"
  paper_specific_compartments <- c("viii", "vwf")
  units <- list(
    time = "h",
    dosing = "U/mL (mass dose divided by plasma volume; state is a concentration, not amount)",
    concentration = "U/mL (both factor VIII and vWF measured on the conventional 1 U/mL == 100% of normal-pool scale)"
  )

  population <- list(
    species       = "human",
    disease_state = "healthy adults with normal endogenous factor VIII and vWF; the same structural model is exercised in the paper to describe hemophilia (reduced factor VIII synthesis), hemophilia carriers, type 1 and type 3 von Willebrand disease (reduced or absent vWF), acute-phase reactions and pregnancy (elevated vWF), and hemophilia / vWD replacement therapy with cryoprecipitate, recombinant factor VIII, or high-purity vWF concentrate",
    notes         = "Parameter values are the final recommended combination reported in the Discussion (p. 300): Svm/kfviii = 0.075 U/mL, kfviii = 1.12 /h, kbviii = 0.028 /h, kequi = 0.2 U/mL, N = 4.67 U factor VIII / U vWF. Fit to mean literature data from Holmberg & Nilsson (subacute/chronic illness), Hermens (hemophilia treated with cryoprecipitate), Fijnvandraat (recombinant factor VIII in hemophilia), Cattaneo (recombinant factor VIII in vWD), Goudemand (high-purity vWF in vWD), and Morfini (type 3 vWD). The paper reports no individual-subject fits; the model is a typical-value mechanism with no inter-individual variability or residual error."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    viii = list(analyte = "Factor VIII", units = NA_character_, specimen = "plasma", verified = FALSE),
    vwf  = list(analyte = "Von Willebrand factor (vWF)", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  ini({
    # Rate constants and structural constants -- Noe 1996 Discussion p. 300.
    # "The kinetic parameter values that, in combination, give the best modeling
    # results are Svm/kfviii = 0.075 U/ml, kfviii = 1.12 h-1, kbviii = 0.028 h-1,
    # kequi = 0.2 U factor VIII/ml, and N = 4.67 U factor VIII/U vWF."
    lkfviii    <- fixed(log(1.12));  label("Elimination rate constant of unbound (free) factor VIII (1/h); half-life ~37 min")               # Noe 1996 Discussion p. 300 (kfviii = 1.12 /h)
    lkbviii    <- fixed(log(0.028)); label("Elimination rate constant of vWF-bound factor VIII; equal to kvWF (1/h); half-life ~24.5 h")     # Noe 1996 Discussion p. 300 (kbviii = 0.028 /h; also equals kvWF by construction, Model section)
    lkd        <- fixed(log(0.2));   label("Dissociation constant of the factor VIII-vWF complex, approximated by kequi (U factor VIII/mL)") # Noe 1996 Simulations p. 292 + Discussion p. 300 (kd = kequi = 0.2 U/mL; kd from in vitro binding studies refs [16,17])
    lnvwf      <- fixed(log(4.67));  label("Proportionality N between vWF concentration and factor VIII binding-site concentration on multimeric vWF (U factor VIII / U vWF)") # Noe 1996 Discussion p. 300 (N = 4.67; multimeric vWF; corresponds to ~1 of every 11 vWF monomers being available for factor VIII binding, p. 300)
    lsvmkfviii <- fixed(log(0.075)); label("Modified-model steady-state parameter Svm/kfviii: plasma factor VIII concentration in the absence of plasma vWF (U/mL)") # Noe 1996 Discussion p. 300 (Svm/kfviii = 0.075 U/mL); also the type 3 vWD asymptote
    # Baseline steady-state concentrations (used to compute endogenous vWF
    # synthesis rate and as ODE initial conditions). Users override these
    # to simulate non-normal populations (see the vignette). Kept on the
    # linear scale so bl_vwf = 0 (type 3 vWD) is expressible.
    bl_vwf     <- fixed(1.0);        label("Baseline plasma vWF concentration used as the vWF steady-state anchor (U/mL); 1.0 for a normal healthy adult, 0 for severe (type 3) vWD") # Noe 1996 Steady State p. 294 ("the concentration of factor VIII equals 1 U/ml when the concentration of vWF equals 1 U/ml, in keeping with conventional practice")
    bl_viii    <- fixed(1.0);        label("Baseline plasma total factor VIII concentration used as the factor VIII initial condition (U/mL); 1.0 for a normal healthy adult, 0.075 (=Svm/kfviii) for severe (type 3) vWD") # Noe 1996 Steady State p. 294 + Discussion p. 300
  })

  model({
    # Un-transform rate constants and structural constants.
    kfviii    <- exp(lkfviii)
    kbviii    <- exp(lkbviii)          # = kvwf by construction (Model section p. 291)
    kd        <- exp(lkd)              # = kequi under the assumption kvwf << koff, kon
    nvwf      <- exp(lnvwf)
    svmkfviii <- exp(lsvmkfviii)

    # Endogenous synthesis rates (concentration units per time).
    # svm is derived from the modified-model parameter Svm/kfviii * kfviii.
    # svwf is derived to hold vWF at its steady-state value bl_vwf, using
    # the assumption that vWF elimination is first-order at rate kvWF = kbviii
    # (Model section p. 291: "The rate constant of elimination of vWF (kvwF)
    # is assumed to be the same for both the bound and the free forms").
    svm  <- svmkfviii * kfviii
    svwf <- kbviii * bl_vwf

    # Instantaneous mass-action binding equilibrium.
    # The paper's kequi is very close to kd (Simulations p. 292), so binding
    # partition between free and vWF-bound factor VIII is set by kd at every
    # instant given total factor VIII (viii) and total vWF (vwf):
    #
    #     [fVIII] * [fSites] = kd * [bVIII]
    #     [fVIII] + [bVIII]  = viii
    #     [fSites] + [bSites] = N * vwf,  with [bSites] = [bVIII]
    #
    # Substituting yields the quadratic
    #     bVIII^2 - bVIII * (viii + N*vwf + kd) + viii * N * vwf = 0
    # whose physical (smaller) root is:
    ssum  <- viii + nvwf * vwf + kd
    disc  <- ssum * ssum - 4.0 * viii * nvwf * vwf
    bviii <- (ssum - sqrt(disc)) / 2.0
    fviii <- viii - bviii

    # Differential equations for total factor VIII and total vWF.
    # Because bviii and vwf are eliminated at kbviii = kvwf, the sum of
    # the paper's four ODEs collapses to two equations for the totals.
    # (Model section p. 291; Discussion p. 300 makes the reduction explicit
    # via the paper's steady-state derivation on p. 294-296.)
    d/dt(viii) <- svm  - kfviii * fviii - kbviii * bviii
    d/dt(vwf)  <- svwf - kbviii * vwf

    # Initial conditions: total factor VIII and total vWF at their reported
    # steady-state values for the target population.
    viii(0) <- bl_viii
    vwf(0)  <- bl_vwf

    # Auxiliary derived quantities (algebraic; carried out at every solver
    # step for reporting and for the vignette's mass-balance checks).
    total_viii <- viii
    total_vwf  <- vwf
    free_viii  <- fviii
    bound_viii <- bviii
  })
}
