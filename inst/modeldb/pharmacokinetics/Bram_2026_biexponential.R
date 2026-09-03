Bram_2026_biexponential <- function() {
  description <- "One-state (no latent compartment) intravenous-bolus model that reproduces bi-exponential plasma kinetics, proposed automatically by the NODE-LASSO workflow of Braem 2026 (Equation 20). Elimination from the single central state is the sum of a conventional first-order term and a transient, dose-proportional flux that decays exponentially with time after dose: dA/dt = -kel * A - D * wdist * exp(-kdist * (t - tD)). The authors prove (Equation 21 and Supporting Information) that its explicit solution is algebraically IDENTICAL to the central compartment of a standard two-compartment model, so the transient flux stands in for the distribution phase without carrying a peripheral state. The data were SIMULATED from a two-compartment model for 50 subjects with 7 concentration measurements over 24 h; no drug, no dose and no units are reported, and the fit identifies the model only up to the ratio D/Vc = 10.1. A nominal unit dose is therefore fixed here (vc = 1/10.1), so a dose of 1 arbitrary unit reproduces the published macro-parameters exactly and any other dose scales the profile linearly, which is exact for this linear model. The companion model Bram_2026_biexponential_2cmt is the two-compartment reference fit to the same data. Because it has no memory of prior doses, this parameterisation is valid for multiple dosing ONLY when each dose falls after the previous distribution phase has completed; Braem 2026 Figure S1 shows the structural misspecification that appears otherwise, and Figure S2 gives an equivalent pseudo-compartment form that repairs it."
  reference <- paste(
    "Braem DS, Steiert B, Steffens B, Pfister M, Koch G.",
    "Automated pharmacometric model development by leveraging",
    "low-dimensional neural ODEs and LASSO regression.",
    "CPT Pharmacometrics Syst Pharmacol. 2026.",
    "doi:10.1002/psp4.70285.",
    "Structure from Equation 20; parameter values from Table S2 of the",
    "Supporting Information (A' = 5.1, p = 0.46, B' = 5.0, k = 0.12)",
    "inverted through the paper's own Equation 21 definitions.",
    sep = " "
  )
  vignette <- "Bram_2026_node_lasso_model_development"
  units <- list(time = "h", dosing = "dose_unit", concentration = "conc_unit")

  compartmentData <- list(
    central = list(analyte = "simulated test compound", units = "dose_unit", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "None (simulated data; no drug and no patients)",
    n_subjects = 50,
    n_studies = 0,
    dose_range = "Single intravenous bolus; the dose is not reported and is fixed at 1 arbitrary unit here.",
    observations = "7 concentration measurements per subject over 24 h (Braem 2026 Section 2.3.2).",
    notes = "Braem 2026 Section 2.3.2: PK data were simulated for 50 subjects from the intravenous two-compartment model of Equation 14 with log-normally distributed individual parameters. The generating parameter values are not reported. This entry is a sensitivity analysis in the source paper, testing whether the NODE-LASSO workflow recovers a usable structure for data arising from a model with a latent state."
  )

  ini({
    # ---------------------------------------------------------------
    # Braem 2026 does not print the estimated parameters of Equation 20
    # directly. Table S2 reports the derived macro-parameters of the
    # explicit solution (Equation 21)
    #     Cc(t) = A' * exp(-kdist * t) + B' * exp(-kel * t)
    # as A' = 5.1, p = 0.46, B' = 5.0, k = 0.12, where the paper's `p`
    # is kdist and its `k` is kel. The remaining two parameters invert
    # exactly from the paper's own definitions
    #     A' = -(D/Vc) * wdist / (kel - kdist)
    #     B' =  (D/Vc) * (kel - kdist + wdist) / (kel - kdist)
    # whose sum telescopes to A' + B' = D/Vc = 10.1, giving
    #     wdist = A' * (kel - kdist) / -(D/Vc) = 0.1716832.
    # The inversion is re-derived and re-checked in the vignette.
    # ---------------------------------------------------------------
    lkel <- log(0.12); label("First-order elimination rate constant (1/h)")  # Table S2, row `k`
    lkdist <- log(0.46); label("Decay rate constant of the transient distribution flux (1/h)")  # Table S2, row `p`
    lwdist <- log(0.1716832); label("Amplitude of the transient distribution flux, per unit dose (1/h)")  # back-solved from Table S2 A' = 5.1, B' = 5.0
    lvc <- log(0.0990099); label("Central volume of distribution (dose_unit per conc_unit)")  # = 1 / (A' + B') = 1 / 10.1, at a nominal unit dose

    # Braem 2026 states that model parameters were assumed log-normally
    # distributed but reports no omega and no residual-error estimate
    # for Equation 20, so this entry is typical-value only. See the
    # vignette Errata.
  })

  model({
    kel <- exp(lkel)
    kdist <- exp(lkdist)
    wdist <- exp(lwdist)
    vc <- exp(lvc)

    # Braem 2026 Equation 20, written with podo() and tad() so the
    # transient flux is referenced to the most recent dose. For the
    # single bolus at t = 0 that the paper fitted, podo(central) is D
    # and tad(central) is t, so this is Equation 20 verbatim; for
    # multiple doses it is the time-after-dose form given in the
    # Supporting Information ("Multiple-dose scenario with the proposed
    # two-compartment model").
    #
    # The compartment argument is LOAD-BEARING. Through the rxUi path
    # that readModelDb() + rxSolve() use, the argument-less podo() and
    # tad() silently evaluate to 0, which drops this whole term and
    # degrades the model to a one-compartment mono-exponential without
    # any warning. Verified on rxode2 5.1.7: podo(central) / tad(central)
    # reproduce Equation 21 to 1.6e-06 relative, bare podo() / tad() do
    # not. The vignette asserts the closed form for exactly this reason.
    d/dt(central) <- -kel * central -
      podo(central) * wdist * exp(-kdist * tad(central))

    Cc <- central / vc
  })
}
