Bram_2026_biexponential_2cmt <- function() {
  description <- "Two-compartment intravenous-bolus reference model of Braem 2026 (Equation 14), fitted to the same simulated bi-exponential dataset as the automatically-proposed one-state model Bram_2026_biexponential and reported alongside it as the conventional-pharmacometrics benchmark (MARE(PMX) in Table 1). The data were SIMULATED for 50 subjects with 7 concentration measurements over 24 h; no drug, no dose and no units are reported. Braem 2026 prints only the macro-parameters of the explicit solution (Table S2: A = 5.1, alpha = 0.46, B = 5.0, beta = 0.12), so the micro-constants here are the standard macro-to-micro inversion of those four values, and the fit identifies the model only up to D/Vc = 10.1. A nominal unit dose is therefore fixed (vc = 1/10.1) so that a dose of 1 arbitrary unit reproduces the published macro-parameters exactly. The point of the pair is that the two models have IDENTICAL explicit solutions for the central compartment despite one of them having no peripheral state."
  reference <- paste(
    "Braem DS, Steiert B, Steffens B, Pfister M, Koch G.",
    "Automated pharmacometric model development by leveraging",
    "low-dimensional neural ODEs and LASSO regression.",
    "CPT Pharmacometrics Syst Pharmacol. 2026.",
    "doi:10.1002/psp4.70285.",
    "Structure from Equation 14; parameter values are the macro-to-micro",
    "inversion of Table S2 of the Supporting Information",
    "(A = 5.1, alpha = 0.46, B = 5.0, beta = 0.12).",
    sep = " "
  )
  vignette <- "Bram_2026_node_lasso_model_development"
  units <- list(time = "h", dosing = "dose_unit", concentration = "conc_unit")

  compartmentData <- list(
    central = list(analyte = "simulated test compound", units = "dose_unit", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "simulated test compound", units = "dose_unit", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "None (simulated data; no drug and no patients)",
    n_subjects = 50,
    n_studies = 0,
    dose_range = "Single intravenous bolus; the dose is not reported and is fixed at 1 arbitrary unit here.",
    observations = "7 concentration measurements per subject over 24 h (Braem 2026 Section 2.3.2).",
    notes = "Braem 2026 Section 2.3.2 and Section 3.2.2: this is the model that generated the data and that was also re-fitted to it as the conventional-pharmacometrics reference."
  )

  ini({
    # ---------------------------------------------------------------
    # Table S2 reports the explicit-solution macro-parameters of
    # Equation 22, Cc(t) = A * exp(-alpha * t) + B * exp(-beta * t),
    # as A = 5.1, alpha = 0.46, B = 5.0, beta = 0.12. The standard
    # inversion for an intravenous bolus into a two-compartment system,
    # at a nominal unit dose D = 1, is
    #     Vc  = D / (A + B)                    = 1 / 10.1 = 0.0990099
    #     k21 = (A*beta + B*alpha) / (A + B)   = 0.2883168 /h
    #     kel = alpha * beta / k21             = 0.1914560 /h
    #     k12 = alpha + beta - kel - k21       = 0.1002271 /h
    # and then CL = kel*Vc, Q = k12*Vc, Vp = Q/k21. Storing CL/Vc/Q/Vp
    # rather than the micro-constants keeps the model solvable through
    # the standard two-compartment machinery. Both invariants check
    # exactly: kel + k12 + k21 = alpha + beta = 0.58 and
    # kel * k21 = alpha * beta = 0.0552. Re-derived in the vignette.
    # ---------------------------------------------------------------
    lcl <- log(0.01895604); label("Clearance (dose_unit per conc_unit per h)")  # = kel * vc, from Table S2
    lvc <- log(0.0990099); label("Central volume of distribution (dose_unit per conc_unit)")  # = 1 / (A + B), at a nominal unit dose
    lq <- log(0.009923478); label("Inter-compartmental clearance (dose_unit per conc_unit per h)")  # = k12 * vc, from Table S2
    lvp <- log(0.03441866); label("Peripheral volume of distribution (dose_unit per conc_unit)")  # = q / k21, from Table S2

    # Braem 2026 reports no omega and no residual-error estimate for
    # Equation 14, so this entry is typical-value only. See the
    # vignette Errata.
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Braem 2026 Equation 14.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
  })
}
