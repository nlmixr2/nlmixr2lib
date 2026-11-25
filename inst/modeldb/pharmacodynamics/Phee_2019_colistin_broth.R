Phee_2019_colistin_broth <- function() {
  description <- "In vitro (Iso-Sensitest broth, no bacteria). First-order loss (degradation and/or binding) of colistin in broth over the 24-h span of a static time-kill experiment, estimated from LC-MS/MS concentrations of colistin alone and with fusidic acid. Fusidic acid concentrations stayed constant and carry no loss term. The state holds the broth concentration itself: the dose is the nominal starting concentration in mg/L. Companion to Phee_2019_colistin_fusidicAcid, whose final CFU model does not use this loss."
  reference <- "Phee LM, Kloprogge F, Morris R, Barrett J, Wareham DW, Standing JF. Pharmacokinetic-pharmacodynamic modelling to investigate in vitro synergy between colistin and fusidic acid against MDR Acinetobacter baumannii. J Antimicrob Chemother. 2019;74(4):961. doi:10.1093/jac/dky524. Structural equations: Methods Eqs 1-2. Parameter estimate: Results, Drug degradation experiments."
  vignette <- "Phee_2019_colistin_fusidicAcid"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "colistin in Iso-Sensitest broth (colistin A + colistin B by LC-MS/MS)",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (cell-free Iso-Sensitest broth)",
    n_subjects = 103L,
    n_studies = 1L,
    disease_state = "Not applicable: drug stability in bacteria-free broth incubated in parallel with the time-kill experiments.",
    model_system = "Simulated samples prepared in parallel to the time-kill assays; 103 aliquots with colistin 0-4 mg/L and fusidic acid 0-64 mg/L, singly or in combination, sampled at 0, 2, 6 and 24 h and assayed by LC-MS/MS after 0.22 um filtration and protein precipitation.",
    dose_range = "Colistin 0-4 mg/L (nominal starting concentration)",
    notes = paste(
      "Fitted in NONMEM 7.3 with ADVAN9 and FOCE-I on natural-log concentrations.",
      "The starting concentration was fixed to the observed baseline with a residual term of the same size as the residual error (Eq 1); the paper does not report that residual error.",
      "Fusidic acid showed no loss and has no parameter; colistin loss was not carried into the final CFU model because it did not improve the fit (P > 0.05)."
    )
  )

  ini({
    lkdeg <- log(0.0175); label("First-order loss rate constant of colistin in broth, theta_k (1/h)") # Results, Drug degradation experiments: 0.0175 1/h (RSE 20.9%)
    # The paper prints no residual error for the concentration model; held at
    # zero so the model simulates the typical profile only.
    propSd <- fixed(0); label("Proportional residual error (fraction), not reported") # not reported (Methods Eq 1 names a proportional error but gives no value)
  })

  model({
    kdeg <- exp(lkdeg)

    # Eq 2: dC/dt = -theta_k * C, the state being the broth concentration.
    d/dt(central) <- -kdeg * central

    Cc <- central
    Cc ~ prop(propSd)
  })
}
