Chen_2022_creatinine <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model for endogenous creatinine in healthy volunteers",
    "with and without ingestion of a cooked-beef meal (Chen 2022), a reanalysis of the Mayersohn 1983",
    "plasma and 24 h urine data. Creatinine is generated at a zero-order rate (CGR), eliminated by",
    "first-order renal clearance into a cumulative urine state, and the meal creatinine load is",
    "absorbed first-order after a lag time. The meal dose is a nominal 180 mg population dose scaled",
    "by an estimated dose-correction factor F1. Plasma starts at the endogenous steady state CGR/CL.",
    sep = " "
  )
  reference <- paste(
    "Chen Z, Chen C, Taubert M, Mayersohn M, Fuhr U.",
    "A population pharmacokinetic model for creatinine with and without ingestion of a cooked meat meal.",
    "Eur J Clin Pharmacol. 2022;78(12):1945-1947.",
    "doi:10.1007/s00228-022-03398-9.",
    "Data from Mayersohn M, Conrad KA, Achari R. Br J Clin Pharmacol. 1983;15(2):227-230.",
    "doi:10.1111/j.1365-2125.1983.tb01490.x.",
    sep = " "
  )
  vignette <- "Chen_2022_creatinine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "creatinine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "creatinine", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "creatinine", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    n_observations = 144L,
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Crossover: a control period without a meal, and a period with ingestion of a cooked (boiled)",
      "beef meal as an exogenous creatinine source. The meal creatinine content was not measured; the",
      "model uses a nominal 180 mg population dose scaled by the estimated F1 (individual estimated",
      "bioavailable dose 293 +/- 62.0 mg, n = 6)."
    ),
    notes = paste(
      "Reanalysis of the Mayersohn 1983 dataset: 133 serial creatinine plasma concentrations and 11",
      "creatinine amounts excreted in urine over 24 h. The letter states that no demographic information",
      "is available in the original publication; the 'male, 31 years, 73 kg' reference subject and the",
      "73 kg implied by 'Vd = 53.9 L (73.8% of the total body weight)' are the only demographic values",
      "it gives, and the authors use 73 kg as the body weight in their total-body-water comparison.",
      "n_subjects = 6 is inferred from the individual dose estimates (n = 6); the 11 urine",
      "amounts are consistent with 6 control-period plus 5 beef-period collections ('n = 5').",
      "Fit with exponential IIV and separate proportional residual errors for plasma and urine;",
      "994 of 1,000 bootstrap runs succeeded."
    )
  )

  ini({
    lka <- log(1.76); label("Absorption rate constant of meal creatinine Ka (1/h)") # Table 1 'Ka (1/h)' = 1.76 (RSE 26.7%; bootstrap median 1.60, 95% CI 0.715-3.34)
    lcl <- log(7.59); label("Creatinine renal clearance CL (L/h)") # Table 1 'CL (L/h)' = 7.59 (RSE 6.4%; bootstrap median 7.59, 95% CI 6.71-8.48). The Results prose quotes 7.92 L/h; see vignette Errata.
    lvc <- log(53.9); label("Creatinine volume of distribution Vd (L)") # Table 1 'V d (L)' = 53.9 (RSE 20.0%; bootstrap median 53.0, 95% CI 31.2-72.0)
    lfdepot <- log(1.57); label("Dose correction factor F1 on the nominal 180 mg meal dose (unitless)") # Table 1 'F1' = 1.57 (RSE 21.1%; bootstrap median 1.60, 95% CI 1.04-2.07)
    lksyn <- log(67.9); label("Creatinine generation rate CGR (mg/h)") # Table 1 'CGR (mg/h)' = 67.9 (RSE 6.9%; bootstrap median 67.9, 95% CI 59.4-76.3)
    ltlag <- log(0.344); label("Absorption lag time of meal creatinine (h)") # Table 1 'Lag time (h)' = 0.344 (RSE 4.7%; bootstrap median 0.344, 95% CI 0.297-0.370)

    # Exponential IIV; Table 1 prints CV% only, converted as omega^2 = log(1 + CV^2).
    # CL and Vd carry no IIV ('Incorporating variability on CL and Vd did not
    # significantly improve the model').
    etalka ~ 0.26683 # Table 1 IIV 'Ka (CV%)' = 55.3 (shrinkage 1.6%); log(1 + 0.553^2) = 0.26683
    etalfdepot ~ 0.058755 # Table 1 IIV 'F1 (CV%)' = 24.6 (shrinkage 0.1%); log(1 + 0.246^2) = 0.058755
    etalksyn ~ 0.0015199 # Table 1 IIV 'CGR (CV%)' = 3.9 (shrinkage 0.1%); log(1 + 0.039^2) = 0.0015199

    propSd <- 0.039; label("Proportional residual error, creatinine plasma concentration (fraction)") # Table 1 'Proportional error, % (plasma)' = 3.9 (shrinkage 7.0%)
    propSd_Aurine <- 0.155; label("Proportional residual error, creatinine amount excreted in urine (fraction)") # Table 1 'Proportional error, % (urine)' = 15.5 (shrinkage 0.2%)
  })

  model({
    # Individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ksyn <- exp(lksyn + etalksyn)
    tlag <- exp(ltlag)

    kel <- cl / vc

    # Zero-order endogenous generation into plasma, first-order meal absorption,
    # and renal elimination accumulating in the urine state.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + ksyn - kel * central
    d/dt(urine) <- kel * central

    f(depot) <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # Plasma creatinine starts at its endogenous steady state, C0 = CGR / CL.
    central(0) <- ksyn / cl * vc

    Cc <- central / vc
    Aurine <- urine

    Cc ~ prop(propSd)
    Aurine ~ prop(propSd_Aurine)
  })
}
