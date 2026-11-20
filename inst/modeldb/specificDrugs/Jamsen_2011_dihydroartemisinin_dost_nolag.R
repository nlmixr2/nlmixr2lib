Jamsen_2011_dihydroartemisinin_dost_nolag <- function() {
  description <- paste(
    "One-compartment population PK model in which the absorption and",
    "elimination rate constants are equal and there is no lag time ('Dost",
    "model without a lag-time', C = D/V * k * t * exp(-k * t)) for",
    "dihydroartemisinin (DHA) after oral artesunate in 32 non-pregnant adult",
    "patients with uncomplicated Plasmodium falciparum malaria on the",
    "Thai-Burmese border and in western Cambodia (Jamsen 2011). One of four",
    "competing structural models used to build robust optimal sampling",
    "designs; the best-fitting model is Jamsen_2011_dihydroartemisinin.",
    "Artesunate is assumed to be completely converted to DHA and the dose is",
    "converted to DHA mass by molecular weight. No covariates."
  )
  reference <- paste(
    "Jamsen KM, Duffull SB, Tarning J, Lindegardh N, White NJ, Simpson JA.",
    "Optimal designs for population pharmacokinetic studies of oral",
    "artesunate in patients with uncomplicated falciparum malaria.",
    "Malar J. 2011;10:181. doi:10.1186/1475-2875-10-181.",
    sep = " "
  )
  vignette <- "Jamsen_2011_dihydroartemisinin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "artesunate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "adults (not reported in this paper)",
    weight_median = "51 kg",
    disease_state = "Uncomplicated Plasmodium falciparum malaria; not pregnant",
    dose_range = "4 mg/kg oral artesunate, single dose; 6-12 venous samples over the first 12 h",
    regions = "Thai-Burmese border (Wang Pha, Thailand) and western Cambodia (Pailin)",
    notes = paste(
      "Intensively sampled subset of the Dondorp 2009 (N Engl J Med",
      "361:455-467) artemisinin-resistance study, re-analysed by Jamsen 2011.",
      "The median weight of 51 kg is from Additional file 1."
    )
  )

  ini({
    lkel <- log(0.95); label("Dost rate constant k, shared by absorption and elimination (1/h)") # Table 1, Dost model, non-pregnant adults, lag-time not included: k = 0.95 /h (95% CI 0.86-1.03)
    lvc <- log(51.3); label("Apparent volume of distribution V/F of DHA (L)") # Table 1, Dost, no lag: V/F = 51.3 L (41.2-61.4)

    # Table 1 footnote: variances printed as %CV = 100*sqrt(variance),
    # covariances as correlation coefficients.
    etalkel + etalvc ~ c(0.043264, -0.07109856, 0.178084) # Table 1: Omega k 20.8 %CV, Omega V/F 42.2 %CV, correlation -0.81 -> -0.81 * 0.208 * 0.422

    propSd <- 0.57; label("Proportional residual error (fraction)") # Table 1: sigma = 0.57 (residual SD, proportional)
  })
  model({
    # Artesunate (MW 384.42 g/mol) -> DHA (MW 284.35 g/mol), complete
    # conversion assumed; f(depot) is the molar conversion, not a
    # bioavailability.
    mw_ratio <- 284.35 / 384.42

    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)

    # Dost (Additional file 1): one rate constant k for both absorption and
    # elimination.
    ka <- kel
    cl <- kel * vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- mw_ratio

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
