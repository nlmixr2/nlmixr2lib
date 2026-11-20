Jamsen_2011_dihydroartemisinin_dost <- function() {
  description <- paste(
    "One-compartment population PK model in which the absorption and",
    "elimination rate constants are equal ('Dost model with a lag-time',",
    "C = D/V * k * (t - tlag) * exp(-k * (t - tlag))) for dihydroartemisinin",
    "(DHA) after oral artesunate in 32 non-pregnant adult patients with",
    "uncomplicated Plasmodium falciparum malaria on the Thai-Burmese border",
    "and in western Cambodia (Jamsen 2011). One of four competing structural",
    "models used to build robust optimal sampling designs; the best-fitting",
    "model is Jamsen_2011_dihydroartemisinin. Artesunate is assumed to be",
    "completely converted to DHA and the dose is converted to DHA mass by",
    "molecular weight. No covariates."
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
    lkel <- log(0.99); label("Dost rate constant k, shared by absorption and elimination (1/h)") # Table 1, Dost model, non-pregnant adults, lag-time included: k = 0.99 /h (95% CI 0.91-1.1)
    lvc <- log(46.3); label("Apparent volume of distribution V/F of DHA (L)") # Table 1, Dost, lag included: V/F = 46.3 L (37.6-55.0)
    ltlag <- fixed(log(0.21)); label("Absorption lag time tlag (h)") # Table 1, Dost, lag included: tlag = 0.21 h (fixed)

    # Table 1 footnote: variances printed as %CV = 100*sqrt(variance),
    # covariances as correlation coefficients.
    etalkel + etalvc ~ c(0.0484, -0.0837716, 0.232324) # Table 1: Omega k 22.0 %CV, Omega V/F 48.2 %CV, correlation -0.79 -> -0.79 * 0.220 * 0.482
    etaltlag ~ 0.012996 # Table 1: Omega tlag 11.4 %CV

    propSd <- 0.47; label("Proportional residual error (fraction)") # Table 1: sigma = 0.47 (residual SD, proportional)
  })
  model({
    # Artesunate (MW 384.42 g/mol) -> DHA (MW 284.35 g/mol), complete
    # conversion assumed; f(depot) is the molar conversion, not a
    # bioavailability.
    mw_ratio <- 284.35 / 384.42

    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag + etaltlag)

    # Dost (Additional file 1): one rate constant k for both absorption and
    # elimination. A depot -> central chain with both rates equal to k
    # integrates to D/V * k * (t - tlag) * exp(-k * (t - tlag)).
    ka <- kel
    cl <- kel * vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- mw_ratio
    alag(depot) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
