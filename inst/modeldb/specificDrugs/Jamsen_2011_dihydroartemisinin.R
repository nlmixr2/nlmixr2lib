Jamsen_2011_dihydroartemisinin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time ('Bateman model with a lag-time') for",
    "dihydroartemisinin (DHA) after oral artesunate in 32 non-pregnant adult",
    "patients with uncomplicated Plasmodium falciparum malaria on the",
    "Thai-Burmese border and in western Cambodia (Jamsen 2011). Artesunate is",
    "assumed to be completely converted to DHA, so only DHA is modelled; the",
    "artesunate dose is converted to DHA mass by molecular weight. This is the",
    "best-fitting of four competing structural models and the one the authors",
    "simulated from; the other three are Jamsen_2011_dihydroartemisinin_nolag,",
    "Jamsen_2011_dihydroartemisinin_dost and",
    "Jamsen_2011_dihydroartemisinin_dost_nolag. No covariates."
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
      "361:455-467) artemisinin-resistance study, re-analysed by Jamsen 2011",
      "(Methods, 'Determining pharmacokinetic models'). The median weight of",
      "51 kg is from Additional file 1. Plasma DHA was measured by LC-MS/MS.",
      "Age, sex and other demographics are given in Dondorp 2009, not here."
    )
  )

  ini({
    lka <- log(0.82); label("Absorption rate constant ka (1/h)") # Table 1, Bateman model, non-pregnant adults, lag-time included: ka = 0.82 /h (95% CI 0.76-0.87)
    lcl <- log(47.5); label("Apparent oral clearance CL/F of DHA (L/h)") # Table 1, Bateman, lag included: CL/F = 47.5 L/h (44.6-50.4)
    lvc <- log(32.1); label("Apparent volume of distribution V/F of DHA (L)") # Table 1, Bateman, lag included: V/F = 32.1 L (24.6-39.6)
    ltlag <- fixed(log(0.21)); label("Absorption lag time tlag (h)") # Table 1, Bateman, lag included: tlag = 0.21 h (fixed)

    # Table 1 footnote: between-subject variances are printed as %CV =
    # 100*sqrt(variance) and covariances as correlation coefficients, so
    # variance = (CV/100)^2 and covariance = r * sd1 * sd2.
    etalcl + etalvc ~ c(0.086436, 0.139656, 0.670761) # Table 1: Omega CL/F 29.4 %CV, Omega V/F 81.9 %CV, correlation 0.58 -> 0.58 * 0.294 * 0.819
    etalka + etaltlag ~ c(0.071289, 0.005073, 0.01) # Table 1: Omega ka 26.7 %CV, Omega tlag 10.0 %CV, correlation 0.19 -> 0.19 * 0.267 * 0.100

    propSd <- 0.41; label("Proportional residual error (fraction)") # Table 1: sigma = 0.41 (residual SD, proportional)
  })
  model({
    # Artesunate (MW 384.42 g/mol) is assumed completely converted to DHA
    # (MW 284.35 g/mol; Methods). f(depot) is that molar conversion, not a
    # bioavailability: CL/F and V/F already carry F.
    mw_ratio <- 284.35 / 384.42

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- mw_ratio
    alag(depot) <- tlag

    # central in mg, vc in L -> mg/L; x1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
