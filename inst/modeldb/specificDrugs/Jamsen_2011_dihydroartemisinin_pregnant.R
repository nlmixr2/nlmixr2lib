Jamsen_2011_dihydroartemisinin_pregnant <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for dihydroartemisinin (DHA) after oral artesunate in",
    "pregnant women with uncomplicated Plasmodium falciparum malaria: the",
    "pregnant-women design and simulation model of Jamsen 2011. Its",
    "structural values are those of McGready 2006 (Eur J Clin Pharmacol",
    "62:367-371), re-expressed for a 50 kg woman; between-subject variability",
    "of ka and the residual error, which McGready 2006 did not report, are",
    "taken from the Jamsen 2011 non-pregnant adult fit, and the BSV of V/F is",
    "reduced to a variance of 1, as Jamsen 2011 did for simulation.",
    "Artesunate is assumed to be completely converted to DHA and the dose is",
    "converted to DHA mass by molecular weight. No covariates."
  )
  reference <- paste(
    "Jamsen KM, Duffull SB, Tarning J, Lindegardh N, White NJ, Simpson JA.",
    "Optimal designs for population pharmacokinetic studies of oral",
    "artesunate in patients with uncomplicated falciparum malaria.",
    "Malar J. 2011;10:181. doi:10.1186/1475-2875-10-181.",
    "Structural values from: McGready R, Stepniewska K, Ward SA, Cho T,",
    "Gilveray G, Looareesuwan S, White NJ, Nosten F. Pharmacokinetics of",
    "dihydroartemisinin following oral artesunate treatment of pregnant women",
    "with acute uncomplicated falciparum malaria. Eur J Clin Pharmacol.",
    "2006;62:367-371. doi:10.1007/s00228-006-0118-y.",
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
    weight_median = "50 kg",
    sex_female_pct = 100,
    disease_state = "Uncomplicated Plasmodium falciparum malaria in pregnancy",
    dose_range = "4 mg/kg oral artesunate",
    regions = "Thai-Burmese border",
    notes = paste(
      "Design population of Jamsen 2011 for pregnant women. The values come",
      "from McGready 2006, 'the only published population PK study of oral",
      "AS in pregnant women with presenting falciparum malaria' (Jamsen 2011",
      "Methods); the number of women and their demographics are reported",
      "there, not here. McGready 2006 reported CL/F in L/kg/h and V/F in",
      "L/kg; Jamsen 2011 Table 1 multiplied them by the median weight of",
      "50 kg."
    )
  )

  ini({
    lka <- log(1.19); label("Absorption rate constant ka (1/h)") # Table 1, Bateman model, pregnant women (from McGready 2006): ka = 1.19 /h (95% CI 0.78-1.60)
    lcl <- log(88.5); label("Apparent oral clearance CL/F of DHA for a 50 kg woman (L/h)") # Table 1, pregnant women: CL/F = 88.5 L/h (60-117); footnote: L/kg/h value x median weight 50 kg
    lvc <- log(232); label("Apparent volume of distribution V/F of DHA for a 50 kg woman (L)") # Table 1, pregnant women: V/F = 232 L (57.0-406); footnote: L/kg value x 50 kg
    ltlag <- log(0.42); label("Absorption lag time tlag (h)") # Table 1, pregnant women: tlag = 0.42 h (0.34-0.50)

    # Table 1 prints the BSV as %CV = 100*sqrt(variance).
    etalka ~ 0.071289 # Results: 'the values of these parameters were set to those of the non-pregnant adults' -> Table 1 non-pregnant Bateman-with-lag Omega ka 26.7 %CV
    etalcl ~ 0.2209 # Table 1, pregnant women: Omega CL/F 47.0 %CV
    etalvc ~ 1 # Table 2 footnote and Results: 'the BSV of V/F was reduced to 1' for design and simulation (McGready 2006 reported 154 %CV, variance 2.37)

    propSd <- 0.41; label("Proportional residual error (fraction)") # Results: sigma not reported by McGready 2006, set to the non-pregnant adult value -> Table 1 Bateman-with-lag sigma = 0.41
  })
  model({
    # Artesunate (MW 384.42 g/mol) -> DHA (MW 284.35 g/mol), complete
    # conversion assumed; f(depot) is the molar conversion, not a
    # bioavailability.
    mw_ratio <- 284.35 / 384.42

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- mw_ratio
    alag(depot) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
