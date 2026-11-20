Jamsen_2011_dihydroartemisinin_children <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for dihydroartemisinin (DHA) after oral artesunate in",
    "infants, children and non-pregnant adults with uncomplicated Plasmodium",
    "falciparum malaria: the adults-and-children design and simulation model",
    "of Jamsen 2011. It is the non-pregnant adult fit",
    "(Jamsen_2011_dihydroartemisinin) extrapolated by allometric scaling of",
    "CL/F (exponent 0.75) and V/F (exponent 1) to a 51 kg reference weight,",
    "with an exponential maturation of CL/F from 1/3.1 of the size-scaled",
    "adult value at birth (first-order rate 0.082 per month of age, the",
    "glomerular-filtration-rate maturation rate of Hayton 2002). The",
    "paediatric parameters are an a-priori extrapolation, not a fit to",
    "paediatric data. Artesunate is assumed to be completely converted to DHA",
    "and the dose is converted to DHA mass by molecular weight."
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

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL/F with exponent 0.75 and of V/F with",
        "exponent 1, referenced to 51 kg, the median weight of the",
        "non-pregnant adults (Additional file 1)."
      ),
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the maturation of CL/F. Additional file 1 writes the",
        "maturation in age in months; the model converts AGE (years) to",
        "months (x 12) before applying the 0.082 per month rate. Maturation",
        "is complete (>99.99%) by about 10 years, so AGE has no effect in",
        "adults."
      ),
      source_name = "age"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "artesunate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    age_range = "infants (<2 years) to adults (>20 years)",
    weight_median = "51 kg (adult reference)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria; not pregnant",
    dose_range = "4 mg/kg oral artesunate",
    notes = paste(
      "Hypothetical design population of Jamsen 2011: 60 patients in four",
      "age groups (<2 years n = 10, 2-10 years n = 10, 11-20 years n = 10,",
      ">20 years n = 30; Table 2), with age/weight combinations drawn from a",
      "linear model built on data from Simpson 2006 (PLoS Med 3:e444). The",
      "adult parameters come from 32 non-pregnant adults (Dondorp 2009",
      "subset); no paediatric PK data entered the model."
    )
  )

  ini({
    lka <- log(0.82); label("Absorption rate constant ka (1/h)") # Table 1, Bateman model, non-pregnant adults, lag-time included: ka = 0.82 /h
    lcl <- log(47.5); label("Apparent oral clearance CL/F of DHA for a mature 51 kg patient (L/h)") # Table 1, Bateman, lag included: CL/F = 47.5 L/h; Additional file 1: CL/Fpop set to the Table 1 value
    lvc <- log(32.1); label("Apparent volume of distribution V/F of DHA for a 51 kg patient (L)") # Table 1, Bateman, lag included: V/F = 32.1 L
    ltlag <- fixed(log(0.21)); label("Absorption lag time tlag (h)") # Table 1, Bateman, lag included: tlag = 0.21 h (fixed)

    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)") # Additional file 1, CL/F equation: (weight_i/weight_pop)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on V/F (unitless)") # Additional file 1, V/F equation: V/Fpop x (weight_i/weight_pop)

    # Additional file 1: CL/F_i = (CL/Fpop / 3.1) * (WT/51)^0.75 * exp(-kmat * age)
    #   + CL/Fpop * (WT/51)^0.75 * (1 - exp(-kmat * age)), age in months.
    # Rearranged: CL/F_i = CL/Fpop * (WT/51)^0.75 * (1 + famp * exp(-kmat * age))
    # with famp = 1/3.1 - 1, the fractional shortfall of CL/F at birth.
    cl_pna_famp <- fixed(1 / 3.1 - 1); label("Fractional difference of CL/F from its size-scaled mature value at birth (fraction)") # Additional file 1: newborn CL/F = CL/Fpop / 3.1 -> famp = 1/3.1 - 1
    lcl_pna_kdes <- fixed(log(0.082)); label("First-order maturation rate of CL/F, kmat (1/month of age)") # Additional file 1: 'kmat was set to 0.082, which assumes CL/F maturates at the same rate that the glomerular filtration rate (GFR) matures'

    # Table 1 footnote: variances printed as %CV = 100*sqrt(variance),
    # covariances as correlation coefficients. The simulation for this design
    # used the non-pregnant adult population estimates (Results).
    etalcl + etalvc ~ c(0.086436, 0.139656, 0.670761) # Table 1: Omega CL/F 29.4 %CV, Omega V/F 81.9 %CV, correlation 0.58
    etalka + etaltlag ~ c(0.071289, 0.005073, 0.01) # Table 1: Omega ka 26.7 %CV, Omega tlag 10.0 %CV, correlation 0.19

    propSd <- 0.41; label("Proportional residual error (fraction)") # Table 1: sigma = 0.41 (residual SD, proportional)
  })
  model({
    # Artesunate (MW 384.42 g/mol) -> DHA (MW 284.35 g/mol), complete
    # conversion assumed; f(depot) is the molar conversion, not a
    # bioavailability.
    mw_ratio <- 284.35 / 384.42

    # Additional file 1 states the maturation in age in months.
    age_months <- AGE * 12
    cl_pna_kdes <- exp(lcl_pna_kdes)
    fmat_cl <- 1 + cl_pna_famp * exp(-cl_pna_kdes * age_months)

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 51)^e_wt_cl * fmat_cl
    vc <- exp(lvc + etalvc) * (WT / 51)^e_wt_vc
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- mw_ratio
    alag(depot) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
