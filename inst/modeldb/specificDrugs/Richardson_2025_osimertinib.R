Richardson_2025_osimertinib <- function() {
  description <- paste(
    "Automated model-search (pyDarwin) one-compartment population PK model",
    "for oral osimertinib in 270 participants from Phase 1 study NCT01802632",
    "(Richardson 2025). This is the top-ranked structure returned by the",
    "paper's automated Bayesian-optimisation model search; it differs from the",
    "manually developed expert structure only by the addition of four transit",
    "compartments to the absorption model, which reduced the variability in",
    "the absorption rate constant. Absorption is first-order through a depot",
    "followed by four transit compartments, with every transfer governed by a",
    "single rate constant (the paper's KA = KTR). Disposition is",
    "one-compartment linear with apparent clearance and volume",
    "(bioavailability fixed to 1, so both are apparent CL/F and V/F). No dose",
    "effect was selected, so the model is dose-linear. IIV is diagonal",
    "log-normal on CL, Vc and the transit rate constant. Residual error is",
    "combined additive plus proportional. Only parent-drug concentrations were",
    "modelled; the metabolite data available in the source study were excluded",
    "because metabolite characterisation is not part of the paper's model",
    "space. Note that this is a machine-search structure, not the manually",
    "developed expert model.",
    sep = " "
  )
  reference <- paste(
    "Richardson S, Irurzun Arana I, Nowojewski A, Zhou D, Leander J,",
    "Tang W, Dearden R, Gibbs M.",
    "A machine learning approach to population pharmacokinetic modelling",
    "automation. Commun Med. 2025;5:325. doi:10.1038/s43856-025-01054-8.",
    "Parameter estimates from Supplementary Table 4; model structure from",
    "Table 4, Supplementary Data 1 and the pyDarwin model-space files",
    "(template.txt, tokens.json) in the Code Availability repository",
    "https://github.com/samjrrr/autopk_synthetic_example.",
    sep = " "
  )
  vignette <- "Richardson_2025_automated_popPK"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot    = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "osimertinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariates were retained in this structure -- the search selected no
  # dose effect for osimertinib (Table 4, Supplementary Data 1), and the model
  # space contains no demographic covariates at all.
  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 270,
    n_studies     = 1,
    disease_state = "not reported (Phase 1 study NCT01802632 in EGFR-mutation-positive advanced non-small-cell lung cancer)",
    dose_range    = "not reported",
    notes         = paste(
      "270 participants after filtering to oral administration and parent-drug",
      "measurements only, contributing 3,766 observation records (Richardson",
      "2025 Table 1 and Supplementary Table 7). Observed concentrations were",
      "modelled on the normal scale (Table 1 'Scale type'). The paper reports",
      "no baseline demographics for this cohort.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- apparent (F fixed to 1). Supplementary Table 4.
    lcl  <- log(13.84)   ; label("Apparent clearance CL/F (L/h)")                        # Suppl Table 4 THETA1 = 13.84 L/h
    lvc  <- log(1101.20) ; label("Apparent central volume V2/F (L)")                     # Suppl Table 4 THETA2 = 1101.20 L
    lktr <- log(1.70)    ; label("Transit / absorption rate constant KA = KTR (1/h)")    # Suppl Table 4 THETA3 = 1.70 1/h

    # Bioavailability was fixed to 1 in the model space (template.txt: F1 = 1),
    # so both disposition parameters above are apparent.
    lfdepot <- fixed(log(1)) ; label("Bioavailability F (unitless)")         # template.txt F1 = 1   ; INPUT1

    # IIV -- diagonal log-normal variances. Supplementary Table 4 reports these
    # as CV%; omega^2 = log(1 + CV^2), the inverse of the
    # CV% = sqrt(exp(omega^2) - 1) * 100 convention verified against the
    # authors' published NONMEM output for the ribociclib model.
    etalcl  ~ 0.247563 # Suppl Table 4 OMEGA(1,1) = 53% CV
    etalvc  ~ 0.239332 # Suppl Table 4 OMEGA(2,2) = 52% CV
    etalktr ~ 0.207339 # Suppl Table 4 OMEGA(3,3) = 48% CV

    # Residual error -- combined additive plus proportional on the linear
    # concentration scale. $SIGMA was 1 FIX and the error was carried in
    # THETA space: W = SQRT(add^2 + (prop * IPRED)^2); Y = IPRED + W * ERR(1).
    addSd  <- 2.61 ; label("Additive residual error (ng/mL)")            # Suppl Table 4 THETA4 = 2.61
    propSd <- 0.17 ; label("Proportional residual error (fraction)")    # Suppl Table 4 THETA5 = 0.17
  })

  model({
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    ktr <- exp(lktr + etalktr)

    kel <- cl / vc

    # Absorption chain: depot -> transit1 -> ... -> transit4 -> central, every
    # transfer governed by the single rate constant ktr (tokens.json
    # TRANSIT_COMP option 3 sets RATE_CONSTANT(...) = KA for all five steps).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(central)  <-  ktr * transit4 - kel * central

    f(depot) <- exp(lfdepot)

    # central is in mg and vc in L, so central / vc is mg/L; the factor of
    # 1000 puts the observation on the ng/mL scale that the additive residual
    # error (2.61) is expressed in.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
