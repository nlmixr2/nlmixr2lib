Richardson_2025_tezepelumab <- function() {
  description <- paste(
    "Automated model-search (pyDarwin) two-compartment population PK model",
    "for subcutaneous tezepelumab in 106 participants pooled from four Phase 1",
    "studies (Richardson 2025). This is the top-ranked structure returned by",
    "the paper's automated Bayesian-optimisation model search; it differs from",
    "the manually developed expert structure primarily in its absorption model",
    "- sequential zero-order input into the depot over a short modelled",
    "duration followed by first-order transfer into the central compartment,",
    "rather than a purely first-order process - and in using a proportional",
    "rather than a combined residual error. Disposition is two-compartment",
    "linear. Subcutaneous bioavailability was fixed to 1 because no",
    "intravenous data were included, so all disposition parameters are",
    "apparent CL/F and V/F. No dose effect was selected, so the model is",
    "dose-linear. IIV is diagonal log-normal and estimated on CL, Vc, Vp, the",
    "absorption rate constant and the zero-order duration, with IIV on Q fixed",
    "at 15% CV. Note that this is a machine-search structure, not the manually",
    "developed expert model.",
    sep = " "
  )
  reference <- paste(
    "Richardson S, Irurzun Arana I, Nowojewski A, Zhou D, Leander J,",
    "Tang W, Dearden R, Gibbs M.",
    "A machine learning approach to population pharmacokinetic modelling",
    "automation. Commun Med. 2025;5:325. doi:10.1038/s43856-025-01054-8.",
    "Parameter estimates from Supplementary Table 6; model structure from",
    "Table 4, Supplementary Data 1 and the pyDarwin model-space files",
    "(template.txt, tokens.json) in the Code Availability repository",
    "https://github.com/samjrrr/autopk_synthetic_example.",
    sep = " "
  )
  vignette <- "Richardson_2025_automated_popPK"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "tezepelumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tezepelumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tezepelumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  # No covariates were retained in this structure -- the search selected no
  # dose effect for tezepelumab (Table 4, Supplementary Data 1), and the model
  # space contains no demographic covariates at all.
  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 106,
    n_studies     = 4,
    disease_state = "not reported (Phase 1 studies; tezepelumab is developed for severe asthma)",
    dose_range    = "not reported",
    notes         = paste(
      "106 participants after filtering to subcutaneous administration only,",
      "contributing 1,477 observation records, pooled from NCT00757042,",
      "NCT00972179, NCT02512900 and NCT01913028 (Richardson 2025 Table 1 and",
      "Supplementary Table 7). The intravenous data available in the source",
      "studies were excluded, which is why subcutaneous bioavailability could",
      "not be identified and was fixed to 1. Observed concentrations were",
      "modelled on the normal scale (Table 1 'Scale type'). The paper reports",
      "no baseline demographics for this cohort.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- apparent (F fixed to 1). Supplementary Table 6.
    lcl  <- log(0.21) ; label("Apparent clearance CL/F (L/day)")                              # Suppl Table 6 THETA1 = 0.21 L/day
    lvc  <- log(2.84) ; label("Apparent central volume V2/F (L)")                             # Suppl Table 6 THETA2 = 2.84 L
    lvp  <- log(3.6)  ; label("Apparent peripheral volume V3/F (L)")                          # Suppl Table 6 THETA6 = 3.6 L
    lq   <- log(6.04) ; label("Apparent inter-compartmental clearance Q/F (L/day)")           # Suppl Table 6 THETA7 = 6.04 L/day
    lka  <- log(0.44) ; label("First-order absorption rate constant KA (1/day)")              # Suppl Table 6 THETA3 = 0.44 /day
    ld1  <- log(0.11) ; label("Duration of zero-order input into the depot D1 (day)")         # Suppl Table 6 THETA8 = 0.11 days

    # Subcutaneous bioavailability was fixed to 1 because no IV data were
    # considered, so every disposition parameter above is apparent.
    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (unitless)") # Suppl Table 6 footnote: The subcutaneous bioavailability was fixed to 1 as no IV data was considered

    # IIV -- diagonal log-normal variances. Supplementary Table 6 reports these
    # as CV%; omega^2 = log(1 + CV^2), the inverse of the
    # CV% = sqrt(exp(omega^2) - 1) * 100 convention verified against the
    # authors' published NONMEM output for the ribociclib model.
    etalcl ~ 0.141586        # Suppl Table 6 OMEGA(1,1) = 39% CV
    etalvc ~ 0.389403        # Suppl Table 6 OMEGA(2,2) = 69% CV
    etalka ~ 0.272771        # Suppl Table 6 OMEGA(3,3) = 56% CV
    etalvp ~ 0.109392        # Suppl Table 6 OMEGA(4,4) = 34% CV
    etalq  ~ fixed(0.022251) # Suppl Table 6 OMEGA(5,5) = 15% (fix); tokens.json OMEGAS_Q 0.0225 FIX
    etald1 ~ 0.514260        # Suppl Table 6 OMEGA(6,6) = 82% CV

    # Residual error -- proportional only. The additive term was fixed to zero
    # (tokens.json ERRORFIX option 1: 0 FIX ; add_error), so the combined
    # form W = SQRT(add^2 + (prop * IPRED)^2) collapses to prop * IPRED.
    addSd  <- fixed(0) ; label("Additive residual error (ug/mL)")  # Suppl Table 6 THETA4 = 0 (fix)
    propSd <- 0.1      ; label("Proportional residual error (fraction)")      # Suppl Table 6 THETA5 = 0.1
  })

  model({
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq  + etalq)
    ka  <- exp(lka + etalka)
    d1  <- exp(ld1 + etald1)

    # Micro-constants, matching the ADVAN5 parameterisation of the source
    # control stream (k20 = CL/V2, k23 = Q/V2, k32 = Q/V3, k12 = KA).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Sequential absorption: the subcutaneous dose enters the depot as a
    # zero-order input over the modelled duration d1, then transfers to the
    # central compartment first-order at ka (tokens.json ZERO_ORDER_ABS option
    # 1). Dose records use rate = -2 so rxode2 takes the duration from
    # dur(depot), matching RATE in the authors' NONMEM data set.
    f(depot)   <- exp(lfdepot)
    dur(depot) <- d1

    # central is in mg and vc in L, so central / vc is mg/L, which equals the
    # ug/mL scale conventionally used for monoclonal-antibody concentrations.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
