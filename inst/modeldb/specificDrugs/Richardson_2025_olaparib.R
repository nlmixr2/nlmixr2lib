Richardson_2025_olaparib <- function() {
  description <- paste(
    "Automated model-search (pyDarwin) two-compartment population PK model",
    "for oral olaparib tablets in 296 participants pooled from five Phase 1",
    "studies (Richardson 2025). This is the top-ranked structure returned by",
    "the paper's automated Bayesian-optimisation model search. It differs",
    "from the manually developed expert structure in its absorption model, in",
    "not estimating IIV on the peripheral volume, and in carrying a dose",
    "effect on the central volume that the authors considered negligible; the",
    "expert model additionally includes a CL autoinhibition mechanism that was",
    "not available in the paper's model space. Absorption is first-order",
    "through a depot followed by four transit compartments, with every",
    "transfer governed by a single rate constant (the paper's KA = KTR).",
    "Disposition is two-compartment linear with apparent clearance and volumes",
    "(bioavailability fixed to 1, so all disposition parameters are apparent",
    "CL/F and V/F). Dose enters as a power-model covariate on Vc, normalised",
    "to a reference dose. IIV is diagonal log-normal and estimated on CL, Vc,",
    "the transit rate constant and Q, with IIV on Vp fixed at 15% CV.",
    "Residual error is proportional only (the additive term was fixed to",
    "zero). Note that this is a machine-search structure, not the manually",
    "developed expert model.",
    sep = " "
  )
  reference <- paste(
    "Richardson S, Irurzun Arana I, Nowojewski A, Zhou D, Leander J,",
    "Tang W, Dearden R, Gibbs M.",
    "A machine learning approach to population pharmacokinetic modelling",
    "automation. Commun Med. 2025;5:325. doi:10.1038/s43856-025-01054-8.",
    "Parameter estimates from Supplementary Table 5; model structure from",
    "Table 4, Supplementary Data 1 and the pyDarwin model-space files",
    "(template.txt, tokens.json) in the Code Availability repository",
    "https://github.com/samjrrr/autopk_synthetic_example.",
    sep = " "
  )
  vignette <- "Richardson_2025_automated_popPK"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "olaparib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "olaparib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olaparib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Administered dose amount of the current dosing record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power model on Vc, normalised to a reference dose:",
        "Vc = TVVc * (DOSE / 100)^e_dose_vc. The authors' NONMEM data set",
        "carried DOSE and a separate DOSEREF column and the model term was",
        "MU(V2) = LOG(TVV2) + LOG(DOSE/DOSEREF) * THETA(cov_V2) (template.txt",
        "and tokens.json in the Code Availability repository). Neither the",
        "numeric value of DOSEREF nor the olaparib dose range is reported",
        "anywhere in the paper or its supplement, so 100 mg is used here as a",
        "deliberately round, arbitrary reference anchor: it is a declared",
        "normalisation constant, not a value read from the source. The",
        "tabulated typical values apply at DOSE = 100 mg. The authors",
        "considered this dose effect negligible because of the small THETA",
        "value (0.17 with 34% relative standard error). See the validation",
        "vignette Assumptions and deviations section - a different reference",
        "dose D rescales Vc by (100 / D)^e_dose_vc and leaves the shape of the",
        "dose dependence unchanged.",
        sep = " "
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 296,
    n_studies     = 5,
    disease_state = "not reported (Phase 1 oncology studies)",
    dose_range    = "not reported",
    notes         = paste(
      "296 participants after filtering to the tablet formulation only,",
      "contributing 7,397 observation records, pooled from NCT01813474,",
      "NCT01921140, NCT01900028, NCT01929603 and NCT00777582 (Richardson 2025",
      "Table 1 and Supplementary Table 7). The capsule-formulation data",
      "available in the source studies were excluded. Observed concentrations",
      "were modelled on the normal scale (Table 1 'Scale type'). The paper",
      "reports no baseline demographics for this cohort.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- apparent (F fixed to 1) typical values at the
    # reference dose. Supplementary Table 5.
    lcl  <- log(5.21)  ; label("Apparent clearance CL/F (L/h)")                                # Suppl Table 5 THETA1 = 5.21 L/h
    lvc  <- log(35.06) ; label("Apparent central volume V2/F at the reference dose (L)")       # Suppl Table 5 THETA2 = 35.06 L
    lvp  <- log(18.77) ; label("Apparent peripheral volume V3/F (L)")                          # Suppl Table 5 THETA7 = 18.77 L
    lq   <- log(0.41)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")             # Suppl Table 5 THETA8 = 0.41 L/h
    lktr <- log(9.22)  ; label("Transit / absorption rate constant KA = KTR (1/h)")            # Suppl Table 5 THETA4 = 9.22 1/h

    # Bioavailability was fixed to 1 in the model space (template.txt: F1 = 1),
    # so every disposition parameter above is apparent.
    lfdepot <- fixed(log(1)) ; label("Bioavailability F (unitless)")               # template.txt F1 = 1   ; INPUT1

    # Dose effect: power model on (DOSE / 100).
    e_dose_vc <- 0.17 ; label("Power exponent on (DOSE / 100) for Vc (unitless)")              # Suppl Table 5 THETA3 = 0.17

    # IIV -- diagonal log-normal variances. Supplementary Table 5 reports these
    # as CV%; omega^2 = log(1 + CV^2), the inverse of the
    # CV% = sqrt(exp(omega^2) - 1) * 100 convention verified against the
    # authors' published NONMEM output for the ribociclib model.
    etalcl  ~ 0.264285        # Suppl Table 5 OMEGA(1,1) = 55% CV
    etalvc  ~ 0.097490        # Suppl Table 5 OMEGA(2,2) = 32% CV
    etalktr ~ 0.307485        # Suppl Table 5 OMEGA(3,3) = 60% CV
    etalvp  ~ fixed(0.022251) # Suppl Table 5 OMEGA(4,4) = 15% (fix); tokens.json OMEGAS_V3 0.0225 FIX
    etalq   ~ 1.528228        # Suppl Table 5 OMEGA(5,5) = 190% CV -- the paper flags this as the one implausibly high OMEGA in the top models (penalty 81.53)

    # Residual error -- proportional only. The additive term was fixed to zero
    # (tokens.json ERRORFIX option 1: 0 FIX ; add_error), so the combined
    # form W = SQRT(add^2 + (prop * IPRED)^2) collapses to prop * IPRED.
    addSd  <- fixed(0) ; label("Additive residual error (ng/mL)")   # Suppl Table 5 THETA5 = 0 (fix)
    propSd <- 0.36     ; label("Proportional residual error (fraction)")       # Suppl Table 5 THETA6 = 0.36
  })

  model({
    # Individual parameters. The dose effect is a power model on the dose
    # normalised to the 100 mg reference anchor, equivalent to the authors'
    # MU(V2) = LOG(TVV2) + LOG(DOSE/DOSEREF) * THETA(cov_V2).
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc) * (DOSE / 100)^e_dose_vc
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq  + etalq)
    ktr <- exp(lktr + etalktr)

    # Micro-constants, matching the ADVAN5 parameterisation of the source
    # control stream (k20 = CL/V2, k23 = Q/V2, k32 = Q/V3).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Absorption chain: depot -> transit1 -> ... -> transit4 -> central, every
    # transfer governed by the single rate constant ktr (tokens.json
    # TRANSIT_COMP option 3 sets RATE_CONSTANT(...) = KA for all five steps).
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(central)     <-  ktr * transit4 + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # central is in mg and vc in L, so central / vc is mg/L; the factor of
    # 1000 puts the observation on the ng/mL scale used throughout this
    # extraction (the additive error is fixed to zero here, so the scale
    # affects only the reported concentration units).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
