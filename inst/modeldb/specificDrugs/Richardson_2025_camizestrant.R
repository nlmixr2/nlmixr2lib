Richardson_2025_camizestrant <- function() {
  description <- paste(
    "Automated model-search (pyDarwin) two-compartment population PK model",
    "for oral camizestrant in 184 participants pooled from Phase 1 study",
    "NCT03616587 (Richardson 2025). This is the top-ranked structure returned",
    "by the paper's automated Bayesian-optimisation model search and is",
    "identical to the manually developed expert structure apart from the",
    "residual-error model. Absorption is first-order through a depot followed",
    "by three transit compartments, with every transfer governed by a single",
    "rate constant (the paper's KA = KTR). Disposition is two-compartment",
    "linear with apparent clearance and volumes (bioavailability fixed to 1,",
    "so all disposition parameters are apparent CL/F and V/F). Dose enters as",
    "a power-model covariate on CL and Vc, normalised to a reference dose. IIV",
    "is diagonal log-normal and estimated on CL, Vc and the transit rate",
    "constant, with IIV on Vp and Q fixed at 15% CV. Residual error is",
    "combined additive plus proportional. Note that this is a machine-search",
    "structure, not the manually developed expert model.",
    sep = " "
  )
  reference <- paste(
    "Richardson S, Irurzun Arana I, Nowojewski A, Zhou D, Leander J,",
    "Tang W, Dearden R, Gibbs M.",
    "A machine learning approach to population pharmacokinetic modelling",
    "automation. Commun Med. 2025;5:325. doi:10.1038/s43856-025-01054-8.",
    "Parameter estimates from Supplementary Table 3; model structure from",
    "Table 4, Supplementary Data 1 and the pyDarwin model-space files",
    "(template.txt, tokens.json) in the Code Availability repository",
    "https://github.com/samjrrr/autopk_synthetic_example.",
    sep = " "
  )
  vignette <- "Richardson_2025_automated_popPK"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "camizestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "camizestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "camizestrant", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "camizestrant", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "camizestrant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "camizestrant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Administered dose amount of the current dosing record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power model on CL and Vc, normalised to a reference dose:",
        "P = TVP * (DOSE / 100)^e_dose_P. The authors' NONMEM data set carried",
        "DOSE and a separate DOSEREF column and the model term was",
        "MU(P) = LOG(TVP) + LOG(DOSE/DOSEREF) * THETA(cov_P) (template.txt and",
        "tokens.json in the Code Availability repository). Neither the numeric",
        "value of DOSEREF nor the camizestrant dose range is reported anywhere",
        "in the paper or its supplement, so 100 mg is used here as a deliberately",
        "round, arbitrary reference anchor: it is a declared normalisation",
        "constant, not a value read from the source. The tabulated typical",
        "values apply at DOSE = 100 mg. See the validation vignette Assumptions",
        "and deviations section - a different reference dose D rescales each",
        "affected typical value by (100 / D)^e_dose_P and leaves the shape of",
        "the dose dependence unchanged.",
        sep = " "
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species     = "human",
    n_subjects  = 184,
    n_studies   = 1,
    disease_state = "not reported (Phase 1 oncology study NCT03616587)",
    dose_range  = "not reported",
    notes       = paste(
      "184 participants after filtering to oral administration, contributing",
      "2,709 observation records (Richardson 2025 Table 1 and Supplementary",
      "Table 7). Observed concentrations were modelled on the log scale",
      "(Table 1 'Scale type'). The paper reports no baseline demographics for",
      "this cohort; the underlying study design is reported separately in the",
      "NCT03616587 primary publication cited by Richardson 2025.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- apparent (F fixed to 1) typical values at the
    # reference dose. Supplementary Table 3.
    lcl  <- log(62.40)  ; label("Apparent clearance CL/F at the reference dose (L/h)")            # Suppl Table 3 THETA1 = 62.40 L/h
    lvc  <- log(857.55) ; label("Apparent central volume V2/F at the reference dose (L)")         # Suppl Table 3 THETA2 = 857.55 L
    lvp  <- log(439.30) ; label("Apparent peripheral volume V3/F (L)")                            # Suppl Table 3 THETA8 = 439.30 L
    lq   <- log(36.48)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")               # Suppl Table 3 THETA9 = 36.48 L/h
    lktr <- log(2.76)   ; label("Transit / absorption rate constant KA = KTR (1/h)")              # Suppl Table 3 THETA5 = 2.76 1/h

    # Bioavailability was fixed to 1 in the model space (template.txt: F1 = 1),
    # so every disposition parameter above is apparent.
    lfdepot <- fixed(log(1)) ; label("Bioavailability F (unitless)")                  # template.txt F1 = 1   ; INPUT1

    # Dose effects: power model on (DOSE / 100).
    e_dose_cl <- -0.17 ; label("Power exponent on (DOSE / 100) for CL (unitless)")                # Suppl Table 3 THETA3 = -0.17
    e_dose_vc <- -0.34 ; label("Power exponent on (DOSE / 100) for Vc (unitless)")                # Suppl Table 3 THETA4 = -0.34

    # IIV -- diagonal log-normal variances. Supplementary Table 3 reports these
    # as CV%; omega^2 = log(1 + CV^2), the inverse of the
    # CV% = sqrt(exp(omega^2) - 1) * 100 convention verified against the
    # authors' published NONMEM output for the ribociclib model.
    etalcl  ~ 0.148420        # Suppl Table 3 OMEGA(1,1) = 40% CV
    etalvc  ~ 0.155378        # Suppl Table 3 OMEGA(2,2) = 41% CV
    etalktr ~ 0.272771        # Suppl Table 3 OMEGA(3,3) = 56% CV
    etalvp  ~ fixed(0.022251) # Suppl Table 3 OMEGA(4,4) = 15% (fix); tokens.json OMEGAS_V3 0.0225 FIX
    etalq   ~ fixed(0.022251) # Suppl Table 3 OMEGA(5,5) = 15% (fix); tokens.json OMEGAS_Q 0.0225 FIX

    # Residual error -- combined additive plus proportional. The camizestrant
    # data were fitted on the log scale, where the search space's combined
    # error W = SQRT(add^2 + (prop * IPRED)^2) becomes W / IPRED; the two forms
    # are the same model to first order, so it is encoded here on the linear
    # concentration scale.
    addSd  <- 1.01 ; label("Additive residual error (ng/mL)")             # Suppl Table 3 THETA6 = 1.01
    propSd <- 0.32 ; label("Proportional residual error (fraction)")      # Suppl Table 3 THETA7 = 0.32
  })

  model({
    # Individual parameters. The dose effect is a power model on the dose
    # normalised to the 100 mg reference anchor, equivalent to the authors'
    # MU(P) = LOG(TVP) + LOG(DOSE/DOSEREF) * THETA(cov_P).
    cl  <- exp(lcl + etalcl) * (DOSE / 100)^e_dose_cl
    vc  <- exp(lvc + etalvc) * (DOSE / 100)^e_dose_vc
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq  + etalq)
    ktr <- exp(lktr + etalktr)

    # Micro-constants, matching the ADVAN5 parameterisation of the source
    # control stream (k20 = CL/V2, k23 = Q/V2, k32 = Q/V3).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Absorption chain: depot -> transit1 -> transit2 -> transit3 -> central,
    # every transfer governed by the single rate constant ktr (tokens.json
    # TRANSIT_COMP option 2 sets RATE_CONSTANT(...) = KA for all four steps).
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # central is in mg and vc in L, so central / vc is mg/L; the factor of
    # 1000 puts the observation on the ng/mL scale that the additive residual
    # error (1.01) is expressed in.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
