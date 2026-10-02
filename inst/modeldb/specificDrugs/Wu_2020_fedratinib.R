Wu_2020_fedratinib <- function() {
  description <- paste(
    "Two-compartment oral PK model of fedratinib with first-order absorption,",
    "fitted naive-pooled to mean plasma profiles after a single 500 mg dose in",
    "healthy volunteers (Wu 2020, Supplemental Material SM3). This is the",
    "compartmental fit the authors used to seed the distribution inputs of",
    "their Simcyp minimal-PBPK DDI model; the PBPK layer itself is not",
    "included. Typical-value model: no IIV and no residual-error magnitude",
    "were reported."
  )
  reference <- "Wu F, Krishna G, Surapaneni S. Physiologically based pharmacokinetic modeling to assess metabolic drug-drug interaction risks and inform the drug label for fedratinib. Cancer Chemother Pharmacol. 2020;86(4):461-473. doi:10.1007/s00280-020-04131-y"
  vignette <- "Wu_2020_fedratinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "fedratinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fedratinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fedratinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_studies = 4L,
    disease_state = "Healthy volunteers",
    dose_range = "Single oral dose of 500 mg fedratinib",
    regions = "Not reported",
    notes = paste(
      "Mean plasma concentration-time profiles pooled from four healthy-volunteer",
      "studies (TDU12620 single ascending dose, BDR12462 tablet-vs-capsule",
      "bioequivalence, FED12258 and ALI13451 food effect), all given single",
      "500 mg doses (Wu 2020 SM3). Three of the four studies are described as",
      "enrolling healthy male subjects (SM1 Table 1); the fourth (BDR12462) does",
      "not state sex. The number of subjects contributing to the mean profiles is",
      "not reported. The fit used the Phoenix NLME 'naive pooled' algorithm, so",
      "no between-subject variability was estimated."
    )
  )

  ini({
    # Wu 2020 Supplemental Material SM3, compartmental PK parameter table
    # (cited as 'Table 3 in Supplemental Material SM 3' in SM2 'Distribution
    # Parameters'). Phoenix NLME v7, naive-pooled, two-compartment, additive
    # error. No standard errors are printed.
    lka <- log(0.221) ; label("First-order absorption rate constant ka (1/h)") # SM3 Table 3: ka = 0.221 1/h
    lcl <- log(27.2)  ; label("Apparent clearance CL/F (L/h)")                  # SM3 Table 3: CL/F = 27.2 L/h
    lvc <- log(107)   ; label("Apparent central volume V/F (L)")                # SM3 Table 3: V/F = 107 L
    lq  <- log(30.8)  ; label("Apparent intercompartmental clearance CL2/F (L/h)") # SM3 Table 3: CL2/F = 30.8 L/h
    lvp <- log(797)   ; label("Apparent peripheral volume V2/F (L)")             # SM3 Table 3: V2/F = 797 L

    # SM3 names an additive residual-error model but prints no magnitude, and
    # the naive-pooled algorithm estimates no IIV, so none is encoded.
    addSd <- fixed(0) ; label("Additive residual SD (ng/mL; ZERO - magnitude not reported in source)") # SM3: 'additive error model', value not reported
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L: mg/L x 1000 = ng/mL (the units of the paper's
    # exposure tables).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
