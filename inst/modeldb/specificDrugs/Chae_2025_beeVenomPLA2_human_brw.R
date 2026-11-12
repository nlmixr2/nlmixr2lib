Chae_2025_beeVenomPLA2_human_brw <- function() {
  description <- paste(
    "Two-compartment PK model with first-order absorption for bee venom",
    "phospholipase A2 (bvPLA2) after subcutaneous new composition bee venom",
    "(NCBV) in a 70 kg human, projected by brain-weight (BrW) corrected",
    "allometry from the mouse, rat and beagle-dog population PK fits of Chae",
    "2025 (Table 3). The projection the paper recommends for the first-in-",
    "human dose (1.3 mg). Typical values only, no between-subject variability",
    "or residual error. The absorption rate constant is not printed for the",
    "human projection and is the rat ka, which reproduces the paper's Table 4",
    "simulation."
  )
  reference <- paste(
    "Chae SU, Min JS, Jo SJ, Lee CB, Park J, Bae SH, Bae SK.",
    "Prediction of first-in-human dose for new composition bee venom based on",
    "allometric scaling and pharmacokinetic modeling approach.",
    "Transl Clin Pharmacol. 2025;33(1):27-39. doi:10.12793/tcp.2025.33.e4",
    sep = " "
  )
  vignette <- "Chae_2025_beeVenomPLA2"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "bvPLA2", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bvPLA2", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bvPLA2", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human (projected by allometry from mouse, rat and beagle dog)",
    n_subjects = 0L,
    n_studies = 0L,
    weight_range = "70 kg (single projected subject)",
    disease_state = "not applicable -- forward projection for a first-in-human study (NCBV is in development for Alzheimer's disease)",
    dose_range = "1.3 mg single SC dose (Table 4; the recommended first-in-human dose)",
    regions = "not applicable",
    notes = paste(
      "Brain-weight corrected allometry across ICR mouse, SD rat and beagle dog",
      "Table 2 parameters (Chae 2025 Methods, Table 3; R2 0.973 for CL/F, the",
      "highest of the three methods). The regression coefficients and the",
      "brain weights used are not printed, so only the 70 kg prediction is",
      "encoded and there is no weight covariate."
    )
  )

  ini({
    # Chae 2025 Table 3, 'BrW' block, predicted values in a 70 kg human.
    # All values are fixed: they are arithmetic consequences of the animal
    # fits, not estimates from human data.
    lcl <- fixed(log(3.25)); label("Apparent clearance CL/F (L/h)") # Table 3 BrW: CL/F = 3.25 L/h (R2 0.973)
    lvc <- fixed(log(4.92)); label("Apparent central volume Vc/F (L)") # Table 3 BrW: Vc/F = 4.92 L (R2 0.993)
    lvp <- fixed(log(17.0)); label("Apparent peripheral volume Vp/F (L)") # Table 3 BrW: Vp/F = 17.0 L (R2 1.00)
    lq <- fixed(log(3.72)); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 BrW: Q/F = 3.72 L/h (R2 0.932)

    # Human ka is not printed in Table 3. The rat Table 2 ka (0.803 1/h)
    # reproduces Table 4 (Cmax 73.4 vs 73.3 ng/mL, Tmax 0.99 vs 1 h, terminal
    # t1/2 7.4 vs 7.40 h after 1.3 mg) -- see the vignette.
    lka <- fixed(log(0.803)); label("Absorption rate constant (1/h)") # Table 2 rat ka = 0.803 1/h; back-solved against Table 4

    # No variability or residual error is reported for the human projection.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)") # not reported for the human projection
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses in mg and vc in L give mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
