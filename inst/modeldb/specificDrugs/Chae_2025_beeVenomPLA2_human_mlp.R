Chae_2025_beeVenomPLA2_human_mlp <- function() {
  description <- paste(
    "Two-compartment PK model with first-order absorption for bee venom",
    "phospholipase A2 (bvPLA2) after subcutaneous new composition bee venom",
    "(NCBV) in a 70 kg human, projected by maximum-lifespan-potential (MLP)",
    "corrected allometry from the mouse, rat and beagle-dog population PK fits",
    "of Chae 2025 (Table 3, as printed). Typical values only, no",
    "between-subject variability or residual error. Note: the printed Table 3",
    "CL/F (0.612 L/h) does not reproduce the paper's own MLP simulation",
    "(Table 4, Figure 5C); see the vignette."
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
    dose_range = "3.5 mg single SC dose (Table 4)",
    regions = "not applicable",
    notes = paste(
      "MLP-corrected allometry across ICR mouse, SD rat and beagle dog Table 2",
      "parameters (Chae 2025 Methods, Table 3; R2 0.904 for CL/F). The",
      "regression coefficients and the lifespans used are not printed, so only",
      "the 70 kg prediction is encoded and there is no weight covariate.",
      "Table 3 and the Results prose give CL/F = 0.612 L/h, but Table 4 lists",
      "CL/F 0.00917 for MLP (read as 9.17 L/h on the scale of the Simple and",
      "BrW entries, consistent with dose/AUC = 3.5/0.378 = 9.26 L/h) and a",
      "terminal t1/2 of 3.35 h. This file carries the Table 3 values as",
      "printed; with them the terminal t1/2 is about 53 h."
    )
  )

  ini({
    # Chae 2025 Table 3, 'MLP' block, predicted values in a 70 kg human,
    # encoded as printed. All values are fixed: they are arithmetic
    # consequences of the animal fits, not estimates from human data.
    lcl <- fixed(log(0.612)); label("Apparent clearance CL/F (L/h)") # Table 3 MLP: CL/F = 0.612 L/h (R2 0.904); Table 4 lists 0.00917 -- see vignette
    lvc <- fixed(log(10.0)); label("Apparent central volume Vc/F (L)") # Table 3 MLP: Vc/F = 10.0 L (R2 0.999)
    lvp <- fixed(log(34.7)); label("Apparent peripheral volume Vp/F (L)") # Table 3 MLP: Vp/F = 34.7 L (R2 0.993)
    lq <- fixed(log(7.59)); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 MLP: Q/F = 7.59 L/h (R2 0.940)

    # Human ka is not printed in Table 3. The rat Table 2 ka (0.803 1/h)
    # reproduces the Simple and BrW Table 4 simulations and is used for all
    # three human projections.
    lka <- fixed(log(0.803)); label("Absorption rate constant (1/h)") # Table 2 rat ka = 0.803 1/h; back-solved against Table 4 Simple/BrW

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
