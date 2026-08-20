Zhang_2024_XZP5610 <- function() {
  description <- "One-compartment oral PK model for the non-steroidal FXR agonist XZP-5610 (NASH) in healthy Chinese adults, forward-predicted from SD rat and beagle dog preclinical PK by allometric scaling (Zhang 2024): first-order absorption from a depot with bioavailability, linear elimination from central. Every parameter is a cross-species prediction used for first-in-human dose selection, not an estimate fitted to human data; the paper reports no IIV and no residual error for it. The companion whole-body PBPK model in the same paper was built in PK-Sim v11.2 and is not reproducible from the published sources (no ODEs, organ volumes, blood flows, or per-tissue partition coefficients are reported)."
  reference <- "Zhang L, Feng F, Wang X, Liang H, Yao X, Liu D. Dose Prediction and Pharmacokinetic Simulation of XZP-5610, a Small Molecule for NASH Therapy, Using Allometric Scaling and Physiologically Based Pharmacokinetic Models. Pharmaceuticals (Basel). 2024 Mar 13;17(3):369. doi:10.3390/ph17030369"
  vignette <- "Zhang_2024_XZP5610"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what molecule each compartment holds, in what units, in what
  # biological matrix. Verified against Zhang 2024 Section 4.3.2 (one-compartment
  # oral PK model) and Table 6.
  compartmentData <- list(
    depot   = list(analyte = "XZP-5610", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "XZP-5610", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  # No covariate was carried into the model. The allometric scaling used a fixed
  # 60 kg adult body weight (Section 4.3.2) rather than a body-weight covariate,
  # and sex was explicitly argued not to be expected to matter in humans
  # (Discussion: "the sex differences in CYP2C9 and CYP3A activity were small").
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    study          = "First-in-human single ascending dose (SAD) trial in healthy Chinese adults; per-cohort n not reported",
    age_range      = "adults (not reported)",
    weight_range   = "60 kg assumed for the whole population (Section 4.3.2: 'The assumption of an average body weight of 60 kg for adult individuals in the Chinese population was employed')",
    sex_female_pct = NA_real_,
    race_ethnicity = "Chinese",
    disease_state  = "healthy volunteers (the therapeutic target indication is non-alcoholic steatohepatitis)",
    dose_range     = "single oral doses 0.15 to 3 mg (Discussion: 'the initial dose of 0.15 mg and maximum dose of 3 mg were used in the SAD trial'); Figure 4 simulates 0.15, 0.5, 1 and 3 mg",
    regions        = "China",
    notes          = paste(
      "Parameters are cross-species forward predictions, not human estimates.",
      "Source preclinical data: 24 SD rats (mean 0.245 kg) and 24 beagle dogs (mean 7.21 kg),",
      "each split into one IV group and three single-dose oral groups (Section 4.3.1, Tables 1 and 2).",
      "Human clearance was taken from the beagle dog only, because rat hepatocyte metabolic",
      "stability differed markedly from human (Section 2.4, Discussion).",
      "The paper's clinical SAD data were used only to check the PBPK predictions (Figure 3);",
      "no parameter here was fitted to human concentrations, so no IIV and no residual error exist."
    )
  )

  ini({
    # All four parameters are allometric-scaling predictions held constant for the
    # first-in-human dose calculation; none was estimated from human data.

    lka <- fixed(log(1.46))
    label("Absorption rate constant (ka, 1/h)")
    # Section 2.4: "the estimated value of Ka was 1.46 h-1 (range: 0.589-2.34 h-1)".
    # Section 4.3.2: mean of the rat and dog one-compartment NONMEM estimates,
    # which bound human Ka; mean(0.589, 2.34) = 1.4645.

    lcl <- fixed(log(8.3))
    label("Intravenous clearance (CL, L/h)")
    # Section 2.4: "the average CLi.v calculated based on data from beagle dogs was
    # used as the predicted human CLi.v, which was 138 mL/min (8.3 L/h)"; also the
    # Table 5 footnote ("Human CLi.v ... which is 8.3 L/h"). Mean of the four dog
    # rows of Table 3 = mean(96.4, 92.3, 181, 184) = 138.4 mL/min. Table 6 lists the
    # weight-normalised form 0.138 L/h/kg (x 60 kg = 8.28 L/h). 8.3 L/h is the value
    # that reproduces the Table 5 HED column exactly (see vignette).

    lvc <- fixed(log(41.8))
    label("Steady-state volume of distribution (Vss, L)")
    # Section 2.4: "the predicted average value (41.8 L) from SD rats and beagle dogs
    # was used as the predicted human Vss". Mean of the three Vss rows of Table 3 =
    # mean(17.1, 45.3, 63.0) = 41.8 L.

    lfdepot <- fixed(log(0.574))
    label("Oral bioavailability (F)")
    # Section 2.4: "the predicted human oral bioavailability was 57.4%
    # (range: 14.7-57.4%)"; Table 5 footnote: "F: the average value from PK study of
    # beagle dog, which is 57.4%". Mean of the six dog oral F values in Table 2 =
    # mean(43.4, 82.9, 51.2, 38.2, 61.4, 67.1) = 57.37%.
  })

  model({
    # Unit bookkeeping: doses are in mg and the states hold mg, so central / vc is
    # mg/L; the paper reports concentrations in ng/mL, and 1 mg/L = 1000 ng/mL.
    mgPerLToNgPerML <- 1000

    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    Cc <- mgPerLToNgPerML * central / vc
  })
}
