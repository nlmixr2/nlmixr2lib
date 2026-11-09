Chae_2025_beeVenomPLA2_human_simple <- function() {
  description <- paste(
    "Two-compartment PK model with first-order absorption for bee venom",
    "phospholipase A2 (bvPLA2) after subcutaneous new composition bee venom",
    "(NCBV) in humans, projected by SIMPLE allometry (Y = a * BW^b) from the",
    "mouse, rat and beagle-dog population PK fits of Chae 2025 (Table 3).",
    "A forward first-in-human projection, not a fit to human data: typical",
    "values only, no between-subject variability or residual error. The",
    "absorption rate constant is not printed for the human projection and is",
    "the rat ka, which reproduces the paper's Table 4 simulation."
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

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference weight 70 kg, the human weight at which Chae 2025 Table 3",
        "reports the projection. Each parameter scales as (WT/70)^b with the",
        "simple-allometry exponent b printed in the Results ('estimated",
        "exponents of CL/F, Vc/F, Vp/F, and Q/F were 0.847, 0.842, 0.729,",
        "and 0.971'); WT = 70 returns Table 3 exactly."
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "bvPLA2", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bvPLA2", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bvPLA2", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human (projected by allometry from mouse, rat and beagle dog)",
    n_subjects = 0L,
    n_studies = 0L,
    weight_range = "70 kg reference",
    disease_state = "not applicable -- forward projection for a first-in-human study (NCBV is in development for Alzheimer's disease)",
    dose_range = "5.5 mg single SC dose (Table 4: the dose matching the mouse efficacious AUC with this projection)",
    regions = "not applicable",
    notes = paste(
      "Simple allometry Y = a * BW^b fitted across ICR mouse, SD rat and",
      "beagle dog Table 2 parameters (Chae 2025 Methods, Table 3; R2 0.781 for",
      "CL/F). Sibling projections using brain weight (BrW) and maximum",
      "lifespan potential (MLP) ship as Chae_2025_beeVenomPLA2_human_brw and",
      "Chae_2025_beeVenomPLA2_human_mlp; the paper recommends BrW."
    )
  )

  ini({
    # Chae 2025 Table 3, 'Simple' block, predicted values in a 70 kg human.
    # All values are fixed: they are arithmetic consequences of the animal
    # fits, not estimates from human data.
    lcl <- fixed(log(12.9)); label("Apparent clearance CL/F at 70 kg (L/h)") # Table 3 Simple: CL/F = 12.9 L/h (R2 0.781)
    lvc <- fixed(log(19.5)); label("Apparent central volume Vc/F at 70 kg (L)") # Table 3 Simple: Vc/F = 19.5 L (R2 1.00)
    lvp <- fixed(log(67.3)); label("Apparent peripheral volume Vp/F at 70 kg (L)") # Table 3 Simple: Vp/F = 67.3 L (R2 0.966)
    lq <- fixed(log(14.7)); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 3 Simple: Q/F = 14.7 L/h (R2 0.862)

    # Human ka is not printed in Table 3. The rat Table 2 ka (0.803 1/h)
    # reproduces Table 4 (Cmax 78.3 vs 78.4 ng/mL, Tmax 0.99 vs 1 h, terminal
    # t1/2 7.4 vs 7.39 h after 5.5 mg) -- see the vignette.
    lka <- fixed(log(0.803)); label("Absorption rate constant (1/h)") # Table 2 rat ka = 0.803 1/h; back-solved against Table 4

    # Simple-allometry exponents b (Results, 'Allometric scaling for the
    # prediction of human PK parameters').
    e_wt_cl <- fixed(0.847); label("Allometric weight exponent on CL/F (unitless)") # Results: CL/F exponent 0.847
    e_wt_vc <- fixed(0.842); label("Allometric weight exponent on Vc/F (unitless)") # Results: Vc/F exponent 0.842
    e_wt_vp <- fixed(0.729); label("Allometric weight exponent on Vp/F (unitless)") # Results: Vp/F exponent 0.729
    e_wt_q <- fixed(0.971); label("Allometric weight exponent on Q/F (unitless)") # Results: Q/F exponent 0.971

    # No variability or residual error is reported for the human projection.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)") # not reported for the human projection
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    q <- exp(lq) * (WT / 70)^e_wt_q

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
