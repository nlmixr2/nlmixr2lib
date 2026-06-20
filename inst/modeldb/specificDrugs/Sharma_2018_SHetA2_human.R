Sharma_2018_SHetA2_human <- function() {
  description <- paste(
    "Allometrically-scaled human (70 kg) projection.",
    "Two-compartment intravenous PK model for SHetA2 (a flexible",
    "heteroarotinoid anti-cancer / chemoprevention drug) with",
    "disposition parameters scaled from preclinical mouse / rat / dog",
    "data via simple allometry (CL = a * BW^b) on a log-log plot",
    "(R^2 = 0.91-0.99 across CL, V1, V2, CLD; Sharma 2018 Fig 5).",
    "Clearance uses the maximum-life-span-potential (MLP) correction",
    "to account for SHetA2's hepatic metabolism (CL_MLP = 17.3 L/h vs",
    "41.0 L/h by simple allometry). The model carries no parametric",
    "oral absorption -- the source paper simulated the oral profile",
    "(Fig 6) externally by linking these disposition parameters to the",
    "Advanced Compartmental Absorption and Transit (ACAT) model in",
    "GastroPlus 9.5, because the preclinical kA values did not",
    "correlate across species and could not be projected to humans. The",
    "predicted human oral bioavailability was 18.8% (range 7.4-42%) at",
    "10 mg/kg, very close to the maximum extent of absorption observed",
    "in preclinical species at doses <100 mg/kg (18.6%). Parameter",
    "values from Sharma 2018 Results (Allometric scaling) and",
    "Prediction of human pharmacokinetics."
  )
  reference <- paste(
    "Sharma A, Benbrook DM, Woo S. (2018).",
    "Pharmacokinetics and interspecies scaling of a novel,",
    "orally-bioavailable anti-cancer drug, SHetA2.",
    "PLoS ONE 13(4):e0194046.",
    "doi:10.1371/journal.pone.0194046.",
    sep = " "
  )
  vignette <- "Sharma_2018_SHetA2"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "human (allometrically-scaled projection; no clinical PK data fit)",
    n_subjects     = 0L,
    n_studies      = 0L,
    age_range      = "(projection; reference adult)",
    weight_median  = "70 kg (interspecies-scaling reference)",
    sex_female_pct = NA,
    disease_state  = "Projection for first-in-human dosing in Phase 0 clinical trial of SHetA2 (IRB #5407, planned at the time of publication)",
    dose_range     = "Proposed starting dose 2 mg/kg PO; escalating toward 10 mg/kg PO (target plasma Cmax 4 uM = 1600 ng/mL) per Sharma 2018 Discussion. The published GastroPlus simulation used a single 10 mg/kg dose in an 85-kg subject under the fasted state for n=100 virtual subjects (Fig 6).",
    regions        = "USA (planned clinical trial)",
    notes          = paste(
      "Disposition parameters scaled to a 70-kg human via simple",
      "allometry P = a * BW^b fit on log-log plots across mouse, rat, and",
      "dog values (Sharma 2018 Fig 5; R^2 = 0.91-0.99 across all four",
      "PK parameters). Predicted CL = 17.3 L/h uses the maximum",
      "life-span-potential (MLP) correction CL_MLP = a * BW^b * MLP",
      "(reported in Results, Allometric scaling) because simple",
      "allometry overpredicts human CL for hepatically-metabolised",
      "small molecules; the simple-allometry value is 41.0 L/h. The",
      "model contains no parametric oral absorption -- in the source",
      "paper, the human oral PK profile (Fig 6) was simulated externally",
      "by linking these disposition parameters to the GastroPlus 9.5",
      "ACAT model (physiologically-based gut absorption with",
      "compartmental pH, transit time, lengths and radii of each GI",
      "compartment, plasma protein binding, stomach volume, hepatic",
      "blood flow, gut enzyme/transporter distribution). The mean",
      "predicted human bioavailability was 18.8% (range 7.4-42%) at 10",
      "mg/kg fasted (Sharma 2018 Results, Prediction of human",
      "pharmacokinetics). Body weight is not a covariate of this model;",
      "to use a different reference body weight, the user should",
      "re-scale CL, V1, V2, CLD outside the model. See Sharma 2018",
      "Allometric scaling section and Results, Prediction of human",
      "pharmacokinetics paragraph."
    )
  )

  ini({
    # --------------------------------------------------------------
    # Structural PK -- Sharma 2018 Results (Allometric scaling and
    # Prediction of human pharmacokinetics) for a 70-kg human.
    # CL uses the MLP-corrected allometric prediction (17.3 L/h);
    # V1, V2, CLD are the simple-allometry projections.
    # No CV% / RSE values are reported by the paper for the
    # allometric projections (they are point predictions from the
    # log-log regressions in Fig 5).
    # --------------------------------------------------------------
    lcl     <- log(17.3)    ; label("Clearance CL (L/h) for 70-kg human (MLP-corrected allometry)")    # Results, Allometric scaling: MLP-CL = 17.3 L/h
    lvc     <- log(36.2)    ; label("Central volume of distribution V1 (L) for 70-kg human")           # Results, Prediction of human pharmacokinetics: V1 = 36.2 L
    lq      <- log(15.2)    ; label("Distributional clearance CLD (L/h) for 70-kg human")              # Results, Prediction of human pharmacokinetics: CLD = 15.2 L/h
    lvp     <- log(68.5)    ; label("Peripheral volume V2 (L) for 70-kg human")                        # Results, Prediction of human pharmacokinetics: V2 = 68.5 L

    # --------------------------------------------------------------
    # Residual error: there is no fitted residual error -- the human
    # profile is a deterministic projection. Encoding a typical 10%
    # proportional SD as a placeholder so the model is a valid
    # nlmixr2 simulation object; downstream consumers can override
    # via zeroRe() for typical-value-only output. NOT a fitted value.
    # --------------------------------------------------------------
    propSd  <- fixed(0.10)  ; label("Proportional residual error (fraction); placeholder, not from paper") # Allometric projection -- no fitted residual error
  })

  model({
    # 1. Individual parameters (typical-value model; no eta).
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)

    # 2. Two-compartment IV disposition only (no parametric oral
    #    absorption -- the source paper simulated the oral profile
    #    externally via GastroPlus 9.5 ACAT, which is not reproducible
    #    inside nlmixr2lib). State amounts in mg.
    d/dt(central)     <- -(cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-   q / vc * central - q / vp * peripheral1

    # 3. Observation: plasma concentration in ng/mL.
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
