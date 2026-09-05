Ma_2023_difloxacin_crucianCarp <- function() {
  description <- paste(
    "Veterinary (crucian carp, Carassius auratus).",
    "One-compartment population PK model with first-order absorption and no",
    "absorption lag for difloxacin after a single 20 mg/kg body-weight oral",
    "(gavage) dose in crucian carp held at 21.3 +/- 1.2 degC. Ma 2023 fitted 25",
    "sparsely sampled fish (three plasma samples each, drawn from a 15-point",
    "0.25-120 h grid) in Phoenix NLME 8.1. Doses, clearance and volume are all",
    "body-weight normalised in the source (dose in mg/kg, tvV in L/kg, tvCL in",
    "L/h/kg), so the model is driven with amt in mg/kg and returns Cc directly",
    "in ug/mL. Final ('covariance') model, Table 3: tvKa 1.18 1/h,",
    "tvV 14.18 L/kg, tvCL 0.20 L/h/kg, with exponential between-fish",
    "variability on all three (omega^2 = 0.69, 1.12 and 0.98) and a Phoenix",
    "log-additive residual error, stdev0 = 0.16 on the log scale, encoded as",
    "lnorm(expSd). Body weight was screened as a covariate on CL/F, V/F and Ka",
    "and rejected (delta -2LL = 3.02, below the 6.64 significance threshold),",
    "so no covariate enters the model. The authors fitted a FULL",
    "variance-covariance omega block over etaKa, etaV and etaCL but published",
    "only the three diagonal variances; the off-diagonals are therefore absent",
    "here and the etas are independent. See the vignette for the consequences.",
    sep = " "
  )
  reference <- paste(
    "Ma K-L, Yang F, Zhang M, Chen J-C, Duan M-H, Li Z-E, Dai Y, Liu Y,",
    "Jin Y-G, Yang F. (2023). Population Pharmacokinetics of Difloxacin in",
    "Crucian Carp (Carassius auratus) after a Single Oral Administration.",
    "Veterinary Sciences 10(7):416. doi:10.3390/vetsci10070416.",
    sep = " "
  )
  vignette <- "Ma_2023_difloxacin_crucianCarp"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states hold difloxacin; the amounts are
  # body-weight normalised (mg/kg) because Ma 2023 reports the dose in
  # mg/kg BW and V/F, CL/F per kg BW. verified = TRUE: the analyte is
  # difloxacin throughout (Section 2.3, oral gavage) and the sampled matrix
  # is plasma (Section 2.3: tail-vein blood centrifuged to plasma).
  compartmentData <- list(
    depot = list(analyte = "difloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "difloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  # Ma 2023 screened body weight and rejected it; nothing else could be
  # screened because sex cannot be determined non-lethally in crucian carp
  # (Discussion, paragraph 2). No covariate is referenced in model(), so the
  # screen is recorded here rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight of the individual fish. Ma 2023 Section 2.5 screened BW by the 'shotgun' method against CL/F, V/F and Ka simultaneously; the change in -2LL was 3.02, below the 6.64 chi-square threshold for p < 0.01, so BW was dropped and no covariate model is reported (Section 3.3 and Discussion paragraph 2). Table 1 lists the per-fish weights (0.21-0.40 kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and rejected; not retained in the final model. Note that the model is nevertheless implicitly weight-normalised: the dose is given in mg/kg BW and V/F and CL/F are reported per kg BW, i.e. linear (exponent 1) weight scaling is built into the units rather than fitted as a covariate. The rejected covariate test was therefore for a departure from that proportionality, not for the presence of any weight effect at all."
    )
  )

  population <- list(
    species            = "crucian carp (Carassius auratus)",
    n_subjects         = 25L,
    n_studies          = 1L,
    weight_range       = "0.21-0.40 kg",
    weight_mean        = "0.298 kg",
    dose_range         = "single 20 mg/kg body weight oral gavage (10 mg/mL difloxacin hydrochloride in 0.9% saline)",
    disease_state      = "healthy",
    regions            = "China (Henan University of Science and Technology, Luoyang)",
    water_temperature  = "21.3 +/- 1.2 degC",
    water_quality      = "pH ~7.4; total ammonia nitrogen <= 0.7 mg/L; nitrite nitrogen < 0.07 mg/L; dissolved oxygen > 7.5 mg/L",
    design             = "30 fish in 6 equal groups of 5; one group served as an undosed control supplying blank plasma, leaving 25 dosed fish. At least 10 days of acclimation before dosing (Section 2.2).",
    sampling           = "Sparse: each dosed fish sampled 3 times from the tail vein. Table 1 shows five repeating schedules of three times each, cycling across the 25 fish: (0.25, 6, 36), (0.5, 8, 48), (1, 12, 72), (2, 16, 96) and (4, 24, 120) h post-dose, so all 15 nominal times are covered by 5 fish each.",
    bioanalysis        = "HPLC-UV at 276 nm; calibration 0.1-5 ug/mL, LOQ 0.1 ug/mL, mean recovery 81.46%, intra-day CV 0.74-7.72%, inter-day CV 4.81-7.66% (Sections 2.4 and 3.1).",
    notes              = "Baseline weights are listed per fish in Table 1; the group description is in Section 2.2. Only one water temperature was studied, and the authors note (Discussion) that difloxacin exposure in carp is strongly temperature dependent, so the parameters here should not be extrapolated to other rearing temperatures."
  )

  ini({
    # =================================================================
    # Structural parameters -- Ma 2023 Table 3, FINAL model
    # ("Final model (Covariance model)"), not the basic model.
    #
    # All three are apparent (per bioavailability) values and are body-
    # weight normalised: V/F is L/kg and CL/F is L/h/kg, matched to a
    # dose expressed in mg/kg BW. No absorption lag was fitted.
    # =================================================================
    lka <- log(1.18)
    label("Absorption rate constant Ka (1/h)")  # Ma 2023 Table 3, final (covariance) model: tvKa = 1.18 1/h (RSE 12.55%, 95% CI 0.88-1.47). Basic model was 1.16.
    lvc <- log(14.18)
    label("Apparent central volume of distribution V/F (L/kg)")  # Ma 2023 Table 3, final (covariance) model: tvV = 14.18 L/kg (RSE 18.88%, 95% CI 8.83-19.53). Basic model was 14.34.
    lcl <- log(0.20)
    label("Apparent clearance CL/F (L/h/kg)")  # Ma 2023 Table 3, final (covariance) model: tvCL = 0.20 L/h/kg (RSE 17.96%, 95% CI 0.13-0.27). Basic model was 0.18. Table 4 implies an unrounded 20/102.05 = 0.196 L/h/kg; the printed 2-dp Table 3 value is used here (see vignette).

    # =================================================================
    # Between-fish (inter-individual) variability
    # -- Ma 2023 Table 3, "Inter-Individual Variation omega^2 (RSE%)"
    #    column, FINAL model row block. The column header states these
    #    are variances (omega^2), not SDs or CV%.
    #
    # Ma 2023 Section 2.5 states the final model used a FULL
    # variance-covariance matrix over (etaV, etaCL, etaKa), but Table 3
    # publishes only the three diagonal variances -- no covariances and
    # no correlations appear anywhere in the paper. Rather than invent
    # off-diagonal values, the etas are left independent here. This is a
    # documented deviation from the published structure; see the
    # vignette's "Assumptions and deviations" section.
    # =================================================================
    etalka ~ 0.69  # Ma 2023 Table 3, final model omega^2 for etaKa = 0.69 (RSE 35.14%); basic model 0.60
    etalvc ~ 1.12  # Ma 2023 Table 3, final model omega^2 for etaV  = 1.12 (RSE 28.54%); basic model 1.16
    etalcl ~ 0.98  # Ma 2023 Table 3, final model omega^2 for etaCL = 0.98 (RSE 27.08%); basic model 0.81

    # =================================================================
    # Residual unexplained variability
    #
    # Ma 2023 Section 2.5: "In both the basic and final population
    # models, the log-additive error model was used". Phoenix NLME's
    # log-additive error model is
    #     CObs = C * exp(eps),  eps ~ N(0, stdev0^2)
    # i.e. an additive residual on ln(concentration), which is exactly
    # nlmixr2's `~ lnorm(expSd)`. stdev0 is therefore a LOG-SCALE SD and
    # is dimensionless; Table 3's "ug/mL" entry in its Units column is a
    # Phoenix reporting artefact (the same column also mis-assigns units
    # across the merged parameter rows) and is not used here.
    # =================================================================
    expSd <- 0.16
    label("Residual SD on the natural-log concentration scale (log-normal; dimensionless)")  # Ma 2023 Table 3, final (covariance) model: stdev0 = 0.16 (RSE 11.06%, 95% CI 0.12-0.19). Basic model was 0.19.
  })

  model({
    # Individual parameters. Exponential (log-normal) between-fish
    # variability, matching Phoenix NLME's default multiplicative eta
    # structure for a "stdev0"/log-additive PK model.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    kel <- cl / vc

    # Ma 2023 fitted a one-compartment model with first-order oral
    # absorption and first-order elimination; no lag time, no peripheral
    # compartment and no estimated bioavailability (V and CL are apparent,
    # i.e. already divided by F).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Doses enter as mg/kg body weight and vc is L/kg, so central/vc is
    # mg/L == ug/mL, the paper's concentration unit.
    Cc <- central / vc

    Cc ~ lnorm(expSd)
  })
}
