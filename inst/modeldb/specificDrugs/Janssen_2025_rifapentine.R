Janssen_2025_rifapentine <- function() {
  description <- "One-compartment population PK model for oral rifapentine in South African adults on once-weekly rifapentine/isoniazid (3HP) tuberculosis preventive therapy (Janssen 2025). Absorption is a Savic transit chain (4 transit compartments, absorption rate equal to the transit rate) whose mean transit time is shortened when the tablets are suspended in water; IOV on bioavailability and MTT is correlated and the MTT IOV is smaller after suspended tablets. Parameters are for a 900 mg dose."
  reference <- paste(
    "Janssen S, Vongjarudech T, Karlsson MO, Upton CM, Garcia-Prats AJ,",
    "Diacon AH, Wiesner L, Sachs T, van der Laan LE, Salazar-Austin N,",
    "Hesseling AC, Svensson EM. Bridging to Paediatric Dosing: Relative",
    "Bioavailability of Suspended Rifapentine and Isoniazid in an Open-Label",
    "Randomized Trial in Adults on Tuberculosis Preventive Therapy.",
    "Clin Pharmacokinet. 2025;64(12):1885.",
    "doi:10.1007/s40262-025-01576-3"
  )
  vignette <- "Janssen_2025_rifapentine_isoniazid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "rifapentine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifapentine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_SUSPENSION = list(
      description = "Suspended-tablet indicator: 1 = tablets quartered, crushed and suspended in 30 mL distilled water immediately before dosing (container rinsed with a further 70 mL); 0 = tablets swallowed whole with 100 mL water",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablets swallowed whole; the typical-value MTT reference in Janssen 2025 Table 2)",
      notes = "Per-dose (per-occasion) indicator: each participant received both forms in a crossover. Shortens the typical MTT by 22.2% and scales the IOV on MTT by (1 - 0.458). The three products (two fixed-dose-combination brands and standalone rifapentine) share one set of parameters in the final model; brand effects on bioavailability were assessed by log-likelihood profiling only.",
      source_name = "FORM (1 = whole tablet, 2 = suspension); FORM_SUSPENSION = FORM - 1"
    ),
    OCC = list(
      description = "Pharmacokinetic occasion (weekly dose number with intensive sampling)",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a; OCC = 1..4 selects the per-occasion IOV etas, any other value applies no IOV",
      notes = "Occasions 1-4 are the first four weekly doses (days 1, 8, 15 and 22), each followed by sampling to 48 h. Doses 5-12 had no PK sampling.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 838L,
    age_range = "18-55 years (inclusion); median 37.5 (IQR 26.5-45.3)",
    weight_range = "40-90 kg (inclusion); median 61.1 kg (IQR 51.1-66)",
    sex_female_pct = 83.3,
    race_ethnicity = c(MixedRace = 91.7, Black = 8.3),
    disease_state = "Adults with a positive tuberculosis infection test and recent household contact with drug-susceptible pulmonary TB, on 3HP tuberculosis preventive therapy; HIV and chronic hepatitis B excluded",
    dose_range = "Rifapentine 900 mg once weekly for 12 weeks (with isoniazid 900 mg), as whole or water-suspended tablets of two fixed-dose combinations (rifapentine 300 mg / isoniazid 300 mg) or standalone 300 mg tablets; given 30 min after a standard 670 kcal breakfast",
    regions = "Cape Town, South Africa (single centre)",
    notes = "Open-label randomized crossover; each participant received two of three products, one whole and one suspended, across the four sampled weekly doses. Demographics from Janssen 2025 Table 1 (FFM median 38.6 kg). BLQ samples (3.0% within 48 h; all at 168 h) were omitted."
  )

  ini({
    lcl <- log(1.03)
    label("Apparent clearance CL/F (L/h)") # Table 2, CL/F = 1.03 (90% CI 0.93-1.11)
    lvc <- log(26.56)
    label("Apparent volume of distribution V/F (L)") # Table 2, V/F = 26.56 (90% CI 24.98-27.96)
    # Table 2 prints MTT as 2.34 with a 90% CI of 2.54-2.97, which does not
    # contain it. The final-model control stream (ESM 2, THETA(3) = 2.75047)
    # sits inside that interval, so the stream value is used.
    lmtt <- log(2.75)
    label("Mean transit time of the absorption chain, whole tablet (h)") # ESM 2 $THETA TVMTT = 2.75047; Table 2 MTT 90% CI 2.54-2.97
    lfdepot <- fixed(log(1))
    label("Reference bioavailability F (unitless)") # Table 2, F = 1 (Fixed)

    e_susp_mtt <- -0.222
    label("Fractional change in MTT for suspended tablets (unitless)") # Table 2, Effect of suspension form on MTT = -22.2% (90% CI -28.7 to -15)
    e_susp_iovmtt <- -0.458
    label("Fractional change in the IOV of MTT for suspended tablets (unitless)") # Table 2, Effect of suspension form on MTT variation = -45.8% (90% CI -57.1 to -34.8)

    # IIV and IOV variances from Table 2 %CV via omega^2 = log(1 + CV^2)
    # (Table 2 footnote c).
    etalcl ~ log(1 + 0.199^2) # Table 2, IIV on CL 19.9% CV
    etalfdepot ~ log(1 + 0.173^2) # Table 2, IIV on F 17.3% CV

    # IOV on F (12.2% CV) and MTT (39.2% CV), correlation 0.37 (Table 2);
    # ESM 2 $OMEGA BLOCK(2) with three BLOCK(2) SAME for occasions 2-4.
    etaiov_lfdepot_1 + etaiov_lmtt_1 ~ c(
      log(1 + 0.122^2),
      0.37 * sqrt(log(1 + 0.122^2) * log(1 + 0.392^2)), log(1 + 0.392^2)
    ) # Table 2, IOV F 12.2%, IOV MTT 39.2%, correlation 0.37
    etaiov_lfdepot_2 + etaiov_lmtt_2 ~ c(
      fixed(log(1 + 0.122^2)),
      fixed(0.37 * sqrt(log(1 + 0.122^2) * log(1 + 0.392^2))), fixed(log(1 + 0.392^2))
    ) # ESM 2 BLOCK(2) SAME, occasion 2
    etaiov_lfdepot_3 + etaiov_lmtt_3 ~ c(
      fixed(log(1 + 0.122^2)),
      fixed(0.37 * sqrt(log(1 + 0.122^2) * log(1 + 0.392^2))), fixed(log(1 + 0.392^2))
    ) # ESM 2 BLOCK(2) SAME, occasion 3
    etaiov_lfdepot_4 + etaiov_lmtt_4 ~ c(
      fixed(log(1 + 0.122^2)),
      fixed(0.37 * sqrt(log(1 + 0.122^2) * log(1 + 0.392^2))), fixed(log(1 + 0.392^2))
    ) # ESM 2 BLOCK(2) SAME, occasion 4

    propSd <- 0.08
    label("Proportional residual error (fraction)") # Table 2, Proportional error = 8%
    addSd <- 1.7
    label("Additive residual error (mg/L)") # Table 2, Additive error = 1.7 mg/L
  })

  model({
    # Occasion indicators for the IOV etas (ESM 2: Q1-Q4 from OCC).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_lfdepot_1 + oc2 * etaiov_lfdepot_2 +
      oc3 * etaiov_lfdepot_3 + oc4 * etaiov_lfdepot_4
    iov_mtt <- oc1 * etaiov_lmtt_1 + oc2 * etaiov_lmtt_2 +
      oc3 * etaiov_lmtt_3 + oc4 * etaiov_lmtt_4

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    kel <- cl / vc

    # ESM 2: MTT = TVMTT * (1 + THETA(4)) * EXP(IOV * (1 + THETA(5))) for a
    # suspension; the IOV deviation is scaled, not the variance.
    mtt <- exp(lmtt + iov_mtt * (1 + e_susp_iovmtt * FORM_SUSPENSION)) *
      (1 + e_susp_mtt * FORM_SUSPENSION)

    # ESM 2 uses the Savic analytical transit input with NN = 4 into a depot
    # absorbed at KA = KTR = (NN + 1) / MTT. The explicit chain below is the
    # same system: five first-order stages (depot, transit1-4) form the
    # Erlang(NN + 1) input, and transit5 is the absorption compartment.
    ktr <- (4 + 1) / mtt
    ka <- ktr

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ka * transit5
    d/dt(central) <- ka * transit5 - kel * central

    f(depot) <- exp(lfdepot + etalfdepot + iov_fdepot)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
