Janssen_2025_isoniazid <- function() {
  description <- "Two-compartment population PK model for oral isoniazid in South African adults on once-weekly rifapentine/isoniazid (3HP) tuberculosis preventive therapy (Janssen 2025). A two-class mixture separates fast/intermediate from slow acetylators on CL/F (supplied here as NAT2_SLOW); CL, Q, V and Vp scale allometrically with fat-free mass (reference 45 kg). Absorption is a Savic transit chain (5 transit compartments, absorption rate equal to the transit rate) whose mean transit time is shortened when the tablets are suspended in water; IOV on bioavailability and MTT, with smaller MTT IOV after suspended tablets."
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
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = "n/a; power scaling (FFM / 45)^0.75 on CL/F and Q/F, (FFM / 45)^1 on V/F and Vp/F",
      notes = "Calculated with the Janmahasatian 2005 formula (Janssen 2025 Methods, reference 12). Cohort median 38.6 kg (IQR 34.5-42.8), Table 1.",
      source_name = "FFM"
    ),
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator indicator: 1 = slow acetylator, 0 = fast or intermediate acetylator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fast/intermediate acetylator; CL/F 23.3 L/h at FFM 45 kg)",
      notes = "No NAT2 genotype was available, so the paper assigned the acetylator class with a two-subpopulation NONMEM $MIX model (ESM 3: MIXNUM 1 = fast/intermediate with probability 0.659, MIXNUM 2 = slow). The model could not separate fast from intermediate acetylators. Supply a known phenotype here, or draw NAT2_SLOW ~ Bernoulli(1 - 0.659) per subject to reproduce the fitted population.",
      source_name = "MIXNUM (NONMEM mixture class); NAT2_SLOW = as.integer(MIXNUM == 2)"
    ),
    FORM_SUSPENSION = list(
      description = "Suspended-tablet indicator: 1 = tablets quartered, crushed and suspended in 30 mL distilled water immediately before dosing (container rinsed with a further 70 mL); 0 = tablets swallowed whole with 100 mL water",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablets swallowed whole; the typical-value MTT reference in Janssen 2025 Table 2)",
      notes = "Per-dose (per-occasion) indicator: each participant received both forms in a crossover. Shortens the typical MTT by 35% and scales the IOV on MTT by (1 - 0.629). The three products (two fixed-dose-combination brands and standalone isoniazid) share one set of parameters in the final model; brand effects on bioavailability were assessed by log-likelihood profiling only.",
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
    n_observations = 707L,
    age_range = "18-55 years (inclusion); median 37.5 (IQR 26.5-45.3)",
    weight_range = "40-90 kg (inclusion); median 61.1 kg (IQR 51.1-66)",
    sex_female_pct = 83.3,
    race_ethnicity = c(MixedRace = 91.7, Black = 8.3),
    disease_state = "Adults with a positive tuberculosis infection test and recent household contact with drug-susceptible pulmonary TB, on 3HP tuberculosis preventive therapy; HIV and chronic hepatitis B excluded",
    dose_range = "Isoniazid 900 mg once weekly for 12 weeks (with rifapentine 900 mg), as whole or water-suspended tablets of two fixed-dose combinations (rifapentine 300 mg / isoniazid 300 mg) or standalone 300 mg tablets; given 30 min after a standard 670 kcal breakfast",
    regions = "Cape Town, South Africa (single centre)",
    notes = "Open-label randomized crossover; each participant received two of three products, one whole and one suspended, across the four sampled weekly doses. Demographics from Janssen 2025 Table 1 (FFM median 38.6 kg, IQR 34.5-42.8). BLQ samples (18.2% within 48 h; all after 48 h) were omitted. Estimated fast/intermediate acetylator fraction 65.9%."
  )

  ini({
    lcl_nonslow <- log(23.3)
    label("Apparent clearance CL/F, fast/intermediate acetylators, FFM 45 kg (L/h)") # Table 2, CL fast/intermediate/F = 23.3 (90% CI 21.3-25.6)
    lcl_slow <- log(9.5)
    label("Apparent clearance CL/F, slow acetylators, FFM 45 kg (L/h)") # Table 2, CL slow/F = 9.5 (90% CI 8.5-10.6)
    lvc <- log(59.6)
    label("Apparent central volume V/F, FFM 45 kg (L)") # Table 2, V/F = 59.6 (90% CI 55.9-63.2)
    lq <- log(1.26)
    label("Apparent intercompartmental clearance Q/F, FFM 45 kg (L/h)") # Table 2, Q/F = 1.26 (90% CI 1.1-1.5)
    lvp <- log(56.7)
    label("Apparent peripheral volume Vp/F, FFM 45 kg (L)") # Table 2, Vp/F = 56.7 (90% CI 48.4-65.9)
    lmtt <- log(1.4)
    label("Mean transit time of the absorption chain, whole tablet (h)") # Table 2, MTT = 1.4 (90% CI 1.3-1.5)
    lfdepot <- fixed(log(1))
    label("Reference bioavailability F (unitless)") # Table 2, F = 1 (Fixed)

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of FFM on CL/F and Q/F (unitless)") # Table 2 footnote b, x = 0.75 for CL and Q
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of FFM on V/F and Vp/F (unitless)") # Table 2 footnote b, x = 1 for V and Vp

    e_susp_mtt <- -0.35
    label("Fractional change in MTT for suspended tablets (unitless)") # Table 2, Effect of suspension form on MTT = -35% (90% CI -40.3 to -29.3)
    e_susp_iovmtt <- -0.629
    label("Fractional change in the IOV of MTT for suspended tablets (unitless)") # Table 2, Effect of suspension form on MTT variation = -62.9% (90% CI -70.8 to -54.5)

    # IIV and IOV variances from Table 2 %CV via omega^2 = log(1 + CV^2)
    # (Table 2 footnote c). The CL eta is shared by both acetylator classes
    # (ESM 3: ETA(1) in both MIXNUM branches).
    etalcl ~ log(1 + 0.174^2) # Table 2, IIV on CL 17.4% CV
    etalfdepot ~ log(1 + 0.181^2) # Table 2, IIV on F 18.1% CV

    # IOV, ESM 3 $OMEGA BLOCK(1) with three BLOCK(1) SAME each for F and MTT.
    etaiov_lfdepot_1 ~ log(1 + 0.071^2) # Table 2, IOV on F 7.1% CV
    etaiov_lfdepot_2 ~ fixed(log(1 + 0.071^2)) # ESM 3 BLOCK(1) SAME, occasion 2
    etaiov_lfdepot_3 ~ fixed(log(1 + 0.071^2)) # ESM 3 BLOCK(1) SAME, occasion 3
    etaiov_lfdepot_4 ~ fixed(log(1 + 0.071^2)) # ESM 3 BLOCK(1) SAME, occasion 4
    etaiov_lmtt_1 ~ log(1 + 0.408^2) # Table 2, IOV on MTT 40.8% CV
    etaiov_lmtt_2 ~ fixed(log(1 + 0.408^2)) # ESM 3 BLOCK(1) SAME, occasion 2
    etaiov_lmtt_3 ~ fixed(log(1 + 0.408^2)) # ESM 3 BLOCK(1) SAME, occasion 3
    etaiov_lmtt_4 ~ fixed(log(1 + 0.408^2)) # ESM 3 BLOCK(1) SAME, occasion 4

    # ESM 3 fits on the log scale, Y = log(C) + SD * EPS with
    # SD = sqrt(THETA(8)^2 + (THETA(9) / C)^2); to first order this is a
    # combined proportional + additive error on the linear scale.
    propSd <- 0.208
    label("Proportional residual error (fraction)") # Table 2, Proportional error = 20.8%
    addSd <- 0.11
    label("Additive residual error (mg/L)") # Table 2, Additive error = 0.11 mg/L
  })

  model({
    # Occasion indicators for the IOV etas (ESM 3: Q1-Q4 from OCC).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_lfdepot_1 + oc2 * etaiov_lfdepot_2 +
      oc3 * etaiov_lfdepot_3 + oc4 * etaiov_lfdepot_4
    iov_mtt <- oc1 * etaiov_lmtt_1 + oc2 * etaiov_lmtt_2 +
      oc3 * etaiov_lmtt_3 + oc4 * etaiov_lmtt_4

    # Acetylator class selects the typical CL/F (ESM 3 $MIX, MIXNUM 1/2).
    lcl_i <- lcl_nonslow * (1 - NAT2_SLOW) + lcl_slow * NAT2_SLOW

    fcl <- (FFM / 45)^e_ffm_cl
    fv <- (FFM / 45)^e_ffm_vc
    cl <- exp(lcl_i + etalcl) * fcl
    vc <- exp(lvc) * fv
    q <- exp(lq) * fcl
    vp <- exp(lvp) * fv
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ESM 3: MTT = TVMTT * (1 + THETA(10)) * EXP(IOV * (1 + THETA(11))) for a
    # suspension; the IOV deviation is scaled, not the variance.
    mtt <- exp(lmtt + iov_mtt * (1 + e_susp_iovmtt * FORM_SUSPENSION)) *
      (1 + e_susp_mtt * FORM_SUSPENSION)

    # ESM 3 uses the Savic analytical transit input with NN = 5 into a depot
    # absorbed at KA = KTR = (NN + 1) / MTT. The explicit chain below is the
    # same system: six first-order stages (depot, transit1-5) form the
    # Erlang(NN + 1) input, and transit6 is the absorption compartment.
    ktr <- (5 + 1) / mtt
    ka <- ktr

    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot + etalfdepot + iov_fdepot)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
