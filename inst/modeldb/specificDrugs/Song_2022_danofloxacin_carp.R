Song_2022_danofloxacin_carp <- function() {
  description <- paste(
    "Veterinary (Yellow River carp, Cyprinus carpio haematopterus).",
    "One-compartment population PK model with first-order absorption and no",
    "absorption lag for danofloxacin after a single 10 mg/kg body-weight oral",
    "(gavage) dose in Yellow River carp held at 23.7 +/- 1.2 degC. Song 2022",
    "fitted 24 sparsely sampled fish (four plasma samples each, drawn from a",
    "16-point 0.25-144 h grid) in Phoenix NLME 8.1 with FOCE-ELS. The dose,",
    "V/F and CL/F are body-weight normalised in the source (dose in mg/kg,",
    "tvV in L/kg, tvCL in L/h/kg), so the model is driven with amt in mg/kg",
    "and returns Cc in ug/mL. Final model (Model 3, Table 5): tvKa 2.48 1/h,",
    "tvV 47.8 L/kg, tvCL 0.694 L/h/kg, exponential between-fish variability on",
    "all three (omega^2 = 3.45, 1.19 and 0.712) and a multiplicative residual",
    "error of 0.339. Model 3 keeps body weight as a power covariate on Ka, V/F",
    "and CL/F, centred on 0.38 kg (Eq. 1), but the three exponents are 'results",
    "not shown' in the paper; they are carried here as placeholders held at 0,",
    "so the model reproduces the published typical fish at every weight. The",
    "authors also fitted correlations between the three etas but did not print",
    "them, so the etas are independent here. See the vignette for both gaps.",
    sep = " "
  )
  reference <- paste(
    "Song Z-W, Yang F, Dai Y, Zhang C-S, Shao H-T, Wang H, Ma K-L, Li Z-E,",
    "Yang F. (2022). Population Pharmacokinetics of Danofloxacin in Yellow",
    "River Carp (Cyprinus carpio haematopterus) After One Single Oral Dose.",
    "Frontiers in Veterinary Science 9:868966. doi:10.3389/fvets.2022.868966.",
    sep = " "
  )
  vignette <- "Song_2022_danofloxacin_carp"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Both states hold danofloxacin (base); the amounts are body-weight
  # normalised (mg/kg) because Song 2022 reports the dose in mg/kg BW and V/F,
  # CL/F per kg BW. verified = TRUE: the analyte is danofloxacin throughout
  # (Drug Administration and Sampling, oral gavage) and the sampled matrix is
  # plasma (tail-vein blood centrifuged to plasma).
  compartmentData <- list(
    depot = list(analyte = "danofloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "danofloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight of the individual fish, recorded before dosing (Table 1 lists every fish).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Song 2022 Eq. 1: Ka_i = tvKa * (wt_i/0.38)^dKadwt * exp(eta_Ka_i), with V and CL 'calculated by the same algorithm'; 0.38 kg is the mean weight of the dosed fish. The three exponents are not reported, so e_wt_ka, e_wt_vc and e_wt_cl are placeholders held at 0 and WT has no effect until a user supplies values. Independently of this power term, the dose (mg/kg) and V/F, CL/F (per kg) are already body-weight normalised.",
      source_name = "wt"
    )
  )

  population <- list(
    species = "Yellow River carp (Cyprinus carpio haematopterus)",
    n_subjects = 24L,
    n_studies = 1L,
    weight_range = "0.25-0.54 kg",
    weight_mean = "0.38 kg",
    dose_range = "single 10 mg/kg body weight oral gavage (danofloxacin mesylate, 5 mg/mL in saline)",
    disease_state = "healthy",
    regions = "China (Henan University of Science and Technology, Luoyang)",
    water_temperature = "23.7 +/- 1.2 degC",
    water_quality = "pH ~7.3; dissolved oxygen > 8 mg/L; ammonia ~0.1 mg/L",
    design = "30 fish; 24 dosed fish in four equal tanks and 6 undosed controls supplying blank plasma. At least 7 days of acclimation before dosing (Animals section).",
    sampling = "Sparse: each dosed fish sampled 4 times from the tail vein. Table 1 cycles four schedules across the 24 fish: (0.25, 4, 16, 72), (0.5, 6, 24, 96), (1, 8, 36, 120) and (2, 12, 48, 144) h post-dose.",
    bioanalysis = "HPLC with fluorescence detection (ex 280 nm, em 450 nm); calibration 0.005-2 ug/mL, LOQ 0.005 ug/mL, mean recovery 94.03%, intra-day CV 0.36-7.36%, inter-day CV 4.74-7.97% (Results, Table 3).",
    notes = "Per-fish weights are in Table 1. Only one water temperature was studied, and the Discussion notes that fluoroquinolone elimination in fish depends on water temperature, so the parameters should not be extrapolated to other rearing temperatures."
  )

  ini({
    # =================================================================
    # Structural parameters -- Song 2022 Table 5, Model 3 ('Covariates
    # body weight and covariance model'), the final model. All three are
    # apparent (per bioavailability) values and are body-weight
    # normalised: V/F is L/kg and CL/F is L/h/kg, matched to a dose in
    # mg/kg BW. No absorption lag was fitted.
    # =================================================================
    lka <- log(2.48)
    label("Absorption rate constant Ka at 0.38 kg (1/h)") # Song 2022 Table 5, Model 3: tvKa = 2.48 1/h (SE 0.123, 95% CI 2.24-2.73)
    lvc <- log(47.8)
    label("Apparent volume of distribution V/F at 0.38 kg (L/kg)") # Song 2022 Table 5, Model 3: tvV = 47.8 L/kg (SE 1.90, 95% CI 44.0-51.6)
    lcl <- log(0.694)
    label("Apparent clearance CL/F at 0.38 kg (L/h/kg)") # Song 2022 Table 5, Model 3: tvCL = 0.694 L/h/kg (SE 0.012, 95% CI 0.669-0.718)

    # =================================================================
    # Body-weight power exponents -- Song 2022 Eq. 1, centred on the
    # 0.38 kg mean weight. The paper describes dKadwt as 'results not
    # shown' and computes V and CL 'by the same algorithm'; no exponent
    # is printed anywhere in the paper and there is no supplement. They
    # are placeholders held at 0 so the published covariate structure is
    # kept and a user with the authors' values can override them.
    # =================================================================
    e_wt_ka <- fixed(0)
    label("Power exponent of body weight on Ka; unreported in source, placeholder (unitless)") # Song 2022 Eq. 1 'dKadwt', value not shown
    e_wt_vc <- fixed(0)
    label("Power exponent of body weight on V/F; unreported in source, placeholder (unitless)") # Song 2022 Eq. 1 text, V 'by the same algorithm', value not shown
    e_wt_cl <- fixed(0)
    label("Power exponent of body weight on CL/F; unreported in source, placeholder (unitless)") # Song 2022 Eq. 1 text, CL 'by the same algorithm', value not shown

    # =================================================================
    # Between-fish variability -- Song 2022 Table 5, Model 3, column
    # 'Interindividual variation omega^2 (CV%)'. The leading number is the
    # variance omega^2; the bracketed number is its relative standard
    # error, not a between-fish CV (see the vignette). Model 3 added a
    # covariance between the etas, but no off-diagonal is printed, so the
    # etas are independent here.
    # =================================================================
    etalka ~ 3.45 # Song 2022 Table 5, Model 3 omega^2 for eta Ka = 3.45 (bracketed 0.203)
    etalvc ~ 1.19 # Song 2022 Table 5, Model 3 omega^2 for eta V = 1.19 (bracketed 8.40)
    etalcl ~ 0.712 # Song 2022 Table 5, Model 3 omega^2 for eta CL = 0.712 (bracketed 4.35)

    # =================================================================
    # Residual error -- Song 2022 Population Pharmacokinetic Modeling:
    # 'Multiplicative was selected as the residual error model'. Phoenix
    # NLME's multiplicative model is C * (1 + eps), eps ~ N(0, stdev0^2),
    # i.e. a proportional error, so stdev0 is a dimensionless fraction;
    # the 'ug/ml' in Table 5's Units column is not used.
    # =================================================================
    propSd <- 0.339
    label("Proportional residual error (fraction)") # Song 2022 Table 5, Model 3: residual error = 0.339 (SE 0.014, 95% CI 0.312-0.367)
  })

  model({
    # Song 2022 Eq. 1: P_i = tvP * (wt_i/0.38)^dPdwt * exp(eta_P_i) for
    # P in Ka, V, CL.
    ka <- exp(lka + e_wt_ka * log(WT / 0.38) + etalka)
    vc <- exp(lvc + e_wt_vc * log(WT / 0.38) + etalvc)
    cl <- exp(lcl + e_wt_cl * log(WT / 0.38) + etalcl)

    kel <- cl / vc

    # One compartment with first-order oral absorption and first-order
    # elimination; V and CL are apparent (already divided by F).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Doses enter as mg/kg body weight and vc is L/kg, so central/vc is
    # mg/L == ug/mL, the paper's concentration unit.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
