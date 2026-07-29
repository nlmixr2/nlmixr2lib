Chirehwa_2017_pyrazinamide <- function() {
  description <- "One-compartment population PK model with Savic-style transit-compartment absorption (NN = 28) for oral pyrazinamide in HIV/TB-coinfected adults on the WHO four-drug fixed-dose combination (Chirehwa 2017); fat-free mass (Janmahasatian formula) drives fixed allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1.0) referenced to a 42 kg subject, and CL/F increases linearly by 14.3% from day 1 to day 29 of treatment, attributed to rifampin-mediated enzyme induction."
  reference <- paste(
    "Chirehwa MT, McIlleron H, Rustomjee R, Mthiyane T, Onyebujoh P, Smith P, Denti P.",
    "Pharmacokinetics of Pyrazinamide and Optimal Dosing Regimens for Drug-Sensitive",
    "and -Resistant Tuberculosis.",
    "Antimicrob Agents Chemother. 2017;61(8):e00490-17.",
    "doi:10.1128/AAC.00490-17.",
    sep = " "
  )
  vignette <- "Chirehwa_2017_pyrazinamide"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass derived from body weight, height, and sex via the Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065) formula.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Drives fixed allometric scaling of CL/F (exponent 0.75)",
        "and V/F (exponent 1.0), both referenced to the cohort-median FFM of 42 kg",
        "(Chirehwa 2017 Table 1 footnote a: 'FFM (kg) 42.2 (28.0-57.6)' and Table 2",
        "footnote b: 'the values reported refer to a subject with FFM of 42 kg ...",
        "CLi = THETA_CL * (FFMi / 42)^0.75 ... Vi = THETA_V * (FFMi / 42)^1').",
        "FFM is computed in the vignette derivation chunk via the Janmahasatian",
        "WHSmax / WHS50 formula: men WHSmax = 42.92 kg/m^2, WHS50 = 30.93 kg/m^2;",
        "women WHSmax = 37.99 kg/m^2, WHS50 = 35.98 kg/m^2;",
        "FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT) with HT in metres and WT in kg",
        "(Chirehwa 2017 Methods 'Model development' paragraph 4).",
        "User data sets supply FFM directly as a covariate column."
      ),
      source_name        = "FFM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 61L,
    n_studies      = 1L,
    n_observations = 1342L,
    age_range      = "18-47 years",
    age_median     = "32 years",
    weight_range   = "34.4-98.7 kg",
    weight_median  = "55.2 kg",
    height_range   = "1.41-1.81 m",
    height_median  = "1.59 m",
    ffm_range      = "28.0-57.6 kg",
    ffm_median     = "42.2 kg",
    sex_female_pct = 54,
    disease_state  = "Adults coinfected with HIV and pulmonary tuberculosis (drug-susceptible TB at trial entry); 67% receiving antiretroviral therapy at PK sampling.",
    dose_range     = paste(
      "Oral pyrazinamide given as part of the WHO four-drug fixed-dose combination tablet",
      "(150 mg rifampin + 75 mg isoniazid + 400 mg pyrazinamide + 275 mg ethambutol per tablet),",
      "with the number of tablets per dose adjusted to body weight band per WHO guidelines.",
      "84% (51/61) of subjects were dosed Monday-Friday only; the remaining 16% were dosed every day.",
      "PK sampling on days 1, 8, 15, and 29 of TB treatment after an overnight fast,",
      "at predose and 1, 2, 4, 6, 8, and 12 h postdose, plus a 12-h pre-dose sample on day 15."
    ),
    regions        = "South Africa.",
    notes          = paste(
      "Baseline demographics from Chirehwa 2017 Table 1; final structural and variability",
      "estimates from Table 2. Lower limit of quantification was 0.2 mg/liter; one BLQ sample",
      "was discarded. Concomitant medications included rifampin (a potent CYP / UGT inducer)",
      "as part of the fixed-dose combination tablet, which the paper proposes as the mechanism",
      "behind the 14.3% increase in CL/F observed from day 1 to day 29 of treatment."
    )
  )

  ini({
    # =========================================================================
    # Structural PK parameters -- Chirehwa 2017 Table 2 'Estimated parameter
    # values from the final model'. CL/F and V/F are reported as the typical
    # values for a subject with FFM = 42 kg (Table 2 footnote b).
    #
    # Pyrazinamide PK was described by a one-compartment open model with
    # transit-compartment absorption (Savic 2007 parameterisation) and
    # first-order elimination ('Structural model and parameter estimates'
    # paragraph 1; 'Model development' paragraph 2).
    # =========================================================================
    lka  <- log(3.54);   label("Absorption rate constant (1/h)")                                       # Chirehwa 2017 Table 2 ka = 3.54 (95% CI 3.0-4.27)
    lcl  <- log(3.35);   label("Apparent oral clearance CL/F at day 1 of treatment, FFM 42 kg (L/h)")  # Chirehwa 2017 Table 2 CL/F day 1 = 3.35 (95% CI 3.11-3.56)
    lvc  <- log(43.2);   label("Apparent central volume of distribution V/F, FFM 42 kg (L)")           # Chirehwa 2017 Table 2 V/F = 43.2 (95% CI 41.5-44.7)
    lmtt <- log(0.542);  label("Mean transit time MTT through the absorption chain (h)")               # Chirehwa 2017 Table 2 MTT = 0.542 (95% CI 0.47-0.61)
    lnn  <- log(28);     label("Number of Savic-style transit compartments NN (continuous, unitless)") # Chirehwa 2017 Table 2 NN = 28 (95% CI 7-52)
    lfdepot <- fixed(log(1));  label("Bioavailability (fixed at 1)")                                   # Chirehwa 2017 Table 2 F = 1 fixed

    # =========================================================================
    # Allometric scaling on fat-free mass with exponents FIXED to the literature
    # values 0.75 (CL/F) and 1.0 (V/F). Chirehwa 2017 Results 'Structural model
    # and parameter estimates' paragraph 4: 'Estimating allometry exponents did
    # not result in a significant improvement in OFV and the estimated values
    # were close, 0.75 for CL and 1 for V, hence the exponents were fixed to
    # these literature values.'
    # =========================================================================
    e_ffm_cl <- fixed(0.75);  label("Allometric exponent of FFM on CL/F (unitless, fixed)")  # Chirehwa 2017 Methods/Results: fixed at the Anderson-Holford literature value
    e_ffm_vc <- fixed(1.0);   label("Allometric exponent of FFM on V/F (unitless, fixed)")   # Chirehwa 2017 Methods/Results: fixed at the Anderson-Holford literature value

    # =========================================================================
    # Linear increase in CL/F across the first month of TB treatment, attributed
    # in the paper's Discussion to rifampin-driven induction of pyrazinamide
    # metabolic enzymes. The relationship is CL(day) = CL_day1 * (1 + dcl_day29_frac
    # * day / 28) where 'day' goes from 0 (= day 1 of treatment) to 28 (= day 29
    # of treatment); this matches Chirehwa 2017 Methods 'Model development'
    # paragraph 3 (equation 2). The fractional increase is reported in Table 2
    # as 'Delta CL/F day 29 (%) = 14.3'.
    # =========================================================================
    dcl_day29 <- 14.3;  label("Linear fractional increase in CL/F from day 1 to day 29 of TB treatment (% of day-1 value)")  # Chirehwa 2017 Table 2 dCL/F day 29 = 14.3% (95% CI 6.0-25.8)

    # =========================================================================
    # IIV (Chirehwa 2017 Table 2 'Between-subject variability' and
    # 'Between-occasion variability' blocks; CV% on the SD scale).
    # nlmixr2lib convention (Bienczak_2016_nevirapine.R, Svensson_2018_
    # bedaquiline.R, Svensson_2016_rifampicin.R): BOV is dropped where a BSV
    # term is reported on the same parameter, and BOV is folded in as a
    # BSV-equivalent where only BOV is reported. omega^2 = log(1 + CV^2).
    #   - CL:  BSV  16.3% kept; BOV 13.3% dropped (see vignette Assumptions)
    #   - F:   BSV  10.7% kept; BOV 11.9% dropped (see vignette Assumptions)
    #   - ka:  no BSV reported; BOV 84.0% folded in as BSV-equivalent
    #   - MTT: no BSV reported; BOV 52.9% folded in as BSV-equivalent
    # =========================================================================
    etalcl     ~ 0.02622  # Chirehwa 2017 Table 2 BSV CL = 16.3% (95% CI 11.0-20.0); omega^2 = log(1 + 0.163^2) = 0.02622
    etalfdepot ~ 0.01138  # Chirehwa 2017 Table 2 BSV F  = 10.7% (95% CI  7.6-13.2); omega^2 = log(1 + 0.107^2) = 0.01138
    etalka     ~ 0.53392  # Chirehwa 2017 Table 2 BOV ka = 84.0% (95% CI 79.7-97.5) folded in as BSV-equivalent; omega^2 = log(1 + 0.84^2) = 0.53392
    etalmtt    ~ 0.24674  # Chirehwa 2017 Table 2 BOV MTT = 52.9% (95% CI 40.1-68.7) folded in as BSV-equivalent; omega^2 = log(1 + 0.529^2) = 0.24674

    # =========================================================================
    # Residual error (Chirehwa 2017 Table 2 'Error' block). Combined additive
    # + proportional structure.
    # =========================================================================
    addSd  <- 1.23;    label("Additive residual SD (mg/L)")                # Chirehwa 2017 Table 2 Additive error = 1.23 mg/liter (95% CI 0.85-1.57)
    propSd <- 0.044;   label("Proportional residual SD (fraction)")        # Chirehwa 2017 Table 2 Coefficient of variation = 4.4% (95% CI 2.8-5.4); 0.044 as a fraction
  })

  model({
    # ----- 1. Individual PK parameters with FFM allometric scaling. -----
    # Reference subject: FFM = 42 kg, day 0 of treatment (= treatment day 1).
    ka    <- exp(lka  + etalka)
    mtt   <- exp(lmtt + etalmtt)
    nn    <- exp(lnn)
    vc    <- exp(lvc) * (FFM / 42)^e_ffm_vc

    # ----- 2. Linear time-on-treatment effect on CL/F. -----
    # Chirehwa 2017 Methods 'Model development' paragraph 3 (equation 2):
    #   CL(day) = CL_day1 * (1 + dcl_day29 / 100 * day / 28)
    # where 'day' is days elapsed from the start of treatment, ranging from 0
    # (the day-1-of-treatment sampling occasion) to 28 (the day-29 sampling
    # occasion). The rxode2 simulation-time variable `t` is in hours; days
    # elapsed = t / 24. The relationship is only validated within the fitted
    # window 0-28 days; predictions beyond day 29 should be interpreted with
    # caution (Chirehwa 2017 Discussion paragraph 5). Users who wish to
    # disable induction (e.g., the paper's MDR-TB simulations, which used
    # day-1 PK because second-line regimens lack rifampin) can set
    # `dcl_day29 = 0` in the `ini()` block before simulating; see the
    # vignette for both day-1 and day-29 simulation recipes.
    day_on_tx     <- t / 24
    cl_time_factor <- 1 + (dcl_day29 / 100) * (day_on_tx / 28)
    cl            <- exp(lcl + etalcl) * (FFM / 42)^e_ffm_cl * cl_time_factor

    # ----- 3. Micro-constants. -----
    ktr <- (nn + 1) / mtt
    kel <- cl / vc

    # ----- 4. One-compartment PK with Savic 2007 analytical transit-chain
    #         input rate. rxode2's built-in `transit(nn, mtt)` function returns
    #         the closed-form gamma-PDF mass-flow rate into depot; setting
    #         bioavailability on depot to 0 (after the lfdepot anchor below)
    #         suppresses the bolus content so the transit-chain rate is the
    #         only input. The standard `f(depot) <- exp(lfdepot)` anchor
    #         is overridden by an explicit f(depot) <- 0 to match the Savic
    #         2007 / Barnett 2018 convention.
    d/dt(depot)   <- transit(nn, mtt) - ka * depot
    d/dt(central) <-                     ka * depot - kel * central

    # ----- 5. Bioavailability. -----
    # Pyrazinamide F is anchored at 1 (Chirehwa 2017 Table 2 'F = 1 fixed').
    # BSV on F enters multiplicatively via etalfdepot. The transit() input
    # rate already integrates the dose mass into depot, so the bolus content
    # is set to 0 to avoid double-counting.
    f(depot) <- 0
    fdepot   <- exp(lfdepot + etalfdepot)

    # ----- 6. Observation: plasma pyrazinamide concentration (mg/L). -----
    # Dose in mg, V/F in L -> central / vc has units of mg/L (= ug/mL), matching
    # the units declared above and Chirehwa 2017's reported assay (lower limit
    # of quantification 0.2 mg/liter; concentrations expressed in mg/liter).
    # The bioavailability multiplier fdepot acts on the observation because
    # f(depot) was set to 0 to suppress the bolus and let transit() drive
    # the input rate; multiplying central / vc by fdepot recovers the
    # F-modulated dose-normalised concentration.
    Cc <- fdepot * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
