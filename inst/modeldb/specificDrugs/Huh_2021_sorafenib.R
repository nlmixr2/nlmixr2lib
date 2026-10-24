Huh_2021_sorafenib <- function() {
  description <- "One-compartment population PK model for single oral doses of sorafenib in healthy Korean men, comparing the nanoparticulated formulation SYO-1644 (100, 150 and 200 mg) with the reference Nexavar 200 mg tablet in a replicated two-period study (Huh 2021). Absorption is a Savic transit-compartment chain (mean transit time MTT, non-integer number of transit compartments NN) feeding an absorption compartment that empties into central at first-order ka; SYO-1644 has a relative bioavailability of 1.62 versus Nexavar. Enterohepatic recirculation follows Jain 2011: a fraction Fent of the drug leaving central is routed to a gallbladder compartment and the rest, (1 - Fent) * CL, is eliminated; the gallbladder releases drug back to the absorption compartment at kEhc gated by a steep Hill-40 switch of time after dose that opens at an early emptying time Tgb1, resets 3 h later and reopens at a late emptying time Tgb2. The central-to-gallbladder rate constant kgb is not reported in the paper; its value here was back-solved by the maintainers from the paper's own Figure 2 Nexavar mean profile. Inter-individual variability on ka, CL and V, inter-occasion variability on CL across the two study periods, and proportional residual error."
  reference <- paste(
    "Huh KY, Hwang S, Park SY, Lim HJ, Jin M, Oh J, Yu KS, Chung JY.",
    "Population Pharmacokinetic Modelling and Simulation to Determine the",
    "Optimal Dose of Nanoparticulated Sorafenib to the Reference Sorafenib.",
    "Pharmaceutics. 2021;13(5):629.",
    "doi:10.3390/pharmaceutics13050629.",
    sep = " "
  )
  vignette <- "Huh_2021_sorafenib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Oral doses enter only the absorption compartment, through transit().
  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder = list(analyte = "sorafenib", units = "mg", specimen = "bile", verified = TRUE)
  )

  covariateData <- list(
    FORM_SORAFENIB_SYO1644 = list(
      description = "Sorafenib formulation indicator: 1 = SYO-1644 (nanoparticulated sorafenib), 0 = Nexavar (reference sorafenib tablet)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies the transit-chain input by the relative bioavailability 1.62 (Table 2, 'Relative bioavailability (SYO-1644 to Nexavar)'). Per-dose-record indicator; constant within a subject in Huh 2021 because each subject received the same assigned treatment in both periods.",
      source_name = "treatment group (SYO-1644 vs Nexavar)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability on CL: 1 = period 1, 2 = period 2 (repeat dose after a two-week washout)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Huh 2021 Methods 2.5: 'The IOV was modeled using the OMEGA BLOCK step in NONMEM'. Two occasions (one per period of the replicated design) with a shared variance, the $OMEGA BLOCK(1) SAME idiom. Records outside the two periods may carry OCC = 0, which zeroes both indicators and removes the IOV contribution.",
      source_name = "period"
    )
  )

  covariatesDataExcluded <- list(
    meal_times = list(
      description = "Lunch and dinner times on the dosing day",
      units = "h",
      type = "continuous",
      notes = "Methods 2.5: 'Meal times (i.e., lunch and dinner) at the dosing date were also included as parameters.' Meals were served about 4 and 10 h after dosing (Methods 2.2). The final model carries no meal-time covariate; the two gallbladder emptying times Tgb1 and Tgb2 are estimated instead and reflect 'the relatively regular meal times in our study' (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "19-50 years (inclusion criterion); mean 30.9 (SD 7.7) years",
    weight_range = "mean 71.0 (SD 8.2) kg; BMI 18.0-27.0 kg/m^2 (inclusion criterion), mean 23.5 (SD 2.4)",
    sex_female_pct = 0,
    race_ethnicity = "Korean (single-centre study at Seoul National University Hospital)",
    disease_state = "Healthy volunteers",
    dose_range = "Single oral doses under fasting of SYO-1644 100, 150 or 200 mg or Nexavar 200 mg, given in period 1 and repeated in period 2 after a two-week washout",
    regions = "South Korea",
    notes = "Randomized, open-label, replicated two-period study (NCT03674060). 32 healthy men enrolled and all included in the population PK analysis; 30 completed both periods (one subject each in the SYO-1644 100 and 150 mg arms withdrew after period 1). Arms: SYO-1644 100 mg (n = 7 completers), 150 mg (n = 7), 200 mg (n = 8), Nexavar 200 mg (n = 8). PK sampling to 168 h post-dose; meals about 4 and 10 h post-dose. Demographics in Results 3.1."
  )

  ini({
    # Absorption: Savic transit chain (Table 2)
    lmtt <- log(3.24); label("Mean transit time of the absorption chain MTT (h)")              # Table 2: MTT = 3.24 h (RSE 5.4%)
    lnn <- log(1.72); label("Number of transit compartments NN (unitless)")                     # Table 2: NN = 1.72 (RSE 7.3%)
    lka <- log(1.6); label("Absorption rate constant ka (1/h)")                                 # Table 2: ka = 1.6 1/h (RSE 16.8%)

    # Disposition (Table 2). Apparent values: sorafenib is oral-only and the
    # Nexavar formulation is the bioavailability reference.
    lcl <- log(3.33); label("Apparent clearance CL/F (L/h)")                                     # Table 2: CL = 3.33 L/h (RSE 9.6%)
    lvc <- log(8.17); label("Apparent central volume V/F (L)")                                   # Table 2: V = 8.17 L (RSE 31.5%)

    # Enterohepatic recirculation (Table 2, Figure 3, Equation 3)
    fent <- 0.215; label("Fraction of central outflow routed to the gallbladder Fent (unitless)") # Table 2: Fent = 0.215 (RSE 13.5%)
    lkbm <- log(7.7); label("Central-to-gallbladder transfer rate constant kgb (1/h)")          # NOT in the paper (Figure 3 arrow 'Fent kgb' carries no value); back-solved by the maintainers to the digitised Figure 2 Nexavar 200 mg mean profile (population-mean fit, all other parameters at Table 2), see vignette
    lkehc <- log(0.212); label("Gallbladder-to-absorption-compartment release rate constant kEhc (1/h)") # Table 2: kEhc = 0.212 1/h (RSE 6.8%)
    ltgb1 <- log(8.14); label("Early gallbladder emptying time Tgb1 (h after dose)")            # Table 2: Tgb1 = 8.14 h (RSE 0.3%)
    ltgb2 <- log(12.2); label("Late gallbladder emptying time Tgb2 (h after dose)")             # Table 2: Tgb2 = 12.2 h (RSE 0.5%)
    hill <- fixed(40); label("Hill exponent of the gallbladder emptying switch Ehc (unitless)") # Results 3.3: 'The slope of the sigmoid function was set to 40'

    # Formulation effect on bioavailability (Table 2)
    e_form_syo1644_fdepot <- log(1.62); label("Log relative bioavailability of SYO-1644 versus Nexavar (unitless; exp = 1.62)") # Table 2: Relative bioavailability (SYO-1644 to Nexavar) = 1.62 (RSE 8.7%)

    # Inter-individual variability. Table 2 reports %CV; the Results text gives
    # the corresponding variances for V and CL as 4.5% and 6.6% (0.045, 0.066),
    # which match omega^2 = log(1 + CV^2) for 21.5% and 26.1%.
    etalka ~ 0.5636 # Table 2: IIV ka = 87.3 %CV; omega^2 = log(1 + 0.873^2)
    etalcl ~ 0.06589 # Table 2: IIV CL = 26.1 %CV; omega^2 = log(1 + 0.261^2) (text: 6.6%)
    etalvc ~ 0.04517 # Table 2: IIV V = 21.5 %CV; omega^2 = log(1 + 0.215^2) (text: 4.5%)

    # Inter-occasion variability on CL across the two study periods, with a
    # shared variance ($OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.008067 # Table 2: IOV CL = 9.0 %CV; omega^2 = log(1 + 0.090^2)
    etaiov_cl_2 ~ fixed(0.008067) # same magnitude as occasion 1 ($OMEGA BLOCK(1) SAME)

    # Residual error: Table 2 reports a proportional term only.
    propSd <- 0.376; label("Proportional residual error (fraction)")                           # Table 2: Proportional residual error = 37.6 %CV (RSE 4.4%)
  })
  model({
    # Occasion indicators; OCC = 0 zeroes the IOV contribution.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    mtt <- exp(lmtt)
    nn <- exp(lnn)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    kbm <- exp(lkbm)
    kehc <- exp(lkehc)
    tgb1 <- exp(ltgb1)
    tgb2 <- exp(ltgb2)
    kel <- cl / vc

    # Relative bioavailability, Nexavar = 1.
    fdepot <- exp(e_form_syo1644_fdepot * FORM_SORAFENIB_SYO1644)

    # Gallbladder emptying switch, Equation 3: a Hill-40 sigmoid of time after
    # dose centred on Tgb1 until Tgb1 + 3 h, when it resets and is centred on
    # Tgb2 thereafter.
    tad_h <- tad()
    if (tad_h < tgb1 + 3) {
      ehc <- tad_h^hill / (tad_h^hill + tgb1^hill)
    } else {
      ehc <- tad_h^hill / (tad_h^hill + tgb2^hill)
    }

    # Savic closed-form transit input (ktr = (NN + 1) / MTT inside transit())
    # into the absorption compartment, which empties into central at ka. The
    # gallbladder returns drug to the absorption compartment (Figure 3).
    # f(depot) <- 0 suppresses the ordinary bolus so the transit input is the
    # only dose pathway; transit() carries the bioavailability itself.
    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot + ehc * kehc * gallbladder
    d/dt(central) <- ka * depot - (1 - fent) * kel * central - fent * kbm * central
    d/dt(gallbladder) <- fent * kbm * central - ehc * kehc * gallbladder
    f(depot) <- 0

    # Dose in mg and volume in L give mg/L; x 1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
