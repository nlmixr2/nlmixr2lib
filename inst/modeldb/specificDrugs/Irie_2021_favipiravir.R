Irie_2021_favipiravir <- function() {
  description <- "One-compartment population PK model for oral favipiravir in hospitalized adults with COVID-19, with dose entered directly into the central compartment (NONMEM ADVAN1, no absorption phase). CL/F is a power function of the last administered dose (dose-dependent nonlinear PK), a multiplicative ratio for time-varying invasive mechanical ventilation, and a power function of body surface area."
  reference <- "Irie K, Nakagawa A, Fujita H, Tamura R, Eto M, Ikesue H, Muroi N, Fukushima S, Tomii K, Hashida T. Population pharmacokinetics of favipiravir in patients with COVID-19. CPT Pharmacometrics Syst Pharmacol. 2021;10(10):1161-1170. doi:10.1002/psp4.12685"
  vignette <- "Irie_2021_favipiravir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "favipiravir", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description = "Last administered favipiravir dose at the time of the record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (use case (a) of the DOSE register entry, but re-evaluated at every dose). Irie 2021 Methods, Covariate analysis: 'FPV dosage (as the last administered dosage at the time)' was incorporated as a time-varying covariate. The deposited control stream (Supplementary Data S2) carries it as data item LAST in micrograms, entering as TVCL = THETA(1) * (600000/LAST)**THETA(3) with THETA(3) = 0.61, i.e. (LAST/600 mg)^-0.61. Set DOSE on every record to the most recent dose given (1600 or 1800 mg on the loading day, 600 or 800 mg thereafter). rxode2 5.1.8 drops a DOSE column placed before 'amt' in the event data, so relocate event columns first.",
      source_name = "LAST"
    ),
    MECH_VENT = list(
      description = "Invasive mechanical ventilation indicator (WHO ordinal clinical status = 7)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on invasive mechanical ventilation)",
      notes = "TIME-VARYING, re-evaluated as ventilation status changes (not time-fixed at admission). Irie 2021 Methods: 'Clinical status ... changed over time, and these factors were incorporated as time-varying covariates'; Results: baseline IMV (BTUB in the control stream) was not significant (Table 2 step 11, p = 0.138) while time-varying IMV was (step 10, p < 0.001). Control-stream data item TUBE, entering as THETA(4)**TUBE. 10 of 39 patients were on IMV at the start of favipiravir, 7 were intubated later and 2 were weaned during treatment. IMV patients received tablet suspension via nasogastric tube.",
      source_name = "TUBE"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL/F normalized to the cohort median 1.72 m^2 (Table 1; control stream '(BSA/1.72)**THETA(5)'). BSA formula not stated in the source. Observed range 1.14-2.20 m^2.",
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 39L,
    n_studies = 1L,
    n_observations = 204L,
    age_range = "27-89 years",
    age_median = "68 years",
    weight_range = "29-100 kg",
    weight_median = "64 kg",
    height_median = "168 cm (range 144-182)",
    bsa_median = "1.72 m^2 (range 1.14-2.20)",
    sex_female_pct = 20.5,
    race_ethnicity = "Not reported (single Japanese centre)",
    disease_state = "Hospitalized COVID-19 (RT-PCR confirmed); WHO ordinal score at start 7 (invasive ventilation) n = 10, 5 (oxygen) n = 24, 4 (no oxygen) n = 5",
    dose_range = "Oral favipiravir tablets (or tablet suspension via nasogastric tube during IMV): 1600 mg b.i.d. on Day 1 then 600 mg b.i.d. (n = 33); 1800 mg b.i.d. on Day 1 then 800 mg b.i.d., switched to 600 mg b.i.d. (n = 6); 5 to 14 days",
    regions = "Japan (Kobe City Medical Center General Hospital)",
    notes = "Retrospective observational study, March-May 2020, using residual serum from routine laboratory samples (median 5 samples per patient, range 1-13; 15 of 219 samples below the 5 ng/mL LOQ excluded). Baseline characteristics per Irie 2021 Table 1."
  )

  ini({
    # Final-model estimates: Irie 2021 Table 4 and the deposited NONMEM control
    # stream (Supplementary Data S2, $THETA/$OMEGA/$SIGMA carry the same values).
    lcl <- log(5.11); label("Apparent clearance CL/F at 600 mg, no IMV, BSA 1.72 m^2 (L/h)") # Table 4 'CL/F' = 5.11 L/h; S2 THETA(1)
    lvc <- log(41.6); label("Apparent volume of distribution V/F (L)") # Table 4 'V/F' = 41.6 (column printed 'L/h', a typo); S2 THETA(2)

    e_dose_cl <- -0.61; label("Power exponent of last administered dose (DOSE/600 mg) on CL/F (unitless)") # Table 4 'Dose on CL/F' = -0.61; S2 THETA(3) = 0.61 in (600000/LAST)**THETA(3)
    e_mech_vent_cl <- 1.71; label("Multiplicative ratio on CL/F during invasive mechanical ventilation (unitless)") # Table 4 'IMV on CL/F' = 1.71; S2 THETA(4)**TUBE
    e_bsa_cl <- 2.22; label("Power exponent of BSA/1.72 on CL/F (unitless)") # Table 4 'BSA on CL/F' = 2.22; S2 THETA(5)

    etalcl ~ 0.355 # Table 4 'omega2 CL/F' = 0.355; S2 $OMEGA, CL = TVCL*EXP(ETA(1))

    propSd <- 0.8654479; label("Proportional residual error (fraction)") # Table 4 'sigma2 proportional error' = 0.749 -> sqrt(0.749)
    addSd <- 0.02764055; label("Additive residual error (ug/mL)") # Table 4 'sigma2 additive error' = 764 (ng/mL)^2 -> sqrt(764) = 27.64 ng/mL = 0.02764 ug/mL
  })

  model({
    # CL/F equation (Irie 2021 Results, Covariate analysis; S2 $PK):
    # CL/F = 5.11 * (Dose/600)^-0.61 * 1.71^IMV * (BSA/1.72)^2.22
    cl <- exp(lcl + etalcl) * (DOSE / 600)^e_dose_cl * e_mech_vent_cl^MECH_VENT * (BSA / 1.72)^e_bsa_cl
    vc <- exp(lvc)

    kel <- cl / vc

    # ADVAN1 TRANS2: the oral dose enters the central compartment directly
    # (no absorption phase was estimated).
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
