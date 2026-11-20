Goeyvaerts_2026_mosnodenvir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral mosnodenvir (JNJ-64281802,",
    "pan-serotype dengue NS3-NS4B interaction inhibitor) in healthy adults,",
    "pooled from three Phase 1 studies and Cohort 1 of the Phase 2a dengue",
    "controlled human infection study (Goeyvaerts 2026). The dose enters a",
    "depot that drains through six transit compartments at rate KTR before",
    "first-order absorption (Ka) into the central compartment, with linear",
    "elimination. Relative bioavailability, Ka and KTR depend on",
    "formulation, food and dose: fed tablet is the F = 1 reference; fasted",
    "solution and fasted tablet carry a dose-dependent F = 1 - DOSE/(400 +",
    "DOSE); fasted capsule F = 0.683; Ka for fasted solution falls with",
    "dose by a power law. BMI (reference 24.7 kg/m^2) is a power covariate",
    "on CL/F, F and Ka. Log-normal IIV on CL/F, Vc/F, F, Ka and KTR, and",
    "inter-occasion variability on F and KTR. This model supplies the",
    "individual steady-state average concentrations that drive the viral",
    "dynamic models Goeyvaerts_2026_mosnodenvir_vdm and",
    "Goeyvaerts_2026_mosnodenvir_vdm_c1."
  )
  reference <- paste(
    "Goeyvaerts N, Woot de Trixhe X, Neyens M, Kakuda TN, Van Wesenbeeck L,",
    "Buelens A, Rasschaert F, Durbin A, Ackaert O. Viral Dynamic",
    "Model-Informed Dose Selection of Mosnodenvir Prophylaxis against DENV-3",
    "in a Controlled Human Infection Model. Clin Pharmacol Ther.",
    "2026;120(5):1379-1388. doi:10.1002/cpt.70439"
  )
  vignette <- "Goeyvaerts_2026_mosnodenvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (BMI/24.7)^e on CL/F, F and Ka; 24.7 kg/m^2 is the median BMI of the PopPK analysis dataset (Supplement Table S2 footnote a).",
      source_name = "BMI"
    ),
    FED = list(
      description = "Fed state at dosing (1 = fed, standard breakfast; 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Only the tablet formulation was studied fed (Supplement Table S1). FED = 1 selects F = 1 (the reference), Ka_fed and KTR_fed; it is not defined for the capsule or solution, which were studied fasted only.",
      source_name = "food"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator (1 = capsule, 0 = tablet or oral solution)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet; FORM_CAPSULE = 0 and FORM_SOLUTION = 0)",
      notes = "Comparator is the tablet; capsule fasted F = 0.683 relative to fed tablet (F = 1), not dose dependent; capsule Ka is shared with the fasted tablet (Supplement Table S2). Capsules were given fasted only (Phase 2a Cohort 1, Phase 1 64281802DNG1001).",
      source_name = "formulation"
    ),
    FORM_SOLUTION = list(
      description = "Oral-solution formulation indicator (1 = oral solution, 0 = tablet or capsule)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet; FORM_CAPSULE = 0 and FORM_SOLUTION = 0)",
      notes = "Comparator is the tablet; the fasted solution shares the dose-dependent F of the fasted tablet and carries its own dose-dependent Ka (Supplement Table S2). Solution was given fasted only (Phase 1 64281802DNG1001 and DNG1006).",
      source_name = "formulation"
    ),
    DOSE_MOSNODENVIR_MG = list(
      description = "Nominal mosnodenvir dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the fasted solution / fasted tablet bioavailability 1 - DOSE/(ID50 + DOSE) (ID50 = 400 mg) and the fasted-solution Ka power term (DOSE/400)^-0.152 (Supplement Table S2 note). Set it to the amount of the dose record being absorbed; for fed-tablet or capsule dosing it has no effect.",
      source_name = "DOSE"
    ),
    OCC = list(
      description = "Occasion (study period) index for inter-occasion variability on F and KTR",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "IOV accounts for absorption variability across the periods of the crossover parts of 64281802DNG1006 (2-period crossover) and 64281802DNG1001 (Supplement, PopPK model). Two occasions are encoded, sharing one variance (NONMEM $OMEGA BLOCK SAME). OCC = 0 (the default for simulation) zeroes the IOV contribution, as in the paper's dose-selection simulations, which ignored IOV.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "mosnodenvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mosnodenvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mosnodenvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 246L,
    n_studies = 4L,
    n_observations = 8352L,
    age_range = "Adults (Phase 2a cohorts 18-53 years; Phase 1 ranges not reported in this paper)",
    weight_range = "Phase 2a Cohort 1: 48.2-107 kg (Table 1); Phase 1 not reported in this paper",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported for the pooled PopPK dataset; Phase 2a Cohort 1 was 48.3% Black, 31.0% White (Table 1).",
    disease_state = "Healthy adult volunteers (Phase 2a participants were dengue-seronegative and later challenged with DENV-3 rDEN3delta30).",
    dose_range = paste(
      "Single doses 50-1200 mg (oral solution, tablet fasted/fed, capsule);",
      "multiple doses 50-560 mg QD (solution), 150/450 mg BID loading then",
      "150-900 mg once or twice weekly (fed tablet), and 600/200, 200/50 or",
      "40/10 mg QD loading/maintenance (fasted capsule)."
    ),
    regions = "Not reported (Johnson & Johnson-sponsored studies; Phase 2a in the United States).",
    notes = paste(
      "Studies 64281802DNG1001, DNG1006, DNG1008 and DNG2002 Cohort 1",
      "(Supplement Table S1); 8,352 plasma samples, LLOQ 0.5-5.0 ng/mL.",
      "The PopPK model was then applied by MAP Bayesian estimation (MAXEVAL = 0)",
      "to Phase 2a Cohort 2, whose data did not enter the estimation."
    )
  )

  ini({
    lcl <- log(2.17)
    label("Apparent clearance CL/F at BMI 24.7 kg/m^2 (L/h)") # Supplement Table S2: CL/F = 2.17 L/h (RSE 2.79%)
    lvc <- log(25.5)
    label("Apparent central volume Vc/F (L)") # Supplement Table S2: Vc/F = 25.5 L (RSE 5.69%)
    lq <- log(8.54)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Supplement Table S2: Q/F = 8.54 L/h (RSE 3.77%)
    lvp <- log(447)
    label("Apparent peripheral volume Vp/F (L)") # Supplement Table S2: Vp/F = 447 L (RSE 2.37%)

    lfdepot <- fixed(log(1))
    label("Relative bioavailability of the fed tablet, reference (unitless)") # Supplement Table S2 / main text Methods: F of capsule fasted is estimated relative to tablet fed, the reference (F = 1)
    lfdepot_capsule <- log(0.683)
    label("Relative bioavailability of the fasted capsule vs fed tablet (unitless)") # Supplement Table S2: F fasted,capsule = 0.683 (RSE 6.63%)
    limax_fdepot <- fixed(log(1))
    label("Maximal dose-dependent reduction of F for fasted solution / tablet, Imax (unitless)") # Supplement Table S2: F fasted,solution,tablet Imax = Fixed 1
    led50 <- fixed(log(400))
    label("Dose giving half-maximal reduction of F for fasted solution / tablet, ID50 (mg)") # Supplement Table S2: F fasted,solution,tablet ID50 = Fixed 400 mg

    lka <- log(0.0473)
    label("Absorption rate constant Ka, fasted tablet or capsule (1/h)") # Supplement Table S2: Ka,fasted,tablet,capsule = 0.0473 1/h (RSE 4.06%)
    lka_solution <- log(0.0417)
    label("Absorption rate constant Ka, fasted oral solution at 400 mg (1/h)") # Supplement Table S2: Ka,fasted,solution,400mg = 0.0417 1/h (RSE 4.63%)
    lka_fed <- log(0.0536)
    label("Absorption rate constant Ka, fed tablet (1/h)") # Supplement Table S2: Ka,fed = 0.0536 1/h (RSE 5.19%)
    e_dose_ka <- -0.152
    label("Power exponent of dose/400 mg on the fasted-solution Ka (unitless)") # Supplement Table S2: Dose on Ka,fasted,solution gamma_D,Ka = -0.152 (RSE 29.1%)

    lktr <- log(6.26)
    label("Transit rate constant KTR, fasted (1/h)") # Supplement Table S2: KTR,fasted = 6.26 1/h (RSE 3.71%)
    lktr_fed <- log(4.19)
    label("Transit rate constant KTR, fed tablet (1/h)") # Supplement Table S2: KTR,fed = 4.19 1/h (RSE 5.70%)

    e_bmi_cl <- -1.69
    label("Power exponent of BMI/24.7 on CL/F (unitless)") # Supplement Table S2: BMI on CL/F gamma_BMI,CL/F = -1.69 (RSE 8.40%)
    e_bmi_fdepot <- -1.58
    label("Power exponent of BMI/24.7 on F (unitless)") # Supplement Table S2: BMI on F gamma_BMI,F = -1.58 (RSE 10.1%)
    e_bmi_ka <- 1.39
    label("Power exponent of BMI/24.7 on Ka (unitless)") # Supplement Table S2: BMI on Ka gamma_BMI,Ka = 1.39 (RSE 14.5%)

    # IIV / IOV variances = log(CV^2 + 1); Supplement Table S2 footnote:
    # CV% = [exp(variance) - 1]^(1/2) * 100.
    etalcl ~ 0.083445 # Supplement Table S2: omega CL/F IIV 29.5 %CV
    etalvc ~ 0.23119 # Supplement Table S2: omega Vc/F IIV 51.0 %CV
    etalfdepot ~ 0.10516 # Supplement Table S2: omega F IIV 33.3 %CV
    etalka ~ 0.12701 # Supplement Table S2: omega Ka IIV 36.8 %CV
    etalktr ~ 0.19499 # Supplement Table S2: omega KTR IIV 46.4 %CV

    # Inter-occasion variability, one indicator-multiplexed eta per occasion
    # sharing one variance (NONMEM $OMEGA BLOCK SAME).
    etaiov_fdepot_1 ~ 0.050245 # Supplement Table S2: omega F IOV 22.7 %CV
    etaiov_fdepot_2 ~ fixed(0.050245) # Supplement Table S2: omega F IOV 22.7 %CV, shared across occasions
    etaiov_ktr_1 ~ 0.21283 # Supplement Table S2: omega KTR IOV 48.7 %CV
    etaiov_ktr_2 ~ fixed(0.21283) # Supplement Table S2: omega KTR IOV 48.7 %CV, shared across occasions

    expSd <- 0.15212
    label("Residual error SD on the log scale (proportional, 15.3 %CV)") # Supplement Table S2: proportional residual error 15.3 %CV; sqrt(log(1 + 0.153^2)) = 0.15212
  })

  model({
    # Covariate terms (Supplement Table S2 footnote a: (BMI/24.7)^gamma)
    bmi_ratio <- BMI / 24.7
    tablet_fasted <- (1 - FORM_CAPSULE) * (1 - FORM_SOLUTION) * (1 - FED)

    # Occasion indicators for IOV; OCC = 0 switches IOV off.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    iov_fdepot <- occ1 * etaiov_fdepot_1 + occ2 * etaiov_fdepot_2
    iov_ktr <- occ1 * etaiov_ktr_1 + occ2 * etaiov_ktr_2

    # Relative bioavailability (Supplement Table S2 note):
    # fasted solution / fasted tablet F = 1 - Imax*DOSE/(ID50 + DOSE);
    # fasted capsule F = 0.683; fed tablet F = 1 (reference).
    imax_fdepot <- exp(limax_fdepot)
    ed50 <- exp(led50)
    f_dose <- 1 -
      imax_fdepot * DOSE_MOSNODENVIR_MG / (ed50 + DOSE_MOSNODENVIR_MG)
    f_form <- FORM_CAPSULE *
      exp(lfdepot_capsule) +
      FORM_SOLUTION * f_dose +
      tablet_fasted * f_dose +
      (1 - FORM_CAPSULE) * (1 - FORM_SOLUTION) * FED * exp(lfdepot)
    fdepot <- f_form * bmi_ratio^e_bmi_fdepot * exp(etalfdepot + iov_fdepot)

    # Absorption rate constant (Supplement Table S2 note):
    # fasted solution Ka = Ka,400mg * (DOSE/400)^gamma_D,Ka.
    ka_form <- FORM_SOLUTION *
      exp(lka_solution) *
      (DOSE_MOSNODENVIR_MG / 400)^e_dose_ka +
      FORM_CAPSULE * exp(lka) +
      tablet_fasted * exp(lka) +
      (1 - FORM_CAPSULE) * (1 - FORM_SOLUTION) * FED * exp(lka_fed)
    ka <- ka_form * bmi_ratio^e_bmi_ka * exp(etalka)

    # Transit rate constant: food effect only (Supplement, PopPK model).
    ktr <- exp((1 - FED) * lktr + FED * lktr_fed + etalktr + iov_ktr)

    cl <- exp(lcl + etalcl) * bmi_ratio^e_bmi_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Supplement Figure S1: Depot -> Transit 1..6 (KTR) -> Central (Ka).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # mg / L -> ng/mL
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
