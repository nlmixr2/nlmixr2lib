Kaikousidis_2023_donepezil_combined_mouse <- function() {
  description <- "Preclinical (mouse). Joint blood + brain PK model for donepezil fitted simultaneously to oral-solution and nasal-film raw data: the oral arm is the three-compartment first-order-absorption model and the intranasal arm the three-way split (lagged zero-order blood and direct-brain inputs plus a Savic transit nose-to-brain chain), sharing the blood elimination, blood-to-brain transfer and brain elimination rate constants, with route-specific apparent volumes (Kaikousidis 2023 Table 4, 'Total Fit')"
  reference <- "Kaikousidis C, Papakyriakopoulou P, Dokoumetzidis A, Valsami G. Donepezil Brain and Blood Pharmacokinetic Modeling after Nasal Film and Oral Solution Administration in Mice. Pharmaceutics. 2023;15(5):1409. doi:10.3390/pharmaceutics15051409"
  vignette <- "Kaikousidis_2023_donepezil_mouse"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "depot_nasal")
  # Route-specific residual magnitudes gated by ROUTE_ORAL in model().
  paper_specific_residual_sds <- c("addSd_oral", "addSd_intranasal", "addSd_Cbrain_oral", "propSd_Cbrain_intranasal")

  covariateData <- list(
    ROUTE_ORAL = list(
      description = "Administration-route indicator: 1 = oral solution by gastric gavage, 0 = intranasal nasal film.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intranasal nasal film)",
      notes = "Per-subject indicator carried on BOTH dose and observation records. Kaikousidis 2023 fits the two routes jointly (Section 3.3, Table 4) with route-specific apparent blood and brain volumes (Vblood_IN / Vbrain_IN vs Vblood / Vbrain_PO), the blood-peripheral exchange present only in the oral sub-model (Figures 2 and 4), and route-specific residual-error models (additive on both matrices for oral; additive on blood and proportional on brain for intranasal). The absorption branch is carried by the dose record's target compartment instead -- oral doses go to `depot`, nasal-film doses to `depot_nasal` -- so ROUTE_ORAL must agree with that choice.",
      source_name = "route of administration (Per Os / IN)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    depot_nasal = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    depot_brain = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "donepezil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6J)",
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = "8 weeks",
    weight_range = "about 25 g (10 mg/kg for a 0.25 mg oral dose; 4 mg/kg for a 0.1 mg nasal film)",
    sex_female_pct = NA_real_,
    disease_state = "healthy",
    dose_range = "single oral dose of 0.25 mg (10 mg/kg, 0.2 mL of a 1.25 mg/mL solution by gastric gavage) or a single intranasal nasal film of 0.1 mg (4 mg/kg)",
    notes = "Destructive (one sample per animal) design in two groups: oral (serum and perfused brain at 5-360 min, n >= 5 per time point) and intranasal (5-120 min, n = 3-5 per time point) (Section 2.1). Raw concentrations of every animal of both groups pooled in one fit (Section 3.3); the authors present it as a robustness check of the separate oral and intranasal fits and did not refine it further. Values below the LOQ (serum 0.043 ug/mL, brain 0.038 ug/mL) were set to LOQ/2. Sex and the exact number of animals are not reported. No between-animal variability is estimated; the model is a typical-value model with residual error only."
  )

  ini({
    # Kaikousidis 2023 Table 4, 'Total Fit - Per Os and IN Raw Data'.
    # Shared by both routes (Section 3.3).
    lkel <- log(0.039); label("Elimination rate constant from blood, paper k20 (1/min)") # Table 4 'k20 (1/min)' 0.039
    lk_central_brain <- log(0.011); label("Blood-to-brain transfer rate constant, paper k23 (1/min)") # Table 4 'k23 (1/min)' 0.011
    lkeff_brain <- log(0.276); label("Lumped brain elimination rate constant, paper k30 (1/min)") # Table 4 'k30 (1/min)' 0.276

    # Oral sub-model (Figure 2).
    lka <- log(0.065); label("First-order oral absorption rate constant, paper k12 (1/min)") # Table 4 'k12 (1/min)' 0.065
    lk12 <- log(0.018); label("Blood-to-peripheral rate constant (oral sub-model), paper k24 (1/min)") # Table 4 'k24 (1/min)' 0.018
    lk21 <- log(0.031); label("Peripheral-to-blood rate constant (oral sub-model), paper k42 (1/min)") # Table 4 'k42 (1/min)' 0.031
    lvc_oral <- fixed(log(0.279)); label("Apparent blood volume of distribution after oral dosing, Vblood (L)") # not listed in Table 4; Table 1 'Vblood (L)' 0.279 as fixed in the oral sub-model (Section 3.1)
    lv_brain_oral <- log(0.015); label("Apparent brain volume of distribution after oral dosing, Vbrain_PO (L)") # Table 4 'Vbrain_PO (L)' 0.015

    # Intranasal sub-model (Figure 4).
    ltlag <- log(3.59); label("Absorption lag time of the direct blood and brain inputs, Tlag (min)") # Table 4 'Tlag (min)' 3.59 (SE 0.004)
    ld1 <- log(11.4); label("Duration of the zero-order direct blood and brain inputs, Tk0 (min)") # Table 4 'Tk0 (min)' 11.4 (SE 0.007)
    lmtt <- log(53.8); label("Mean transit time of the delayed nose-to-brain chain, Mtt (min)") # Table 4 'Mtt (min)' 53.8
    lktr <- log(0.205); label("Transit rate constant, ktr (1/min)") # Table 4 'ktr (1/min)' 0.205
    lvc_intranasal <- log(0.024); label("Apparent blood volume of distribution after intranasal dosing, Vblood_IN (L)") # Table 4 'Vblood_IN (L)' 0.024
    lv_brain_intranasal <- log(0.005); label("Apparent brain volume of distribution after intranasal dosing, Vbrain_IN (L)") # Table 4 'Vbrain_IN (L)' 0.005
    logitfrel <- log(0.557 / (1 - 0.557)); label("Logit of the dose fraction absorbed directly into blood, P1 (unitless)") # Table 4 'P1' 0.557
    logitfrel2 <- log(0.306 / (1 - 0.306)); label("Logit of the share of the non-blood dose absorbed directly into brain, P2 (unitless)") # Table 4 'P2' 0.306
    lka_brain <- log(0.127); label("First-order absorption rate constant from the delayed nose-to-brain compartment into brain, ka (1/min)") # Table 4 'ka (1/min)' 0.127

    # Residual error (Section 3.3): additive on both matrices for oral data;
    # additive on blood and proportional on brain for intranasal data. The
    # Table 4 rows are a1, a2, a3, b1; the footnote's b1/b2 labels do not
    # match those rows, so they are assigned by the text and by agreement
    # with the separate IN raw fit (Table 3: a1 0.082, b2 0.238). Read as
    # NONMEM variances; see the vignette.
    addSd_oral <- sqrt(0.017); label("Additive residual error on serum concentration, oral data (ug/mL)") # Table 4 'a1' 0.017 (variance)
    addSd_Cbrain_oral <- sqrt(0.004); label("Additive residual error on brain concentration, oral data (ug/mL)") # Table 4 'a2' 0.004 (variance)
    addSd_intranasal <- sqrt(0.084); label("Additive residual error on serum concentration, intranasal data (ug/mL)") # Table 4 'a3' 0.084 (variance)
    propSd_Cbrain_intranasal <- sqrt(0.236); label("Proportional residual error on brain concentration, intranasal data (fraction)") # Table 4 'b1' 0.236 (variance)
  })
  model({
    kel <- exp(lkel)
    k_central_brain <- exp(lk_central_brain)
    keff_brain <- exp(lkeff_brain)
    ka <- exp(lka)
    # Blood-peripheral exchange belongs to the oral sub-model only
    # (Figures 2 and 4; Section 3.3 joins 'the two previous individual
    # models'). After an intranasal dose peripheral1 stays empty.
    k12 <- exp(lk12) * ROUTE_ORAL
    k21 <- exp(lk21)
    vc <- exp(lvc_oral * ROUTE_ORAL + lvc_intranasal * (1 - ROUTE_ORAL))
    v_brain <- exp(lv_brain_oral * ROUTE_ORAL + lv_brain_intranasal * (1 - ROUTE_ORAL))
    tlag <- exp(ltlag)
    d1 <- exp(ld1)
    mtt <- exp(lmtt)
    ktr <- exp(lktr)
    frel <- expit(logitfrel)
    frel2 <- expit(logitfrel2)
    ka_brain <- exp(lka_brain)

    # Intranasal dose split (Section 3.2): blood P1, direct brain
    # (1 - P1) * P2, delayed brain (1 - P1) * (1 - P2).
    f_blood <- frel
    f_brain <- (1 - frel) * frel2
    f_transit <- (1 - frel) * (1 - frel2)
    # Savic transit with ktr = (N + 1) / Mtt (Section 3.2).
    ntr <- mtt * ktr - 1

    # Nasal-film doses are recorded on depot_nasal and released by an
    # explicit zero-order input over [tlag, tlag + d1] and an explicit Savic
    # gamma-density input; the equivalent f()/alag()/dur() form fails to
    # solve in rxode2 5.1.8 (see Kaikousidis_2023_donepezil_intranasal_meandata_mouse).
    # podo()/tad() refer to the most recent nasal-film dose.
    t_in <- tad(depot_nasal)
    dose_in <- podo(depot_nasal)
    zo_in <- 0
    if (t_in > tlag && t_in <= tlag + d1) zo_in <- dose_in / d1
    transit_in <- 0
    if (t_in > 0) transit_in <- f_transit * exp(log(dose_in) + log(ktr) + ntr * log(ktr * t_in) - ktr * t_in - lgamma(ntr + 1))

    d/dt(depot) <- -ka * depot
    d/dt(depot_nasal) <- -(f_blood + f_brain) * zo_in - transit_in
    d/dt(depot_brain) <- transit_in - ka_brain * depot_brain
    d/dt(central) <- ka * depot + f_blood * zo_in - (kel + k12 + k_central_brain) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(brain) <- f_brain * zo_in + ka_brain * depot_brain + k_central_brain * central - keff_brain * brain

    Cc <- central / vc
    Cbrain <- brain / v_brain
    ruvAdd <- addSd_oral * ROUTE_ORAL + addSd_intranasal * (1 - ROUTE_ORAL)
    ruvAdd_Cbrain <- addSd_Cbrain_oral * ROUTE_ORAL
    ruvProp_Cbrain <- propSd_Cbrain_intranasal * (1 - ROUTE_ORAL)
    Cc ~ add(ruvAdd)
    Cbrain ~ add(ruvAdd_Cbrain) + prop(ruvProp_Cbrain)
  })
}
