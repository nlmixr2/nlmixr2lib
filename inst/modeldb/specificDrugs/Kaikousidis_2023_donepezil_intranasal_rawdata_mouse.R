Kaikousidis_2023_donepezil_intranasal_rawdata_mouse <- function() {
  description <- "Preclinical (mouse). Blood + brain PK model for donepezil after a nasal film: the dose is split three ways into a lagged zero-order input to blood, a lagged zero-order direct nose-to-brain input, and a delayed nose-to-brain input through a Savic transit chain and a first-order brain absorption compartment; blood-to-brain transfer and brain elimination fixed from the oral raw-data fit. Fitted by naive pooling to the RAW per-animal serum and brain concentrations (Kaikousidis 2023 Table 3, 'IN - Raw Data')"
  reference <- "Kaikousidis C, Papakyriakopoulou P, Dokoumetzidis A, Valsami G. Donepezil Brain and Blood Pharmacokinetic Modeling after Nasal Film and Oral Solution Administration in Mice. Pharmaceutics. 2023;15(5):1409. doi:10.3390/pharmaceutics15051409"
  vignette <- "Kaikousidis_2023_donepezil_mouse"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")
  dosing <- "depot_nasal"

  covariateData <- list()

  compartmentData <- list(
    depot_nasal = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    depot_brain = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "serum", verified = TRUE),
    brain = list(analyte = "donepezil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6J)",
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = "8 weeks",
    weight_range = "about 25 g (4 mg/kg for a 0.1 mg nasal film; 10 mg/kg for a 0.25 mg oral dose)",
    sex_female_pct = NA_real_,
    disease_state = "healthy",
    dose_range = "single intranasal nasal film containing 0.1 mg donepezil (4 mg/kg)",
    notes = "Destructive (one sample per animal) design: serum and perfused brain sampled at 5, 10, 15, 30, 45, 60, 90 and 120 min after the nasal film, n = 3-5 animals per time point (Section 2.1). This file is the fit to the RAW concentrations of every animal pooled as one typical subject (Section 2.2); its mean-data sibling is Kaikousidis_2023_donepezil_intranasal_meandata_mouse. Values below the LOQ (serum 0.043 ug/mL, brain 0.038 ug/mL) were set to LOQ/2. Sex and the exact number of animals are not reported. No between-animal variability is estimated; the model is a typical-value model with residual error only."
  )

  ini({
    # Kaikousidis 2023 Table 3, 'IN - Raw Data' block.
    ltlag <- log(3.72); label("Absorption lag time of the direct blood and brain inputs, Tlag (min)") # Table 3 IN Raw 'Tlag (min)' 3.72
    ld1 <- log(11.6); label("Duration of the zero-order direct blood and brain inputs, Tk0 (min)") # Table 3 IN Raw 'Tk0 (min)' 11.6
    lmtt <- log(60.9); label("Mean transit time of the delayed nose-to-brain chain, Mtt (min)") # Table 3 IN Raw 'Mtt (min)' 60.9
    lktr <- log(0.183); label("Transit rate constant, ktr (1/min)") # Table 3 IN Raw 'ktr (1/min)' 0.183
    lvc <- log(0.027); label("Apparent blood volume of distribution after intranasal dosing, V1 (L)") # Table 3 IN Raw 'V1' 0.027 (footnote: Vblood; litres per Discussion)
    lv_brain <- log(0.003); label("Apparent brain volume of distribution after intranasal dosing, V2 (L)") # Table 3 IN Raw 'V2' 0.003 (footnote: Vbrain; litres per Discussion)
    logitfrel <- log(0.683 / (1 - 0.683)); label("Logit of the dose fraction absorbed directly into blood, P1 (unitless)") # Table 3 IN Raw 'P1' 0.683 = Section 3.2 'P2' in F2 = P2
    logitfrel2 <- log(0.276 / (1 - 0.276)); label("Logit of the share of the non-blood dose absorbed directly into brain, P2 (unitless)") # Table 3 IN Raw 'P2' 0.276 = Section 3.2 'P3' in F3 = (1 - P2) * P3
    lka_brain <- log(0.566); label("First-order absorption rate constant from the delayed nose-to-brain compartment into brain, ka (1/min)") # Table 3 IN Raw 'ka (1/min)' 0.566
    lkel <- log(0.05); label("Elimination rate constant from blood, paper k20 (1/min)") # Table 3 IN Raw 'k20 (1/min)' 0.05
    lk_central_brain <- fixed(log(0.008)); label("Blood-to-brain transfer rate constant, paper k23 (1/min)") # fixed to Table 2 Per Os Raw 'k23' 0.008 (Section 3.2)
    lkeff_brain <- fixed(log(0.369)); label("Lumped brain elimination rate constant, paper k30 (1/min)") # fixed to Table 2 Per Os Raw 'k30' 0.369 (Section 3.2)

    # Residual error: additive on blood, proportional on brain (Section 3.2).
    # Table 3 prints this block under a 'Per Os - Raw Data' heading, but it
    # sits in the IN table and its a1/b2 structure is the IN structure the
    # text describes. Read as NONMEM variances; see the vignette.
    addSd <- sqrt(0.082); label("Additive residual error on serum concentration (ug/mL)") # Table 3 (IN) Raw 'a1' 0.082 (variance)
    propSd_Cbrain <- sqrt(0.238); label("Proportional residual error on brain concentration (fraction)") # Table 3 (IN) Raw 'b2' 0.238 (variance)
  })
  model({
    tlag <- exp(ltlag)
    d1 <- exp(ld1)
    mtt <- exp(lmtt)
    ktr <- exp(lktr)
    vc <- exp(lvc)
    v_brain <- exp(lv_brain)
    frel <- expit(logitfrel)
    frel2 <- expit(logitfrel2)
    ka_brain <- exp(lka_brain)
    kel <- exp(lkel)
    k_central_brain <- exp(lk_central_brain)
    keff_brain <- exp(lkeff_brain)

    # Section 3.2: F2 = P2 (blood), F3 = (1 - P2) * P3 (direct brain),
    # F4 = (1 - P2) * (1 - P3) (delayed brain), with the text's P2 / P3
    # printed as P1 / P2 in Table 3.
    f_blood <- frel
    f_brain <- (1 - frel) * frel2
    f_transit <- (1 - frel) * (1 - frel2)

    # Savic transit: Section 3.2 reports N = 10.14 for this fit, which equals
    # ktr * Mtt - 1, i.e. ktr = (N + 1) / Mtt.
    ntr <- mtt * ktr - 1

    # The whole nasal-film dose is recorded on depot_nasal and released by a
    # zero-order input over [tlag, tlag + d1] shared by blood and direct
    # brain, plus the Savic gamma-density input into depot_brain. Written
    # out because the equivalent f()/alag()/dur() form fails to solve in
    # rxode2 5.1.8; see Kaikousidis_2023_donepezil_intranasal_meandata_mouse.
    # podo()/tad() refer to the most recent dose.
    t_in <- tad(depot_nasal)
    dose_in <- podo(depot_nasal)
    zo_in <- 0
    if (t_in > tlag && t_in <= tlag + d1) zo_in <- dose_in / d1
    transit_in <- 0
    if (t_in > 0) transit_in <- f_transit * exp(log(dose_in) + log(ktr) + ntr * log(ktr * t_in) - ktr * t_in - lgamma(ntr + 1))

    # Figure 4. No peripheral compartment after intranasal dosing.
    d/dt(depot_nasal) <- -(f_blood + f_brain) * zo_in - transit_in
    d/dt(depot_brain) <- transit_in - ka_brain * depot_brain
    d/dt(central) <- f_blood * zo_in - (kel + k_central_brain) * central
    d/dt(brain) <- f_brain * zo_in + ka_brain * depot_brain + k_central_brain * central - keff_brain * brain

    Cc <- central / vc
    Cbrain <- brain / v_brain
    Cc ~ add(addSd)
    Cbrain ~ prop(propSd_Cbrain)
  })
}
