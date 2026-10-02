Krzyzanski_2021_dexamethasone <- function() {
  description <- "Two-compartment population PK model for dexamethasone after single 6 mg doses of dexamethasone phosphate given intramuscularly or orally to healthy nonpregnant Indian women, with separate first-order absorption from an IM and an oral depot. Clearances and volumes are apparent (divided by the IM bioavailability FIM); oral bioavailability is relative to IM (Fr = FPO / FIM). No covariates."
  reference <- "Krzyzanski W, Milad MA, Jobe AH, Peppard T, Bies RR, Jusko WJ. Population pharmacokinetic modeling of intramuscular and oral dexamethasone and betamethasone in Indian women. J Pharmacokinet Pharmacodyn. 2021;48(2):261-272. doi:10.1007/s10928-020-09730-z"
  vignette <- "Krzyzanski_2021_corticosteroids"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  dosing <- c("depot1", "depot2")

  covariateData <- list()

  compartmentData <- list(
    depot1 = list(analyte = "dexamethasone", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "dexamethasone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dexamethasone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dexamethasone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 36,
    n_studies = 1,
    age_range = "22-39 years",
    weight_range = "47.0-68.7 kg",
    weight_mean = "56.8 kg",
    bmi_range = "20.6-25.0 kg/m^2",
    height_range = "144-167 cm",
    sex_female_pct = 100,
    race_ethnicity = c(Indian = 100),
    disease_state = "healthy nonpregnant women of reproductive age",
    dose_range = "6 mg dexamethasone (free-alcohol equivalent) as dexamethasone phosphate, single dose IM (solution) or PO (0.5 mg tablets)",
    regions = "India (Bangalore)",
    notes = "Open-label, randomized, two-period partial crossover in 48 healthy women (NCT03668860; Krzyzanski 2021 Tables 1-2). Dexamethasone was given to the 36 women of sequences AB, BA, CD, DC, ED and DE: 12 received DEX-P IM (treatment A) and 24 received DEX-P PO (treatment D). 578 dexamethasone concentrations, 103 below the 0.1 ng/mL LLOQ (handled by Beal M3). The homogeneous ages, body weights and ethnicity precluded covariate analysis."
  )

  ini({
    lcl <- log(9.29)
    label("Apparent clearance CL/FIM (L/h)") # Table 3, CL/FIM = 9.29 L/h (RSE 4.4%)
    lvc <- log(51.3)
    label("Apparent central volume Vp/FIM (L)") # Table 3, Vp/FIM = 51.3 L (RSE 4.4%)
    lq <- log(0.538)
    label("Apparent distributional clearance CLD/FIM (L/h)") # Table 3, CLD/FIM = 0.538 L/h (RSE 4.1%)
    lvp <- log(5.06)
    label("Apparent peripheral volume VT/FIM (L)") # Table 3, VT/FIM = 5.06 L (RSE 4.7%)
    lka_im <- log(0.460)
    label("First-order IM absorption rate constant kaIM (1/h)") # Table 3, kaIM = 0.460 1/h (RSE 8.8%)
    lka_oral <- log(0.936)
    label("First-order oral absorption rate constant kaPO (1/h)") # Table 3, kaPO = 0.936 1/h (RSE 15.2%)
    lfdepot_oral <- log(1.04)
    label("Oral bioavailability relative to IM, Fr = FPO/FIM (fraction)") # Table 3, Fr = 1.04 (RSE 5.3%); Eq 7

    # IIV (log-normal, Eq 10); Table 3 variance column with %CV in parentheses.
    # IIV on Vp/FIM was estimated close to 0 and fixed at 0 (Table 3 '0**';
    # Results). It is omitted rather than written as fixed(0) because a
    # zero-variance diagonal makes OMEGA singular for rxSolve's sampler.
    # No IIV was estimated on Fr, CLD/FIM or VT/FIM (Table 3 'NA').
    etalcl ~ 0.0265 # Table 3, omega^2 CL/FIM = 0.0265 (16.4% CV)
    etalka_im ~ 0.0633 # Table 3, omega^2 kaIM = 0.0633 (25.6% CV)
    etalka_oral ~ 0.395 # Table 3, omega^2 kaPO = 0.395 (69.6% CV)

    # Residual error: log C = log Cp + e, e ~ N(0, sigma^2) (Eq 11), i.e.
    # additive on the log scale = lnorm; SD = sqrt(0.0455) = 0.2133.
    expSd <- sqrt(0.0455)
    label("Residual error SD on the log scale (unitless)") # Table 3, sigma^2 = 0.0455 (RSE 18.9%)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka_im <- exp(lka_im + etalka_im)
    ka_oral <- exp(lka_oral + etalka_oral)
    fr <- exp(lfdepot_oral)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # depot1 = IM dexamethasone phosphate (Eq 3, reference route FIM);
    # depot2 = oral dexamethasone phosphate (Eq 4).
    d/dt(depot1) <- -ka_im * depot1
    d/dt(depot2) <- -ka_oral * depot2
    d/dt(central) <- ka_im * depot1 + ka_oral * depot2 - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Bioavailabilities relative to IM (Eq 6-7); FIM is the reference (1).
    f(depot2) <- fr

    # Dose in mg, volumes in L: mg/L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
