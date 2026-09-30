Krzyzanski_2021_betamethasone <- function() {
  description <- "Two-compartment population PK model for betamethasone after single 6 mg doses of betamethasone phosphate given intramuscularly or orally, or of a 1:1 betamethasone phosphate/acetate IM suspension (Celestone), to healthy nonpregnant Indian women. Three parallel first-order absorption depots: IM phosphate, oral phosphate, and a slow IM acetate depot (flip-flop terminal phase). Clearances and volumes are apparent (divided by the IM bioavailability FIM); oral and acetate bioavailabilities are relative to IM phosphate. No covariates."
  reference <- "Krzyzanski W, Milad MA, Jobe AH, Peppard T, Bies RR, Jusko WJ. Population pharmacokinetic modeling of intramuscular and oral dexamethasone and betamethasone in Indian women. J Pharmacokinet Pharmacodyn. 2021;48(2):261-272. doi:10.1007/s10928-020-09730-z"
  vignette <- "Krzyzanski_2021_corticosteroids"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  dosing <- c("depot1", "depot2", "depot3")

  covariateData <- list()

  compartmentData <- list(
    depot1 = list(analyte = "betamethasone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "betamethasone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "betamethasone", units = "mg", specimen = "plasma", verified = TRUE),
    depot2 = list(analyte = "betamethasone", units = "mg", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "betamethasone", units = "mg", specimen = "administration site", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    age_range = "22-39 years",
    weight_range = "47.0-68.7 kg",
    weight_mean = "56.8 kg",
    bmi_range = "20.6-25.0 kg/m^2",
    height_range = "144-167 cm",
    sex_female_pct = 100,
    race_ethnicity = c(Indian = 100),
    disease_state = "healthy nonpregnant women of reproductive age",
    dose_range = "6 mg betamethasone (free-alcohol equivalent) single dose: betamethasone phosphate IM (solution) or PO (0.5 mg tablets), or IM Celestone (3 mg phosphate + 3 mg acetate suspension)",
    regions = "India (Bangalore)",
    notes = "Open-label, randomized, two-period partial crossover in 48 healthy women (NCT03668860; Krzyzanski 2021 Tables 1-2). Every woman received betamethasone at least once: BET-P IM (treatment B, n = 12), BET-PA IM (treatment C, n = 24) and BET-P PO (treatment E, n = 24); sequences CE and EC received two betamethasone doses 10 days apart. 949 betamethasone concentrations, 19 below the 0.1 ng/mL LLOQ (handled by Beal M3). The homogeneous ages, body weights and ethnicity precluded covariate analysis."
  )

  ini({
    lcl <- log(5.95)
    label("Apparent clearance CL/FIM (L/h)") # Table 4, CL/FIM = 5.95 L/h (RSE 8.6%); supplement $THETA 1
    lvc <- log(67.5)
    label("Apparent central volume Vp/FIM (L)") # Table 4, Vp/FIM = 67.5 L (RSE 3.2%); supplement $THETA 2
    lq <- log(0.173)
    label("Apparent distributional clearance CLD/FIM (L/h)") # Table 4, CLD/FIM = 0.173 L/h (no RSE printed); supplement $THETA 8
    lvp <- log(4.94)
    label("Apparent peripheral volume VT/FIM (L)") # Table 4, VT/FIM = 4.94 L (no RSE printed); supplement $THETA 9
    lka_im <- log(0.971)
    label("First-order IM absorption rate constant, phosphate, kaIM (1/h)") # Table 4, kaIM = 0.971 1/h (RSE 3.6%); supplement $THETA 3
    lka_im_acetate <- log(0.00638)
    label("First-order IM absorption/hydrolysis rate constant, acetate, kaIMa (1/h)") # Table 4, kaIMa = 0.00638 1/h (RSE 14.7%); supplement $THETA 4
    lka_oral <- log(1.21)
    label("First-order oral absorption rate constant kaPO (1/h)") # Table 4, kaPO = 1.21 1/h (RSE 1.5%); supplement $THETA 5
    lfdepot_oral <- log(0.935)
    label("Oral bioavailability relative to IM phosphate, Fr = FPO/FIM (fraction)") # Table 4, Fr = 0.935 (RSE 9.7%); supplement $THETA 6
    lfdepot_im_acetate <- log(0.819)
    label("IM acetate bioavailability relative to IM phosphate, Fra = FIMa/FIM (fraction)") # Table 4, Fra = 0.819 (RSE 9.1%); supplement $THETA 7

    # IIV (log-normal, Eq 10; supplement $PK THETA*EXP(ETA)). Table 4 variance
    # column with %CV in parentheses. No IIV on CLD/FIM or VT/FIM.
    etalcl + etalvc ~ c(0.0210, 0.0155, 0.0188) # Table 4, omega^2 CL/FIM = 0.0210 (14.6% CV), Cov = 0.0155 (r = 0.78), omega^2 Vp/FIM = 0.0188 (13.8% CV); supplement $OMEGA BLOCK(2)
    etalka_im ~ 0.0441 # Table 4, omega^2 kaIM = 0.0441 (21.2% CV)
    etalka_im_acetate ~ 0.147 # Table 4, omega^2 kaIMa = 0.147 (39.8% CV)
    etalka_oral ~ 0.241 # Table 4, omega^2 kaPO = 0.241 (52.2% CV)
    etalfdepot_oral ~ 0.0182 # Table 4, omega^2 Fr = 0.0182 (13.6% CV)
    etalfdepot_im_acetate ~ 0.00773 # Table 4, omega^2 Fra = 0.00773 (8.8% CV)

    # Residual error: log C = log Cp + e, e ~ N(0, sigma^2) (Eq 11; supplement
    # $ERROR Y = LOG(A(2)/V) + EPS(1)), i.e. lnorm; SD = sqrt(0.0211) = 0.1453.
    expSd <- sqrt(0.0211)
    label("Residual error SD on the log scale (unitless)") # Table 4, sigma^2 = 0.0211; supplement $SIGMA
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka_im <- exp(lka_im + etalka_im)
    ka_im_acetate <- exp(lka_im_acetate + etalka_im_acetate)
    ka_oral <- exp(lka_oral + etalka_oral)
    fr <- exp(lfdepot_oral + etalfdepot_oral)
    fra <- exp(lfdepot_im_acetate + etalfdepot_im_acetate)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # depot1 = IM betamethasone phosphate (Eq 3; supplement COMP IM, reference
    # route FIM); depot2 = oral betamethasone phosphate (Eq 4; COMP PO);
    # depot3 = IM betamethasone acetate (Eq 5; COMP IMA). A Celestone (BET-PA)
    # dose is split: half into depot1 and half into depot3.
    d/dt(depot1) <- -ka_im * depot1
    d/dt(central) <- ka_im * depot1 + ka_oral * depot2 + ka_im_acetate * depot3 - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(depot2) <- -ka_oral * depot2
    d/dt(depot3) <- -ka_im_acetate * depot3

    # Bioavailabilities relative to IM phosphate (Eq 6-7; supplement F4 = FR,
    # F5 = FRA); FIM is the reference (1).
    f(depot2) <- fr
    f(depot3) <- fra

    # Dose in mg, volumes in L: mg/L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
