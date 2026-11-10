Lee_2024_colistin_plasma <- function() {
  description <- "One-compartment population PK model with Weibull (time-dependent first-order) absorption for colistin in plasma after a 30-min jet nebulization of low-dose colistimethate sodium (CMS) in ventilated adults with ventilator-associated pneumonia caused by carbapenem-resistant Acinetobacter baumannii (Lee 2024). Fit to plasma concentrations only; the companion epithelial-lining-fluid model is Lee_2024_colistin_elf. Clearance and volume are apparent (CL/fm, V/fm) and the nebulized dose carries a fixed 9% bioavailability. No covariates."
  reference <- "Lee DH, Kim SY, Kim YK, Jung SY, Jang JH, Jang HJ, Lee JH, Jang SM, Lipman J. Intrapulmonary and Systemic Pharmacokinetics of Colistin Following Nebulization of Low-Dose Colistimethate Sodium in Patients with Ventilator-Associated Pneumonia Caused by Carbapenem-Resistant Acinetobacter baumannii. Antibiotics (Basel). 2024;13(3):258. doi:10.3390/antibiotics13030258"
  vignette <- "Lee_2024_colistin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The dose is colistin base activity (CBA): Lee 2024 Figure 1 caption gives
  # one nebulization as 0.9375 MIU CMS = 75 mg CMS = 28.125 mg CBA. Figure 4
  # prints fAUC24 = 2,182 ug/L*h for '75 mg qid', which equals
  # 0.09 * 0.34 * 4 * 28.125 mg / CL.
  compartmentData <- list(
    depot = list(analyte = "colistin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "colistin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    n_observations = 100,
    age_range = "mean 68.9 years (SD 12.6)",
    weight_range = "mean 64.9 kg (SD 15.6)",
    sex_female_pct = 35,
    race_ethnicity = "Korean (single centre, Busan, Republic of Korea)",
    disease_state = "Mechanically ventilated ICU adults with ventilator-associated pneumonia caused by carbapenem-resistant Acinetobacter baumannii; 4/20 on ECMO, 2/20 on CRRT, 3/20 in shock.",
    dose_range = "Nebulized CMS 0.9375 MIU (75 mg CMS, 28.125 mg CBA) over 30 min by jet nebulizer every 6 h (3.75 MIU/day)",
    regions = "Republic of Korea",
    notes = "Baseline demographics from Lee 2024 Table 1 (13 male / 7 female). 100 plasma concentrations at pre-dose, 1, 2, 3 and 5 h after the start of a nebulization. The pre-dose concentrations were zero (Figure 1), i.e. the sampled dose was the first. NONMEM 7.5.0, FOCEI. Covariates screened (demographics, comorbidities, ECMO, CRRT, laboratory values, renal function) but none retained."
  )

  ini({
    # Lee 2024 Table 2, 'Plasma' columns. CL and V are apparent values divided
    # by fm, the unidentifiable fraction of CMS converted to colistin. CL is in
    # mL/min in the paper and is converted to L/h in model().
    lcl <- log(26.3); label("Apparent colistin clearance from plasma, CL/fm (mL/min)") # Table 2 plasma: CL/fm = 26.3 mL/min (RSE 27.3%)
    lvc <- log(67.4); label("Apparent colistin plasma volume of distribution, V/fm (L)") # Table 2 plasma: V/fm = 67.4 L (RSE 22.2%)

    # Weibull absorption, WB = KAMAX * (1 - exp(-(RA * tad)^GAM)); see model().
    lkamax <- log(1.20); label("Weibull-absorption maximum rate constant, KAMAX (1/h)") # Table 2 plasma: KAMAX = 1.20 (RSE 17.2%)
    lra <- fixed(log(2.70)); label("Weibull-absorption rate-scaling parameter, RA (1/h)") # Table 2 plasma: RA = 2.70, footnote '*, fixed value'
    lgam1 <- fixed(log(1.79)); label("Weibull-absorption shape parameter, GAM (unitless)") # Table 2 plasma: GAM = 1.79, footnote '*, fixed value'

    # Nebulized-dose bioavailability. Section 4.5: 'The bioavailability of
    # colistin was assumed to be 9% [24]' (Boisson 2014). It is applied inside
    # the fitted model, not only in the fAUC post-processing: the Figure 4
    # typical-value simulations reproduce the observed Figure 1 medians only
    # with F = 0.09.
    lfdepot <- fixed(log(0.09)); label("Bioavailability of the nebulized CMS dose (fraction)") # Section 4.5 and Discussion: bioavailability assumed 9%

    # IIV. The Table 2 percentages are sqrt(omega^2) * 100 (the PsN
    # convention), so omega^2 = (pct/100)^2; the companion ELF model file
    # records the Figure S3 VPC check that settles this. For plasma the same
    # reading reproduces the simulated median of the Figure S4 VPC at 3 h
    # (25.4 ug/L, digitised 25.4) where log(1 + CV^2) gives 29.8.
    # GAM IIV is reported as '0 *' (held at zero) and so carries no eta.
    etalcl ~ 3.8416 # Table 2 plasma: IIV CL/fm = 196%; omega^2 = 1.96^2
    etalvc ~ 1.2321 # Table 2 plasma: IIV V/fm = 111%; omega^2 = 1.11^2
    etalkamax ~ 0.413449 # Table 2 plasma: IIV KAMAX = 64.3%; omega^2 = 0.643^2
    etalra ~ fixed(0.6889) # Table 2 plasma: IIV RA = 83.0%, held constant per the table footnote; omega^2 = 0.83^2

    propSd <- 0.185; label("Proportional residual error (fraction)") # Table 2 plasma: proportional error = 18.5% (RSE 11.3%)
  })

  model({
    # Paper units are CL in mL/min; convert to L/h so that a dose in mg gives
    # concentrations in mg/L on an hourly clock.
    cl <- exp(lcl + etalcl) * 60 / 1000
    vc <- exp(lvc + etalvc)
    kamax <- exp(lkamax + etalkamax)
    ra <- exp(lra + etalra)
    gam1 <- exp(lgam1)

    kel <- cl / vc

    # Lee 2024 Eq. (2), WB = KAMAX * (1 - exp(-(RA * time)^GAM)), with time
    # taken as time after the start of the current nebulization. The printed
    # MathML groups the exponent as (exp(-RA * time))^GAM, which collapses to
    # a single rate RA * GAM; the standard Weibull form of the cited
    # isavuconazole model (ref. 15, Desai 2016) is used, as in the ELF model.
    # tad(depot) is NA before the first dose, which the guard maps to zero.
    wb <- 0.0
    if (tad(depot) >= 0.0) wb <- kamax * (1 - exp(-(ra * tad(depot))^gam1))

    # Lee 2024 Eqs. (1) and (3).
    d/dt(depot) <- -wb * depot
    d/dt(central) <- wb * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Total plasma colistin concentration (mg/L). The paper's unbound plasma
    # concentration assumes fu = 0.34 (Mohamed 2012) and is not part of the model.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
