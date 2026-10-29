Lee_2024_colistin_elf <- function() {
  description <- "One-compartment population PK model with Weibull (time-dependent first-order) absorption for colistin in epithelial lining fluid (ELF) after a 30-min jet nebulization of low-dose colistimethate sodium (CMS) in ventilated adults with ventilator-associated pneumonia caused by carbapenem-resistant Acinetobacter baumannii (Lee 2024). Fit to urea-corrected bronchoalveolar-lavage ELF concentrations only; the companion plasma model is Lee_2024_colistin_plasma. Clearance and volume are apparent (CL/fm, V/fm) and the nebulized dose carries a fixed 9% bioavailability. No covariates."
  reference <- "Lee DH, Kim SY, Kim YK, Jung SY, Jang JH, Jang HJ, Lee JH, Jang SM, Lipman J. Intrapulmonary and Systemic Pharmacokinetics of Colistin Following Nebulization of Low-Dose Colistimethate Sodium in Patients with Ventilator-Associated Pneumonia Caused by Carbapenem-Resistant Acinetobacter baumannii. Antibiotics (Basel). 2024;13(3):258. doi:10.3390/antibiotics13030258"
  vignette <- "Lee_2024_colistin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The dose is colistin base activity (CBA): Lee 2024 Figure 1 caption gives
  # one nebulization as 0.9375 MIU CMS = 75 mg CMS = 28.125 mg CBA. The CBA
  # convention is pinned by Figure 4, whose printed fAUC24 for '75 mg qid'
  # (31,192 ug/L*h) equals 0.09 * 0.10 * 4 * 28.125 mg / CL exactly.
  compartmentData <- list(
    depot = list(analyte = "colistin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "colistin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    n_observations = 80,
    age_range = "mean 68.9 years (SD 12.6)",
    weight_range = "mean 64.9 kg (SD 15.6)",
    sex_female_pct = 35,
    race_ethnicity = "Korean (single centre, Busan, Republic of Korea)",
    disease_state = "Mechanically ventilated ICU adults with ventilator-associated pneumonia caused by carbapenem-resistant Acinetobacter baumannii; 4/20 on ECMO, 2/20 on CRRT, 3/20 in shock.",
    dose_range = "Nebulized CMS 0.9375 MIU (75 mg CMS, 28.125 mg CBA) over 30 min by jet nebulizer every 6 h (3.75 MIU/day)",
    regions = "Republic of Korea",
    notes = "Baseline demographics from Lee 2024 Table 1 (13 male / 7 female). 80 ELF concentrations from bronchoalveolar lavage at pre-dose, 1, 3 and 5 h after the start of a nebulization; ELF concentration = BAL concentration x (plasma urea / BAL urea). The pre-dose concentrations were zero (Figure 1), i.e. the sampled dose was the first. NONMEM 7.5.0, FOCEI. Covariates screened (demographics, comorbidities, ECMO, CRRT, laboratory values, renal function) but none retained."
  )

  ini({
    # Lee 2024 Table 2, 'Epithelial Lining Fluid' columns. CL and V are
    # apparent values divided by fm, the unidentifiable fraction of CMS
    # converted to colistin. CL is in mL/min and V in mL in the paper; both
    # are converted to L/h and L in model().
    lcl <- log(0.541); label("Apparent colistin clearance from ELF, CL/fm (mL/min)") # Table 2 ELF: CL/fm = 0.541 mL/min (RSE 23.2%)
    lvc <- log(23.0); label("Apparent colistin ELF volume of distribution, V/fm (mL)") # Table 2 ELF: V/fm = 23.0 mL (RSE 28.9%)

    # Weibull absorption, WB = KAMAX * (1 - exp(-(RA * tad)^GAM)); see model().
    lkamax <- log(0.531); label("Weibull-absorption maximum rate constant, KAMAX (1/h)") # Table 2 ELF: KAMAX = 0.531 (RSE 16.1%)
    lra <- fixed(log(2.16)); label("Weibull-absorption rate-scaling parameter, RA (1/h)") # Table 2 ELF: RA = 2.16, footnote '*, fixed value'
    lgam1 <- fixed(log(2.92)); label("Weibull-absorption shape parameter, GAM (unitless)") # Table 2 ELF: GAM = 2.92, footnote '*, fixed value'

    # Nebulized-dose bioavailability. Section 4.5: 'The bioavailability of
    # colistin was assumed to be 9% [24]' (Boisson 2014). It is applied inside
    # the fitted model, not only in the fAUC post-processing: the Figure 4
    # typical-value simulations reproduce the observed Figure 1 medians only
    # with F = 0.09 (with F = 1 every prediction is 11-fold too high).
    lfdepot <- fixed(log(0.09)); label("Bioavailability of the nebulized CMS dose (fraction)") # Section 4.5 and Discussion: bioavailability assumed 9%

    # IIV. Lee 2024 reports the variability as a percentage under
    # theta_i = theta * exp(eta_i). The percentage is sqrt(omega^2) * 100 (the
    # PsN convention), so omega^2 = (pct/100)^2: under that reading the model
    # reproduces the simulated 10th/50th/90th percentiles of the Figure S3 VPC
    # (1 h: 1,600 / 10,700 / 41,700 ug/L against 1,580 / ~10,500 / 41,100
    # digitised), whereas log(1 + CV^2) puts the 10th percentile 60% too high.
    # GAM IIV is reported as '0 *' (held at zero) and so carries no eta.
    etalcl ~ 1.5376 # Table 2 ELF: IIV CL/fm = 124%; omega^2 = 1.24^2
    etalvc ~ 0.356409 # Table 2 ELF: IIV V/fm = 59.7%; omega^2 = 0.597^2
    etalkamax ~ 0.710649 # Table 2 ELF: IIV KAMAX = 84.3%; omega^2 = 0.843^2
    etalra ~ fixed(0.6889) # Table 2 ELF: IIV RA = 83.0%, held constant per the table footnote; omega^2 = 0.83^2

    propSd <- 0.271; label("Proportional residual error (fraction)") # Table 2 ELF: proportional error = 27.1% (RSE 18.6%)
  })

  model({
    # Paper units are CL in mL/min and V in mL; convert to L/h and L so that a
    # dose in mg gives concentrations in mg/L on an hourly clock.
    cl <- exp(lcl + etalcl) * 60 / 1000
    vc <- exp(lvc + etalvc) / 1000
    kamax <- exp(lkamax + etalkamax)
    ra <- exp(lra + etalra)
    gam1 <- exp(lgam1)

    kel <- cl / vc

    # Lee 2024 Eq. (2), WB = KAMAX * (1 - exp(-(RA * time)^GAM)), with time
    # taken as time after the start of the current nebulization. The printed
    # MathML groups the exponent as (exp(-RA * time))^GAM, which collapses to
    # a single rate RA * GAM and leaves the two fixed values meaningless; the
    # standard Weibull form of the cited isavuconazole model (ref. 15, Desai
    # 2016) reproduces the Figure 4 peak height and timing (see vignette).
    # tad(depot) is NA before the first dose, which the guard maps to zero.
    wb <- 0.0
    if (tad(depot) >= 0.0) wb <- kamax * (1 - exp(-(ra * tad(depot))^gam1))

    # Lee 2024 Eqs. (1) and (3).
    d/dt(depot) <- -wb * depot
    d/dt(central) <- wb * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Total colistin concentration in ELF (mg/L). The paper's unbound ELF
    # concentration assumes fu = 0.10 (Gkoufa 2022) and is not part of the model.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
