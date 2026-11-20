Purohit_2026_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model for oral linezolid in plasma and",
    "lumbar cerebrospinal fluid (CSF) in Ugandan adults with HIV-associated",
    "tuberculous meningitis given standard- (10 mg/kg) or high-dose (35",
    "mg/kg) rifampicin (ALTER trial, Purohit 2026). First-order absorption",
    "(ka 0.40 1/h) and linear clearance (6.01 L/h; V 27.7 L). A daily",
    "rifampicin dose of 1650 or 2100 mg raises linezolid clearance by 34.2%",
    "relative to 300-1200 mg. CSF is a Sheiner-style effect compartment",
    "holding a concentration, equilibrating with plasma at kp 0.383 1/h",
    "toward a CSF-to-plasma partition coefficient of 0.535. Between-subject",
    "variability on clearance (17.7% CV) and on kp (92.5% CV). Residual",
    "error is combined proportional (33.3%) plus additive (0.26 mg/L) for",
    "plasma and proportional (31.6%) for CSF."
  )
  reference <- paste(
    "Purohit D, van Wijk R, Kafeero P, Kibengo F, Zimmerman M, D'artois V,",
    "Chow FC, Savic RM (2026). Impact of high-dose rifampicin on linezolid",
    "pharmacokinetics in tuberculous meningitis. Open Forum Infect Dis",
    "13(4):ofag154. doi:10.1093/ofid/ofag154",
    sep = " "
  )
  vignette <- "Purohit_2026_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state is an exception to the usual "states hold an amount"
  # rule. The Discussion describes CSF "as levels within a hypothetical
  # effect compartment" characterised by a rate (kp) and an extent (PC) of
  # entry (Figure 1), i.e. dC_CSF/dt = kp * (PC * C_plasma - C_CSF), so the
  # state carries mg/L and drug entering the CSF is not removed from plasma.
  compartmentData <- list(
    depot = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "linezolid",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    csf = list(
      analyte = "linezolid",
      units = "mg/L",
      specimen = "CSF",
      verified = TRUE
    )
  )

  covariateData <- list(
    DOSE_RIFAMPICIN_MG = list(
      description = paste(
        "Administered daily rifampicin dose (mg/day). Enters the model only",
        "through the indicator DOSE_RIFAMPICIN_MG >= 1650, which selects the",
        "34.2% higher linezolid clearance."
      ),
      units = "mg (per day)",
      type = "continuous",
      reference_category = paste(
        "Daily rifampicin doses of 300, 450, 600 or 1200 mg (the reference",
        "group; Results 'Participants receiving rifampicin at doses of 300,",
        "450, 600, or 1200 mg were considered the reference group')."
      ),
      notes = paste(
        "Weight-band dosing (Supplementary Table 1): standard-dose arms",
        "received 300 / 450 / 600 / 750 mg for 30-37 / 38-54 / 55-70 / >=71",
        "kg and high-dose arms 1200 / 1650 / 2100 / 2550 mg for the same",
        "bands, so the lightest high-dose band (1200 mg) falls in the",
        "reference group. The paper analysed the dose as a categorical",
        "variable (Discussion, limitations) and reports only the 1650-or-2100",
        "mg category effect; any threshold between 1200 and 1650 mg",
        "reproduces it for the studied doses. Doses of 750 and 2550 mg were",
        "not observed (all participants weighed 34-69 kg) and are assigned",
        "to the reference and high categories respectively by the >= 1650 mg",
        "threshold, which is an extrapolation."
      ),
      source_name = "absolute dose of rifampicin"
    )
  )

  population <- list(
    n_subjects = 18,
    n_studies = 1,
    species = "human",
    age_range = "median 38.5 years (range 26-72; IQR 33.3-42.8)",
    weight_range = "median 50 kg (range 34-69) on the first PK day",
    sex_female_pct = 50,
    race_ethnicity = "Ugandan; not further reported",
    disease_state = paste(
      "Definite (5), probable (12) or possible (1) tuberculous meningitis;",
      "all participants living with HIV (6 on antiretroviral therapy).",
      "MRC grade I 27.8%, II 50%, III 22.2%."
    ),
    dose_range = paste(
      "Linezolid 1200 mg once daily (n = 12) or 600 mg twice daily (n = 6)",
      "orally for the first 4 weeks, with rifampicin 10 mg/kg (n = 10) or",
      "35 mg/kg (n = 8) daily by weight band, plus isoniazid,",
      "pyrazinamide, ethambutol and dexamethasone."
    ),
    regions = "Uganda (Masaka Regional Referral Hospital)",
    renal_function = "Serum creatinine median 0.72 mg/dL (IQR 0.56-0.75)",
    notes = paste(
      "Linezolid substudy of the phase II randomised open-label ALTER trial",
      "(NCT04021121). 73 plasma and 30 CSF concentrations: intensive plasma",
      "sampling on day 1 (0, 2, 4, 8 h) and one sparse plasma sample on day",
      "14, with one lumbar CSF sample per visit randomised to an early (0-2",
      "h), medium (2-6 h) or late (6-8 h) window (5, 6 and 7 participants).",
      "Day-28 samples were excluded. Two of 20 linezolid-randomised",
      "participants with detectable pre-dose linezolid were excluded.",
      "Demographics: Table 1 and Results; individual concentrations:",
      "Supplementary Table 2."
    )
  )

  ini({
    lka <- log(0.40)
    label("Log of the first-order absorption rate constant (1/h)") # Table 2 ka 0.40 (.19-.61) [RSE 27.2%]
    lcl <- log(6.01)
    label("Log of apparent clearance with a 300-1200 mg/day rifampicin dose (L/h)") # Table 2 CL 6.01 (5.01-7.01) [RSE 8.5%]
    lvc <- log(27.70)
    label("Log of apparent volume of distribution (L)") # Table 2 V 27.70 (16.27-38.77) [RSE 20.4%]
    lfdepot <- fixed(log(1))
    label("Log of relative oral bioavailability (unitless)") # Table 2 F relative bioavailability 1 (fixed)

    lke0 <- log(0.383)
    label("Log of the plasma-to-CSF transfer (equilibration) rate constant kp (1/h)") # Table 2 kp 0.383 (.282-.484) [RSE 13.5%]
    lppc <- log(0.535)
    label("Log of the CSF-to-plasma partition coefficient PC (unitless)") # Table 2 PC 0.535 (.458-.612) [RSE 7.4%]

    e_dose_rifampicin_cl <- 0.342
    label("Fractional increase in CL with a rifampicin dose of 1650 or 2100 mg/day (unitless)") # Table 2 effect of 1650 or 2100 mg rifampicin on CL 34.2% (2.4-66.4) [RSE 47.4%]

    # Table 2 footnote: 'Interindividual variability is expressed as percent
    # coefficient of variation'; converted with omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.177^2) # Table 2 IIV on CL 17.7% CV [RSE 30.0%]
    etalke0 ~ log(1 + 0.925^2) # Table 2 IIV on kp 92.5% CV [RSE 24.7%]

    propSd <- 0.333
    label("Proportional residual error for plasma linezolid (fraction)") # Table 2 plasma proportional error 33.3% [RSE 14.0%]
    addSd <- 0.26
    label("Additive residual error for plasma linezolid (mg/L)") # Table 2 plasma additive error 0.26 ug/mL [RSE 57.6%]
    propSd_Ccsf <- 0.316
    label("Proportional residual error for CSF linezolid (fraction)") # Table 2 CSF proportional error 31.6% [RSE 14.5%]
  })

  model({
    # Rifampicin dose category: 1650 or 2100 mg/day vs 300-1200 mg/day
    # (Results, Pharmacokinetic Model; Table 2).
    rif_high <- (DOSE_RIFAMPICIN_MG >= 1650)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 + e_dose_rifampicin_cl * rif_high)
    vc <- exp(lvc)
    ke0 <- exp(lke0 + etalke0)
    ppc <- exp(lppc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc
    # CSF effect compartment holding a concentration (Figure 1, Discussion).
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    f(depot) <- exp(lfdepot)

    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
