vanEsdonk_2018_pregabalin <- function() {
  description <- "One-compartment population pharmacokinetic model for a single 300 mg oral dose of pregabalin in healthy adults, with first-order absorption after a lag time, linear elimination, and fixed allometric body-weight scaling on apparent clearance (exponent 0.75) and apparent volume of distribution (exponent 1) centred on 70 kg (van Esdonk 2018)"
  reference <- paste(
    "van Esdonk MJ, Lindeman I, Okkerse P, de Kam ML, Groeneveld GJ, Stevens J. (2018).",
    "Population pharmacokinetic/pharmacodynamic analysis of nociceptive pain models",
    "following an oral pregabalin dose administration to healthy subjects.",
    "CPT Pharmacometrics Syst Pharmacol 7(9):573-580. doi:10.1002/psp4.12318."
  )
  vignette <- "vanEsdonk_2018_pregabalin_pain_models"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling centred on 70 kg, with both exponents fixed at their theory-based values rather than estimated (van Esdonk 2018 Covariate analysis: 'allometric scaling (centered around 70 kg) was tested on the volume of distribution (Vd; exponent = 1) and clearance (CL; exponent = 0.75)'). Inclusion gave dOFV = -22.4. Cohort weight 68.0 kg (SD 8.22), range 54.25-77.50 kg (Table 1).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pregabalin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pregabalin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "19-25 years",
    age_mean = "21.75 years (SD 1.61)",
    weight_range = "54.25-77.50 kg",
    weight_mean = "68.0 kg (SD 8.22)",
    height_mean = "176 cm (SD 8.54); range 163.5-192.5 cm",
    bmi_mean = "21.89 kg/m^2 (SD 1.60); range 19.4-24.9",
    ffm_mean = "50.26 kg (SD 9.95); range 36.62-63.26 (Janmahasatian equation)",
    serum_creatinine_mean = "82.19 umol/L (SD 12.95); range 52-99",
    gfr_mean = "112.7 mL/min (SD 18.18); range 79-149 (Cockcroft-Gault)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers",
    dose_range = "Single 300 mg oral dose",
    regions = "The Netherlands (Centre for Human Drug Research, Leiden)",
    notes = "Part 2 (oral analgesics) of a two-part, four-way randomised placebo-controlled crossover study (Okkerse 2017, Br J Clin Pharmacol 83:976-990); subjects received imipramine, pregabalin, ibuprofen and placebo on separate visits with a 1-week washout. Baseline demographics from van Esdonk 2018 Table 1; the 16 subjects were evenly split between men and women (n = 8/8). PK sampling predose and at 0.5, 1, 2, 3, 4, 5, 6, 8 and 10 h after dosing; 136 of 144 planned samples were above the 20 ug/L lower limit of quantification and used for model building."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Reference subject: 70 kg.
    # All values are van Esdonk 2018 Table 2, 'Population parameters (RSE)'
    # column. The bootstrap column of the same table is reported as
    # similar throughout and is not used here.
    # Vd and CL are APPARENT (oral study, no intravenous arm), so F is not
    # separately identifiable and is left at its implicit value of 1.
    # ------------------------------------------------------------------
    lka <- log(6.07)
    label("Log of apparent first-order absorption rate constant (1/h)")
    # Table 2: 'k a (/hour)' = 6.07 [RSE 42.3%]

    ltlag <- log(0.495)
    label("Log of absorption lag time (h)")
    # Table 2: 'Lag time (hour)' = 0.495 [RSE 0.39%]

    lvc <- log(31.1)
    label("Log of apparent central volume of distribution Vd/F at 70 kg (L)")
    # Table 2: 'V d /70 kg (L)' = 31.1 [RSE 3.13%]

    lcl <- log(4.5)
    label("Log of apparent clearance CL/F at 70 kg (L/h)")
    # Table 2: 'CL/70 kg (L/hour)' = 4.5 [RSE 2.53%]

    # Allometric exponents are theory-based constants the paper imposed
    # rather than estimated (no RSE is reported for either), so both are
    # encoded with fixed().
    e_wt_cl <- fixed(0.75)
    label("Allometric body-weight exponent on apparent clearance (unitless)")
    # Covariate analysis: 'clearance (CL; exponent = 0.75)'

    e_wt_vc <- fixed(1)
    label("Allometric body-weight exponent on apparent volume of distribution (unitless)")
    # Covariate analysis: 'volume of distribution (V d ; exponent = 1)'

    # ------------------------------------------------------------------
    # Inter-individual variability. The paper states IIV was implemented
    # 'from a ln-normal distribution', and Table 2 reports the omega^2
    # values below together with a %CV column. The two are mutually
    # consistent under CV = sqrt(exp(omega^2) - 1): 2.6 -> 353%,
    # 7.09e-5 -> 0.842%, 0.0101 -> 10.1%, 0.00672 -> 8.21%, matching every
    # printed CV. The omega^2 values are therefore log-scale variances.
    # ------------------------------------------------------------------
    etalka ~ 2.6 # Table 2: 'omega 2 k a' = 2.6 (CV 353%, shrinkage 8%)
    etaltlag ~ 7.09e-5 # Table 2: 'omega 2 lag time' = 7.09E-5 (CV 0.842%, shrinkage 36%)
    etalvc ~ 0.0101 # Table 2: 'omega 2 V d / F' = 0.0101 (CV 10.1%, shrinkage 23%)
    etalcl ~ 0.00672 # Table 2: 'omega 2 CL/ F' = 0.00672 (CV 8.21%, shrinkage 21%)

    # ------------------------------------------------------------------
    # Residual error. Table 2 reports sigma^2 = 0.0146 as a VARIANCE (the
    # parenthesised 14% is shrinkage, per the row header 'Residual error
    # (shrinkage)', not a CV); nlmixr2 takes the standard deviation.
    # sqrt(0.0146) = 0.120830.
    # ------------------------------------------------------------------
    propSd <- 0.120830
    label("Proportional residual error (fraction)")
    # Table 2: 'sigma 2 proportional' = 0.0146 (shrinkage 14%)
  })

  model({
    # 1. Individual parameters with allometric weight scaling on a 70 kg
    #    reference subject.
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment ODE system with first-order absorption
    #    (Figure 1a: DEPOT -> Central -> elimination).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Absorption lag time. The paper found a lag time (dOFV = -27)
    #    superior to both a plain first-order and a transit-compartment
    #    absorption model.
    alag(depot) <- tlag

    # 5. Observation. Dose in mg and vc in L give mg/L (== ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
